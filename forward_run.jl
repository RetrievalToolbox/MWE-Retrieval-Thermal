using RetrievalToolbox
const RE = RetrievalToolbox

using ArgParse
using JLD2
using HDF5
using Random
using Unitful
using YAML

include("forward_model.jl")


function parse_commandline()

    s = ArgParseSettings(description = "Produces forward run radiances")

    @add_arg_table! s begin

        "-i"
            help = "Path to windows YAML file"
            arg_type = String
            required = true

        "-o"
            help = "Output file name"
            arg_type = String
            required = true

    end

    return parse_args(ARGS, s)

end




function main()

    # Parse arguments
    args = parse_commandline()

    # Read retrieval window YAML 
    yml = YAML.load_file(args["i"])


    #=
        Program-level constants
        =======================

        We set values here, for ease and convenience. More involved applications
        should probably set these as either command-line options or passed via
        a configuration file.
    =#

    # Gas VMR unit lookup, so that each gas may use different units
    vmr_unit_lookup = Dict(
        "CO2" => u"ppm",
        "H2O" => Unitful.NoUnits,
    )

    # Surface pressure
    surface_pressure = uparse(yml["surface_pressure"])

    #=
        Figure out how many spectral samples our instrument will have
        =============================================================

        This is mostly taken from the YAML config file. The spectral sampling provides
        the spacing between detector samples. The instrument wavelength range is fixed
        here but could be changed by the user. Users MUST also make the same change in
        the `inverse_run.jl` script, otherwise the inconsistencies will break the 
        retrieval system.

    =#

    instrument_wavelength_start = 4000.0u"nm"
    instrument_wavelength_end = 5000.0u"nm"
    dλ = yml["instrument"]["spectral_sampling"] * uparse(yml["instrument"]["unit"])
    Nsamples = round((instrument_wavelength_end - instrument_wavelength_start) / dλ) |> Int


    #=
        Create spectral windows, spectroscopy and gases!
        ================================================
    =#

    spectral_window_list = RE.SpectralWindow[] # List of spectral windows
    ABSCO_dict = Dict{String, RE.AbstractSpectroscopy}() # Dict to hold spectroscopy
    gas_dict = Dict{String, RE.GasAbsorber}() # Dict to hold gases
    dispersion_dict = Dict{RE.SpectralWindow, RE.AbstractDispersion}() # To hold dispersions

    # Load the spectroscopy tables
    for gas_name in yml["gases"]
        absco_fname = "$(gas_name)_spectroscopy.jld2"

        @info "Loading spectroscopy $(absco_fname)"
        my_spectroscopy = jldopen(absco_fname)["spectroscopy"]

        if haskey(ABSCO_dict, gas_name)
            @warn "Gas $(gas) already in ABSCO dict!"
        end

        ABSCO_dict[gas_name] = my_spectroscopy
    end

    for (window_name, wdict) in yml["spectral_windows"]

        @info "Creating spectral window objects: $(window_name)"

        if haskey(wdict, "wavelength_unit")
            ww_key = "wavelength"
        elseif haskey(wdict, "wavenumber_unit")
            ww_key = "wavenumber" 
        end

        # Extract the needed wavelength or wavenumber boundaries
        ww_min = wdict["$(ww_key)_min"]
        ww_max = wdict["$(ww_key)_max"]
        ww_ref = wdict["$(ww_key)_ref"]
        ww_buffer = wdict["$(ww_key)_buffer"]
        ww_unit = uparse(wdict["$(ww_key)_unit"])

    
        #= Create the spectral window
           (this function creates a spectral window that uses the exact high-res
            grid as the chosen spectroscopy, which is convenient for the calculation
            of optical properties - no spectral interpolation needed there..)
        =#

        swin = RE.spectralwindow_from_ABSCO(
            window_name,
            ww_min,
            ww_max,
            ww_ref,
            ww_buffer,
            ABSCO_dict[yml["gases"][1]], # just use the first ABSCO
            ww_unit
        )

        # Add to list
        push!(spectral_window_list, swin)

    end

    # Create dispersion object
    # (characterizes the spectral sampling of the instrument, and
    #  determines how model radiances at wavelengths are mapped
    #  onto instrument spectral samples)
    for swin in spectral_window_list

        @info "Building dispersion for: $(swin)"
        this_dispersion = RE.SimplePolynomialDispersion(
            [instrument_wavelength_start, dλ],
            0:Nsamples - 1,
            swin
        )

        dispersion_dict[swin] = this_dispersion

    end

    # Create the ISRFs
    ISRF_dict = Dict{RE.AbstractDispersion, RE.AbstractISRF}()

    # Make all ISRFs the same for now..
    for (swin, disp) in dispersion_dict

        ISRF_dict[disp] = RE.GaussISRF(
            yml["instrument"]["FWHM"], # Gauss ISRF FWHM
            uparse(yml["instrument"]["unit"]) # wavelength unit
        )
        
    end

    #=
        Set the contents of the prior atmosphere
    =#
    @info "Creating and setting atmosphere"

    N_RT_lev = 72;

    # Create an Earth atmosphere object, in this case a US midwest summer one
    atmosphere = RE.create_example_atmosphere("US-midwest-summer", N_RT_lev);   

    # Set the retrieval / RT grid
    # (unit conversion done automatically)
    #RE.ingest!(
    #    atmosphere, :pressure_levels,
    #    create_UoL_pressure_grid(1000.0u"hPa", 80.0u"Pa"; N_total=N_RT_lev)
    #)

    RE.ingest!(
        atmosphere, :pressure_levels,
        atmosphere.met_pressure * atmosphere.met_pressure_unit
    )

    # Collect and create required gases
    for gas_name in yml["gases"]

        @info "Creating gas $(gas_name)"
        
        gas = create_example_gas_profile(
            "US-midwest-summer",
            gas_name,
            ABSCO_dict[gas_name],
            atmosphere.pressure_levels * atmosphere.pressure_unit
        )

        #=
        gas = RE.GasAbsorber(
            gas_name,
            ABSCO_dict[gas_name],
            zeros(N_RT_lev),
            vmr_unit_lookup[gas_name]
        )
        =#

        gas_dict[gas_name] = gas

    end


    for (gas_name, gas) in gas_dict
        @info "Adding $(gas) to atmosphere"
        push!(atmosphere.atm_elements, gas)
    end

    # We take the H2O VMR from the atmosphere's specific humidity
    # (1) sample specific humidity profile at retrieval pressure levels
    sh_profile = RE.atmospheric_profile_interpolator_linear(
        atmosphere.met_pressure,
        atmosphere.specific_humidity,
        atmosphere.pressure_levels
    )

    # (2) convert specific humidity to H2O volume mixing ratio
    h2o_profile = RE.specific_humidity_to_H2O_VMR.(sh_profile)
    # (3) insert into gas object
    gas_dict["H2O"].vmr_levels[:] .= h2o_profile

    # OPTIONAL / ALTERNATIVE
    # Set CO2 profile to some value instead of using the provided example atmosphere..
    # gas_dict["CO2"].vmr_levels[:] .= 405.0 # remember, this is in ppm

    @info "Adding thermal emission of atmosphere"
    push!(atmosphere.atm_elements, RE.ThermalAtmosphereIsotropic())

    @info "Adding thermal emission of surface"
    push!(atmosphere.atm_elements, RE.ThermalSurfaceIsotropic(300.0, u"K"))
    

    @info "Create buffer objects"
    # Need max value of all spectroscopy wavelengths
    N1 = 2 * (map(x -> length(x.ww), values(ABSCO_dict)) |> maximum)
    # Need max value of all dispersion sample lengths
    N2 = 2 * (map(x -> length(x.detector_samples), values(dispersion_dict)) |> maximum)

    instrument_buffer = RE.InstrumentBuffer(
        zeros(Float64, N1),
        zeros(Float64, N1),
        zeros(Float64, N2),
    )

    rt_buffer = RE.ScalarRTBuffer(
        dispersion_dict, # Already a SpectralWindows -> Dispersion dictionary
        RE.ScalarRadiance(Float64, N2), # Hold the radiance - we use ScalarRadiance because we don't need polarization
        nothing, # holds Jacobians, but we do not need them in the forward run..
        Dict(swin => zeros(Int, 0) for swin in spectral_window_list), # Hold the detector indices
        u"W/m^2/sr/µm" # Radiance units for the forward model, these should match the L1 observation units
    )
    
    buffer = RE.EarthAtmosphereBuffer(
        RE.RetrievalStateVector([]), # This is empty for now
        spectral_window_list,
        [(:Lambert, 2) for swin in spectral_window_list], # Use 2 coefficients for surface polynomial
        atmosphere,
        Dict(swin => RE.NoSolarModel() for swin in spectral_window_list),
        [:XRTM for swin in spectral_window_list],
        RE.ScalarRadiance, # Use ScalarRadiance for high-res radiance calculations
        rt_buffer,
        instrument_buffer,
        N_RT_lev, # The number of retrieval or RT pressure levels
        atmosphere.N_met_level, # The number of meteorological pressure levels, as given by the atmospheric inputs
        Float64 # The chosen Float data type (e.g. Float16, Float32, Float64)
    )

    # Set the instrument viewing angles
    buffer.scene.observer.viewing_zenith = yml["viewing_zenith"]
    buffer.scene.observer.viewing_azimuth = yml["viewing_azimuth"]

    # Add location details. We must create a new EarthLocation object, and then
    # insert it into the EarthScene object.

    # .. unit-aware handling of elevation
    elevation = uparse(yml["elevation"])

    new_location = RE.EarthLocation(
        yml["longitude"],
        yml["latitude"],
        ustrip(elevation),
        unit(elevation)
    )

    buffer.scene.location = new_location

    @info "Location:"
    @info "$(buffer.scene.location)"

    
    @info "Observer:"
    @info "$(buffer.scene.observer)"


    # Set the surface reflectance to some value
    for (swin, surface) in buffer.scene.surfaces
        surface.kernels[1].coefficients[1] = 0.25 # value at reference wavelength
        surface.kernels[1].coefficients[2] = 0.0001 # slope per wavelength
    end

    #=
        Settings for the XRTM RT solver
`       
        We use the `eig_bvp` solver, since it is one of the few that
        supports (thermal) emission from both the atmosphere and surface.
    =#

    rt_options = Dict()
    rt_options["solvers"] = ["eig_bvp"]
    rt_options["streams"] = 2
    rt_options["add"] = true
    rt_options["options"] = [
        "output_at_levels", 
        "source_thermal", # Use a thermal source
        #"source_solar", # Use a solar source -> not needed when wavelength > 3µm or so
        "psa", # pseudo-spherical approximation, curved Earth for incoming solar beam
        "sfi", # source function integration
    ]

    for (swin, rt) in buffer.rt
        push!(rt.model_options, rt_options)
    end

    # Set up is finished!
    # Now call the forward model to compute radiances!

    @info "Running the forward model!"

    success_fm = forward_model!(
        RE.ForwardModelStateVector();
        buffer=buffer,
        ISRF_dict=ISRF_dict
    )

    if !success_fm
        @error "Forward model did not execute correctly! Aborting."
        exit()
    end

    # Calculate random noise
    # Set the seed here, so the results are
    # reproducible, at least for the same run..
    rng = Xoshiro(12345)

    noise_dict = Dict{RE.SpectralWindow, Vector{Float64}}()

    for swin in spectral_window_list

        # Our noise model is:
        #   noise = √(a*I + b)

        noise_a = yml["spectral_windows"][swin.window_name]["noise_a"]
        noise_b = yml["spectral_windows"][swin.window_name]["noise_b"]

        # This is a simple SNR model, where noise is a normal distribution
        # sample with sigma = signal / SNR. You may want to replace this
        # by an actual noise model.
        signal = rt_buffer.radiance.I[rt_buffer.indices[swin]]

        noise = @. sqrt(signal * noise_a + noise_b)
        noise .*= randn(rng, Float64, (length(noise)))

        # Add noise to radiance
        rt_buffer.radiance.I[rt_buffer.indices[swin]] .+= noise

        # Store to be written out later, so you could do "radiance - noise" to
        # get the raw noise-free radiance..
        noise_dict[swin] = noise

    end



    @info "Writing radiances into file!"

    # Save out the radiances:
    h5open(args["o"], "w") do fid

        for swin in spectral_window_list

            @info "Writing out radiances for $(swin)"
            disp = dispersion_dict[swin]
            
            g = create_group(fid, swin.window_name)
            
            # Write radiance
            output_radiance = zeros(length(disp.detector_samples)) * NaN
            output_radiance[disp.index] .= rt_buffer.radiance.I[rt_buffer.indices[swin]]
            
            dset = create_dataset(g, "radiance", Float64, (length(output_radiance),))
            write(dset, output_radiance)

            # Write noise
            output_noise = zeros(length(disp.detector_samples)) * NaN
            output_noise[disp.index] .= noise_dict[swin]

            dset = create_dataset(g, "noise", Float64, (length(output_noise),))
            write(dset, output_noise)

            # Write spectral
            output_ww = zeros(length(disp.detector_samples)) * NaN
            output_ww[disp.index] .= disp.ww

            dset = create_dataset(g, "wavelength", Float64, (length(output_ww),))
            write(dset, output_ww)

        end

    end


end 


main()