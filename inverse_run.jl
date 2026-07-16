using RetrievalToolbox
const RE = RetrievalToolbox

using ArgParse
using JLD2
using HDF5
using LinearAlgebra
using Statistics
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

        "-m"
            help = "Path to synthetic measurement file"
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
        a configuration file. Some of these values you may want to set differently
        compared to the forward run, or keep them the same.
    =#

    # Number of pressure levels for the radiative transfer
    N_RT_lev = 72

    # Gas VMR unit lookup
    vmr_unit_lookup = Dict(
        "CO2" => u"ppm",
        "H2O" => Unitful.NoUnits,
    )

    # Surface pressure
    surface_pressure = 1000.0u"hPa"

    # Figure out how many spectral samples our instrument will have
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
    # (characterizes the spectral sampling of the instrument)
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
        
        #=
        gas = RE.GasAbsorber(
            gas_name,
            ABSCO_dict[gas_name],
            zeros(N_RT_lev),
            vmr_unit_lookup[gas_name]
        )
        =#
        gas = create_example_gas_profile(
            "US-midwest-summer",
            gas_name,
            ABSCO_dict[gas_name],
            atmosphere.pressure_levels * atmosphere.pressure_unit
        )        

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

    # Set CO2 profile to some value
    # NOTE! YOU MAY WANT TO SET THIS TO A DIFFERENT VALUE COMPARED TO THE FORWARD RUN
    # gas_dict["CO2"].vmr_levels[:] .= 405.0 # remember, this is in ppm

    @info "Adding thermal emission of atmosphere"
    push!(atmosphere.atm_elements, RE.ThermalAtmosphereIsotropic())

    # NOTE! YOU MAY WANT TO SET THIS TO A DIFFERENT VALUE COMPARED TO THE FORWARD RUN
    @info "Adding thermal emission of surface"
    push!(atmosphere.atm_elements, RE.ThermalSurfaceIsotropic(300.0, u"K"))
    

    #=

        Setting up the state vector
        ===========================

        Here we decide which elements of the forward model we want to adjust,
        as well as the details (prior value, prior covariance, ..) of each
        state vector element.

        Adjust this section manually, or modify to read parameters from the
        configuration file.
    =#


    # Adjusts the temperature profile by adding a constant offset. We
    # thus retrieve the 
    DeltaT_profile_SVE = RE.TemperatureOffsetSVE(
        u"K", # unit
        0.0, # first guess
        0.0, # prior value
        10.0, # prior covariance
    )

    surfaceT_SVE = RE.SurfaceTemperatureSVE(
        u"K",
        300.0, # first guess
        300.0, # prior value
        10.0, # prior covariance
    )


    # Gas SVEs for every present gas:
    gas_SVEs = Vector{RE.GasLevelScalingFactorSVE}()
    for gas_name in yml["gases"]

        this_gas = RE.get_gas_from_name(atmosphere, gas_name)

        this_gas_SVE = RE.GasLevelScalingFactorSVE(
            1, # Start level index
            N_RT_lev, # End level index
            this_gas,
            Unitful.NoUnits,
            1.0, # first guess
            1.0, # prior value
            1.0 # prior covariance
        )

        push!(gas_SVEs, this_gas_SVE)

    end

    # Retrieve the surface reflectance..
    surf_SVEs = Vector{RE.BRDFPolynomialSVE}()
    for swin in spectral_window_list

        for order in [0,1]

            # Set the first guess according to coefficient order
            # NOTE! YOU MAY WANT TO SET THIS TO A DIFFERENT VALUE COMPARED TO THE FORWARD RUN
            if order == 0
                fg = 0.25
            else
                fg = 0.0001
            end

            this_surf_SVE = RE.BRDFPolynomialSVE(
                swin,
                RE.LambertianPolynomialKernel,
                order,
                u"nm",
                fg, # first guess
                fg, # prior = first guess here
                1.0^(-2*order)
            )

            push!(surf_SVEs, this_surf_SVE)
        
        end

    end    



    # Collect all SVEs and construct the state vector
    state_vector = RE.RetrievalStateVector([
        DeltaT_profile_SVE,
        surfaceT_SVE,
        gas_SVEs..., # unpack this list and add..
        surf_SVEs..., # unpack this list and add..
    ])



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
        Dict(sve => RE.ScalarRadiance(Float64, N2) for sve in state_vector.state_vector_elements), # holds Jacobians
        Dict(swin => zeros(Int, 0) for swin in spectral_window_list), # Hold the detector indices
        u"W/m^2/sr/µm" # Radiance units for the forward model, these should match the L1 observation units
    )
    
    buffer = RE.EarthAtmosphereBuffer(
        state_vector,
        spectral_window_list,
        [(:Lambert, 2) for swin in spectral_window_list],
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

    @info "Observer:"
    @info "$(buffer.scene.observer)"    

    # Set the surface reflectance to some value
    for (swin, surface) in buffer.scene.surfaces
        surface.kernels[1].coefficients[1] = 0.25 # value at reference wavelength
        surface.kernels[1].coefficients[2] = 0.0 # slope per wavelength
    end

    rt_options = Dict()
    rt_options["solvers"] = ["eig_bvp"]
    rt_options["streams"] = 2
    rt_options["add"] = true
    rt_options["options"] = [
        "output_at_levels",
        "calc_derivs", # CALCULATE JACOBIANS!
        "source_thermal", # Use a thermal source
        # "source_solar", # Use a solar source -> not needed when wavelength > 3µm or so
        "psa", # pseudo-spherical approximation, curved Earth for incoming solar beam
        "sfi", # source function integration
    ]

    for (swin, rt) in buffer.rt
        push!(rt.model_options, rt_options)
    end

    # Read in the radiance from the synthetic radiance file
    # created by the forward run:

    meas_dict = Dict{RE.AbstractDispersion, Vector{Float64}}()
    noise_dict = Dict{RE.AbstractDispersion, Vector{Float64}}()
    
    @info "Reading synthetic measurements"
    RE.calculate_indices!(buffer)

    h5open(args["m"], "r") do fid

        for (swin, disp) in dispersion_dict
            meas_dict[disp] = zeros(length(disp.detector_samples))
            meas_dict[disp][disp.index] = fid["$(swin.window_name)/radiance"][:]

            noise_dict[disp] = zeros(length(disp.detector_samples))
            noise_dict[disp][disp.index] = fid["$(swin.window_name)/noise"][:]
        end

    end

    # The following keyword-argument list is needed to run iterations'
    # of the inverse solver. Expand this, as well as the forward_model!
    # function signature if your forward model needs more or different data!
    fm_kwargs = (
        buffer=buffer,
        ISRF_dict=ISRF_dict,
    )


    # Create solver, this object managers the linear algebra and
    # is needed to find a solution.
    @info "Creating IMAP solver object"
    solver = RE.IMAPSolver(
        forward_model!,
        state_vector,
        Diagonal(RE.get_prior_covariance(state_vector)), # Prior covariance matrix - just use diagonal
        20, # max number of iterations
        1, # dsigma scale
        dispersion_dict,
        rt_buffer.indices,
        rt_buffer.radiance,
        rt_buffer.jacobians,
        meas_dict, # the measurement (full spectrometer)
        noise_dict # noise (full spectrometer)
    )

    iter_result = true
    converged = false

    # ITERATE!
    while !(RE.check_convergence(solver)) # Loop until converged

        @info "Doing iteration.."
        # Evaluate and produce next iteration state vector!
        iter_result = RE.next_iteration!(solver; fm_kwargs)
        
        RE.print_state_vector_update(state_vector)

        if !(iter_result)
            @error "Error in the forward model evaluation! Exiting."
            exit()
        end

    end

    # Update the atmosphere object
    RE.atmosphere_statevector_update!(
        buffer.scene.atmosphere,
        state_vector
    )

    # Update the atmosphere elements (this must be done separately)
    RE.atmosphere_element_statevector_update!(
        buffer.scene.atmosphere.atm_elements,
        state_vector
    )

    # Print out the fit quality
    chi2_all = RE.calculate_chi2(solver)
    @info "Fit quality:"
    for (swin, chi2) in chi2_all
        @info "$(swin): $(chi2)"
    end 

    # Print out the posterior
    RE.print_posterior(solver)


    # Store the radiances as output
    @info "Writing radiances into file!"

    # Save out the fit:
    h5open(args["o"], "w") do fid

        for swin in spectral_window_list

            @info "Writing out radiances for $(swin)"
            disp = dispersion_dict[swin]
            
            g = create_group(fid, swin.window_name)
            
            modeled = RE.get_modeled(solver, swin; view=false) # take a copy of the model radiance
            dset = create_dataset(g, "modeled", Float64, (length(modeled),))
            write(dset, modeled)

            measured = RE.get_measured(solver, swin; view=false) # take a copy of the model radiance
            dset = create_dataset(g, "measured", Float64, (length(measured),))
            write(dset, measured)

            dset = create_dataset(g, "wavelength", Float64, (length(disp.ww),))
            write(dset, disp.ww)

        end

    end

end

main()