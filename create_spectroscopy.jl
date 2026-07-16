using RetrievalToolbox

using ArgParse
using JLD2
using Unitful
using YAML


function parse_commandline()

    s = ArgParseSettings(description = "Creates spectroscopy files")

    @add_arg_table! s begin

        "-i"
            help = "Path to windows YAML file"
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

    if !haskey(yml, "spectral_windows")
        @error "YAML file must contain `spectral_windows`!"
        return false
    end

    #=
    # For now, let's just allow single-band situations only
    if length(keys(yml["spectral_windows"])) != 1
        @error "The YAML file must only contain exactly one spectral window!"
        @error "Multi-window situations are not supported yet!"
        return false
    end
    =#

    use_wavelengths = false

    # Loop through every gas
    for gas in yml["gases"]

        ww_tuples = []

        # Loop through windows
        for (window, wdict) in yml["spectral_windows"]

            @info "Processing spectroscopy: $(window)"

            if haskey(wdict, "wavelength_unit")
                ww_key = "wavelength"
                use_wavelengths = true
            elseif haskey(wdict, "wavenumber_unit")
                ww_key = "wavenumber" 
            end

            # Extract the needed wavelength or wavenumber boundaries
            ww_min = wdict["$(ww_key)_min"]
            ww_max = wdict["$(ww_key)_max"]
            ww_buffer = wdict["$(ww_key)_buffer"]
            ww_unit = uparse(wdict["$(ww_key)_unit"])

            push!(ww_tuples, 
                (
                    (ww_min - 2 * ww_buffer) * ww_unit,
                    (ww_max + 2 * ww_buffer) * ww_unit
                )
            )

        end


        #@info "$(gas) from $(ww_min * ww_unit) to $(ww_max * ww_unit)"

        # Set output file name
        outname = "$(gas)_spectroscopy.jld2"

        # Calculate spectroscopy
        spectroscopy = RetrievalToolbox.create_ABSCO_from_HITRAN(
            gas,
            [ww_tuples...],
            h2o_broadening=false,
            wavelength_output=use_wavelengths
        )

        @info "Saving as $(outname)"
        jldsave("$(outname)"; spectroscopy)

    end


end


main()