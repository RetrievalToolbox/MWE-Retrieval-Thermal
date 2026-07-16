function forward_model!(
    state_vector::RE.AbstractStateVector;
    buffer::RE.EarthAtmosphereBuffer,
    ISRF_dict::Dict{RE.AbstractDispersion, RE.AbstractISRF}
)

    # Handy short-cuts
    rt_buf = buffer.rt_buf
    inst_buf = buffer.inst_buf

    # Update the dispersions object based on the state vector
    # (You must make sure that the dispersion state vector elements
    #  are bound to the used dispersion!)
    # Dispersions are stored in the RT buffer as `SpectralWindow`->Dispersion dicts
    if state_vector isa RE.RetrievalStateVector
        for (swin, disp) in rt_buf.dispersion
            RE.update_dispersion!(disp, state_vector)
        end
    end

    # Produce the correct indices (mappings between dispersion and buffer positions)
    RE.calculate_indices!(buffer)

    # If we have a surface state vector element, plug in the current value
    if state_vector isa RE.RetrievalStateVector
        for (swin, rt) in buffer.rt
            RE.surfaces_statevector_update!(rt.scene, state_vector)
        end
    end

    # Ingest retrieved H2O into specific humidity!
    # RE.update_specific_humidity_from_H2O!(buffer.scene.atmosphere)

    # Update the atmosphere with state vector contents (e.g. temperature scale)
    if state_vector isa RE.RetrievalStateVector
        RE.atmosphere_statevector_update!(buffer.scene.atmosphere, state_vector)
    end

    # Given the current statevector - modify any atmospheric element
    # according to whatever the state vector element dictates.
    # (e.g. scale the VMR of a gas according to the state vector)
    if state_vector isa RE.RetrievalStateVector
        RE.atmosphere_element_statevector_update!(
            buffer.scene.atmosphere.atm_elements, 
            state_vector
        )
    end

    # Update the solar scaler field in the RT container
    for rt in values(buffer.rt)
        RE.solar_scaler_statevector_update!(rt, state_vector)
    end

    # Calculate z and g profiles
    RE.calculate_altitude_and_gravity!(buffer.scene)

    # Calculate optical properties for every spectral window
    for (swin, rt) in buffer.rt
        RE.calculate_earth_optical_properties!(
            rt.optical_properties,
            buffer.scene,
            state_vector
        )
    end

    # RT for every spectral window
    for (swin, rt) in buffer.rt
        
        # Clear out RT results
        RE.clear!(rt)

        # Calculate the TOA radiance with XRTM
        RE.calculate_radiances_and_jacobians!(rt)

    end


    # Apply the ISRF to the hi-res spectrum to obtain detector-level radiance
    for (swin, rt) in buffer.rt

        disp = rt_buf.dispersion[swin]

        # Let's make sure we actually have the correct ISRF in the dict
        if !haskey(ISRF_dict, disp)
            @error "ISRF dict has no entry for dispersion: $(disp)"
            @error "Exiting forward model!"
            return false
        end

        success = RE.apply_isrf_to_spectrum!(
            inst_buf, # Convolution buffer
            ISRF_dict[disp], # ISRF table
            disp, # Dispersion
            rt.hires_radiance.I, # Radiance
            swin # Spectral window
        )

        if !(success)
            @debug "Application of ISRF on radiances FAILED for window $(swin)."
            return false
        else
            @debug "Application of ISRF on radiances successful for window $(swin)."
        end

        # Store in RT buffer, but account for unit differences!
        @views rt_buf.radiance.I[rt_buf.indices[swin]] = inst_buf.low_res_output[disp.index] *
            rt.radiance_unit / rt_buf.radiance_unit

    end

    # Apply the ISRF to the hires Jacobians
    if state_vector isa RE.RetrievalStateVector
        for (i, sve) in enumerate(state_vector.state_vector_elements)
            if RE.calculate_jacobian_before_isrf(sve)

                for (swin, rt) in buffer.rt

                    disp = rt_buf.dispersion[swin]

                    success = RE.apply_isrf_to_spectrum!(
                        inst_buf,
                        ISRF_dict[disp],
                        disp,
                        rt.hires_jacobians[sve].I,
                        swin
                    )

                    if !success
                        @debug "Application of ISRF on Jacobian $(i) for $(sve) FAILED."
                        return false
                    else
                        @debug "Application of ISRF on Jacobian $(i) for $(sve) successful."
                    end

                    @views rt_buf.jacobians[sve].I[rt_buf.indices[swin]] = inst_buf.low_res_output[disp.index] *
                        rt.radiance_unit / rt_buf.radiance_unit

                end
            end
        end    
    end
    


    # Important!
    # Re-set the atmosphere to its prior state!
    if state_vector isa RE.RetrievalStateVector
        RE.atmosphere_element_statevector_rollback!(buffer.scene.atmosphere.atm_elements, state_vector)
        RE.atmosphere_statevector_rollback!(buffer.scene.atmosphere, state_vector)
    end

    @info "FM done!"

    # Must always return true if everything worked!
    return true

end