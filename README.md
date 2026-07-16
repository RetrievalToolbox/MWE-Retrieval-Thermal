# Minimum working example of a thermal CO2 and H2O retrieval

This example application is a fully self-contained demo that illustrates how a retrieval algorithm can be built with [RetrievalToolbox](https://github.com/US-GHG-Center/RetrievalToolbox.jl).

## Introduction

This code example serves as a bare-bones implementation of a trace gas retrieval for CO$_2$ and H$_2$O in thermal wavelength region between 4 µm and 5 µm. It contains two seperate scripts: one that runs the forward model in "forward mode", to produce synthetic observations; and a second script that runs a retrieval on those to adjust a small number of parameters and finds the best fit to the generated synthetic observations.

The observer is set as a space-based observer looking down at some user-specified viewing angle. Due to the wavelength ranges considered here, there is no solar model needed (emission spectroscopy) and thus also no solar geometry. No atmospheric scattering due to aerosols is considered either.

Error checking and input sanitization is kept to a minimum to allow users to look at the most necessary architecture of a RetrievalToolbox retrieval set-up.


## Requirements and Installation

RetrievalToolbox along with the XRTM library must be able to run on your computer, which necessitates a working Julia installation. For details see the documentation of RetrievalToolbox [here](https://github.com/US-GHG-Center/RetrievalToolbox.jl#installation), and pay attention to how XRTM must be compiled [here](https://github.com/US-GHG-Center/RetrievalToolbox.jl#building-xrtm-and-making-retrievaltoolbox-aware-of-its-location). RetrievalToolbox's internal RT solver does not (yet) support thermal emission from surface or the atmosphere.

To install the dependencies, run following line inside the cloned directory:

    julia --project=. -e 'using Pkg; Pkg.instantiate();'

## Set up the window file

Either directly edit `example_windows.yml` YAML file or make a copy and maniuplate the copy. Choose the bounds of the spectral window(s), their reference wavelengths or the spectral window names etc. We recommend not making the spectral window too large, as the underlying spectral grid will grow and the RT calculations may take too long. The high-resolution spectral grid is based on a constant spacing in wavenumbers (0.01 cm-1), so the spacing will not be regular in wavelength space.

For a first attempt, we recommend not changing the windows file at all, but the entries are easily understood:

    spectral_windows:

        window1:
            wavelength_unit: nm
            wavelength_min: 4350.0
            wavelength_max: 4500.0
            wavelength_ref: 4500.0
            wavelength_buffer: 2.0
            SNR: 100

        window2:
            ...

    gases:
        - CO2
        - H2O

    instrument:
        spectral_sampling: 0.02
        FWHM: 0.25
        unit: nm

    viewing_zenith: 15.0
    viewing_azimuth: 0.0


Users may add additional spectral windows under the corresponding section, as long as they have a unique name and the required sub-entries are set.

## Creating spectroscopy files

First, we need to have the spectroscopy tables ready. To create those, run the `create_spectroscopy.jl` script, by pointing it to the window file:

    julia --project=. create_spectroscopy.jl -i example_windows.yml

Change the command line argument `-i` if you want to use a different window file. This launches the calculation of spectroscopy tables based on the HITRAN database, and using the [HITRAN.jl](https://github.com/TacHawkes/HITRAN.jl) package. For every gas present in the windows, a seperate spectroscopy table will be created.

## Running in forward mode

Once the spectroscopy is created, we must run the forward script which generates the synthetic observations:

    julia --project=. forward_run.jl -i example_windows.yml -o synth_obs.h5

The `-o` commandline argument defines the location of the newly created file which contains the radiances of the synthetic observation.

## Running the retrieval

After synthetic observations are created, we can run the inverse script, which adjusts a small number of parameters to fit the (synthetic) observations:

    julia --project=. inverse_run.jl -i example_windows.yml -m synth_obs.h5 -o fit_results.h5

If no edits are made to the files, the retrieval will do very little, as the true values of the forward run are used as the initial values, and thus the fit will converge after two forward model evaluations. The program will print out some details about the posterior state vector and then quit. The radiance fits are stored out into the `fit_results.h5` file (or whichever filename the user supplied via the `-o` argument).

The un-modified script should re-produce the following. The atmospheric temperature offset should be ~0 K, the surface temperature 300 K. Both gas scaling factors for CO$_2$ and H$_2$O should be ~1.0, and the surface albedo coefficients 0.25 and 0.001.

Should users change the contents of the driving YAML configuration, both forward and inverse must be run for changes to take effect, unless a mis-matched situation is desired. If the spectral windows themselves are changed (extent, number, ..), the spectroscopy creation may also need to be re-run.

## Known Issues

- Currently, H2O (water vapor) spectroscopy tables in RetrievalToolbox, as generated via the HITRAN.jl interface, is buggy. This is a RetrievalToolbox-internal bug. For a simple retrieval exercise and means of learning RetrievalToolbox, this is still workable since both forward and inverse parts use the same spectroscopy. Alternatively, users may modify the code to read in trustworthy spectroscopy tables instead.