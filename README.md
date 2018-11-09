# EMRI Kludge Suite

**Version 0.3.1**

This is a C/C++ suite that allows kludge waveforms for extreme-mass-ratio inspirals (EMRIs) to be generated with shared settings and parameters. The three waveforms included in the suite are the augmented analytic kludge (AAK) [1,2], the analytic kludge (AK) [3], and the numerical kludge (NK) [4]. EMRI Kludge Suite is part of the Black Hole Perturbation Toolkit; visit http://bhptoolkit.org for more information.

The GSL and FFTW libraries are required for compilation. Clean up any previous installation first with `make clean`. Running `make` will create the following executables in the folder `./bin`:

- `AK_Waveform`
- `NK_Waveform`
- `AAK_Waveform`
- `AAK_Phase`

The file `./examples/SetPar_Waveform` is a template for a formatted settings/parameters file that contains the output file path, various toggles, and the EMRI parameters. More details are provided in the template file itself.

As an example, running `bin/AAK_Waveform examples/SetPar_Waveform` will generate an AAK waveform with default settings and parameters. Three files will be created in `./bin`:

- `example_wave.dat` contains waveform data (t, h_I, h_II)
- `example_traj.dat` contains inspiral trajectory data (t, p/M, e, iota, E, L_z, Q)
- `example_info.txt` contains additional information such as signal-to-noise ratio and waveform timing

The other template file `./examples/SetPar_Phase` contains default settings and parameters for the executable `./bin/AAK_Phase`, which computes the evolution of the radial, polar and azimuthal phases in the AAK (without amplitude information). These phases are fast to generate and can be downsampled significantly; their time derivatives are the Kerr fundamental frequencies (see [2] for their explicit relation to the AK frequencies).

Running `bin/AAK_Phase examples/SetPar_Phase` will create two files in `./bin`:

- `example_wave.dat` contains phase data (t, phase_r, phase_theta, phase_phi, omega_r, omega_theta, omega_phi, eccentricity)
- `example_info.txt` contains timing information

Python support is also available for the AAK waveform and phases. The `AAKwrapper` module is installed by running `python setup.py install`; see the file `./examples/AAKdemo.py` for example usage.

NEW IN VERSION 0.3.1: Added eccentricity output to phase executable.

NEW IN VERSION 0.3.0: Executable for fast phase/frequency generation; Python wrapper and demo.

Please check https://github.com/alvincjk/EMRI_Kludge_Suite for any version updates.

&mdash; Alvin Chua, Nov 2018

## Work in progress

- Suite utilities to compute waveform matches, Fisher matrices, etc. may be included in the future

## List of (important) known bugs

- The LISA response functions h_I/h_II for the AAK/AK and the NK do not match up, likely due to differing conventions when implementing the Doppler shift
- The NK may produce NaNs at times that coincide with specific fractions of the LISA orbital period (when using integer values of dt), or near plunge; this needs to be fixed, but in the meantime a band-aid solution for dealing with isolated NaNs is included
- All waveforms may have difficulties with zero or epsilon values for certain parameters such as spin, eccentricity and inclination

## Authors

**Alvin Chua**  
Jet Propulsion Laboratory  
`alvin.j.chua@jpl.nasa.gov`

**Jonathan Gair**  
School of Mathematics, University of Edinburgh  
`j.gair@ed.ac.uk`

The EMRI Kludge Suite is also based on code written by Leor Barack (for the AK) and Scott Hughes (for the NK). The Python wrapper for the AAK is provided by Michele Vallisneri.

## References

[1] A. J. K. Chua & J. R. Gair. Improved analytic extreme-mass-ratio inspiral model for scoping out eLISA data analysis. *Class. Quantum Grav.* 32:232002, 2015.

[2] A. J. K. Chua, C. J. Moore & J. R. Gair. Augmented kludge waveforms for detecting extreme-mass-ratio inspirals. *Physical Review D* 96:044005, 2017.

[3] L. Barack & C. Cutler. LISA capture sources: Approximate waveforms, signal-to-noise ratios, and parameter estimation accuracy. *Physical Review D* 69:082005, 2004.

[4] S. Babak, H. Fang, J. R. Gair, K. Glampedakis & S. A. Hughes. "Kludge" gravitational waveforms for a test-body orbiting a Kerr black hole. *Physical Review D* 75:024005, 2007.
