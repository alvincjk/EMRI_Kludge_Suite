# EMRI Kludge Suite

**Version 0.2.0**

This is a C/C++ suite that allows kludge waveforms for extreme-mass-ratio inspirals (EMRIs) to be generated with shared settings and parameters. The three waveforms included in the suite are the augmented analytic kludge (AAK) [1,2], the analytic kludge (AK) [3], and the numerical kludge (NK) [4]. EMRI Kludge Suite is part of the Black Hole Perturbation Toolkit; visit http://bhptoolkit.org for more information.

The GSL and FFTW libraries are required for compilation. Running `make` will create the corresponding executables in the folder `./bin`:

- `AAK_Waveform`
- `AK_Waveform`
- `NK_Waveform`

The file `./SetPar_Template` is a template for a formatted settings/parameters file that contains the output file path, various toggles, and the EMRI parameters. More details are provided in the template file itself.

As an example, running

`bin/AAK_Waveform SetPar_Template`

will generate an AAK waveform with default settings and parameters. Three files will be created in `./bin`:

- `example_wave.dat` contains waveform data (t, h_I, h_II)
- `example_traj.dat` contains inspiral trajectory data (t, p/M, e, iota, E, L_z, Q)
- `example_info.txt` contains additional information such as signal-to-noise ratio and waveform timing

NEW IN VERSION 0.2.0: Improvements to AAK, e.g. automated backward integration for plunging orbits; adaptive fitting duration; better speed and robustness.

Please check https://github.com/alvincjk/EMRI_Kludge_Suite for any version updates.

&mdash; Alvin Chua, Feb 2018

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

The EMRI Kludge Suite is also based on code written by Leor Barack (for the AK) and Scott Hughes (for the NK).

## References

[1] A. J. K. Chua & J. R. Gair. Improved analytic extreme-mass-ratio inspiral model for scoping out eLISA data analysis. *Class. Quantum Grav.* 32:232002, 2015.

[2] A. J. K. Chua, C. J. Moore & J. R. Gair. Augmented kludge waveforms for detecting extreme-mass-ratio inspirals. *Physical Review D* 96:044005, 2017.

[3] L. Barack & C. Cutler. LISA capture sources: Approximate waveforms, signal-to-noise ratios, and parameter estimation accuracy. *Physical Review D* 69:082005, 2004.

[4] S. Babak, H. Fang, J. R. Gair, K. Glampedakis & S. A. Hughes. "Kludge" gravitational waveforms for a test-body orbiting a Kerr black hole. *Physical Review D* 75:024005, 2007.
