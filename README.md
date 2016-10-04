# EMRI Kludge Suite

**Version 0.1.2**

This is a C/C++ suite that allows kludge waveforms for extreme-mass-ratio inspirals (EMRIs) to be generated with shared settings and parameters. The three waveforms included in the suite are the augmented analytic kludge (AAK) [1], the analytic kludge (AK) [2], and the numerical kludge (NK) [3].

Running `make` will create the corresponding executables in the folder `./bin`:

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

NEW IN VERSION 0.1.1: The AAK implementation now uses a 10-point quartic polynomial fit, which provides high overlaps with the NK over longer durations [4]. The time duration spanned by the fitting points may be specified by the new parameter T_fit (default value is two weeks).

NEW IN VERSION 0.1.2: A fast approximate calculation of the last stable orbit has been implemented for the AAK, allowing the waveform to automatically include plunge [4]. The trajectory file is truncated at the last stable orbit, while a one-sided Planck-taper window is used to zero the waveform smoothly over 10 additional orbits.

Please check https://github.com/alvincjk/EMRI_Kludge_Suite for any version updates.

&mdash; Alvin Chua, Oct 2016

## Work in progress

- Suite utilities to compute waveform matches, Fisher matrices, etc. may be included in the future

## List of (important) known bugs

- A constant time offset for the AAK/AK is required to ensure far-field agreement with the NK at t=0; this is likely due to differing conventions, and needs to be fixed
- The option to output h_+/h_x instead of the LISA response functions h_I/h_II does not provide trustworthy waveforms (yet) &mdash; use at own risk
- The NK may produce NaNs at times that coincide with specific fractions of the LISA orbital period (when using integer values of dt), or near plunge; this needs to be fixed, but in the meantime a band-aid solution for dealing with isolated NaNs is included
- All waveforms may have difficulties with zero or epsilon values for certain parameters such as spin, eccentricity and inclination

## Authors

**Alvin Chua**  
Institute of Astronomy, University of Cambridge  
`ajkc3@ast.cam.ac.uk`

**Jonathan Gair**  
School of Mathematics, University of Edinburgh  
`j.gair@ed.ac.uk`

The EMRI Kludge Suite is also based on code written by Leor Barack (for the AK) and Scott Hughes (for the NK).

## References

[1] A. J. K. Chua & J. R. Gair. Improved analytic extreme-mass-ratio inspiral model for scoping out eLISA data analysis. *Class. Quantum Grav.* 32:232002, 2015.

[2] L. Barack & C. Cutler. LISA capture sources: Approximate waveforms, signal-to-noise ratios, and parameter estimation accuracy. *Physical Review D* 69:082005, 2004.

[3] S. Babak, H. Fang, J. R. Gair, K. Glampedakis & S. A. Hughes. "Kludge" gravitational waveforms for a test-body orbiting a Kerr black hole. *Physical Review D* 75:024005, 2007.

[4] A. J. K. Chua, C. J. Moore & J. R. Gair. The fast and the fiducial: Improving kludge waveforms for extreme-mass-ratio inspirals. *In prep.*
