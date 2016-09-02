# EMRI Kludge Suite

**Version 0.1.0**

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

Please check https://github.com/alvincjk/EMRI_Kludge_Suite for any version updates.

&mdash; Alvin Chua, Sep 2016

## Work in progress

- A fast approximate calculation of the last stable orbit, and a smooth truncation at that point, will be included in the next version of the AAK
- Higher-order polynomial fits for the AAK are being investigated, and will also be included if they provide significantly improved performance
- Suite utilities to compute waveform matches, Fisher matrices, etc. may be included in the future

## List of (important) known bugs

- A constant time offset for the AAK/AK is required to ensure far-field agreement with the NK at t=0; this is likely due to differing conventions, and needs to be fixed
- The option to output h_+/h_x instead of the LISA response functions h_I/h_II does not provide trustworthy waveforms (yet) &mdash; use at own risk
- Integer values of dt for the NK may produce NaNs at times that coincide with specific fractions of the LISA orbital period
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
