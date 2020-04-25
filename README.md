# PyTracts
![](./examples/movie.gif)

### Tractography pipeline:
* Tractography pipeline assuming DWI image had already been preprocessed.
1. Generate Fractional anisotrophy (FA) image for further group analysis (for more information, see [Mrtrix3's dwi2tensor](https://mrtrix.readthedocs.io/en/latest/reference/commands/dwi2tensor.html)).
2. Estimate tissue response function to be used by the deconvolution algorithm (for more information, see [Mrtrix3 documentation](https://mrtrix.readthedocs.io/en/latest/constrained_spherical_deconvolution/response_function_estimation.html))
3. Estimate fiber orientation distributions using spherical deconvolution (for more informatin, see [Mrtrix3's dwi2fod](https://mrtrix.readthedocs.io/en/latest/reference/commands/dwi2fod.html))
4. Generate streamlines using iFOD2 algorithm (implamented by [Mrtrix3's tckgen](https://mrtrix.readthedocs.io/en/latest/reference/commands/tckgen.html)).
>Mrtrix3's documentation regarding the algortihm:
>>" iFOD2 (default): Second-order Integration over Fiber Orientation Distributions. A probabilistic algorithm that takes as input a Fiber Orientation Distribution (FOD) image represented in the Spherical Harmonic (SH) basis. Candidate streamline paths (based on short curved “arcs”) are drawn, and the underlying (trilinear-interpolated) FOD amplitudes along those arcs are sampled. A streamline is more probable to follow a path where the FOD amplitudes along that path are large; but it may also rarely traverse orientations where the FOD amplitudes are small, as long as the amplitude remains above the FOD amplitude threshold along the entire path."
* note that the algorithm used may change in futuer versions, as different algorithms perform better for different acusitions and data.
5. Converting streamlines from .tck format to the more widely used .trk format.

* note that the pipeline assumes that all images and data are organized according to recently published PyPrep packages's output.
for moer information, see [PyPrep](https://github.com/TheLabbingProject/PyPrep)

# Tracts_processing module:
Perform a complete pipeline aimed at producing streamlines from a preprocessed DWI image.

# Generate_Tracts_with_Mrtrix3:
Use mostly Mrtrix3's functions and tools to produce streamlines from a preprocessed DWI image.

> Initiate the streamlines generator class with either a specific subject or None = all subjects (default)
* Arguments:
    * mother_dir {Path} -- [Path to a PyPrep derivatives (output) directory. (should contain0 "mother_dir/sub-xx")]
* Keyword Arguments:
    * subj {str} - ["sub-xx" got specific subject or None for all subjects] (default: {None})

## Usage:
```
>>> from NetworksPlasticity.code import Tracts_processing
>>> from pathlib import Path

>>> tracts = Tracts_processing.Generate_Tracts_with_Mrtrix3(
        mother_dir = Path("/path/to/your/derivatives_dir"),
        subj = "sub-01",
        )
```
* note that the {mother_dir} directory must contain subdirectories with names such as "sub-XX"
```
>>> tracts.run()
```
### Output:
* derivatives/
    * sub-01/
        * dwi/
        * fmap/
        * anat/
        * scripts/
        * atlases/
        * func/
        * tractography/ (--> newly created directory, containing all tractography-related files)
            * FOD_wm.mif
            * response_wm.txt
            * seeds.tsv
            * FOD_gm.mif
            * dti.mif
            * fa.mif
            * tractogram.tck
            * FOD_csf.mif
            * tractogram.trk
            * response_gm.txt
            * response_csf.txt
    * sub-02/
        * dwi/
        * fmap/
        * anat/
        * scripts/
        * atlases/
        * func/
        * tractography/
    * .
    * .
    * .

