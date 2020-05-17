GENERATE_FA = """Fractional anisotrophy image generator.
    Working directory: {subj_dir}
    Inputs:
        - Preprocessed DWI: {dwi_name}
        - Mask: {mask_name}
    Outputs:
        - Output-containing directory: {tracts_dir}
        - Fractional anisotrophy (FA) file: {FA_name}
    """

GENERATE_RESPONSES = """Tissue response function estimator.
    Working directory: {subj_dir}
    Inputs:
        - Preprocessed DWI: {dwi_name}
        - Mask: {mask_name}
    Outputs:
        - Outputs-containing directory: {tracts_dir}
        - Tissue response function files: ["response_wm.txt", "response_gm.txt", "response_csf.txt"]
    """

CALCULATE_FIBER_ORIENTATION = """Fibre orientation distributions estimator.
    Working directory: {subj_dir}
    Inputs:
        - Preprocessed DWI: {dwi_name}
        - Mask: {mask_name}
    Outputs:
        - Outputs-containing directory: {tracts_dir}
        - Fibre orientation distributions (FOD) files: ["FOD_wm.mif", "FOD_gm.txt", "FOD_csf.txt"]
    """

GENERATE_TCK = """Tracts generator parameters:
    working directory: {subj_dir}
    Inputs:
        - White matter fiber orientation distributions file: {fod_wm_name} (at {tracts_dir})
        - Five-tissue-type: {seg_5tt_name} (at "dwi/Mrtrix_prep")
        - number of analyzed tracts: 350000
    Outputs:
        - Outputs-contianing directory: {tracts_dir}
        - Tractogram (.tck) file: {tractogram_name}
        - Seeds file: seeds.csv
    """

CONVERT_TCK_2_TRK = """tck to trk converter.
    Working directory: {subj_dir}
    Inputs:
        - Preprocessed DWI image: {dwi_name} (at "dwi" subdirectory)
        - tck format streamlines: {tck_name} (at "tractography" subdirectory)
    Outputs:
        - trk format streamlines: {trk_name}
    """

TRACTS_WITH_MRTRIX3 = """Mrtrix3-based tractography generator
    Working ("mother") dir: {mother_dir}
    Subjects: {subjects}
    """

MRTRIX_PRINT_START = """Initializing tracts processing for {subj}...
    Initial input files were extracted from subject`s directory at: {folder_name}
        Preprocessed DWI image: {dwi_name} (at "dwi" subdirectory)
        Brain mask: {dwi_mask} (at "dwi" subdirectory)
        Five-tissue-type (5TT): {five_tissue} (at "anat" subdirectory)
    Output files file will be located at "{tracts_dir}" subdirectory under subject's directory."""

GENERATE_WHITE_MATTER = """"Resampling subject's high-resolution segmented image to DWI space...
    Inputs:
        - High-resolution segmented image: {labels_file}
        - DWI image: {dwi_file}
    Outputs:
        - nibabel.Nifti1Image of resampled labels images in DWI space."""

LOAD_DATA = """Loading data from given file names to fit Dipy necessities.
    Inputs:
        - Path to DWI file: {dwi_file}
        - Path to .bvec file: {bvec}
        - Path to .bval file: {bval}
    Outputs:
        - DWI image's data, affine (as numpy.ndarrays), and nib.Nifti1Image
        - Dipy's gradient table
"""

CREATE_MODEL = """Create Constrained Spherical Deconvolution (CSD) and Constant Solid Angle (CSA) models to fit the data.
    Inputs:
        - data: DWI nump.ndarray to model.
        - gtab: Dipy's gradient table.
        - white_matter: mask of white matter.
        - sh_order,roi_radius,fa_thr,relatie_peak_threshold,min_seperation_angle: numerous parameter for model generation.
        {sh},{roi},{fa},{relative},{min_separation} accordingly.
    Outputs:
        - CSD and CSA models that fit the data.
"""

GENERATE_STREAMLINES = """Generate tractogram ("streamlines") based on fitted CSD model.
    Inputs:
        - Directory in which the tractogram will be saved: {folder}.
        - Fitted CSD model.
        - DWI image.
        - Fiber tracking (reconstruction) method: {recon}.
        - Stopping criterion.
        - Seed mask.
        - Image's affine matrix.
        - Tractogram's file name: {fname}.
        - Sphere type, maximum fiber's angle, step size: {sphere},{angle},{step} accordingly.
    Outputs:
        - Tractogram file located at {tracts_loc}
"""
