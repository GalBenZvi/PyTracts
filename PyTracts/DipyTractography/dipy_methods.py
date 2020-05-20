from dipy.core.gradients import gradient_table
from dipy.data import get_fnames
from dipy.io.gradients import read_bvals_bvecs
from dipy.io.image import load_nifti, load_nifti_data
from pathlib import Path
import nibabel as nib
from nipype.interfaces import fsl
from nilearn.image import resample_to_img
import numpy as np

FSLOUTTYPE = ".nii.gz"


def load_dwi_files(folder_name: Path):
    """
    Load relevant files for subject's preprocessed derivatives folder
    Arguments:
        folder_name {Path} -- [Path to subjects]

    Returns:
        dwi_file {Path} -- [Path to subject's preprocessed DWI image]
        bvec_file {Path} -- [Path to DWI's .bvec file]
        bval_file {Path} -- [Path to DWI's .bval file]
        labels_file {Path} -- [Path to subject's segmentation file (at structural space)]
    """
    dwi_folder = folder_name / "dwi"
    for file in dwi_folder.iterdir():
        if file.suffix == ".bvec":
            bvec_file = file
        elif "preprocessed_biascorr.nii.gz" in str(file):
            dwi_file = file
    labels_file = folder_name / "anat" / "prep.anat" / "T1_fast_seg.nii.gz"
    bval_file = bvec_file.with_suffix(".bval")
    return dwi_file, labels_file, bvec_file, bval_file


def resample_to_dwi(labels_file: Path, dwi_file: Path, out_file: Path):
    flt = fsl.FLIRT()
    flt.inputs.in_file = labels_file
    flt.inputs.reference = dwi_file
    flt.inputs.out_file = out_file
    if not out_file.exists():
        flt.run()
    resampled_img = round_seg(nib.load(out_file))
    return resampled_img


def round_seg(img):
    """
    Round image's data (to keep segmentation values)
    Arguments:
        img {[nibabel.Nifti1Image]} -- [image to be rounded]

    Returns:
        [new_img] -- [rounded image]
    """
    orig_data = img.get_fdata()
    new_data = np.round(orig_data)
    new_img = nib.Nifti1Image(new_data.astype(int), img.affine)
    return new_img
