from PyTracts.DipyTractography.dipy_methods import resample_to_dwi, FSLOUTTYPE
from nipype.interfaces import fsl
from logs import messages
from pathlib import Path
import nibabel as nib
import json
import os


class GenerateGrayAndWhite:
    def __init__(
        self,
        wm_fname: Path,
        gm_fname: Path,
        dwi_fname: Path,
    ):
        self.wm_fname = wm_fname
        self.gm_fname = gm_fname
        self.dwi_fname = dwi_fname

    def __str__(self):
        str_to_print = messages.GENERATE_WHITE_MATTER(
            labels_file=self.labels_fname, dwi_file=self.dwi_fname
        )
        return str_to_print

    def register_masks_to_epi(self):
        registered = []
        for mask, out_file in zip(
            [self.wm_fname, self.gm_fname],
            [f"white_matter{FSLOUTTYPE}", f"gray_matter{FSLOUTTYPE}"],
        ):
            out_file = self.dwi_fname.parent / out_file
            registered.append(out_file)
            flt = fsl.FLIRT()
            flt.inputs.in_file = mask
            flt.inputs.reference = self.dwi_fname
            flt.inputs.bins = 256
            flt.inputs.searchr_x = [-90, 90]
            flt.inputs.searchr_y = [-90, 90]
            flt.inputs.searchr_z = [-90, 90]
            flt.inputs.dof = 12
            flt.inputs.interp = "trilinear"
            flt.inputs.out_file = out_file
            if not out_file.exists():
                print(flt.cmdline)
                flt.run()
        wm_registered, gm_registered = registered
        return wm_registered, gm_registered

    def extract_white_and_gray(self, wm_fname: Path, gm_fname: Path):
        white_matter = nib.load(wm_fname).get_fdata()
        white_matter[white_matter != 0] = 1
        gray_matter = nib.load(gm_fname).get_fdata()
        gray_matter[gray_matter != 0] = 1
        return white_matter, gray_matter

    def run(self):
        wm_fname, gm_fname = self.register_masks_to_epi()
        white_matter, gray_matter = self.extract_white_and_gray(wm_fname, gm_fname)
        return white_matter, gray_matter
