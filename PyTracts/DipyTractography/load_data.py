from dipy.io.image import load_nifti, load_nifti_data
from dipy.io.gradients import read_bvals_bvecs
from dipy.core.gradients import gradient_table
from logs import messages
from pathlib import Path


class LoadData:
    def __init__(self, dwi_file: Path, bvec_fname: Path, bval_fname: Path):
        self.dwi_file = dwi_file
        self.bvec_fname = bvec_fname
        self.bval_fname = bval_fname

    def __str__(self):
        str_to_print = messages.LOAD_DATA(
            dwi_file=self.dwi_file, bvec=self.bvec_fname, bval=self.bval_fname
        )
        return str_to_print

    def load_dwi(self):
        data, affine, hardi_img = load_nifti(self.dwi_file, return_img=True)
        return data, affine, hardi_img

    def load_gtab(self):
        bvals, bvecs = read_bvals_bvecs(str(self.bval_fname), str(self.bvec_fname))
        gtab = gradient_table(bvals, bvecs)
        return gtab

    def run(self):
        data, affine, hardi_img = self.load_dwi()
        gtab = self.load_gtab()
        return data, affine, hardi_img, gtab
