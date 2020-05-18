import PyTracts.DipyTractography.dipy_methods as dipy_method
from logs import messages
from pathlib import Path
import nibabel as nib


class GenerateGrayAndWhite:
    def __init__(
        self, labels_fname: Path, dwi_fname: Path, white_label: int, gray_label: int
    ):
        self.labels_fname = labels_fname
        self.dwi_fname = dwi_fname
        self.white_label = white_label
        self.gray_label = gray_label

    def __str__(self):
        str_to_print = messages.GENERATE_WHITE_MATTER(
            labels_file=self.labels_fname, dwi_file=self.dwi_fname
        )
        return str_to_print

    def resample_labels(self):
        resampled_img = dipy_method.resample_to_dwi(self.labels_fname, self.dwi_fname)
        return resampled_img

    def extract_white_and_gray(self, resampled_img: nib.Nifti1Image):
        data = resampled_img.get_fdata()
        white_matter = data == self.white_label
        gray_matter = data == self.gray_label
        return white_matter, gray_matter

    def run(self):
        resampled_img = self.resample_labels()
        white_matter, gray_matter = self.extract_white_and_gray(resampled_img)
        return white_matter, gray_matter
