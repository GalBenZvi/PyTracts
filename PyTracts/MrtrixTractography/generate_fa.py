from pathlib import Path
from PyTracts.utils.utillities import check_dir_existence, check_files_existence
import PyTracts.MrtrixTractography.mrtrix_methods as mrt_methods
from logs import messages
import logging


class GenerateFA:
    def __init__(self, dwi_file: Path, mask_file: Path, tracts_dir: Path):
        """
        Class to generate FA and DTI images from preprocessed dwi image and mask.
        Arguments:
            tracts_dir {Path} -- [output directory, where fa.mif and dti.if will be produced]
            dwi_file {Path} -- [path to preprocessed dwi image]
            mask_file {Path} -- [path to dwi mask image]
        """
        self.tract_dir = tracts_dir
        check_dir_existence(self.tract_dir)
        self.dwi_file = dwi_file
        self.mask_file = mask_file
        self.dti_file = tracts_dir / "dti.mif"
        self.fa_file = tracts_dir / "fa.mif"
        self.exist = check_files_existence([self.fa_file, self.dti_file])

    def __str__(self):
        str_to_print = messages.GENERATE_FA.format(
            subj_dir=self.tract_dir.parent,
            dwi_name=self.dwi_file.name,
            mask_name=self.mask_file.name,
            tracts_dir=self.tract_dir,
            FA_name=self.fa_file.name,
        )
        return str_to_print

    def generate_fa(self):
        fa_file = mrt_methods.fit_tensors(
            self.dwi_file, self.mask_file, self.dti_file, self.fa_file
        )
        return fa_file

    def run(self):
        if not self.exist:
            logging.info("Generating FA image for group-level analysis")
            self.fa_file = self.generate_fa()
        else:
            logging.warning(
                "Already generated FA image for group-level analysis, continuing..."
            )
        return self.fa_file
