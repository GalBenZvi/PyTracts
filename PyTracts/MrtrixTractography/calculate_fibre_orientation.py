from pathlib import Path
from PyTracts.utils.utillities import check_dir_existence, check_files_existence
import PyTracts.MrtrixTractography.mrtrix_methods as mrt_methods
from logs import messages
import logging


class CalculateFibreOrientation:
    def __init__(
        self, dwi_file: Path, dwi_mask: Path, response_dict: dict, tracts_dir: Path
    ):
        """
        Estimation of fibre orientation distributions
        Arguments:
            dwi_file {Path} -- [path to preprocessed dwi file]
            response_dict {dict} -- [dictionary with "wm","gm","csf" as keys, and paths to corresponding response_{tissue}.txt files as values]

        Returns:
            [type] -- [description]
        """
        self.dwi_file = dwi_file
        self.tracts_dir = tracts_dir
        self.mask = dwi_mask
        self.response_dict = response_dict
        self.fod_dict = dict()
        for tissue in ["wm", "gm", "csf"]:
            self.fod_dict[tissue] = tracts_dir / f"FOD_{tissue}.mif"
        self.exist = check_files_existence(list(self.fod_dict.values()))

    def __str__(self):
        str_to_print = messages.CALCULATE_FIBER_ORIENTATION.format(
            subj_dir=self.tracts_dir.parent,
            dwi_name=self.dwi_file.name,
            mask_name=self.mask.name,
            tracts_dir=self.tracts_dir,
        )
        return str_to_print

    def fibre_orientation(self):
        mrt_methods.calculate_fibre_orientation(
            self.dwi_file, self.mask, self.response_dict, self.fod_dict
        )

    def run(self):
        if not self.exist:
            logging.info("Estimating Fibre Orientation Distributions")
            self.fibre_orientation()
        else:
            logging.warning(
                "Already estimated fibre orientation distributions, continuing..."
            )
        return self.fod_dict
