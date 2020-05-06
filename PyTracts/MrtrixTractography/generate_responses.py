from pathlib import Path
from PyTracts.utils.utillities import check_dir_existence, check_files_existence
import PyTracts.MrtrixTractography.mrtrix_methods as mrt_methods
from logs import messages
import logging


class GenerateResponses:
    def __init__(self, dwi_file: Path, mask: Path, tracts_dir: Path):
        """
        Estimating tissue response functions for spherical deconvolution
        Arguments:
            dwi_file {Path} -- [path to preprocessed dwi image]
            mask {Path} -- [path to dwi mask image]
            tracts_dir {Path} -- [path to output directory, where response_*tissue*.txt file will be produced]

        Returns:
            [type] -- [description]
        """
        self.dwi_file = dwi_file
        self.mask = mask
        self.tracts_dir = tracts_dir
        self.response_wm, self.response_gm, self.response_csf = [
            tracts_dir / f
            for f in ["response_wm.txt", "response_gm.txt", "response_csf.txt"]
        ]
        self.exist = check_files_existence(
            [self.response_wm, self.response_gm, self.response_csf]
        )

    def __str__(self):
        str_to_print = messages.GENERATE_RESPONSES.format(
            subj_dir=self.tracts_dir.parent,
            dwi_name=self.dwi_file.name,
            mask_name=self.mask.name,
            tracts_dir=self.tracts_dir,
        )
        return str_to_print

    def generate_response(self):
        response_wm, response_gm, response_csf = mrt_methods.generate_response(
            self.dwi_file, self.mask, self.tracts_dir
        )
        return response_wm, response_gm, response_csf

    def run(self):
        if not self.exist:
            logging.info(
                "Estimating tissue response functions for spherical deconvolution"
            )
            (
                self.response_wm,
                self.response_gm,
                self.response_csf,
            ) = self.generate_response()
        else:
            logging.warning(
                "Already estimated tissue response functions, continuing..."
            )
        response_dict = dict()
        for key, val in zip(
            ["wm", "gm", "csf"], [self.response_wm, self.response_gm, self.response_csf]
        ):
            response_dict[key] = val
        return response_dict
