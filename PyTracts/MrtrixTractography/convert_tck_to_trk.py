from pathlib import Path
from PyTracts.utils.utillities import check_dir_existence, check_files_existence
import PyTracts.MrtrixTractography.mrtrix_methods as mrt_methods
from logs import messages
import logging


class ConvertTck2Trk:
    def __init__(self, tck_stream: Path, dwi_nii: Path, trk_stream: Path):
        self.tck_stream = tck_stream
        self.dwi = dwi_nii
        self.trk_stream = trk_stream
        self.exist = check_files_existence([self.trk_stream])

    def __str__(self):
        str_to_print = messages.CONVERT_TCK_2_TRK.format(
            subj_dir=self.dwi.parent.parent,
            dwi_name=self.dwi.name,
            tck_name=self.tck_stream.name,
            trk_name=self.trk_stream.name,
        )
        return str_to_print

    def convert(self):
        trk_file = mrt_methods.convert_tck_to_trk(
            self.tck_stream, self.dwi, self.trk_stream
        )

    def run(self):
        if not self.exist:
            logging.info("Converting tractography file from .tck format to .trk")
            self.convert()
        else:
            logging.warning("Given input for .trk tractography file already exists.")
        return self.trk_stream
