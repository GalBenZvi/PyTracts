from pathlib import Path
from PyTracts.utils.utillities import check_dir_existence, check_files_existence
import PyTracts.MrtrixTractography.mrtrix_methods as mrt_methods
from logs import messages
import logging


class GenerateTck:
    def __init__(self, fod_dict: dict, seg_5tt: Path, tract_dir: Path):
        """
        Generate .tck tracts file using Mrtrix3's iFOD2 algorithm
        Arguments:
            fod_dict {dict} -- [dictionary with "wm","gm","csf" as keys, and paths to corresponding FOD_{tissue}.mif files as values]
            tract_dir {Path} -- [Path to directory containing all tracts-processing-related files]
            seg_5tt {Path} -- [Path to 5-tissue-type.mif file]

        Returns:
            [type] -- [description]
        """
        self.fod_wm = fod_dict["wm"]
        self.seg_5tt = seg_5tt
        self.tractogram = tract_dir / "tractogram.tck"
        self.exist = check_files_existence([self.tractogram])

    def __str__(self):
        str_to_print = messages.GENERATE_TCK.format(
            subj_dir=self.tractogram.parent.parent,
            fod_wm_name=self.fod_wm.name,
            tracts_dir=self.tractogram.parent,
            seg_5tt_name=self.seg_5tt.name,
            tractogram_name=self.tractogram.name,
        )
        return str_to_print

    def generate_tracts(self):
        tractogram = mrt_methods.generate_tracts(
            self.fod_wm, self.tractogram, self.seg_5tt
        )
        logging.info(f"Generated tractogram.tck file at {Path(tractogram).parent}")

    def run(self):
        if not self.exist:
            logging.info("Generating tractogram.tck file using iFOD2 algorithm...")
            self.generate_tracts()
        else:
            logging.warning(
                "Already generated tractogram.\n to recreate it, please remove the currently existing .tck file"
            )
        return self.tractogram
