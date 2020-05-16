from dipy.tracking import utils
import numpy as np


class CalculateSeeds:
    def __init__(
        self, gray_mask: np.ndarray, affine: np.ndarray, density: list = [2, 2, 2]
    ):
        self.gray_mask = gray_mask
        self.affine = affine
        self.density = density

    def calculate_seeds(self):
        seeds = utils.seeds_from_mask(self.gray_mask, self.affine, density=self.density)
        return seeds
