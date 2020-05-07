from pathlib import Path

FSLOUTTYPE = ".nii.gz"


def check_files_existence(files: list):
    """
    Checks the existence of a list of files (all should exist)
    Arguments:
        file {Path} -- [description]

    Returns:
        [type] -- [description]
    """
    exist = all(Path(f).exists() for f in files)
    return exist


def check_dir_existence(directory: Path):
    if not directory.exists():
        print(f"Creating directory: {directory}")
        directory.mkdir()
