from pathlib import Path
import shutil

def generate_path(parent_path: str, run_specific_dirname: str, target_path: str, output_type: str, create_dir: bool = False) -> str:
    """_summary_

    Args:
        parent_path (str): parent directories
        run_specific_dirname (str): hmm_database_name argument, basically all files from a run are all assigned inside this dir
        target_path (str): sufix path
        output_type (str): type of the file
        create_dir (bool, optional): if it is to create a path with its parents or just return the full path. Defaults to False.

    Returns:
        str: full path
    """
    if create_dir:
        Path(parent_path + run_specific_dirname + target_path + "." + output_type).mkdir(parents = True, exist_ok = True)


def dir_generator_from_list(list_paths: list):
    """Given a list of paths, create those directories of still not present, with all parent directories

    Args:
        list_paths (list): list of paths to be created
    """
    for path in list_paths:
        Path(path).mkdir(parents = True, exist_ok = True)


def dir_remover(parent_path_list: list, run_specific_dirname: str) -> None:
    """Given a list of parent directories, removes all files from a given run inside those parents

    Args:
        parent_path_list (list): list of parent pathds
        run_specific_dirname (str): hmm_database_name argument, basically all files from a run are all assigned inside this dir
    """
    for path in parent_path_list:
        path_to_remove = Path(path)
        path_to_remove = path_to_remove / run_specific_dirname
        if path_to_remove.exists() and path_to_remove.is_dir():
            shutil.rmtree(path_to_remove)


def check_results_directory(output: str) -> str:
    """Automatically creats the path where output should appear. Checks if folder already exists or not in the 
    execution path
    Args:
        output (str): Name for the output folder
    Returns:
        str: Path for the output folder
    """
    Path(output).mkdir(exist_ok=True, parents=True)