from command_run import run_command
import pandas as pd


def run_KMA(input_db: str, output_db: str, meta_input: str, meta_out: str, threads: int, paired_end: bool = False, second_input: str = None) -> str:
    """Runs KMA with in the command line and returns the output file

    Args:
        input_DB (str): Reference nucleotide database
        output_DB (str): Name for the indexed referente database
        meta_input (str): Metagenome sample
        meta_out (str): Output name
        threads (int): Number of threads
        paired_end (bool, optional): Set to true of a second input (paired end) is given. Defaults to False.
        second_input (str, optional): Paired end input file. Defaults to None.

    Raises:
        ValueError: Raises a Value error if paired_end is True but no file is given

    Returns:
        str: Path of the output file
    """
    run_command(f'kma`index`-i`{input_db}`-o`{output_db}', sep = "`")
    run_command(f'kma`-i`{meta_input}`-o`{meta_out}`-t_db`{output_db}`-t`{threads}`-1t1`-mem_mode`-ef', sep = "`")
    if paired_end and second_input == None:
        raise ValueError("When paired end is flaged, a second input file is mandatory")
    elif paired_end and second_input != None:
        run_command(f'kma`-i`{meta_input}`{second_input}`{meta_out}`-t_db`{output_db}`-t`{threads}`-1t1`-mem_mode`-ef', sep = "`")
    return meta_out

def kma_parser(input_file: str) -> pd.DataFrame:
    """Process the output of the KMA run for the desired information

    Args:
        input_file (str): Path for the KMA output file

    Returns:
        pd.Dataframe: A pandas Dataframe with relevant information
    """
    df = pd.read_csv(input_file, sep = "\t", decimal = ",")
    df1 = df.loc[(pd.to_numeric(df["Template_Identity"]) >= 80) & 
                 (pd.to_numeric(df["Template_Coverage"]) >= 60)]
    return df1[["#Template", "Template_Identity", "Template_Coverage", "q_value", "p_value"]]

def get_hit_sequences(dataframe: pd.DataFrame, to_list: bool = False) -> list:
    """Given a Dataframe with info from KMA run, return a list of the hit IDs

    Args:
        dataframe (pd.DataFrame): Dataframe from kma_parser
        to_list (bool, optional): Flag True to return a list as output. Defaults to False.

    Returns:
        list: list of hit IDs from KMA
    """
    if to_list:
        return dataframe["#Template"].tolist()
    else:
        return dataframe["#Template"]
