import pandas as pd
import re
from command_run import run_command
from mparty_util import save_as_tsv
from seq_download import get_fasta_sequences2


def run_UPIMAPI(query: str, outpath: str, upi_database: str, threads: int) -> str:
    """Function to run UPIMAPI by bash with the given arguments from M-PARTY.

    Args:
        query (str): Database path to be searched on by DIAMOND.
        outpath (str): Output directory name with filename.
        upi_database (str): Sequence list to be searched against the database.
        threads (int): Number of threads.

    Returns:
        str: Path to the final TSV file.
    """
    run_command(f'upimapi.py`-i`{query}`-o`{outpath}`--database`{upi_database}`-t`{threads}', sep = "`")
    return outpath + "/UPIMAPI_results.tsv"


def UPIMAPI_parser(filepath: str):
    UPIMAPI_outfile = pd.read_csv(filepath, sep="\t")
    return UPIMAPI_outfile


def UPIMAPI_iter_per_sim(dataframe: pd.DataFrame, search: bool = False, expansion: bool = False, cutoff: float = None) -> dict or list:
    """Given a pandas DataFrame, return a dictionary with a list of sequences form the iteration of the sequence similarity between queries and database sequences.

    Args:
        dataframe (pd.DataFrame): A pandas dataframe with diamond documented columns names as header.

    Returns:
        dict: A dictionary where the keys are intervals of sequence similarity, and values are lists of UniProtKB queries.
    """
    # selecionar colunas com perc. identity juntamente com os IDs das sequencias
    seq_id = dataframe[["qseqid", "sseqid", "pident"]]
    if expansion:
        target_enzymes = {}
        for perc in range(60, 91, 5):
            for _, seq in seq_id.iterrows():
                if seq["pident"] >= perc:
                    ident = re.findall("\|.*\|", seq["qseqid"])
                    ident = re.sub("\|", "", ident[0])
                    if perc not in target_enzymes.keys():
                        target_enzymes[perc] = [ident]
                    else:
                        target_enzymes[perc].append(ident)
        return target_enzymes
    else:
        target_enzymes = []
        if cutoff:
            for _, seq in seq_id.iterrows():
                if seq["pident"] >= cutoff:
                    if search:
                        ident = re.findall("\|.*\|", seq["qseqid"])
                        ident = re.sub("\|", "", ident[0])
                    else:
                        ident = seq["sseqid"]
                    target_enzymes.append(ident)
        else:
            for _, seq in seq_id.iterrows():
                if search:
                    ident = re.findall("\|.*\|", seq["qseqid"])
                    ident = re.sub("\|", "", ident[0])
                else:
                    ident = seq["sseqid"]
                target_enzymes.append(ident)
        return [*set(target_enzymes)]       


def sigasiga(listinha, outpath):
    df = pd.DataFrame(listinha, columns=["IDs"])
    df.to_csv(outpath, index=False, sep = "\t")
