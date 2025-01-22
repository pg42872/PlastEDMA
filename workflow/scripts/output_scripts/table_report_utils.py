import pandas as pd
from hmm_process import *

def write_excel(outdir: str, dataframe: pd.DataFrame, ids_per_model: dict) -> None:
    writer = pd.ExcelWriter(outdir + "report_table.xlsx", engine = "openpyxl")
    dataframe.to_excel(writer, sheet_name = "Table_Report", index = 0)
    df1 = pd.DataFrame.from_dict(ids_per_model, orient = "index")
    df1.to_excel(writer, sheet_name = "Model_Sequences")
    writer.close()

def create_summary_dict(dataframe) -> dict:
    return {
        "models": [model for model in get_models_names(dataframe, to_list = True, only_relevant = True)], 
        "querys": get_match_ids(dataframe, to_list = True, only_relevant = True),
        "bit_scores": get_bit_scores(dataframe, to_list = True, only_relevant = True),
        "e_values": get_e_values(dataframe, to_list = True, only_relevant = True)
        }

def write_tsv(outdir, table_name, dataframe) -> None:
    dataframe.to_csv(outdir + table_name, sep = "\t")

def write_csv(outdir, table_name, dataframe) -> None:
    dataframe.to_csv(outdir + table_name)

def check_output(type: str, outdir: str, table_name: str, dataframe: pd.DataFrame, ids_per_model: dict) -> None:
    if type == "tsv":
        write_tsv(outdir, table_name, dataframe)
    elif type == "csv":
        write_csv(outdir, table_name, dataframe)
    elif type == "excel":
        write_excel(outdir, dataframe, ids_per_model)
    else:
        raise TypeError(f'Specified table format {type} is not available. Read documentation for --output_type')
