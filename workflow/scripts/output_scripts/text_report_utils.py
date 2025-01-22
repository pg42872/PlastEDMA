def write_var_file(variable_to_update: str, value) -> dict:
    """_summary_

    Args:
        variable_to_update (str): key refering to the cariable to be updated
        value (str or int): value to be assigned to the 

    Returns:
        dict: _description_
    """
    variables = {"number_init_hmms": 0, 
            "query_names": "", 
            "number_validated_hmms": 0, 
            "bit_threshold": 0, 
            "eval_threshold": 0, 
            "inputed_seqs": 0, 
            "unique_seqs": 0}
    variables[variable_to_update] = value
    return variables


def write_text_report(config, path: str, args: dict, variables: dict):

    if config["hmm_validation"] == True:
        with open(path + "text_report.txt", "w") as f:
            f.write(f"M-PARTY hits report:\n \
                    \nFrom a total number of {variables.number_init_hmms} HMM profiles initially considered, only {len(variables.query_names)} where considered"
                    f"for the final report.\nUser defined validation to true with {args.negative_db} database, from which resulted"
                    f"in {variables.number_validated_hmms}. After annotation, another filtering process was performed considering the values from bit score and E-value from the HMM search run, \n"
                    f"in which the considered bit score threshold was {variables.bit_threshold} and E-value was {variables.eval_threshold}.\n"
                    f"Also, {len(variables.inputed_seqs)} initial query sequences where inputed form {args.input} file, from which {len(variables.unique_seqs)}\n"
                    f"out of these {len(variables.inputed_seqs)} were considered to have a hit against the HMM database.")
            f.close()
    else:
        with open(path + "text_report.txt", "w") as f:
            f.write(f"M-PARTY hits report:\n \
                    \nFrom a total number of {variables.number_init_hmms} HMM profiles initially considered, only {len(variables.query_names)} where considered"
                    "for the final report.\nFiltering process was performed considering the values from bit score and E-value from the HMM search run, \n"
                    f"in which the considered bit score threshold was {variables.bit_threshold} and E-value was {variables.eval_threshold}.\n"
                    f"Also, {len(variables.inputed_seqs)} initial query sequences where inputed form {args.input} file, from which {len(variables.unique_seqs)}\n"
                    f"out of these {len(variables.inputed_seqs)} were considered to have a hit against the HMM database.")
            f.close()