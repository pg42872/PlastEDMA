import yaml
import json
from workflow.pathing_utils.fixed_paths import PathManager

def get_arguments(args: dict, sequences: list, output_dir: str) -> dict:
    """Converts the arguments given by the CLI to a dictionary"""
    arguments = {"seqids": sequences,
                "database": args.database,
                "input_file": None if sequences == [] else args.input.split("/")[-1],
                "input_file_db_const": args.input_seqs_db_const,
                "consensus": args.consensus,
                "KEGG_ID": args.kegg,
                "InterPro_ID": args.interpro,
                "hmm_database_name": args.hmm_db_name,
                "alignment_method": args.align_method.lower(),
                "msa_aligner": args.aligner,
                "input_type": None if args.input == None else args.input_type,
                "metagenomic": True if args.input_type == "metagenome" else False,
                "hmm_validation": args.hmm_validation,
                "expansion": args.expansion,
                "concat_models": args.concat_hmm_models,
                "output_directory": output_dir,
                "out_table_format": args.output_type,
                "hmmsearch_out_type": args.hmms_output_type,
                "threads": args.threads,
                "workflow": args.workflow,
                "thresholds": [*range(60, 91, 5)] if args.expansion else False,
                "verbose": args.verbose,
                "overwrite": args.overwrite
                }
    return arguments


def check_input_arguments(args: dict, verbose: bool, kma_res: bool) -> bool:
    if args.hmm_validation == True and args.workflow == "annotation" and args.input == None:
        if verbose:
            print("No input file detected. Proceding to validation")
        return False
    
    elif args.workflow == "database_construction" and args.input == None and args.kegg == None and args.interpro == None and args.input_seqs_db_const == None:
        if verbose:
            print("No input file detected. Proceding to model construction")
        return False
    
    elif args.input_type == "metagenome" and kma_res == False:
        return False
    
    else: return True

def check_config(args: dict):
    """Function that runs at the begining of every pipeline to check if the necessary arguments were given
    Args:
        args (dict): Arguments passed through the CLI
    """
    if args.config_file != None and (args.input != None or args.kegg != None or args.interpro != None):
        raise ValueError("config file cannot be given with other arguments")
    elif args.hmm_db_name is None:
        raise TypeError("Missing hmm database name! Make sure --hmm_db_name option is filled")
    elif args.workflow == "database_construction" and args.input_seqs_db_const is None and args.kegg is None and args.interpro is None:
        raise TypeError("Missing input sequences to build HMM database")
    if args.interpro is not None:
        # for interpro is only possible to run for aminoacids and so for HMM and not KMA and raw metagenomes
        if args.input_type == "metagenome":
            raise ValueError("Metagenomic samples cannot be annalyzed with proteins as database")
        elif args.interpro[0].startswith("IPR") and len(args.interpro) > 1:
            raise ValueError("Give only 1 InterPro ID (IPR******)")
        elif not args.interpro[0].startswith("IPR") and not args.interpro[0].startswith("A"):
            raise ValueError("Must input and IPR ID or protein ID from InterPro starting with 'A'")
    
def write_yaml_json(config_type: str, out_dir: str, args_dict: dict, to_output: bool):
    if to_output:
        config_filename = "parameters"
    else:
        out_dir = PathManager.config_path
        config_filename = "config"
    if config_type == "yaml":
        with open(f'{out_dir}/{config_filename}.yaml', "w") as file:
            yaml.dump(args_dict, file)
            file.close()
    else:
        with open(f'{out_dir}/{config_filename}.json', "w") as file:
            document = json.dumps(args_dict)
            file.write(document)
            file.close()