#!/usr/bin/env python
# run tool main script without indicating python

"""
M-PARTY - Mining Protein dAtasets for Target EnzYmes

by José Freitas

Set 2024
"""

import sys
import shutil
# sys.path.insert(0, f'{"/".join(sys.path[0].split("/")[:-1])}/share')
sys.path.append(f'{sys.path[0]}/workflow/scripts')
sys.path.append(f'{sys.path[0]}/workflow/pathing_utils')
# sys.path.append(f'{sys.path[0]}/M-PARTY')
import os
from pathlib import Path
import time
import yaml
import json
import re
import pandas as pd
from tqdm import tqdm 
import snakemake
import itertools
import threading

from workflow.pathing_utils.cli_args import get_parser, process_arguments
from hmmsearch_run import run_hmmsearch
from hmm_process import *
from hmm_vali import concat_final_model, file_generator, exec_testing, hmm_filtration, remove_fp_models, make_paths_dic, delete_inter_files
import UPIMAPI_parser
from seq_download import get_fasta_sequences
from CDHIT_seq_download import fasta_retriever_from_cdhit
import CDHIT_parser
from mparty_util import build_upi_query_db, threshold2clusters, get_tsv_files, save_as_tsv, concat_code_hmm, compress_fasta, return_fasta_content, check_id, ask_for_overwrite
import BLAST_parser
import DIAMOND_parser
from command_run import run_tcoffee, run_hmmbuild, run_hmmemit, concat_fasta
from InterPro_retriever import get_IP_sequences
from KEGG_retriever import get_kegg_genes
from KMA_parser import run_KMA, kma_parser, get_hit_sequences
from config.process_arguments import get_arguments, check_input_arguments, check_config, write_yaml_json
import output_scripts.table_report_utils as table_report_utils
import output_scripts.text_report_utils as text_report_utils
from workflow.pathing_utils.fixed_paths import PathManager, declare_fixed_paths
from workflow.pathing_utils.path_generator import dir_generator_from_list, generate_path, dir_remover, check_results_directory


def read_config(filename: str) -> tuple:
    config_type = filename.split(".")[-1]
    if config_type == "yaml":
        with open(filename) as stream:
            try:
                config_file = yaml.safe_load(stream)
                stream.close()
            except yaml.YAMLError as exc:
                print(exc)
    elif config_type == "json":
        with open(filename) as stream:
            try:
                config_file == json.load(stream)
                stream.close()
            except json.decoder.JSONDecodeError as exc:
                print(exc)
    else:
        quit("Config file must be in .yaml or .json format! Get an example config file from config/ folder.")
    return config_file, config_type


def clean_sequence_ids(line: str, remove_excess_id: bool, ip: bool, kegg: bool, kma_res: bool) -> str:
    """Function that receives a string and cleans it based on predefined patterns

    Args:
        line (str): line string
        remove_excess_id (bool): Decide wether to remove the excess part of UniProt IDs
        ip (bool): Set to True if sequences from filename were retrieved from InterPro, which has a specific nomenclature
        for the FASTA entries
        kegg (bool): Set to True if sequences from filenames were retrieved from KEGG, which has a specific nomenclature
        for the FASTA entries.
        kma_res (bool): Set to True if this function is set to run for the processing of KMA results

    Returns:
        str: Cleaned string
    """
    if kegg:
        return re.search(r">(\S+)", line).group(1)
    elif ip:
        return re.search(r">([^|]+)\|", line).group(1)
    elif kma_res:
        return line.replace(">", "").strip()
    else: 
        if not remove_excess_id:
            return line.split(" ")[0][1:]
        else:
            try:
                return re.search(r"\|(.*)\|", line).group(1)
            except Exception:
                identi = line.split(" ")[0]
                return identi.replace(">", "")


def parse_fasta(filename: str, remove_excess_id: bool = True, ip: bool = False, kegg: bool = False, kma_res: bool = False, verbose: bool = False) -> list:
    """Given a FASTA file, returns the IDs from all sequences in that file.
    If file not present, program will be quited and TypeError message raised.

    Args:
        filename (str): Name of FASTA file.
        remove_excess_id (bool, optional): Decide wether to remove the excess part of UniProt IDs. Defaults to True.
        ip (bool, optional): Set to True if sequences from filename were retrieved from InterPro, which has a specific nomenclature
        for the FASTA entries.
        kegg (bool, optional): Set to True if sequences from filenames were retrieved from KEGG, which has a specific nomenclature
        for the FASTA entries.
        kma_res (bool, optional): Set to True if this function is set to run for the processing of KMA results. Defaults to False.
        verbose (bool, optional): Set to True to print aditional messages of wath is happening. Defaults to False.

    Returns:
        list: A list containing IDs from all sequences
    """
    uniq_ids = []
    if check_input_arguments(args, verbose=verbose,kma_res=kma_res) == False:
        return uniq_ids
    else:
        try:
            with open(filename, "r") as handlefile:
                try:
                    for line in handlefile:
                        if line.startswith(">"):
                            uniq_ids.append(clean_sequence_ids(line, remove_excess_id, ip, kegg, kma_res))
                    if verbose:
                        print(f'Input file {filename} detected and sequence IDs retrieved\n')
                        time.sleep(1)
                except Exception as exc:
                    print(exc)
                    quit("File must be in FASTA format.")
        except TypeError:
            raise TypeError("Missing input file! Make sure -i option is filled")
        return uniq_ids


def write_config(input_file: str, out_dir: str, to_output: bool  = False):
    """Given a input file, output directory, and a name to assign to the new config file, write that same config file
    accordingly to the given arguments

    Args:
        input_file (str): Name for the input FASTA file
        out_dir (str): Name for the output directory where result shall be directed
        config_filename (str): Name to be given to the new config file
    """
    if args.workflow == "database_construction" and args.input == None:
        seq_ids = []
    if args.input != None:
        file_stats = os.stat(input_file)
        if file_stats.st_size / (1024 * 1024) > 400:
            seq_ids = "too_big"
        else:
            seq_ids = parse_fasta(input_file)
    if args.hmm_validation and args.workflow != "database_construction" and args.workflow != "both" and args.input == None:
        out_dir = None
    else:
        check_results_directory(out_dir)
        arguments = get_arguments(args, seq_ids, out_dir)

    try:
        Path(PathManager.config_path).mkdir(parents = True, exist_ok = True)
        write_yaml_json(config_type="yaml", out_dir=out_dir, args_dict=arguments, to_output=to_output)
    except Exception as exc:
        raise Exception(exc)


def file_generator(path: str, full_path: bool = False):
    """Function that yield the name of all and only files inside a directory in the given path, for iteration purposes
    Args:
        path (str): Path for the folder to be analyzed

    Yields:
        str: Yield the name of each file inside the given directory
    """
    for file in os.listdir(path):
        if os.path.isfile(os.path.join(path, file)):
            if full_path:
                yield os.path.join(path, file)
            else:
                yield file


def table_report(dataframe: pd.DataFrame, path: str, type_format: str, db_name: str):
    """Saves a table in a user specified format, with the processed and filtered information from the 
    hmmsearch execution with the HMMs against the query sequences.

    Args:
        dataframe (pd.DataFrame): Dataframe with only the relevant information from hmmsearch execution 
        for all hmm from all threshold ranges.
        path (str): output path.
        type_format (str): Specify the output format.
        db_name(str): Name of the databases name.

    Raises:
        TypeError: Raises TypeError error if user gives an unsupported output format.
    """
    summary_dic = table_report_utils.create_summary_dict(dataframe=dataframe)
    
    if args.expansion:
        indexes = dataframe.index.values.tolist()
        for i in range(len(indexes)):
            summary_dic["models"][i] = f'{indexes[i]}_{summary_dic["models"][i]}'
    
    df = pd.DataFrame.from_dict(summary_dic)
    list_ids_permodel = {}
    if not args.expansion:
            mother_seqs = f'{sys.path[0]}/resources/Data/FASTA/{db_name}/CDHIT/clusters/'
            for model in tqdm(list(set(summary_dic["models"])), desc = "Tracebacking model's sequences", unit = "model"):
                for file in file_generator(mother_seqs):
                    if file.split(".")[0] == model:
                        if model not in list_ids_permodel:
                            list_ids_permodel[model] = parse_fasta(os.path.join(mother_seqs, file))
                            break
    else:
        mother_seqs = f'{sys.path[0]}/resources/Data/FASTA/{db_name}/CDHIT/'
        for val in summary_dic["models"]:
            thresh = val.split("_")[0]
            model = val.split("_")[-1]
            for folder in os.listdir(mother_seqs):
                if os.path.isdir(os.path.join(mother_seqs, folder)) and folder == thresh:
                    for file in os.listdir(os.path.join(mother_seqs, folder)):
                        if file.endswith(".fasta") and model == file.split(".")[0]:
                            key = thresh + "_" + model
                            if key not in list_ids_permodel:
                                list_ids_permodel[key] = parse_fasta(mother_seqs + folder + "/" + file)
                            # else:
                            #     list_ids_permodel[key].append(parse_fasta(mother_seqs + folder + "/" + file, meta_gen = True if args.input_type == "metagenome" else False))
    
    table_name = "report_table." + type_format
    table_report_utils.check_output(type=type_format, outdir=path, table_name=table_name, dataframe=df, ids_per_model=list_ids_permodel)


def text_report(dataframe: pd.DataFrame, path: str, bit_threshold: float, eval_threshold: float, vali: bool = False, kma: bool = False):
    """Write the final report as .txt file, with a summary of the results from the annotation 
    performed with hmmsearch. Starts by calculating the number of in-built HMM profiles, and gives an insight of the 
    filtration thresholds.

    Args:
        dataframe (pd.DataFrame): Dataframe with only the relevant information from hmmsearch execution 
        for all hmm from all threshold ranges.
        path (str): output path.
    """
    # number of initial HMM profiles
    number_init_hmms, number_validated_hmms = 0, 0
    for dir in os.listdir(PathManager.hmm_database_path):
        if os.path.isdir(os.path.join(PathManager.hmm_database_path, dir)):
            for _ in os.listdir(os.path.join(PathManager.hmm_database_path, dir)):
                number_init_hmms += 1

    if vali:
        for dir in os.listdir(PathManager.validated_hmm_dir):
            if os.path.isdir(os.path.join(PathManager.validated_hmm_dir, dir)):
                for _ in os.listdir(os.path.join(PathManager.validated_hmm_dir, dir)):
                    number_validated_hmms += 1

    # get the IDs from all hits after quality check
    query_names = get_match_ids(dataframe, to_list = True, only_relevant = True)

    # get number of hits given for each sequence
    number_hits_perseq = get_number_hits_perseq(query_names)

    # get the unique sequences
    unique_seqs = get_unique_hits(query_names)
    inputed_seqs = config["seqids"]
    
    variables = text_report_utils.write_var_file()
    text_report_utils.write_text_report(config, path, args, variables)


def get_number_hits_perseq(hit_ids_list: list) -> dict:
    """Given a list of sequences IDs from the hits against the hmm models from hmmsearch, counts the number of each ID.

    Args:
        hit_ids_list (list): List of sequence IDs.

    Returns:
        dict: Dictionary containing each ID as key and the respective number of occurrences as value.
    """
    counter = {}
    for i in hit_ids_list:
        counter[i] = counter.get(i, 0) + 1
    return counter


def get_unique_hits(hit_ids_list: list) -> list:
    """Given a list of sequence IDs from the hits against the hmm models from hmmsearch, return a new list with only the unique elements.

    Args:
        hit_IDs_list (list): List of sequence IDs.

    Returns:
        list: List with only a single occurrence of each ID.
    """
    unique_ids_list = []
    for x in hit_ids_list:
        if x not in unique_ids_list:
            unique_ids_list.append(x)
    return unique_ids_list


def get_aligned_seqs(config, hit_ids_list: list, path: str, inputed_seqs: str, kma_alignfile: str = None):
    """Writes an ouput Fasta file with the sequences from the input files that had a hit in hmmsearch 
    annotation against the hmm models.

    Args:
        hit_ids_list (list): list of IDs that hit.
        path (str): ouput path.
        inputed_seqs (str): name of the initial input file.
    """
    # returns a list the sequences that hit against the models (only one entry)
    unique_ids = get_unique_hits(hit_ids_list)

    if config["seqids"] == "too_big":
        check_id(inputed_seqs, path, unique_ids)
    
    else:
        with open(path + "aligned.fasta", "w") as wf:
            if args.input_type == "metagenome":
                input_ids = parse_fasta(kma_alignfile, remove_excess_id = False, kma_res = True)
            else:
                input_ids = parse_fasta(inputed_seqs, remove_excess_id = False)
            
            if args.input_type == "metagenome":
                inp_seqs = kma_alignfile
            else: 
                inp_seqs = inputed_seqs
            with open(inp_seqs, "r") as rf:
                lines = rf.readlines()
                for x in unique_ids:
                    if x in input_ids:
                        iterador = iter(lines)
                        linha = next(iterador)
                        while linha is not None:
                            if x not in linha:
                                linha = next(iterador, None)
                                continue
                            elif x in linha:
                                wf.write(linha)
                                linha = next(iterador, None)
                                while linha is not None and not linha.startswith(">"):
                                    wf.write(linha)
                                    linha = next(iterador, None)
                            elif x not in linha and linha.startswith(">"):
                                break
                            linha = next(iterador, None)
                    else:
                        continue
            rf.close()
        wf.close()


def generate_output_files(dataframe: pd.DataFrame, 
                          hit_ids_list: list, 
                          inputed_seqs: str,
                          config_file,
                          bit_threshold: float = None, 
                          eval_threshold: float = None, 
                          kma: bool = False,
                          kma_alignfile: str = None):
    """Function that initializes the output files creation simultaneously, for now, only two files are generated:
    report and aligned sequences.
    Path will always be the output folder defined by the user when running tool in CLI, so no pat argument is required.

    Args:
        dataframe (pd.DataFrame): Dataframe with only the relevant information from hmmsearch execution.
        hit_ids_list (list): list of Uniprot IDs that hit.
        inputed_seqs (str): name of the initial input file.
    """
    out_folder = args.output + "/"
    if kma:
        get_aligned_seqs(config_file, hit_ids_list, out_folder, inputed_seqs, kma_alignfile = kma_alignfile)
        dataframe.to_excel(f'{out_folder}report_table.xlsx', sheet_name = "Table_Report", index = 0)
    else:
        table_report(dataframe, out_folder, args.output_type, args.hmm_db_name)
        if args.report_text:
            if args.hmm_validation:
                text_report(dataframe, out_folder, bit_threshold, eval_threshold, vali = True)
            else:
                text_report(dataframe, out_folder, bit_threshold, eval_threshold)
        get_aligned_seqs(config_file, hit_ids_list, out_folder, inputed_seqs)
    if args.display_config:
        write_config(args.input, args.output, to_output=True)


def database_construction(config):
    """Pipeline for the database construction workflow

    Args:
        config (file): The parsed config file object

    Raises:
        ValueError: If input type is metagenome at the same time as the interpro flag is given with an ID
    """
    print("HMM database construction workflow from user input started...\n")
    time.sleep(1)
    if os.path.exists(os.path.join(sys.path[0], f'resources/Data/FASTA/{args.hmm_db_name}/')):
        if args.overwrite:
            try:
                dir_remover(['resources/Data/FASTA', 'resources/Alignments', 'resources/Data/HMMs'], args.hmm_db_name)
                if args.verbose:
                    print(f"Deleted previously created files from {args.hmm_db_name}\n")
            except Exception as exc:
                print(exc)
        else:
            overwrite = ask_for_overwrite(args.hmm_db_name, verbose=args.verbose)
            if overwrite:
                try:
                    dir_remover(['resources/Data/FASTA', 'resources/Alignments', 'resources/Data/HMMs'], args.hmm_db_name)
                    if args.verbose:
                        print(f"Deleted previously created files from {args.hmm_db_name}\n")
                except Exception as exc:
                    raise(exc)
            else:
                print("Database for that name is already present. If you wish to create a new database,\neither overwrite the existant or give a different HMM database name.")
                quit("M-PARTY has finished its execution")

    time.sleep(2)

    if args.expansion:
        expand_base_sequences(config=config)

    else:
        # make necessary directories
        dir_generator_from_list([PathManager.tcoffee_path, PathManager.cdhit_path / "clusters", PathManager.hmm_database_path])
        if args.kegg:
            # if given ID is Kegg Orthology
            if args.kegg[0].startswith("K"):
                if args.input_type_db_const == "nucleic":
                    kegg_sequences = get_kegg_genes(PathManager.fasta_type_dir / Path(args.kegg[0]).with_suffix(".fasta"), 
                                                    type_seq = "nuc",
                                                    ko = args.kegg, 
                                                    verbose = args.verbose)
                else:
                    kegg_sequences = get_kegg_genes(PathManager.fasta_type_dir / Path(args.kegg[0]).with_suffix(".fasta"), 
                                                    ko = args.kegg, 
                                                    verbose = args.verbose)

            # If given ID is an E.C. number
            else:
                if args.input_type_db_const == "nucleic":
                    kegg_sequences = get_kegg_genes(PathManager.fasta_type_dir / Path(args.kegg[0]).with_suffix(".fasta"), 
                                                    type_seq = "nuc", 
                                                    ec_number = args.kegg, 
                                                    verbose = args.verbose)
                else:
                    kegg_sequences = get_kegg_genes(PathManager.fasta_type_dir / Path(args.kegg[0]).with_suffix(".fasta"), 
                                                    ec_number = args.kegg, 
                                                    verbose = args.verbose)

            # Only build HMMs if input is protein or nucleic
            if args.input_type != "metagenome":
                build_hmms_from_seqs(kegg_sequences)

        if args.interpro:
            # for interpro is only possible to run for aminoacids and so for HMM and not KMA and raw metagenomes
            if args.input_type == "metagenome":
                raise ValueError("Metagenomic samples cannot be annalyzed with proteins as database")

            # if given ID is a InterProt ID
            elif args.interpro[0].startswith("IPR") and len(args.interpro) == 1:
                filename = args.interpro[0]
                inp_seqs = get_IP_sequences(PathManager.fasta_type_dir / Path(filename).with_suffix(".fasta"), 
                                            interpro_ID = args.interpro, 
                                            reviewed = args.curated, 
                                            verbose = args.verbose)

            # if given ID is a list of proteins from InterProt
            elif args.interpro[0].startswith("A"):
                filename = args.interpro[0]
                inp_seqs = get_IP_sequences(PathManager.fasta_type_dir / Path(filename).with_suffix(".fasta"), 
                                            protein = args.interpro, 
                                            verbose = args.verbose)

            # Start HMM construction
            build_hmms_from_seqs(inp_seqs, "InP", ident_perc=0.8)

        # if a FASTA file with interest proteins/nucleiotides is given
        if args.input_seqs_db_const:
            # Will not build HMMs if input is a metagenome
            if args.input_type == "metagenome":
                shutil.copyfile(args.input_seqs_db_const, PathManager.fasta_type_dir)

            else:
                # Start HMM construction
                shutil.copyfile(args.input_seqs_db_const, PathManager.fasta_type_dir / args.input_seqs_db_const.split("/")[-1].split(".")[0])
                build_hmms_from_seqs(sequences=args.input_seqs_db_const,
                                type_seq=args.input_type_db_const, 
                                from_database=args.input_seqs_db_const.split("/")[-1].split(".")[0])

        # remove files wrongly going to the root dir
        files = [f for f in os.listdir('.') if os.path.isfile(f)]
        for file in files:
            if file.endswith(".dnd"):
                delete_inter_files(file)

    if args.hmm_validation:
        validate_hmm(config=config)


def expand_base_sequences(config):
    Path("resources/Data/FASTA/DataBases").mkdir(parents = True, exist_ok = True)
    Path(f'resources/Data/Tables/{args.hmm_db_name}').mkdir(parents = True, exist_ok = True)
    query_db = build_upi_query_db("resources/Data/FASTA/DataBases", config = config, verbose = config["verbose"])

    if config["alignment_method"] == "diamond":
        ### FASTA to DMND
        diamond_file = build_diamond_DB(query_db, "resources/Data/FASTA/", verbose = config["verbose"])  # ver a cena do overwrite para estes passos
        Path(f'resources/Alignments/{args.hmm_db_name}/BLAST/diamond_output/').mkdir(parents = True, exist_ok = True)
        aligned_tsv = run_DIAMOND(args.input_seqs_db_const, f'resources/Alignments/{args.hmm_db_name}/{config["alignment_method"].upper()}/diamond_output/out.tsv', diamond_file, args.threads)
        handle = DIAMOND_parser(aligned_tsv)
        dic_enzymes = DIAMOND_iter_per_sim(handle)
        if config["verbose"]:
            print(f'Saving IDs from the ranges of {config["thresholds"]} percentages of similarity.\n')
        save_as_tsv(dic_enzymes, f'resources/Data/Tables/{args.hmm_db_name}/DIAMOND_results_per_sim.tsv')

    elif config["alignment_method"] == "upimapi":
        # aligned_TSV = run_UPIMAPI(query_DB, f'resources/Alignments/{args.hmm_db_name}/{config["alignment_method"].upper()}/upimapi_results', args.input_seqs_db_const, args.threads)
        aligned_tsv = f'resources/Alignments/{args.hmm_db_name}/{config["alignment_method"].upper()}/upimapi_results/UPIMAPI_results.tsv'
        handle = UPIMAPI_parser(aligned_tsv)
        dic_enzymes = UPIMAPI_iter_per_sim(handle)
        if config["verbose"]:
            print(f'Saving IDs for the minimum cutoff values of {config["thresholds"]} percentages of similarity.\n')
        save_as_tsv(dic_enzymes, f'resources/Data/Tables/{args.hmm_db_name}/UPIMAPI_results_per_sim.tsv')

    elif config["alignment_method"] == "blast":
        # blastdb_file = build_blast_DB(query_DB, "resources/Data/FASTA/DataBases/BLAST", args.input_type_db_const, verbose = config["verbose"])
        Path(f'resources/Alignments/{args.hmm_db_name}/BLAST/BLAST_results').mkdir(parents = True, exist_ok = True)
        # run_BLAST(args.input_seqs_db_const, f'resources/Alignments/{args.hmm_db_name}/BLAST/BLAST_results/test.tsv', blastdb_file, 8)
        aligned_tsv = f'resources/Alignments/{args.hmm_db_name}/BLAST/BLAST_results/test.tsv'
        handle = BLAST_parser(aligned_tsv)
        dic_enzymes = BLAST_iter_per_sim(handle)
        if config["verbose"]:
            print(f'Saving IDs from the ranges of {config["thresholds"]} percentages of similarity.\n')
        Path(f'resources/Data/Tables/{args.hmm_db_name}/').mkdir(parents = True, exist_ok = True)
        save_as_tsv(dic_enzymes, f'resources/Data/Tables/{args.hmm_db_name}/BLAST_results_per_sim.tsv')

    else:
        raise ValueError("--align_method flag only ranges from 'diamond', 'upimapi' or 'blast'. Chose one from the list.")

    Path(f'resources/Data/FASTA/{args.hmm_db_name}/{config["alignment_method"].upper()}/').mkdir(parents = True, exist_ok = True)
    Path(f'resources/Data/Tables/{args.hmm_db_name}/CDHIT_clusters/').mkdir(parents = True, exist_ok = True)
    for thresh in config["thresholds"]:
        if config["verbose"]:
            print(f'Retrieving sequences from {thresh} range\n')
        try:
            get_fasta_sequences(f'resources/Data/Tables/{args.hmm_db_name}/{config["alignment_method"].upper()}_results_per_sim.tsv', f'resources/Data/FASTA/{args.hmm_db_name}/{config["alignment_method"].upper()}/{thresh}.fasta')
        except Exception as exc:
            print(exc)
            raise FileNotFoundError(f'resources/Data/Tables/{config["alignment_method"].upper()} not found.')
        ### run CDHIT
        if config["verbose"]:
            print(f'CDHIT run for {thresh} range\n')
            Path(f'resources/Data/FASTA/{args.hmm_db_name}/CDHIT/{thresh}/').mkdir(parents = True, exist_ok = True)
        try:
            run_CDHIT(f'resources/Data/FASTA/{args.hmm_db_name}/{config["alignment_method"].upper()}/{thresh}.fasta', f'resources/Data/FASTA/{args.hmm_db_name}/CDHIT/cd-hit_after_{config["alignment_method"]}_{thresh}.fasta', 8)
            handle = cdhit_parser(f'resources/Data/FASTA/{args.hmm_db_name}/CDHIT/cd-hit_after_{config["alignment_method"]}_{thresh}.fasta.clstr')
            handle2 = counter(handle, tsv_ready = True, remove_duplicates = True)
            save_as_tsv(handle2, f'resources/Data/Tables/{args.hmm_db_name}/CDHIT_clusters/cdhit_clusters_{thresh}_after{config["alignment_method"]}.tsv')

            if config["verbose"]:
                print("Retrieving sequences divided by clusters from CDHIT\n")
            fasta_retriever_from_cdhit(f'resources/Data/Tables/{args.hmm_db_name}/CDHIT_clusters/cdhit_clusters_{thresh}_after{config["alignment_method"]}.tsv', 
                                        f'resources/Data/FASTA/{args.hmm_db_name}/CDHIT/{thresh}')
        except Exception as exc:
            print(exc)
            if config["verbose"]:
                print(f'[WARNING] Minimum cutoff of {thresh} of similarity not detected.\n')
            time.sleep(2)
            continue

    ### add cluster per threshol to config
    os.remove("config/config.yaml")

    # 
    if config_format == "yaml":
        files = get_tsv_files(config)
        threshandclust = threshold2clusters(files)
        print(threshandclust)
        for thresh, cluster in threshandclust.items():
            for c in range(len(cluster)):
                cluster[c] = str(cluster[c])
            config[thresh] = cluster
        newthresh = []
        for thresh in config["thresholds"]:
            if thresh not in threshandclust.keys():
                continue
            else:
                newthresh.append(thresh)
        config["thresholds"] = newthresh

    with open("config/config.yaml", "w") as dump_file:
        yaml.dump(config, dump_file)
        dump_file.close()

    snakemake.main(
        f'-s {args.snakefile} --printshellcmds --cores {config["threads"]} --configfile config/{args.config_file}'
        f'{" --unlock" if args.unlock else ""}')

    files = [f for f in os.listdir('.') if os.path.isfile(f)]
    for file in files:
        if file.endswith(".dnd"):
            delete_inter_files(file)

    print("HMM database created!")
    time.sleep(2)


def build_hmms_from_seqs(sequences: list, 
                    type_seq: str = "AA", 
                    from_database: str = "KEGG", 
                    ident_perc: float = 0.7):
    """Sequence that condenses the code responsible for building the hmm from the sequences received from the user.

    Args:
        sequences (list): a list of sequences. Required
        type_seq (str, optional): the type of the given list of sequences. Defaults to "AA"
        from_database (str, optional): The database the sequences come from. Defaults to "KEGG"
        ident_perc (float, optional): The identity threshold to be used for CDHIT clustering. Defaults to 0.7
    """
    if from_database == "KEGG" or from_database == "InP":
        filename = sequences.stem
    else:
        filename = sequences.split("/")[-1].split(".")[0]
    
    CDHIT_parser.run_CDHIT(sequences, PathManager.cdhit_path / Path(filename).with_suffix(".fasta"), args.threads, identperc=ident_perc, type_seq=type_seq)
    if from_database == "KEGG":
        seqs = CDHIT_parser.cdhit_parser(PathManager.cdhit_path / Path(filename).with_suffix(".fasta.clstr"), kegg = True)
    elif from_database == "InP":
        seqs = CDHIT_parser.cdhit_parser(PathManager.cdhit_path / Path(filename).with_suffix(".fasta.clstr"), ip = True)
    input_ids = parse_fasta(sequences, kegg = True, verbose = args.verbose)
    CDHIT_parser.get_clustered_sequences(seqs, PathManager.cdhit_path / "clusters", sequences, input_ids, from_database)

    for file in os.listdir(PathManager.cdhit_path / "clusters"):
        try:
            if args.input_type_db_const == "nucleic":
                run_tcoffee(PathManager.cdhit_path / "clusters" / file, 
                            PathManager.tcoffee_path / Path(file.split(".")[0]).with_suffix(".clustal_aln"), 
                            type_seq = "DNA")
            elif args.input_type_db_const == "protein":
                run_tcoffee(PathManager.cdhit_path / "clusters" / file,
                            PathManager.tcoffee_path / Path(file.split(".")[0]).with_suffix(".clustal_aln"))
        except Exception as exc:
            print(exc)
            if args.verbose:
                print(f'[WARNING] T-COFFEE for file {file} not working')
            continue

    for msa in os.listdir(PathManager.tcoffee_path):
        run_hmmbuild(PathManager.tcoffee_path / msa, PathManager.hmm_database_path / Path(msa.split(".")[0]).with_suffix(".hmm"))
        if args.consensus:
            run_hmmemit(PathManager.hmm_database_path /Path(msa.split(".")[0]).with_suffix(".hmm"), 
                    PathManager.consensus_path / Path(msa.split(".")[0]).with_suffix(".fasta"))
            concat_fasta(PathManager.consensus_path, PathManager.consensus_path / "consensus")

    concat_code_hmm(args.hmm_db_name, from_database + "_model")


def validate_hmm(config):
    print("Starting HMM validation procedures...")
    time.sleep(2)

    pathing = make_paths_dic(args.hmm_db_name)
    exec_testing(thresholds = config["thresholds"], path_dictionary = pathing, database = args.negative_db)
    to_remove = hmm_filtration(pathing)
    remove_fp_models(to_remove, pathing)
    concat_final_model(pathing)
    time.sleep(2)
    print("M-PARTY has concluded model validation! Will now switch to the newlly created models (in the validated_HMM folder\n")


def annotation():
    pass


def main_pipeline(args):
    ### Clean database ###
    if args.clean:
        paths = ["Data/FASTA/", "Alignments/", "Data/HMMs/"]
        for path in paths:
            path = f'resources/{path}{args.hmm_db_name}/'
            if Path(path).exists():
                shutil.rmtree(path)
                if args.verbose:
                    print(f'Deleted path {path}')
        return

    done = False

    if args.verbose:
        def animate():
            for c in itertools.cycle(['|', '/', '-', '\\']):
                if done:
                    break
                sys.stdout.write('\rParsing input sequences IDs ' + c)
                sys.stdout.flush()
                time.sleep(0.1)
            time.sleep(0.5)
            sys.stdout.write('\nDone!\n')

        t = threading.Thread(target=animate)
        t.start()

    if args.config_file is not None:
        config, config_format = read_config(args.config_file)
    else:
        write_config(args.input, args.output)
        config, config_format = read_config(f"{PathManager.config_path}/config.yaml")

    done = True
    time.sleep(1)

    hmmsearch_results_path = sys.path[0].replace("\\", "/") + "/results/" + args.hmm_db_name + "/HMMsearch_results/"
    Path(hmmsearch_results_path).mkdir(parents = True, exist_ok = True)

    st = time.time()

    ### VALIDATION ###
    # first only runs for if user flags --validation alone without input sequences, will validate the models inside database only
    if args.hmm_validation and args.workflow != "database_construction" and args.workflow != "both" and args.input == None:

        validate_hmm(config=config)

    ### ANNOTATION ###
    # runs if input sequences are given
    if args.workflow == "annotation" and args.input is not None:

        print("Annotation workflow started...\n")
        time.sleep(2)

        if args.hmm_validation:

            if not os.path.exists(PathManager.validated_hmm_dir):
                validate_hmm(config=config)
            else:
                print("Validated HMM already up, proceding to annotation...\n")
                time.sleep(1)

        # if a metagenome is given, runs KMA
        if args.input_type == "metagenome":
            Path(f'resources/Data/Tables/{args.hmm_db_name}/kma_hits/').mkdir(parents = True, exist_ok = True)
            Path(f'resources/Data/FASTA/DataBases/{args.hmm_db_name}/kma_db/KEGG_cons/').mkdir(parents = True, exist_ok = True)
            if args.consensus:
                kma_out = run_KMA(f'resources/Data/FASTA/{args.hmm_db_name}/Consensus/consensus.fasta', f'resources/Data/FASTA/DataBases/{args.hmm_db_name}/kma_db/',
                    args.input, f'resources/Data/Tables/{args.hmm_db_name}/kma_hits/{args.input.split(".")[0]}', threads = args.threads)
            else:
                for file in os.listdir(f'resources/Data/FASTA/{args.hmm_db_name}'):
                    if os.path.isfile(os.path.join(f'resources/Data/FASTA/{args.hmm_db_name}', file)):
                        kma_out = run_KMA(f'resources/Data/FASTA/{args.hmm_db_name}/{file}', f'resources/Data/FASTA/DataBases/{args.hmm_db_name}/kma_db/',
                            args.input, f'resources/Data/Tables/{args.hmm_db_name}/kma_hits/{args.input.split("/")[-1].split(".")[0]}', threads = args.threads)
                        
            df = kma_parser(kma_out + ".res")
            hit_seqs = get_hit_sequences(df, to_list = True)
            generate_output_files(df, hit_seqs, kma_out, config, kma = True, kma_alignfile = kma_out + ".fsa")

        # if input file is not a metagenome
        else:
            if args.hmm_validation:
                for hmm_file in file_generator(PathManager.validated_hmm_dir, full_path = True):
                    run_hmmsearch(args.input, hmm_file,
                                hmmsearch_results_path + "search_" + config["input_file"].split("/")[-1].split(".")[0] +
                                "_" + hmm_file.split("/")[-1].split(".")[0] + "." + args.hmms_output_type, verbose = args.verbose, eval = 0.00001,
                                out_type = args.hmms_output_type)
            else:
            # if models have been concatenated
                if args.concat_hmm_models:
                    print(PathManager.hmm_database_path)
                    for hmm_file in file_generator(PathManager.hmm_database_path + "concat_model/", full_path = True):
                        if os.path.exists(hmmsearch_results_path + "search_" + args.input.split("/")[-1].split(".")[0] +
                                "_" + hmm_file.split("/")[-1].split(".")[0] + "." + args.hmms_output_type):
                            os.remove(hmmsearch_results_path + "search_" + args.input.split("/")[-1].split(".")[0] +
                                "_" + hmm_file.split("/")[-1].split(".")[0] + "." + args.hmms_output_type)
                            run_hmmsearch(args.input, hmm_file, 
                                hmmsearch_results_path + "search_" + args.input.split("/")[-1].split(".")[0] +
                                "_" + hmm_file.split("/")[-1].split(".")[0] + "." + args.hmms_output_type, verbose = args.verbose, eval = 0.00001,
                                out_type = args.hmms_output_type)
                        else:
                            run_hmmsearch(args.input, hmm_file, 
                                hmmsearch_results_path + "search_" + args.input.split("/")[-1].split(".")[0] +
                                "_" + hmm_file.split("/")[-1].split(".")[0] + "." + args.hmms_output_type, verbose = args.verbose, eval = 0.00001,
                                out_type = args.hmms_output_type)
                else:
                    p = os.listdir(PathManager.hmm_database_path)
                    for thresh in p:
                        path = os.path.join(PathManager.hmm_database_path, thresh)
                        Path(path).mkdir(parents = True, exist_ok = True)
                        for hmm_file in file_generator(path, full_path = True):
                            run_hmmsearch(args.input, hmm_file, 
                                path + "search_" + config["input_file"].split("/")[-1].split(".")[0] +
                                "_" + hmm_file.split("/")[-1].split(".")[0] + "." + args.hmms_output_type, verbose = args.verbose, eval = 0.00001,
                                out_type = args.hmms_output_type)
                        concat_hmmsearch_results(path, hmmsearch_results_path)

            if args.expansion:
                lista_dataframes = dict.fromkeys(config["thresholds"])
                for file in file_generator(hmmsearch_results_path):
                    thresh = file.split("_")[-1].split(".")[0]
                    lista_dataframes[thresh] = read_hmmsearch_table(hmmsearch_results_path + file)

                final_df = concat_df_byrow(df_dict = lista_dataframes)
                rel_df = relevant_info_df(final_df)
                quality_df, bs_thresh, eval_thresh = quality_check(rel_df, give_params = True)
                hited_seqs = get_match_ids(quality_df, to_list = True, only_relevant = True)
                generate_output_files(quality_df, hited_seqs, args.input, config, bs_thresh, eval_thresh)

            else:
                for file in file_generator(hmmsearch_results_path):
                    if args.input.split("/")[-1].split(".")[0] in file:
                        dataframe = read_hmmsearch_table(hmmsearch_results_path + file)
                rel_df = relevant_info_df(dataframe)
                quality_df, bs_thresh, eval_thresh = quality_check(rel_df, give_params = True)
                hited_seqs = get_match_ids(quality_df, to_list = True, only_relevant = True)
                generate_output_files(quality_df, hited_seqs, args.input, config, bs_thresh, eval_thresh)

    ### DATABASE CONSTRUCTION ###
    elif args.workflow == "database_construction":
        database_construction(config=config)

    ### DATABASE CONSTRUCTION + ANNOTATION ###
    elif args.workflow == "both":

        database_construction(config=config)

        print("Annotation workflow started...\n")
        time.sleep(2)

        if args.input_type == "metagenome":
            Path(f'resources/Data/Tables/{args.hmm_db_name}/kma_hits/').mkdir(parents = True, exist_ok = True)
            Path(f'resources/Data/FASTA/DataBases/{args.hmm_db_name}/kma_db/KEGG_cons/').mkdir(parents = True, exist_ok = True)
            if args.consensus:
                kma_out = run_KMA(f'resources/Data/FASTA/{args.hmm_db_name}/Consensus/consensus.fasta', f'resources/Data/FASTA/DataBases/{args.hmm_db_name}/kma_db/',
                    args.input, f'resources/Data/Tables/{args.hmm_db_name}/kma_hits/{args.input.split(".")[0]}', threads = args.threads)
            else:
                for file in os.listdir(f'resources/Data/FASTA/{args.hmm_db_name}'):
                    if os.path.isfile(os.path.join(f'resources/Data/FASTA/{args.hmm_db_name}', file)):
                        kma_out = run_KMA(f'resources/Data/FASTA/{args.hmm_db_name}/{file}', f'resources/Data/FASTA/DataBases/{args.hmm_db_name}/kma_db/',
                            args.input, f'resources/Data/Tables/{args.hmm_db_name}/kma_hits/{args.input.split("/")[-1].split(".")[0]}', threads = args.threads)
                        
            df = kma_parser(kma_out + ".res")
            hit_seqs = get_hit_sequences(df, to_list = True)
            generate_output_files(df, hit_seqs, kma_out, config, kma = True, kma_alignfile = kma_out + ".fsa")

        # if input file is not a metagenome
        else:
            if args.hmm_validation:
                for hmm_file in file_generator(PathManager.validated_hmm_dir, full_path = True):
                    run_hmmsearch(args.input, hmm_file,
                                hmmsearch_results_path + "search_" + config["input_file"].split("/")[-1].split(".")[0] +
                                "_" + hmm_file.split("/")[-1].split(".")[0] + "." + args.hmms_output_type, verbose = args.verbose, eval = 0.00001,
                                out_type = args.hmms_output_type)
            else:
            # se os modelos estiverem concatenados
                if args.concat_hmm_models:
                    # pass
                    for hmm_file in file_generator(PathManager.hmm_database_path + "concat_model/", full_path = True):
                        if os.path.exists(hmmsearch_results_path + "search_" + args.input.split("/")[-1].split(".")[0] +
                                "_" + hmm_file.split("/")[-1].split(".")[0] + "." + args.hmms_output_type):
                            os.remove(hmmsearch_results_path + "search_" + args.input.split("/")[-1].split(".")[0] +
                                "_" + hmm_file.split("/")[-1].split(".")[0] + "." + args.hmms_output_type)
                            run_hmmsearch(args.input, hmm_file, 
                                hmmsearch_results_path + "search_" + args.input.split("/")[-1].split(".")[0] +
                                "_" + hmm_file.split("/")[-1].split(".")[0] + "." + args.hmms_output_type, verbose = args.verbose, eval = 0.00001,
                                out_type = args.hmms_output_type)
                        else:
                            run_hmmsearch(args.input, hmm_file, 
                                hmmsearch_results_path + "search_" + args.input.split("/")[-1].split(".")[0] +
                                "_" + hmm_file.split("/")[-1].split(".")[0] + "." + args.hmms_output_type, verbose = args.verbose, eval = 0.00001,
                                out_type = args.hmms_output_type)
                else:
                    p = os.listdir(PathManager.hmm_database_path)
                    for thresh in p:
                        path = os.path.join(PathManager.hmm_database_path, thresh)
                        Path(path).mkdir(parents = True, exist_ok = True)
                        for hmm_file in file_generator(path, full_path = True):
                            run_hmmsearch(args.input, hmm_file, 
                                path + "search_" + config["input_file"].split("/")[-1].split(".")[0] +
                                "_" + hmm_file.split("/")[-1].split(".")[0] + "." + args.hmms_output_type, verbose = args.verbose, eval = 0.00001,
                                out_type = args.hmms_output_type)
                        concat_hmmsearch_results(path, hmmsearch_results_path)

            if args.expansion:
                lista_dataframes = dict.fromkeys(config["thresholds"])
                for file in file_generator(hmmsearch_results_path):
                    thresh = file.split("_")[-1].split(".")[0]
                    lista_dataframes[thresh] = read_hmmsearch_table(hmmsearch_results_path + file)

                final_df = concat_df_byrow(df_dict = lista_dataframes)
                rel_df = relevant_info_df(final_df)
                quality_df, bs_thresh, eval_thresh = quality_check(rel_df, give_params = True)
                hited_seqs = get_match_ids(quality_df, to_list = True, only_relevant = True)
                generate_output_files(quality_df, hited_seqs, args.input, config, bs_thresh, eval_thresh)

            else:
                for file in file_generator(hmmsearch_results_path):
                    if args.input.split("/")[-1].split(".")[0] in file:
                        dataframe = read_hmmsearch_table(hmmsearch_results_path + file)
                rel_df = relevant_info_df(dataframe)
                quality_df, bs_thresh, eval_thresh = quality_check(rel_df, give_params = True)
                hited_seqs = get_match_ids(quality_df, to_list = True, only_relevant = True)
                generate_output_files(quality_df, hited_seqs, args.input, config, bs_thresh, eval_thresh)

    elif args.workflow != "annotation" and args.workflow != "database_construction" and args.workflow != "both":
        raise ValueError("-w worflow flag only ranges from 'annotation', 'database_construction' or 'both'. Chose one from the list.")


    et = time.time()
    elapsed_time = et - st
    elapsed_time = elapsed_time * 1000
    minutes_time = (elapsed_time * 1000) / 60
    print(f'Execution time: {elapsed_time:.4f} milliseconds and {minutes_time:.2f} minutes!')
    if args.workflow == "database_construction":
        if args.consensus:
            print("Consensus sequences generated!")
        else:
            print("HMMs generated!")
            print("Next step should be annotation with the just created HMM database.\nIf you need further guidance, refer to M-PARTY documentation in GitHub")
    else:
        print(f'M-PARTY has stoped running! Results are displayed in the {args.output} folder :)')
    print("Thank you for using M-PARTY! ")


if __name__ == "__main__":
    # get CLI arguments
    parser = get_parser()
    args = parser.parse_args()
    print(args)

    process_arguments(args)

    # check arguments
    check_config(args=args)

    # initialize paths
    declare_fixed_paths(args)

    # start pipeline
    main_pipeline(args)