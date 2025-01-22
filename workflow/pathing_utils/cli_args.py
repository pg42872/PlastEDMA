import argparse
import os

def get_parser():
    version = "1.0.1"

    parser = argparse.ArgumentParser(description="M-PARTY's main script")
    parser.add_argument("-i", "--input", help = "input FASTA file containing\
                        a list of protein sequences to be analysed")
    parser.add_argument("--input_seqs_db_const", help = "input a FASTA file with a set of sequences from which the user \
                        wants to create the HMM database from scratch.")
    parser.add_argument("-db", "--database", help = "FASTA database to run against the also user inputted sequences. DIAMOND \
                        is performed in order to expand the data and build the models. M-PARTY has no in-built database for this \
                        matter. If flag is given, download of the default database will start and model built from that. Defaults to UniProt DataBase.",
                        default = "UniProt")
    parser.add_argument("--hmm_db_name", help = "name to be assigned to the hmm database to be created. Its recomended to give a name that \
                        that describes the family or other characteristic of the given sequences. Be carefull as what name to use, as this will \
                        define what HMMs will be used for the search")
    parser.add_argument("-it", "--input_type", default = "protein", help = "specifies the nature of the sequences in the input file between \
                        'protein', 'nucleic' or 'metagenome'. Defaults to 'protein'")
    parser.add_argument("--input_type_db_const", help = "specifies the nature of the input sequences for the database construction between \
                        'nucleic' and 'protein'. Defaults to 'protein'.", default = "protein")
    parser.add_argument("--consensus", default = False, action = "store_true", help = "call to build consensus sequences when building the database, \
                        in order to run KMA against raw metagenomes")
    parser.add_argument("-o", "--output", default = "MPARTY_results", help = "name for the output directory. Defaults to 'MPARTY_results'")
    parser.add_argument("--output_type", default = "tsv", help = "choose report table outpt format from 'tsv', 'csv' or 'excel'. Defaults to 'tsv'")
    parser.add_argument("-rt", "--report_text", default = False, action = "store_true", help = "decides whether to produce or not a friendly report in \
                        txt format with easy to read information")
    parser.add_argument("--hmms_output_type", default = "tsv", help = "chose output type of hmmsearch run from 'out', 'tsv' or 'pfam' format. Defaults to 'tsv'")
    parser.add_argument("--hmm_validation", default = False, action = "store_true", help = "decides whether to perform models validation and filtration with \
                        the 'leave-one-out' cross validation methods. Call to set to True. Defaults to False")
    parser.add_argument("-p", "--produce_inter_tables", default = False, action = "store_true", help = "call if user wants to save intermediate\
                        tables as parseale .csv files (tables from hmmsearch results processing)")
    parser.add_argument("--negative_db", help = "path to a user defined negative control database. Default use of human gut microbiome")
    parser.add_argument("-s", "--snakefile", help = "user defined snakemake workflow Snakefile. Defaults to '/workflow/Snakefile",
                        default = "workflow/Snakefile")
    parser.add_argument("-ex", "--expansion", default = False, action = "store_true", help = "Decides wheter to expand the interest dataset. Defaults to False.")
    parser.add_argument("--kegg", help = "input KEGG ID(s) to download respective sequences, in order to build a pHMM based on those", nargs = "+")
    parser.add_argument("--interpro", help = "input InterPro ID(s) to download the respective sequences, in order to build a pHMM based on those", nargs = "+")
    parser.add_argument("--curated", default = False, action = "store_true", help = "call to only retrieve reviewed sequences from InterPro")
    parser.add_argument("-t", "--threads", type = int, help = "number of threads for Snakemake to use. Defaults to max number of available logical CPUs.",
                        default = os.cpu_count())
    parser.add_argument("--align_method", default = "upimapi", help = "chose the alignment method for the initial sequences database expansion, between\
                        'diamond', 'blast' and 'upimapi'. Defaults to 'upimapi'")
    parser.add_argument("--aligner", default = "tcoffee", help = "chose the aligner program to perform the multiple sequence alignment for the models\
                        between 'tcoffee' and 'muscle'. Defaults to 'tcoffee'.")
    parser.add_argument("-hm", "--hmm_models", type=str, help = "path to a directory containing HMM models previously created by the user. By default, M-PARTY\
                        does not have any in-built HMMs, so the user always needs to either create a database with the database construction workflow or \
                        inputing them this way.")
    parser.add_argument("--concat_hmm_models", action = "store_false", default = True, help = "call to not concatenate HMM models into a single file. Defaults to True")
    parser.add_argument("--unlock", action = "store_true", default = False, help = "could be required after forced workflow termination")
    parser.add_argument("-w", "--workflow", default = "annotation", help = 'defines the workflow to follow,\
                        between "annotation", "database_construction" and "both". Latter keyword makes the database construction\
                        first and posterior annotation. Defaults to "annotation"')
    parser.add_argument("-c", "--config_file", help = "user defined config file. Only recommended for\
                        advanced users.", default = None)
    parser.add_argument("--clean", default = False, action = "store_true", help = "could be required after running tool multiple times and files inside \
                        databases start to mix up. Defaults to False")
    parser.add_argument("--overwrite", action = "store_true", default = False, help = "Call to overwrite inputted files. Defaults to False")
    parser.add_argument("--verbose", action = "store_true", default = False, help = "Call so M-PARTY display more messaging")
    parser.add_argument("--display_config", default = False, action = "store_true", 
                        help = "declare to output the written config file together with results. Useful in case of debug")
    parser.add_argument("-v", "--version", action = "version", version = "M-PARTY {}".format(version))

    return parser

def process_arguments(args):
    if args.clean and args.hmm_db_name is None:
        raise ValueError("You need to provide the name of the folder to be deleted")