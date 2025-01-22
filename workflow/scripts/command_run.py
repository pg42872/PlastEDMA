from subprocess import run
import sys
import os


def run_command(bash_command, output='', mode='w', sep=' ', print_message=True, verbose=True):
    if print_message:
        print(f"{bash_command.replace(sep, ' ')}{' > ' + output if output != '' else ''}")
    if output == '':
        run(bash_command.split(sep), stdout=sys.stdout if verbose else None, check=True)
    else:
        with open(output, mode) as output_file:
            run(bash_command.split(sep), stdout=output_file)

def docker_run_tcoffee(volume, input_file, output_type, output_name):
    run_command(f'docker`run`--rm`-v`{volume}`pegi3s/tcoffee:latest`t_coffee`/data/{input_file}`-run_name`/data/{output_name}`-output`{output_type}', sep="`")

def docker_run_hmmbuild(volume, input_file, output_file):
    run_command(f'docker`run`--rm`-v`{volume}`biocontainers/hmmer:v3.2.1dfsg-1-deb_cv1`hmmbuild`{output_file}`{input_file}', sep="`")

def docker_run_hmmsearch(volume, hmm_file, db, output_file):
    run_command(f'docker`run`--rm`-v`{volume}`biocontainers/hmmer:v3.2.1dfsg-1-deb_cv1`hmmsearch`{hmm_file}`{db}`>`{output_file}', sep="`")

def run_tcoffee(input: str, output: str, type_seq: str = "PROTEIN"):
    run_command(f't_coffee`{input}`-output`clustalw_aln`-outfile`{output}`-type`{type_seq}`-n_core`4', sep = "`")

def run_hmmbuild(input: str, output: str):
    run_command(f'hmmbuild`{output}`{input}', sep = "`")

def run_hmmemit(input: str, output: str):
    run_command(f'hmmemit`-o`{output}`{input}', sep = "`")

def concat_hmm(input_path: str, output_path: str):
    run_command(f'cat`{input_path}*.hmm`>`{output_path}.hmm', sep = "`")

def concat_fasta(input_path: str, output_path: str):
    run_command(f'cat`{input_path}*.fasta`>`{output_path}.fasta', sep = "`")