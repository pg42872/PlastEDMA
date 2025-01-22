from glob import glob
import sys
from itertools import product
from clint.textui import progress
import pandas as pd
import requests
import shutil
from command_run import run_command
import os
import time
from hmm_vali import delete_inter_files
from pathlib import Path
import fileinput
from tqdm import tqdm
from requests.exceptions import HTTPError
from bs4 import BeautifulSoup


def get_clusters(tsv_file: str) -> list:
    df = pd.read_csv(tsv_file, sep="\t", index_col=0)
    return list(df.index.values)


def get_number_clusters(tsv_file: str) -> int:
    df = pd.read_csv(tsv_file, sep="\t", index_col=0)
    return len(df.index)


# vai buscar aos .tsv criados antes, e não fica dependente dos fasta que vão ser criados
def get_tsv_files(config_file: str) -> dict:
	files = {threshold: glob(f'resources/Data/Tables/{config_file["hmm_database_name"]}/CDHIT_clusters/cdhit_clusters_{threshold}_after{config_file["alignment_method"].upper()}.tsv') for threshold in config_file["thresholds"]}
	return files


def threshold2clusters(file_dic: dict) -> dict:
	threshold2clusters, threshold2clustersclean = {}, {}
	for thresh, path in file_dic.items():
		try:
			threshold2clusters[thresh] = get_clusters(path[0])
		except Exception:
			print(f'No clusters were found for {thresh} threshold.')
			continue
	for k, v in threshold2clusters.items():
		if v != []:
			threshold2clustersclean[k] = v
	return threshold2clustersclean


def clusters_in_list(dic: dict) -> list:
	lista_clusters = [v for k, v in dic.items()]
	return lista_clusters


def get_all_clusters(config_file: str) -> tuple:
# fazer uma lista de listas com todos os clusters, por ordem de threshold
	big_list_clusters = [v for k, v in threshold2clusters.items()]
	max_clusters = max([max(x) for x in big_list_clusters])
	all_clusters = [str(i) for i in range(0, max_clusters+1)]
	return big_list_clusters, all_clusters, max_clusters


# função vai fazer todas as combinações entre thresholds e clusters correspondentes
def util(lista_thresholds, lista_de_listas_clusters):
	autorized_combs = []
	for threshold in range(len(lista_thresholds)):
		for cluster in lista_de_listas_clusters[threshold]:
			combinacao = (lista_thresholds[threshold], str(cluster))
			autorized_combs.append(combinacao)
	autorized_combs_frozen = {frozenset(t) for t in autorized_combs}
	return autorized_combs_frozen


# função que vai fazer o produto entre todos, e so devolve os desejados
def match_threshold_W_cluster(combinador, desired_combs) -> tuple:
    def match_threshold_W_cluster(*args, **kwargs):
        for combo in combinador(*args, **kwargs):
            if frozenset(combo) in desired_combs:
                yield combo
    return match_threshold_W_cluster


def cat_hmms_input(wildcards):
	return expand("resources/Data/HMMs/After_tcoffee_UPI/{threshold}/{cluster}.hmm", threshold=wildcards, cluster=threshold2clusters[wildcards])


def get_target_db(config):
	return config["hmm_database_name"]


def get_upi_querydb(config):
	return config["database"]


def save_as_tsv(dic: dict, out_path: str):
    int_df = pd.DataFrame.from_dict(dic, orient="index")
    int_df.to_csv(out_path, sep="\t")


def get_output_dir(path: str, config: str, hmm: bool = False) -> str:
	c = path.split("_")
	if hmm:
		ind = c.index("HMMs")
	else:
		try:
			ind = c.index("FASTA")
		except Exception:
			ind = c.index("Tables")
	c.insert(ind + 1, config["hmm_database_name"])
	return "/".join(c)


def ask_for_overwrite(path: str, verbose: bool = False) -> bool:
	"""Function will ask if the user wants to overwrite or not an existing file.

	Args:
		path (str): Path to the file to overwrite or not.
		verbose (bool, optional): Decides to print aditional information. Defaults to False.

	Raises:
		ValueError: If user insists in not giving the required input between 'y' and 'n' program will cease.

	Returns:
		bool: True if it is to overwrite. False if not.
	"""
	overwrite = input(f'[WARNING] {path} already exists - overwrite?'
                    	f' [y/n] ({path}) ')
	count = 0
	while overwrite not in ["y", "n", "Y", "N"]:
		overwrite = input("Enter 'y' (overwrite) or 'n' (cancel).")
		count += 1
		if count == 10:
			raise ValueError("Too many tries. Try again...\n")
	if overwrite.lower() == "n":
		if verbose:
			print(f'Not overwriting {path}. Using previous file.\n')
			time.sleep(0.5)
		return False
	elif overwrite.lower() == "y":
		if verbose:
			print("Proceding to overwrite present file...\n")
		print("[TIP] Next time specify --overwrite = True\n")
		return True


def download_with_progress_bar(url: str, database_folder: str):
	"""Function that builds a progress bar given an url.

	Args:
		url (str): link to get the data.
		database_folder (str): path to the folder.
	"""
	r = requests.get(url, stream=True)
	path = f'{database_folder}/{url.split("/")[-1]}'
	with open(path, "wb") as wf:
		total_length = int(r.headers.get('content-length'))
		for chunk in progress.bar(r.iter_content(chunk_size=1024), expected_size=(total_length/1024) + 1): 
			if chunk:
				wf.write(chunk)
				wf.flush()
	wf.close()


def download_uniprot(database_folder: str):
	"""will download and read the content of the compressed output, without actually decompresing.

	Args:
		database_folder (str): path to the folder.
	"""
	for url in [
	"https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz", 
	"https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_trembl.fasta.gz",
	"https://ftp.uniprot.org/pub/databases/uniprot/relnotes.txt"]:
		download_with_progress_bar(url, database_folder)
	run_command(f'zcat {database_folder}/{url[0].split("/")[-1]} {database_folder}/{url[1].split("/")[-1]} > {database_folder}/uniprot.fasta')


def build_upi_query_db(database_folder: str, config: str = None, verbose: bool = False) -> str:
	"""Function that will download the database from uniprot to a specified folder.

	Args:
		database_folder (str): name for the output folder.
		config (str, optional): path to the config file. Defaults to None.
		verbose (bool, optional): flag to print what is going on. Defaults to False.

	Raises:
		TypeError: If database to be downloaded is not unipror, swissprot or a FASTA file.

	Returns:
		str: Path for the new database.
	"""
	# database = "uniprot"
	database = get_upi_querydb(config)
	if database.lower() == "uniprot":
		if not os.path.exists(database_folder + "/uniprot.fasta"):
			print(f'Download of {database} database started...\n')
			download_uniprot(database_folder)
			if verbose:
				print("Done\n")
				time.sleep(0.5)
		else:
			overw = ask_for_overwrite(database_folder + "/uniprot.fasta.gz", verbose = verbose)
			if overw:
				os.remove(database_folder + "/uniprot.fasta")
				print(f'Download of {database} database started...\n')
				download_uniprot(database_folder)
				if verbose:
					print("Done\n")
					time.sleep(0.5)
			else:
				print("UniProt database already present. Proceding...\n")
				time.sleep(0.5)
		return database_folder + "/uniprot.fasta"
	elif database.lower() == "swissprot":
		if not os.path.exists(database_folder + "/uniprot_sprot.fasta.gz") and os.path.exists(database_folder + "/uniprot_sprot.fasta"):
			overw = ask_for_overwrite(database_folder + "/uniprot_sprot.fasta", verbose = verbose)
			if overw:
				os.remove(database_folder + "/uniprot_sprot.fasta")
				print(f'Download of {database} database started...\n')
				download_with_progress_bar("https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz", database_folder)
				run_command(f'gunzip -v {database_folder}/uniprot_sprot.fasta.gz')
				if verbose:
					print("Done\n")
					time.sleep(0.5)
			return database_folder + "/uniprot_sprot.fasta"
		if not os.path.exists(database_folder + "/uniprot_sprot.fasta.gz"):
			print(f'Download of {database} database started...\n')
			download_with_progress_bar("https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz", database_folder)
			run_command(f'gunzip -v {database_folder}/uniprot_sprot.fasta.gz')
			if verbose:
				print("Done\n")
				time.sleep(0.5)
		else:
			overw = ask_for_overwrite(database_folder + "/uniprot_sprot.fasta.gz", verbose = verbose)
			if overw:
				os.remove(database_folder + "/uniprot_sprot.fasta.gz")
				print(f'Download of {database} database started...\n')
				download_with_progress_bar("https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz", database_folder)
				run_command(f'gunzip -v {database_folder}/uniprot_sprot.fasta.gz')
				if verbose:
					print("Done\n")
					time.sleep(0.5)
			else:
				run_command(f'gunzip -v {database_folder}/uniprot_sprot.fasta.gz')
				print("SwissProt database already present. Proceding...\n")
				time.sleep(0.5)
			return database_folder + "/uniprot_sprot.fasta"
	elif database.split(".")[-1] == "fasta":
		if not os.path.exists(database_folder + "/" + database.split("/")[-1]):
			shutil.move(database, database_folder)
			if verbose:
				print(f'Inputed database moved to {database_folder}')
				time.sleep(0.5)
		else:
			overw = ask_for_overwrite(database_folder + "/" + database.split("/")[-1], verbose = verbose)
			if overw:
				os.remove(database_folder + "/" + database.split("/")[-1])
				shutil.move(database, database_folder)
				if verbose:
					print(f'Inputed database moved to {database_folder}')
					time.sleep(0.5)
			else:
				if verbose:
					print("Inputed database already present. Proceding...\n")
					time.sleep(0.5)
		return database_folder + "/" + database.split("/")[-1]
	else:
		raise TypeError("--database given parameter is not accepted. Chose between 'uniprot', 'swissprot' or a path to a FASTA file of protein sequences.")


def concat_code_hmm(db_name: str, model_name: str):
	"""Will concat the created HMMs.

	Args:
		db_name (str): name given by the user to his database.
		model_name (str): name to give to the concatenated model.
	"""
	Path(f'resources/Data/HMMs/{db_name}/concat_model/').mkdir(parents = True, exist_ok = True)
	with open(f'resources/Data/HMMs/{db_name}/concat_model/{model_name}.hmm', "w") as wf:
		for hmm in os.listdir(f'resources/Data/HMMs/{db_name}/'):
			if os.path.isfile(os.path.join(f'resources/Data/HMMs/{db_name}/', hmm)):
				with open(os.path.join(f'resources/Data/HMMs/{db_name}/', hmm), "r") as rf:
					lines = rf.readlines()
					wf.writelines(lines)
				rf.close()
	wf.close()


def delete_previous_same_run(dir_path: str):
	shutil.rmtree(dir_path)


def check_id(filepath: str, outpath: str, id_list: list):
	"""Checks the existance a number of IDs inside a list in a given file, and writes the respective found sequences in a output file.

	Args:
		filepath (str): path to the file to be checked on.
		outpath (str): path to the file to be writen, which will be added "aligned.fasta".
		id_list (list): list of IDs to find in the "filepath".
	"""
	num_lines = sum(1 for _ in open(filepath, "r"))
	sequence = ""
	in_seq = False
	with open(outpath + "aligned.fasta", "w") as wf:
		with open(filepath, "r") as rf:
			for line in tqdm(rf, desc = "Searching for matched sequences in input file (this migth take a while)", total = num_lines, unit = "B", unit_scale = True):
				if line.startswith(">") and in_seq == False:
					for id in id_list:
						if id in line:
							in_seq = True
							sequence += line
							break
				elif line.startswith(">") and in_seq == True:
					in_seq = False
					wf.write(sequence)
					sequence = ""
				elif not line.startswith(">") and in_seq == True:
					sequence += line
		rf.close()
	wf.close()


def compress_fasta(filepath: str) -> str:
	"""Given a large file, concatenate the sequences lines, so each line is a FASTA identifier followed by "|||" and the full sequence.

	Args:
		filepath (str): file to be concatenated.

	Returns:
		str: the same sime file, now with much less lines.
	"""
	separator = "|||"
	for line in fileinput.input(filepath, inplace=True):
		if line.startswith(">"):
			line = line.replace("\n", separator)
			sys.stdout.write("\n" + line)
		else:
			sys.stdout.write(line.strip())
	return filepath


def return_fasta_content(filepath: str, outpath: str, identifier: list = None):
	"""Given a compressed FASTA from compress_fasta(), return the ID and the sequence of each given identifier, in fasta format
	for a output file.

	Args:
		filepath (str): compressed file.
		outpath (str): path to the file to be writen, which will be added "aligned.fasta".
		identifier (list, optional): list of IDs to be checked. Defaults to None.
	"""
	with open(filepath, "r") as rf:
		with open(outpath + "aligned.fasta", "w") as wf:
			# uma linha do ficheiro de input em memória
			for line in rf:
				for ident in identifier:
					# se o id que deu hit estiver nessa linha
					if ident in line:
						# divide-se pelo separador
						new = line.split("|||")
						wf.write(new[0] + "\n")
						# passamos a formato de fasta
						fasta_form = []
						if len(new[1]) > 60:
							for i in range(0, len(new[1]), 60):
								fasta_form.append(new[1][i: i + 60] + "\n")
							fasta_form.append("\n")
							for x in fasta_form:
								wf.write(x)
						else:
							wf.write(new[1])


def get_soup(url: str, status: int = 200):
	"""Simplified function to retrieve the content of a given URL in HTML format with the BeatifulSoup package

	Args:
		url (str): requested URL
		status (int, optional): Defaults to 200.

	Returns:
		object: Return a BeautifulSoup object to be extracted posteriorly
	"""
	try:
		response = requests.get(url)
	except HTTPError as http_err:
		print(f'HTTP error occurred: {http_err}')
	except Exception as err:
		print(f'Other error occurred: {err}')
	soup = BeautifulSoup(response.text, "html.parser")
	return soup


def retry(tries: int, url: str):
	"""Due to server overload, it's necessary to perform some tries and to connect to the URL.

	Args:
		tries (int): number of tries.
		url (str): requested URL

	Returns:
		object: requests object with the URL's information
	"""
	connected = False
	i = 0
	while not connected and i < tries:
		try:
			response = requests.get(url, timeout = 10)
			connected = True
		except Exception as e:
			print(e)
			i += 1
			time.sleep(2)
	return response		
