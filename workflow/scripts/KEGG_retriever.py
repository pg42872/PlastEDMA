import requests
from requests.exceptions import HTTPError
import re
from tqdm import tqdm
import math
import time
from mparty_util import get_soup, retry


def find_between(string: str, first: str, last: str):
    """Function that, given a string, returns the content between two points.

    Args:
        string (str): whole string.
        first (str): first character
        last (str): last character

    Returns:
        str: String between first and last.
    """
    try:
        start = string.index(first) + len(first)
        end = string.index(last, start)
        return string[start:end]
    except ValueError:
        return ""


def get_kegg_genes(filepath: str, type_seq: str = "AA", ec_number: str = None, ko: str = None, verbose: bool = False):
    """Function that will connect to KEGG and find information based on what it's asked for, in distinct databases

    Args:
        filepath (str): 
        type_seq (str, optional): _description_. Defaults to "AA".
        ec_number (str, optional): _description_. Defaults to None.
        ko (str, optional): _description_. Defaults to None.
        verbose (bool, optional): _description_. Defaults to False.

    Raises:
        ValueError: _description_
        ValueError: _description_

    Returns:
        _type_: _description_
    """
    if ko == None and ec_number == None:
        raise ValueError("Either an E.C. number or a KO must be given")
    elif ko != None and ec_number != None:
        raise ValueError("Only give eith a E.C. number or KO")
    elif ko:
        ko_urlist = []
        for k in ko:
            ko_urlist.append(f'https://rest.kegg.jp/get/{k}')

            if verbose:
                print(f'Downloading {k} orthology genes')
        for link in ko_urlist:
            get_kegg_kosequences(filepath, link, korec = "ko", type_seq = type_seq, verbose=verbose)
    elif ec_number:
        ec_urlist = []
        for ec in ec_number:
            ec_urlist.append(f'https://rest.kegg.jp/get/{ec}')
            if verbose:
                print(f'Downloading {ec} genes')
        for link in ec_urlist:
            get_kegg_kosequences(filepath, link, korec = "ec", type_seq = type_seq, verbose=verbose)
    return filepath


def get_kegg_kosequences(filepath: str, url: str, korec: str = None, type_seq: str = "AA", verbose: bool = False, tries: int = 5):
    """Function that will retrieve the genes (nucleic or protein) for a given KO or EC from KEGG directly

    Args:
        filepath (str): path to the file where the output will be writen
        url (str): requested URL
        korec (str, optional): argument that is passed depending on user input, either "KO" or "EC". Defaults to None.
        type_seq (str, optional): type of sequence to be retrieved from KEGG databases, either "AA" or "NUC". Defaults to "AA".
        verbose (bool, optional): Decides to print more information. Defaults to False.
        tries (int, optional): number of tries to connect to KEGG. Defaults to 5.
    """
    soup = get_soup(url)
    genes = []
    start_line = False
    end_line = False
    with open(filepath, "a") as wf:
        for line in soup.get_text().split('\n'):
            if line.startswith("GENES"):
                start_line = True
                end_line = False
                gene_line = line.split(None, 1)[1]
                genes.extend(gene_line.split(';'))
            elif start_line and not end_line:
                if line.startswith(' '):
                    gene_line = line.strip()
                    genes.extend(gene_line.split(';'))
                else:
                    end_line = True
        # case of EC numbers with no acess to api genes
        # print(genes)
        if genes == []:
            if verbose:
                print(f'Sequences from {url} not found. Trying in RefGene')
            wf.close()
            try:
                get_kegg_refgene_sequences(filepath, url, korec, type_seq = type_seq, verbose = verbose)
            except Exception as exc:
                # print(exc)
                print(f"[WARNING] Sequence for {url} not found")
        else:
            genes2 = []
            for entry in genes:
                indiv = entry.split(" ")
                indiv[0] = indiv[0].lower()
                if len(indiv) <= 2:
                    genes2.append("".join(indiv))
                else:
                    for x in range(1, len(indiv)):
                        genes2.append(indiv[0] + indiv[x])
            for i in tqdm(range(len(genes2)), desc = f'Downloading genes for {url.split("/")[-1]}', unit = "sequence", disable=True if not verbose else False):
                genes2[i] = re.sub("\\(.*\\)", "", genes2[i])
                if type_seq == "AA":
                    url2 = f'https://rest.kegg.jp/get/{genes2[i]}/aaseq'
                else:
                    url2 = f'https://rest.kegg.jp/get/{genes2[i]}/ntseq'
                response = retry(tries, url2)
                wf.write(response.text)
            if verbose:
                print(f'Checking the existance of other genes for {url.split("/")[-1]} in RefSeq')
            get_kegg_refgene_sequences(filepath, url, korec, type_seq, verbose)
        wf.close()


def get_kegg_refgene_sequences(filepath: str, url: str, korec: str = None, type_seq: str = "AA", verbose: bool = False, tries: int = 5):
    """Function called to find more genes in the same run in other databases from KEGG

    Args:
        filepath (str): path to the file where the output will be writen
        url (str): resquest URL
        korec (str, optional): argument that is passed depending on user input, either "KO" or "EC". Defaults to None.
        type_seq (str, optional): type of sequence to be retrieved from KEGG databases, either "AA" or "NUC". Defaults to "AA".
        verbose (bool, optional): Decides to print more information. Defaults to False.
        tries (int, optional): number of tries to connect to KEGG. Defaults to 5.
    """
    url_list = []
    if korec == "ec":
        ec = url.split("/")[-1]
        url = f"https://www.genome.jp/dbget-bin/get_linkdb?-t+refgene+ec:{ec}"
        soup = get_soup(url)
        try:
            number_hits = soup.find(string = re.compile("Hits")) 
            number_pages = math.ceil(int(re.findall("[0-9]{2,}", number_hits)[0]) / 1000)
        except:
            return
        if number_pages < 2:
            url_list.append(url)
        else:
            for i in range(1, number_pages + 1):
                url_list.append(f"https://www.genome.jp/dbget-bin/get_linkdb?-t+refgene+-p+{i}+ec:{ec}")
    else:
        ko = url.split("/")[-1]
        url = f"https://www.genome.jp/dbget-bin/get_linkdb?-t+refgene+ko:{ko}"
        soup = get_soup(url)
        try:
            number_hits = soup.find(string = re.compile("Hits"))
            number_pages = math.ceil(int(re.findall("[0-9]{2,}", number_hits)[0]) / 1000)
        except:
            return
        if number_pages < 2:
            url_list.append(url)
        else:
            for i in range(1, number_pages + 1):
                url_list.append(f"https://www.genome.jp/dbget-bin/get_linkdb?-t+refgene+-p+{i}+ko:{ko}")
    with open(filepath, "a") as wf:
        not_found = 0
        for url2 in url_list:
            soup = get_soup(url2)
            # Extract the gene IDs from the HTML
            gene_ids = [a.text for a in soup.find_all("a")]
            # Use the KEGG API to get the sequence information for a given list of gene IDs
            fasta = ""
            for gene_id in tqdm(gene_ids, desc = f'Downloading genes pag.{url_list.index(url2) + 1}/{len(url_list)} in RefGene', position = 0, leave = True, unit = "sequence", disable=True if not verbose else False):
                if gene_id.startswith("RG"):
                    if type_seq == "AA":
                        url3 = f"https://www.genome.jp/entry/-f+-n+a+{gene_id}"
                    else:
                        url3 = f"https://www.genome.jp/entry/-f+-n+n+{gene_id}"
                    response = retry(tries, url3)
                    if response.status_code == 403:
                        if verbose:
                            print("Waiting for server decongestion\n")
                        time.sleep(181)
                        response = requests.get(url3)
                    text = response.text
                    # Add the sequence information to the fasta string
                    try:
                        start = text.index("-->&gt;") + len("-->&gt;")
                        end = text.index("</pre>", start)
                        fasta += ">" + text[start:end]
                    except Exception as e:
                        # print(e)
                        # print(f"[WARNING] Sequence for {gene_id} not found")
                        not_found += 1
                        continue
            wf.write(fasta)
        if verbose:
            print(f'{not_found} IDs not found')
        wf.close()

# ec_number = ["3.1.1.101", "1.14.12.15", "1.18.1.-", "1.3.1.53"]
# ko = ["K18251", "K18252", "K18253", "K18254"]