import numpy as np
from Bio import Entrez
from tqdm import tqdm

def calc_sparsity(matrix):
    bool_matrix = np.array(matrix, dtype=bool)
    return 1 - (matrix.sum() / np.array(matrix.shape).prod())


def fetch_protein_names(accession_numbers):
    Entrez.email = "j.leusch@tum.de"
    handle = Entrez.efetch(db="protein", id=accession_numbers, rettype="fasta")
    records = handle.readlines()
    names = []
    for record in tqdm(records):
        if record.startswith(">"):
            names.append(record.strip().lstrip(">"))
    return names