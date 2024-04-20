import numpy as np

from Bio import Entrez
from tqdm import tqdm
from scipy.stats import median_abs_deviation

def calc_sparsity(matrix):
    bool_matrix = np.array(matrix, dtype=bool)
    return 1 - (bool_matrix.sum() / np.array(bool_matrix.shape).prod())


def fetch_protein_names(accession_numbers):
    Entrez.email = "j.leusch@tum.de"
    handle = Entrez.efetch(db="protein", id=accession_numbers, rettype="fasta")
    records = handle.readlines()
    names = []
    for record in tqdm(records):
        if record.startswith(">"):
            names.append(record.strip().lstrip(">"))
    return names

def is_outlier(adata, metric: str, nmads: int):
    M = adata.obs[metric]
    outlier = (M < np.median(M) - nmads * median_abs_deviation(M)) | (
        np.median(M) + nmads * median_abs_deviation(M) < M
    )
    return outlier