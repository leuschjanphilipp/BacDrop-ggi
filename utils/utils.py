import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

from Bio import Entrez
from scipy.stats import median_abs_deviation

def calc_sparsity(matrix):
    bool_matrix = np.array(matrix, dtype=bool)
    return 1 - (bool_matrix.sum() / np.array(bool_matrix.shape).prod())


def fetch_protein_names(accession_numbers):
    Entrez.email = "j.leusch@tum.de"
    handle = Entrez.efetch(db="protein", id=accession_numbers, rettype="fasta")
    records = handle.readlines()
    names = []
    for record in records:
        if record.startswith(">"):
            names.append(record.strip().lstrip(">"))
    return names

def is_outlier(adata, metric: str, nmads: int):
    M = adata.obs[metric]
    outlier = (M < np.median(M) - nmads * median_abs_deviation(M)) | (
        np.median(M) + nmads * median_abs_deviation(M) < M
    )
    return outlier

def plot_lambda1_path(dict):
    lambda1 = dict["reg_params"]["lambda1"]
    lambda1_range = dict["modelselect_params"]["lambda1_range"]
    eBIC = dict["modelselect_stats"]["BIC"][0.1]

    sns.lineplot(x=lambda1_range, y=eBIC.squeeze())
    plt.axvline(x=lambda1, ls="--", color="C3", label=f"lambda1={np.round(lambda1, 7)}")
    plt.xscale("log")
    plt.xlabel("lambda")
    plt.ylabel("BIC")
    plt.gca().invert_xaxis()
    plt.legend()
    plt.show()

def plot_sparsity_path(dict):
    lambda1 = dict["reg_params"]["lambda1"]
    lambda1_range = dict["modelselect_params"]["lambda1_range"]
    sparsity = dict["modelselect_stats"]["SP"]

    sns.lineplot(x=lambda1_range, y=sparsity.squeeze())
    plt.axvline(x=lambda1, ls="--", color="C3", label=f"selected lambda1={np.round(lambda1, 7)}")
    plt.xscale("log")
    plt.xlabel("lambda")
    plt.ylabel("Sparsity")
    plt.gca().invert_xaxis()
    plt.legend()
    plt.show()