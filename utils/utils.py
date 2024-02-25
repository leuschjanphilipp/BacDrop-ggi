import numpy as np

def calc_sparsity(matrix):
    bool_matrix = np.array(matrix, dtype=bool)
    return matrix.sum() / np.array(matrix.shape).prod()