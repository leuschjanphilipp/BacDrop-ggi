import numpy as np

def calc_sparsity(matrix):
    bool_matrix = np.array(matrix, dtype=bool)
    return 1 - (matrix.sum() / np.array(matrix.shape).prod())