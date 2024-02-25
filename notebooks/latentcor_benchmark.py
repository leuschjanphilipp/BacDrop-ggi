import time
import sys
import scanpy as sc
import numpy as np

import warnings
from statsmodels.tools.sm_exceptions import IterationLimitWarning
warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.simplefilter(action='ignore', category=IterationLimitWarning)

from latentcor import latentcor

if __name__ == "__main__":

    print('Benchmark the influence of use_nearPD argument on speed of latentcore with increasing gene and cell counts')
    n_exp = int(sys.argv[1])
    step_size = int(sys.argv[2])
    start = int(sys.argv[3])
    print(f'Running {n_exp} experiments with increasing step size of {step_size} starting at {start}')
    print('Loading data...')
    ec_adata = sc.read_h5ad(filename='../data/preprocessed/psa_prepro.h5ad')
    print('Finish!')

    times_false = []
    times_true = []

    for i in range(start, min(n_exp*step_size + start, ec_adata.shape[1]), step_size):
        
        print(f'Running trial {i/step_size} with {i} HVGs')
        sc.pp.highly_variable_genes(ec_adata, n_top_genes=i)
        ec_adata_hvg = ec_adata[:, ec_adata.var['highly_variable']].copy()

        print('Shape pre filter: ', ec_adata_hvg.shape)
        #sc.pp.filter_genes(ec_adata_hvg, min_cells=2)
        #sc.pp.filter_cells(ec_adata_hvg, min_genes=1)
        print('Shape post filter: ', ec_adata_hvg.shape) 

        tps_hvg = ['tru' for i in range(ec_adata_hvg.shape[1])]

        start_false = time.time()
        est_pd_false = latentcor(ec_adata_hvg.layers['counts'].A, tps=tps_hvg, method='approx', use_nearPD=False, nu=0.001, showplot=False)
        end_false = time.time()
        print('Sanity check: ', est_pd_false['R'][0][0])
        time_false = np.round(end_false-start_false, 2)
        times_false.append(time_false)
        print('time for use_nearPD=False:', time_false)

        start_true = time.time()
        est_pd_true = latentcor(ec_adata_hvg.layers['counts'].A, tps=tps_hvg, method='approx', use_nearPD=True,  nu=0.001, showplot=False)
        end_true = time.time()
        print('Sanity check: ', est_pd_true['R'][0][0])
        time_true = np.round(end_true-start_true, 2)
        times_true.append(time_true)
        print('time for use_nearPD=True: ', time_true)

        print('Current flase times: ', times_false)
        print('Current true times : ', times_true)

        print('saving correlation matricies')
        np.save(f'../data/benchmark_latentcor/filter_HVG_{i}_false.npy', est_pd_false)
        np.save(f'../data/benchmark_latentcor/filter_HVG_{i}_true.npy', est_pd_true)

        del(ec_adata_hvg, tps_hvg, start_false, end_false, start_true, end_true, est_pd_false, est_pd_true)

        print('__________________________________')

    np.save(f'../data/benchmark_latentcor/filter_false_times.npy', times_false)
    np.save(f'../data/benchmark_latentcor/filter_true_times.npy', times_true)
    print('Finish benchmark.')