import pickle
import numpy as np
import pandas as pd
import seaborn as sns
import networkx as nx

import matplotlib.pyplot as plt

from gglasso.problem import glasso_problem
from gglasso.helper.basic_linalg import adjacency_matrix

class gg_lasso_network_analysis():
    def __init__(self, N, file_path=None, load_data=False, estimate=None):
        self.N = N
        self.S = None
        self.P = None
        self.G = None

        if load_data:
            self.file_path = file_path
            self.estimate = None
            self.load_data()
        else:
            self.estimate = estimate

    def load_data(self):
        with open(self.file_path, 'rb') as f:
            self.estimate = pickle.load(f)
    
    def create_problem(self, S_col=None):
        if S_col:
            self.S = self.estimate[S_col].to_numpy()
        else:
            self.S = self.estimate
        self.P = glasso_problem(self.S, self.N, reg_params = {'lambda1': 0.05}, latent = False, do_scaling = False)
        print(self.P)

    def model_selection(self, lambda1_range=np.logspace(0, -3, 30), method='eBIC', gamma = 0.1):

        self.P.model_selection(modelselect_params = {'lambda1_range': lambda1_range}, method = method, gamma = gamma)
        # regularization parameters are set to the best ones found during model selection
        print(self.P.reg_params)
    
    def plot_graph_and_percision_matrix(self, fig_size=(10,10), node_size=100, font_size=9, node_color="peru", edge_color="peru", font_color='white', with_labels=True, seed=42):
        sol = self.P.solution.precision_
        self.P.solution.calc_adjacency(t = 1e-4)

        fig, axs = plt.subplots(2, figsize=fig_size)

        self.G = nx.from_numpy_array(self.P.solution.adjacency_)
        pos = nx.drawing.layout.spring_layout(self.G, seed=seed)

        nx.draw_networkx(self.G, pos=pos, node_size=node_size, node_color=node_color, edge_color=edge_color, \
                         font_size=font_size, font_color=font_color, with_labels=with_labels, ax=axs[0])
        axs[0].axis('off')
        axs[0].set_title("Recovered graph")

        sns.heatmap(sol, cmap="coolwarm", vmin=-0.5, vmax=0.5, linewidth=.5, square=True, cbar=False, \
                    xticklabels=[], yticklabels=[], ax=axs[1])
        axs[1].set_title("Recovered precision matrix")
