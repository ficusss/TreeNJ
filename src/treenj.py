import numpy as np
import pandas as pd

from Bio.Phylo import BaseTree


class TreeNJ:
    def __init__(self, dist_matrix=None):
        self.dist_matrix = dist_matrix
        self.tree = None
        self._nodes = None
        self._q_matrix = None
        self._d_matrix = None
        
    def set_distance_matrix(self, dist_matrix):
        if type(dist_matrix) is np.ndarray:
            dist_matrix = pd.DataFrame(dist_matrix)
        
        self.tree = None
        self.dist_matrix = dist_matrix
        
    def fit(self):
        self._nodes = {n: BaseTree.Clade(None, n)
                        for n in self.dist_matrix.columns}
        self._d_matrix = self.dist_matrix
        
        
        TreeNJ._update_q_matrix(self)
        raw_min, col_min = TreeNJ._get_pos_min_from_q_matrix(self)
        range_raw = ( self._d_matrix[raw_min][col_min] / 2 
               + (self._d_matrix[raw_min].sum() - self._d_matrix[col_min].sum())
               / (2 * (self._d_matrix.shape[0] - 2)))
        range_col = self._d_matrix[raw_min][col_min] - range_raw

        
    
    def get_tree(self):
        return self.tree
    
    def _update_q_matrix(self):
        n = self._d_matrix.shape[0]
        q_matrix = np.zeros((n, n))
        for i in self._d_matrix.columns:
            for j in self._d_matrix.index:
                if i == j:
                    q_matrix[i, j] = 0
                else:
                    q_matrix[i, j] = ((n - 2) * self._d_matrix[i][j]
                                      -self._d_matrix[i].sum()
                                      -self._d_matrix[j].sum())
        
        self._q_matrix = q_matrix
    
    def _get_pos_min_from_q_matrix(self):
        tmp = np.argmin(self._q_matrix)
        raw = tmp // self._d_matrix.shape[0]
        column = tmp % self._d_matrix.shape[0]

        return(raw, column)
    
    def _get_new_dist_matrix(self, raw, column):
        new_dist_matrix = self._d_matrix.drop([raw, column], axis=0)
        new_dist_matrix = new_dist_matrix.drop([raw, column], axis=1)
        
        return(new_dist_matrix)
        
    def _update_nodes(self, n1, n2, d1, d2, new_n):
        tmp1 = self._nodes.pop(n1)
        tmp2 = self._nodes.pop(n2)
        
        tmp1.branch_length = float(d1)
        tmp2.branch_length = float(d2)
        
        self._nodes[new_n] = BaseTree.Clade(None, new_n, [tmp1, tmp2])

    
    
