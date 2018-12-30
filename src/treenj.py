import numpy as np
import pandas as pd


class TreeNJ:
    def __init__(self, dist_matrix=None):
        self.dist_matrix = dist_matrix
        self.tree = None
        self._q_matrix = None
    
    def set_distance_matrix(self, dist_matrix):
        if type(dist_matrix) is np.ndarray:
            dist_matrix = pd.DataFrame(dist_matrix)
        
        self.tree = None
        self.dist_matrix = dist_matrix
        
    def fit(self):
        pass
    
    def get_tree(self):
        return self.tree
    
    def _get_q_matrix(self):
        pass
    
    def
    
