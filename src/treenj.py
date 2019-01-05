import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from Bio.Phylo import BaseTree, draw_graphviz


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
        assert(self.dist_matrix is not None)
        self.tree = None
        
        self._nodes = {n: BaseTree.Clade(None, str(n))
                        for n in self.dist_matrix.columns}
        self._d_matrix = self.dist_matrix
        
        while self._d_matrix.shape[0] > 2:
            
            self._update_q_matrix()
            
            raw_min, col_min = self._get_pos_min_from_q_matrix()
            range_raw, range_col = self._get_dist_for_neighborhood(raw_min, col_min)
            new_name = "{}{}".format(str(raw_min), str(col_min))
            
            self._update_nodes(raw_min, col_min, range_raw, range_col, new_name)
            
            new_dist_matrix = self._get_new_dist_matrix(raw_min, col_min)
            dists_node = [self._get_dist_for_nodes(raw_min, col_min, index)
                        for index in new_dist_matrix.index]
            new_item = pd.Series(dists_node, name=new_name, index = new_dist_matrix.index)
            new_dist_matrix = new_dist_matrix.append(new_item).T
            new_item = new_item.append(pd.Series(0, name=new_name, index=[new_name]))
            new_dist_matrix = new_dist_matrix.append(new_item)
            self._d_matrix = new_dist_matrix
            
        assert(len(self._nodes) == 2)
        name1 = self._d_matrix.index[0]
        name2 = self._d_matrix.index[1]
        node1 = self._nodes.pop(name1)
        node2 = self._nodes.pop(name2)
        node1.branch_length = self._d_matrix[name1][name2]
        node2.clades.append(node1)
        self.tree = BaseTree.Tree(node2, rooted=False)
        

    def get_tree(self):
        return (self.tree)
    
    def draw(self, filename=None):
        if (self.tree is None):
            self.fit()
        
        draw_graphviz(self.tree)
        
        if filename:
            plt.savefig(filename)
            
    
    def _get_dist_for_nodes(self, r, c, i):
        return ((self._d_matrix[r][i] + self._d_matrix[c][i] - self._d_matrix[r][c]) / 2)
    
    def _get_dist_for_neighborhood(self, r, c):
        range_r = ( self._d_matrix[r][c] / 2 
               + (self._d_matrix[r].sum() - self._d_matrix[c].sum())
               / (2 * (self._d_matrix.shape[0] - 2)))
        range_c = self._d_matrix[r][c] - range_r
        
        return (range_r, range_c)
    
    
    def _update_q_matrix(self):
        n = self._d_matrix.shape[0]
        q_matrix = np.zeros((n, n))
        for ii, i in enumerate(self._d_matrix.columns):
            for jj, j in enumerate(self._d_matrix.index):
                if i == j:
                    q_matrix[ii, jj] = 0
                else:
                    q_matrix[ii, jj] = ((n - 2) * self._d_matrix[i][j]
                                      -self._d_matrix[i].sum()
                                      -self._d_matrix[j].sum())
        
        self._q_matrix = q_matrix
    
    def _get_pos_min_from_q_matrix(self):
        tmp = np.argmin(self._q_matrix)
        raw = tmp // self._d_matrix.shape[0]
        column = tmp % self._d_matrix.shape[0]
        
        r = self._d_matrix.index[raw]
        c = self._d_matrix.index[column]
        return(r, c)
    
    def _get_new_dist_matrix(self, r, c):
        new_dist_matrix = self._d_matrix.drop([r, c], axis=0)
        new_dist_matrix = new_dist_matrix.drop([r, c], axis=1)
        
        return(new_dist_matrix)
        
    def _update_nodes(self, n1, n2, d1, d2, new_n):
        tmp1 = self._nodes.pop(n1)
        tmp2 = self._nodes.pop(n2)
        
        tmp1.branch_length = float(d1)
        tmp2.branch_length = float(d2)
        
        self._nodes[new_n] = BaseTree.Clade(None, new_n, [tmp1, tmp2])


q = 2
    
if q == 1:
    
    dist_matrix = np.array([[  0,   1,   2,   3, 1.5],
                            [  1,   0, 2.3, 4.1,   1],
                            [  2, 2.3,   0,   3,   5],
                            [  3, 4.1,   3,   0,   4],
                            [1.5,   1,   5,   4,   0]])
                    
    tmp = TreeNJ()
    tmp.set_distance_matrix(dist_matrix)
    tmp.fit()
    print(tmp.get_tree())
    tmp.draw()


if q == 2:
    
    dist_matrix = np.array([[  0,   5,   9,   9,   8],
                            [  5,   0,  10,  10,   9],
                            [  9,  10,   0,   8,   7],
                            [  9,  10,   8,   0,   3],
                            [  8,   9,   7,   3,   0]])
                    
    dist_matrix = pd.DataFrame(dist_matrix, columns=['a', 'b', 'c', 'd', 'e'],
                               index=['a', 'b', 'c', 'd', 'e'])
    
    tmp.set_distance_matrix(dist_matrix)
    tmp.fit()
    print(tmp.get_tree())
    tmp.draw("figure.png")

