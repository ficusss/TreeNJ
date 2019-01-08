import unittest
import pandas as pd
#import numpy as np
from networkx import is_isomorphic
from Bio import Phylo
from treenj import TreeNJ
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor, _DistanceMatrix


class BaseTestTreeNJ(unittest.TestCase):
    
    def test_set_np_matrix(self):
        pass
    
    def test_set_pd_matrix(self):
        dist_matrix = pd.read_csv("data/wiki_tree.csv", index_col=0)
        tree = TreeNJ()
        self.assertTrue(tree.set_distance_matrix(dist_matrix))

    def test_set_not_zero_diagonal_matrix(self):
        dist_matrix = pd.read_csv("data/no_zero_diagonal.csv", index_col=0)
        tree = TreeNJ()
        self.assertFalse(tree.set_distance_matrix(dist_matrix))
    
    def test_set_not_square_matrix(self):
        dist_matrix = pd.read_csv("data/no_square.csv", index_col=0)
        tree = TreeNJ()
        self.assertFalse(tree.set_distance_matrix(dist_matrix))

    def test_empty(self):
        tree = TreeNJ()
        self.assertFalse(tree.fit())

    def test_correct_work(self):    
        dist_matrix = pd.read_csv("data/wiki_tree.csv", index_col=0)
        tree = TreeNJ(dist_matrix)
        tree.fit()
        self.assertTrue(tree.get_tree())
    
    def test_correct_res(self):
        
        dist_matrix = pd.read_csv("data/wiki_tree.csv", index_col=0)
        tree1 = TreeNJ(dist_matrix)
        tree1.fit()
        
        dist_matrix = _DistanceMatrix(names=['a', 'b', 'c', 'd', 'e'], 
                                      matrix=[[0],
                                              [5, 0],
                                              [9, 10, 0],
                                              [9, 10, 8, 0],
                                              [8,  9, 7, 3, 0]])
        constructor = DistanceTreeConstructor()
        tree2 = constructor.nj(dist_matrix)
        
        self.assertTrue(is_isomorphic(Phylo.to_networkx(tree2).to_undirected(),
                                       Phylo.to_networkx(tree1.get_tree()).to_undirected()))


if __name__ == '__main__':
    unittest.main()
