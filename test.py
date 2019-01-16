import unittest
import pandas as pd

from networkx import is_isomorphic
from Bio import Phylo
from treenj import TreeNJ
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor, _DistanceMatrix


class BaseTestTreeNJ(unittest.TestCase):
    def setUp(self):
        self.tree = TreeNJ()

    def test_set_pd_matrix(self):
        dist_matrix = pd.read_csv("data/wiki_tree.csv", index_col=0)
        self.assertTrue(self.tree.set_distance_matrix(dist_matrix))

    def test_set_not_zero_diagonal_matrix(self):
        dist_matrix = pd.read_csv("data/no_zero_diagonal.csv", index_col=0)
        self.assertFalse(self.tree.set_distance_matrix(dist_matrix))
    
    def test_set_not_square_matrix(self):
        dist_matrix = pd.read_csv("data/no_square.csv", index_col=0)
        self.assertFalse(self.tree.set_distance_matrix(dist_matrix))

    def test_empty(self):
        self.assertFalse(self.tree.fit())

    def test_correct_work(self):    
        dist_matrix = pd.read_csv("data/wiki_tree.csv", index_col=0)
        self.tree.set_distance_matrix(dist_matrix)
        self.tree.fit()
        self.assertTrue(self.tree.get_tree())
    
    def test_correct_res(self):
        
        dist_matrix = pd.read_csv("data/wiki_tree.csv", index_col=0)
        self.tree.set_distance_matrix(dist_matrix)
        self.tree.fit()
        
        dist_matrix = _DistanceMatrix(names=['a', 'b', 'c', 'd', 'e'], 
                                      matrix=[[0],
                                              [5, 0],
                                              [9, 10, 0],
                                              [9, 10, 8, 0],
                                              [8,  9, 7, 3, 0]])
        constructor = DistanceTreeConstructor()
        lib_tree = constructor.nj(dist_matrix)
        
        self.assertTrue(is_isomorphic(Phylo.to_networkx(lib_tree).to_undirected(),
                                       Phylo.to_networkx(self.tree.get_tree()).to_undirected()))


if __name__ == '__main__':
    unittest.main()
