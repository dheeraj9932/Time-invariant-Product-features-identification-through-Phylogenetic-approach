import Bio
from Bio import AlignIO
from Bio import Phylo
from fastdtw import fastdtw
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
from Bio.Phylo.TreeConstruction import *

#DNA sequence substitution matrices
#['identity', 'blastn', 'trans', 'benner6', 'benner22', 'benner74', 'blosum100', 'blosum30', 'blosum35', 'blosum40', 'blosum45', 'blosum50', 'blosum55', 'blosum60', 'blosum62', 'blosum65',
# 'blosum70', 'blosum75', 'blosum80', 'blosum85', 'blosum90', 'blosum95', 'feng', 'fitch', 'genetic', 'gonnet', 'grant', 'ident', 'johnson', 'levin', 'mclach', 'miyata', 'nwsgappep', 'pam120',
# 'pam180', 'pam250', 'pam30', 'pam300', 'pam60', 'pam90', 'rao', 'risler', 'structure']

#Protein models
#['benner6', 'benner22', 'benner74', 'blosum100', 'blosum30', 'blosum35', 'blosum40', 'blosum45', 'blosum50', 'blosum55', 'blosum60', 'blosum62', 'blosum65', 'blosum70', 'blosum75', 'blosum80',
# 'blosum85', 'blosum90', 'blosum95', 'feng', 'fitch', 'genetic', 'gonnet', 'grant', 'ident', 'johnson', 'levin', 'mclach', 'miyata', 'nwsgappep', 'pam120', 'pam180', 'pam250', 'pam30', 'pam300',
# 'pam60', 'pam90', 'rao', 'risler', 'structure']

# aln = AlignIO.read('data/swir_converted.phy', 'phylip')
import numpy as np
from numpy import genfromtxt

my_data_swir = genfromtxt('data/SWIR.AllData_avg - Copy.csv', delimiter=',')
m1_y = my_data_swir[1][2:-1]
m1_g = my_data_swir[2][2:-1]
m1_r = my_data_swir[3][2:-1]
m1_b = my_data_swir[4][2:-1]

m2_y = my_data_swir[5][2:-1]
m2_r = my_data_swir[6][2:-1]
m2_b = my_data_swir[7][2:-1]
m2_g = my_data_swir[8][2:-1]

m3_r = my_data_swir[9][2:-1]
m3_y = my_data_swir[10][2:-1]
m3_b = my_data_swir[11][2:-1]
m3_g = my_data_swir[12][2:-1]

m4_y = my_data_swir[13][2:-1]
m4_g = my_data_swir[14][2:-1]
m4_r = my_data_swir[15][2:-1]
m4_b = my_data_swir[16][2:-1]

# aln = AlignIO.read('/home/dheeraj/CRNT/ovgu/fraunhofer_work/biopython/data/swir_converted.phy', 'phylip')
aln = AlignIO.read('/home/dheeraj/CRNT/ovgu/fraunhofer_work/biopython/data/material_pair_trees/m1/swir_converted.phy', 'phylip')

# print(aln)
DNA_models  = ['identity', 'blastn', 'trans', 'benner6', 'benner22', 'benner74', 'blosum100', 'blosum30', 'blosum35', 'blosum40', 'blosum45', 'blosum50', 'blosum55', 'blosum60', 'blosum62', 'blosum65',
               'blosum70', 'blosum75', 'blosum80', 'blosum85', 'blosum90', 'blosum95', 'feng', 'fitch', 'genetic', 'gonnet', 'grant', 'ident', 'johnson', 'levin', 'mclach', 'miyata', 'nwsgappep', 'pam120',
               'pam180', 'pam250', 'pam30', 'pam300', 'pam60', 'pam90', 'rao', 'risler', 'structure']
protein_models = ['benner6', 'benner22', 'benner74', 'blosum100', 'blosum30', 'blosum35', 'blosum40', 'blosum45', 'blosum50', 'blosum55', 'blosum60', 'blosum62', 'blosum65', 'blosum70', 'blosum75',
                  'blosum80', 'blosum85', 'blosum90', 'blosum95', 'feng', 'fitch', 'genetic', 'gonnet', 'grant', 'ident', 'johnson', 'levin', 'mclach', 'miyata', 'nwsgappep', 'pam120', 'pam180',
                  'pam250', 'pam30', 'pam300', 'pam60', 'pam90', 'rao', 'risler', 'structure']

# protein_models = ['blosum90']
for i in range(len(protein_models)):
    i = 0
    try:
        calculator = DistanceCalculator(protein_models[i])
        # print(calculator)
        # print(type(calculator))
        lst = []
        distance, path = fastdtw(m1_y, m1_g);lst.append(distance);   print('yellow, gray ', str(distance));
        distance, path = fastdtw(m1_y, m1_r);lst.append(distance);   print('yellow, red ', str(distance));
        distance, path = fastdtw(m1_y, m1_b);lst.append(distance);   print('yellow, blue ', str(distance));
        distance, path = fastdtw(m1_g, m1_r);lst.append(distance);   print('gray, red ', str(distance));
        distance, path = fastdtw(m1_g, m1_b);lst.append(distance);   print('gray, blue ', str(distance));
        distance, path = fastdtw(m1_r, m1_b);lst.append(distance);   print('red, blue ', str(distance));
        # print(lst)
        # lost = [[0, lst[0], lst[1], lst[2]],
        #         [lst[0], 0, lst[3], lst[4]],
        #         [lst[1], lst[3], 0, lst[5]],
        #         [lst[2], lst[4], lst[5], 0]]
        # print(np.asarray(lost))
        lst = (lst - min(lst))/(max(lst)-min(lst))
        print(lst)
        lost = [[0,0,0,0],
               [lst[0],0,0,0],
               [lst[1],lst[3],0,0],
               [lst[2],lst[4],lst[5],0]]
        print(np.asarray(lost))

        dm = calculator.get_distance(aln)
        print('\nDistance Matrix\n===================')
        print(dm)
        print(type(dm))
        constructor = DistanceTreeConstructor()
        # tree = constructor.nj(dm)
        tree = constructor.upgma(dm)
        # Phylo.draw(tree)
        # print('\nPhylogenetic Tree\n===================')
        # Phylo.draw_ascii(tree)

        # aln = Bio.AlignIO.read(open('data/test_1.phy'), 'phylip')
        # print(aln)
        starting_tree = tree
        # starting_tree = Bio.Phylo.read('data/nj.tre', 'newick')
        # Bio.Phylo.draw(starting_tree)
        print(type(starting_tree))
        scorer = Bio.Phylo.TreeConstruction.ParsimonyScorer()
        searcher = Bio.Phylo.TreeConstruction.NNITreeSearcher(scorer)
        constructor = Bio.Phylo.TreeConstruction.ParsimonyTreeConstructor(searcher, starting_tree)
        pars_tree = constructor.build_tree(aln)
        print(pars_tree)
        print(protein_models[i])
        Bio.Phylo.draw(pars_tree)
        # print('\nPhylogenetic Tree\n===================')
        # Phylo.draw_ascii(pars_tree) # Print the phylogenetic tree in the terminal
        break
    except:
        print(protein_models[i] + " DID NOT WORK")
        pass

########################################