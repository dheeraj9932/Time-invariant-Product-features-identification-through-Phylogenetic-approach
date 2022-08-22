import Bio
from Bio import AlignIO
from Bio import Phylo
from fastdtw import fastdtw
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
from Bio.Phylo.TreeConstruction import *

# DNA_models  = ['identity', 'blastn', 'trans', 'benner6', 'benner22', 'benner74', 'blosum100', 'blosum30', 'blosum35', 'blosum40', 'blosum45', 'blosum50', 'blosum55', 'blosum60', 'blosum62', 'blosum65',
#                'blosum70', 'blosum75', 'blosum80', 'blosum85', 'blosum90', 'blosum95', 'feng', 'fitch', 'genetic', 'gonnet', 'grant', 'ident', 'johnson', 'levin', 'mclach', 'miyata', 'nwsgappep', 'pam120',
#                'pam180', 'pam250', 'pam30', 'pam300', 'pam60', 'pam90', 'rao', 'risler', 'structure']
# protein_models = ['benner6', 'benner22', 'benner74', 'blosum100', 'blosum30', 'blosum35', 'blosum40', 'blosum45', 'blosum50', 'blosum55', 'blosum60', 'blosum62', 'blosum65', 'blosum70', 'blosum75',
#                   'blosum80', 'blosum85', 'blosum90', 'blosum95', 'feng', 'fitch', 'genetic', 'gonnet', 'grant', 'ident', 'johnson', 'levin', 'mclach', 'miyata', 'nwsgappep', 'pam120', 'pam180',
#                   'pam250', 'pam30', 'pam300', 'pam60', 'pam90', 'rao', 'risler', 'structure']
protein_models = ['benner6']
def phylogeny_tree(alignment):
    aln = alignment
    for i in range(len(protein_models)):
        try:
            calculator = DistanceCalculator(protein_models[i])
            dm = calculator.get_distance(aln)
            # print('\nDistance Matrix\n===================')
            # print(dm)
            constructor = DistanceTreeConstructor()
            starting_tree = constructor.nj(dm)
            # starting_tree = constructor.upgma(dm)
            scorer = Bio.Phylo.TreeConstruction.ParsimonyScorer()
            searcher = Bio.Phylo.TreeConstruction.NNITreeSearcher(scorer)
            constructor = Bio.Phylo.TreeConstruction.ParsimonyTreeConstructor(searcher, starting_tree)
            parsimonous_tree = constructor.build_tree(aln)
            # print(parsimonous_tree)
            print(protein_models[i])
            Bio.Phylo.draw(parsimonous_tree)
            # print('\nPhylogenetic Tree\n===================')
            # Phylo.draw_ascii(parsimonous_tree) # Print the phylogenetic tree in the terminal
            # break
        except:
            print(protein_models[i] + " DID NOT WORK")
            pass

########################################
# alignment = AlignIO.read('data/normalised_sequence_swir.phy', 'phylip')
file_name = 'txt_resp_phy_files/normalised_sequence_swir_50.phy'
alignment = AlignIO.read('data/' + file_name, 'phylip')
phylogeny_tree(alignment)