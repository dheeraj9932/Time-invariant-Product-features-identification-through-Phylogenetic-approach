import Bio
from Bio import AlignIO
from Bio import Phylo
import numpy as np
from Bio.Phylo.Consensus import _BitString
from fastdtw import fastdtw
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
from Bio.Phylo.TreeConstruction import *
import networkx
from networkx.drawing import nx_agraph
import pylab
networkx.graphviz_layout = nx_agraph.graphviz_layout

# DNA_models  = ['identity', 'blastn', 'trans', 'benner6', 'benner22', 'benner74', 'blosum100', 'blosum30', 'blosum35', 'blosum40', 'blosum45', 'blosum50', 'blosum55', 'blosum60', 'blosum62', 'blosum65',
#                'blosum70', 'blosum75', 'blosum80', 'blosum85', 'blosum90', 'blosum95', 'feng', 'fitch', 'genetic', 'gonnet', 'grant', 'ident', 'johnson', 'levin', 'mclach', 'miyata', 'nwsgappep', 'pam120',
#                'pam180', 'pam250', 'pam30', 'pam300', 'pam60', 'pam90', 'rao', 'risler', 'structure']
# protein_models = ['benner6', 'benner22', 'benner74', 'blosum100', 'blosum30', 'blosum35', 'blosum40', 'blosum45', 'blosum50', 'blosum55', 'blosum60', 'blosum62', 'blosum65', 'blosum70', 'blosum75',
#                   'blosum80', 'blosum85', 'blosum90', 'blosum95', 'feng', 'fitch', 'genetic', 'gonnet', 'grant', 'ident', 'johnson', 'levin', 'mclach', 'miyata', 'nwsgappep', 'pam120', 'pam180',
#                   'pam250', 'pam30', 'pam300', 'pam60', 'pam90', 'rao', 'risler', 'structure']
# protein_models = ['benner6_10','benner6_20','benner6_26','benner6_30','benner6_40','benner6_52','benner6_53','benner6_62','benner6_70','benner6_80','benner6_93']
protein_models = ['benner6_93']
def _bitstrs(tree):
    bitstrs = set()
    term_names = [term.name for term in tree.get_terminals()]
    term_names.sort()
    for clade in tree.get_nonterminals():
        clade_term_names = [term.name for term in clade.get_terminals()]
        boolvals = [name in clade_term_names for name in term_names]
        bitstr = _BitString("".join(map(str, map(int, boolvals))))
        bitstrs.add(bitstr)
    return bitstrs
def compare(tree1, tree2):
    term_names1 = [term.name for term in tree1.get_terminals()]
    term_names2 = [term.name for term in tree2.get_terminals()]
    # false if terminals are not the same
    if set(term_names1) != set(term_names2):
        return False
    # true if _BitStrings are the same
    if _bitstrs(tree1) == _bitstrs(tree2):
        return True
    else:
        return False

def to_distance_matrix(tree):
    """Create a distance matrix (NumPy array) from clades/branches in tree.
    # https://biopython.org/wiki/Phylo_cookbook
    A cell (i,j) in the array is the length of the branch between allclades[i]
    and allclades[j], if a branch exists, otherwise infinity.

    Returns a tuple of (allclades, distance_matrix) where allclades is a list of
    clades and distance_matrix is a NumPy 2D array.
    """
    allclades = list(tree.find_clades(order="level"))
    lookup = {}
    for i, elem in enumerate(allclades):
        lookup[elem] = i
    distmat = np.repeat(np.inf, len(allclades) ** 2)
    distmat.shape = (len(allclades), len(allclades))
    for parent in tree.find_clades(terminal=False, order="level"):
        for child in parent.clades:
            if child.branch_length:
                distmat[lookup[parent], lookup[child]] = child.branch_length
    if not tree.rooted:
        distmat += distmat.transpose
    np.savetxt("distace_mat_from_tree.csv", distmat, delimiter=",")
    return allclades, np.matrix(distmat)

def terminal_neighbor_dists(self):
    """Return a list of distances between adjacent terminals."""

    def generate_pairs(self):
        pairs = itertools.tee(self)
        pairs[1].__next__()
        return zip(pairs[0], pairs[1])

    return [self.distance(*i) for i in generate_pairs(self.find_clades(terminal = True))]

def phylogeny_tree(alignment):
    aln = alignment
    for i in range(len(protein_models)):
        try:
            calculator = DistanceCalculator(protein_models[i])
            dm = calculator.get_distance(aln)
            print('\nDistance Matrix\n===================')
            print((dm))
            constructor = DistanceTreeConstructor()
            starting_tree = constructor.nj(dm)

            # Phylo.write(starting_tree, "example_tree_format_1.xml", "phyloxml")
            # unrooted_tree = Phylo.read("example_tree_format_1.xml", "phyloxml")
            # Phylo.draw_graphviz(unrooted_tree, node_size=0)
            # pylab.show()
            print(type(starting_tree))
            Phylo.write(starting_tree,'data/unrooted_trees/'+newick_name+'_starting_tree.newick',format='newick')
            # Phylo.write(starting_tree,'trial.xml',format='phyloxml')
            scorer = Bio.Phylo.TreeConstruction.ParsimonyScorer()
            searcher = Bio.Phylo.TreeConstruction.NNITreeSearcher(scorer)
            constructor = Bio.Phylo.TreeConstruction.ParsimonyTreeConstructor(searcher, starting_tree)
            parsimonous_tree = constructor.build_tree(aln)
            Phylo.write(parsimonous_tree,'data/unrooted_trees/'+newick_name+'_parsimonous_tree.newick', format='newick')
            print(protein_models[i], ' respective score for tree is ', scorer.get_score(parsimonous_tree, aln))
            # pars_tree, distance_mat = to_distance_matrix(parsimonous_tree)
            # Phylo.write(parsimonous_tree, "example_tree_format.xml", "phyloxml")
            # Bio.Phylo.draw(parsimonous_tree)
            print(type(parsimonous_tree))

            # Bio.Phylo.draw_graphviz(parsimonous_tree)
            # unrooted_tree = Phylo.read("example_tree_format.xml", "phyloxml")
            # Phylo.draw_graphviz(unrooted_tree, node_size = 0)
            # pylab.show()


        except:
            print(protein_models[i] + " DID NOT WORK")
            continue
    return parsimonous_tree

########################################0
newick_name = "segment5of5"
file_name = 'txt_resp_phy_files/2022_05_segment_wise_trees/5_segments/local_normalised_sequence_segment_normalised_swir_5_of_5_segments.phy'; alignment = AlignIO.read('data/' + file_name, 'phylip'); parsimonous_tree9 = phylogeny_tree(alignment)
print(terminal_neighbor_dists(parsimonous_tree9))



# file_name = 'txt_resp_phy_files/global_normalised_sequence_swir_80.phy'; alignment = AlignIO.read('data/' + file_name, 'phylip'); parsimonous_tree8 = phylogeny_tree(alignment)
# file_name = 'txt_resp_phy_files/global_normalised_sequence_swir_70.phy'; alignment = AlignIO.read('data/' + file_name, 'phylip'); parsimonous_tree7 = phylogeny_tree(alignment)
# file_name = 'txt_resp_phy_files/global_normalised_sequence_swir_62.phy'; alignment = AlignIO.read('data/' + file_name, 'phylip'); parsimonous_tree6 = phylogeny_tree(alignment)
# file_name = 'txt_resp_phy_files/global_normalised_sequence_swir_53.phy'; alignment = AlignIO.read('data/' + file_name, 'phylip'); parsimonous_tree5 = phylogeny_tree(alignment)
# file_name = 'txt_resp_phy_files/global_normalised_sequence_swir_40.phy'; alignment = AlignIO.read('data/' + file_name, 'phylip'); parsimonous_tree4 = phylogeny_tree(alignment)
# file_name = 'txt_resp_phy_files/global_normalised_sequence_swir_30.phy'; alignment = AlignIO.read('data/' + file_name, 'phylip'); parsimonous_tree3 = phylogeny_tree(alignment)
# file_name = 'txt_resp_phy_files/global_normalised_sequence_swir_20.phy'; alignment = AlignIO.read('data/' + file_name, 'phylip'); parsimonous_tree2 = phylogeny_tree(alignment)
# file_name = 'txt_resp_phy_files/global_normalised_sequence_swir_10.phy'; alignment = AlignIO.read('data/' + file_name, 'phylip'); parsimonous_tree1 = phylogeny_tree(alignment)
# print('Are the trees similar in structure?  '+str(compare(parsimonous_tree1, parsimonous_tree2)))
# print('#############################################################################################################')
# print('\n')
# file_name = 'txt_resp_phy_files/global_normalised_sequence_swir_52.phy'
# alignment = AlignIO.read('data/' + file_name, 'phylip')
# parsimonous_tree3 = phylogeny_tree(alignment)
# print('#############################################################################################################')
# print('Are the trees similar in structure?  '+str(compare(parsimonous_tree1, parsimonous_tree3)))
# print('#############################################################################################################')
# print('\n')
# file_name = 'txt_resp_phy_files/global_normalised_sequence_swir_80.phy'
# alignment = AlignIO.read('data/' + file_name, 'phylip')
# parsimonous_tree4 = phylogeny_tree(alignment)
# print('#############################################################################################################')
# print('Are the trees similar in structure?  '+str(compare(parsimonous_tree1, parsimonous_tree1)))
# print('#############################################################################################################')
