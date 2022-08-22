import PyQt5
from ete3 import Tree, TreeStyle, NodeStyle, faces, AttrFace, CircleFace
def layout(node):
    if node.is_leaf():
        N = AttrFace("name", fsize=2)
        faces.add_face_to_node(N, node, 0, position="aligned")
def cluster_check(get_leaves):
    k = []
    for i in range(len(get_leaves)):
        material = (str(get_leaves[i])[3:5])
        k.append(material)
        if len(set(k))> 1:
            return False
        else:
            continue
    return True
def get_example_tree():
    # Set dashed blue lines in all leaves
    nst1 = NodeStyle()
    nst1["bgcolor"] = "DarkGray"
    nst1["fgcolor"] = "OrangeRed"
    nst1["size"] = 10
    # nst1["shape"] = 'square'
    nst2 = NodeStyle()
    nst2["fgcolor"] = "Magenta"
    nst2["size"] = 10
    nst2["bgcolor"] = "Moccasin"

    nst3 = NodeStyle()
    nst3["fgcolor"] = "LimeGreen"
    nst3["size"] = 10
    nst3["bgcolor"] = "DarkSeaGreen"

    nst4 = NodeStyle()
    nst4["fgcolor"] = "MediumBlue"
    nst4["size"] = 10
    nst4["bgcolor"] = "Khaki"

    t = Tree(filename+".newick", format=1)

    for n in t.traverse():
        n.dist = 0

    n1 = t.get_common_ancestor("m1_gelb", "m1_grau", "m1_rot", "m1_blau")
    n1.set_style(nst1)
    n2 = t.get_common_ancestor("m2_gelb", "m2_grau", "m2_rot", "m2_blau")
    n2.set_style(nst2)
    n3 = t.get_common_ancestor("m3_gelb", "m3_grau", "m3_rot", "m3_blau")
    n3.set_style(nst3)
    n4 = t.get_common_ancestor("m4_gelb", "m4_grau", "m4_rot", "m4_blau")
    n4.set_style(nst4)

    ts = TreeStyle()
    ts.layout_fn = layout
    ts.show_leaf_name = False
    ts.show_branch_support  = False

    ts.mode = "c"
    ts.root_opening_factor = 0 #default = 1
    good_cluster_bool = False
    if cluster_check(n1.get_leaves()) and cluster_check(n2.get_leaves()) and cluster_check(n3.get_leaves()) and cluster_check(n4.get_leaves()):
        good_cluster_bool = True
    return t, ts, good_cluster_bool
if __name__ == "__main__":
    filename = "data/unrooted_trees/segment1of4_parsimonous_tree"
    t, ts, good_cluster_bool = get_example_tree()
    t.show(tree_style=ts)
    print(good_cluster_bool)
