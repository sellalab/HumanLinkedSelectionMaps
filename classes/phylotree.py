__author__ = 'davidmurphy'


import numpy as np


tstring = '(((((hg19:0.0042187,panTro4:0.00882663):4.88763e-05,gorGor3:0.00753179):0.00989961,ponAbe2:0.0161739):0.0084068,nomLeu3:0.0258706):0.0175353,(((rheMac3:0.00376621,macFas5:0.00409417):0.00508014,papHam1:0.00797024):0.00901611,chlSab1:0.0181186):0.0221685);'


def parse_tree(tree_string):
    popen = 0
    pclos = 0
    # this dict holds spec:dist pairs for nodes in the tree
    tree_dict = {}

    # these flags indicate what is being recorded
    rec_spec = False
    rec_dist = False

    # these string variables hold the current name and distance while scanning
    cur_spec = ''
    cur_dist = ''

    # the list keeps track of the current pairs at the current node
    cur_pair = []

    # scan the tree string one character at a time and create tree data struct
    for char in tree_string:
        # scan past outer parentheses
        if char == '(':
            rec_spec = True
            rec_dist = False
            cur_spec = ''
            cur_dist = ''
            cur_pair = []
            popen += 1
        # complete current node
        elif char == ')':
            pclos += 1
            tree_dict[cur_spec] = float(cur_dist)
            # the tree string is finished, no further processing required
            if pclos == popen:
                break
            cur_pair.append(cur_spec)
            label = '{}_{}'.format(*cur_pair)
            cur_pair = []
            cur_spec = label
        # stop recording species, begin recording distance
        elif char == ':':
            rec_spec = False
            rec_dist = True
            cur_dist = ''
        # at node in tree, record current spec:dist pair
        elif char == ',':
            cur_pair.append(cur_spec)
            tree_dict[cur_spec] = float(cur_dist)
            cur_spec = ''
            cur_dist = ''
            rec_dist = False
            rec_spec = True
        # record current char to species string
        elif rec_spec:
            cur_spec += char
        # record current char to distance string
        elif rec_dist:
            cur_dist += char
        elif char == ';':
            break

    return tree_dict


# TODO: create a genuine tree structure while parsing
def make_tree(tree_string):
    popen = 0
    pclos = 0
    # this dict holds spec:dist pairs for nodes in the tree
    tree_dict = {}

    # these flags indicate what is being recorded
    rec_spec = False
    rec_dist = False

    # these string variables hold the current name and distance while scanning
    cur_spec = ''
    cur_dist = ''

    # the list keeps track of the current pairs at the current node
    cur_pair = []

    # scan the tree string one character at a time and create tree data struct
    for char in tree_string:
        # scan past outer parentheses
        if char == '(':
            rec_spec = True
            rec_dist = False
            cur_spec = ''
            cur_dist = ''
            cur_pair = []
            popen += 1
        # complete current node
        elif char == ')':
            pclos += 1
            tree_dict[cur_spec] = float(cur_dist)
            # the tree string is finished, no further processing required
            if pclos == popen:
                break
            cur_pair.append(cur_spec)
            label = '{}_{}'.format(*cur_pair)
            cur_pair = []
            cur_spec = label
        # stop recording species, begin recording distance
        elif char == ':':
            rec_spec = False
            rec_dist = True
            cur_dist = ''
        # at node in tree, record current spec:dist pair
        elif char == ',':
            cur_pair.append(cur_spec)
            tree_dict[cur_spec] = float(cur_dist)
            cur_spec = ''
            cur_dist = ''
            rec_dist = False
            rec_spec = True
        # record current char to species string
        elif rec_spec:
            cur_spec += char
        # record current char to distance string
        elif rec_dist:
            cur_dist += char
        elif char == ';':
            break

    return tree_dict


def parse_exptotsub(xts_file, n_cells=16):
    """parse the results of expected total substitution matrices for MA"""
    # read lines into a single list, cut off the header
    with open(xts_file, 'r') as f:
        lines = f.readlines()[8:]

    # save matrix to dict hashed by node label
    mat_dict = {}

    # create temp variables for each matrix
    cur_mat = []
    cur_lbl = ''

    # extract each matrix and its label from the lines
    for line in lines:
        # split line on whitespace
        line = line.split()

        # skip blank lines
        if len(line) == 0:
            continue

        # get node number and name
        if line[0] == 'Branch':
            cur_lbl = line[-1].split("'")[1]

        # record the matrix
        if len(line) == n_cells + 1:
            cur_mat.append(map(float, line[1:]))

        # when current 16x16 matrix is complete, store in dict
        if len(cur_mat) == n_cells:
            mat_dict[cur_lbl] = np.array(cur_mat)
            cur_mat = []

    return mat_dict


def load_tree(tree_file):
    """get parse branch lengths from tree line in input tree file"""
    with open(tree_file, 'r') as f:
        for line in f:
            if line.startswith('TREE'):
                # return parse_tree(line)
                return line.replace('TREE: ', '')
    raise IOError('{} is not a tree file'.format(tree_file))


class TreeDict:
    def __init__(self, tree_file):
        self.tstring = load_tree(tree_file)
        self.tdict = parse_tree(self.tstring)

    @property
    def branches(self):
        return self.tdict.keys()

    @property
    def lengths(self):
        return self.tdict.values()

    def spec_len(self, sp):
        """get the length of all branches including the species"""
        return sum(self.tdict[b] for b in self.branches if sp in b)

    def group_len(self, grp):
        """get the length of all branches including species in the list"""
        tot = 0
        for sp in grp:
            tot += self.spec_len(sp)
        return tot


class SpeciesPair:
    pass


class Leaf:
    """class to represent leaves on phylogenetic tree"""
    def __init__(self, name, branch, node):
        self.name = name
        self.branch = branch
        self.node = node

    def root_dist(self):
        pass


class Node:
    """class to represent nodes on phylogenetic tree"""
    def __init__(self, branch_1, branch_2, parent=None):
        """intialize node and assign its two branches"""
        self.branch_1 = branch_1
        self.branch_2 = branch_2
        self.parent = parent

    @property
    def root(self):
        return self.parent is None


class PhyloTree:
    """class to represent a phylogenetic tree"""
    def __init__(self, tree_string):
        """read a formatted tree string and create the tree"""
        self.root = None
        self.nodes = []
        self.leaves = []

        # parse the tree string
        for group in tree_string.split(','):
            # remove paranthesis
            group = group.strip('(').strip(')')
            # split on ":" separating branches
            sub_group = group.split(':')
            # split into leaf name, branch length
            if len(sub_group) == 2:
                leaf = sub_group[0]
                blen = sub_group[1]
            # split into leaf name, inner branch, outer branch lengths
            elif len(sub_group) == 3:
                pass