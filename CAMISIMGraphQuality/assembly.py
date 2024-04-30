"""Parses assemblies from various different algorithms."""
import os
import numpy as np
import pandas as pd
import parse_seq as parse_seq
import tqdm
from scipy import sparse as sp


class Assembly:
    def __init__(self, assembler, assembly_file):

        self.assembler = assembler
        # get the file extension
        _, ext = os.path.splitext(assembly_file)

        # Needed to construct partial evaluations when some assemblers have not finished.
        if not os.path.exists(assembly_file):
            print(f'Evaluation of {self.assembler} failed. Assembly file does not exist', flush=True)
            self.assembler = None
            return

        # parse the assembly file into an adjacency matrix and the lis of contig
        if ext == '.gfa':
            self.adjM, self.contigs = self.gfa2adjM(assembly_file)
        elif ext == '.fastg':
            self.adjM, self.contigs = self.fastg2adjM(assembly_file)
        elif ext == '.fasta' or ext == '.fa':
            self.adjM, self.contigs = self.fasta2adjM(assembly_file) # return empty nodes X nodes matrix
        else:
            raise Exception(f'{assembly_file} has extension {ext}, which is invalid: .fastg, .fasta/.fa, .gfa, .sif only')

    def __repr__(self):
        return f'Assembly(assembler={self.assembler})'

    def gfa2adjM(self, fl):

        nodes = [e for e in parse_seq.parse_gfa(open(fl)) if e.type == parse_seq.GFATypes.S]
        edges = [e for e in parse_seq.parse_gfa(open(fl)) if e.type == parse_seq.GFATypes.L]

        edges = [(e.l_nid, e.r_nid) for e in edges]
        node_ids = set([e.nid for e in nodes])
        contigs = [(n.nid, n.seq) for n in nodes]

        return self.build_adjM(node_ids, edges), contigs

    def sif2adjM(self, fl, sif_delimiter = ' - '):
        with open(fl) as sif:
            nodes = set()
            edges = list()

            for line in sif:
                u, v = line.strip().split(sif_delimiter) if sif_delimiter in line else (line.strip(), None)
                nodes.add(u)
                if v is not None:
                    edges.append((u,v))
                    nodes.add(v)

        return self.build_adjM(nodes, edges)

    def fasta2adjM(self, fl):

        nodes = [e for e in parse_seq.parse(open(fl), parse_seq.Fasta)]

        node_ids = set([e.hdr for e in nodes])
        # empty list for edges
        edges = []
        contigs = [(n.hdr, n.seq) for n in nodes]

        return self.build_adjM(node_ids, edges), contigs

    def rc(self, s):
        '''
        :param s: sequence string
        :return: reverse complement of s
        '''

        translate = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
        return "".join([translate[e] for e in s])[::-1]

    def fastg2adjM(self, fl):

        nodes = [e for e in parse_seq.parse(open(fl), parse_seq.Fastg)]
        contig = list()
        edge_set = set()
        node_set = set()

        for rec in nodes:

            # add the out node to the contig set.
            # in fastg, each sequence can be present twice as it's reverse complement
            out_node, in_node = rec.out_node, rec.in_node
            if out_node[-1] == "'":
                out_node = out_node[:-1]  # Remove ' prefix
                rec.seq = self.rc(rec.seq)
            if out_node not in node_set:
                contig.append((out_node, rec.seq))
                node_set.add(out_node)

            # if there is an in_node, then the fastg header is recording an edge, so
            # update our edge set.
            if in_node is not None:

                # in fastg, the "in_node" is infact a comma separated list of multiple
                # possible adjacent in nodes, so iterate through one by one
                for adj_node in in_node.split(','):
                    if adj_node[-1] == "'":
                        adj_node = adj_node[:-1]
                    if adj_node[0] == '~':
                        adj_node = adj_node[1:]
                    edge_already_present = (out_node, adj_node) in edge_set or (adj_node, out_node) in edge_set
                    if not edge_already_present:
                        edge_set.add((out_node, adj_node))
        return self.build_adjM(node_set, edge_set), contig

    def build_adjM(self, nodes, edges):

        nodes = sorted(list(set(nodes)))
        node_to_index = {v: i for i, v in enumerate(nodes)}
        num_nodes = len(nodes)
        #adjM = np.full(
        #    (num_nodes, num_nodes), False, dtype=bool
        #)
        adjM = sp.lil_matrix((num_nodes, num_nodes), dtype=bool)
        n_missed = 0
        set_nodes = set(nodes) 
        for u, v in tqdm.tqdm(edges):
            if u not in set_nodes or v not in set_nodes:
                n_missed += 1
                continue
            adjM[node_to_index[u], node_to_index[v]] = True
            adjM[node_to_index[v], node_to_index[u]] = True

        adjM_df = pd.DataFrame.sparse.from_spmatrix(adjM, index=nodes, columns=nodes)
        print(f'missed edges:', n_missed)

        return adjM_df

    #def build_adjM(self, nodes, edges):

    #    nodes = sorted(list(set(nodes)))
    #    node_to_index = {v: i for i, v in enumerate(nodes)}
    #    num_nodes = len(nodes)
    #
    #    sparse_dtype = pd.SparseDtype(bool, fill_value=False)
    #    adjM = pd.DataFrame(False, index=range(num_nodes), columns=range(num_nodes), dtype=sparse_dtype)
    #    n_missed = 0
    #    for u, v in edges:
    #        if u not in nodes or v not in nodes:
    #            n_missed += 1
    #            continue
    #        #adjM.loc[node_to_index[u], node_to_index[v]] = True
    #        #adjM.loc[node_to_index[v], node_to_index[u]] = True
    #        adjM.update(pd.DataFrame([[True]], index=[node_to_index[u]], columns=[node_to_index[v]]))
    #        adjM.update(pd.DataFrame([[True]], index=[node_to_index[v]], columns=[node_to_index[u]]))
    #    adjM.index = nodes
    #    adjM.columns=nodes
    #    print(f'missed edges:', n_missed)

    #    return adjM
