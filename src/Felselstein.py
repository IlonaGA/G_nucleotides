import networkx as nx
import matplotlib.pyplot as plt
from networkx.drawing.nx_pydot import graphviz_layout
import pydot
import numpy as np
from itertools import product
from collections import deque

class Felselstein:

    def __init__(self, graph, sequence_length, nucls):

        # Initial graph
        self.graph = graph

        # Length on aligned sequences
        self.sequence_length = sequence_length

        # Alphabet
        self.nucls = nucls
        self.nucl_to_ind_dict = {nucl: i for i, nucl in enumerate(nucls)}

        # Trained model parameters
        self.v = ...
        self.w = ...


    def Felselstein_step(self, graph, k_node, a_nucl, position_nucl):
        
        # Get index of predesscor nucleotide
        a_nucl_ind = self.nucl_to_ind_dict[a_nucl]
        
        # If node is leaf then fill probability as 1 or 0
        if len(list(self.graph.successors(k_node))) == 0:
            if graph.nodes[k_node]['sequence'][position_nucl] == a_nucl:
                graph.nodes[k_node]['likelihood_matrix'][position_nucl, a_nucl_ind] = 1
            else:
                graph.nodes[k_node]['likelihood_matrix'][position_nucl, a_nucl_ind] = 0

        else:

            # Set likelihood to 0 
            graph.nodes[k_node]['likelihood_matrix'][position_nucl, a_nucl_ind] = 0
            
            # Get successor nodes
            l_successor, r_successor = graph.successors(k_node)
            
            # Get successor lengths
            l_time = graph.edges[k_node, l_successor]['weight']
            r_time = graph.edges[k_node, r_successor]['weight']

            # Sum probabilities of all possible pairs 
            for l_nucl, r_nucl in product(self.nucls, self.nucls):
                l_nucl_ind = self.nucl_to_ind_dict[l_nucl]
                r_nucl_ind = self.nucl_to_ind_dict[r_nucl]

                # Somehow use `self.transition_probability`
                likelihood = ...

                graph.nodes[k_node]['likelihood_matrix'][position_nucl, a_nucl_ind] += likelihood

    def Felselstein(self, graph):
        
        # Get head to run BFS
        head = [node for node in graph.nodes if len(list(graph.predecessors(node))) == 0][0]

        # Run BFS to get node order
        nodes_order = list(nx.bfs_layers(graph, head))[::-1]

        # Go through all positions
        for seq_pos in range(self.sequence_length):
            
            # Go through all nodes in correct order
            for layer in nodes_order:
                for node in layer:

                    # Go through all possible nucleotides
                    for nucl in self.nucls:

                        # Evaluate probability of having 
                        # nucleotide `nucl` in the node `node` in position `seq_pos`
                        self.Felselstein_step(graph, node, nucl, seq_pos)

    def g(self, x):
        ''' Formula 15 '''
        return 2 * x / (1 - np.exp(-2 * x))

    def e(self, seq):
        ''' Formula 30 '''
        v_w_sum = 0
        for i in range(len(seq)):
            nucl_i_ind = self.nucl_to_ind_dict[seq[i]]
            nucl_j_ind = self.nucl_to_ind_dict[seq[j]]

            v_w_sum += self.v[i, nucl_i_ind]

            for j in range(i + 1, len(seq)):
                v_w_sum += self.w[i, j, nucl_i_ind, nucl_j_ind]

        return v_w_sum


    def rate(self, nucl_1, nucl_2):
        ...

    def p_rate(self, seq, ind, nucl):
        ''' Formula 29 '''
        seq_a = seq[:ind] + nucl + seq[ind + 1:]
        return self.rate(seq[ind], nucl) * self.g((self.e(seq_a) - self.e(seq)) / 2)


    def transition_probability(self, seq_successor, seq_predecessor, length):
        ''' Formula 28 '''
        log_left_product = 0
        log_right_product = 0

        for i, (nucl_successor, nucl_predecessor) in enumerate(zip(seq_successor, seq_predecessor)):
            if nucl_successor != nucl_predecessor:
                log_left_product += np.log(self.p_rate(seq_predecessor, i, nucl_successor)) 
                log_left_product += np.log(length)

            else:
                summond = 0
                
                for nucl in self.nucls:
                    if nucl == nucl_predecessor:
                        continue

                    # TODO: Rewrite as logs
                    summond += self.p_rate(seq_predecessor, i, nucl) * length

                log_right_product += np.log(1 - summond)



