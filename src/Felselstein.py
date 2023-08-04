import networkx as nx
import matplotlib.pyplot as plt
from networkx.drawing.nx_pydot import graphviz_layout
import numpy as np
from itertools import product
from collections import deque
import torch

from transition_matrix_eq_approx import transition_matrix

class Felselstein:

    def __init__(self, graph, sequence_length, nucls):

        # Initial graph
        self.graph = graph

        # Length on aligned sequences
        self.sequence_length = sequence_length

        # Alphabet
        self.nucls = nucls

        # Trained model parameters
        self.v = torch.rand([2, 4], dtype=torch.float64, requires_grad=True)
        self.w = torch.rand([4, 4], dtype=torch.float64, requires_grad=True)
        self.mutation_rate = 1e-6
        self.u = np.full((4, 4), self.mutation_rate)
        for i in range(4):
            self.u[i, i] = 1 - 3 * self.mutation_rate
        self.u = torch.from_numpy(self.u)
        self.Neff = torch.full([], 1000.0, requires_grad=True, dtype=torch.float64)

        # Create transition marix
        self.T_M = transition_matrix(self.v, self.w, self.u, self.Neff)

        # Create likelihood matrices
        for node, data in self.graph.nodes(data=True):
            data['likelihood_matrix'] = torch.zeros((sequence_length, len(nucls)))


    def felselstein_step(self, k_node, a_nucl, position_nucl):

        # If node is leaf then fill probability as 1 or 0
        if len(list(self.graph.successors(k_node))) == 0:
            if self.graph.nodes[k_node]['sequence'][position_nucl] == a_nucl:
                self.graph.nodes[k_node]['likelihood_matrix'][position_nucl, a_nucl] = 1
            
            else:
                self.graph.nodes[k_node]['likelihood_matrix'][position_nucl, a_nucl] = 0

        else:

            # Set likelihood to 0
            self.graph.nodes[k_node]['likelihood_matrix'][position_nucl, a_nucl] = 0

            # Get successor nodes
            l_successor, r_successor = self.graph.successors(k_node)

            # Get successor lengths
            l_time = self.graph.edges[k_node, l_successor]['weight']
            r_time = self.graph.edges[k_node, r_successor]['weight']

            # Sum probabilities of all possible pairs 
            for l_nucl, r_nucl in product(self.nucls, self.nucls):
            
                likelihood = self.T_M[a_nucl, l_nucl] * l_time * \
                             self.graph.nodes[l_successor]['likelihood_matrix'][position_nucl, l_nucl] * \
                             self.T_M[a_nucl, r_nucl] * r_time * \
                             self.graph.nodes[r_successor]['likelihood_matrix'][position_nucl, r_nucl]

                self.graph.nodes[k_node]['likelihood_matrix'][position_nucl, a_nucl] += likelihood

    def run(self):

        # Get head to run BFS
        head = [node for node in self.graph.nodes if len(list(self.graph.predecessors(node))) == 0][0]

        # Run BFS to get node order
        nodes_order = list(nx.bfs_layers(self.graph, head))[::-1]

        # Go through all positions
        for seq_pos in range(self.sequence_length):

            # Go through all nodes in correct order
            for layer in nodes_order:
                for node in layer:

                    # Go through all possible nucleotides
                    for nucl in self.nucls:

                        # Evaluate probability of having 
                        # nucleotide `nucl` in the node `node` in position `seq_pos`
                        self.felselstein_step(node, nucl, seq_pos)


