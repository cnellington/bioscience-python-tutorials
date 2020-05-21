"""
Caleb Ellington
2020.05.22
"""

import numpy as np
import random


def fasta_to_data(fasta_filepath):
    """
    Returns a singular sequence from fasta file

    :param fasta_filepath: input filepath
    :return: uppercase sequence
    """
    with open(fasta_filepath) as sequence_file:
        text = sequence_file.read().upper()
        sequence = (text.split('\n')[1]).strip('\n')
        return sequence

def blosum_to_data(blosum_filepath):
    """
    Returns a scoring matrix from .txt file

    :param blosum_filepath: path to blosum62 matrix (for protein alignment)
    :return: nested dictionary scoring matrix
    """
    with open(blosum_filepath) as blosum_file:
        header = blosum_file.readline().upper().strip().replace('  ', ' ')
        species = header.split(' ')
        score = {}
        for line in blosum_file.readlines():
            row = line.upper().strip().replace('  ', ' ').split(' ')
            species2 = row[0]
            species2_score = {}
            for i in range(len(species)):
                species2_score[species[i]] = int(row[i+1])
            score[species2] = species2_score
        return score


class Alignment:
    """
    Sequence alignment tool
    Requires a scoring matrix
    """

    def __init__(self, scoring_matrix, gap_cost):
        self.scoring_matrix = scoring_matrix
        self.gap_cost = gap_cost

    def get_alignment(self, seq1, seq2):
        # Returns alignment string and score
        score, path_matrix = self._alignment_helper(seq1, seq2)
        align_top = ""
        matches = ""
        align_bottom = ""
        i = len(seq1) - 1
        j = len(seq2) - 1
        while (i, j) != (-1, -1):
            path = path_matrix[i+1][j+1]
            species1 = seq1[i]
            species2 = seq2[j]
            match = " "
            if path == 0:
                i -= 1
                j -= 1
                if species1 == species2:
                    match = "|"
            elif path == 1:
                species1 = "-"
                j -= 1
            else:
                species2 = "-"
                i -= 1
            align_top = species1 + align_top
            align_bottom = species2 + align_bottom
            matches = match + matches
        return f"{align_top}\n{matches}\n{align_bottom}\nscore: {score}"

    def _alignment_helper(self, seq1, seq2):
        # Sets scoring matrix
        alignment_matrix, path_matrix = self.run_needleman_wunsch(seq1, seq2)
        alignment_score = alignment_matrix[-1][-1]
        return alignment_score, path_matrix

    def run_needleman_wunsch(self, seq1, seq2):
        # Returns scoring matrix based on input sequences
        """
        [0, 1]
        [2, X]
        """
        M = len(seq1)
        N = len(seq2)
        alignment_matrix = np.zeros((M+1, N+1))
        path_matrix = np.zeros((M+1, N+1))
        for i in range(1, M+1):
            alignment_matrix[i][0] = self.gap_cost * i
            path_matrix[i][0] = 2
        for j in range(1, N+1):
            alignment_matrix[0][j] = self.gap_cost * j
            path_matrix[0][j] = 1
        for i in range(1, M+1):
            for j in range(1, N+1):
                species1 = seq1[i-1]
                species2 = seq2[j-1]
                diag_score = alignment_matrix[i-1][j-1] + self.scoring_matrix[species1][species2]
                right_score = alignment_matrix[i-1][j] + self.gap_cost
                down_score = alignment_matrix[i][j-1] + self.gap_cost
                score = max(diag_score, max(right_score, down_score))
                alignment_matrix[i][j] = score
                if score == diag_score:
                    path_matrix[i][j] = 0
                elif score == down_score:
                    path_matrix[i][j] = 1
                else:
                    path_matrix[i][j] = 2
        return alignment_matrix, path_matrix
    
    def get_p_value(self, seq1, seq2, n=1e4):
        score, _ = self._alignment_helper(seq1, seq2)
        k = 0
        for i in range(int(n)):
            seq1_permutation = ''.join(random.sample(seq1, len(seq1)))
            score_permutation, _ = self._alignment_helper(seq1_permutation, seq2)
            if score_permutation >= score:
                k += 1
        return (k + 1) / (n + 1)
