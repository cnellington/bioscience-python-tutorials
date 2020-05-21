"""
Caleb Ellington
2020/05/20
University of Washington
Department of Bioengineering
"""

import numpy as np
import random


class Alignment:
    """
    Runs various alignment types from 2 input text files. Amino acids are considered
    by defaul unless another scoring matrix is given.
    """

    def __init__(self, file1, file2, scorefile="scoring_matrices/blosum62.txt", gap_cost=-4, verbose=False):
        """
        Initializes the alignment model

        :param file1: first sequence file
        :param file2: second sequence file
        :param scorefile: scoring matrix (BLOSUM62 by default)
        :param gap_cost: cost of gaps in the alignment
        :param verbose: run with flags for alignment completion and updates
        """
        self._file1 = file1
        self._file2 = file2
        self._scorefile = scorefile
        self._verbose = verbose
        self.gap_cost = gap_cost
        self.score_matrix = None
        self.path_matrix = None
        self._load()

    def _load(self):
        with open(self._file1, 'r') as sequence_file:
            self.seq1 = sequence_file.read().upper()
        with open(self._file2, 'r') as sequence_file:
            self.seq2 = sequence_file.read().upper()
        with open(self._scorefile, 'r') as score_file:
            header = score_file.readline().upper().strip().replace('  ', ' ')
            species = header.split(' ')
            score = {}
            for line in score_file.readlines():
                row = line.upper().strip().replace('  ', ' ').split(' ')
                species2 = row[0]
                species2_score = {}
                for i in range(len(species)):
                    species2_score[species[i]] = int(row[i+1])
                score[species2] = species2_score
            self.score = score
        if self._verbose:
            print('Finished loading successfully')

    def get_alignment(self):
        """
        Returns the alignment string from two sequence inputs.
        Requires that the needleman-wunsch algorithm has been run.

        :return: Alignment of self.seq1 and self.seq2
        """
        align1 = ""
        matches = ""
        align2 = ""
        # Start at the bottom right corner
        i = len(self.seq1) - 1
        j = len(self.seq2) - 1
        # Work our way back until we hit the top left corner
        while (i, j) != (-1, -1):
            path = self.path_matrix[i+1][j+1]
            species1 = self.seq1[i]
            species2 = self.seq2[j]
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
            align1 = species1 + align1
            align2 = species2 + align2
            matches = match + matches
        return f"{align1}\n{matches}\n{align2}"

    def get_local_alignment(self):
        """
        Returns the best local alignment after running smith-waterman.
        Requires that the smith-waterman algorithm has been run.

        :return: Local alignment from self.seq1 and self.seq2
        """
        align1 = ""
        matches = ""
        align2 = ""
        # Find max score indices
        max_coords = np.where(self.score_matrix == np.amax(self.score_matrix))
        (i, j) = list(zip(max_coords[0], max_coords[1]))[0]
        # Traceback until a 0 is encountered, build alignment as we go
        while self.score_matrix[i][j] != 0:
            path = self.path_matrix[i][j]
            species1 = self.seq1[i-1]
            species2 = self.seq2[j-1]
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
            align1 = species1 + align1
            align2 = species2 + align2
            matches = match + matches
        return f"Alignment at seq1 index {i}, seq2 index {j}\n{align1}\n{matches}\n{align2}"

    def get_nw_score(self):
        """
        Gets the score from a needleman-wunsch alignment.
        :return: Needleman-wunsch global alignment score
        """
        return self.score_matrix[-1][-1]

    def get_sw_score(self):
        """
        Gets the score from a smith-waterman alignment
        :return: Smith-waterman local alignment score
        """
        return np.amax(self.score_matrix)

    def get_score_matrix(self):
        """
        Gets the alignment scoring matrix
        :return: alignment scoring matrix of floats
        """
        return self.score_matrix

    def get_path_matrix(self):
        """
        Gets the alignment path matrix
        :return: alignment path matrix, directed from the bottom right
                 upward, considering each point X as [[0, 2], [1, X]]
                 for backtracing left and upward from the bottom right.
        """
        return self.path_matrix

    def get_p_value(self, n=1e4):
        """
        Gets the p-value for a needleman-wunsch global alignment
        :param n: permutations to test (affects p-value fidelity)
        :return: p-value
        """
        k = 0.0
        seq1_score = self.get_nw_score()
        seq1 = self.seq1
        for i in range(int(n)):
            self.seq1 = ''.join(random.sample(self.seq1, len(self.seq1)))
            self.run_needleman_wunsch()
            k_score = self.get_nw_score()
            if k_score >= seq1_score:
                k += 1
        return (k+1)/(n+1)


    def run_needleman_wunsch(self):
        """
        Runs the needleman-wunsch alignment algorithm on self._seq1 and self._seq2.
        After completion, alignment metrics are accessible by other class methods.
        """
        if self._verbose:
            print(f"seq1: {self.seq1}\nseq2: {self.seq2}")
        M = len(self.seq1)
        N = len(self.seq2)
        self.score_matrix = np.zeros((M+1, N+1))
        self.path_matrix = np.zeros((M+1, N+1))
        for i in range(1, M+1):
            self.score_matrix[i][0] = self.gap_cost * i
            self.path_matrix[i][0] = 2
        for j in range(1, N+1):
            self.score_matrix[0][j] = self.gap_cost * j
            self.path_matrix[0][j] = 1
        for i in range(1, M+1):
            for j in range(1, N+1):
                species1 = self.seq1[i-1]
                species2 = self.seq2[j-1]
                diag_score = self.score_matrix[i-1][j-1] + self.score[species1][species2]
                right_score = self.score_matrix[i-1][j] + self.gap_cost
                down_score = self.score_matrix[i][j-1] + self.gap_cost
                score = max(diag_score, max(right_score, down_score))
                self.score_matrix[i][j] = score
                # Catalog path using a 2x2 array representation [[0,2],[1,3]]
                if score == diag_score:
                    self.path_matrix[i][j] = 0
                elif score == down_score:
                    self.path_matrix[i][j] = 1
                else:
                    self.path_matrix[i][j] = 2
        if self._verbose:
            print("Finished running needleman-wunsch alignment")

    def run_smith_waterman(self):
        """
        Runs the smith-waterman alignment algorithm on self.seq1 and self.seq2.
        After completion, alignment metrics are accessible by other class methods.
        """
        if self._verbose:
            print(f"seq1: {self.seq1}\nseq2: {self.seq2}")
        M = len(self.seq1)
        N = len(self.seq2)
        self.score_matrix = np.zeros((M + 1, N + 1))
        self.path_matrix = np.zeros((M + 1, N + 1))
        for i in range(1, M + 1):
            for j in range(1, N + 1):
                species1 = self.seq1[i - 1]
                species2 = self.seq2[j - 1]
                diag_score = self.score_matrix[i - 1][j - 1] + self.score[species1][species2]
                right_score = self.score_matrix[i - 1][j] + self.gap_cost
                down_score = self.score_matrix[i][j - 1] + self.gap_cost
                score = max(0, max(diag_score, max(right_score, down_score)))
                self.score_matrix[i][j] = score
                # Catalog path using a 2x2 array representation [[0,2],[1,3]]
                if score == diag_score:
                    self.path_matrix[i][j] = 0
                elif score == down_score:
                    self.path_matrix[i][j] = 1
                else:
                    self.path_matrix[i][j] = 2
        if self._verbose:
            print("Finished running smith-waterman alignment")


