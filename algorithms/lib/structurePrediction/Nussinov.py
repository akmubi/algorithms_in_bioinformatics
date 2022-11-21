#!/usr/bin/python3
# Copyright 2015 Joachim Wolff
# Programming Course: Algorithms in Bioinformatics
# Tutors: Robert Kleinkauf, Omer Alkhnbashi
# Winter semester 2014/2015
#
# Chair of Bioinformatics
# Department of Computer Science
# Faculty of Engineering
# Albert-Ludwig-University Freiburg im Breisgau
#
# Nussinov algorithm

# TODO: understand & refactor
class Nussinov():
    """
        The algorithm of Nussinov is a RNA secondary structure folding
        algorithm. It was developed by Ruth Nussinov et al. and was published
        in 1978:
        Nussinov, Ruth, et al. "Algorithms for loop matchings."
        SIAM Journal on Applied mathematics 35.1 (1978): 68-82.
        http://rci.rutgers.edu/~piecze/GriggsNussinovKleitmanPieczenik.pdf
    """

    def __init__(self, rnaSequence):
        """
            rnaSequence: The RNA sequence for which the folding should be
                         computed.
        """
        self.sequence = rnaSequence
        self.pairedBases = {}
        self.matrix = [[]]

    def computeMatrix(self):
        """
            This function computes the matrix which the Nussinov-algorithm is
            based on.
        """
        length = len(self.sequence)
        self.matrix = [[0] * (length + 1)] * length

        i = 2
        while i <= len(self.sequence):
            k = i
            j = 0
            while j <= (len(self.sequence)-2) and k <= (len(self.sequence)):
                self.computeMatrixCell(j, k)
                j += 1
                k += 1
            i += 1

    def computeMatrixCell(self, i, j):
        """
            This function computes the value for every cell of the matrix for
            the Nussinov-algorithm.
            i: First index of cell of the Nussinov-matrix
            j: Second index of cell of the Nussinov-matrix
            Every cell is the maximum of:
                            |       N_(i, j-1)
            N_(i,j) = max   |max i <= k < j N_(i, k-1) + N_(k+1, j-1) + 1
                            |       S_k and S_j are complementary
        """
        self.matrix[i][j-1]
        maximumValue = [0,0,0]
        k = i
        while i <= k and k < j:
            if self.areComplementary(self.sequence[k], self.sequence[j-1]):
                pairingValue = self.matrix[i][k-1] + self.matrix[k+1][j-1] + 1
                if maximumValue[2] < pairingValue:
                    maximumValue[0] = k
                    maximumValue[1] = j
                    maximumValue[2] = pairingValue
            k += 1
        self.matrix[i][j] = max(self.matrix[i][j-1], maximumValue[2])

    def areComplementary(self, charA, charB):
        """
            Returns True if two RNA nucleotides are complementary, False
            otherwise. Nucleotides are complemetary if there are 'A' and 'U' or
            'C' and 'G'.
            characterA: First nucleotide
            characterB: Second nucleotide
        """
        if (charA == 'A' and charB == 'U') or (charA == 'U' and charB == 'A'):
            return True
        elif (charA == 'C' and charB == 'G') or (charA == 'G' and charB == 'C'):
            return True
        return False

    def traceback(self, i, j):
        """
            Computes the traceback for the Nussinov-algorithm.
            i: First index of cell of the Nussinov-matrix
            j: Second index of cell of the Nussinov-matrix
        """
        if j <= i:
            return

        elif self.matrix[i][j] == self.matrix[i][j-1]:
            self.traceback(i, j-1)
            return

        else:
            k = i
            while i <= k and k < j:
                if self.areComplementary(self.sequence[k-1], self.sequence[j-1]):

                    if self.matrix[i][j] == self.matrix[i][k-1] + self.matrix[k][j-1] + 1:
                        self.pairedBases[k] = j
                        self.traceback(i, k-1)
                        self.traceback(k, j -1)
                        return
                k += 1

    def execute(self):
        """
            To compute the Nussinov-algorithm execute this method. It returns a
            dictionary with the paired bases.
        """
        self.computeMatrix()
        self.traceback(0, len(self.sequence))
        print(self.pairedBases)
        print(str(len(self.pairedBases)))
        return self.pairedBases
