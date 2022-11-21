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
# Sum of pairs algorithm
from lib.helper.PairwiseAlignmentHelper import PairwiseAlignmentHelper as helper

class SumOfPairs():
    """
        This class computes the Sum-of-pairs algorithm by Carillo and Lipman:
        Carrillo, Humberto, and David Lipman.
        "The multiple sequence alignment problem in biology."
        SIAM Journal on Applied Mathematics 48.5 (1988): 1073-1082.
        http://www.academia.edu/download/30855770/Articulo03.pdf
    """

    def __similarity(self, a, b):
        return helper.pam250(a, b)

    def __init__(self, sequences):
        """
            To initialize a object of the SumOfPairs class please define a list
            with the multiple sequence alignment and a similarity score method
            which is defined in class PairwiseAlignmentHelper of package helper.
        """
        self.sequences = sequences

    def execute(self):
        """
            Run this method to compute the sum of pairs scoring for multiple
            alignment.
        """
        score_value = 0
        for i in range(0, len(self.sequences)):
            for j in range(i + 1, len(self.sequences)):
                score_value += self.__score(self.sequences[i], self.sequences[j])
        return score_value

    def __score(self, seqA, seqB):
        """
            Returns the pairwise alignment for seqA and seqB.
        """
        scoreValue = 0
        maxLen = max(len(seqA), len(seqB))
        for i in range(maxLen):
            if i < len(seqA) and i < len(seqB):
                scoreValue += self.__similarity(seqA[i], seqB[i])
            elif i < len(seqA):
                scoreValue += self.__similarity(seqA[i], '-')
            elif i < len(seqB):
                scoreValue += self.__similarity('-', seqB[i])

        return scoreValue
