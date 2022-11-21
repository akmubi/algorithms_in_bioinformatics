#!/usr/bin/python3
# Copyright 2014 Joachim Wolff
# Programming Course: Algorithms in Bioinformatics
# Tutors: Robert Kleinkauf, Omer Alkhnbashi
# Winter semester 2014/2015
#
# Chair of Bioinformatics
# Department of Computer Science
# Faculty of Engineering
# Albert-Ludwig-University Freiburg im Breisgau
class PairwiseAlignmentHelper():
    """
        Class to support the pairwise alignment algorithms
        Needleman-Wunsch and Gotoh.
    """

    # needleman-wunsch
    left     = 0
    up       = 1
    diagonal = 2

    # gotoh
    diagonalD = 0
    dotQ      = 1
    dotP      = 2
    upD       = 3
    upP       = 4
    leftD     = 5
    leftQ     = 6
    matrixIndexD = 0
    matrixIndexP = 1
    matrixIndexQ = 2

    @staticmethod
    def weightFunctionDifference(a ,b):
        if a == b:
            return 0
        else:
            return 1

    @staticmethod
    def weightFunctionDifference2(a, b):
        # gap cost: ('G', -) or (-, 'A') -> 0
        if a == '-' or b == '-':
            return 0
        # match: ('A', 'A') -> 1
        elif a == b:
            return 1
        # mismatch: ('T', 'G') -> 0
        else:
            return 0

    @staticmethod
    def gapCost(x):
        """
            Returns a gap cost of g(x) = 2 + k.
        """
        return 2 + x
