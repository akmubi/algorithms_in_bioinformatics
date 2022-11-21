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
#
from helper import PairwiseAlignmentHelper as helper

class MultipleAlignmentHelper():
    noGap = 0
    gapA = 1
    gapB = 2
    gapC = 3
    gapAB = 4
    gapBC = 5
    gapAC = 6

    @staticmethod
    def weightFunctionDifference(a, b, c):
        """
            Weight function with 0 if a==b==c, 1 if a==b, a==c or b==c, 2 else.
        """
        # if a == b == c:
        #     return 0
        # elif a == b or b == c or a == c:
        #     return 1
        # else:
        #     return 2

        if a == '-' or b == '-' or c == '-':
            if a == b != '-' or b == c != '-' or a == c != '-':
                return 0
            return 0
        elif a == b == c:
            return 2
        elif a == b or b == c or a == c:
            return 1
        else:
            return 0

    def createDataForUpgmaWpgma(self, sequences):
        """
            Preprocessing of the sequences for the upgm/wpgm algorithm.
        """
        distances = {}
        seqToIdMap = {}
        seqToLengthMap = {}

        for id, sequence in enumerate(sequences):
            seqToIdMap[sequence] = id
            seqToLengthMap[id] = len(sequence)

        for i in range(len(sequences)):
            for j in range(i + 1, len(sequences)):
                score = 0
                maxLength = max( len(sequence[i]), len(sequences[j]) )
                for k in range(maxLength):
                    a = sequences[i][k] if k < len(sequences[i]) else '-'
                    b = sequences[j][k] if k < len(sequences[j]) else '-'
                    score += helper.weightFunctionDifference(a, b)

                distances[f'{i} {j}'] = score

        return [
            distances,
            seqToIdMap,
            seqToLengthMap
        ]




