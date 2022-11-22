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
# from pairwiseAlignmentHelper import PairwiseAlignmentHelper as helper
from lib.helper.PairwiseAlignmentHelper import PairwiseAlignmentHelper as helper

class MultipleAlignmentHelper():
    noGap = 0
    gapA = 1
    gapB = 2
    gapC = 3
    gapAB = 4
    gapBC = 5
    gapAC = 6

    def generateScoreFunction(match, mismatch, gapCost, partialMatch):
        def scoreFunction(a, b, c):
            if a == '-' or b == '-' or c == '-':
                return gapCost
            elif a == b == c:
                return match
            elif a == b or b == c or a == c:
                return partialMatch
            else:
                return mismatch

        return scoreFunction

    @staticmethod
    def nw3DefaultScoreFunction(a, b, c):
        if a == '-' or b == '-' or c == '-':
            return 0
        elif a == b == c:
            return 6
        elif a == b or b == c or a == c:
            return 1
        else:
            return -6

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
                    score += helper.gotohScoreFunction(a, b)

                distances[f'{i} {j}'] = score

        return [
            distances,
            seqToIdMap,
            seqToLengthMap
        ]




