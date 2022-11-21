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
from helper import PairwiseAlignmentHelper as pah

class MultipleAlignmentHelper():
    noGap = 0
    gapA = 1
    gapB = 2
    gapC = 3
    gapAB = 4
    gapBC = 5
    gapAC = 6

    def weightFunctionDifference(self, a, b, c):
        """
            Weight function with 0 if a==b==c, 1 if a==b, a==c or b==c, 2 else.
        """
        if a == b == c:
            return 0
        elif a == b or b == c or a == c:
            return 1
        else:
            return 2

    def createDataForUpgmaWpgma(self, sequences):
        """
            Preprocessing of the sequences for the upgm/wpgm algorithm.
        """
        differenceDictionary = {}
        sequenceToIdMapping = {}
        sequenceToLengthMapping = {}

        for id, sequence in enumerate(sequences):
            sequenceToIdMapping[sequence] = id
            sequenceToLengthMapping[id] = len(sequence)

        for i in range(len(sequences)):
            for j in range(i + 1, len(sequences)):
                score = 0
                for k in range(0, max(len(sequence[i])), len(sequences[j])):
                    a, b = '-', '-'
                    if k < len(sequences[i]):
                        a = sequences[i][k]
                    if k < len(sequences[j]):
                        b = sequences[j][k]
                    score += pah.weightFunctionDifference(a, b)

                differenceDictionary[f'{i} {j}'] = score

        return [
            differenceDictionary,
            sequenceToIdMapping,
            sequenceToLengthMapping
        ]




