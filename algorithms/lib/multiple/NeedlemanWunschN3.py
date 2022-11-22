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
from lib.helper.MultipleAlignmentHelper import MultipleAlignmentHelper as helper

class NeedlemanWunschN3():
    """
        This class computes the Needleman-Wunsch algorithm with three sequences.
    """

    def __init__(self,
                 seqA,
                 seqB,
                 seqC,
                 scoreFunction=helper.nw3DefaultScoreFunction,
                 maxSolutions=-1,
                 echo=False):
        """
            Initalize all variables and methods needed to compute
            the Needleman-Wunsch algorithm with three sequences.
        """
        self.__echo = echo
        self.__score = scoreFunction
        self.seqA = seqA
        self.seqB = seqB
        self.seqC = seqC
        self.matrix = [[[]]]
        self.traceStack = [[]]
        self.traceIndex = 0
        self.traceIndices = [[]]
        self.maxSolutions = maxSolutions

    def __computeMatrix(self, seqA, seqB, seqC):
        """
            Computes the matrix which is needed by the Needleman-Wunsch
            algorithm for three sequences.
        """
        n, m, o = len(seqA) + 1, len(seqB) + 1, len(seqC) + 1
        self.matrix = [[[0 for _ in range(o)] for _ in range(m)] for _ in range(n)]

        # initalize matrix
        for i in range(1, n):
            self.matrix[i][0][0] = self.matrix[i - 1][0][0] + \
                              self.__score('-', '-', seqA[i - 1])

        for i in range(1, m):
            self.matrix[0][i][0] = self.matrix[0][i - 1][0] + \
                              self.__score('-', '-', seqB[i - 1])

        for i in range(1, o):
            self.matrix[0][0][i] = self.matrix[0][0][i - 1] + \
                              self.__score('-', '-', seqC[i - 1])

        for i in range(1, n):
            for j in range(1, m):
                self.matrix[i][j][0] = self.matrix[i - 1][j - 1][0] + \
                                  self.__score(seqA[i - 1], seqB[j - 1], '-')

        for i in range(1, n):
            for k in range(1, o):
                self.matrix[i][0][k] = self.matrix[i - 1][0][k - 1] + \
                                  self.__score(seqA[i - 1], '-', seqC[k - 1])

        for j in range(1, m):
            for k in range(1, o):
                self.matrix[0][j][k] = self.matrix[0][j - 1][k - 1] + \
                                  self.__score('-', seqB[j - 1], seqC[k - 1])

        for i in range(1, n):
            for j in range(1, m):
                for k in range(1, o):
                    self.matrix[i][j][k] = self.__computeMaximum(
                        self.matrix,
                        i, j, k,
                        seqA[i - 1], seqB[j - 1], seqC[k - 1]
                    )

        return self.matrix

    def __computeMaximum(self, matrix, i, j, k, charA, charB, charC):
        """
            Compute the maximum for a given cell of the matrix.
            The maximum is choosen of the following values:
                D(i-1, j-1, k-1) + w(a_i-1, b_j-1, c_k-1)
                D(i, j-1, k-1) + w(a_i, b_j-1, c_k-1)
                D(i-1, j, k-1) + w(a_i-1, b_j, c_k-1)
                D(i-1, j-1, k) + w(a_i-1, b_j-1, c_k)
                D(i, j, k-1) + w(a_i, b_j, c_k-1)
                D(i-1, j, k) + w(a_i-1, b_j, c_k)
                D(i, j-1, k) + w(a_i, b_j-1, c_k)
        """

        # no gap
        noGap = matrix[i - 1][j - 1][k - 1] + self.__score(charA, charB, charC)

        # one gap
        gapA = matrix[i][j - 1][k - 1] + self.__score('-', charB, charC)
        gapB = matrix[i - 1][j][k - 1] + self.__score(charA, '-', charC)
        gapC = matrix[i - 1][j - 1][k] + self.__score(charA, charB, '-')

        # two gaps
        gapAB = matrix[i][j][k - 1] + self.__score('-', '-', charC)
        gapBC = matrix[i - 1][j][k] + self.__score(charA, '-', '-')
        gapAC = matrix[i][j - 1][k] + self.__score('-', charB, '-')

        return max(noGap, gapA, gapB, gapC, gapAB, gapBC, gapAC)

    def __traceback(self, seqA, seqB, seqC, matrix):
        """
            Computes the traceback for the Needleman-Wunsch n=3 matrix.
        """
        self.traceStack = [[]]
        self.traceIndex = 0
        self.traceIndices = [[
            len(matrix) - 1,
            len(matrix[0]) - 1,
            len(matrix[0][0]) - 1,
        ]]

        done = False
        solutionCount = 0

        while not done:
            solutionCount += 1
            i = self.traceIndices[self.traceIndex][0]
            j = self.traceIndices[self.traceIndex][1]
            k = self.traceIndices[self.traceIndex][2]

            appendTraceback = self.traceStack[self.traceIndex].append

            while i > 0 or j > 0 or k > 0:
                traceSplit = False
                pathVariableI = i
                pathVariableJ = j
                pathVariableK = k

                # no gap
                if i > 0 and j > 0 and k > 0:
                    if matrix[i][j][k] == matrix[i - 1][j - 1][k - 1] + \
                            self.__score(seqA[i - 1], seqB[j - 1], seqC[k - 1]):
                        appendTraceback(helper.noGap)
                        pathVariableI -= 1
                        pathVariableJ -= 1
                        pathVariableK -= 1
                        traceSplit = True

                # a gap in sequence a
                if j > 0 and k > 0:
                    if matrix[i][j][k] == matrix[i][j - 1][k - 1] + \
                            self.__score('-', seqB[j - 1], seqC[k - 1]):
                        if not traceSplit:
                            appendTraceback(helper.gapA)
                            pathVariableJ -= 1
                            pathVariableK -= 1
                            traceSplit = True

                        elif [i, j - 1, k - 1] not in self.traceIndices:
                            self.__split([i, j - 1, k - 1], helper.gapA)

                # a gap in sequence b
                if i > 0 and k > 0:
                    if matrix[i][j][k] == matrix[i - 1][j][k - 1] + \
                            self.__score(seqA[i - 1], '-', seqC[k - 1]):
                        if not traceSplit:
                            appendTraceback(helper.gapB)
                            pathVariableI -= 1
                            pathVariableK -= 1
                            traceSplit = True

                        elif [i - 1, j, k - 1] not in self.traceIndices:
                            self.__split([i - 1, j, k - 1], helper.gapB)

                # a gap in sequence c
                if i > 0 and j > 0:
                    if matrix[i][j][k] == matrix[i - 1][j - 1][k] + \
                            + self.__score(seqA[i - 1], seqB[j - 1], '-'):
                        if not traceSplit:
                            appendTraceback(helper.gapC)
                            pathVariableI -= 1
                            pathVariableJ -= 1
                            traceSplit = True

                        elif [i - 1, j - 1, k] not in self.traceIndices:
                            self.__split([i - 1, j - 1, k], helper.gapC)

                # a gap in sequence a and b
                if k > 0:
                    if matrix[i][j][k] == matrix[i][j][k - 1] + \
                            self.__score('-', '-', self.seqC[k - 1]):
                        if not traceSplit:
                            appendTraceback(helper.gapAB)
                            pathVariableK -= 1
                            traceSplit = True

                        elif [i, j, k - 1] not in self.traceIndices:
                            self.__split([i, j, k - 1], helper.gapAB)

                # a gap in sequence a and c
                if j > 0:
                    if matrix[i][j][k] == matrix[i][j - 1][k] + \
                            self.__score('-', seqB[j - 1], '-'):
                        if not traceSplit:
                            appendTraceback(helper.gapAC)
                            pathVariableJ -= 1
                            traceSplit = True

                        elif [i, j - 1, k] not in self.traceIndices:
                            self.__split([i, j - 1, k], helper.gapAC)

                # a gap in sequence b and c
                if i > 0:
                    if matrix[i][j][k] == matrix[i - 1][j][k] + \
                            self.__score(seqA[i - 1], '-', '-'):
                        if not traceSplit:
                            appendTraceback(helper.gapBC)
                            pathVariableI -= 1

                        elif [i - 1, j, k] not in self.traceIndices:
                            self.__split([i - 1, j, k], helper.gapBC)

                i = pathVariableI
                j = pathVariableJ
                k = pathVariableK

            self.traceIndex += 1

            if self.__echo:
                print(f'\rTraces processed: {self.traceIndex:05}/{len(self.traceIndices):05}', end='')

            # done if we have traversed all traces
            if self.traceIndex >= len(self.traceIndices):
                done = True

            # or got exact amount of solutions
            if self.maxSolutions != -1 and solutionCount >= self.maxSolutions:
                done = True

        if self.__echo:
            print()

        resultTraces = None
        if self.maxSolutions == -1:
            resultTraces = self.traceStack
        else:
            resultTraces = self.traceStack[:self.maxSolutions]

        computedAlignment = []
        for trace in resultTraces:
            computedAlignment.append(self.__buildAlignment(trace, seqA, seqB, seqC))

        return computedAlignment

    def __split(self, index, gapSymbol):
        """
            Splits the actual traceback path into two paths.
            index:     The index values for the next cell of the path.
            gapSymbol: A symbol for the computed step for the path.
        """
        self.traceStack.append(self.traceStack[self.traceIndex][0:-1])
        self.traceStack[len(self.traceStack) - 1].append(gapSymbol)
        self.traceIndices.append(index)

    def __buildAlignment(self, traceStack, seqA, seqB, seqC):
        """
            Builds the alignment for one traceback path.
        """
        i, j, k = 0, 0, 0
        alignmentA, alignmentB, alignmentC = '', '', ''

        for direction in traceStack[::-1]:
            if direction == helper.noGap:
                alignmentA += seqA[i]
                alignmentB += seqB[j]
                alignmentC += seqC[k]
                i += 1
                j += 1
                k += 1
            elif direction == helper.gapA:
                alignmentA += '-'
                alignmentB += seqB[j]
                alignmentC += seqC[k]
                j += 1
                k += 1
            elif direction == helper.gapB:
                alignmentA += seqA[i]
                alignmentB += '-'
                alignmentC += seqC[k]
                i += 1
                k += 1
            elif direction == helper.gapC:
                alignmentA += seqA[i]
                alignmentB += seqB[j]
                alignmentC += '-'
                i += 1
                j += 1
            elif direction == helper.gapAB:
                alignmentA += '-'
                alignmentB += '-'
                alignmentC += seqC[k]
                k += 1
            elif direction == helper.gapAC:
                alignmentA += '-'
                alignmentB += seqB[j]
                alignmentC += '-'
                j += 1
            elif direction == helper.gapBC:
                alignmentA += seqA[i]
                alignmentB += '-'
                alignmentC += '-'
                i += 1

        while i < len(seqA):
            alignmentA += seqA[i]
            i += 1

        while j < len(seqB):
            alignmentB += seqB[j]
            j += 1

        while k < len(seqC):
            alignmentC += seqC[k]
            k += 1

        return [alignmentA, alignmentB, alignmentC]

    def compute(self):
        """
            Method to start the computation of the Needleman-Wunsch algorithm
            with three sequences. It returns the computed alignment.
        """

        return self.__traceback(
            self.seqA,
            self.seqB,
            self.seqC,
            self.__computeMatrix(
                self.seqA,
                self.seqB,
                self.seqC
            )
        )
