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
# Needleman-Wunsch algorithm
import sys
from helper import PairwiseAlignmentHelper as helper

class NeedlemanWunsch():
    """
        This class holds methods which are needed to compute the pairwise
        alignment algorithm from Saul Needleman and Christian Wunsch, published in 1970:
        Needleman, Saul B.; and Wunsch, Christian D. (1070).
        A general method applicable to search for similarities in the aminoacid
        sequence of two proteins. Journal of Molecular Biology 48 (3): 443-53
        http://www.cise.ufl.edu/class/cis4930sp09rab/00052.pdf
    """

    def __init__(self, seqA, seqB, maxSolutions=-1):
        """
            Initalize all variables and methods needed to compute
            the Needleman-Wunsch algorithm.
            seqA: A string with the first DNA sequence.
            seqB: A string with the second DNA sequence.
        """
        self.seqA = seqA
        self.seqB = seqB
        self.maxSolutions = maxSolutions

    def __score(self, a, b):
        return helper.weightFunctionDifference2(a, b)

    def __computeMatrix(self, seqA, seqB):
        """
            Initalize and computes the values for the Needleman-Wunsch matrix.
                seqA: A string with the first DNA sequence.
                seqB: A string with the second DNA sequence.
        """
        matrix = [[0 for i in range(len(seqB) + 1)] for j in range(len(seqA) + 1)]

        # initalize matrix
        for i in range(1, len(seqA) + 1):
            matrix[i][0] = matrix[i - 1][0] + self.__score(seqA[i - 1], '-')
        for i in range(1, len(seqB) + 1):
            matrix[0][i] = matrix[0][i - 1] + self.__score('-', seqB[i - 1])

        for i in range(1, len(seqA) + 1):
            for j in range(1, len(seqB) + 1):
                matrix[i][j] = self.__computeMaximum(
                    seqA[i - 1],
                    seqB[j - 1],
                    matrix[i][j - 1],
                    matrix[i - 1][j],
                    matrix[i - 1][j - 1]
                )
        return matrix

    def __computeMaximum(self, charA, charB, predLeft, predUp, predDiagonal):
        """
            Computes the maximum of a given cell for the Needleman-Wunsch matrix.
                charA:        The character in sequence A at position i.
                charB:        The character in sequence B at position j.
                predLeft:     The value i, j-1 in the matrix.
                predUp:       The value i-1, j in the matrix.
                predDiagonal: The value i-1, j-1 in the matrix.
        """

        costUp = predUp + self.__score(charA, '-')
        costDiagonal = predDiagonal + self.__score(charA, charB)
        costLeft = predLeft + self.__score('-', charB)
        return max(costUp, costDiagonal, costLeft)

    def __traceback(self, seqA, seqB, matrix):
        """
            Computes the traceback for the Needleman-Wunsch matrix.
                seqA:   A string with the first DNA sequence.
                seqB:   A string with the second DNA sequence.
                matrix: The computed matrix for the two sequences.
        """
        m = len(matrix)
        n = len(matrix[0])

        traces = [[]]
        traceIndices = [[m - 1, n - 1]]
        currentTrace = 0
        done = False
        solutionCount = 0
        l = 0 # current trace 2
        # traversedCount = 0
        appendTrace = traces.append
        appendIndex = traceIndices.append

        while not done:
            solutionCount += 1
            i = traceIndices[currentTrace][0]
            j = traceIndices[currentTrace][1]

            appendTraceback = traces[currentTrace].append
            # start traceback from cell [m-1][n-1]
            while i > 0 or j > 0:
                pathVariableI = i
                pathVariableJ = j
                # if we encounter two possible paths,
                # we need to split current traceback into two branches
                splitTraceback = False

                # if min cost is from left cell
                if j > 0 and matrix[i][j] == matrix[i][j - 1] + self.__score('-', seqB[j-1]):
                    # add left arrow
                    appendTraceback(helper.left)
                    # go into left cell
                    pathVariableJ -= 1
                    splitTraceback = True

                # if min cost is from upper cell
                if i > 0 and matrix[i][j] == matrix[i - 1][j] + self.__score(seqA[i-1], '-'):
                    # add trace only if it doesn't exists yet
                    if not splitTraceback:
                        # add up arrow
                        appendTraceback(helper.up)
                        # go into upper cell
                        pathVariableI -= 1
                        splitTraceback = True
                    elif [i - 1, j] not in traceIndices:
                        # create alternative traceback (branch)
                        # example:
                        # original sequence: [[up, diagonal, left]]
                        # add alternative:
                        # -> [[up, diagonal, left], [up, diagonal]]
                        appendTrace(traces[currentTrace][0:-1])
                        # -> [[up, diagonal, left], [up, diagonal, up]]
                        traces[len(traces) - 1].append(helper.up)
                        # save start coords of alternative trace
                        appendIndex([i - 1, j])

                # if min cost if from diagonal cell
                if (i > 0 and j > 0) and matrix[i][j] == matrix[i - 1][j - 1] + self.__score(seqA[i - 1], seqB[j - 1]):
                    # add trace only if it doesn't exists yet
                    if not splitTraceback:
                        # add diagonal arrow
                        appendTraceback(helper.diagonal)
                        # go into diagonal cell
                        pathVariableI -= 1
                        pathVariableJ -= 1
                    elif [i - 1, j - 1] not in traceIndices:
                        # create alternative traceback (branch)
                        appendTrace(traces[currentTrace][0:-1])
                        traces[len(traces) - 1].append(helper.diagonal)
                        # save start coords of alternative trace
                        appendIndex([i - 1, j - 1])
                i = pathVariableI
                j = pathVariableJ

            currentTrace += 1

            print(f'\rComputing traces: {currentTrace:05}/{len(traceIndices):05}', end='')

            # done if we have traversed all traces
            if currentTrace >= len(traceIndices):
                done = True
            # or got exact amount of solutions
            if self.maxSolutions != -1 and solutionCount >= self.maxSolutions:
                done = True

        print()

        computedAlignment = []
        resultTraces = None
        if self.maxSolutions == -1:
            resultTraces = traces
        else:
            resultTraces = traces[:self.maxSolutions]

        for trace in resultTraces:
            computedAlignment.append(self.__buildAlignment(trace, seqA, seqB))
        return computedAlignment

    def __buildAlignment(self, traceStack, seqA, seqB):
        """
            Builds the alignment for one traceback path.
                traceStack: The computed tracebackpath as a list = []
                seqA:       A string with the first DNA sequence.
                seqB:       A string with the second DNA sequence.
        """
        i, j = 0, 0
        alignmentA, alignmentB = '', ''

        # iterate backwards
        for direction in traceStack[::-1]:
            if direction == helper.left:
                alignmentA += '-'
                alignmentB += seqB[j]
                j += 1
            elif direction == helper.up:
                alignmentA += seqA[i]
                alignmentB += '-'
                i += 1
            elif direction == helper.diagonal:
                alignmentA += seqA[i]
                alignmentB += seqB[j]
                i += 1
                j += 1
            else:
                print(f'unknown element encountered in trace stack: {direction}')
                sys.exit()

        # append rest of sequences
        while i < len(seqA):
            alignmentA += seqA[i]
            i += 1

        while j < len(seqB):
            alignmentB += seqB[j]
            j += 1

        return [alignmentA, alignmentB]

    def compute(self):
        """
            Method to execute the Needleman-Wunsch algorithm.
                sequences:      A list with two strings which represents the DNA sequences.
        """

        return self.__traceback(
                    self.seqA,
                    self.seqB,
                    self.__computeMatrix(self.seqA, self.seqB)
                )
