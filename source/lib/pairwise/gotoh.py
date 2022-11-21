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
# Gotoh algorithm
import math
from helper import PairwiseAlignmentHelper as helper

class Gotoh():
    """
        This class holds methods which are needed to compute the pairwise
        alignment algorithm from Osamu Gotoh, published in 1982:
        Osamu Gotoh (1982). "An improved algorithm for matching biological sequences".
        Journal of molecular biology 162: 705.
        https://www.cs.umd.edu/class/spring2003/cmsc838t/papers/gotoh1982.pdf
    """

    def __score(self, a, b):
        return helper.weightFunctionDifference(a, b)

    def __cost(self, x):
        return helper.gapCost(x)

    def __init__(self, seqA, seqB):
        """
            Initalize all variables and methods needed to compute the Gotoh algorithm.
                seqA:      A string with the first DNA sequence.
                seqB:      A string with the second DNA sequence.
        """
        self.seqA = seqA
        self.seqB = seqB
        self.beta = self.__cost(1) - self.__cost(0)
        self.i = 0
        self.j = 0
        self.traceStack = [[]]
        self.currentTrace = 0
        self.traceIndices = [[]]

    def __computeMatrices(self, seqA, seqB, beta):
        """
            Initalize the three matricies needed for the Gotoh-Algorithm.
            The sequences A and B, the weight function and the gap costs have
            to be defined by the creation of the object of this class.
        """
        n, m = len(seqA) + 1, len(seqB) + 1
        matrixD = [[0] * m] * n
        matrixP = [[0] * m] * n
        matrixQ = [[0] * m] * n

        # initalize matrix
        for i in range(1, n):
            matrixD[i][0] = self.__cost(i)
            matrixP[i][0] = math.nan
            matrixQ[i][0] = math.inf
        for i in range(1, m):
            matrixD[0][i] = self.__cost(i)
            matrixP[0][i] = math.inf
            matrixQ[0][i] = math.nan

        for i in range(1, n):
            for j in range(1, m):
                matrixP[i][j] = self.__computeP(
                    matrixD[i - 1][j],
                    matrixP[i - 1][j],
                    beta
                )
                matrixQ[i][j] = self.__computeQ(
                    matrixD[i][j - 1],
                    matrixQ[i][j - 1],
                    beta
                )
                matrixD[i][j] = self.__computeD(
                    matrixD[i - 1][j - 1],
                    matrixP[i][j],
                    matrixQ[i][j],
                    seqA[i - 1],
                    seqB[j - 1]
                )

        return [matrixD, matrixP, matrixQ]

    def __computeP(self, valD, valP, beta):
        """
            Compute the values for matrix P.
                This is the minimum value of:
                    matrix D of cell (i-1, j) + gap costs
                    and
                    matrix P of cell (i-1, j) + 1
                valD: The value from matrix D of cell i-1, j.
                valP: The value from matrix P of cell i-1, j.
                beta: The beta value from the gap costs.
        """
        return min(valD + self.__cost(1), valP + beta)

    def __computeQ(self, valD, valQ, beta):
        """
            Compute the values for matrix Q.
                This is the minimum value of:
                    matrix D of cell (i, j-1) + gap costs
                    and
                    matrix Q of cell (i, j-1) + 1
                valD: The value from matrix D of cell i, j-1.
                valQ: The value from matrix Q of cell i, j-1.
                beta: The beta value from the gap costs.
        """
        return min(valD + self.__cost(1), valQ + beta)

    def __computeD(self, valD, valP, valQ, charA, charB):
        """
            Compute the values for matrix D.
                This is the minimum value of:
                    matrix D of cell (i-1, j-1) + w(a,b)
                    and
                    matrix P of cell (i, j)
                    and
                    matrix Q of cell (i, j)
                valD:  The value from matrix D of cell i-1, j-1.
                valP:  The value from matrix P of cell i, j.
                valQ:  The value from matrix Q of cell i, j.
                charA: The character in sequence A at position i.
                charB: The character in sequence B at position j.
        """
        return min(valP, min(valQ, valD + self.__score(charA, charB)))

    def __traceback(self, seqA, seqB, matrices):
        """
            Computes the traceback for the Gotoh algorithm.
        """
        self.j = len(matrices[0][0]) - 1
        self.i = len(matrices[0]) - 1
        self.currentTrace = 0
        self.traceIndices[self.currentTrace] = [
            self.i,
            self.j,
            helper.matrixIndexD
        ]
        done = False
        while not done:
            currentIndex = self.traceIndices[self.currentTrace]
            while self.i > 0 or self.j > 0:
                if currentIndex[2] == helper.matrixIndexD:
                    self.__tracebackD(matrices)

                elif currentIndex[2] == helper.matrixIndexP:
                    self.__tracebackP(matrices)

                elif currentIndex[2] == helper.matrixIndexQ:
                    self.__tracebackQ(matrices)

                self.i = currentIndex[0]
                self.j = currentIndex[1]
            done = True

            for index in self.traceIndices:
                if index[0] > 0 or index[1] > 0:
                    self.currentTrace = index
                    done = False
                    break
            self.i = currentIndex[0]
            self.j = currentIndex[1]

        computedAlignment = []
        for trace in self.traceStack:
            computedAlignment.append(self.__buildAlignment(trace, seqA, seqB))

    def __tracebackD(self, matrices):
        """
            Computes the traceback for a cell of the matrix D.
        """
        a = self.seqA[self.i - 1]
        b = self.seqB[self.j - 1]
        traceSplit = False
        i = pathVariableI = self.i
        j = pathVariableJ = self.j

        matrixD = matrices[helper.matrixIndexD]
        matrixQ = matrices[helper.matrixIndexQ]
        matrixP = matrices[helper.matrixIndexP]

        if j > 0 and i > 0:
            if matrixD[i][j] == \
               matrixD[i-1][j-1] + self.__score(a,b):
                self.traceStack[self.currentTrace].append(helper.diagonalD)
                pathVariableI -= 1
                pathVariableJ -= 1
                traceSplit = True

            if matrixD[i][j] == matrixQ[i][j]:
                if not traceSplit:
                    self.traceStack[self.currentTrace].append(helper.dotQ)
                    self.traceIndices[self.currentTrace][2] = helper.matrixIndexQ
                    traceSplit = True
                elif [i, j] not in self.traceIndices:
                    self.traceStack.append(self.traceStack[self.currentTrace][0:-1])
                    self.traceStack[len(self.traceStack) - 1].append(helper.dotQ)
                    self.traceIndices.append([i, j, helper.matrixIndexQ])

            if matrixD[i][j] == matrixP[i][j]:
                if not traceSplit:
                    self.traceStack[self.currentTrace].append(helper.dotP)
                    self.traceIndices[self.currentTrace][2] = helper.matrixIndexP
                elif [i, j] not in self.traceIndices:
                    self.traceStack.append(self.traceStack[self.currentTrace][0:-1])
                    self.traceStack[len(self.traceStack)-1].append(helper.dotP)
                    self.traceIndices.append([i, j, helper.matrixIndexP])

        if i == 0:
            self.traceStack[self.currentTrace].append(helper.leftD)
            pathVariableJ -= 1

        if j == 0:
            self.traceStack[self.currentTrace].append(helper.upD)
            pathVariableI -= 1

        if i <= 0 or pathVariableI <= 0:
            pathVariableI = 0

        if j <= 0 or pathVariableJ <= 0:
            pathVariableJ = 0

        self.traceIndices[self.currentTrace][0] = pathVariableI
        self.traceIndices[self.currentTrace][1] = pathVariableJ

    def __tracebackP(self, matrices):
        """
            Computes the traceback for a cell of the matrix P
        """
        traceSplit = False
        matrixP = matrices[helper.matrixIndexP]
        matrixD = matrices[helper.matrixIndexD]

        i, j = self.i, self.j

        if i > 0:
            if matrixP[i][j] == matrixD[i - 1][j] + self.__cost(1):
                self.traceStack[self.currentTrace].append(helper.upD)
                self.traceIndices[self.currentTrace][0] -= 1
                self.traceIndices[self.currentTrace][2] = helper.matrixIndexD
                traceSplit = True

            if matrixP[i][j] == matrixP[i - 1][j] + self.beta:
                if not traceSplit:
                    self.traceStack[self.currentTrace].append(helper.upP)
                    self.traceIndices[self.currentTrace][0] -= 1
                    self.traceIndices[self.currentTrace][2] = helper.matrixIndexP

                elif [i - 1, j] not in self.traceIndices:
                    self.traceStack.append(self.traceStack[self.currentTrace][0:-1])
                    self.traceStack[len(self.traceStack)-1].append(helper.upP)
                    self.traceIndices.append([i - 1, j, helper.matrixIndexP])

    def __tracebackQ(self, matrices):
        """
            Computes the traceback for a cell of the matrix Q
        """
        traceSplit = False
        matrixQ = matrices[helper.matrixIndexQ]
        matrixD = matrices[helper.matrixIndexD]

        i, j = self.i, self.j

        if j > 0:
            if matrixQ[i][j] == matrixD[i][j-1] + self.__cost(1):
                self.traceStack[self.currentTrace].append(helper.leftD)
                self.traceIndices[self.currentTrace][1] -= 1
                self.traceIndices[self.currentTrace][2] = helper.matrixIndexD
                traceSplit = True

            if matrixQ[i][j] == matrixQ[i][j-1] + self.beta:
                if not traceSplit:
                    self.traceStack[self.currentTrace].append(helper.leftQ)
                    self.traceIndices[self.currentTrace][1] -= 1
                    self.traceIndices[self.currentTrace][2] = helper.matrixIndexQ

                elif [i, j - 1] not in self.traceIndices:
                    self.traceStack.append(self.traceStack[self.currentTrace][0:-1])
                    self.traceStack[len(self.traceStack)-1].append(helper.leftQ)
                    self.traceIndices.append([i , j - 1, helper.matrixIndexQ])

    def __buildAlignment(self, traceStack, seqA, seqB):
        """
            A method to compute the alignment of a given traceback of the Gotoh algorithm.
                tracebackStack: The computed traceback path for one alignment as a list.
        """
        i, j = 0, 0
        alignmentA, alignmentB = '', ''

        for direction in traceStack[::-1]:
            if direction == helper.leftQ or direction == helper.leftD:
                alignmentA += '-'
                alignmentB += seqB[j]
                j += 1
            elif direction == helper.upP or direction == helper.upD:
                alignmentA += seqA[i]
                alignmentB += '-'
                i += 1
            elif direction == helper.diagonalD:
                alignmentA += seqA[i]
                alignmentB += seqB[j]
                i += 1
                j += 1

        while i < len(seqA):
            alignmentA += seqA[i]
            i += 1

        while j < len(seqB):
            alignmentB += seqB[j]
            j += 1

        return [alignmentA, alignmentB]

    def compute(self):
        """
            Method to start the computation of the Gotoh algorithm.
        """

        return self.__traceback(
            self.seqA,
            self.seqB,
            self.__computeMatrices(
                self.seqA,
                self.seqB,
                self.beta
            )
        )

