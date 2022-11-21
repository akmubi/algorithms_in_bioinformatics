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
import math
from pairwise import NeedlemanWunsch
from helper import PairwiseAlignmentHelper as helper
from multiple import UpgmaWpgma
import random

# TODO: understand
class FengDoolittle():
    """
        This class computes the Feng-Doolittle algorithm by Da-Fei Feng and Russell F. Doolittle:
        Feng, Da-Fei, and Russell F. Doolittle.
        "Progressive sequence alignment as a prerequisitetto correct phylogenetic trees."
        Journal of molecular evolution 25.4 (1987): 351-360.
        http://dna.bio.puc.cl/cardex/papersbio252/Grupo06-2013.pdf
    """

    def score(self, a, b):
        return helper.pam250(a, b)

    def __init__(self, sequences):
        """
            To initialize an object of class FengDoolittle you have to define.
        """
        self.sequences = sequences
        self.alignments = []
        self.alignmentToIndexMapping = {}
        self.sequenceToIndexMapping = {}
        self.distanceDict = {}
        self.newickTree = ''
        self.orderToAlign = []

    def computeAlignments(self):
        """
            This function computes all pairwise alignments between every
            sequence with the Needleman-Wunsch algorithm.
        """
        nw = NeedlemanWunsch(self.sequences[i], self.sequences[j], 1)
        alignmentsAppend = self.alignments.append
        for i in range(0, len(self.sequences)):
            for j in range(i + 1, len(self.sequences)):
                alignmentsAppend([nw.compute(), i, j])

    def __computeDistanceDict(self):
        """
            This function computes the distance between every alignment.
            The distances are used to generate a phylogenetic tree.
        """
        for alignment in self.alignments:
            key = f'{alignment[1]} {alignment[2]}'
            self.distanceDict[key] = self.similarityToDistance(alignment[0])

    def similarityToDistance(self, alignment):
        """
            Computes from the given similarity the distance measure.
        """
        sMax = self.similarity(alignment[0], alignment[0]) + self.similarity(alignment[1], alignment[1])
        sMax /= 2

        shuffle = []
        for a in [list(alignment[0]), list(alignment[1])]:
            shuffle.append(''.join(random.shuffle(a)))

        sRand = self.similarity(shuffle[0], shuffle[1])
        if sMax == sRand:
            sRand = sRand - 0.0001
        else:
            sEff = (self.similarity(alignment[0], alignment[1]) - sRand) / float(sMax - sRand)
        if sEff <= 0.0:
            return 1

        return -math.log10(sEff)

    def similarity(self, a, b):
        """
            Returns the similarity of two sequences a and b with the similarity
            score defined at the initialization.
        """
        similarity = 0
        for i in range(len(a)):
            similarity += self.score(a[i], b[i])
        return similarity

    def buildTree(self):
        """
            This function computes the phylogenetic tree with UPGMA and stores
            it in the Newick-Tree format.
        """
        upgma = UpgmaWpgma(self.distanceDict, len(self.sequences))
        upgma.computeClustering()
        self.newickTree = upgma.getNewickTree()

    def buildMultipleAlignment(self, firstGroup, secondGroup):
        """
            This function returns which is the best pairwise alignment out of
            all alignments of firstGroup and secondGroup.
        """
        highestScore = 0
        optimalAlignment = []
        for i in firstGroup:
            for j in secondGroup:
                nw = NeedlemanWunsch(i[0], j[0], 1)
                alignment = nw.compute()
                score = self.similarity(alignment[0], alignment[1])
                if highestScore < score:
                    highestScore = score
                    optimalAlignment = [alignment[0], alignment[1], i[1], j[1]]
        return optimalAlignment

    def computeOrderOfSequencesToAlign(self):
        """
            This function computes out of the phylogenetic tree in which order
            the sequences are aligned.
        """
        indexBegin = 0
        indexEnd = len(self.newickTree)
        while indexEnd != -1:
            indexBegin = self.newickTree.rfind('(', indexBegin, indexEnd)
            if indexBegin == -1:
                break
            i = indexBegin + 1
            stack = 0
            while stack >= 0 and i < len(self.newickTree):
                if self.newickTree[i] == '(':
                    stack += 1
                elif self.newickTree[i] == ')':
                    stack -= 1
                i += 1
            indexEnd = i

            firstGroup = ''
            secondGroup = ''
            substring = self.newickTree[indexBegin:indexEnd]
            if substring[1] != '(':
                firstGroupIndex = substring.find(',')
                firstGroup = substring[0:firstGroupIndex].strip(',')
                secondGroup = substring[firstGroupIndex:-1].strip(',')
            else:
                k = 1
                stack = 0
                while k < len(substring):
                    if substring[k] == '(':
                        stack += 1
                    elif substring[k] == ')':
                        stack -= 1
                    k += 1
                    if stack <= 0:
                        break
                firstGroup = substring[0:k].strip(',')
                secondGroup = substring[k:-1].strip(',')
            group0List = firstGroup.split(',')
            group1List = secondGroup.split(',')
            firstList = []
            secondList = []
            for j in group0List:
                firstList.append(int(j.strip('(').strip(')').strip(',')))
            for j in group1List:
                secondList.append(int(j.strip('(').strip(')').strip(',')))

            self.orderToAlign.append(sorted([sorted(firstList), sorted(secondList)]))
            indexEnd = indexBegin
            indexBegin = 0

    def computeMultipleAlignment(self):
        """
            This function returns the multiple sequence alignment.
        """
        self.computeAlignments()
        self.__computeDistanceDict()
        self.buildTree()
        self.computeOrderOfSequencesToAlign()
        i = 0
        indexAlignments = {}
        # create index to algnment realation
        while i < len(self.orderToAlign):
            if len(self.orderToAlign[i][0]) == 1 and len(self.orderToAlign[i][1]):
                for j in self.alignments:
                    if (j[1] == self.orderToAlign[i][0][0] and j[2] == self.orderToAlign[i][1][0]):
                        indexAlignments[self.orderToAlign[i][0][0]] = j[0][0]
                        indexAlignments[self.orderToAlign[i][1][0]] = j[0][1]
                        break
                    elif(j[1] == self.orderToAlign[i][1][0] and j[2] == self.orderToAlign[i][0][0]):
                        indexAlignments[self.orderToAlign[i][0][0]] = j[0][1]
                        indexAlignments[self.orderToAlign[i][1][0]] = j[0][0]
                        break
            elif len(self.orderToAlign[i][0]) == 1:
                indexAlignments[self.orderToAlign[i][0][0]] = self.sequences[self.orderToAlign[i][0][0]]
            elif len(self.orderToAlign[i][1]) == 1:
                try:
                    indexAlignments[self.orderToAlign[i][1][0]] = self.sequences[self.orderToAlign[i][1][0]]
                except:
                    print('Exception!')
                    print(f'i: {i}')
                    print(f'OrderToAlign: {self.orderToAlign}')
                    print(f'orderAlign: {self.orderToAlign[i][1][0]}')
                    print(self.sequences)
            i += 1

        for i in self.orderToAlign:
            # one sequence with one sequence
            if len(i[0]) == 1 and len(i[1]):
                indexAlignments[i[0][0]] = indexAlignments[i[0][0]].replace('-', 'X')
                indexAlignments[i[1][0]] = indexAlignments[i[1][0]].replace('-', 'X')
            # one sequence with one group
            # two groups
            else:
                firstGroup = []
                secondGroup = []
                for j in i[0]:
                    firstGroup.append([indexAlignments[j], j])

                for j in i[1]:
                    secondGroup.append([indexAlignments[j],j])

                pairwiseAlignment = self.buildMultipleAlignment(firstGroup, secondGroup)
                indexAlignments[pairwiseAlignment[2]] = pairwiseAlignment[0].replace('-', 'X')
                indexAlignments[pairwiseAlignment[3]] = pairwiseAlignment[1].replace('-', 'X')

                for j in i[0]:
                    nw = NeedlemanWunsch(pairwiseAlignment[0], indexAlignments[j], 1)
                    indexAlignments[j] = nw.compute()[1]

                for j in i[1]:
                    nw = NeedlemanWunsch(pairwiseAlignment[1], indexAlignments[j], 1)
                    indexAlignments[j] = nw.compute()[1]

                for j in indexAlignments:
                    indexAlignments[j] = indexAlignments[j].replace("-", "X")

        return indexAlignments
