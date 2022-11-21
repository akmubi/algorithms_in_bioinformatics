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

class UpgmaWpgma():
    """
        Upgma/Wpgma is a clustering method to generate phylogenetic trees.
    """

    def __init__(self,
                 distances,
                 nodeCount,
                 useUpgma=True,
                 sequenceSizeMapping={}):
        """
            To initalize a object of this class, please define the following:
            distances:           A dictionary with the distance between two
                                 sequences. Should have the form
                                 'Key0 Key1':distance. The key0 and key1 have
                                 to be integers.
            nodeCount:           The number of sequences.
            useUpgma:            If True, the upgma weighting is used,
                                 if False, wpgma.
            sequenceSizeMapping: Only necessary if wpgma is executed. It
                                 defines the size of each sequence. Should
                                 have the form: 'Key:len(sequence)'
        """
        self.distances = distances
        self.mapping = {}
        self.nodeCount = nodeCount
        self.numNodes = nodeCount
        self.useUpgma = useUpgma
        self.seqToSize = sequenceSizeMapping
        self.edgeWeight = {}

    def computeClustering(self):
        """
            This function computes the clustering to get the phylogenetic tree.
        """
        done = False
        j = 0
        while not done:
            j += 1
            minKey, minValue = self.computeMinDistance()
            nodes = minKey.split(' ') # ['key0', 'key1']

            if len(nodes) > 1:
                self.mapping[minKey] = self.nodeCount
                self.computeEdgeWeight(minValue, nodes)

                if minKey in self.distances:
                    del self.distances[minKey]

                for i in range(0, self.nodeCount + 1):
                    keys = self.keyInDict(f'{nodes[0]} {i}', f'{nodes[1]} {i}')
                    if keys[0] != '':
                        # new distance
                        self.distances[f'{i} {self.nodeCount}'] = \
                            self.computeNewDistance(
                                self.distances[keys[0]],
                                self.distances[keys[1]],
                                nodes[0],
                                nodes[1]
                            )
                        if not self.useUpgma:
                            self.seqToSize[self.nodeCount] = \
                                self.seqToSize[int(nodes[0])] + \
                                self.seqToSize[int(nodes[1])]
                        del self.distances[keys[0]]
                        del self.distances[keys[1]]
                self.nodeCount += 1
            else:
                done = True

    def keyInDict(self, key1, key2):
        """
            Returns True if the given keys are in the distance dictionary,
            False otherwise.
            key1: The first key value.
            key2: The second key value.
        """
        rkey1 = key1[::-1]
        rkey2 = key2[::-1]

        if key1 in self.distances and key2 in self.distances:
            return [key1, key2]

        elif rkey1 in self.distances and key2 in self.distances:
            return [rkey1, key2]

        elif key1 in self.distances and rkey2 in self.distances:
            return [key1, rkey2]

        elif rkey1 in self.distances and rkey2 in self.distances:
            return [rkey1, rkey2]

        return ['', '']

    def computeMinDistance(self):
        """
            Returns the next two clusters for merging.
        """
        minKey   = ''       # 'key0 key1'
        minValue = math.inf # distance

        for key, value in self.distances.items():
            if value < minValue:
                minValue = value
                minKey = key

        return minKey, minValue

    def computeNewDistance(self, distanceAX, distanceBX, indexA, indexB):
        """
            Returns the new distance between the new merged cluster and
            another cluster.
            distanceAX: The old distance between cluster a and x.
            distanceBX: The old distance between cluster b and x.
            indexA:     The index of a.
            indexB:     The index of b.
        """
        if self.useUpgma:
            return self.upgmaDistance(distanceAX, distanceBX)
        else:
            return self.wpgmaDistance(
                distanceAX,
                distanceBX,
                self.seqToSize[int(indexA)],
                self.seqToSize[int(indexB)]
            )

    def upgmaDistance(self, distanceAX, distanceBX):
        """
            Returns the upgma-distance between the new merged cluster a and
            another cluster x.
            distanceAX: The old distance between cluster a and x.
            distanceBX: The old distance between cluster b and x.
        """
        return (distanceAX + distanceBX) / 2

    def wpgmaDistance(self, distanceAX, distanceBX, lengthA, lengthB):
        """
            Returns the wpgma-distance between the new merged cluster a and
            another cluster x.
            distanceAX: The old distance between cluster a and x.
            distanceBX: The old distance between cluster b and x.
            lengthA:    The index of a.
            lengthB:    The index of b.
        """
        return (lengthA * distanceAX + lengthB * distanceBX) / (lengthA + lengthB)

    def computeEdgeWeight(self, weight, nodes):
        """
            This method computes the new edge weight for a new cluster.
            weight: The edge weight equal to the distance of the to merged
                    clusters.
            nodes:  A list containing the indices of the two merged clusters.
        """
        node0, node1 = int(nodes[0]), int(nodes[1])

        if node0 < self.numNodes and node1 < self.numNodes:
            self.edgeWeight[self.nodeCount] = [
                weight / 2.0,
                weight / 2.0,
            ]

        elif node0 < self.numNodes:
            self.edgeWeight[self.nodeCount] = [
                weight / 2.0,
                weight / 2.0 - self.edgeWeight[node1][1],
            ]

        elif node1 < self.numNodes:
            self.edgeWeight[self.nodeCount] = [
                weight / 2.0 - self.edgeWeight[node0][1],
                weight / 2.0,
            ]

        else:
            self.edgeWeight[self.nodeCount] = [
                weight / 2.0 - self.edgeWeight[node0][1],
                weight / 2.0 - self.edgeWeight[node1][1],
            ]

    def getNewickTree(self, widthEdgeWeights=False):
        """
            Returns the computed cluster in the Newick tree format.
            widthEdgeWeights: If True, edge weights are part of the output,
                              if False, not.
        """
        # expectedValue = {"2 3": 5, "0 1": 7, "4 5": 6, "6 7": 8}
        newickDict = dict([[v, k] for k, v in self.mapping.items()])
        if widthEdgeWeights:
            for i in newickDict:
                if i in self.edgeWeight:
                    nodesWithWeights = newickDict[i].split(' ')
                    nodesWithWeights[0] = nodesWithWeights[0].strip(' ')
                    nodesWithWeights[0] += f':{self.edgeWeight[i][1]}'
                    nodesWithWeights[1] = nodesWithWeights[1].strip(' ')
                    nodesWithWeights[1] += f':{self.edgeWeight[i][0]}'
                    newickDict[i] = f'{nodesWithWeights[0]} {nodesWithWeights[1]}'
            self.mapping = dict([[v, k] for k, v in newickDict.items()])

        for i in self.mapping:
            index = -1
            leadingSequence = True
            for j in newickDict:
                stringToFind = f' {self.mapping[i]}'# + ''
                if newickDict[j].find(stringToFind) != -1:
                    index = j
                    leadingSequence = False
                    break

                stringToFind = f'{self.mapping[i]} '
                if newickDict[j].find(stringToFind) != -1:
                    index = j
                    leadingSequence = True
                    break

                if widthEdgeWeights:
                    stringToFind = f'{self.mapping[i]}:'
                else:
                    stringToFind = f'{self.mapping[i]},'

                if newickDict[j].find(stringToFind) != -1:
                    index = j
                    leadingSequence = True
                    break

                stringToFind = f',{self.mapping[i]}'
                if newickDict[j].find(stringToFind) != -1:
                    index = j
                    leadingSequence = False
                    break

            if index != -1:
                key = int(stringToFind.strip().strip(',').strip(':'))
                value = newickDict[key].replace(' ', ',')
                if leadingSequence:
                    stringToReplace = f'({value}):'
                else:
                    stringToReplace = f',({value})'

                newickDict[index] = \
                    newickDict[index].replace(stringToFind, stringToReplace) \
                                     .replace(',,', ',')
                del newickDict[key]

        for i in newickDict:
            return f'({newickDict[i]})'
