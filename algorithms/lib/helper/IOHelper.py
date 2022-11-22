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
import os

class IOHelper():
    """
        Helper class for reading an writing files in different formats.
    """
    @staticmethod
    def readFastaFile(inputFileName):
        """
            Reads a given fasta file and returns it as a array.
            inputFileName: The path (relative or absolute) to the input fasta
                           file.
        """
        sequences = []
        if not os.path.exists(inputFileName):
            return sequences

        with open(inputFileName, 'r') as inputFile:
            sequence = ''
            for line in inputFile:
                if line.startswith('>'):
                    if len(sequence) > 0:
                        sequences.append(sequence)
                        sequence = ''
                    continue
                sequence += line.strip('\n')
        sequences.append(sequence)

        return sequences

    @staticmethod
    def writeFastaFile(alignments, outputFileName):
        """
            Writes a the given sequences to a file in the fasta format.
            sequences:      All computed alignemnts.
                            A list of lists with two elements: [[,],...,[,]].
            outputFileName: The path (relative or absolut) and the output file
                            name. e.g.: "/path/to/file" or "file" to write it
                            in the local directory.
        """
        if not outputFileName.endswith('.fas'):
            outputFileName += str('.fas')

        with open(outputFileName, 'w') as out:
            for i, alignmentSequences in enumerate(alignments):
                for j, sequence in enumerate(alignmentSequences):
                    out.write(f'>Alignment {i} sequence {j}\n')
                    out.write(f'{sequence}\n')

    @staticmethod
    def writeGraphMLFile(clusteredNodesDict, outputFileName):
        """
            Writes a tree computed by the UpgmaWpgma class in graphML-format to
            specified outputFileName.
        """
        if not outputFileName.endswith('.graphml'):
            outputFileName += str('.graphml')

        with open(outputFileName, 'w') as out:
            out.write(
'''<?xml version="1.0" encoding="UTF-8"?>
<graphml xmlns="http://graphml.graphdrawing.org/xmlns"
    xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"
    xsi:schemaLocation="http://graphml.graphdrawing.org/xmlns/1.0/graphml.xsd">
  <graph id="G" edgedefault="undirected">
''')

            for key in clusteredNodesDict:
                nodes = key.split(' ')
                out.write(f'\t\t<node id=\"{nodes[0]}\"/>\n')
                out.write(f'\t\t<node id=\"{nodes[1]}\"/>\n')
                out.write(f'\t\t<node id=\"{str(clusteredNodesDict[key])}\"/>\n')

            for i, key in enumerate(clusteredNodesDict):
                nodes = key.split(' ')
                out.write(f'\t\t<edge id=\"{i * 2}\" source=\"{nodes[0]}\" target=\"{clusteredNodesDict[key]}\"/>\n')
                out.write(f'\t\t<edge id=\"{i * 2 + 1}\" source=\"{nodes[0]}\" target=\"{clusteredNodesDict[key]}\"/>\n')

            out.write('\t</graph>\n</graphml>')

        # for i in clusteredNodesDictionary:
        #     nodes = i.split(" ")
        #     fileToWrite.write("\t\t<node id=\"" + nodes[0] + "\"/>\n")
        #     fileToWrite.write("\t\t<node id=\"" + nodes[1] + "\"/>\n")
        #     fileToWrite.write("\t\t<node id=\"" + str(clusteredNodesDictionary[i]) + "\"/>\n")
        # j = 0
        # for i in clusteredNodesDictionary:
        #     nodes = i.split(" ")
        #     fileToWrite.write("\t\t<edge id=\"" + str(j) + "\" source=\"" + nodes[0] + "\" target=\""+ str(clusteredNodesDictionary[i]) + "\"/>\n")
        #     j += 1
        #     fileToWrite.write("\t\t<edge id=\"" + str(j) + "\" source=\"" + nodes[1] + "\" target=\""+ str(clusteredNodesDictionary[i]) + "\"/>\n")
        #     j += 1

        # fileToWrite.write("\t</graph>\n</graphml>")
        # fileToWrite.close()

    @staticmethod
    def writeRnaDotBracketNotation(sequence, pairedBases, outputFileName):
        """
            Writes a given RNA sequence and the computed matching bases in
            dot-bracket notation to the file outputFileName.
        """
        stack = {}
        for i in range(len(sequence)):
            if i in pairedBases:
                stack[i], stack[pairedBases[i]] = '(', ')'
            elif i not in stack:
                stack[i] = '.'

        with open(outputFileName, 'w') as out:
            out.write(sequence + '\n')
            for key in sorted(stack):
                out.write(stack[key])

        return stack

    @staticmethod
    def writeNewickTree(newickTree, outputFileName):
        with open(outputFileName, 'w') as out:
            out.write(newickTree)
