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
# main class
import os
import sys
import argparse

from lib.helper.IOHelper import IOHelper as io
from lib.helper.MultipleAlignmentHelper import MultipleAlignmentHelper as mah
from lib.pairwise.NeedlemanWunsch import NeedlemanWunsch as nw
from lib.pairwise.Gotoh import Gotoh
from lib.multiple.NeedlemanWunschN3 import NeedlemanWunschN3 as NW3
from lib.multiple.UpgmaWpgma import UpgmaWpgma
from lib.multiple.FengDoolittle import FengDoolittle
from lib.multiple.SumOfPairs import SumOfPairs
from lib.structurePrediction.Nussinov import Nussinov

def main():
    """
        Method to parse the arguments and start the defined algorithms.
    """

    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-a', '--algorithm',
        choices=[
            'nw', 'gotoh',
            'nw3','fengDoolittle', 'sumOfPairs',
            'upgma', 'wpgma', 'nussinov'
        ],
        required=True,
        help='''
Define which algorithm should be executed.
Options are: 'nw' for the algorithm of Needleman and Wunsch,
             'gotoh' for the algorithm of Osamu Gotoh,
             'nw3' for the Needleman-Wunsch algorithm with three sequences,
             'fengDoolittle' for the heuristic multiple sequence alignment
                             algorithm by Da-Fei Feng and Russell F. Doolittle,
             'sumOfPairs' for the scoring of a multiple sequence alignment by
                          Humberto Carrillo and David Lipman.
             'upgma' or 'wpgma' is a clustering method to generate phylogenetic
                                trees,
             'nussinov' for the RNA secondary structure prediction algorithm
                        by Ruth Nussinov.
''')

    parser.add_argument(
        '-f', '--inputFile',
        dest='inputFile',
        help='''
Define the file in which the input sequences are defined.
It have to be in fasta-format.
''')

    parser.add_argument(
        '-o', '--outputFile',
        help='''
Define in which file the output should be written. If not defined, it is
written to 'outputFile.fas' in the local directory.
''')

    # parser.add_argument("-gc", "--gapCost", dest="gapCost",
    #                     help="Name of a gap function definde in class PairwiseAligmentHelper.")

    parser.add_argument(
        '--numSolutions',
        dest='numSolutions',
        help='''
Define the number of optimal solutions the Needleman-Wunsch algorithm
should compute.
''')

    parser.add_argument(
        '--outputFormat',
        dest='outputFormat',
        choices=['graphML', 'newickTree'],
        help='''
Define the output format of the output file. This function is only parsed if
you choose 'upgma' or 'wpgma' as an algorithm. Default is 'Newick tree'.
''')

    # parser.add_argument("--scoring", dest="similarityScore",
    #                     help="Name of a similarity score defined in class PairwiseAligmentHelper. Per default "
    #                          "\"pam\" and \"blosum\" (pam250 and blosum62) are implemented. Feel free to extend, you can find the "
    #                          "file \"PairwiseAligmentHelper.py\" in lib/helper. If this option is not defined, the pam250 matrix is choosen.")

    parser.add_argument(
        '--gapPenalty',
        dest='gapPenalty',
        help='''
Define a gap penalty. Default for pam is 8 and blosum 6.
''')

    args = parser.parse_args()

    outputFile = ''

    if args.outputFile:
        outputFile = args.outputFile

    # if args.similarityScore:
    #     weightFunction = args.similarityScore

    sequences = getSequencesFromFile(args.inputFile)

    if len(sequences) > 1:
        # pairwise alignment
        if args.algorithm == 'nw':
            if outputFile == '':
                outputFile = 'needlemanWunsch.fas'

            numSolutions = -1
            if args.numSolutions:
                numSolutions = int(args.numSolutions)
            needlemanWunsch(
                sequences[0:2],
                outputFile = outputFile,
                maxSolutions=numSolutions
            )

        elif args.algorithm == 'gotoh':
            if outputFile == '':
                outputFile = 'gotoh.fas'
            gotoh(sequences[0:2], outputFile = outputFile)

        # multiple alignment
        elif args.algorithm == 'upgma' or args.algorithm == 'wpgma':
            newickTree = True
            if args.outputFormat == 'graphML':
                newickTree = False
            if outputFile == '':
                if args.algorithm == 'upgma':
                    outputFile = 'upgma'
                else:
                    outputFile = 'wpgma'
            upgmaWpgma(args.algorithm == "upgma", sequences, outputFile, newickTree)

        elif args.algorithm == 'fengDoolittle':
            if outputFile == '':
                outputFile = 'fengDoolittle.fas'
            fengDoolittle(sequences, outputFile)

        elif args.algorithm == 'sumOfPairs':
            sumOfPairs(sequences)

        elif args.algorithm == 'nw3':
            if not (len(sequences) == 3):
                print('Wrong number of input sequences.')
                print('Needleman-Wunsch n=3 requires exactly three sequences')
                print(f'{len(sequences)} sequences are given.')
                sys.exit()

            numSolutions = -1
            if args.numSolutions:
                numSolutions = int(args.numSolutions)

            if outputFile == '':
               outputFile = 'nw3.fas'
            needlemanWunschN3(
                sequences[0:3],
                outputFile = outputFile,
                maxSolutions=numSolutions
            )

    elif len(sequences) == 1:
        # structure prediction
        if args.algorithm == 'nussinov':
            if outputFile == '':
                outputFile = 'nussinov.dotBracket'
            nussinov(sequences[0:1], outputFile)
        else:
            print('You have defined only one input sequence.')
            print(f'But the algorithm \'{args.algorithm}\' requires two')

    else:
        print('No sequences are defined in input file.')
        sys.exit(0)

def getSequencesFromFile(inputFile):
    """
        Parse the input file to get the sequences.
        Returns the sequences as an array.
        inputFile: A fasta format file with the input sequences.
    """
    return io.readFastaFile(inputFile)

def needlemanWunsch(sequences, outputFile, maxSolutions=-1):
    """
        Executes the Needleman-Wunsch algorithm with a default score function
        defined as: a == b -> 0 and a !=b --> 1. Stores the alignments per
        default in file needlemanWunsch.fas.
        outputFile:         The path to the output file.
    """
    print('\nThe following sequences are given:')
    for sequence in sequences:
        print(sequence)

    print('\nComputing solution...\n')
    nwObj = nw(sequences[0], sequences[1], maxSolutions)
    result = nwObj.compute()

    print(f'Number of optimal solutions: {len(result)}')
    print(f'Score: {nwObj.matrix[-1][-1]}')
    print(f'One solution is:\n{result[0][0]}\n{result[0][1]}')

    io.writeFastaFile(result, outputFile)
    print('\nFor more solutions look in the file \'needlemanWunsch.fas\' '\
          'in the bin directory.')

def gotoh(sequences, outputFile='gotoh.fas'):
    """
        Executes the Gotoh algorithm with a default score function defined
        as: a == b -> 0 and a !=b --> 1 and a cost function defined as:
        g(x) = 2 + k. Stores the alignments per default in file gotoh.fas.
        outputFile: The path to the output file.
    """
    print('\nThe following sequences are given:')
    for i in sequences:
        print(i)

    print('\nComputing solution...\n')
    gotoh = Gotoh(sequences[0], sequences[1])
    result = gotoh.compute()

    score = max(
        gotoh.matrices[0][-1][-1],
        gotoh.matrices[1][-1][-1],
        gotoh.matrices[2][-1][-1],
    )
    print(f'Number of solutions: {len(result)}')
    print(f'Score: {-score}')
    print(f'One solution is:\n{result[0][0]}\n{result[0][1]}')

    io.writeFastaFile(result, outputFile)
    print('\nFor more solutions look in the file \'gotoh.fas\' in the bin directory.')

def needlemanWunschN3(sequences, outputFile, maxSolutions=-1):
    """
        Executes the Needleman-Wunsch algorithm with three sequences
    """
    print('\nThe following sequences are given:')
    for sequence in sequences:
        print(sequence)

    print('\nComputing solution...\n')
    nw3 = NW3(sequences[0], sequences[1], sequences[2], maxSolutions)
    result = nw3.compute()

    print(f'Score: {nw3.matrix[-1][-1][-1]}')
    print(f'Number of optimal solutions: {len(result)}')
    print(f'One solution is:\n{result[0][0]}\n{result[0][1]}\n{result[0][2]}')

    io.writeFastaFile(result, outputFile)
    print('\nFor more solutions look in the file \'nw3.fas\' in the bin directory.')


def upgmaWpgma(upgmaWpgma, sequences, outputFile, fileFormat):
    """
        Executes the a phylogenetic clustering with a upgm or wpgm weighting.
        sequences:  All defined input sequences as a list.
        outputFile: The name of the output file
        fileFormat: The file format of the output file
    """
    print('The following sequences are given:')
    for sequence in sequences:
        print(sequence)
    print('Computing clustering...')

    data = mah.createDataForUpgmaWpgma(sequences)

    if upgmaWpgma:
        upgma = UpgmaWpgma(data[0], len(data[1]))
        upgma.computeClustering()

        if not fileFormat:
            outputFile += '.graphML'
            io.writeGraphMLFile(upgma.mapping, outputFile)
            print(f'Clustering written as graphML file: {os.path.abspath(outputFile)}')

        else:
            outputFile += '.newickTree'
            cluster = upgma.getNewickTree(with_edge_weights=True)

            io.writeNewickTree(cluster, outputFile)
            print(f'Computed upgma cluster: {cluster}')
            print(f'The clustering was also written to: {os.path.abspath(outputFile)}')

    else:
        wpgma = UpgmaWpgma(data[0], len(data[1]), False, data[2])
        wpgma.computeClustering()

        if not fileFormat:
            outputFile += '.graphML'

            io.writeGraphMLFile(wpgma.mapping, outputFile)
            print(f'Clustering written as graphML file: {os.path.abspath(outputFile)}')

        else:
            outputFile += '.newickTree'
            cluster = wpgma.getNewickTree(with_edge_weights=True)

            io.writeNewickTree(cluster, outputFile)
            print(f'Computed wpgma cluster: {cluster}')
            print(f'The clustering was also written to: {os.path.abspath(outputFile)}')


def nussinov(sequence, outputFile):
    """
        Executes the RNA-folding algorithm from Nussinov.
        sequence:   The RNA-sequnce as a list.
        outputFile: The name of the output file.
    """
    print('\nThe following sequence is given:')
    print(f'{sequence[0]}\n')

    nussinov = Nussinov(sequence[0])
    nussinov.execute()

    print('\nDot-bracket: ')
    io.writeRnaDotBracketNotation(sequence[0], nussinov.pairedBases, outputFile)
    print(f'The result was also written to: {os.path.abspath(outputFile)}')

def sumOfPairs(sequences):
    """
        This method scores a multiple sequence alignment with
        the sum of pairs algorithm.
        sequences:       The multiple sequence alignment.
    """
    print('The following sequences are given:')
    for sequence in sequences:
        print(sequence)

    sop = SumOfPairs(sequences)
    print(f'Sum-of-pairs scoring: {sop.execute()}')

def fengDoolittle(sequences, outputFile):
    """
        Executes the heuristic multiple sequence alignment by Feng and Doolittle.
        sequences:       All input sequnces to align.
        outputFile:      The output file name.
    """
    fd = FengDoolittle(sequences)
    alignmentDict = fd.computeMultipleAlignment()
    alignment = [[]]

    for value in alignmentDict.values():
        alignment[0].append(value)

    io.writeFastaFile(alignment, outputFile)

    print('Input sequences:\n')
    for sequence in sequences:
        print(sequence)

    print('\nAlignment:')
    for value in alignmentDict.values():
        print(value)

    print(sumOfPairs(alignment[0]))

if __name__ == '__main__':
        main()
