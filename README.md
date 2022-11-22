# Algorithms In Bioinformatics
To run the algorithms execute the file "main.py".

## Parameters

#### Help
  -h, --help

  Show this help message and exit

#### Algorithms

  -a {nw,gotoh,nw3,fengDoolittle,sumOfPairs,upgma,wpgma,nussinov},

  --algorithm {nw,gotoh,nw3,fengDoolittle,sumOfPairs,upgma,wpgma,nussinov}

  Define which algorithm should be executed. Options are:

  * 'nw' for the algorithm of Needleman and Wunsch.
  * 'gotoh' for the algorithm of Osamu Gotoh.
  * 'nw3' for the Needleman-Wunsch algorithm with three sequences.
  * 'fengDoolittle' for the heuristic multiple sequence alignment algorithm by Da-Fei Feng and Russell F. Doolittle.
  * 'sumOfPairs' for the scoring of a multiple sequence alignment by Humberto Carrillo and David Lipman.
  * 'upgma' or 'wpgma' is a clustering method to generate pylogenetic trees.
  * 'nussinov' for the RNA secondary structure prediction algorithm by Ruth
  Nussinov.

#### Input file

  -f INPUTFILE, --inputFile INPUTFILE

  Define the file in which the input sequences are defined. It have to be in fasta-format.

#### Input directory

  -d INPUTDIR, --inputdir INPUTDIR

  Define the directory in which the input sequence files are located. Files must be in fasta-format.

#### Output file

  -o OUTPUTFILE, --outputFile OUTPUTFILE

  Define in which file the output should be written. If
  not defined, it is written to "outputFile.fas" in the
  local directory.

#### Number of bases

  -b NUMBASES, --num_bases NUMBASES

  Define the number of bases to read for each input sequence.

#### Number of solutions

  --numberOfSolutions NUMBEROFSOLUTIONS

  Define the number of optimal solutions the Needleman-Wunsch algorithm should compute.

#### Output format

  --outputFormat {graphML,newickTree}

  Define the output format of the output file. This function is only parsed if you choose 'upgma' or 'wpgma' as an algorithm. Default is Newick tree.

#### Match score

  -m MATCH, --match MATCH

  Score for matching two nucleotides. Used in alignment.

#### Mismatch score

  -mm MISMATCH, --mismatch MISMATCH

  Score for mismatching two nucleotides. Used in alignment.

#### Gap cost

  -gc GAPCOST, --gap_cost GAPCOST

  Score for opening a gap in sequence. Is't used in Gotoh (use --gap_opening and --gap_enlargement instead).

#### Partial match score

  -pm PARTIALMATCH, --partial_match PARTIALMATCH

  Score for partial matching several nucleotides.

#### Gap opening (Gotoh)

  -go GAPOPENING, --gap_opening GAPOPENING

  Gap opening coefficient. Used in Gotoh.

#### Gap enlargement (Gotoh)
  -ge GAPENLARGEMENT, --gap_enlargement GAPENLARGEMENT

  Gap enlargement coefficient. Used in Gotoh.
