# -*- coding: utf-8 -*-

# Global imports
import argparse as ap # To parse all the command-line arguments in a more fancy way

# Local imports
from SequenceOperator import *

# Script information
__author__ = "Sergi Rodà"
__license__ = "MIT"
__version__ = "1.0.1"
__maintainer__ = "Sergi Rodà"
__email__ = "sergi.rodallordes@bsc.es"

def parseArgs():
    """
    Parse arguments from command-line

    RETURNS
    -------
    input : string
              path of the input sequence file
    output : string
              path of the output sequence file
    operation : string
              path of the input sequence file
    trimleft : integer
              number of bases trimmed at the left
    trimright : integer
              number of bases trimmed at the left
    adaptor : string
              DNA adaptor sequence
    input2 : string
              path of the second input sequence file (alignment operation)
    """

    parser = ap.ArgumentParser()
    optional = parser._action_groups.pop()
    required = parser.add_argument_group('required arguments')

    required.add_argument("-i", "--input", required=True, metavar="FILE",
                          type=str, help="path of input file")
    required.add_argument("-o", "--output", required=True, metavar="FILE",
                          type=str, help="path of output file")
    required.add_argument("-O", "--operation", required=True, metavar="STRING",
                          type=str, help="operation flag (rc, trim, adaptor-removal, alignment)")
    optional.add_argument("-TL", "--trimleft", metavar="INTEGER",
                          type=int, help="number of bases trimmed at the left")
    optional.add_argument("-TR", "--trimright", metavar="INTEGER",
                          type=int, help="number of bases trimmed at the right")    
    optional.add_argument("-A", "--adaptor", metavar="STRING",
                          type=str, help="DNA adaptor sequence")
    optional.add_argument("-i2", "--input2", metavar="FILE",
                          type=str, help="path of 2nd input file")
    parser._action_groups.append(optional)
    args = parser.parse_args()

    Input_filename = args.input
    Output_filename = args.output
    Operation = args.operation
    TL, TR, Adaptor,Input_filename2 = None, None, None, ""

    if Operation.lower() == "trim":
      TL = args.trimleft
      TR = args.trimright
    
    elif Operation.lower() == "adaptor-removal" or Operation.lower() == "ar":
      Adaptor = args.adaptor

    elif Operation.lower() == "alignment" or Operation.lower() == "aln":
      Input_filename2 = args.input2

    return Input_filename, Output_filename, Operation, TL, TR, Adaptor, Input_filename2

def main():
    """
    Main function

    It is called when this script is the main program called by the interpreter
    """

    # Parse command-line arguments
    Input_filename, Output_filename, Operation, TL, TR, Adaptor, Input_filename2 = parseArgs()

    fasta_q = SequenceOperator(Input_filename, Output_filename, TL, TR, Adaptor, Input_filename2)

    # Execute the operation specified in the flag
    if Operation.lower() == "rc":
      fasta_q.ReverseComplementFasta()
    elif Operation.lower() == "trim":
      fasta_q.TrimFasta()
    elif Operation.lower() == "adaptor-removal" or Operation.lower() == "ar":
      fasta_q.AdaptorRemovalFasta()
    elif Operation.lower() == "alignment" or Operation.lower() == "aln":
      fasta_q.Alignment()
    else:
      print ("It doesn't exist such operation. The existing operations are rc(reverse-complement),trim "
                   "and adaptor-removal")


if __name__ == "__main__":
    """Call the main function"""
    main()
