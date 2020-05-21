"""
Caleb Ellington
2020.05.22
"""

import argparse
from os.path import isdir, isfile

from services.helpers import *


def main():
    parser = argparse.ArgumentParser(description="Protein Alignment CLI")

    parser.add_argument(
        "sequence_file_1",
        action="store",
        help="Location of the first .fasta file to align"
    )

    parser.add_argument(
        "sequence_file_2",
        action="store",
        help="Location of the first .fasta file to align"
    )

    parser.add_argument(
        "-p",
        action="store_true",
        help="Include p-value measurement of significance"
    )

    args = parser.parse_args()
    if not (isfile(args.sequence_file_1) and isfile(args.sequence_file_2)):
        print("Invalid sequence path")
        exit(1)

    sequence1 = fasta_to_data(args.sequence_file_1)
    sequence2 = fasta_to_data(args.sequence_file_2)
    blosum = blosum_to_data('documents/blosum62.txt')
    align = Alignment(blosum, -8)
    print(align.get_alignment(sequence1, sequence2))
    if args.p:
        print(f"P-val: {round(align.get_p_value(sequence1, sequence2, n=1e2), 4)}")


if __name__ == '__main__':
    main()
