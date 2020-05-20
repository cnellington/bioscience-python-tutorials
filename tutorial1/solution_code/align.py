"""
Caleb Ellington
2020/05/20
University of Washington
Department of Bioengineering
"""

import argparse
from os.path import isfile
from model import Alignment


"""
Usage: align.py <sequence1.fasta> <sequence2.fasta> [-p] [-v] [-s]
"""
def main():
    parser = argparse.ArgumentParser(description="Protein/DNA alignment")

    parser.add_argument(
        "seqfile1",
        action="store",
        help="Location of the first sequence file to align",
    )

    parser.add_argument(
        "seqfile2",
        action="store",
        help="Location of the second sequence file to align",
    )

    parser.add_argument(
        "--verbose",
        "-v",
        action="store_true",
        default=False,
        help="Display score and path matrices with resulting alignment",
    )

    parser.add_argument(
        "--p-value",
        "-p",
        action="store_true",
        default=False,
        help="Return alignment with p-value",
    )

    parser.add_argument(
        "--smith-waterman",
        "-s",
        action="store_true",
        default=False,
        help="Use Smith-Waterman local alignment instead of Needleman-Wunsch global alignment",
    )

    args = parser.parse_args()
    if not (isfile(args.seqfile1) and isfile(args.seqfile2)):
        print("Invalid sequence path")

    alignment = Alignment(args.seqfile1, args.seqfile2, gap_cost=-8)
    if args.smith_waterman:
        alignment.run_smith_waterman()
        print(alignment.get_local_alignment())
        print(f"score: {alignment.get_sw_score()}")
    else:
        alignment.run_needleman_wunsch()
        print(alignment.get_alignment())
        print(f"score: {alignment.get_nw_score()}")

    if args.p_value:
        n = 1e3
        print(f"P-value: {alignment.get_p_value(n=n)} (N={n})")

    if args.verbose:
        print("\nScore matrix (using BLOSUM62): ")
        print(alignment.get_score_matrix())
        print("\nPath traceback matrix [[0,2],[1,X]]: ")
        print(alignment.get_path_matrix())
        print()


if __name__ == '__main__':
    main()
