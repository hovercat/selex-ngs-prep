#!/usr/bin/env python
import os
import sys
import pandas as pd
import argparse

alphabet_DNA = "ACGT"
alphabet_RNA = "ACGU"
alphabet_PROT = "ABCDEFGHIKLMNPQRSTUVWYZ"


# Takes a fasta file and a length
def analyse_round_composition(in_file, alphabet, aptamer_n_region, percentages, collapse):
    nt_dict = {nt_position: {c: 0 for c in alphabet} for nt_position in range(0, aptamer_n_region)}
    try:
        in_stream = None
        if in_file == "-" or in_file is None:
            in_stream = sys.stdin
        else:
            in_stream = open(in_file, "r")

        for desc in in_stream:
            seq = in_stream.readline().rstrip('\n')
            if len(seq) == aptamer_n_region:
                for i in range(0, aptamer_n_region):
                    nt_dict[i][seq[i]] += 1
    finally:
        in_stream.close()

    df_composition = pd.DataFrame(nt_dict)

    if collapse:
        df_composition = df_composition.sum(axis=1)
        if percentages:
            df_composition = df_composition / df_composition.sum()
    elif percentages:
        df_composition = df_composition / df_composition.sum(axis=0)[0]

    return df_composition


def parse_args(args):
    # ======= Parameter validation =======
    in_file = args.in_file
    if in_file != "-" and not os.path.exists(in_file):
        print("Input file {} does not seem to be at the specified path.".format(in_file))
        print("Halting.")
        return 1

    force = args.force
    out_file = args.out_csv
    if out_file != "-" and not force and os.path.exists(out_file):
        print("Output file {} already exists. Add --force to overwrite the output file".format(
            out_file))
        print("Halting.")
        return 1

    n_region = args.n_region

    alphabet = args.alphabet
    dna = args.DNA
    rna = args.RNA
    prot = args.PROT
    if not ((bool(dna) ^ bool(rna) ^ bool(prot)) or not (bool(dna) or bool(rna) or bool(prot))):
        print("\'--DNA\', \'--RNA\' and \'--PROT\' are exclusive. Only one argument is allowed.")
        print("Halting.")
        return 1
    if dna:
        alphabet = alphabet_DNA
    if rna:
        alphabet = alphabet_RNA
    if prot:
        alphabet = alphabet_PROT

    collapse = args.general_composition

    return in_file, out_file, n_region, alphabet, args.percent, collapse


def analyse_selex_round(in_file, out_file, n_region, alphabet, percentages, collapse):  # this is wonky
    # ======= Round composition =======
    df_composition = analyse_round_composition(in_file, alphabet, n_region, percentages, collapse)

    # ======= Printing =======
    try:
        out_stream = None
        if out_file == "-" or out_file is None:
            out_stream = sys.stdout
        else:
            out_stream = open(out_file, "w")

        if collapse:
            round_name = os.path.splitext(os.path.basename(in_file))[0]
            df_composition = pd.DataFrame({round_name: df_composition})
            df_composition = df_composition.transpose()
            df_composition.to_csv(out_stream, sep="\t", index_label="round_name")
        else:
            df_composition = df_composition.transpose()
            df_composition.to_csv(out_stream, sep="\t", index_label="nt_position")
    finally:
        out_stream.close()


args_parser = argparse.ArgumentParser(
    prog="SELEXntcomposition",
    description="",
    formatter_class=argparse.RawTextHelpFormatter
)

args_parser.add_argument("-i", "--in-file", type=str, default="-", help="Input file in FASTA format. Default: stdin")
args_parser.add_argument("-o", "--out-csv", type=str, default="-", help="Output csv file location. Default: stdout")
args_parser.add_argument("--force", action="store_true", default=False,
                         help="Overwrite output files if they already exists")
args_parser.add_argument("--DNA", action="store_true", default=False, help="Shortcut for --alphabet ACGT.")
args_parser.add_argument("--RNA", action="store_true", default=False, help="Shortcut for --alphabet ACGU")
args_parser.add_argument("--PROT", action="store_true", default=False,
                         help="Shortcut for --alphabet ABCDEFGHIKLMNPQRSTUVWYZ")
args_parser.add_argument("-a", "--alphabet", type=str, default="ACGT",
                         help="Define FASTA letter alphabet. Default is DNA (ACGT).")
args_parser.add_argument("-n", "--n-region", type=int, required=True,
                         help="All sequences of length n are considered. Any other sequence is discarded.")
args_parser.add_argument("-p", "--percent", action="store_true",
                         help="If present the values will be divided by the number of sequences.")
args_parser.add_argument("-g", "--general-composition", action="store_true", help="If set, the round's composition "
                                                                                  "will be printed out non-position "
                                                                                  "dependent.")


def main():
    args = args_parser.parse_args()
    in_file, out_file, n_region, alphabet, percent, collapse = parse_args(args)
    analyse_selex_round(in_file, out_file, n_region, alphabet, percent, collapse)


main()
