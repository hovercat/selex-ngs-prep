#!/usr/bin/env python

import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-i", type=str)
parser.add_argument("-o", type=str)
parser.add_argument("--min", type=int)
parser.add_argument("--mino", type=str)
parser.add_argument("--max", type=int)
parser.add_argument("--maxo", type=str)

args = parser.parse_args()

in_bounds = 0
too_short = 0
too_long = 0

with open(args.i, "r") as input_fastq, open(args.o, "w") as output_fastq, open(args.mino, "w") as output_fastq_short, open(args.maxo, "w") as output_fastq_long:
    while True:
        fastq_id = input_fastq.readline()
        fastq_seq = input_fastq.readline()
        fastq_desc = input_fastq.readline()
        fastq_quality = input_fastq.readline()

        if fastq_id is None or fastq_id == "":
            break

        seq = fastq_seq.strip()
        if len(seq) < args.min:
            current_file = output_fastq_short
            too_short += 1
        elif len(seq) > args.max:
            current_file = output_fastq_long
            too_long += 1
        else:
            current_file = output_fastq
            in_bounds += 1

        current_file.write(fastq_id)
        current_file.write(fastq_seq)
        current_file.write(fastq_desc)
        current_file.write(fastq_quality)
        current_file.flush()

print("Sequences, correct length:\t{}".format(in_bounds))
print("Sequences, too short:\t\t{}".format(too_short))
print("Sequences, too long:\t\t{}".format(too_long))