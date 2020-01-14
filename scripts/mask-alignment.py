"""
Mask initial bases from alignment FASTA
"""
import argparse
import Bio
import Bio.SeqIO
from Bio.Seq import Seq


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Mask initial bases from alignment FASTA",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--alignment", required=True, help="FASTA file of aligment")
    parser.add_argument("--mask-length", required=True, help="number of bases to mask from start")
    parser.add_argument("--output", required=True, help="FASTA file of output aligment")
    args = parser.parse_args()

    length = int(args.mask_length)

    with open(args.output, 'w') as outfile:
        for record in Bio.SeqIO.parse(args.alignment, 'fasta'):
            seq = str(record.seq)
            start = "N" * length
            end = seq[length:]
            record.seq = Seq(start + end)
            Bio.SeqIO.write(record, outfile, 'fasta')
