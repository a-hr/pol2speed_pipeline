#!/usr/bin/env python
import sys


def main(infile: str, outfile: str) -> None:
    with open(infile, "r") as f, open(outfile, "w") as bed:
        # skip header
        f.readline()

        while (l := f.readline()):
            exon_start, exon_end, rank, chr, strand, gene_id, tr_id = l.split(";")
            strand = "+\n" if strand == "1" else "-\n"
            bed_str = "\t".join(
                ["chr" + chr.strip('\n'), exon_start, exon_end, f"{tr_id}: {rank}", "100", strand]
                )
            bed.writelines(bed_str)


if __name__ == "__main__":
    in_file = sys.argv[1]
    out_file = sys.argv[2]
    main(in_file, out_file)
