#!/usr/bin/env python
import sys


def main(infile: str, outfile: str) -> None:
    with open(infile, "r") as f:
        with open(outfile, "w") as bed:
            # skip header
            f.readline()

            while True:
                l = f.readline()
                if not l:
                    break
                
                exon_start, exon_end, rank, chr, strand, gene_id, tr_id = l.split(";")
                strand = "+" if strand == "1" else "-"
                tr_id = tr_id.strip("\n")
                elements = [chr, exon_start, exon_end, f"{tr_id}: {rank}", "100", strand]
                elements = [e.strip("\n") for e in elements]
                bed_str = "\t".join(elements) + "\n"
                bed.writelines(bed_str)


if __name__ == "__main__":
    in_file = sys.argv[1]
    out_file = sys.argv[2]
    main(in_file, out_file)
