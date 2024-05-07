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
                
                tr_id, int_start, int_end, int_len, rank, intr_interval, chr_, strand = l.split("\t")
                strand = "+\n" if strand == "1\n" else "-\n"
                bed_str = "\t".join(
                    [chr_, int_start, int_end, f"{tr_id}: {rank}", "100", strand]
                    ).replace('"', '')
                bed.writelines(bed_str)


if __name__ == "__main__":
    in_file = sys.argv[1]
    out_file = sys.argv[2]
    main(in_file, out_file)
