#!/usr/bin/env python
import sys

def main(file: str) -> None:
    outname = file.split("/")[-1].split(".")[0] + ".bed"
    with open(file, "r") as f, open(outname, "w") as bed:
        # skip header
        f.readline()
        while (l := f.readline()):
            tr_id, int_start, int_end, int_len, intr_interval, chr_, rank, strand = l.split(";")
            strand = "+\n" if strand == "1\n" else "-\n"
            bed_str = "\t".join(
                [chr_, int_start, int_end, f"{tr_id}: {rank}", "100", strand]
                ).replace('"', '')
            bed.writelines(bed_str)


if __name__ == "__main__":
    files = sys.argv[1:]
    for file in files:
        main(file)
