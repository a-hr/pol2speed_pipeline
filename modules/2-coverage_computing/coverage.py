#!/usr/bin/env python

import sys
from pathlib import Path

import numpy as np
import deeptools.countReadsPerBin as crpb


def get_regions(bed):
    regions = Path(bed).read_text().splitlines()
    regions = [region.split("\t") for region in regions]
    return regions

def get_coverage(bams, regions):
    cr = crpb.CountReadsPerBin(bams, binLength=1000, stepSize=1000, numberOfProcessors=12)
    coverages_by_region = {}
    i = 0
    for chr_, start, end, enst, _, strand in regions:
        enst = enst.replace(": ", "-")
        region_name = f"{enst}_{chr_}_{start}_{end}_{strand}"
        coverages_by_region[region_name] = cr.count_reads_in_region(chr_, int(start), int(end)) # nbin x nbams
        i += 1

    return coverages_by_region


def main(bams_path, bed, output_dir):
    bams = list(Path(bams_path).glob("*.bam"))
    header = [bam.name for bam in bams]

    
    # save bam files to trace columns in the output matrix
    with open("bam_files.txt", "w") as f:
        for i, bam in enumerate(bams):
            f.write(f"{i}\t{bam}\n")

    # get target regions from bed file
    regions = get_regions(bed)

    # compute coverages
    print(f"Computing coverage for {len(bams)} bam files")
    coverages_by_region = get_coverage(bams, regions)

    # save coverages to a csv each
    print(f"Saving coverages to {output_dir}")
    for region, coverage in coverages_by_region.items():
        np.savetxt(str(output_dir / f"{region}.csv"), coverage, delimiter=",", header=header)

    print("Done computing coverage!")
    

if __name__ == "__main__":
    bams_path = sys.argv[1]
    bed = sys.argv[2]
    output_dir = Path(sys.argv[3]) / "coverages"

    output_dir.mkdir(exist_ok=True, parents=True)

    main(bams_path, bed, output_dir)
