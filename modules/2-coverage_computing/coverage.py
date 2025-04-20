#!/usr/bin/env python

import sys
from pathlib import Path

from tqdm import tqdm

import numpy as np
import deeptools.countReadsPerBin as crpb


def get_regions(bed):
    regions = Path(bed).read_text().splitlines()
    regions = [region.split("\t") for region in regions]
    return regions


def get_coverage(bams, regions, od, header):
    cr = crpb.CountReadsPerBin(
        bams, binLength=1000, stepSize=1000, numberOfProcessors=10
    )

    for chr_, start, end, enst, _, strand in tqdm(
        regions, file=sys.stdout, colour="green"
    ):
        enst = enst.replace(": ", "-")
        region = f"{enst.replace(".", "-")}_chr{chr_}_{start}_{end}_{strand}"

        # nbin x nbams
        result = cr.count_reads_in_region(f"chr{chr_}", int(start), int(end))[0]
        np.savetxt(str(od / f"{region}.csv"), result, delimiter="\t", header="\t".join(header))


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
    get_coverage(bams, regions, output_dir, header)

    print("Done computing coverage!")


if __name__ == "__main__":
    bams_path = sys.argv[1]
    bed = sys.argv[2]
    output_dir = Path(sys.argv[3]) / "coverages"

    output_dir.mkdir(exist_ok=True, parents=True)

    main(bams_path, bed, output_dir)
