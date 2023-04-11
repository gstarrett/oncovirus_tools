#!/usr/bin/env python3
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("g", type=str, help="The virus bedgraph file")
parser.add_argument("s", type=str, help="the bamindexstat file")
args = parser.parse_args()

# this can be simplified by just importing the bam file and then running bedtools and bamindexstat then processing those outputs

nums = list(range(1, 23))
chroms = ["chr" + str(num) for num in nums]
chroms.extend(["chrX", "chrY", "chrM", "X", "Y", "M", "MT"])
chroms.extend(str(num) for num in nums)

huReads = 0

path = args.s.split("/")
ext = path[-1].split(".")
lengths = {}
# stat file for human reads
virDict = {}

with open(args.g) as covFH:
    for line in covFH:
        f = line.split()
        if str(f[0]) in virDict:
            if int(f[3])>0:
                virDict[str(f[0])][0] += int(f[3]) * (int(f[2]) - int(f[1]))
                virDict[str(f[0])][1] += int(f[2]) - int(f[1])
        else:
            if int(f[3])>0:
                virDict.update({str(f[0]): [int(f[3]) * (int(f[2]) - int(f[1])), int(f[2]) - int(f[1]), 0]})

with open(args.s) as statFH:
    for line in statFH:
        f = line.split()
        if str(f[0]) in chroms:
            huReads += int(f[4])
        elif str(f[0]) != 'NoCoordinateCount=':
            lengths.update({str(f[0]): int(f[2])})
            if str(f[0]) in virDict:
                virDict[str(f[0])][2] = int(f[4])

# mcpyv.cov file for coverage of virus

for key in virDict:
    weightedCov = virDict[key][0]
    bases = virDict[key][1]
    reads= virDict[key][2]
    if (weightedCov > 0):
        avgCov = weightedCov/bases
    else:
        avgCov = 0
    normCp = (avgCov) / (huReads/1000)
    # viral genome coverage per 1000 human reads
    print("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(key, ext[0], normCp, huReads, bases, avgCov, reads, lengths[key]))
