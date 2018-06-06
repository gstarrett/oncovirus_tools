#!/usr/bin/env python3
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("s", type=str, help="the bamindexstat file")
parser.add_argument("g", type=str, help="The virus bedgraph file")
args = parser.parse_args()

# this can be simplified by just importing the bam file and then running bedtools and bamindexstat then processing those outputs

chroms  = list(range(1, 22))
chroms.extend(["X","Y","MT"])

huReads = 0

path = args.s.split("/")
ext = path[-1].split(".")

# stat file for human reads
with open(args.s) as statFH:
  for line in statFH:
    f = line.split()
    if f[0] in chroms:
        huReads += int(f[4])

cov = 0
bases = 0
# mcpyv.cov file for coverage of virus
with open(args.g) as covFH:
  for line in covFH:
    f = line.split()
    cov += int(f[3])
    bases += int(f[2]) - int(f[1])

if (bases > 0):
    avgCov = cov/bases
else:
    avgCov = 0
normCp = (avgCov) / (huReads/1000)
# viral genome coverage per 1000 human reads
print(ext[0], normCp, huReads, bases, avgCov, cov)
