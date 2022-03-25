import sys, os
import argparse

parser = argparse.ArgumentParser(description="Merges EnzymeMatch outputfiles")
parser.add_argument("-i", metavar="inputFiles", help="EnzymeMatch output files to merge", nargs="+", required=True)
parser.add_argument("-o", metavar="outFile", help="Output file", required=True)
parser.add_argument("--overwrite", help="Allow overwrite", action="store_true")
args = parser.parse_args()

for filename in args.i:
    if not os.path.isfile(filename):
        sys.stderr.write(f"File {filename} does not exist.\n")
        sys.exit(1)

if os.path.isfile(args.o):
    if not args.overwrite:
        sys.stderr.write(f"File {args.o} already exist, provide --overwrite flag to allow overwrite.\n")
        sys.exit(1)
        
#parse databases
remember_lines = {}
for fileNr, filename in enumerate(args.i):
    with open(filename, "r") as f:
        for line in f.readlines():
            pdbID = line.split("\t")[0]
            if pdbID == "PDB_ID":
                if fileNr == 0:
                    header = line
                else:
                    if not header == line:
                        sys.stderr.write("Headers of file differ, this cannot be merged. Exiting...\n")
                        sys.exit(1)
            elif pdbID in remember_lines.keys():
                sys.stdout.write(f"Entry {pdbID} from file {args.i[fileNr]} is skipped, since this is a duplicate.\n")
            else:
                remember_lines[pdbID] = line

#write merged database
with open(args.o, "w+") as fileOut:
   fileOut.write(header)
   for entry in remember_lines.keys():
       fileOut.write(remember_lines[entry])

sys.stdout.write(f"Merged files {args.i} written to {args.o}\n")

