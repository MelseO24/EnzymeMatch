import sys, os, csv
try:
    import pypdb  #needs v2.0
    if pypdb.get_info('7lzm') == None:
        sys.stderr.write("FATAL: pypdb does not work properly. Note pypdb version 2.0 is required.\n")
        sys.exit(1)
except:
    sys.stderr.write("FATAL: import pypdb failed. Try to install pypdb: 'pip install pypdb'\n")
    sys.exit(1)
import argparse
fout = sys.stdout

parser = argparse.ArgumentParser(description="Adds PDB annotation and UniProtID to EnzymeMatch output.")
parser.add_argument("-i", help="EnzymeMatch output file", required=True)
parser.add_argument("-o", help="name of output database", required=True)
parser.add_argument("--overwrite", help="Allow overwrite", action="store_true")
args = parser.parse_args()

fileInName = args.i 
new_db = args.o 
overwrite = args.overwrite
seperator = ' '

if not os.path.isfile(fileInName):
    sys.stderr.write(f"File {file} does not exist.\n")
    sys.exit(1)
elif os.path.isfile(new_db) and not overwrite:
    sys.stderr.write(f"File {new_db} exists, use --overwrite flag if you want to overwrite.\n")
    sys.exit(1)

fout.write("Conversion started, this may take a while...\n")
with open(fileInName, "r") as fileIn:
    csvreader = csv.reader(fileIn, delimiter="\t")
    with open(new_db, "w+") as outfile:
        csvwriter = csv.writer(outfile, delimiter="\t")
        for lineNr, csvLine in enumerate(csvreader):
            line = csvLine
            if line[0] == "PDB_ID":
                index_insert = line.index("Ligand_serial_number")
                csvwriter.writerow(line[:index_insert] +
                                   ["pdbx descriptor","PDB annotation"] + 
                                   line[index_insert:])
            elif lineNr == 0:
                sys.stderr.write("FATAL: First line of input database should start with column titles, starting with 'PDB_ID'")
                sys.exit(1)
            else:
                entryExists = True
                try:
                    annotation = pypdb.get_info(line[0])['struct']['title']
                except:
                    annotation = "unknown"
                    pdbx = "unknown"
                    entryExists = False

                if entryExists:
                    try:
                        pdbx = pypdb.get_info(line[0])['struct']['pdbx_descriptor'].split(",")[0]
                    except:
                        if 'pdbx_decriptor' in pypdb.get_info(line[0])['struct'].keys():
                            pdbx = pypdb.get_info(line[0])['struct']['pdbx_descriptor']
                        else:
                            pdbx = "unknown"

                csvwriter.writerow(line[:index_insert] + 
                                   [pdbx, annotation] +
                                   line[index_insert:])

fout.write(f"Annotated file written to {new_db}\n")

