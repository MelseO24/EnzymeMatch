    # EnzymeMatch
    # Copyright (C) 2022  Okke Melse
    #
    # This program is free software: you can redistribute it and/or modify
    # it under the terms of the GNU General Public License as published by
    # the Free Software Foundation, either version 3 of the License, or
    # (at your option) any later version.
    #
    # This program is distributed in the hope that it will be useful,
    # but WITHOUT ANY WARRANTY; without even the implied warranty of
    # MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    # GNU General Public License for more details.
    #
    # You should have received a copy of the GNU General Public License
    # along with this program.  If not, see <https://www.gnu.org/licenses/>.


"""
This script updates the BioLip database.
Make sure you have a file "updates.txt" in your working directory, listing all previously performed BioLip updates
"""

import sys, os, re, argparse
import urllib.request
from datetime import datetime
import tarfile
fout = sys.stdout


overview_url = "https://zhanglab.ccmb.med.umich.edu/BioLiP/weekly.html"
download_url_base = "https://zhanglab.ccmb.med.umich.edu/BioLiP/weekly"
download_url_anno_before2013 = "https://zhanglab.ccmb.med.umich.edu/BioLiP/download/BioLiP_nr.tar.bz2"   # DO NOT CHANGE!
download_url_rec_before2013 = "https://zhanglab.ccmb.med.umich.edu/BioLiP/download/receptor_nr.tar.bz2"  # DO NOT CHANGE!
download_url_lig_before2013 = "https://zhanglab.ccmb.med.umich.edu/BioLiP/download/ligand_nr.tar.bz2"    # DO NOT CHANGE!

updated_weeks = []

parser = argparse.ArgumentParser(description="Script updating BioLiP database")
parser.add_argument("-ano", metavar="folder", help="Folder to save annotations files in, default = annotations", default="annotations")
parser.add_argument("-append", metavar="BioLiP DB", help="Main BioLiP database to add updates in")
parser.add_argument("--wholeDB", help="Download entire BioLiP database, including the database before 2013--03-06", action="store_true")
args = parser.parse_args()

folder_ano = args.ano
folder_rec = "receptor_nr"
folder_lig = "ligand_nr"
append_db = args.append
wholeDB = args.wholeDB

if not os.path.isdir(folder_ano):
    os.mkdir(folder_ano)


if wholeDB:
    if append_db:
        sys.stderr.write("cannot combine --wholeDB and -append flag.\n")
        sys.exit(1)
    elif os.path.isdir(folder_rec):
        sys.stderr.write(f"Folder '{folder_rec}' found. This folder should not be present when the --wholeDB flag is provided, since this is created by the script.\n")
        sys.exit(1)
    elif os.path.isdir(folder_lig):
        sys.stderr.write(f"Folder '{folder_lig}' found. This folder should not be present when the --wholeDB flag is provided, since this is created by the script.\n")
        sys.exit(1)
    elif os.path.isfile("updates.txt"):
        sys.stderr.write("'updates.txt' found, please remove this file, since this script creates this file when --wholeDB flag provided\n")
        sys.exit(1)
    else:
        os.mkdir(folder_rec)
        os.mkdir(folder_lig)
else:
    if not os.path.isdir(folder_rec):
        sys.stderr.write(f"Folder '{folder_rec}' not found. Please make sure that the folder containing the receptor pdb files are in a folder named {folder_rec}\n")
        sys.exit(1)
    elif not os.path.isdir(folder_lig):
        sys.stderr.write(f"Folder '{folder_lig}' not found. Please make sure that the folder containing the ligand pdb files are in a folder named {folder_lig}\n")
        sys.exit(1)
    elif not append_db:
        sys.stderr.write("-append flag is required, when --wholeDB flag is not provided.\n")
        sys.exit(1)

def decompress(archive_file, out_path=None):
    print(f"Extracting {archive_file}...")
    tar = tarfile.open(archive_file, "r:bz2")
    if out_path:
        tar.extractall(out_path)
    else:
        tar.extractall()
    tar.close()
    os.remove(archive_file)


def downloadWeek(date, append = True):
    ano_filename = "BioLiP_" + date + "_nr.txt"
    rec_filename = "receptor_" + date + "_nr.tar.bz2"
    lig_filename = "ligand_" + date + "_nr.tar.bz2"

    print(f"Downloading {ano_filename}...")
    urllib.request.urlretrieve(os.path.join(download_url_base, ano_filename),
                               os.path.join(os.curdir, folder_ano, ano_filename))
    with open(append_db, 'a') as ano_file:
        with open(os.path.join(os.curdir, folder_ano, ano_filename), 'r') as update_file:
            for line in update_file:
                ano_file.write(line)

    print(f"Downloading {rec_filename}...")
    urllib.request.urlretrieve(os.path.join(download_url_base, rec_filename),
                               os.path.join(os.curdir, folder_rec, rec_filename))
    decompress(os.path.join(os.curdir, folder_rec, rec_filename))

    print(f"Downloading {lig_filename}...")
    urllib.request.urlretrieve(os.path.join(download_url_base, lig_filename),
                               os.path.join(os.curdir, folder_lig, lig_filename))
    decompress(os.path.join(os.curdir, folder_lig, lig_filename))



## MAIN ##
if wholeDB:
   updateFile =  open("updates.txt", "w+")
   updateFile.write("File listing all weekly updates already performed:\n")
   updateFile.close()

   print(f"Downloading annotations before 2013-03-06...")
   urllib.request.urlretrieve(download_url_anno_before2013, os.path.join(os.curdir, folder_ano, "BioLiP_nr.tar.bz2"))
   decompress(os.path.join(os.curdir, folder_ano, "BioLiP_nr.tar.bz2"), folder_ano)
   append_db = "BioLiP_database_updated.txt"
   with open(append_db, 'w+') as ano_file:
       with open(os.path.join(os.curdir, folder_ano, "BioLiP_2013-03-6_nr.txt"), 'r') as update_file:
           for line in update_file:
               ano_file.write(line)

   print(f"Downloading receptors before 2013-03-06, this may take a while...")
   urllib.request.urlretrieve(download_url_rec_before2013, os.path.join(os.curdir, folder_rec, "receptor_nr.tar.bz2"))
   decompress(os.path.join(os.curdir, folder_rec, "receptor_nr.tar.bz2"))

   print(f"Downloading ligands before 2013-03-06, this may take a while...")
   urllib.request.urlretrieve(download_url_lig_before2013, os.path.join(os.curdir, folder_lig, "ligand_nr.tar.bz2"))
   decompress(os.path.join(os.curdir, folder_lig, "ligand_nr.tar.bz2"))

elif not os.path.exists("updates.txt"):
    sys.stderr.write("FATAL: File updates.txt missing\n")
    sys.exit(1)
else:
    # read performed updated weekly updates
    print("Parsing updates.txt...")
    with open("updates.txt", "r") as updateFile:
        for line in updateFile.readlines():
            updated_weeks.append(line.strip())


# find which updates need to be performed
weeks_to_download = []
with open("updates.txt", "a") as updateFile:
    urllib.request.urlretrieve(overview_url, os.path.join(os.curdir,"weekly.html"))
    with open("weekly.html", "r") as weekly:
        table_found = False

        for weeklyline in weekly.readlines():

            if weeklyline.startswith("<div align=\"justify\">"):
                table_found = True  # ignore header links
    
            if table_found:
                if weeklyline.startswith("<tr><td>"):
                    date = re.search(r'(\d+-\d+-\d+)', weeklyline)
                    # check if this week is already updated
                    if not date.group(1) in updated_weeks:
                        weeks_to_download.append(date.group(1))


# Sort dates from oldest date to newest date
weeks_to_download.sort(key=lambda date: datetime.strptime(date, '%Y-%m-%d'))

with open("updates.txt", "a") as updateFile:
    for week in weeks_to_download:
        print(f"Updating week {week}")
        downloadWeek(week)
        updateFile.write(f"{week}")
        updateFile.write("\n")
        updateFile.flush()
        
        
