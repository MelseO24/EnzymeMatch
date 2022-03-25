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

import fnmatch
import os
import sys

from parsers.parse_BioLipIterator import BioLipIterator
from utils import utilities
from utils import globalParameters
from utils.utilities import raiseError

fout = sys.stdout
util = utilities


def matchResidue(parm_general, parm_residue):
    """
    :type parm_general: parsers.parse_inputFile.MatchGlobalParameter
    :type parm_residue: parsers.parse_inputFile.MatchResidueParameter
    :param parm_residue: Object containing input parameters of residue mode
    matches according to the parsed parameter file in parse input
    """

    fout.write("Starting Matchresidue...\n")

    # create output directory
    out_main_path = os.path.join(util.globalParameters.workdir, parm_residue.MR_directory_out)
    if not os.path.exists(out_main_path):
        os.makedirs(out_main_path)

    # open summary file, write the specified allowed residues
    writer_out_summary = util.openFile(os.path.join(out_main_path, parm_residue.SummaryOut), "w+")
    util.printInputInFile(writer_out_summary, parm_general, parm_residue)
    writer_out_summary.write("Residue Matching:\n")
    for position, allowed_residues in enumerate(parm_general.allowed_residues, 1):
        writer_out_summary.write(f"Allowed residues at position {position}: ")
        for allowed_residue in allowed_residues:
            writer_out_summary.write(allowed_residue + " ")
        writer_out_summary.write("\n")

    # initiate biolip iterator
    biolipDB = BioLipIterator(parm_residue.get_DBfile())
    biolipDBparser = iter(biolipDB)

    # BindingSites & CatalyticSite in one loop
    BiS_writer_out = util.openFile(os.path.join(out_main_path, parm_residue.BiSOut), 'w+')
    BiS_matchcounter = 0
    BiS_writer_out.write(utilities.header)

    Cat_writer_out = util.openFile(os.path.join(out_main_path, parm_residue.CatOut), 'w+')
    Cat_matchcounter = 0
    Cat_writer_out.write(utilities.header)

    # residues_Check only contains list of positions, where the search has to be done
    # if the backbone interaction is turned off, all positions are taken.
    # If the backbone interaction is turned on,
    # all the positions with "H_don" and "H_acc" are not included in this list
    # because all amino acids have hydrogen donors and acceptors in their backbone.
    if parm_general.backbone_interaction == 0:
        residues_Check = parm_general.allowed_residues
        if len(residues_Check) == 0:
            raiseError("FATAL: No allowed amino acids provided for ResidueMatch\n", True)
    else:
        residues_Check = []
        for position, residueList in enumerate(parm_general.allowed_residues):
            if parm_general.aminoacid_inputs[position] not in ["H_don", "H_acc"]:
                residues_Check.append(residueList)

        if len(parm_general.aminoacid_inputs) == 0:
            raiseError("FATAL: No allowed amino acids provided for ResidueMatch\n", True)



    #Using iterator
    entry_count = 0
    for entry in biolipDBparser:
        if globalParameters.debug > 0:
            sys.stdout.write(f"Entry being analyzed: {entry['entry_name']}\n")
        # Counting the entries counted until now
        entry_count += 1
        BiS_matched_intpoints = parm_residue.number_of_interaction_points
        Cat_matched_intpoints = parm_residue.number_of_interaction_points
        # Saving the positions that couldn't be found (For the error message)
        Bis_Erroneous = []
        Cat_Erroneous = []
        if parm_residue.EC is not []:
            biolipEC = biolipDB.get_EC().strip(".").split(".")  # Strip out the points from the biolip entry in the end
            skipEntry = False
            if len(biolipEC) < len(parm_residue.EC):  # if EC given in the biolip is shorter than the EC from input, skip
                skipEntry = True
            else:  # subcategory-wise comparison of the EC entries
                for i in range(len(parm_residue.EC)):
                    if parm_residue.EC[i] != biolipEC[i]:
                        skipEntry = True
                        break
            if skipEntry:
                if globalParameters.debug > 0:
                    sys.stdout.write(
                        f"Skipped, unwanted EC value. Wanted: {parm_residue.EC}, Found: {biolipDB.get_EC()}\n")
                continue  # break out of loop if entry enzyme is not in the right class
        # the BiS_matched_intpoints and Cat_match_intpoints are initialized with the maximum number of interaction points.
        # per position, where no amino acid can be assigned, the point number decreases
        # and if it is less than the minimum number of allowed residues, the entry is removed

       #First: the BiS or Cat entry in the database should at least have min_num_of_interaction_points AA
        if len(biolipDB.get_BisResidues()) < parm_residue.min_num_of_interaction_points:
            BiS_matched_intpoints = -1
        if len(biolipDB.get_CatResidues()) < parm_residue.min_num_of_interaction_points:
            Cat_matched_intpoints = -1

        for checklist in residues_Check:
            if not any([fnmatch.filter(biolipDB.get_BisResidues(), AA) for AA in checklist]):
                Bis_Erroneous.append(checklist)
                BiS_matched_intpoints -= 1  # none of the AA in a checklist was matched for the binding site
            if not any([fnmatch.filter(biolipDB.get_CatResidues(), AA) for AA in checklist]):
                Cat_Erroneous.append(checklist)
                Cat_matched_intpoints -= 1  # none of the AA in a checklist was matched for the catalytic site
        if BiS_matched_intpoints < parm_residue.min_num_of_interaction_points:
            BiS_matched_intpoints = False
        if Cat_matched_intpoints < parm_residue.min_num_of_interaction_points:
            Cat_matched_intpoints = False
        if globalParameters.debug > 0:
            print(f"Sites of the entry: {list(set([Amino[0] for Amino in biolipDB.get_BisResidues()]))}")
            if not BiS_matched_intpoints:
                print(f"Erroneous Site for BiS match: {Bis_Erroneous}")
            if not Cat_matched_intpoints:
                print(f"Erroneous Site for Cat match: {Cat_Erroneous}")
            print(f"BiS matched: {BiS_matched_intpoints}\tCat matched: {Cat_matched_intpoints}\n")
        if BiS_matched_intpoints:
            BiS_writer_out.write(biolipDB.get_fullLine())
            BiS_matchcounter += 1
        if Cat_matched_intpoints:
            Cat_writer_out.write(biolipDB.get_fullLine())
            Cat_matchcounter += 1

    BiS_writer_out.close()
    Cat_writer_out.close()

    try:
        fout.write(f"{entry_count} entries were analyzed:\n")
        fout.write(f"All entries ({BiS_matchcounter}) which fit residue combination in binding site are"
                   " written in " + parm_residue.BiSOut + "\n")
        writer_out_summary.write(f"\nNumber of entries analyzed: {entry_count}\n")
    except NameError:
        entry_count = 0
        raiseError(f"WARNING: {entry_count} entries were analyzed:\n"
                   "database might be empty", False)

    # summary
    writer_out_summary.write(f"Number of entries in BindingSite: {BiS_matchcounter:d}\n")
    fout.write(f"All entries ({Cat_matchcounter}) which fit residue combination in catalytic site are "
               "written in " + parm_residue.CatOut + "\n")
    writer_out_summary.write(f"Number of entries in Catalytic site: {Cat_matchcounter:d}\n")

