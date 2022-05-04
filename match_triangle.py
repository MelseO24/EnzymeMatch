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

from itertools import product, combinations, repeat
from math import inf
import os, re, sys
import numpy as np
from utils.pdbutils import get_nrAtoms
import utils.pdbutils
import utils.globalParameters
import Bio.PDB as biopdb
from AQD import calcBindingsiteDistances
import utils.globalParameters
import utils.pdbutils
from parsers.parse_BioLipIterator import BioLipIterator
from utils import globalParameters
from utils.utilities import CentralAtomName, AccurateCentralAtomName, Acc_positions, raiseError, calculatedAngles
import utils.utilities as util
import multiprocessing as mp

# Suppress inconsistent pdb warnings
import warnings
from Bio.PDB.StructureBuilder import PDBConstructionWarning
warnings.simplefilter('ignore', PDBConstructionWarning)

fout = sys.stdout
pdbparser = biopdb.PDBParser()

def matchTriangle(parm_general, parm_triangle):
    """
    :type parm_general: parsers.parse_inputFile.MatchGlobalParameter
    :type parm_triangle: parsers.parse_inputFile.MatchTriangleParameter
    """
    if parm_triangle.download_format == "mmtf":
        try:
            import Bio.PDB.mmtf
        except:
            sys.stderr.write("FATAL: Could not import the mmtf-python package. Either change the format in the input file "
                             "to PDB, or install the mmtf-python package in your python environment. See 'requirements.txt' "
                             "for the required version.\n")
            sys.exit(1)

    # create output directory
    out_main_path = os.path.join(util.globalParameters.workdir, parm_triangle.MT_directory_out)
    if not os.path.exists(out_main_path):
        util.makeDirectory(out_main_path)

    #create download directory
    if parm_triangle.download_mode != 0:
        pdb_download_path = os.path.join(util.globalParameters.workdir, parm_triangle.download_dir)
        if os.path.exists(pdb_download_path):
            util.makeDirectory(pdb_download_path)

    # create output files
    output = util.openFile(os.path.join(out_main_path, parm_triangle.outFile), 'w+')
    summary_output = util.openFile(os.path.join(out_main_path, parm_triangle.summaryFile), 'w+')
    biolipDB = BioLipIterator(parm_triangle.input)
    biolipDBparser = iter(biolipDB)

    #Writing headers
    if parm_triangle.scoring == 0:
        output.write(utils.utilities.header)
    else:
        headerlist = list()
        headerlist.append(utils.utilities.header.split("\t")[0])
        headerlist.append("Matched_points_count")
        if parm_general.backbone_interaction != 0:
            headerlist.append("Backbone_interaction_count")
        headerlist.append("Score")
        headerlist.append(",".join(f"pos{i+1}" for i in range(parm_triangle.number_of_interaction_points)))
        headerlist += utils.utilities.header.split("\t")[1:]
        util.writeLines(output, headerlist)

    # calculate distances for ligand if conformational search is required, i.e. user didn't specify distances
    parm_triangle = calcBindingsiteDistances(parm_triangle, parm_general)
    fout.write("Starting MatchTriangle...\n")
    util.printInputInFile(summary_output, parm_general, parm_triangle)
    summary_output.write("Triangle Matching:\n")
    for position, allowed_residues in enumerate(parm_general.allowed_residues, 1):
        summary_output.write(f"Allowed residues at position {position}: ")
        for allowed_residue in allowed_residues:
            summary_output.write(allowed_residue + " ")
        summary_output.write("\n")
    summary_output.flush()

    if parm_triangle.numprocs == 1:
        fout.write("TriangleMatch runs in single-core performance.\n")
        match_counter = 0
        for entry_counter, entry in enumerate(biolipDBparser):
            result = is_trianglematch(entry, parm_general, parm_triangle, None)
            if result:
                match_counter += 1
                util.writeLines(output, result)
                output.flush()
        output.close()
    else:
        # Set up multiprocessing
        output.close()
        manager = mp.Manager()
        queue = manager.Queue()
        pool = mp.Pool(parm_triangle.numprocs + 1)
        fout.write(f"TriangleMatch runs on {parm_triangle.numprocs} processes/CPU cores.\n")
        pool.apply_async(listener, (queue,os.path.join(out_main_path, parm_triangle.outFile)))  # listener thread with writing rights to write matched entries
        results = pool.starmap(is_trianglematch, zip(biolipDBparser,
                                           repeat(parm_general), repeat(parm_triangle), repeat(queue)))
        queue.put("kill")
        pool.close()
        entry_counter = len(results) -1
        match_counter = results.count(True)

    try:
        fout.write(f"{entry_counter} entries were analyzed by TriangleMatch.\n")
        fout.write(f"{match_counter} entries matched to the input query.\n")
        fout.write(f"All matched entries are saved in {parm_triangle.outFile}\n")
        summary_output.write(f"\nNumber of entries analyzed: {entry_counter}\n")
    except NameError:
        entry_counter = 0
        raiseError(f"WARNING: {entry_counter} entries were analyzed:\n"
                   "database might be empty", False)

    # summary
    summary_output.write(f"Number of matched entries: {match_counter}\n")
    fout.write("Triangle matching done...\n")

def listener(queue, fileName):
    """
    Listenes to queue (i.e. is_trianglematch() output) and writes to file when something adds queue
    :param queue: line to print (type: list)
    """
    entry_counter = 0
    output = util.openFile(fileName, 'a')
    while 1:
        message = queue.get()
        if message == "kill":
            output.close()
            break
        entry_counter += 1
        util.writeLines(output, message)
        output.flush()

def is_trianglematch(entry, parm_general, parm_triangle, queue):
    """
    Main function to check if the binding site analyzed matched the query interaction point graph
    :param entry: entry from BioLipPDB
    :param queue: queue to put output to in case of match, which should be written by a listener thread
    :type entry: dict
    :type parm_general: parsers.parse_inputFile.MatchGlobalParameter
    :type parm_triangle: parsers.parse_inputFile.MatchTriangleParameter
    :return: False if no match, otherwise line to print to outputfile (type: list)
    """

    # set some initials
    # Verbose if the verbosity is not turned off
    verbosity = util.globalParameters.verbosity > 1

    if parm_triangle.download_mode > 0:
        pdbDownloader = biopdb.PDBList(verbose=verbosity, obsolete_pdb="abc")

    # now start
    if entry["fullLine"].startswith(utils.utilities.header):
        return False
    if globalParameters.debug > 0:
        print(f"Entry being analyzed: {entry['entry_name']}")
    ligand_path = os.path.join(parm_triangle.ligand_folder, entry["LigPDBName"])
    receptor_path = os.path.join(parm_triangle.receptor_folder, entry["ProtPDBName"])
    pdb_protein = get_protein_structure(ligand_path, receptor_path, parm_triangle, entry,
                                        pdbDownloader if parm_triangle.download_mode > 0 else None)
    if pdb_protein is None:
        if globalParameters.debug > 0:
            sys.stderr.write(f"PDB not found, or problem while parsing: {entry['ProtPDBName']}\n")
        return False
    # get the chain of the residues
    chain = entry["protchain"]

    # create list of binding site residues (BiS_residues: ['F18, 'G23', ..])
    if parm_triangle.match_side == "bis":
        BiS_residues = list(entry["BisResidues"])
    else:
        BiS_residues = list(entry["CatResidues"])

    # create list of binding site residues per vertex
    # BiS_residues_at_vertex = [[D2_CG, G523_CA, G548_CA],[G7_CA, G18_CA], ...] (with list index as vertex nr)
    BiS_residues_at_vertex = potentialResiduesAtPositions(BiS_residues, parm_triangle, parm_general)

    # find whether there residues fit within the query distances
    residueDistances = residuesWithinQueryDistances(parm_general, parm_triangle, BiS_residues_at_vertex, pdb_protein, chain, entry)
    # Either boolean or an int, it wouldn't be None, if it returned right scores.
    if residueDistances:
        # match_counter += 1
        if parm_triangle.scoring == 0:
            outputlist = entry["fullLine"].split("\t")
        else:
            outputlist = list()
            # The entry name
            outputlist.append(entry["fullLine"].split("\t")[0])
            # The number of positions in the triangle
            outputlist.append(str(len(residueDistances[1])))
            if parm_general.backbone_interaction == 1:
                # ...of them the backbone interaction points
                outputlist.append(str(sum(1 for i in residueDistances[1] if i[1].split("_")[1] in ["N", "O"])))
            # The score of this entry
            outputlist.append(str(round(residueDistances[0],3)))
            mappingPositions = dict()
            for pos, res in residueDistances[1]:
                mappingPositions[int(pos)] = res
            interaction_point_output = ""
            for i in range(parm_triangle.number_of_interaction_points):
                if i in mappingPositions:
                    interaction_point_output += mappingPositions[i]
                interaction_point_output += ","
            outputlist.append(interaction_point_output[:len(interaction_point_output)-1])
            outputlist.append("\t".join(entry["fullLine"].split("\t")[1:]))
        if queue:
            queue.put(outputlist)
            return True
        else:
            return outputlist
    return False

####################################3
def get_protein_structure(ligand_path, receptor_path, parm_triangle, entry, pdbDownloader=None,):
    """
    Checking if the ligand file and the protein file exist. If the download parameter is turned on, the file will be downloaded too
    :type entry: dict
    :param ligand_path: path of the ligand folder
    :param receptor_path: path of the receptor folder
    :param parm_triangle: triangle parameters
    :param entry: entry from biolipiterator
    :param pdbDownloader: pdb downloader
    :return:protein structure file. None is returned if the protein structure is inappropriate or not suitable
    """
    if parm_triangle.download_format == "mmtf":
        mmtfparser = biopdb.mmtf.MMTFParser
    # if one of the biolip files aren't there: download pdb from pdb database
    if not os.path.isfile(ligand_path) or not os.path.isfile(receptor_path):
        # returning error message for ligand file, if not given
        if not os.path.isfile(ligand_path):
            if util.globalParameters.verbosity > 1:
                sys.stderr.write(
                    f"WARNING: Ligand file {os.path.join(parm_triangle.ligand_folder, entry['LigPDBName'])} does not exist.\n")
            if utils.globalParameters.strict_mode:
                sys.stderr.write(f"Program terminated, since strict mode is switched on, and Ligand file"
                                 f" {os.path.join(parm_triangle.ligand_folder, entry['LigPDBName'])} does not exist "
                                 f"Either switch download function on, or download the respective PDB file manually\n")
                sys.exit(1)
        # returning error message for receptor file, if not given
        if not os.path.isfile(receptor_path):
            if util.globalParameters.verbosity > 1:
                sys.stderr.write(
                    f"WARNING: Receptor file {os.path.join(parm_triangle.receptor_folder, entry['ProtPDBName'])} does not exist.\n")
            if utils.globalParameters.strict_mode:
                sys.stderr.write(f"Program terminated, since strict mode is switched on, and Receptor file"
                                 f" {os.path.join(parm_triangle.receptor_folder, entry['ProtPDBName'])} does not exist. "
                                 f"Either switch download function on, or download the respective PDB file manually.\n")
                sys.exit(1)
        # skipping if download function is turned off
        if parm_triangle.download_mode == 0:
            if util.globalParameters.verbosity > 1:
                sys.stderr.write(f"This entry: '{entry['entry_name']}' will be skipped.\n")
            return None
        else:
            # Downloading file
            if util.globalParameters.verbosity > 1:
                sys.stdout.write(f"This entry: '{entry['entry_name']}' will be downloaded.\n")
            try:
                # The download path for the download file
                download_path = os.path.join(util.globalParameters.workdir, parm_triangle.download_dir)
                # Won't download multiple times, only downloading the file if it is not given
                if not os.path.isfile(os.path.join(download_path, entry['pdbDownloadName'] + ".pdb")) and parm_triangle.download_format == "pdb":
                    # Downloading the pdbfile from pdb database
                    pdbDownloader.retrieve_pdb_file(entry['pdbDownloadName'], pdir=download_path, file_format="pdb",
                                                    overwrite=True)
                    # Renaming the pdb file, because it is downloaded in .ent format
                    os.rename(os.path.join(download_path, "pdb" + entry['pdbDownloadName'] + ".ent"),
                              os.path.join(download_path, entry['pdbDownloadName'] + ".pdb"))
                    # getting protein structure
                    pdb_protein = pdbparser.get_structure("protein", os.path.join(download_path,
                                                                                  entry['pdbDownloadName'] + ".pdb"))
                # mmtf file option also checks for the .pdb files
                elif not os.path.isfile(os.path.join(download_path, entry['pdbDownloadName'] + ".mmtf")) and \
                          not os.path.isfile(os.path.join(download_path, entry['pdbDownloadName'] + ".pdb")) and\
                                parm_triangle.download_format == "mmtf":
                    # Downloading the mmtf file form pdb database
                    pdbDownloader.retrieve_pdb_file(entry['pdbDownloadName'], pdir=download_path,
                                                    file_format="mmtf", overwrite=True)
                    # getting protein structure
                    pdb_protein = mmtfparser.get_structure(
                        os.path.join(download_path, entry['pdbDownloadName'] + ".mmtf"))
                else:
                    # if the pdb file is found, it's read, if not, the mmtf file is read
                    if os.path.isfile(os.path.join(download_path, entry['pdbDownloadName'] + ".pdb")):
                        # if the file is found
                        if util.globalParameters.verbosity > 1:
                            sys.stdout.write(
                                f"Found downloaded pdb file for entry {entry['entry_name']}, using this file...\n")
                        pdb_protein = pdbparser.get_structure("protein", os.path.join(download_path,
                                                                                      entry['pdbDownloadName'] + ".pdb"))
                    else:
                        # if the file is found
                        if util.globalParameters.verbosity > 1:
                            sys.stdout.write(
                                f"Found downloaded mmtf file for entry {entry['entry_name']}, using this file...\n")
                        pdb_protein = mmtfparser.get_structure(os.path.join(download_path, entry['pdbDownloadName'] + ".mmtf"))
            except:
                # errormessage for failed downloads
                sys.stderr.write(f"WARNING: downloading file {entry['pdbDownloadName']} failed.\n"
                                 f"This entry: '{entry['entry_name']}' will be skipped.\n")
                if utils.globalParameters.strict_mode:
                    sys.stderr.write(f"Program terminated, since strict mode is switched on, and downloading file"
                                     f" {entry['pdbDownloadName']} failed.\n")
                    sys.exit(1)
                return None
            # Reference to the biolip paper https://zhanglab.dcmb.med.umich.edu/papers/2013_1.pdf
            # Page 2, Materials and methods, step 2:
            # "To avoid conflicts with the existing three-letter code, we name k-mer, DNA/RNA and peptide ligands as UUU, NUC and III, respectively."
            # Skipped anyways, because a multiple ligand ids imply that there are more than 2 residues anyway.
            try:
                if entry['LigID'] not in ['UUU','NUC','III']:
                    # Searching for ligand sizes
                    # filtering out the chain
                    receptor_Chain = None
                    ligand_Chain = None
                    # filtering out the residue in order to check for its size
                    Wanted_Residue = None
                    # Getting the wanted chain
                    for chain in pdb_protein.get_chains():
                        if chain.get_id() == entry['protchain']:
                            receptor_Chain = chain
                        if chain.get_id() == entry['LigChain']:
                            ligand_Chain = chain
                    # break if the wanted chain is not present
                    if receptor_Chain is None:
                        if util.globalParameters.verbosity > 1:
                            sys.stderr.write(
                                f"WARNING: No {entry['protchain']} chain found for file {entry['entry_name']}\n"
                                f"This entry: '{entry['pdbDownloadName']}' will be skipped\n")
                        if utils.globalParameters.strict_mode:
                            sys.stderr.write(f"Program terminated, since strict mode is switched on, and no"
                                             f" {entry['protchain']} chain found for file {entry['pdbDownloadName']}.\n")
                            sys.exit(1)
                        return None
                    if ligand_Chain is None:
                        if util.globalParameters.verbosity > 1:
                            sys.stderr.write(
                                f"WARNING: No {entry['LigChain']} chain found for file {entry['pdbDownloadName']}\n"
                                f"This entry: '{entry['pdbDownloadName']}' will be skipped\n")
                        if utils.globalParameters.strict_mode:
                            sys.stderr.write(f"Program terminated, since strict mode is switched on, and no"
                                             f" {entry['LigChain']} chain found for file {entry['pdbDownloadName']}.\n")
                            sys.exit(1)
                        return None
                    # Get the wanted residue
                    for residues in ligand_Chain.get_residues():
                        if residues.get_resname().strip() == entry['LigID']:
                            Wanted_Residue = residues
                            break
                    if Wanted_Residue is None:
                        if util.globalParameters.verbosity > 1:
                            sys.stderr.write(
                                f"WARNING: No {entry['LigID']} ligand found for file {entry['entry_name']}\n"
                                f"This entry: '{entry['entry_name']}' will be skipped\n")
                        if utils.globalParameters.strict_mode:
                            sys.stderr.write(f"Program terminated, since strict mode is switched on, and "
                                             f"no {entry['LigID']} ligand found for file {entry['entry_name']}.\n")
                            sys.exit(1)
                        return None
                    # Break if the ligand size is too small
                    if len(Wanted_Residue.get_unpacked_list()) < 4:
                        return None
            except:
                # errormessage for failed parsing
                sys.stderr.write(f"WARNING: parsing file {entry['pdbDownloadName']} failed.\n"
                                 f"This entry: '{entry['entry_name']}' will be skipped.\n")
                if utils.globalParameters.strict_mode:
                    sys.stderr.write(f"Program terminated, since strict mode is switched on, and "
                                     f"parsing file {entry['pdbDownloadName']} failed.\n")
                    sys.exit(1)
                return None
    else:
        pdb_ligand = pdbparser.get_structure("ligand", os.path.join(parm_triangle.ligand_folder, entry['LigPDBName']))
        # Only consider binding sites of ligands >3 atoms
        if get_nrAtoms(pdb_ligand) < 4:
            return None
        pdb_protein = pdbparser.get_structure("protein",
                                              os.path.join(parm_triangle.receptor_folder, entry['ProtPDBName']))
    return pdb_protein

def potentialResiduesAtPositions(enzyme_residues, param_triangle, param_general):
    """
    checks for all residues if they are allowed at each position, and returns a list with the allowed
    residues at each position.
    list of residues for each vertex position: [['N39','D60'],['D60','A19','G23'],...]
    :type enzyme_residues: List[str,str,..]
    :type param_triangle: parsers.parse_inputFile.MatchTriangleParameter
    :type param_general: parsers.parse_inputFile.MatchGlobalParameter
    :return list with allowed matching residues in each position
    :rtype: List[List[str,str,..],..]
    """
    # initiating a list of empty lists, saving possible positions for every interaction points
    residues_at_vertex = [[] for _ in range(param_triangle.number_of_interaction_points)]

    # the positions, lists of the possible residues are iterated
    # If the backbone interaction is turned on, the "N" locations and the "O" locations of ALL residues
    # in enzyme_residues for positions that had "H_don" and "H_acc" as query as they are hydrogen donors/acceptors
    for position, allowed_residues in enumerate(param_general.allowed_residues):
        # Here, only the backbone interaction points are taken into account
        if param_general.aminoacid_inputs[position] in ["H_don", "H_acc"] and param_general.backbone_interaction == 1:
            for residue in enzyme_residues:
                if param_general.aminoacid_inputs[position] == "H_acc":
                    residues_at_vertex[position].append(residue + "_" + "O")
                else:
                    residues_at_vertex[position].append(residue + "_" + "N")
        # Here, only the non-backbone interaction points are taken into account
        for allowed_residue in allowed_residues:
            for residue in enzyme_residues:
                # binding/catalytic site (enzyme) residue should be an allowed residue by the given parameters
                if allowed_residue.strip('*') in re.findall(r'\D', residue):
                    # Depending on the accurate positions parameter, different dictionary is used, which is also read differently.
                    # Look inside utils.utilities.py for the dictionary
                    if param_triangle.accurate_mode == 0: # less accurate positions
                        residues_at_vertex[position].append(residue + "_" + util.CentralAtomName[residue[0]])
                    else: # more accurate positions
                        for element in util.AccurateCentralAtomName[residue[0]]:
                            residues_at_vertex[position].append(residue + "_" + element)

    return residues_at_vertex


# triangle matching measurement functions:
def calcMeanAndVariance(v, vertex_labels=None, parm_triangle=None):
    """
    Calculates mean and variance of input list, if vertex_labels given, then an offset of 2A is added to the variance
    for every non H-acc and H-don in the labels (since they are centered on the ligand atom rather than the proten)
    :param v: list of values to calculate variance and mean
    :param vertex_labels:the two positions between which the distances were measured (set of 2 ints)
    :type parm_triangle: parsers.parse_inputFile.MatchTriangleParameter
    :return: list of mean and variance from input list
    """
    var = np.var(v)
    mean = np.mean(v)
    if vertex_labels:
        if not parm_triangle:
            sys.stderr.write("Internal error: alcMeanAndVariance needs parm_triangle, if vertex_labels are given.\n")
            sys.exit(1)
        ##add offset 2 Angstrom for each non-Hacc/Hadon entry (since center is on ligand, not on protein)
        for position in vertex_labels:
            if parm_triangle.auto_query_types[position] not in ["H_acc", "H_don"]:
                if globalParameters.debug > 0:
                    fout.write(f"Offset for distance {vertex_labels} was increased with 2, because position {position} "
                               f"is {parm_triangle.auto_query_types[position]}.\n")
                var += 2
    return mean, var

def residuesWithinQueryDistances(param_general, param_triangle, residues_at_vertex, protein_pdb, chain, entry):
    """
    creates combinations of residues at each vertex, which are essentially possible triangles (polygons)
    and finds whether these potential triangle are within the query distances
    :type param_triangle: parsers.parse_inputFile.MatchTriangleParameter
    :param residues_at_vertex: [[D36,D39,D55],[D36,D39,D55,N123],[G81,G129]]
    :type protein_pdb: Bio.PDB.Structure.Structure
    :param chain: protein chain
    :return: float or boolean.
    If the scoring function is turned on, the full score is returned, when it's off, it only returns a boolean value
    """
    # dictionary of already calculated distances. Will be updated in every loop.
    distance_dict = {}
    max_match = param_triangle.number_of_interaction_points
    min_match = param_triangle.min_num_of_interaction_points
    score = param_triangle.scoring
    # saving the residues with enumerated positions
    positionVertex = list(enumerate(residues_at_vertex))
    if score != 0:
        # Iterates over the possible number interaction points
        for pos_num in range(max_match, min_match-1, -1):
            best_score = getTriangleScore(param_general, param_triangle, protein_pdb, chain, entry, positionVertex,
                                          pos_num, distance_dict)
            if best_score:
                return best_score
    else:
        # Only smallest size need to be calculated when the score is turned off
        triangleScore = getTriangleScore(param_general, param_triangle, protein_pdb, chain, entry, positionVertex,
                                         min_match, distance_dict)
        if triangleScore:
            return True
    return False

def getTriangleScore(param_general, param_triangle, protein_pdb, chain, entry, positionVertex, pos_num, distance_dict):
    # set to infinity, as ANY valid score would be under this value, and the smaller the score value, the better.
    best_score = inf
    best_triangle = None
    # Getting interaction points (list of all possible points in the position) of selected size
    # (e.g. (pos0,pos2,pos4), when one chooses 3 points out of 5)
    for positions in combinations(positionVertex, pos_num):
        # The position numbers have to be saved, as the position information should determine if the distance between
        # two points are valid or not
        positionNumbers = [x[0] for x in positions]
        positions = [x[1] for x in positions]
        # Getting actual triangles out of selected interaction points
        for triangle in product(*positions):
            # duplicate check
            if len(list(set(triangle))) != pos_num:
                continue
            else:
                triangle = [(positionNumbers[i], triangle[i]) for i in range(pos_num)]
                triangleScore = triangleWithinQueryDistances(triangle, param_triangle, protein_pdb, chain,
                                                             entry, distance_dict)
                if triangleScore:
                    if param_triangle.scoring != 0:
                        if param_general.backbone_interaction == 1:
                            for residue in triangle:
                                if residue[1].split("_")[1] in ["N", "O"]:
                                    triangleScore += param_triangle.backbone_penalty # adding backbone penalty
                        if best_score > triangleScore: # Best score is updated upon finding better scores
                            best_score = triangleScore
                            best_triangle = triangle
                    else:
                        return True
    if best_score == inf:
        return False
    else:
        return best_score, best_triangle

def triangleWithinQueryDistances(triangle, param_triangle, protein_pdb, chain, entry, distance_dict):
    """
    calculates distances of / calls isInRange on all vertex combinations of a potential triangle
    :param entry: the entry name
    :param distance_dict: the dictionary of already calculated distances
    :param triangle: list of residues of a potential triangle=[D20, D32, G80]
    :type protein_pdb: Bio.PDB.Structure.Structure
    :param chain: protein chain
    :type param_triangle: parsers.parse_inputFile.MatchTriangleParameter
    :return: float or boolean.
    If the scoring function is turned on, the full score is returned, when it's off, it only returns a boolean value
    """
    # enumerate(triangle)=[(0,D20_POS), (1,D32_POS), (2,G80_POS)]
    # combinations(...,2) = [[(0,D20_POS), (1,D32_POS)], [...], ...]
    # score_sum adds up all the distance scores calculated
    score_sum = 0
    _nr_edges = 0
    for residue_combinations in combinations(triangle, 2):
        if frozenset([residue_combinations[0], residue_combinations[1]]) in distance_dict:
            score = distance_dict[frozenset([residue_combinations[0], residue_combinations[1]])]
        else:
            score = isResiduePairInRange(residue_combinations, param_triangle, protein_pdb, chain, entry)
            distance_dict[frozenset([residue_combinations[0], residue_combinations[1]])] = score
        if not score:
            if globalParameters.debug > 1:
                print(f"residues which are checked on their inter-distance: {triangle}")
            return False
        else:
            score_sum += score
            _nr_edges += 1
    if globalParameters.debug > 1:
        print(f"residues which are checked on their inter-distance: {triangle}")
    # all the distance scores aren't None. This triangle is valid
    if param_triangle.scoring == 0:
        return True
    else:
        return score_sum/_nr_edges

def isResiduePairInRange(residue_combinations, param_triangle, protein_pdb, chain, entry):
    """
    calculates for the given residue combination whether they are close enough
    :param residue_combinations: combinations of two vertices and their residue e.g. [(0,D20), (1,D32)]
    :type param_triangle: parsers.parse_inputFile.MatchTriangleParameter
    :type protein_pdb: Bio.PDB.Structure.Structure
    :param chain: protein chain
    :param entry: protein entry name
    :return: boolean or float, containing the score
    """
    residue1 = residue_combinations[0][1]
    residue2 = residue_combinations[1][1]
    vertex_combination = frozenset((residue_combinations[0][0], residue_combinations[1][0]))
    # calculating the allowed distance between positions, given from the input
    allowed_distance = param_triangle.distances[vertex_combination]
    allowed_offset = param_triangle.offset[vertex_combination]

    try:
        real_distance = getResDistance(residue1, residue2, protein_pdb, chain, entry, param_triangle.accurate_mode)
    except KeyError as error:
        key = error.args[0]
        if (param_triangle.accurate_mode == 0 and key in CentralAtomName.values()) or \
                (param_triangle.accurate_mode == 1 and key in AccurateCentralAtomName.values()) or \
                (param_triangle.accurate_mode == 1 and key in Acc_positions) or \
                        key in ['N', 'O', 'rotation', 'midpoint']:
            if globalParameters.verbosity > 1:
                raiseError(f"WARNING: Triangle matching could not be performed in this entry, probably due to "
                           f" non-resolved residue or side chain:"
                           f" entry: {entry['entry_name']} res1: {residue1}, res2: {residue2}.\n", False)
        else:
            raise KeyError

        real_distance = inf
    if globalParameters.debug > 1:
        print(f"[{allowed_distance - allowed_offset} < {real_distance} < {allowed_distance+allowed_offset}]")
        print(param_triangle.offset)
    if allowed_distance - allowed_offset < real_distance < allowed_distance + allowed_offset:
        # This distance is valid, if the scoring is not needed, only a boolean is returned
        if param_triangle.scoring == 0:
            return True
        else:
            return abs(allowed_distance-real_distance)
    return False

def getResDistance(res1, res2, pdb_protein, chain, entry, accurate):
    """
    Matches 'representative/central' atomnames and
    caculates distance via pdbReader module
    :type entry: str
    :type pdb_protein: Bio.PDB.Structure.Structure
    :type chain: str
    :type res1: str
    :type res2: str
    """

    #split the residue ID and atom name
    res_pos1 = res1.split("_")
    res_pos2 = res2.split("_")

    # get the positions of the residues: re for integers
    residue1_nr = int(re.findall(r'-?[0-9]+', res_pos1[0])[0])
    residue2_nr = int(re.findall(r'-?[0-9]+', res_pos2[0])[0])

    # get central atom and if applicable, insertion code: re for alphabet
    central_atom_insertion_code1 = re.findall(r'[a-zA-Z]', res_pos1[0])
    central_atom_insertion_code2 = re.findall(r'[a-zA-Z]', res_pos2[0])

    # inspection check for insertion code of residue
    if len(central_atom_insertion_code1) == 1:
        insertion1 = " "
    else:
        try:
            insertion1 = central_atom_insertion_code1[1]
        except:
            raiseError(f"Likely BioLiP entry {entry['entry_name']} has a problem.\n", False)
            if globalParameters.strict_mode:
                sys.exit(1)
            return inf

    if len(central_atom_insertion_code2) == 1:
        insertion2 = " "
    else:
        try:
            insertion2 = central_atom_insertion_code2[1]
        except:
            raiseError(f"Likely BioLiP entry {entry['entry_name']} has a problem.\n", False)
            if globalParameters.strict_mode:
                sys.exit(1)
            return inf

    if not accurate:
        # Calls the Atom Object
        atom1 = pdb_protein[0][chain][(" ", residue1_nr, insertion1)][res_pos1[1]]
        atom2 = pdb_protein[0][chain][(" ", residue2_nr, insertion2)][res_pos2[1]]
        #if globalParameters.debug > 0:
        #    print(f"receptor pdb: {biolipDB.ProtPDBName(entryName)}")
        #    print(f"residue: {res1}")
        #    print(f"Residue {residue1_nr}")
        #    print(f"chain {chain}")
        #    print(f"central atom1 {central_atom1}")

        # Returns the euclidean distance
        return atom1-atom2
    else:
        try:
            positionOne = getVectorOfPositions(res_pos1, pdb_protein, chain, residue1_nr, insertion1)
        except:  # As it can be seen above, the inf score is treated as a void score. It means that the distance is invalid
            return inf
        try:
            positionTwo = getVectorOfPositions(res_pos2, pdb_protein, chain, residue2_nr, insertion2)
        except:  # As it can be seen above, the inf score is treated as a void score. It means that the distance is invalid
            return inf
        # For a rotation-rotation distance calculation, one needs to find a minimum of all possible n * m combinations
        # Which requires a double-loop
        # For a roataion-point distance calculation, one only needs to find a minimum of n combinations
        # For a point point distance calculation, just the euclidean distance is returned
        # If a position is a list, it means it was a rotation and need comparing between every possible points.
        if not type(positionOne) is list and not type(positionTwo) is list:
            return (positionTwo-positionOne).norm()
        elif not type(positionOne) is list:
            return min([(positionOne-positionTwoElement).norm() for positionTwoElement in positionTwo])
        elif not type(positionTwo) is list:
            return min([(positionTwo - positionOneElement).norm() for positionOneElement in positionOne])
        else: # Double loop, comparing every pairs with each other, calculating the minimum distance
            return min([(positionOneElement - positionTwoElement).norm() for positionTwoElement in positionTwo for positionOneElement in positionOne])

def getVectorOfPositions(residue_position, pdb_protein, chain, residue1_nr, insertion):
    """
    Calculating positions of the possible interaction points when the accurate postions are turned on.
    :param residue_position: the position of the residue
    :param pdb_protein: the structure file of the protein
    :param chain: chain identifier
    :param residue1_nr: residue number
    :param insertion: if an insertion is given
    :return: if rotation, list of points in 3D space, else, a point in 3D space.
    """
    # For a rotation, Rodrigues' rotation formula is used.
    # A rotation axis, k, is calculated from the base of the rotating residue and the chain
    # The rotated vector, v, is calculated from the base of the rotating residue to one of the tips
    # The function is then v_rot = v * cos(x) + (k x v) * sin(x) + k * (k * v) * (1 - cos(x))
    # The first term scaled the vector down on a plane, second skews it towards the new rotational position
    # the third re-adds the height. We assign various angles in x and return a list
    # Graphical explination in the manuals
    if residue_position[1] == "rotation":
        point1 = pdb_protein[0][chain][(" ", residue1_nr, insertion)][residue_position[2]].get_vector()
        point2 = pdb_protein[0][chain][(" ", residue1_nr, insertion)][residue_position[3]].get_vector()
        point3 = pdb_protein[0][chain][(" ", residue1_nr, insertion)][residue_position[4]].get_vector()
        kVector = (point2 - point1).normalized()
        vVector = point3 - point2
        """
        EVALUATION SCRIPT! WILL BE DELETED BEFORE MERGE
        out = [point2 + ((vVector ** np.cos(angle)) + ((kVector ** vVector) ** np.sin(angle)) + ((kVector ** (kVector * vVector)) ** (1 - np.cos(angle)))) for angle in
                utils.utilities.calculatedAngles]
        for elementt in out:
            posSaver = elementt.get_array()
            gg = []
            for poss in posSaver:
                gg.append(str(round(poss, 3)).rjust(8))
            print(f"ATOM  {str(globalParameters.counternumber).rjust(5)}  N   DUM A{str(globalParameters.counternumber).rjust(4)}    {gg[0]}{gg[1]}{gg[2]}")
            globalParameters.counternumber += 1
        return out
        """
        return [point2 + ((vVector ** np.cos(angle)) + ((kVector ** vVector) ** np.sin(angle)) + ((kVector ** (kVector * vVector)) ** (1 - np.cos(angle)))) for angle in
                calculatedAngles]
    # For a midpoint calculation, we necessarily add all the points and average it, getting a position in the middle
    elif residue_position[1] == "midpoint":
        sumPoints = pdb_protein[0][chain][(" ", residue1_nr, insertion)][residue_position[2]].get_vector()
        for point in residue_position[3:]:
            sumPoints += pdb_protein[0][chain][(" ", residue1_nr, insertion)][point].get_vector()
        """
        EVALUATION SCRIPT! WILL BE DELETED BEFORE MERGE
        out = sumPoints / (len(residue_position) - 2)
        posSaver = out.get_array()
        gg = []
        for poss in posSaver:
            gg.append(str(round(poss, 3)).rjust(8))
        print(f"ATOM  {str(globalParameters.counternumber).rjust(5)}  N   DUM A{str(globalParameters.counternumber).rjust(4)}    {gg[0]}{gg[1]}{gg[2]}")
        globalParameters.counternumber += 1
        return out
        """
        return sumPoints / (len(residue_position) - 2)
    else:
        return pdb_protein[0][chain][(" ", residue1_nr, insertion)][residue_position[1]].get_vector()

