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

import sys, re, os
from utils import globalParameters


fout = sys.stdout
ferr = sys.stderr

CentralAtomName = {
    "R": "CZ",  # Arginine
    "K": "NZ",  # Lysine
    "H": "CE1",  # Histidine
    "C": "SG",  # Cysteine
    "Y": "OH",  # Tyrosine
    "S": "OG",  # Serine
    "G": "CA",  # Glycine
    "N": "CG",  # Asparagine
    "D": "CG",  # Aspartate
    "Q": "CD",  # Glutamine
    "E": "CD",  # Glutamate
    "T": "OG1",  # Threonine
    "U": "SEG",  # Selenocysteine: couldnt find anything
    "P": "CG",  # Proline
    "A": "CB",  # Alanine
    "L": "CG",  # Leucine
    "I": "CB",  # Isoleucine
    "M": "SD",  # Methionine
    "F": "CG",  # Phenylalanine
    "W": "NE1",  # Tryptophan
    "V": "CB",  # Valine
}

# Dictionary used when accurate positions are turned on.
# rotation: Residues that rotate. "rotation_X1_X2_X3" means: X1 and X2 to calculate axis vector, X2 and X3 to calculate vector to rotate
# midpoint: getting midpoints from two or three points in a ring. Two is needed in a hexagon ring, three in pentagon
AccurateCentralAtomName = {
    "R": ["NE", "rotation_NE_CZ_NH2"],  # Arginine
    "K": ["NZ"],  # Lysine
    "H": ["midpoint_CG_ND1_CE1_NE2_CD2"],  # Histidine
    "C": ["SG"],  # Cysteine
    "Y": ["OH", "midpoint_CG_CD1_CE1_CZ_CE2_CD2"],  # Tyrosine
    "S": ["OG"],  # Serine
    "G": ["CA"],  # Glycine
    "N": ["rotation_CB_CG_OD1"],  # Asparagine
    "D": ["rotation_CB_CG_OD1"],  # Aspartate
    "Q": ["rotation_CG_CD_OE1"],  # Glutamine
    "E": ["rotation_CG_CD_OE1"],  # Glutamate
    "T": ["OG1"],  # Threonine
    "U": ["SEG"],  # Selenocysteine: couldnt find anything
    "P": ["CG"],  # Proline
    "A": ["CB"],  # Alanine
    "L": ["rotation_CB_CG_CD1"],  # Leucine
    "I": ["CG2", "CG1", "CD1"],  # Isoleucine
    "M": ["SD"],  # Methionine
    "F": ["midpoint_CG_CD1_CE1_CZ_CE2_CD2"],  # Phenylalanine
    "W": ["midpoint_CG_CD1_NE1_CE2_CD2", "midpoint_CD2_CE2_CZ2_CH2_CZ3_CE3"],  # Tryptophan
    "V": ["rotation_CA_CB_CG1"],  # Valine
}

# List of all the valid positions. Is used to check if the entries are valid in MatchTriangle
# In matchTriangle, the program throws error messages, containing the erroneous sites.
# To be able to catch such sites, all the positions that are taken into account are saved
backbone_positions = ['N', 'O']
Acc_positions = ['CA', 'CB', 'CD', 'CD1', 'CD2', 'CE1', 'CE2', 'CE3', 'CG', 'CH2', 'CG1', 'CG2', 'CZ', 'CZ2', 'CZ3',
                 'ND1', 'NE', 'NE1', 'NE2', 'NH2', 'NZ', 'OD1', 'OE1', 'OG', 'OG1', 'OH', 'SD', 'SEG', 'SG']


# All the angles in 5 degrees step, change the step value in order to be more exact/ less exact
calculatedAngles = [i for i in range(0, 359, 5)]

header = "PDB_ID\tPDB_chain\tResolution\tBinding_site_number_code\tLigand_ID\tLigand_chain\tLigand_serial_number\tBiS_residues\tBis_residues\tCat_residues\tCat_residues\tEC\tGO\tBinding_affinity_literature\tBinding_affinity_MOAD\tBinding_affinity_PDBbind-CN\tBinding_affinity_BindingDB\tUniprotID\tPubMedID\tReceptor_sequence\n"

# File handling functions
def openFile(path, mode):
    """
    Uses Python's open function to open a file, in the mode provided.
    First checks if file exists via utilities.doesFileExists(), and terminiates if file exists and no overwrite allowed
    :param path: str
    :param mode: a mode accepted by Python's open function, eg. w+, r, a, a+ etc.
    :return: writer object
    """
    if "r" in mode:
        if not _doesFileExist(path):
            raise FileNotFoundError
    else:
        _isFileWritable(path,mode)
    try:
        return open(path,mode)
    except PermissionError:
         raiseError(f"FATAL: Insufficient permission to access or create {path}\n", True)

def makeDirectory(directory_path):
    try:
        if not os.path.exists(directory_path):
            os.makedirs(directory_path)
    except PermissionError:
        raiseError(f"FATAL: Insufficient permission to access or create {path}\n", True)

def _isFileWritable(path,mode):
    """
    IMPORTANT: Generally not necessary to communicate with directly, since this function ins intrinsically used by
    the utilities.openFile() function before returning the write object!

    Checks if file exists, terminates execution if file does exist and overwrite is switched off.
    NOTE: utilities function "openFile()" uses this function automatically to check if file exists!
    :param path: path to file to check if file exists
    """
    if mode == "r":
        if not globalParameters.overwrite and os.path.isfile(path):
            raiseError(
                f"WARNING: File {os.path.basename(path)} already exists.\nProvide -O flag to allow overwrite.\n", True)
    elif mode == "a":
        if not os.path.isfile(path):
            raiseError(f"WARNING: file {path} cannot be accessed.\n", True)

def writeLines(fileobj, elements):
    """
    Writing list of elements into a file
    :param fileobj: file to be written to. REQUIRES AN open(file) OBJECT
    :param elements: elements to be written
    :type elements: List[str, str]
    """
    fileobj.write('\t'.join([elem.strip('\n') for elem in elements]) + '\n')

def _doesFileExist(path):
    return os.path.isfile(path)


def printInputInFile(fileobject, parm_general, parm_mode_specific):
    """
    Prints relevant input settings in summary file.
    :type fileobject: file
    :param fileobject: fileobject which can be written to (created with open(xx))
    :param parm_general: object containing general options
    :param parm_mode_specific: object of either residue or triangle parameters
    :type parm_general: parsers.parse_inputFile.MatchGlobalParameter
    :type parm_mode_specific: parsers.parse_inputFile.MatchResidueParameter or parsers.parse_inputFile.MatchTriangleParameter
    """
    fileobject.write("Input parameters:\n\n")
    input_data = parm_general.allAttributes()
    fileobject.write("General parameters:\n")
    for attribute in input_data.keys():
        fileobject.write(f"{attribute} : {str(input_data[attribute])}\n")
    input_data = parm_mode_specific.allAttributes()
    fileobject.write("-----\n")
    fileobject.write("Mode specific parameters:\n")
    for attribute in input_data.keys():
        fileobject.write(f"{attribute} : {str(input_data[attribute])}\n")
    fileobject.write("\n")
    fileobject.flush()
    return 0


def raiseError(message, StrictMode):
    """
    :type message: str
    :type StrictMode: bool
    :param message: Message to write in stdout.err
    :param StrictMode: if it should exit or not
    :return:
    """
    sys.stderr.write(message)
    if StrictMode:
        exit(1)


# Utilities for Input parsing
def parseInputLine(line):
    """
    parses raw input line. i.e. removes all comments (after semicolon), and splits after first colon.
    :type line: str
    :param line: raw unstripped line from input file
    :return var1: input option, var2: setting
    :rtype: (str, str)
    """
    line_clean = line.split(";")[0].strip()  # remove comments
    input_option = line_clean.split(":")[0].strip()
    setting = re.sub('^.*?:', '', line_clean).strip()
    return input_option, setting


class NonResolvedResidueError(Exception):
    def __init__(self, *args):
        if args:
            self.message = args[0]
        else:
            self.message = None

    def __str__(self):
        if self.message:
            return self.message
        else:
            "NonResolvedResidueError has been raised"
