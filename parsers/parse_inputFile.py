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

import sys
import os
import itertools
import utils.globalParameters
from parsers.parse_ResidueGroups import parseResidueGroups
from utils.utilities import raiseError, _doesFileExist

from utils.utilities import ferr

fout = sys.stdout

"""
HowTo add a new input keyword:
In order to add a new keyword, follow these 4 steps:
1) add keyword and set_function in the respective class (MatchGlobalParameter, MatchResidueParameter, MatchTriangleParameter)
2) add keyword in switcher_errormsg, with respecting object key and input keyword
3) add keyword in corresponding helper_function (parseGeneral/parseResidue/parseTriangle)
4) add checks in <class>_isValid function
"""

switcher_errormsg = {
    # general parameters
    "allowed_residues": ">Aminoacids",
    "aminoacid_inputs": ">Aminoacids",
    "verbosity": ">verbosity",
    "debug": ">debug",
    "strict_mode": ">strict_mode",
    "backbone_interaction": ">backbone_interaction",
    "min_num_of_interaction_points":">min_interactionpoints_num",
    "MatchSide":">MC_MatchSide",

    # residue matching
    "DBfile": ">DBfile",
    "SummaryOut": ">MR_summaryfile",
    "MR_directory_out": ">MR_outdir",
    "number_of_interaction_points":">interactionpoints_num",
    "BiSOut": ">MR_bis_outfile",
    "CatOut": ">MR_cat_outfile",

    "EC": ">MR_EC",
          
    # triangle matching
    "input": ">MT_input",
    "ligand_folder": ">MT_Ligand_folder_path",
    "receptor_folder": ">MT_Receptor_folder_path",
    "summaryFile": ">MT_summaryfile",
    "outFile": ">MT_outfile",
    "MT_directory_out": ">MT_outdir",
    "offset": ">MT_offset_fixed",
    "distances": ">MT_distances",
    "auto_query_input": ">MT_AQD_ligand_input",
    "auto_query_input_format": ">MT_AQD_ligand_inputformat",
    "auto_query_atoms": ">MT_AQD_atomids",
    "auto_query_types": ">MT_AQD_types",
    "auto_query_nr_conformations": ">MT_AQD_conf_num",
    "auto_query_offset_factor": ">MT_AQD_offsetfactor",
    "scoring":">MT_scoring",
    "download_mode":">MT_download_mode",
    "download_dir":">MT_download_dir",
    "download_format":">MT_download_format",
    "backbone_penalty":">MT_backbone_penalty",
    "auto_query_min_offset":">MT_AQD_min_offset",
    "accurate_mode":">MT_accurate_mode",
    "match_side":">MT_MatchSide",
    "numprocs":">MT_numprocs"
}

def isValidAminoAcid(amino_acid):
    valid_amino_acids = ["R", "H", "K", "D", "E", "S", "T", "N", "Q", "C", "U", "G",
                         "P", "A", "V", "I", "L", "M", "F", "Y", "W"]
    # Check if the aminoacid is valid. Whitespaces are ignored.
    if amino_acid.strip() in valid_amino_acids:
        return True
    else:
        raiseError(f"Invalid input: '{amino_acid.strip()}' is not a valid amino acid.\n"
                   f"Choose from {valid_amino_acids}\n", True)


def parseGeneral(parser_object, options, only_return_options=False):
    """
    Helper function to set options in parameter object General
    :param parser_object: object of MatchGlobalParameter class
    :param options: option to set
    :param only_return_options: If True, all possible options will be returned, nothing will be set.
    :return: If only_return_switcher=False, the function to set the requested option is returned and
        the arguments need to be provided manually
    """
    switcher_general = {
        ">Aminoacids": parser_object.setAllowedResidues,
        ">debug": parser_object.setDebug,
        ">ResidueGroup_input": parser_object.setResidueGroups,
        ">verbosity": parser_object.setVerbosity,
        ">strict_mode": parser_object.setStrictMode,
        ">backbone_interaction":parser_object.setBackboneInteraction,
        "Aminoacids_OrigInput": parser_object.setAminoacidInputs,
        ">MC_MatchSide":parser_object.setMatchSide
    }
    if only_return_options:
        return switcher_general.keys()
    else:
        func = switcher_general.get(options)
        return func


def parseResidue(parser_object, options, only_return_options=False):
    """
    Helper function to set options in parameter object Residue
    :type parser_object: MatchResidueParameter
    :param parser_object: object of MatchResidueParameter class
    :param options: option to set
    :param only_return_options: If True, all possible options will be returned, nothing will be set.
    :return: If only_return_switcher=False, the function to set the requested option is returned and
     the arguments need to be provided manually
    """
    switcher_residue = {
        ">DBfile": parser_object.setDBfile,
        ">MR_summaryfile": parser_object.setSummaryOut,
        ">MR_bis_outfile": parser_object.setBiSOut,
        ">MR_cat_outfile": parser_object.setCatOut,
        ">interactionpoints_num": parser_object.setInteractionPoints,
        ">MR_EC": parser_object.setEC,
        ">MR_outdir": parser_object.setDirectory,
        ">min_interactionpoints_num":parser_object.minNumInteractionPoints
        }
    if only_return_options:
        return switcher_residue.keys()
    else:
        func = switcher_residue.get(options)
        return func


def parseTriangle(parser_object, options, only_return_options=False):
    """
    Helper function to set options in parameter object Triangle
                 # By setting :type to the object, PyCharm knows which class it is
    :type parser_object: MatchTriangleParameter
    :param parser_object:  input_residue: object of MatchPolygonParameter class
    :param options: option to set
    :param only_return_options: If True, all possible options will be returned, nothing will be set.
    :return: If only_return_switcher=False, function to set requested option is returned, and
    arguments need to be provided manually
    """
    switcher_triangle = {
        ">MT_input": parser_object.setInputPath,
        ">MT_Ligand_folder_path": parser_object.setLigandPath,
        ">MT_Receptor_folder_path": parser_object.setReceptorPath,
        ">MT_summaryfile": parser_object.setSummaryPath,
        ">MT_offset_fixed": parser_object.setDistanceOffset,
        ">interactionpoints_num": parser_object.setNumberOfInteractionPoints,
        ">MT_outfile": parser_object.setOutPath,
        ">MT_distances": parser_object.setDistancesViaInputFile,
        ">MT_outdir": parser_object.setDirectory,
        ">MT_AQD_ligand_input": parser_object.setAutoQueryInput,
        ">MT_AQD_ligand_inputformat": parser_object.setAutoQueryInputFormat,
        ">MT_AQD_atomids": parser_object.setAutoQueryAtoms,
        ">MT_AQD_types": parser_object.setAutoQueryTypes,
        ">MT_AQD_conf_num": parser_object.setAutoQueryNrConformations,
        ">MT_AQD_offsetfactor": parser_object.setAutoQueryOffsetFactor,
        ">MT_scoring":parser_object.setScoreFunction,
        ">MT_download_mode":parser_object.setDownloadParam,
        ">MT_download_dir":parser_object.setDownloadDir,
        ">MT_download_format":parser_object.setDownloadFormat,
        ">MT_backbone_penalty":parser_object.setBackbonePenalty,
        ">min_interactionpoints_num":parser_object.minNumInteractionPoints,
        ">MT_AQD_min_offset":parser_object.setAQDMinOffset,
        ">MT_accurate_mode":parser_object.setAccPosMode,
        ">MT_MatchSide":parser_object.setMatchSide,
        ">MT_numprocs":parser_object.setNrProcs
    }
    if only_return_options:
        return switcher_triangle.keys()
    else:
        func = switcher_triangle.get(options)
        return func


def parseAminoAcids(parm_general, settings, option=">Aminoacids"):
    """
    sets aminoacids in the residueparse parser
    :type parm_general: MatchGlobalParameter
    :type settings" list[str]
    :param parm_general: MatchGlobalParameter object
    :param settings: list of all aminoacids/settings to set
    :param option: setting to set, default=">AminoAcids"
    """

    for AA in settings:
        for amino_acid in AA.split(','):
            isValidAminoAcid(amino_acid)
        setting = ("*,".join(AA.replace(" ","").split(',')) + '*').split(',')
        parseGeneral(parm_general, option)(setting)


def printMissingVars(errors):
    """
    List errors and stops execution of the script.
    :param errors:
    """
    for missing_var in errors:
        missing_var = switcher_errormsg.get(missing_var)
        ferr.write(missing_var)
        ferr.write("\n")
    exit(1)


class MatchGlobalParameter:
    def __init__(self):
        self.allowed_residues = []
        self._residue_groups = parseResidueGroups(os.path.abspath(utils.globalParameters.path_to_main +
                                                                   "/example/input_resgroups.in"))
        self.debug = None
        self.strict_mode = False
        self.verbosity = 1
        self.backbone_interaction = None
        self.aminoacid_inputs= []
        self.MatchSide = None

    def setMatchSide(self, param):
        self.MatchSide = param.lower()

    def setAminoacidInputs(self, aminoAcid):
        if isinstance(aminoAcid, list):
            self.aminoacid_inputs = aminoAcid
        else:
            self.aminoacid_inputs.append(aminoAcid)

    def setAllowedResidues(self, allowed_residues):
        self.allowed_residues.append(allowed_residues)

    def setResidueGroups(self, path_residue_groups):
        parsedGroups = parseResidueGroups(path_residue_groups)
        for group in parsedGroups.keys():
            if group in self._residue_groups.keys():
                ferr.write(f"WARNING: default residuegroup '{group}' overwritten.\n")
            self._residue_groups[group] = parsedGroups[group]

    def setBackboneInteraction(self, setting):
        try:
            setting = int(setting)
        except:
            raiseError(f"FATAL: {switcher_errormsg['backbone_interaction']} only takes 0 (off) or 1 (on) as argument. \n", True)
        if setting == 0:
            self.backbone_interaction = 0
        elif setting == 1:
            self.backbone_interaction = 1
        else:
            raiseError(
                f"FATAL: {switcher_errormsg['backbone_interaction']} only takes 0 (off) or 1 (on) as argument. \n", True)

    def setDebug(self, setting):
        try:
            setting = int(setting)
        except:
            raiseError(f"FATAL: {switcher_errormsg['debug']} only takes 0 (off), 1 (low) or 2 (high) as argument.\n", True)

        if setting in [0,1,2]:
            self.debug = setting  # just for remembering, such that it is printed with allAttributes()
            utils.globalParameters.setDebug(setting)
            if setting > 0:
                fout.write(f"Debugging mode switched onto {setting}...\n")
        else:
            raiseError(f"FATAL: {switcher_errormsg['debug']} only takes 0 (off), 1 (low) or 2 (high) as argument.\n", True)

    def setStrictMode(self, setting):
        try:
            setting = int(setting)
        except:
            raiseError(f"FATAL: {switcher_errormsg['strict_mode']} only takes 0 (off) or 1 (on) as argument.\n", True)

        if setting == 0:
            self.strict_mode = False  # just for remembering, such that it is printed with allAttributes()
            utils.globalParameters.setStrictMode(False)
        elif setting == 1:
            self.strict_mode = True  # just for remembering, such that it is printed with allAttributes()
            utils.globalParameters.setStrictMode(True)
            fout.write("Strict mode switched on...\n")
        else:
            raiseError(f"FATAL: {switcher_errormsg['strict_mode']} only takes 0 (off) or 1 (on) as argument.\n", True)

    def setVerbosity(self, verbosity_value):
        try:
            verbosity_value = int(verbosity_value)
        except:
            raiseError(f"Invalid input: {verbosity_value} is not a valid input for {switcher_errormsg['verbosity']}. \n"
                       f"Allowed values: 1 (low verbosity) or 2 (high verbosity).\n", True)

        if not verbosity_value in [1, 2]:
            raiseError(f"Invalid input: {verbosity_value} is not a valid input for {switcher_errormsg['verbosity']}. \n"
                       f"Allowed values: 1 (low verbosity) or 2 (high verbosity).\n", True)
        else:
            self.verbosity = verbosity_value    # just for remembering, such that it is printed with allAttributes()
            utils.globalParameters.setVerbosity(verbosity_value)

    def allAttributes(self):
        """"
        Returns all attributes of the object, including values
        :return dictionary containing attributes and their values
        :rtype: Dict[str: str]
        """
        attributes_and_values = {}
        for attribute in vars(self).keys():
            if not attribute.startswith("_"):  # do not return internal attributes
                attributes_and_values[switcher_errormsg[attribute]] = vars(self)[attribute]
        return attributes_and_values

    def isValid(self, parm_triangle):
        """
        Checks if all required parameters are set to run residue matching
        :param parm_triangle: parm_triangle object, in order to check if automatic query design switched on
        :type parm_triangle: parsers.parse_inputFile.MatchTriangleParameter
        :return: True if valid, or throws exception if not
        """
        # check for unset variables
        errors = []
        for key in vars(self).keys():
            if key == "MatchSide":
                continue
            if (vars(self)[key] is None) or (vars(self)[key] == []):
                errors.append(key)


        # MatchSide type can only be given when the complete mode is turned on,
        # because in MatchResidue, we trivially don't need any input files
        # and in MatchTriangle, the input file has to be specified anyway.
        if utils.globalParameters.complete_mode and self.MatchSide is None:
            raiseError(f"FATAL: {switcher_errormsg['MatchSide']} has to be given in complete Mode\n"
                          f"Available options: bis, cat\n", True)
        elif not utils.globalParameters.complete_mode and self.MatchSide is not None:
            raiseError(f"FATAL: {switcher_errormsg['MatchSide']} cannot be given when not in complete Mode\n", True)

        if utils.globalParameters.complete_mode:
            if self.MatchSide not in ['cat', 'bis']:
                raiseError(f"FATAL: {switcher_errormsg['MatchSide']} is invalid\n"
                          f"Available options: bis, cat\n", True)

        if not errors == []:
            if not (("allowed_residues" in errors) and parm_triangle.auto_query_input):
                raiseError("FATAL: The following variables are essential: \n", False)
                printMissingVars(errors)

        # if backbone interaction is turned on, but any single amino acid query is given, exit
        if self.backbone_interaction == 1 and not all([i.strip() in list(self._residue_groups.keys()) for i in self.aminoacid_inputs]):
            raiseError("FATAL: Only residue groups can be given when backbone interaction is turned on\n", True)

class MatchResidueParameter:
    def __init__(self):
        """
        self.DBfile = file containing BioLip database
        self.SummaryOut = summary output file name
        self.BiSOut = binding site output file name
        self.CatOut = catalytic site output file name
        self.number_of_interactionPoints = number of interaction points
        self.EC = "Filtering for specific Enzyme Classes"
        """
        self._errors = []
        self.DBfile = None
        self.MR_directory_out = None
        self.SummaryOut = "MatchResidue.summary"
        self.BiSOut = "BiS_matches.out"
        self.CatOut = "Cat_matches.out"
        self.number_of_interaction_points = None
        self.EC = []
        self.min_num_of_interaction_points = None

    # setter for parameter attributes
    def setDBfile(self, fname):
        self.DBfile = fname

    def setSummaryOut(self, summary_output):
        self.SummaryOut = summary_output

    def setBiSOut(self, binding_site_output):
        self.BiSOut = binding_site_output

    def setCatOut(self, catalytic_site_output):
        self.CatOut = catalytic_site_output

    def setInteractionPoints(self, interaction_points: int):
        try:
            self.number_of_interaction_points = int(interaction_points)
        except ValueError:
            raiseError(f"Invalid input: {interaction_points} is not a valid input for {switcher_errormsg['number_of_interaction_points']}. \n"
                       f"Allowed values: Integers\n", True)
        except:
            raiseError(f"FATAL: An unknown error has occured during parsing of the input.\n"
                             f"Argument: {switcher_errormsg['number_of_interaction_points']}\n", True)
        if self.min_num_of_interaction_points is None:
            self.min_num_of_interaction_points = int(interaction_points)

    def minNumInteractionPoints(self, min_interaction_points):
        try:
            self.min_num_of_interaction_points = int(min_interaction_points)
        except ValueError:
            raiseError(f"Invalid input: {min_interaction_points} is not a valid input for {switcher_errormsg['min_num_of_interaction_points']}. \n"
                       f"Allowed values: Integers\n", True)
        except:
            raiseError(f"FATAL: An unknown error has occured during parsing of the input.\n"
                             f"Argument: {switcher_errormsg['min_num_of_interaction_points']}\n", True)
        if int(min_interaction_points) < 2:
            raiseError(
                f"Invalid input: {switcher_errormsg['min_num_of_interaction_points']} has to be greater than or equal to 2\n", True)

    def setEC(self, EC):
        EClist = EC.split(".")
        for num in EClist:
            if num == "":
                raiseError(f"FATAL: The EC number has the false format\n"
                           f"Right format: point separated integers without points in the end\n", True)
        self.EC = EC.split(".")

    def setDirectory(self, directory):
        self.MR_directory_out = directory

    def get_DBfile(self):
        if _doesFileExist(self.DBfile):
            return self.DBfile
        elif _doesFileExist(os.path.join(utils.globalParameters.workdir, self.DBfile)):
            return os.path.join(utils.globalParameters.workdir, self.DBfile)
        else:
            raiseError(f"FATAL: Database {self.DBfile} not found.\n", True)

    def allAttributes(self):
        """"
        Returns all attributes of the object, including values
        :return dictionary containing attributes and their values
        :rtype: Dict[str: str]
        """
        attributes_and_values = {}
        for attribute in vars(self).keys():
            if not attribute.startswith("_"):  # do not return internal attributes
                attributes_and_values[switcher_errormsg[attribute]] = vars(self)[attribute]
        return attributes_and_values

    def isValid(self, parm_general, parm_triangle=None, matchcomplete=False):
        """
        Checks if all required parameters are set to run residue matching
        :param matchcomplete: True if matchcomplete is running, because >interactionpoints_num is not required if aqd is used
        :param parm_triangle: Only for matchComplete, to check if AQD is used (needs to be provided when matchcomplete=True!)
        :type parm_general: MatchGlobalParameter
        :type parm_triangle: MatchTriangleParameter
        :type matchcomplete: bool
        :return: True if valid, throws exception if not
        """
        # check for unset variables
        for key in vars(self).keys():
            if vars(self)[key] is None:
                #number_of_interaction_points not requried if matchcomplete runs with AQD
                if (key == "number_of_interaction_points") and (matchcomplete) and (parm_triangle.auto_query_input):
                    continue
                else:
                    self._errors.append(key)

        if self.min_num_of_interaction_points is not None and self.number_of_interaction_points is not None and (self.min_num_of_interaction_points > self.number_of_interaction_points):
            raiseError(f"FATAL: {switcher_errormsg['min_num_of_interaction_points']} is greater than {switcher_errormsg['number_of_interaction_points']}\n", True)

        if self._errors:
            raiseError("FATAL: The following variables are essential: \n", False)
            printMissingVars(self._errors)

        if not self.number_of_interaction_points == len(parm_general.allowed_residues):
            # number_of_interaction_points not requried if matchcomplete runs with AQD
            if not (matchcomplete and parm_triangle.auto_query_input):
                raiseError("FATAL: Illegal input: The amount of '>Aminoacids' entries does not match '>interaction points'.\n", True)

class MatchTriangleParameter:

    def __init__(self):
        """
        self.input = binding site/cat site input file path
        self.ligand_folder = ligand folder file path
        self.receptor_folder = receptor folder file path
        self.number_of_interaction_points = # of binding number_of_interaction_points
        self.offset = distance offset for the distances between interaction_points, {frozenset(pos 1, pos 2):offset}
        self.out = output file name
        self.distances = allowed distance between 2 positions {frozenset(pos 1, pos 2):distance}
        self.score_function = 0:off 1:on(distance based scoring)
        self.download_param = 0:off 1:Downloading without deletion 2:Downloading with deletion
        self.download_dir = pdb download directory
        self.download_format = either "pdb" or "mmtf"
        self.match_side = either "bis" or "cat"
        self.numprocs = number of processors for multiprocessing
        """
        self._errors = []
        self.input = None
        self.ligand_folder = None
        self.receptor_folder = None
        self.summaryFile = "MatchTriangle.summary"
        self.number_of_interaction_points = None
        self.offset = {}
        self.outFile = "MatchTriangle.out"
        self.distances = {}
        self.MT_directory_out = None
        self.auto_query_input = False
        self.auto_query_input_format = False
        self.auto_query_atoms = []
        self.auto_query_types = []
        self.auto_query_nr_conformations = 100
        self.auto_query_offset_factor = 1.0
        self.auto_query_min_offset = 1.0
        self.scoring = 0
        self.download_mode = 0
        self.download_dir = ""
        self.download_format = ""
        self.backbone_penalty = 0.0
        self.min_num_of_interaction_points = None
        self.accurate_mode = 0
        self.match_side = "bis"
        self.numprocs = os.cpu_count()

    def __str__(self):
        return f"distances: {str(self.distances)}"

    def setAQDMinOffset(self, minOffset):
        try:
            self.auto_query_min_offset = float(minOffset)
        except:
            raiseError(
                f"Illegal input: '{switcher_errormsg['auto_query_min_offset']}' needs to be a floating point number.\n", True)
        if self.auto_query_min_offset < 0:
            raiseError(
                f"Illegal input: '{switcher_errormsg['auto_query_min_offset']}' cannot be negative.\n", True)

    def setAccPosMode(self, mode):
        try:
            self.accurate_mode = int(mode)
            if self.accurate_mode not in [0,1]:
                raise Exception
        except:
            raiseError(f"Illegal input: '{switcher_errormsg['accurate_mode']}' only takes 0 (off) and 1 (on) as argument\n", True)

    def setBackbonePenalty(self, penaltyScore):
        try:
            penaltyScore = float(penaltyScore)
        except:
            raiseError(f"Illegal input: '{switcher_errormsg['backbone_penalty']}' needs to be a floating point number.\n", True)
        self.backbone_penalty = penaltyScore

    def setDownloadDir(self, download_path):
        """
        sets the download path for the pdb files
        :param download_path:
        """
        self.download_dir = download_path

    def setDownloadFormat(self, download_format):
        """
        determines the download format
        :param download_format: pdb or mmtf
        """
        if download_format.strip() not in ["mmtf", "pdb"]:
            raiseError(f"Illegal input: file format {download_format} is not supported\n"
                       f"Please enter one of the following: mmtf, pdb\n", True)
        self.download_format = download_format.strip()

    def setInputPath(self, input_path):
        input_file = os.path.join(utils.globalParameters.workdir, input_path)
        if os.path.isfile(input_file):
            self.input = input_file
        # if complete mode is being run, the file does not exist yet and will be created during the matchResidue part.
        elif not utils.globalParameters.complete_mode:
            raiseError(f"Illegal input: file {input_file} does not exist.\n", True)
        else:
            self.input = input_file

    def setLigandPath(self, ligand_folder_path):
        _path_to_ligand = os.path.join(utils.globalParameters.workdir, ligand_folder_path)
        if os.path.isdir(_path_to_ligand):
            self.ligand_folder = _path_to_ligand
        else:
            raiseError(f"Illegal input: folder {_path_to_ligand} does not exist or is unreachable.\n", True)

    def setReceptorPath(self, receptor_folder_path):
        _path_to_receptor = os.path.join(utils.globalParameters.workdir, receptor_folder_path)
        if os.path.isdir(_path_to_receptor):
            self.receptor_folder = _path_to_receptor
        else:
            raiseError(f"Illegal input: folder {_path_to_receptor} does not exist or is unreachable.\n", True)

    def setSummaryPath(self, filename):
        self.summaryFile = filename

    def setOutPath(self, filename):
        self.outFile = filename

    def setNumberOfInteractionPoints(self, interaction_points: int):
        try:
            self.number_of_interaction_points = int(interaction_points)
        except ValueError:
            raiseError(
                f"Invalid input: {interaction_points} is not a valid input for {switcher_errormsg['number_of_interaction_points']}. \n"
                f"Allowed values: Integers\n", True)
        except:
            raiseError(f"FATAL: An unknown error has occured during parsing of the input.\n"
                             f"Argument: {switcher_errormsg['number_of_interaction_points']}", True)
        if self.min_num_of_interaction_points is None:
            self.min_num_of_interaction_points = int(interaction_points)

    def minNumInteractionPoints(self, min_interaction_points):
        try:
            self.min_num_of_interaction_points = int(min_interaction_points)
        except ValueError:
            raiseError(
                f"Invalid input: {min_interaction_points} is not a valid input for {switcher_errormsg['min_num_of_interaction_points']}. \n"
                f"Allowed values: Integers\n", True)
        except:
            raiseError(f"FATAL: An unknown error has occured during parsing of the input.\n"
                             f"Argument: {switcher_errormsg['min_num_of_interaction_points']}", True)
        if int(min_interaction_points) < 2:
            raiseError(
                f"Invalid input: {switcher_errormsg['min_num_of_interaction_points']} has to be greater than or equal to 2\n", True)

    def setDistancesViaInputFile(self, distances):
        """
        Set distances parsed via input file.
        in case the distance between specific distance pairs needs to be set manually,
         use setDistancesManually() instead.
        :param distances: string of input file
        :type distances: string
        """
        try:
            float(distances.split("=")[1].strip())
        except:
            raiseError(f"Illegal input: '{switcher_errormsg['distances']}' needs to be a floating point number.\n", True)
        interaction_points_no = distances.split("=")[0].strip().split("-")
        interaction_points_no = frozenset([int(i)-1 for i in interaction_points_no])
        self.distances[interaction_points_no] = float(distances.split("=")[1].strip())

    def setDistancesManually(self, distance, interaction_point1, interaction_point2):
        """"
        Set distance between 2 specific positions manually.
        If parsing of input file is needed, use setDistancesViaInputFile() instead.
        :param distance: value of optimal distance between 2 positions
        :param interaction_point1: number of first position (counting starts at 0)
        :param interaction_point2: number of second position (counting starts at 0)
        #type offset: float
        :type interaction_point1: int
        :type interaction_point2: int
        """
        if len(self.distances.keys()) == 0:
            """Distances dictionary still has to be initialized, set all values to 0"""
            for vertex_combination in itertools.combinations(list(range(self.number_of_interaction_points)), 2):
                self.distances[frozenset(list(vertex_combination))] = 0

        # Actual parsing to object distances variable
        self.distances[frozenset(list((interaction_point1, interaction_point2)))] = distance

    def setDistanceOffset(self, offset, interaction_point1=0, interaction_point2=0):
        """"
        Set offset values to object variable.
        If positions1 and positions2 == 0, then this offset will be used for all position pairs (which is the case
        in normal input)
        If positions1 and positions2 are provided, this offset will only be set between these 2 positions. (which
        is the case for automatic query design.
        :param offset: value of max allowed offset (on distance) between 2 positions
        :param interaction_point1: number of first position (counting starts at 0)
        :param interaction_point2: number of second position (counting starts at 0)
        :type offset: float
        :type interaction_point1: int
        :type interaction_point2: int
        """
        try:
            offset = float(offset)
        except:
            raiseError(f"Illegal input: '{switcher_errormsg['offset']}' needs to be a floating point number.\n", True)

        if not self.number_of_interaction_points:
            raiseError(f"Illegal input: '{switcher_errormsg['number_of_interaction_points']}' should be provided before "
                       f"'{switcher_errormsg['offset']}' in the input file.\n", True)
        elif len(self.offset.keys()) == 0:
            """Offset dictionary still has to be initialized, set all values to 0"""
            for vertex_combination in itertools.combinations(list(range(self.number_of_interaction_points)), 2):
                self.offset[frozenset(list(vertex_combination))] = offset
        if (interaction_point1 == 0 and interaction_point2 == 0):
            #set all values to 0
            for vertex_combination in itertools.combinations(list(range(self.number_of_interaction_points)), 2):
                self.offset[frozenset(list(vertex_combination))] = offset
        else:
            self.offset[frozenset(list((interaction_point1, interaction_point2)))] = offset

    def setAutoQueryInput(self, ligand_input_path):
        _path_to_input = os.path.join(utils.globalParameters.workdir, ligand_input_path)
        if os.path.isfile(_path_to_input):
            self.auto_query_input = _path_to_input
        else:
            raiseError(f"Illegal input: folder {_path_to_input} does not exist or is unreachable.\n", True)
        self.auto_query_input = _path_to_input

    def setAutoQueryInputFormat(self, ligand_format):
        allowed_formats = ["mol2", "sdf"]
        if ligand_format not in allowed_formats:
            raiseError(f"Illegal input: {switcher_errormsg['auto_query_input_format']} only accepts following input: {allowed_formats}\n", True)
        else:
            self.auto_query_input_format = ligand_format

    def setAutoQueryAtoms(self, atoms):
        """
        :param atoms: string of atomIDs which should interact with the binding site
        :type atoms: string
        """
        try:
            self.auto_query_atoms = [int(atomID) for atomID in atoms.split(",")]
        except:
            raiseError(f"Illegal input: '{switcher_errormsg['auto_query_atoms']}''"
                       f" only takes list of integers as input.\n", True)

    def setAutoQueryTypes(self, types):
        """
        :param types: list of interaction types of atoms defined in self.auto_query_atoms
        :type types: string
        """
        self.auto_query_types = [interaction_type.strip() for interaction_type in types.split(",")]

    def setAutoQueryNrConformations(self, nr_conformations):
        """
        :param nr_conformations: number of conformations
        :type nr_conformations: int
        """
        try:
            int(nr_conformations)
            if int(nr_conformations) < 50:
                raiseError(f"WARNING: {switcher_errormsg['auto_query_nr_conformations']} is lower than 50, be aware "
                           f"that low values may result in inaccurate queries, and thus results. "
                           f"Consider increasing this value.\n", False)
            self.auto_query_nr_conformations = int(nr_conformations)
        except:
            raiseError(f"Illegal input: '{switcher_errormsg['auto_query_nr_conformations']}''"
                       f" only takes integer as input.\n", True)

    def setAutoQueryOffsetFactor(self, factor):
        """
        :param factor: factor to multiply distance variance with (from ligand conf. search), to retrieve allowed offset
        :type factor: float
        """
        try:
            factor = float(factor)
        except:
            raiseError(f"Illegal input: '{switcher_errormsg['auto_query_offset_factor']}'"
                       f" needs to be a floating point number.\n", True)

        self.auto_query_offset_factor = factor

    def setDirectory(self, directory):
        self.MT_directory_out = directory

    def allAttributes(self):
        """"
        Returns all attributes of the object, including values
        :return dictionary containing attributes and their values
        :rtype: Dict[str: str]
        """
        attributes_and_values = {}
        for attribute in vars(self).keys():
            if not attribute.startswith("_"):  # do not return internal attributes
                # Rewrite the positions for the summary output
                if attribute in ["distances", "offset"]:
                    outputDict = {}
                    for positionsets in vars(self)[attribute]:
                        rightPos = [str(position + 1) for position in list(positionsets)]
                        outputDict["position " + "-".join(rightPos)] = vars(self)[attribute][positionsets]
                    attributes_and_values[switcher_errormsg[attribute]] = outputDict
                else:
                    attributes_and_values[switcher_errormsg[attribute]] = vars(self)[attribute]
        return attributes_and_values

    def setScoreFunction(self, ScoreMode):
        """
        :type ScoreMode: int
        :param ScoreMode: chosen mode of the scoring function, if 0, turned off
        """
        try:
            ScoreMode = int(ScoreMode)
        except:
            raiseError(f"FATAL: {switcher_errormsg['scoring']} only takes 0 (off) or 1 (Distance-based) as argument.\n", True)

        if ScoreMode == 0:
            self.scoring = 0
        elif ScoreMode == 1:
            self.scoring = 1
            fout.write("Scoring mode: Distance based...\n")
        else:
            raiseError(f"FATAL: {switcher_errormsg['scoring']} only takes 0 (off) or 1 (Distance-based) as argument.\n", True)

    def setDownloadParam(self, downloadMode):
        """
        :type downloadMode: int
        :param downloadMode: chosen mode of the download mode
        """
        try:
            downloadMode = int(downloadMode)
        except:
            raiseError(
                f"FATAL: {switcher_errormsg['download_mode']} only takes 0 (off), 1 (on) as argument\n", True)
        if downloadMode == 0:
            self.download_mode = 0
        elif downloadMode == 1:
            self.download_mode = 1
        else:
            raiseError(
                f"FATAL: {switcher_errormsg['download_mode']} only takes 0 (off), 1 (on) as argument\n", True)

    def setMatchSide(self, side):
        if side in ["bis", "cat"]:
            self.match_side = side
        else:
            raiseError("Illegal input: >MatchSide only accepts the following keywords: 'bis', 'cat'\n", True)

    def setNrProcs(self, procs):
        try:
            self.numprocs = int(procs)
        except:
            raiseError(f"Illegal input: '{switcher_errormsg['numprocs']}' needs to be an integer.\n", True)


    def isValid(self, parm_general, aqd_standalone=False):
        """
        Checks if all required parameters are set to run triangle matching
        :param parm_general: parm_general object to retrieve allowed residues (>Aminoacids input)
        :param aqd_standalone: True if mode=AQD, since less needs to be set
        :type parm_general: parsers.parse_inputFile.MatchGlobalParameter
        :type aqd_standalone, bool
        :return: True if valid, throws exception if not
        """
        # check for unset variables
        # don't check these, since they will be checked later in this isValid() function
        variables_no_check = ['auto_query_offset_factor', 'auto_query_min_offset']
        if aqd_standalone:
            aqd_no_check = ['input', 'ligand_folder', 'receptor_folder', 'scoring',
                            'download_mode', 'min_num_of_interaction_points', 'accurate_mode']
            variables_no_check = variables_no_check + aqd_no_check

        for key in vars(self).keys():
            if vars(self)[key] is None:
                if key not in variables_no_check:
                    self._errors.append(key)

        if self._errors:
            raiseError("Illegal input: The following variables are essential: \n", False)
            printMissingVars(self._errors)

        # If download parameter = 0, no download dir can be given. If download parameter != 0, download dir has to be given
        if not aqd_standalone:
            if self.download_mode == 0 and (self.download_dir is not "" or self.download_format is not ""):
                raiseError(f"Illegal input: '{switcher_errormsg['download_dir']}' and '{switcher_errormsg['download_format']}' "
                           f"shouldn't be given when download function is turned off\n", True)
            elif self.download_mode != 0 and (self.download_dir is "" or self.download_format is ""):
                raiseError(f"Illegal input: '{switcher_errormsg['download_dir']}' and '{switcher_errormsg['download_format']}' "
                           f"should be given when download function is turned on\n", True)

        # penalty can only be given when the backbone interaction is turned on
        if self.backbone_penalty != 0.0 and (parm_general.backbone_interaction == 0 or self.scoring == 0):
            raiseError(f"Illegal input: '{switcher_errormsg['backbone_penalty']}' shouldn't be given when "
                       f"'{switcher_errormsg['backbone_interaction']}' is turned off"
                       f" or '{switcher_errormsg['scoring']}' is turned off\n", True)

        # when any automatic query design input given, all of them need to be provided
        if not self.auto_query_input:
            if self.auto_query_atoms:
                raiseError(f"Illegal input: '{switcher_errormsg['auto_query_atoms']}' provided, which requires "
                           f"'{switcher_errormsg['auto_query_input']}' input as well.\n", True)
            if self.auto_query_types:
                raiseError(f"Illegal input: '{switcher_errormsg['auto_query_types']}' provided, which requires "
                           f"'{switcher_errormsg['auto_query_input']}' input as well.\n", True)
            if self.auto_query_input_format:
                raiseError(f"Illegal input: '{switcher_errormsg['auto_query_input_format']}'"
                           f" is provided, which requires "
                           f"'{switcher_errormsg['auto_query_input']}' as well.\n", True)
            if aqd_standalone:
                raiseError(f"Illegal input: '{switcher_errormsg['auto_query_input']}' needs to be provided when run in 'AQD' mode.\n", True)
        else:
            if not self.auto_query_atoms:
                raiseError(f"Illegal input: '{switcher_errormsg['auto_query_input']}' provided, which requires "
                           f"'{switcher_errormsg['auto_query_atoms']}' input as well.\n", True)
            if not self.auto_query_types:
                raiseError(f"Illegal input: '{switcher_errormsg['auto_query_input']}' provided, which requires "
                           f"'{switcher_errormsg['auto_query_types']}' input as well.\n", True)
            if not self.auto_query_offset_factor:
                raiseError(f"Illegal input: '{switcher_errormsg['auto_query_input']}' is provided, which requires "
                           f"'{switcher_errormsg['auto_query_offset_factor']}' as well.\n", True)
            if not self.auto_query_input_format:
                raiseError(f"Illegal input: '{switcher_errormsg['auto_query_input']}' is provided, which requires "
                           f"'{switcher_errormsg['auto_query_input_format']}' as well.\n", True)
            if parm_general.allowed_residues:
                raiseError(f"Illegal input: if '{switcher_errormsg['auto_query_input']}'"
                           f" provided, no '>Aminoacids' are "
                           "allowed, since this is created internally by the automatic query design.\n", True)
            if not self.auto_query_min_offset:
                if self.auto_query_min_offset == 0:
                    raiseError(f"WARNING: It is not advised to used the AQD minimum offset of 0. "
                               f"Consider using a higher value.\n", False)
                else:
                    raiseError(f"Illegal input: '{switcher_errormsg['auto_query_input']}' is provided, which requires "
                               f"'{switcher_errormsg['auto_query_min_offset']}' as well.\n", True)
            if self.distances:
                raiseError("Illegal input: both '>MT_distances' and '>automatic query design ligand input' provided.\n", True)


        # check if query types are valid
        if self.auto_query_types:
            valid_query_types = ["H_acc", "H_don", "pos_charged", "neg_charged", "hydrophobic",
                                 "aromatic_ring", "hydrophobic_ring"]
            for query_type in self.auto_query_types:
                if query_type not in valid_query_types:
                    raiseError(f"{query_type} is not a valid input for '{switcher_errormsg['auto_query_types']}'. "
                               f"Choose from {valid_query_types}.\n", True)

        if self.auto_query_atoms:
            # check if number of provided ligand atoms equals number of provided types
            if not len(self.auto_query_atoms) == len(self.auto_query_types):
                raiseError(f"Illegal input: number of entries in '{switcher_errormsg['auto_query_atoms']}' should be "
                           f"equal to '{switcher_errormsg['auto_query_types']}'.", True)
                sys.exit(1)
            # check if given number_of_interaction_points match number of given types
            if not self.number_of_interaction_points == len(self.auto_query_atoms) == len(self.auto_query_types):
                if not aqd_standalone:
                    raiseError(f"Illegal input: The following input arguments have to contain an equal amount of entries:"
                               f"'{switcher_errormsg['number_of_interaction_points']}', "
                               f"''{switcher_errormsg['auto_query_atoms']}, "
                               f"'{switcher_errormsg['auto_query_types']}'.\n", True)
            if not len(self.offset.keys()) == 0:
                raiseError(f"Illegal input: '{switcher_errormsg['offset']}' "
                           f"cannot be combined with automatic query design.\n", True)

        if not self.auto_query_input:
            if not self.number_of_interaction_points == len(parm_general.allowed_residues):
                raiseError(f"Illegal input: The amount of '{switcher_errormsg['allowed_residues']}' entries ({len(parm_general.allowed_residues)}) does not "
                           f"match '{switcher_errormsg['number_of_interaction_points']} ({self.number_of_interaction_points})'.\n", True)
            if not len(self.distances.keys()) == ((self.number_of_interaction_points - 1) * (self.number_of_interaction_points) / 2):
                if len(self.distances.keys()) > ((self.number_of_interaction_points - 1) * (self.number_of_interaction_points) / 2):
                    raiseError(f"Illegal input: There are more than {self.number_of_interaction_points}"
                               f" '>MT_distances' provided.\n", True)
                else:
                    raiseError(f"Illegal input: There are less than {self.number_of_interaction_points}"
                               f" '>MT_distances' provided.\n", True)
            for vertex_combination in itertools.combinations(list(range(self.number_of_interaction_points)), 2):
                if not frozenset(list(vertex_combination)) in self.distances.keys():
                    raiseError(f"Illegal input: '>MT_distances: {vertex_combination[0] + 1}-{vertex_combination[1] + 1} is missing.\n", True)
            if len(self.offset.keys()) == 0:
                raiseError(f"Illegal input: '{switcher_errormsg['offset']}' is missing.\n", True)

        if self.min_num_of_interaction_points > self.number_of_interaction_points:
            raiseError(f"FATAL: {switcher_errormsg['min_num_of_interaction_points']} is greater than {switcher_errormsg['number_of_interaction_points']}\n", True)

