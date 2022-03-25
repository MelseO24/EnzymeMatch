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
import os.path
import utils.utilities
from utils import utilities
from utils import globalParameters
from utils.utilities import raiseError
from parsers.parse_ResidueGroups import parseResidueGroups
from parsers.parse_inputFile import MatchResidueParameter, MatchTriangleParameter, MatchGlobalParameter,\
    parseGeneral, parseResidue, parseTriangle, parseAminoAcids

fout = sys.stdout
ferr = sys.stderr


def parseParameter(filepath):
    """
    Parses input file, mode dependent on mode given in input file (match_residue or match_triangle)
    If running triangle matching, please note that the distance list should be given in a very specific order:
    distances = [pos12, pos 13, ..., pos 1n, pos23,...,pos2n,...,pos n-1 n]
    :rtype (str, MatchGlobalParameter, MatchResidueParameter, MatchTriangleParameter)
    """

    global residue_groups
    residue_groups = {}
    input_general: MatchGlobalParameter = MatchGlobalParameter()
    input_residue: MatchResidueParameter = MatchResidueParameter()
    input_triangle: MatchTriangleParameter = MatchTriangleParameter()
    all_input_parameters = {}

    fpar = utilities.openFile(filepath, 'r')

    # Locally parse all options in list all_input_parameters
    for line in fpar.readlines():
        if not (line.strip().startswith(";")) and (not len(line.strip()) == 0):
            # allow comment lines in input file, and ignore empty lines

            input_option, setting = utils.utilities.parseInputLine(line)

            if input_option == ">ResidueGroup_input":
                residue_groups = parseResidueGroups(setting)
            if input_option not in all_input_parameters.keys():  # check if input is given more than once
                if input_option == ">Aminoacids":
                    all_input_parameters["Aminoacids_OrigInput"] = [setting]
                    if setting in getResidueGroups():
                        all_input_parameters[input_option] = [getResidueGroupResidues(setting)]
                    else:
                        all_input_parameters[input_option] = [setting]
                elif input_option == ">MT_distances":
                    all_input_parameters[input_option] = [setting]
                else:
                    all_input_parameters[input_option] = setting
            elif input_option == ">Aminoacids":
                all_input_parameters["Aminoacids_OrigInput"].append(setting)
                if setting in getResidueGroups():
                    all_input_parameters[input_option].append(getResidueGroupResidues(setting))
                else:
                    all_input_parameters[input_option].append(setting)
            elif input_option == ">MT_distances":
                all_input_parameters[input_option].append(setting)
            else:
                ferr.write(f"FATAL: {input_option} found twice in input file.\n")
                sys.exit(1)
    #print(all_input_parameters) # can't use if globalParameters.debug, since debugging mode is set later

    if ">setCWD" not in all_input_parameters.keys():
        raiseError("WARNING: Working directory not provided in input file, location of input file taken as CWD.\n", False)
        globalParameters.setWorkdir(os.path.dirname(os.path.realpath(filepath)))
    else:
        globalParameters.workdir = all_input_parameters[">setCWD"]

    try:
        if "complete" in all_input_parameters["mode"]:
            # Has to be set to complete mode, because the BiS.out file does not exist yet for the triangleMatch to work
            globalParameters.setCompleteMode()
    except:
        raiseError("FATAL: No specified mode in parameter file.\n", True)

    # Set options in corresponding parameter object
    if ">MT_pickle" in all_input_parameters.keys():
        try:
            import pickle
            if not os.path.isfile(all_input_parameters[">MT_pickle"]):
                sys.stderr.write(f"FATAL: File {all_input_parameters['>MT_pickle']} does not exist.\n")
            else:
                pickleFile = utils.utilities.openFile(all_input_parameters[">MT_pickle"], "rb")
                input_general, input_triangle = pickle.load(pickleFile)

                if ">Aminoacids" in all_input_parameters.keys():
                    # If aminoacids are provided in inputfile, reset the settings from pickle to be able to overwrite
                    input_general.allowed_residues = []
                    input_general.aminoacid_inputs = []
        except:
            raiseError("FATAL: Failed to import 'pickle' module. Either do not provide '>MT_pickle, or install pickle python module.\n", True)

    for key in all_input_parameters.keys():
        # General options
        if key in parseGeneral(input_general, None, only_return_options=True):
            if key == ">Aminoacids":
                parseAminoAcids(input_general, all_input_parameters[key], option=">Aminoacids")
            else:
                parseGeneral(input_general, key)(all_input_parameters[key])

        # ResidueMatch options
        if key in parseResidue(input_residue, None, only_return_options=True):
            parseResidue(input_residue, key)(all_input_parameters[key])
        # TriangleMatch options
        if key in parseTriangle(input_triangle, None, only_return_options=True):
            if key == ">MT_distances":
                for distance in all_input_parameters[key]:
                    parseTriangle(input_triangle, key)(distance)
            else:
                parseTriangle(input_triangle, key)(all_input_parameters[key])
        if key not in (list(parseGeneral(input_general, None, only_return_options=True))
                       + list(parseResidue(input_residue, None, only_return_options=True))
                       + list(parseTriangle(input_triangle, None, only_return_options=True)))\
                and key not in ["mode", ">setCWD", ">MT_pickle"]:
            raiseError(f"FATAL: Option '{key}' unknown.\n", True)

    ferr.flush()
    # find out which mode to run in
    if "mode" not in all_input_parameters.keys():
        raiseError("FATAL: No specified mode in parameter file.\n", True)
    elif "residue" in all_input_parameters["mode"]:
        input_residue.isValid(input_general)
        input_general.isValid(input_triangle)
        # Only return general parameters and residue parameters
        return "MatchResidue", input_general, input_residue, None
    elif "triangle" in all_input_parameters["mode"]:
        input_triangle.isValid(input_general)
        input_general.isValid(input_triangle)
        # Only return general parameters and triangle parameters
        return "MatchTriangle", input_general, None, input_triangle
    elif "complete" in all_input_parameters["mode"]:
         input_residue.isValid(input_general, parm_triangle=input_triangle, matchcomplete=True)
         # Check if the MatchTriangle outputs have to be saved into the same folder as the MatchResidue
         if input_triangle.MT_directory_out is None:
             raiseError("WARNING: No output directory for MatchTriangle given. Saving outputs in the residue directory\n", False)
             input_triangle.setDirectory(input_residue.MR_directory_out)
         # the output of the MatchResidue consists of a Bis file and a Cat file.
         # The file used for MatchTriangle is determined by the parameter MatchFile in general param
         if input_triangle.input is None:
             if input_general.MatchSide == "bis":
                input_triangle.setInputPath(os.path.join(input_residue.MR_directory_out, input_residue.BiSOut))
             elif input_general.MatchSide == "cat":
                input_triangle.setInputPath(os.path.join(input_residue.MR_directory_out, input_residue.CatOut))
                input_triangle.setMatchSide("cat")
             else:
                 raiseError(f"FATAL: Unknown error on choosing MatchSide\n", True)
         else:
             raiseError("FATAL: Illegal input: Parameter >MT_input not needed in complete mode.\n", True)

         input_triangle.isValid(input_general)
         input_general.isValid(input_triangle)
         # Return both residue and triangle parameters
         return "complete", input_general, input_residue, input_triangle
    elif "aqd" in all_input_parameters["mode"].lower():
        input_triangle.isValid(input_general, aqd_standalone=True)
        input_general.isValid(input_triangle)
        # Only return general parameters and triangle parameters
        return "aqd", input_general, None, input_triangle
    else:
        raiseError("FATAL: No valid mode specified in parameter file. Choose from 'residue', 'triangle', 'complete' or 'AQD'.\n", True)


def getResidueGroups():
    if residue_groups:
        return residue_groups.keys()
    else:
        raiseError("FATAL: Illegal input: '>ResidueGroup_input' has to be provided at the beginning of the input file, "
                   "directly after the mode specification.\n", True)


def getResidueGroupResidues(res_group):
    if residue_groups:
        if res_group in residue_groups.keys():
            return residue_groups[res_group]
        else:
            raiseError(f"FATAL: Illegal input: requested residue group {res_group} is not known. "
                       f"Add this residue group to the file given in '>ResidueGroup input'.\n", True)
    else:
        raiseError("FATAL: Illegal input: '>ResidueGroup_input' has to be provided at the beginning of the input file, "
                   "directly after the mode specification.\n", True)