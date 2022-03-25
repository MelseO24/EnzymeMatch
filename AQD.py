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
import match_triangle
import utils.utilities
from utils import globalParameters
from parsers.parse_inputFile import parseAminoAcids
from utils.utilities import raiseError
from utils.AQDUtils import nrOfBonds,isDoubleBoundHeavyAtom,isHydrogen,isHacceptor,isHdonor,containsHydrogens,calcLinearOffset,calcLonepairOffset,calcSphericOffset, calcRingCenterOffset
from rdkit import Chem
from rdkit.Chem import AllChem
import numpy as np
import itertools
fout = sys.stdout
global aqdStandAlone
aqdStandAlone = False

def runAQDstandalone(parm_triangle, parm_general):
    """
    Runs AQD stand-alone, writes summary file.
    :type parm_triangle: parsers.parse_inputFile.MatchTriangleParameter
    :type parm_general: parsers.parse_inputFile.MatchGlobalParameter
    :return: parm_triangle with calculated distances between interaction points
    """
    global aqdStandAlone
    aqdStandAlone = True
    out_main_path = os.path.join(utils.utilities.globalParameters.workdir, parm_triangle.MT_directory_out)
    if not os.path.exists(out_main_path):
        utils.utilities.makeDirectory(out_main_path)
    summary_output = utils.utilities.openFile(os.path.join(out_main_path, parm_triangle.summaryFile), 'w+')
    utils.utilities.printInputInFile(summary_output, parm_general, parm_triangle)\

    parm_triangle = calcBindingsiteDistances(parm_triangle, parm_general)

    summary_output.write("AQD standalone:\nNew input after AQD:\n")
    utils.utilities.printInputInFile(summary_output, parm_general, parm_triangle)
    for position, allowed_residues in enumerate(parm_general.allowed_residues, 1):
        summary_output.write(f"Allowed residues at position {position}: ")
        for allowed_residue in allowed_residues:
            summary_output.write(allowed_residue + " ")
        summary_output.write("\n")

    #Save parm_triangle as pickle. to allow easy follow-up MatchTriangle; note: all aqd settings are set to default, otherwise matchtriangle mode will explain
    parm_triangle.auto_query_input = False
    parm_triangle.auto_query_types = []
    parm_triangle.auto_query_atoms = []
    parm_triangle.auto_query_input_format = False
    parm_triangle.auto_query_offset_factor = None
    parm_triangle.auto_query_min_offset = None
    try:
        import pickle
        with open("AQDoutput.pickle", "wb") as pickleFile:
            pickle.dump([parm_general, parm_triangle], pickleFile)
        fout.write("TriangleMatching parameters as well as general parameters saved as pickle to: AQDoutput.pickle\n")
    except:
        sys.stderr.write("WARNING: Failed to import 'pickle' module. Therefore no pickle saved.\n")


def calcBindingsiteDistances(parm_triangle, parm_general):
    """
    Main function of conf_search, performing conformational search and returning distances and varience between binding
    site interaction points
    More detailed:
    If no distances provided by user (and AQD activated), then conformational search of ligand is performed
    and distances between interaction groups are predicted. If not, input parm_triangle is returned
    :type parm_triangle: parsers.parse_inputFile.MatchTriangleParameter
    :return: parm_triangle with calculated distances between interaction points
    """
    global StrictMode
    StrictMode = parm_general.strict_mode
    if parm_triangle.auto_query_input:
        fout.write("Automatic query design activated...\n")
        # calculates ligand conformations
        ligand_conformations = ligandConfSearch(parm_triangle)
        # calculates the distances, return dictionary (same as orig input)
        mean_distance_var, label_mean_var = calcBisInteractionPointsDistances(parm_triangle, ligand_conformations)
        _setDistanceAndOffset(parm_triangle, mean_distance_var, label_mean_var) # set distances in parm_triangle to what has been calculated
        setProteinAllowedResidues(parm_triangle, parm_general)
        fout.write("Automatic query design finished...\n")
    return parm_triangle

def ligandConfSearch(parm_triangle):
    """
    given ligand coordinates (mol file),
    write "pdb" file containing coordinates of different conformations using rdkit
    :param parm_triangle: triangle parameters to retrieve path to ligand .mol2 file and number of ligand conformations to generate
    :type parm_triangle: parsers.parse_inputFile.MatchTriangleParameter
    :return: conformers of input molecule
    """
    if globalParameters.debug > 0:
        conformations_out = utils.utilities.openFile(os.path.join(globalParameters.workdir, "confSearch_output.pdb"), "w+")

    if not (parm_triangle.auto_query_input.split(".")[-1]) == parm_triangle.auto_query_input_format:
        raiseError(f"FATAL: The ligand input file {parm_triangle.auto_query_input} is not of"
                         f" {parm_triangle.auto_query_input_format} format.\n", True)

    try:
        if parm_triangle.auto_query_input_format == "mol2":
            ligand_mol = Chem.MolFromMol2File(parm_triangle.auto_query_input, removeHs=False)  # type: Chem.rdchem.Mol
        elif parm_triangle.auto_query_input_format == "sdf":
            ligand_input_supplier = Chem.SDMolSupplier(parm_triangle.auto_query_input, removeHs=False)
            ligand_mol = next(ligand_input_supplier) # type: Chem.rdchem.Mol
    except OSError as error:
        if str(error).startswith("Bad input file"):
            raiseError(f"File {parm_triangle.auto_query_input} not found.\n", True)
        else:
            raise
    if not ligand_mol:
        if parm_triangle.auto_query_input_format == "mol2":
            raiseError(f"\nCould not parse {parm_triangle.auto_query_input}, check your file, "
                             f"did you use Corina atom types? Other atom types are not (fully) supported. Perhaps try to use a '.sdf' file.\n", True)
        else:
            raiseError(f"\nCould not parse {parm_triangle.auto_query_input}, check your file. Perhaps try to use a '.mol2' file.\n", True)

    containsHydrogens(ligand_mol)  # check if user provided protonated structure

    etkdg_ps = AllChem.ETKDG()
    etkdg_ps.randomSeed = 1234
    etkdg_ps.clearConfs = True  # overwrites ligand_mol initial conformation
    fout.write(f"{parm_triangle.auto_query_nr_conformations} ligand conformations are generated...\n")
    AllChem.EmbedMultipleConfs(ligand_mol, parm_triangle.auto_query_nr_conformations, etkdg_ps)  # actual conformational search, creating ligand conformations


    if globalParameters.debug > 0:
        fout.write("Generated ligand conformations written to 'confSearch_output.pdb'.\n")
        conformations_out.write(Chem.MolToPDBBlock(ligand_mol))
    return ligand_mol


def calcBisInteractionPointsDistances(parm_triangle, mol_conformations):
    """
    (previously calcMeanAndVar())
    given ligand conformations, calculate position of protein (bis) interaction points and
    calculate mean of distances and variances between bis interaction points
    :type mol_conformations: Chem.rdchem.Mol
    :type parm_triangle: parsers.parse_inputFile.MatchTriangleParameter
    :return: list([avgDist01, varience01], [avgDist02, varience02],..) and
             a list listing between which atoms the distances are measured: list((atomID1,atomID2), (atomID1,atomID2))
    """
    query_atoms = dict(zip(                                  # dict(atomID1: type, atomID2: type)
        [atom-1 for atom in parm_triangle.auto_query_atoms], # rdkit starts counting at 0
        parm_triangle.auto_query_types))
    if globalParameters.debug > 0:
        print(f"query atoms: {query_atoms}")
    # allconformers_bis_int_points contains for each conformation the coordinates of the binding site interaction points [[conf1(X,Y,Z)], [conf2(X,Y,Z]..]
    allconformers_bis_int_points = [[] for _ in range(len(mol_conformations.GetConformers()))]

    if globalParameters.debug > 0:
        for atom in mol_conformations.GetAtoms():
            print(f"Atom: {atom.GetIdx()+1},"
                  f" atomicNr: {atom.GetAtomicNum()},"
                  f" NrBonds: {nrOfBonds(atom)},"
                  f" Formal charge: {atom.GetFormalCharge()},"
                  f" Aromatic: {atom.GetIsAromatic()},"
                  f" in a ring: {atom.IsInRing()}")

    # Check if all provided atomIDs exist
    for atomID in parm_triangle.auto_query_atoms:
        if not atomID in [atom.GetIdx()+1 for atom in mol_conformations.GetAtoms()]:
            raiseError(f"FATAL: Ligand atom {atomID} not present in {parm_triangle.auto_query_input} input ligand."
                             f" Please provide correct atomIDs.\n", True)

    # for each conformation, get the bis_int_points (characterized by query_atoms) coordinates
    global conf_index
    for conf_index, conformation in enumerate(mol_conformations.GetConformers()):
        if globalParameters.debug > 0:
            fout.write(f"Ligand conformer {conf_index + 1}:\n")
        for atom_index in query_atoms.keys():
            if isHydrogen(mol_conformations.GetAtomWithIdx(atom_index)):
                raiseError(f"FATAL: Ligand atom interaction point with "
                                 f"index {mol_conformations.GetAtomWithIdx(atom_index).GetIdx()+1} is a hydrogen, which is not allowed."
                                 f"Provide the connected heavy-atom instead. In case of H-don, we automatically identify the donated hydrogen internally.\n", True)
            else:
                allconformers_bis_int_points[conf_index].append(predictBiSInteractionPoints(conformation,
                                                                                atom_index, query_atoms[atom_index]))

    # Write predicted binding site interaction points to file, for debugging mode
    if globalParameters.debug > 0:
        bis_file=utils.utilities.openFile("predicted_BiS_int.pdb", "w+")
        bis_file.write("HEADER    Predicted binding site interaction points\n")
        for nr_res, conf in enumerate(allconformers_bis_int_points):
            bis_file.write(f"MODEL {nr_res+1}\n")
            for nr, int_site in enumerate(conf):
                bis_file.write(f"ATOM   {float(nr+1):4.0f}  C   DUM A {float(nr_res+1):3.0f}     {int_site[0]:7.3f} {int_site[1]:7.3f} {int_site[2]:7.3f}\n")
            bis_file.write(f"TER\nENDMDL\n")
        fout.write("Predicted binding site interaction points written to 'predicted_BiS_int.pdb'.\n")


    # use list of coordinates of atom_index to calculate mean and variance
    allconformers_distances = [] # list of all required distances in all conformations
                        # [[conf1(dist01,dist02,dist12,...)], [conf2(dist01,dist02,dist12,...)]..]
    allconformers_distance_mean_and_var = [] # mean and varience for every distance
    vertex_combination_label = [] # labels which distances are given in allconformers_distance_mean_and_var, e.g [(0,1),(0,2),(1,2),..]
    for bis_int_points in allconformers_bis_int_points:
        _tmp_conf_labels = _calcDist(bis_int_points)
        allconformers_distances.append(_tmp_conf_labels[0])
        vertex_combination_label.append(_tmp_conf_labels[1])

    allconformers_distances = np.array(allconformers_distances)
    for distance_nr, distance_of_combination in enumerate(np.transpose(allconformers_distances)):
        allconformers_distance_mean_and_var.append(match_triangle.calcMeanAndVariance(distance_of_combination,
                                                                                      vertex_labels= vertex_combination_label[0][distance_nr],
                                                                                      parm_triangle=parm_triangle))

    if globalParameters.debug > 0:
        for combination, distance_mean_var in zip(vertex_combination_label[0], allconformers_distance_mean_and_var):
            fout.write(f"Average distance between atom {parm_triangle.auto_query_atoms[combination[0]]} and "
                       f"{parm_triangle.auto_query_atoms[combination[1]]}: {distance_mean_var[0]:.3f} "
                       f"+- {distance_mean_var[1]:.3f} Angstrom.\n")

    # only vertex_combination_label[0] since all entries (labels) are identical
    return allconformers_distance_mean_and_var, vertex_combination_label[0]


def predictBiSInteractionPoints(conformer, atom_index, atom_type):
    """
    Add offset on ligand interaction point coordinates to predict coordinates of binding site interaction point
    Only performed for H-don and H-acc, for others unaltered input coordinates are returned.
    :param atom: rdkit Atom object of the atom of interest (i.e. ligand interaction point)
    :param coords: xyz coordinates of ligand interaction point (atom)
    :param int_type: interaction type from list input_resgroups.in
    :param atom_index: atomID (starting at 0)
    :param atom_type: atom interaction type (H_don, H_acc)
    :type atom: Chem.rdchem.Atom
    :type coords: list(x,y,z)
    :type atom_index: int
    :type atom_type: str
    :return: new coordinates of protein ligand interaction point
    """
    atom = conformer.GetOwningMol().GetAtomWithIdx(atom_index)
    if atom_type == "H_don":
        if not isHdonor(atom):
            if conf_index == 0:
                raiseError(f"WARNING: Atom {atom_index + 1} is not a hydrogen donor. "
                             f"Please provide correct atom interaction type.\n", StrictMode)
                fout.write("Since strict mode is off, trying to continue... It is recommended to run this again with"
                           " 'debug: 1', and check 'predicted_BiS_int.pdb' for predicted binding site interaction points.\n")
        hydrogen_neighbors = []
        for neighbor in atom.GetNeighbors():
            if isHydrogen(neighbor):
                hydrogen_neighbors.append(neighbor)

        if len(hydrogen_neighbors) == 1:
            """ in case of a single hydrogen neighbors, add offset in direction of heavy-atom -> hydrogen """
            return calcLinearOffset(conformer, atom_index, "to_hydrogen")
        else:
            """ in case of multiple hydrogen neighbors, add offset in C-heavy_atom direction"""
            return calcLinearOffset(conformer, atom_index, "from_heavyatom")

    elif atom_type == "H_acc":
        if not isHacceptor(atom.GetOwningMol(), atom.GetIdx()):
            if conf_index == 0:
                raiseError(f"WARNING: Atom {atom_index + 1} is not a hydrogen acceptor. "
                             f"Please provide correct atom interaction type.\n", StrictMode)
                fout.write("Since strict mode is off, trying to continue... It is recommended to run this again with"
                           " 'debug: 1', and check 'predicted_BiS_int.pdb' for predicted binding site interaction points.\n")
        elif isDoubleBoundHeavyAtom(atom):
            return calcLinearOffset(conformer, atom_index, "from_heavyatom")
        else:
            return calcLonepairOffset(conformer, atom_index)

    elif atom_type in ["pos_charged", "neg_charged"]:
        """spheric offset"""
        if (atom.GetFormalCharge() >= 0) and atom_type == "neg_charged":
            if conf_index == 0:
                raiseError(f"WARNING: Atom {atom_index + 1} is defined as neg_charged, but formal charge on this atom is {atom.GetFormalCharge()}.\n", StrictMode)
                fout.write("Since strict mode is off, trying to continue... It is recommended to run this again with"
                           " 'debug: 1', and check 'predicted_BiS_int.pdb' for predicted binding site interaction points.\n")
        elif (atom.GetFormalCharge() <= 0) and atom_type == "pos_charged":
            if conf_index == 0:
                raiseError(f"WARNING: Atom {atom_index + 1} is defined as pos_charged, but formal charge on this atom is {atom.GetFormalCharge()}.\n", StrictMode)
                fout.write("Since strict mode is off, trying to continue... It is recommended to run this again with"
                           " 'debug: 1', and check 'predicted_BiS_int.pdb' for predicted binding site interaction points.\n")
        return calcSphericOffset(conformer, atom_index)

    elif atom_type == "hydrophobic":
        """spheric offset"""
        return calcSphericOffset(conformer, atom_index)

    elif atom_type in ["aromatic_ring", "hydrophobic_ring"]:
        ring_info = conformer.GetOwningMol().GetRingInfo()
        ring_found = False
        for ring_atoms in ring_info.AtomRings():
            if atom_index in ring_atoms:
                if ring_found:
                    raiseError(f"FATAL: Atom {atom_index + 1} is defined as 'aromatic ring' or 'hydrophobic ring',"
                                     f" but this atom was found in multiple rings.\n"
                                     f"Please provide an atom which is only present in a single ring.\n", True)
                else:
                    ring_found = True
        if not ring_found:
            raiseError(f"FATAL: Atom {atom_index + 1} is defined as 'aromatic_ring', but no aromatic ring containing this atom was found.\n", True)
        else:
            return calcRingCenterOffset(conformer, atom_index, ring_atoms)


def _calcDist(bis_int_points):
    """
    calculate the distances between given atoms in triangle, while keeping track of the positions
    :param bis_int_points: list of bis interaction points to calculate distances between
    :type bis_int_points: List[List[float,float,float],..]
    :return: distance
    :rtype: List[float,...], List[float,...]
    """
    # enumerate(bis_int_points) = [(0, [x,y,z]), (1, [x,y,z]), (2, [x,y,z])]
    # itertools.combinations(...) = [((0, [x,y,z]), (1, [x,y,z]) ), ... ]
    distances = []
    vertex_combination_label = []
    for atom_combination in itertools.combinations(enumerate(bis_int_points), 2):
        # atom_combination = [ ((0, [x,y,z]), (1, [x,y,z])) , ...]
        vertex_combination = (atom_combination[0][0], atom_combination[1][0])
        distance = np.linalg.norm([atom_combination[0][1][axis] - atom_combination[1][1][axis] for axis in range(3)])
        distances.append(distance)
        vertex_combination_label.append(vertex_combination)
    return distances, vertex_combination_label


def _setDistanceAndOffset(parm_triangle, mean_and_var, vertex_combination_label):
    """
    sets distances and offsets in parameter class using the calculated mean and variance of the conformational search results,
    :param mean_and_var: (position combination, (mean, var))
    :param vertex_combination_label: labels which distances are given in mean_and_var, e.g [(0,1),(0,2),(1,2),..]
    :type parm_triangle: parsers.parse_inputFile.MatchTriangleParameter
    :type mean_and_var: List[Tuple[float,float],...]
    :type vertex_combination_label: List[Tuple[int,int],...]
    """
    distances = {}
    for index, entry in enumerate(mean_and_var):
        position_no = vertex_combination_label[index]
        parm_triangle.setDistancesManually(entry[0], position_no[0], position_no[1])
        parm_triangle.setDistanceOffset(parm_triangle.auto_query_offset_factor * max(entry[1], parm_triangle.auto_query_min_offset),
                                        position_no[0], position_no[1])
    return 0


def setProteinAllowedResidues(parm_triangle, parm_general):
    """
    Takes ligand atom type (from parm_triangle object) and converts towards protein interaction types, which
    are saved in parm_general.allowed_residues
    :type parm_triangle: parsers.parse_inputFile.MatchTriangleParameter
    :type parm_general: parsers.parse_inputFile.MatchGlobalParameter
    """
    convert_ligandInttype_to_proteinInttype = {
        "H_acc": "H_don",
        "H_don": "H_acc",
        "pos_charged": "neg_charged",
        "neg_charged": "pos_charged",
        "aromatic_ring": "aromatic",
        "hydrophobic_ring": "hydrophobic",
        "hydrophobic": "hydrophobic"
    }
    if len(parm_general.aminoacid_inputs) == 0: #if this was already performed, e.g. complete_aqd mode: skip
        protein_types = []
        for ligand_type in parm_triangle.auto_query_types:
            protein_types.append(parm_general._residue_groups[convert_ligandInttype_to_proteinInttype[ligand_type]])
            parm_general.setAminoacidInputs(convert_ligandInttype_to_proteinInttype[ligand_type])
        parseAminoAcids(parm_general, protein_types)
    return parm_general
