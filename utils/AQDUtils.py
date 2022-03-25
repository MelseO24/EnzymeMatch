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
import numpy as np
from utils import globalParameters
from utils.utilities import raiseError
from rdkit import Chem
from rdkit.Chem import Lipinski
from scipy.spatial.transform import Rotation as scipy_rot
fout = sys.stdout

# definition of functional groups:
fg_carbonyl = Chem.MolFromSmarts("[CX3]=[OX1]")
fg_phoshponyl = Chem.MolFromSmarts("[P*]=[OX1]")
fg_sulfonyl = Chem.MolFromSmarts("[S*]=[OX1]")

def nrOfBonds(atom):
    return atom.GetDegree()

def containsHydrogens(molecule):
    """
    Checks if hydrogens are present in input molecule
    :param molecule: molecule to check for presence of hydrogens
    :type molecule: Chem.rdchem.Mol
    :return: bool
    """
    H_found = False
    for atom in molecule.GetAtoms():
        if isHydrogen(atom):
            H_found = True
            break
    if not H_found:
        raiseError("WARNING: input ligand does not contain any hydrogens, please check if they the ligand "
                         "is correctly protonated\n", False)
        return False
    else:
        return True

def isHydrogen(atom):
    """
    Returns if atom is a hydrogen by checking it's atomic number
    :param atom: atom of interest
    :type atom: Chem.rdchem.Atom
    :return: bool
    """
    if atom.GetAtomicNum() == 1:
        return True
    else:
        return False

def isDoubleBoundHeavyAtom(atom):
    """
    Checks if the atom given in atomID is the heavy atom of a carbonyl group, or P=O or S=O
    :param atom: rdkit Atom object of the atom of interest (i.e. ligand interaction point)
    :type atom: Chem.rdchem.Atom
    :return: bool if atom is the carbonyl oxygen
    """
    if (atom.GetAtomicNum() == 8) and (len(atom.GetNeighbors()) == 1):
        #Carbonyl does not always work with SMARTS (especially when in aromatic rings), thus manually:
        if len(atom.GetBonds()) == 1 \
                and atom.GetBonds()[0].GetBondType() == Chem.rdchem.BondType.DOUBLE \
                and atom.GetNeighbors()[0].GetAtomicNum() == 6:
            return True
        #if {atom.GetIdx(), atom.GetNeighbors()[0].GetIdx()} in \
        #        [set(x) for x in atom.GetOwningMol().GetSubstructMatches(fg_carbonyl)]:  # sets since they are unordered
        #    return True
        if {atom.GetIdx(), atom.GetNeighbors()[0].GetIdx()} in \
                [set(x) for x in atom.GetOwningMol().GetSubstructMatches(fg_phoshponyl)]:  # sets since they are unordered
            return True
        elif {atom.GetIdx(), atom.GetNeighbors()[0].GetIdx()} in \
                [set(x) for x in atom.GetOwningMol().GetSubstructMatches(fg_sulfonyl)]:  # sets since they are unordered
            return True
    return False

def isHacceptor(molecule, atomID):
    return atomID in [x[0] for x in Lipinski._HAcceptors(molecule)]

def isHdonor(atom):
    """
    Returns if the atom is a hydrogen bond donor.
    Note: does explicitly NOT consider tautomerization!
    :param atom: atom which should be checked
    :type atom: Chem.rdchem.Atom
    :return: bool
    """
    if atom.GetIdx() in [x[0] for x in Lipinski._HDonors(atom.GetOwningMol())]:
        # Do not consider tautomerization, thus needs to be connected to hydrogen
        for neighbor in atom.GetNeighbors():
            if isHydrogen(neighbor):
                return True
    return False

def calcLinearOffset(conformer, heavyatom_index, int_type):
    """
    Calculates where to place the offset on the hydrogen donor ([O/N/S]-H bond, thus linear) to predict position of the
    binding site interaction point.
    :param conformer: Chem.rdchem.Conformer
    :param heavyatom_index: atomID of heavy atom (starting from 0)
    :param int_type: controls if offset is build on H_don or carbonyl, allowed arguments: ["to_hydrogen", "from_heavyatom"]
    :return: coordinates of binding site interaction point (i.e. coords with offset)
    :rtype: List[float, float, float]
    """
    molecule = conformer.GetOwningMol()
    heavy_atom = molecule.GetAtomWithIdx(heavyatom_index)
    neighbors = heavy_atom.GetNeighbors()
    if int_type not in ["to_hydrogen","from_heavyatom"]:
        raiseError("calcLinearOffset function, argument 'int_type' only takes arguments: ['to_hydrogen', 'from_heavyatom']\n", True)

    if int_type == "to_hydrogen":
        """Offset is placed in direction of heavy_atom -> hydrogen bond"""
        hydrogen_neighbors = []
        # Search for neighboring Hydrogens ('to donate')
        for neighbor in neighbors:
            if isHydrogen(neighbor):
                hydrogen_neighbors.append(neighbor)
        if len(hydrogen_neighbors) == 1:
            linear_offset = _addLinearOffset(conformer.GetAtomPosition(heavyatom_index),
                                             conformer.GetAtomPosition(hydrogen_neighbors[0].GetIdx()),
                                             1.9,
                                             int_type)
                                              # distance from  Boehm(Journal of Computer-Aided Molecular Design, 6 (1992) 593406.
        else:
            raiseError(
                f"Multiple hydrogens bound to heavy_atom {heavy_atom.GetIdx() + 1}, this cannot be handled yet...\n", True)

    if int_type == "from_heavyatom":
        """Offset is placed in direction of non-hydrogen-neighbor -> heavy-atom, since multiple hydrogen atoms are present and only
        one offset is allowed"""
        non_hydrogen_neighbors = [neighbor for neighbor in neighbors if not isHydrogen(neighbor)]
        if not len(non_hydrogen_neighbors) == 1:
            raiseError(f"FATAL: Invalid atomtype for atom {heavyatom_index + 1}.\n", True)
        else:
            linear_offset = _addLinearOffset(conformer.GetAtomPosition(heavyatom_index),
                                             conformer.GetAtomPosition(non_hydrogen_neighbors[0].GetIdx()),
                                             1.9,
                                             int_type)
                                              # distance from  Boehm-Journal of Computer-Aided Molecular Design, 6 (1992) 593406.

    if globalParameters.debug > 0:
        fout.write(f"For ligand atom {heavyatom_index+1} ({int_type}), the binding site "
                   f"interaction point is predicted at coordinates:"
                   f" [{linear_offset[0]:.3f}, {linear_offset[1]:.3f}, {linear_offset[2]:.3f}]\n")
    return linear_offset

def calcLonepairOffset(conformer, heavyatom_index):
    """
    Calculates where to place the offset on the H-acc C-O-C, or C=N-C, in the direction of the lone pair of the heavy-atom
    to predict the position of the binding site interaction point.
    :param conformer: Chem.rdchem.Conformer
    :param heavyatom_index: atomID of heavy heavy_atom (starting from 0)
    :return: coordinates of binding site interaction point (i.e. coords with offset)
    :rtype: List[float, float, float]
    """
    molecule = conformer.GetOwningMol()
    heavy_atom = molecule.GetAtomWithIdx(heavyatom_index)
    neighbors = heavy_atom.GetNeighbors()

    if not len(neighbors) == 2:
        raiseError(f"Illegal: Atom {heavy_atom.GetIdx() + 1} does not have 2 neighbors, but defined as a non-carbonyl "
                         f"hydrogen acceptor.\n", True)
    else:
        lonepair_offset = _addLonepairOffest(conformer.GetAtomPosition(heavyatom_index),
                           conformer.GetAtomPosition(neighbors[0].GetIdx()),
                           conformer.GetAtomPosition(neighbors[1].GetIdx()),
                           1.9)
        if globalParameters.debug > 0:
            fout.write(f"For ligand atom {heavyatom_index+1} (H_acc), the binding site "
                   f"interaction point is predicted at coordinates:"
                   f" [{lonepair_offset[0]:.3f}, {lonepair_offset[1]:.3f}, {lonepair_offset[2]:.3f}]\n")
    return lonepair_offset

def calcSphericOffset(conformer, atom_index):
    """
    Calculates a spheric offset around the provided atom. The center will be on the atom itself
    :param conformer: Chem.rdchem.Conformer
    :param atom_index: atomID of atom (starting from 0)
    :return: coordinates of binding site interaction point (i.e. coords with offset)
    :rtype: List[float, float, float]
    """
    atom = conformer.GetOwningMol().GetAtomWithIdx(atom_index)
    spheric_offset = conformer.GetAtomPosition(atom_index) # center of atom, since spheric offset

    if globalParameters.debug > 0:
        fout.write(f"For ligand atom {atom_index + 1} (spheric offset), the binding site "
               f"interaction point is predicted at coordinates:"
               f" [{spheric_offset[0]:.3f}, {spheric_offset[1]:.3f}, {spheric_offset[2]:.3f}]\n")
    return list(spheric_offset)

def calcRingCenterOffset(conformer, atom_index, ring_atoms):
    """
    Calculates the center of the aromatic ring
    :param conformer: Chem.rdchem.Conformer
    :param atom_index: selected atomID by user
    :param ring_atoms: atomIDs of atoms forming a ring
    :type conformer: Chem.rdchem.Conformer
    :type atom_index: int
    :type ring_atoms: List[int,...]
    :return: coordinates of binding site interaction point (i.e. center of aromatic ring)
    :rtype: List[float, float, float]
    """
    ring_center = []
    for axis in range(3):
        ring_center.append(np.average([conformer.GetAtomPosition(ring_atom)[axis] for ring_atom in ring_atoms]))

    if globalParameters.debug > 0:
        fout.write(f"For ligand atom {atom_index + 1} (aromatic ring), the binding site "
               f"interaction point (i.e. ring center) is predicted at coordinates:"
               f" [{ring_center[0]:.3f}, {ring_center[1]:.3f}, {ring_center[2]:.3f}]\n")
    return ring_center

def _addLinearOffset(coords_heavyatom, coords_linked_atom, offset_length, interaction_type):
    """
    Adds offset in direction of vector coords_heavyatom -> coords_hydrogen of length offset_length.
    Intended for internal use only.
    :param coords_heavyatom: coordinates of heavy atom
    :param coords_linked_atom: coordinates of hydrogen atom in case of H_don, or carbonyl carbon for H_acc carbonyl
    :param offset_length: distance from hydrogen to create offset
    :param interaction_type: H_don or carbonyl (ie. a H_acc)
    :type coords_heavyatom: List[float, float, float]
    :type coords_linked_atom: List[float, float, float]
    :type offset_length: float
    :type interaction_type: str
    :return: new coordinates of binding site interaction point (i.e. with offset)
    :rtype: List[float, float, float]
    """
    init_vector = [coords_linked_atom[axis] - coords_heavyatom[axis] for axis in range(3)]
    offset_vector = _convertVectorLength(init_vector, offset_length)
    if interaction_type == "to_hydrogen":
        return [coords_linked_atom[axis] + offset_vector[axis] for axis in range(3)]
    elif interaction_type == "from_heavyatom":
        # sign flip since direction of offset vector is turned (C=O) is opposite bond then (O-H)
        return [coords_heavyatom[axis] - offset_vector[axis] for axis in range(3)]

def _addLonepairOffest(coords_heavyatom, coords_neighbor1, coords_neighbor2, offset_length):
    """
    Adds offset in direction of vector coords_heavyatom -> lone pair of length offset_length.
    Intended for internal use only.
    :param coords_heavyatom: coordinates of heavy atom (the one with lone-pair)
    :param coords_neighbor1: coordinates of the first neighbor of the heavy-atom
    :param coords_neighbor2: coordinates of the second neighbor of the heavy-atom
    :param offset_length: distance from heavy-atom to create offset
    :type coords_heavyatom: List[float,float,float]
    :type coords_neighbor1: List[float,float,float]
    :type coords_neighbor2: List[float,float,float]
    :return: new coordinates of binding site interaction point (i.e. with offset)
    :rtype: List[float,float,float]
    """
    vector1 = [coords_neighbor1[axis] - coords_heavyatom[axis] for axis in range(3)]  # vector between heavy-atom and neighbor 1
    vector2 = [coords_neighbor2[axis] - coords_heavyatom[axis] for axis in range(3)]  # vector between heavy-atom and neighbor 2

    # calculate angle beteween vector 1 and vector 2:
    # using cos(teta) = (a.b)/(|a||b|)  (note, (a.b) is dot-product, and |a|,|b| respect. are vector lengths!)
    angle_v1v2_rads = np.arccos(np.dot(vector1, vector2) / (_calcVectorLength(vector1)*_calcVectorLength(vector2)))

    # calculate axis of rotation (i.e. axis orthogonal on vector1 and vector2, using cross product: v1xv2
    axis_of_rotation = np.cross(vector1, vector2)
    # length of rotation axis is the amount (in radians) of rotation, thus needs to be rescaled
    axis_of_rotation = _convertVectorLength(list(axis_of_rotation), ((2*np.pi-angle_v1v2_rads)/2))
    lonepair_rotation = scipy_rot.from_rotvec(axis_of_rotation)   # scipy rotation function object
    lonepair_vector = lonepair_rotation.apply(vector2)            # the vector from heavy-atom to lone-pair, i.e. the vector to place offset on
    lonepair_vector = _convertVectorLength(lonepair_vector, offset_length)
    return [coords_heavyatom[axis] + lonepair_vector[axis] for axis in range(3)]

def _calcVectorLength(vector):
    """
    Calculates length of input vector
    :param vector: input vector
    :type vector: List[float, float, float]
    :return: length of vector
    :rtype: float
    """
    return np.linalg.norm(vector)

def _convertVectorLength(vector, target_length):
    """
    Returns a rescaled vector with target_length
    :param vector: vector to resize
    :param target_length: target length of vector
    :type vector: List[float, float, float]
    :type target_length: float
    :return: rescaled vector
    :rtype: List[float, float, float]
    """
    init_length = _calcVectorLength(vector)
    scale_factor = target_length / init_length
    return [vector[axis]*scale_factor for axis in range(3)]