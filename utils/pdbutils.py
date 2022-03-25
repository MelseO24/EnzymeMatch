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

def get_nrAtoms(entity):
    """
    Returns number of atoms in structure or residue
    :type entity: Bio.PDB.Entity.Entity or Bio.PDB.Structure.Structure
    :param entity: Residue or Structure object of Bio.PDB
    :return: number of atoms in residue or structure
    :rtype int
    """
    if entity.get_level() == "S":      # structure given as input
        return len([atom for atom in entity.get_atoms()])
    elif entity.get_level() == "R":    # residue given as input
        return len([atom for atom in entity])

def get_nrResidues(entity):
    """
    Returns number of atoms in entity
    :type entity: Bio.PDB.Entity.Entity or Bio.PDB.Structure.Structure
    :param entity: Chain or Structure object of Bio.PDB
    :return: number of residues in chain or structure
    :rtype int
    """
    if entity.get_level() == "S":      # structure given as input
        return len([residue for residue in entity.get_residues()])
    elif entity.get_level() == "C":    # chain given as input
        return len([residue for residue in entity])

def get_nrModels(structure):
    """
    :type structure: Bio.PDB.Structure.Structure
    :param structure: Bio.PDB structure
    :return: number of models in structure
    :rtype int
    """
    return len([model for model in structure.get_models()])

def pdbiterator(structure):
    """
    Iterates over all atoms in the PDB file, providing chainID, resID, atomID
    :type structure: Bio.PDB.Structure.Structure
    :param structure: Bio.PDB structure
    :return: generator
    """
    for atom in structure.get_atoms():
        yield {"chainID": atom.parent.parent.get_id(), "resName": atom.parent.get_resname(), "atomID": atom.get_name(),
               "coordinates": atom.get_coord(), "resID": atom.parent.get_id()[1]
               }

#def get_coordinates(structure, residueNr, atomNr):

# print(pdb_ligand[0]["A"][0]["H_PO4", 453, " "].get_name())
# for atom in pdb_ligand.get_atoms():
#    print(atom.get_name())
#    print(atom.get_coord())
#    print(atom.get_full_id())
