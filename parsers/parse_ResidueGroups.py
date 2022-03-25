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

import sys, os
from utils import utilities
import utils.globalParameters

def parseResidueGroups(filepath):
    """
    Parses input file, mode dependent on mode given in input file (match_residue or match_triangle)
    :rtype Dict[str: str,...]  # {'H_don':'R,K,W', 'H_acc':'R,H,K,D...', ...}
    """

    groups = {}
    if not os.path.isfile(filepath):
        sys.stderr.write(f"FATAL: File {filepath} does not exist.\n")
        sys.exit(1)
    fres_groups = utilities.openFile(filepath, 'r')

    # Locally parse all options in list all_input_parameters
    for line in fres_groups.readlines():
        if not (line.strip().startswith(";")) and (not len(line.strip()) == 0):
            # allow comment lines in input file, and ignore empty lines
            residue_attribute, residues = utils.utilities.parseInputLine(line)
            # cosmetics: attribute names without >
            groups[residue_attribute.replace('>', '')] = residues
    return groups
