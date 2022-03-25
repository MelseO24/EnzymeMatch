# EnzymeMatch

    Copyright (C) 2022  Okke Melse

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
    
**Modes:**

- CompleteMatch: Complete run of EnzymeMatch, consecutively conducting AQD, ResidueMatch and TriangleMatch.
- Automatic Query Design (AQD): Stand-alone Automatic Query Design mode, i.e. automatic input generation for ResidueMatch and TriangleMatch based on input of substrate (.mol2 and .sdf supported).
- ResidueMatch: Match enzymes from BioLiP database containing complement residue types with respect to the target substrate.
- TriangleMatch: Geometrical match of interaction point pattern between enzyme and substrate. This mode should ideally run on the output of ResidueMatch.
