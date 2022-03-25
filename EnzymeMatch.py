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


import os
import sys
from datetime import datetime

fout = sys.stdout
ferr = sys.stderr

def runEnzymeMatch(input_file):
    """
    Runs EnzymeMatch, based on the mode defined by the user
    # :param run_mode: string defining mode, either: residue_match/triangle_match/both
    :param input_file: path to input file
    """
    from parsers.inputParser_main import parseParameter
    from match_triangle import matchTriangle
    from match_residue import matchResidue
    from AQD import runAQDstandalone, setProteinAllowedResidues

    time_start = datetime.now()
    fout.write(f"Execution started at {time_start.strftime('%d/%m/%Y %H:%M:%S')}\n")

    if not os.path.isfile(input_file):
        ferr.write(f"FATAL: {input_file} does not exist.\n")
        sys.exit(1)
    fout.write("Starting parsing parameters...\n")
    run_mode, gen_options, residue_options, triangle_options = parseParameter(input_file)
    fout.write("Finished parsing parameters...\n\n")
    fout.write(f"EnzymeMatch was started in mode: {run_mode}\n")

    if run_mode.lower() == "matchresidue":
        matchResidue(gen_options, residue_options)
    elif run_mode.lower() == "matchtriangle":
        matchTriangle(gen_options, triangle_options)
    elif run_mode.lower() == "complete":
        # complete mode protocols time needed for matchResidue and matchTriangle separately
        if triangle_options.auto_query_input:
            setProteinAllowedResidues(triangle_options, gen_options) #set >Aminoacids
        matchResidue(gen_options, residue_options)
        time_now = datetime.now()
        intermediate_runtime = round((time_now-time_start).total_seconds())
        print(f"Runtime: {intermediate_runtime} seconds\n")
        matchTriangle(gen_options, triangle_options)
    elif run_mode.lower() == "aqd":
        runAQDstandalone(triangle_options, gen_options)
    else:
        ferr.write("FATAL: something went wrong, no mode identified...\nChoose from: 'matchresidue', 'matchtriangle', 'complete', 'AQD'\n")
        sys.exit(1)
    fout.write(f"EnzymeMatch finished without errors in {run_mode} mode.\n")
    time_end = datetime.now()
    fout.write(f"Execution ended at {time_end.strftime('%d/%m/%Y %H:%M:%S')}\n")
    runtime = round((time_end - time_start).total_seconds())
    print(f"Runtime: {runtime} seconds\n")


if __name__ == "__main__":
    import argparse
    from utils import globalParameters

    parser = argparse.ArgumentParser(description="BioLiP, find binding sites binding your ligand.")
    parser.add_argument('-i', help="Input file", required=True)
    parser.add_argument("-O", help="allow overwrite of output file", default=False, action="store_true")
    args = parser.parse_args()

    globalParameters.overwrite = args.O
    runEnzymeMatch(input_file=args.i)
