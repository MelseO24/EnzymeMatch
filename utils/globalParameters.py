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

# File to save all GLOBAL PARAMETERS
# Only for parameters which are often used in many modules, and which control the normal behavior of the program
import EnzymeMatch
import os.path
import sys
fout = sys.stdout

overwrite = False
path_to_main = os.path.dirname(EnzymeMatch.__file__)
workdir = None
debug = 0
strict_mode = False
verbosity = 1
parsing_counter = 0
# For use in the parsing
complete_mode = False

def setCompleteMode():
    global complete_mode
    complete_mode = True

def setWorkdir(wdir):
    if os.path.isdir(wdir):
        global workdir
        workdir = wdir
    else:
        sys.stderr.write(f"FATAL: {wdir} is not a valid directory.\n")
        exit(1)

def setDebug(setting):
    global debug
    debug = setting

def setVerbosity(verbosity_value):
    global verbosity
    verbosity = verbosity_value

def setStrictMode(setting):
    global strict_mode
    strict_mode = setting