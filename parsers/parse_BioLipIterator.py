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

import re
import utils.globalParameters
from utils.utilities import raiseError



class BioLipIterator(object):
    """
    Parses the Annotations acquired from the BioLip database. (https://zhanglab.ccmb.med.umich.edu/BioLiP/download.html)
    In every iteration, each line of the BioLip annotation file is read and saved until the next iteration.
    This iterator helps reading the file, while not having the save the whole file into the memory during run.
    The information read is used to navigate and test the proteins of eligibility.
    """
    def __init__(self, path):
        """

        :param path: The path of the biolipdata file
        :type path: String
        """
        self.path = path
        self.protchain = None
        self.LigID = None
        self.LigChain = None
        self.LigSerialNR = None
        self.BisResidues = []
        self.CatResidues = []
        self.EC = None
        self.fullLine = None

    def __iter__(self):
        with open(self.path) as f:
            for line in f.readlines():
                try:
                    self.entry_name = (str(line.split('\t')[0]) + str(line.split('\t')[1]) + "_" + str(line.split('\t')[3]))
                except IndexError:
                    raiseError(f"FATAL: An error occurred in during parsing of the BioLip database.\n"
                                     f"Please check if this entry is tab-separated, and if all BioLip entries are provided.\n"
                                     f"Line: {line}", True)
                except:
                    raiseError(f"FATAL: An unknown error has occurred during parsing of the BioLiP database.\n"
                                     f"Line: {line}", True)
                self.protchain = str(line.split('\t')[1])
                self.LigID = str(line.split('\t')[4])
                self.LigChain = str(line.split('\t')[5])
                self.LigSerialNR = str(line.split('\t')[6])
                self.BisResidues = line.split('\t')[7].split()
                self.CatResidues = re.split('[\s;]+', line.split('\t')[9])
                self.CatResidues = [i for i in self.CatResidues if i]  # remove empty strings
                self.EC = line.split('\t')[11]
                self.fullLine = line
                self.pdbDownloadName = str(line.split('\t')[0])
                yield {"entry_name": self.entry_name,
                       "protchain": self.protchain,
                       "LigID": self.LigID,
                       "LigChain": self.LigChain,
                       "LigSerialNR" : self.LigSerialNR,
                       "BisResidues" : self.BisResidues,
                       "CatResidues" : self.CatResidues,
                       "EC" : self.EC,
                       "fullLine" : self.fullLine,
                       "pdbDownloadName" : self.pdbDownloadName,
                       "ProtPDBName" : self.entry_name.split('_')[0] + ".pdb",
                       "LigPDBName" : self.entry_name.split('_')[0][:-1] + "_" + self.LigID + "_" + self.LigChain + "_" + self.LigSerialNR + ".pdb"}

    def get_pdbDownloadName(self):
        return self.pdbDownloadName

    def get_entryNames(self):
        return self.entry_name

    def get_LigID(self):
        return self.LigID

    def get_LigChain(self):
        return self.LigChain

    def get_ProtChain(self):
        return self.protchain

    def get_LigSerialNr(self):
        return self.LigSerialNR

    def get_BisResidues(self):
        return self.BisResidues

    def get_CatResidues(self):
        return self.CatResidues

    def get_fullLine(self):
        return self.fullLine

    def ProtPDBName(self):
        return self.entry_name.split('_')[0] + ".pdb"

    def get_EC(self):
        return self.EC

    def LigPDBName(self):
        self.pdbname = self.entry_name.split('_')[0][:-1]
        return self.pdbname + "_" + self.LigID + "_" + self.LigChain + "_" + self.LigSerialNR + ".pdb"