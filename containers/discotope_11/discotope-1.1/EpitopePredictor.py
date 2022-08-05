
"""
This class will perform the actual prediction of the
antibody binding potential of each amino-acid in
the protein.
"""

__author__ = 'Nicholas Gauthier'
__version__ = '1.1'
__date__ = '2006/10/06'


from CAlphaResidue import CAlphaResidue

from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.NeighborSearch import NeighborSearch 


"""
GLOBAL
"""
pdbParser = PDBParser()


"""
METHODS
"""
oddsRatios = { 
	'ASP': 0.691,
	'GLU': 0.346,
	'ASN': 1.242,
	'SER': -0.145,
	'GLN': 1.082,
	'GLY': 0.189,
	'LYS': 1.136,
	'THR': -0.233,
	'ARG': 1.180, 
	'PRO': 1.164,
	'HIS': 1.098,
	'CYS': -3.519,
	'ALA': -1.522,
	'TYR': 0.030,
	'VAL': -1.474,
	'MET': 0.273,
	'ILE': -0.713,
	'PHE': -1.147,
	'LEU': -1.836,
	'TRP': -0.064
}
		   
def getLogOddsRatio(residue):
	return oddsRatios[residue.get_resname()]

# this method checks whether or not the residue has a C-Alpha atom
# and that it is not a calcium atom of a hetero atom
def isValidCAlphaResidue(residue):
	if residue.has_id("CA") and residue.get_id()[0].strip() == "":
		return True
	return False

class EpitopePredictor:
	# visibleName: a name the user can see (best is pdb id)
	# pdbEntry: should be a StringIO object containing the
	#    pdb file.  
	# chainsToUse: should be an array of the chains to predict
	#    predict on...or omit if prediction should be on the
	#    entire file
	def __init__(self, visibleName, pdbEntry, chains = [], verbose = False):
		self.verbose = verbose
		self.visibleName = visibleName
		self.pdbEntry = pdbEntry
		self.chainsToUse = []
		for chain in chains:
			self.chainsToUse.append(chain.upper()) #put all chains into uppercase

	def getVisibleName(self):
		return self.visibleName

	def getPdbEntry(self):
		return self.pdbEntry
	
	def getChainsToUse(self):
		return self.chainsToUse
	
	def getChainsToUseFormatted(self):
		if len(self.chainsToUse) == 0:
			return "Looking at All Chains"
		toReturn = ""
		for chain in self.chainsToUse:
			if toReturn != "":
				toReturn += ", "
			toReturn += chain
		return toReturn
	
	
	def getShouldExamineChain(self, chainid):
		if len(self.chainsToUse) == 0:
			return True
		for chain in self.chainsToUse:
			if chain.upper() == chainid.upper(): #ignore case
				return True
		return False

	# print the output as verbose if the verbose flag is set
	# to true.  otherwise don't do anything
	def printVerbose(self, str):
		if self.verbose == True:
			print str

	def predictAntibodyEpitopes(self, contactDistance, threshold):

		structure_id = self.visibleName
		structure = pdbParser.get_structure(structure_id, self.pdbEntry)

		residue_list = []
		ca_list = []
		
		allResiduesHash = {}

		# iterate over structures/models/chains/residues
		self.printVerbose("Found %i models in the structure" % len(structure.get_list()))
		
		# take the first model not matter what.  This will get rid of the
		# NMR case where there is more then one model
		model = structure.get_list()[0]
		self.printVerbose("\tA) Found " + str(len(model.get_list())) + " chains in model '" + str(model) + "'")
		for chain in model.get_list():
			if self.getShouldExamineChain(chain.get_id()):
				self.printVerbose("\t\t-Examining chain '"+str(chain)+"'...it contains "+str(len(chain.get_list()))+" residues")
				numberCA = 0

				for i in range(0, len(chain.get_list())):
					residue = chain.get_list()[i]
					if isValidCAlphaResidue(residue):
						ca_list.append(residue["CA"])
						
						self.printVerbose("\t\t\t\tResidue "+str(residue.get_id()[1])+" has a C-Alpha")
						
						#caluculate log-odds-ratio
						numberNeighbors = 0
						propSmooth = 0
						for j in range(i-4, i+5):
							if j >= 0 and j < len(chain.get_list()):
								neighborRes = chain.get_list()[j]
								if isValidCAlphaResidue(neighborRes):
									propSmooth = propSmooth + getLogOddsRatio(neighborRes)
									numberNeighbors = numberNeighbors + 1
						propSmooth = propSmooth / numberNeighbors
						caResidue = CAlphaResidue(residue, propSmooth)

						residue_list.append(caResidue)
						allResiduesHash[residue] = caResidue
						numberCA = numberCA + 1
					else:
						self.printVerbose("\t\t\t\tResidue "+str(residue.get_id()[1])+" doesn't have a C-Alpha (not predicting on)")
						
				self.printVerbose("\t\t\t-"+str(numberCA)+" residues in chain "+str(chain)+" are C-Alphas")
			else:
				self.printVerbose("\t\t-Not examining chain '"+str(chain)+"' (user choose only other chains)")
				
		self.printVerbose("In total "+str(len(ca_list))+" C-Alpha atoms were found")
		
		if len(ca_list) < 1:
			print "DiscoTope predictions for '" + self.visibleName + "'."
			print "\tLooking only at Chain: ", self.getChainsToUseFormatted()
			print "\tContact Distance = %.3f Angstroms\n\tThreshold = %.3f\n" % (contactDistance, threshold)
			print "Error! No C-Alpha residues found in chain(s) " + self.getChainsToUseFormatted()
		else:
			# Search for neighbors and add to each residue as a neighbor
			neighborSearch = NeighborSearch(ca_list)

			mycount = 0
			for residue in residue_list:
				closeResidues = neighborSearch.search(residue.getResidue()["CA"].coord, contactDistance, "R")
				nameChain = str(residue.getResidue().get_id()[1]) + " in chain " + str(residue.getResidue().get_parent().get_id())
				self.printVerbose("Found "+str(len(closeResidues))+" residues within "+str(contactDistance)+" Angstroms of "+nameChain)
				for neighbor in closeResidues:
					try:
						cAlphaResidue = allResiduesHash[neighbor]
						residue.addNeighbor(cAlphaResidue)
					except KeyError:
						#its ok, it just wasn't necessarly a c-alpha residue
						pass

			# output results
			print "DiscoTope predictions for '" + self.visibleName + "'."
			print "\tLooking only at Chain: ", self.getChainsToUseFormatted()
			print "\tContact Distance = %.3f Angstroms\n\tThreshold = %.3f\n" % (contactDistance, threshold)
			#print "chain\tseq #\ttype\tcontact.num\tprox.sum\toverall.score\n"
			numberResiduesPredicted = 0
			numberResidues = 0
			for residue in residue_list:
				toPrint = str(residue.getResidue().get_parent().get_id()) + "\t"
				toPrint += str(residue.getResidue().get_id()[1]) + str(residue.getResidue().get_id()[2]).strip() + "\t"
				toPrint += residue.getResidue().get_resname() + "\t"
				toPrint += str(len(residue.getCAlphaNeighbors())) + "\t"
				toPrint += "%.3f" % residue.getProxSum() + "\t"

				overallscore = residue.getOverallScore()
				toPrint += "%.3f" % overallscore
				if overallscore >= threshold:
					toPrint += "\t<=B"
					numberResiduesPredicted += 1
				print toPrint
			print "\nIdentified "+str(numberResiduesPredicted)+" B-Cell epitope residues out of "+str(len(residue_list))+" total residues"
			
