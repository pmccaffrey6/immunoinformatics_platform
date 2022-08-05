#! /usr/local/python/bin/python

"""
This class represents a C-alpha residue
"""

__author__ = 'Nicholas Gauthier'
__version__ = '1.1'
__date__ = '2006/10/10'

class CAlphaResidue:
	
	def __init__(self, residue, propSmooth):
		self.residue = residue
		self.neighbors = []
		self.propSmooth = propSmooth
		self.overallScore = 0

	def residueName(self):
		return self.residue.get_resname()
		
	def getResidue(self):
		return self.residue
		
	def getProxSum(self):
		proxSum = 0
		for neighbor in self.neighbors:
			proxSum = proxSum + neighbor.propSmooth
		return proxSum


	# the overall score is defined as:
	#    prox.sum -0.5*contact number
	def getOverallScore(self):
		#overallScore = self.proxSum - (len(self.getCAlphaNeighbors())/2) 
		
		# The original program always rounds UP the odd number of contacts (ie. if
		# the contact number is 13 and 13/2 = 6.5, the original program subtracted
		# 7 from the proxSum.  This may be undesirable, but the initial release will
		# work the same, rounding up those numbers...
		contactsWithWeight = len(self.getCAlphaNeighbors())/2.0
		#print "overall:", int((len(self.getCAlphaNeighbors())/2.0)+.00001)
		overallScore = self.getProxSum() - contactsWithWeight
		return overallScore
		
		
	def addNeighbor(self, neighborResidue):
		self.neighbors.append(neighborResidue)
	
	def setNeighbors(self, newNeighbors):
		self.neighbors = newNeighbors
	
	
	def getCAlphaNeighbors(self):
		toReturn = []
		for neighbor in self.neighbors:
			if(neighbor.getResidue().has_id("CA")):
				toReturn.append(neighbor)
		return toReturn
