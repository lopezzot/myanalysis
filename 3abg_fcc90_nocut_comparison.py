#!/usr/bin/env python

import sys
import ROOT
import includeme

if len(sys.argv) < 2:
	print " Usage: python examples/analysis.py input.root output.root"
	sys.exit(1)

ROOT.gSystem.Load("libDelphes")

try:
	ROOT.gInterpreter.Declare('#include "classes/SortableObject.h"')
	ROOT.gInterpreter.Declare('#include "classes/DelphesClasses.h"')
	ROOT.gInterpreter.Declare('#include "external/ExRootAnalysis/ExRootTreeReader.h"')
except:
	pass

inputFile = sys.argv[1]
outputFile = sys.argv[2]

# Create chain of root trees
chain = ROOT.TChain("Delphes")
chain.Add(inputFile)

# Create object of class ExRootTreeReader
treeReader = ROOT.ExRootTreeReader(chain)
numberOfEntries = treeReader.GetEntries()

# Get pointers to branches used in this analysis
branchPhoton = treeReader.UseBranch("Photon")

testmass = 80
#ecm
ecm = 45.8*2 #GeV

#ptot
#ptot = ROOT.TLorentzVector(0.,0.,0.,91.6)
#print "event*MALP"
# Loop over all events
for entry in range(0, treeReader.GetEntries()):
	# Load selected branches with data from specified event
	treeReader.ReadEntry(entry)

	#print "Event: "+str(entry)

	if branchPhoton.GetEntries() == 3:
		
		threephotons_vec = [] #array with 3 TLorentzVector from 3 photons

		#print "# photons: "+str(branchPhoton.GetEntries())
		
		for counter, photon in enumerate(branchPhoton):
			Photon_TLV = photon.P4()
			#print "entry "+str(entry)
			#print "photon "+str(counter)+" energy "+str(Photon_TLV.E())    
			#print "photon "+str(counter)+" pt "+str(photon.PT)
			#print "photon "+str(counter)+" eta "+str(photon.Eta)
			#print "photon "+str(counter)+" phi "+str(photon.Phi)
		
			threephotons_vec.append(Photon_TLV)
		
		#find two photons 
		ipalp1, ipalp2, imind, egtest = includeme.ass_3l(threephotons_vec, ecm, testmass)
		
		vec_alpphotons = threephotons_vec[ipalp1]+threephotons_vec[ipalp2]
		MALP = vec_alpphotons.M()
		#print str(entry)+"*"+str(vec_alpphotons.M())
		if MALP < 40.:
			print "entry "+str(entry)+" MALP "+str(MALP)		
			for counter, photon in enumerate(branchPhoton):
				Photon_TLV = photon.P4()
				
				print "photon "+str(counter)+" energy "+str(Photon_TLV.E())    
				print "photon "+str(counter)+" pt "+str(photon.PT)
				print "photon "+str(counter)+" eta "+str(photon.Eta)
				print "photon "+str(counter)+" phi "+str(photon.Phi)
			
	#else:
		#print str(entry)+"*"+"-1"
