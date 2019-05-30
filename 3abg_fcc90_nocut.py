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
branchTrack = treeReader.UseBranch("Track")

#testmasses to scan
#testmasses = [10, 20, 30, 40, 50, 60, 70, 80, 90]
testmasses = [80]
#ecm
ecm = 45.8*2 #GeV

#ptot
#ptot = ROOT.TLorentzVector(0.,0.,0.,91.6)

#recreate root output file
out_root = ROOT.TFile(outputFile,"RECREATE")	

#Book histograms
histPhotonEnergy = ROOT.TH1F("E_photons", "Energy [GeV]", 150, 0.0, 100.0)
histALPMass2D = ROOT.TH2F("Mass_2D", "Mass [GeV]", 150, 0.0, 100.0, 90, 10., 90.)

#loop over testmasses to scan
for counter, testmass in enumerate(testmasses):

	#Book histograms inside loop
	histALPMass = ROOT.TH1F("M_alp_"+str(testmass), "Mass [GeV]", 150, 0.0, 100.0) #mass reconstructed with 

	# Loop over all events
	for entry in range(0, numberOfEntries):
		# Load selected branches with data from specified event
		treeReader.ReadEntry(entry)

		threephotons_vec = [] #array with 3 TLorentzVector from 3 photons

		#use only events with 3 photons
		if branchPhoton.GetEntries() != 3:
			continue
		if branchTrack.GetEntries() != 0:
			continue

		#add in histo all photon energies
		for i in range(0, branchPhoton.GetEntries()):

			photon = branchPhoton.At(i)
			Photon_TLV = photon.P4()                  #TLorentzVector
			if counter == 0:
				histPhotonEnergy.Fill(Photon_TLV.E())
	
			threephotons_vec.append(Photon_TLV)
		
		photons_vec_added = threephotons_vec[0]+threephotons_vec[1]+threephotons_vec[2]
		if photons_vec_added.E() < 70.:
			print "continued"
			continue

		#find two photons 
		ipalp1, ipalp2, imind, egtest = includeme.ass_3l(threephotons_vec, ecm, testmass)
		
		if entry % 100 == 0:
			print testmass, entry, ipalp1, ipalp2, imind, egtest
			
		vec_alpphotons = threephotons_vec[ipalp1]+threephotons_vec[ipalp2]
		histALPMass.Fill(vec_alpphotons.M()) #fill mass of alp
		histALPMass2D.Fill(egtest, testmass)

	histALPMass.Write()

histPhotonEnergy.Write()
histALPMass2D.Write()

