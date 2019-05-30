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

testmass = raw_input("insert testmass: ")

# Book histograms
histPhotonEnergy = ROOT.TH1F("E_photons", "Energy [GeV]", 150, 0.0, 100.0)
histALPPhotonEnergy = ROOT.TH1F("E_alp_photon", "Energy [GeV]", 150, 0.0, 100.0)
histALPMass = ROOT.TH1F("M_alp", "Mass [GeV]", 150, 0.0, 100.0)
histZPhotonEnergy = ROOT.TH1F("E_photon", "Energy [GeV]", 150, 0.0, 100.0)

#ptot = ROOT.TLorentzVector(0.,0.,0.,90.)
# Loop over all events

for entry in range(0, numberOfEntries):
	# Load selected branches with data from specified event
	treeReader.ReadEntry(entry)

	threephotons_vec = [] #array with 3 TLorentzVector from 3 photons

	#use only events with 3 photons
	if branchPhoton.GetEntries() != 3:
		continue

	#add in histo all photon energies
	for i in range(0,branchPhoton.GetEntries()):

		photon = branchPhoton.At(i)
		Photon_TLV = photon.P4()                  #TLorentzVector
		histPhotonEnergy.Fill(Photon_TLV.E())
	
		threephotons_vec.append(Photon_TLV)

	#find two photons 
	ipalp1, ipalp2, imind, egtest = includeme.ass_3l(threephotons_vec, 91.6, float(testmass))
	print ipalp1, ipalp2, imind, egtest

	histALPPhotonEnergy.Fill(threephotons_vec[ipalp1].E()) #fill energy of first photon from alp
	histALPPhotonEnergy.Fill(threephotons_vec[ipalp2].E()) #fill energy of second photon from alp

	histZPhotonEnergy.Fill(threephotons_vec[imind].E()) #fill energy of photon from Z
	alp_vec = threephotons_vec[ipalp1]+threephotons_vec[ipalp2]
	histALPMass.Fill(alp_vec.M()) #fill mass of alp

out_root = ROOT.TFile(outputFile,"RECREATE")	
histPhotonEnergy.Write()
histALPMass.Write()
histALPPhotonEnergy.Write()
histZPhotonEnergy.Write()

