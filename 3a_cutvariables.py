#!/usr/bin/env python
#code to show histo of background and signal samples in order to understand how to apply cuts
#to be used with /home/delphesresults/3abg_fcc90_nocut for background samples
#to be used with /home/delphesresults/alpscan for signal samples
#created on 8/5/2019

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

BgFile = sys.argv[1]
SgFile = sys.argv[2]
outputFile = sys.argv[3]

# Create chain of root trees
Bgchain = ROOT.TChain("Delphes")
Bgchain.Add(BgFile)

Sgchain = ROOT.TChain("Delphes")
Sgchain.Add(SgFile)

# Create object of class ExRootTreeReader
BgtreeReader = ROOT.ExRootTreeReader(Bgchain)
BgnumberOfEntries = BgtreeReader.GetEntries()

SgtreeReader = ROOT.ExRootTreeReader(Sgchain)
SgnumberOfEntries = SgtreeReader.GetEntries()

print str(BgnumberOfEntries)+" background events\n"
print str(SgnumberOfEntries)+" signal events\n"

# Get pointers to branches used in this analysis
BgbranchPhoton = BgtreeReader.UseBranch("Photon")
BgbranchTrack = BgtreeReader.UseBranch("Track")

SgbranchPhoton = SgtreeReader.UseBranch("Photon")
SgbranchTrack = SgtreeReader.UseBranch("Track")

testmass = 80

#ecm
ecm = 45.8*2 #GeV

#ptot
#ptot = ROOT.TLorentzVector(0.,0.,0.,91.6)

#recreate root output file
out_root = ROOT.TFile(outputFile,"RECREATE")	

#Book histograms
histBgepho2epho1 = ROOT.TH1F("Bg_Epho2/Epho1", "Bg_Epho2/Epho1", 150, 0.0, 1.1)
histSgepho2epho1 = ROOT.TH1F("Sg_Epho2/Epho1", "Sg_Epho2/Epho1", 150, 0.0, 1.1)

histBgdeltar = ROOT.TH1F("Bg_DeltaR", "Bg_DeltaR", 100, 0., 5.)
histSgdeltar = ROOT.TH1F("Sg_DeltaR", "Sg_DeltaR", 100, 0., 5.)

histBgMALPcut = ROOT.TH1F("Bg_MALPcut", "Bg", 100, 0., 50.)
histSgMALPcut = ROOT.TH1F("Sg_MALPcut", "Sg", 100, 0., 50.)
# Loop over signal events
for entry in range(0, SgnumberOfEntries):
	# Load selected branches with data from specified event
	SgtreeReader.ReadEntry(entry)

	threephotons_vec = [] #array with 3 TLorentzVector from 3 photons

	#use only events with 3 photons and not tracks
	if SgbranchPhoton.GetEntries() != 3:
		continue
	if SgbranchTrack.GetEntries() != 0:
		continue

	#add lorentzvetor photons in vector
	for i in range(0, SgbranchPhoton.GetEntries()):

		photon = SgbranchPhoton.At(i)
		Photon_TLV = photon.P4()                  #TLorentzVector
		threephotons_vec.append(Photon_TLV)
		
	#find two photons 
	ipalp1, ipalp2, imind, egtest = includeme.ass_3l(threephotons_vec, ecm, testmass)
		
	#assing ipalp1 to the photon with max energy between ipalp1 and ipalp2
	if threephotons_vec[ipalp1].E() < threephotons_vec[ipalp2].E():
		ipalp1_temporary = ipalp2
		ipalp2 = ipalp1
		ipalp1 = ipalp1_temporary
		del ipalp1_temporary

	if entry % 100 == 0:
		print testmass, entry, ipalp1, ipalp2, imind, egtest
		
	epho2epho1 = threephotons_vec[ipalp2].E()/threephotons_vec[ipalp1].E()	
	histSgepho2epho1.Fill(epho2epho1)
		
	Sg_DeltaR = threephotons_vec[ipalp1].DeltaR(threephotons_vec[ipalp2])
	histSgdeltar.Fill(Sg_DeltaR)

	MALP = (threephotons_vec[ipalp1]+threephotons_vec[ipalp2]).M()
	epho3 = threephotons_vec[imind].E()
	sigmaepho3 = 0.11*(epho3**0.5)+0.01*epho3
	sigmaalp = 1.057
	ephotest = (ecm*ecm-testmass*testmass)/2./ecm

	MALPcut = ((MALP-testmass)**2/sigmaalp**2+(epho3-ephotest)**2/sigmaepho3**2)**0.5
	
	histSgMALPcut.Fill(MALPcut)

histSgepho2epho1.Write()
histSgdeltar.Write()
histSgMALPcut.Write()
# Loop over background events
for entry in range(0, BgnumberOfEntries):
	# Load selected branches with data from specified event
	BgtreeReader.ReadEntry(entry)

	threephotons_vec = [] #array with 3 TLorentzVector from 3 photons

	#use only events with 3 photons and not tracks
	if BgbranchPhoton.GetEntries() != 3:
		continue
	if BgbranchTrack.GetEntries() != 0:
		continue

	#add lorentzvetor photons in vector
	for i in range(0, BgbranchPhoton.GetEntries()):

		photon = BgbranchPhoton.At(i)
		Photon_TLV = photon.P4()                  #TLorentzVector
		threephotons_vec.append(Photon_TLV)
		
	#find two photons 
	ipalp1, ipalp2, imind, egtest = includeme.ass_3l(threephotons_vec, ecm, testmass)
		

	#assing ipalp1 to the photon with max energy between ipalp1 and ipalp2
	if threephotons_vec[ipalp1].E() < threephotons_vec[ipalp2].E():
		ipalp1_temporary = ipalp2
		ipalp2 = ipalp1
		ipalp1 = ipalp1_temporary
		del ipalp1_temporary
 	
	if entry % 100 == 0:
		print testmass, entry, ipalp1, ipalp2, imind, egtest
		
	epho2epho1 = threephotons_vec[ipalp2].E()/threephotons_vec[ipalp1].E()	
	histBgepho2epho1.Fill(epho2epho1)

	Bg_DeltaR = threephotons_vec[ipalp1].DeltaR(threephotons_vec[ipalp2])
	histBgdeltar.Fill(Bg_DeltaR)

	MALP = (threephotons_vec[ipalp1]+threephotons_vec[ipalp2]).M()
	epho3 = threephotons_vec[imind].E()
	sigmaepho3 = 0.11*(epho3**0.5)+0.01*epho3
	sigmaalp = 1.057
	ephotest = (ecm*ecm-testmass*testmass)/2./ecm

	MALPcut = ((MALP-testmass)**2/sigmaalp**2+(epho3-ephotest)**2/sigmaepho3**2)**0.5
	histBgMALPcut.Fill(MALPcut)


histBgepho2epho1.Write()
histBgdeltar.Write()
histBgMALPcut.Write()

out_root.Close()
