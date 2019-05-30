#!/usr/bin/env python
#code to show histo of background and signal samples in order to understand how to apply cuts
#to be used with /home/delphesresults/3abg_fcc90_nocut for background samples
#to be used with /home/delphesresults/alpscan for signal samples
#created on 10/5/2019
#from 3a_cutvariables, here I implement the cut on the Mass ad first cut

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

# Get pointers to branches used in this analysis
BgbranchPhoton = BgtreeReader.UseBranch("Photon")
BgbranchTrack = BgtreeReader.UseBranch("Track")

SgbranchPhoton = SgtreeReader.UseBranch("Photon")
SgbranchTrack = SgtreeReader.UseBranch("Track")

testmass = 80

#ecm
ecm = 45.8*2 #GeV

#cross-section to rescale graphs
#xsec scanalp90_80 eta=10 = 3.294275e-06 pb
xsec_scanalp90_80_eta10 = 3.294275 #ab #from lorenzo/alp/lorenzo/xsectocyy.py
xsec_3abg90_eta10 = 0.070153e06 #ab from polesell/work/alpbg/fcc3a/Events/3abf_90_1
luminosity = 150 #ab^-1

#counter for statistics afer cuts
Sg_3a_counter = 0
Sg_3a_MALP2_counter = 0

Bg_3a_counter = 0
Bg_3a_MALP2_counter = 0

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

histBgEphoDR = ROOT.TH2F("Bg", "Bg", 100, 0.65, 1.1, 100, 2., 6.)
histSgEphoDR = ROOT.TH2F("Sg", "Sg", 100, 0.65, 1.1, 100, 2., 6.)

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
	Sg_3a_counter = Sg_3a_counter+1

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
	
	MALP = (threephotons_vec[ipalp1]+threephotons_vec[ipalp2]).M()
	epho3 = threephotons_vec[imind].E()
	sigmaepho3 = 0.11*(epho3**0.5)+0.01*epho3
	sigmaalp = 1.057 #from sigma of reconstructed MALP
	ephotest = (ecm*ecm-testmass*testmass)/2./ecm

	MALPcut = ((MALP-testmass)**2/sigmaalp**2+(epho3-ephotest)**2/sigmaepho3**2)**0.5

	if MALPcut > 2.:
		continue
	Sg_3a_MALP2_counter = Sg_3a_MALP2_counter+1

	epho2epho1 = threephotons_vec[ipalp2].E()/threephotons_vec[ipalp1].E()	
	histSgepho2epho1.Fill(epho2epho1, xsec_scanalp90_80_eta10*1000)
		
	Sg_DeltaR = threephotons_vec[ipalp1].DeltaR(threephotons_vec[ipalp2])
	histSgdeltar.Fill(Sg_DeltaR)

	histSgMALPcut.Fill(MALPcut)

	if Sg_DeltaR > 2.5:
		histSgEphoDR.Fill(epho2epho1, Sg_DeltaR)

histSgepho2epho1.Write()
histSgdeltar.Write()
histSgMALPcut.Write()
histSgEphoDR.Write()

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
	Bg_3a_counter = Bg_3a_counter+1

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
		
	MALP = (threephotons_vec[ipalp1]+threephotons_vec[ipalp2]).M()
	epho3 = threephotons_vec[imind].E()
	sigmaepho3 = 0.11*(epho3**0.5)+0.01*epho3
	sigmaalp = 1.057
	ephotest = (ecm*ecm-testmass*testmass)/2./ecm

	MALPcut = ((MALP-testmass)**2/sigmaalp**2+(epho3-ephotest)**2/sigmaepho3**2)**0.5
	
	if MALPcut > 2.:
		continue
	Bg_3a_MALP2_counter = Bg_3a_MALP2_counter + 1

	epho2epho1 = threephotons_vec[ipalp2].E()/threephotons_vec[ipalp1].E()	
	histBgepho2epho1.Fill(epho2epho1, xsec_3abg90_eta10)

	Bg_DeltaR = threephotons_vec[ipalp1].DeltaR(threephotons_vec[ipalp2])
	histBgdeltar.Fill(Bg_DeltaR)

	histBgMALPcut.Fill(MALPcut)

	if Bg_DeltaR > 2.5:
		histBgEphoDR.Fill(epho2epho1, Bg_DeltaR)

histBgepho2epho1.Write()
histBgdeltar.Write()
histBgMALPcut.Write()
histBgEphoDR.Write()

out_root.Close()

print str(BgnumberOfEntries)+" background events\n"
print str(SgnumberOfEntries)+" signal events\n"
print str(Bg_3a_counter)+" background events 3a\n"
print str(Sg_3a_counter)+" signal events 3a\n"
print str(Sg_3a_MALP2_counter)+" signal event 3a MALP<2\n"
print str(Bg_3a_MALP2_counter)+" background event 3a MALP<2\n"
