#!/usr/bin/env python
#code to show 2d histo of s/sqrt(b)
#to be used with /home/delphesresults/3abg_fcc90_nocut for background samples
#to be used with /home/delphesresults/alpscan for signal samples
#created on 15/5/2019
#from 3a_cutMALPvariables

import sys
import ROOT
import includeme
from array import array
from ROOT import TGraph2D

if len(sys.argv) < 3:
	print " Usage: python examples/analysis.py inputBg.root inputSg.root output.root testmass"
	sys.exit(1)

ROOT.gSystem.Load("libDelphes")

try:
	ROOT.gInterpreter.Declare('#include "classes/SortableObject.h"')
	ROOT.gInterpreter.Declare('#include "classes/DelphesClasses.h"')
	ROOT.gInterpreter.Declare('#include "external/ExRootAnalysis/ExRootTreeReader.h"')
except:
	pass

BgFile = sys.argv[1] #background file
SgFile = sys.argv[2] #signal file
outputFile = sys.argv[3]
test_testmass = sys.argv[4]

# Create chain of root trees
Bgchain = ROOT.TChain("Delphes")
Bgchain.Add(BgFile)

Sgchain = ROOT.TChain("Delphes")
Sgchain.Add(SgFile)

# Create object of class ExRootTreeReader
BgtreeReader = ROOT.ExRootTreeReader(Bgchain)
global BgnumberOfEntries
BgnumberOfEntries = BgtreeReader.GetEntries()
print "bg entries: "+str(BgnumberOfEntries)

SgtreeReader = ROOT.ExRootTreeReader(Sgchain)
global SgnumberOfEntries
SgnumberOfEntries = SgtreeReader.GetEntries()
print "sg entries: "+str(SgnumberOfEntries)

# Get pointers to branches used in this analysis
global BgbranchPhoton
global BgbranchTrack
BgbranchPhoton = BgtreeReader.UseBranch("Photon")
BgbranchTrack = BgtreeReader.UseBranch("Track")

global SgbranchPhoton
global SgbranchTrack
SgbranchPhoton = SgtreeReader.UseBranch("Photon")
SgbranchTrack = SgtreeReader.UseBranch("Track")
SgbrancParticle = SgtreeReader.UseBranch("Particle")

#testmass for ass3l()
global testmass
#testmass = 80
testmass = float(test_testmass)
print testmass

#ecm
global ecm
ecm = 45.8*2 #GeV

#cross-section to rescale graphs
#xsec scanalp90_80 eta=10 = 3.294275e-06 pb
global xsec_scanalp90_1_eta10
global xsec_3abg90_eta10
global luminosity
#couplings theoretical from xsecttocyy.py eta10 ecm=90
#etaa=10.
#[0.00010584473200797015, 0.00011195767438563123, 0.00012333718369766136, 0.00014287565957617757, 0.0001768077913388226, 0.00024087176728341438, 0.0003875681881371373, 0.0008997139223959913, 0.016127450112823835]

#xsec10 = [0.0002380292, 0.0002127458, 0.0001752995, 0.0001306328, 8.530332e-05, 4.596179e-05, 1.775303e-05, 3.294275e-06]
xsec_scanalp90_1_eta10 = 242.5564 #ab #from lorenzo/alp/lorenzo/xsectocyy.py
xsec_3abg90_eta10 = 5.6361e06 #ab #from polesell/work/alpbg/fcc3a/Events/3abf_90_1
luminosity = 150 #ab^-1

out_root = ROOT.TFile(outputFile,"RECREATE")	

#Book histograms
histSgsqrtBEphoDR = ROOT.TH2F("Sg/sqrt(B)", "Sg/sqrt(B)", 100, 0.0, 2.0, 100, 0., 10.)
 
def funcs(cutDR, cutEpho):
	s_counter = 0
	
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

		MALP = (threephotons_vec[ipalp1]+threephotons_vec[ipalp2]).M()
		epho1 = threephotons_vec[ipalp1].E()
		sigmaepho1 =(0.1*0.1*epho1+0.01*0.01*epho1*epho1)**0.5
		epho2 = threephotons_vec[ipalp2].E()
		sigmaepho2 =(0.1*0.1*epho2+0.01*0.01*epho2*epho2)**0.5
		epho3 = threephotons_vec[imind].E()
		sigmaepho3 =(0.1*0.1*epho3+0.01*0.01*epho3*epho3)**0.5
		#sigmaalp = 1.057 #from sigma of reconstructed MALP
		sigmaalp=MALP*0.5*((sigmaepho1/epho1)*(sigmaepho1/epho1)+(sigmaepho2/epho2)*(sigmaepho2/epho2))**0.5

		ephotest = (ecm*ecm-testmass*testmass)/2./ecm

		MALPcut = ((MALP-testmass)**2/(sigmaalp**2)+(epho3-ephotest)**2/(sigmaepho3**2))**0.5
		epho2epho1 = threephotons_vec[ipalp2].E()/threephotons_vec[ipalp1].E()	
		
		Sg_DeltaR = threephotons_vec[ipalp1].DeltaR(threephotons_vec[ipalp2])
			
		if MALPcut < 1.5 and ((Sg_DeltaR < 0.6 and epho2epho1<0.2) or (Sg_DeltaR<0.3 and epho2epho1<0.2)):
			s_counter = s_counter+1
		
	#print s_counter
	return s_counter

def funcb(cutDR, cutEpho):
	b_counter = 0
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
		
		MALP = (threephotons_vec[ipalp1]+threephotons_vec[ipalp2]).M()
		epho1 = threephotons_vec[ipalp1].E()
		sigmaepho1 =(0.1*0.1*epho1+0.01*0.01*epho1*epho1)**0.5
		epho2 = threephotons_vec[ipalp2].E()
		sigmaepho2 =(0.1*0.1*epho2+0.01*0.01*epho2*epho2)**0.5
		epho3 = threephotons_vec[imind].E()
		sigmaepho3 =(0.1*0.1*epho3+0.01*0.01*epho3*epho3)**0.5
		#sigmaalp = 1.057 #from sigma of reconstructed MALP
		sigmaalp=MALP*0.5*((sigmaepho1/epho1)*(sigmaepho1/epho1)+(sigmaepho2/epho2)*(sigmaepho2/epho2))**0.5
		
		ephotest = (ecm*ecm-testmass*testmass)/2./ecm

		MALPcut = ((MALP-testmass)**2/(sigmaalp**2)+(epho3-ephotest)**2/(sigmaepho3**2))**0.5
		
		epho2epho1 = threephotons_vec[ipalp2].E()/threephotons_vec[ipalp1].E()	
	
		Bg_DeltaR = threephotons_vec[ipalp1].DeltaR(threephotons_vec[ipalp2])
	
		if MALPcut < 1.5 and ((Bg_DeltaR < 0.6 and epho2epho1<0.2) or (Bg_DeltaR<0.3 and epho2epho1<0.2)):
			b_counter = b_counter+1
		
	#print b_counter
	return b_counter

def funcssqrtb(cutDR, cutEpho):
	s = funcs(cutDR, cutEpho)
	b = funcb(cutDR, cutEpho)
	s1 = s*xsec_scanalp90_1_eta10*luminosity/SgnumberOfEntries
	b1 = b*xsec_3abg90_eta10*luminosity/BgnumberOfEntries
	print str(cutDR)+"  "+str(cutEpho)+"  "+str(s)+"  "+str(b)+"  "+str(s1)+"  "+str(b1)+"  "+str(s1/(b1**0.5))+"\n"
	return s,b,s1/(b1**0.5)  

def funcshisto():
	#Book histograms
	histSgepho2epho1 = ROOT.TH1F("Sg_Epho2/Epho1", "Sg_Epho2/Epho1", 150, 0.0, 1.1)
	histSgdeltar = ROOT.TH1F("Sg_DeltaR", "Sg_DeltaR", 20, 0., 0.2)
	histSgMALP = ROOT.TH1F("Sg_MALP", "Sg", 100, 0., 4.)
	histSgMALPcut = ROOT.TH1F("Sg_MALPcut", "Sg", 100, 0., 50.)
	histSgEphoDR = ROOT.TH2F("Sg", "Sg", 100, 0.0, 1.1, 100, 0., 6.)
	histSgEphoMalp = ROOT.TH2F("Sg_2","Sg_2", 100, -4, 4, 100, -5, 5)
	histSgEphoDR_afterMALPcut = ROOT.TH2F("Sg_afterMALPcut", "Sg_afterMALPcut", 100, 0.0, 1.1, 100, 0., 6.)
	histSgEpho3 = ROOT.TH1F("Sg_Epho3", "Sg_epho3", 100, 0., 100)
	histepho3egtest = ROOT.TH1F("Sg_epho3egtest", "Sg_epho3egtest", 100, -6, 6.)
	histTruthSgdeltar = ROOT.TH1F("Truth_Sg_DeltaR", "Truth_Sg_DeltaR", 100, 0., 2.)
	
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

		truthphoton1 = SgbrancParticle.At(7)
		truthphoton2 = SgbrancParticle.At(8)
		if truthphoton1.M1 != 5 or truthphoton2.M1 != 5:
			continue
		
		truthP1 = truthphoton1.P4()
		truthP2 = truthphoton2.P4()
		truthDR = truthP1.DeltaR(truthP2)
		histTruthSgdeltar.Fill(truthDR)
		
		#find two photons 
		ipalp1, ipalp2, imind, egtest = includeme.ass_3l(threephotons_vec, ecm, testmass)
		
		#assing ipalp1 to the photon with max energy between ipalp1 and ipalp2
		if threephotons_vec[ipalp1].E() < threephotons_vec[ipalp2].E():
			ipalp1_temporary = ipalp2
			ipalp2 = ipalp1
			ipalp1 = ipalp1_temporary
			del ipalp1_temporary

		MALP = (threephotons_vec[ipalp1]+threephotons_vec[ipalp2]).M()
		#print MALP
		epho1 = threephotons_vec[ipalp1].E()
		sigmaepho1 =(0.11*0.11*epho1+0.01*0.01*epho1*epho1)**0.5
		epho2 = threephotons_vec[ipalp2].E()
		sigmaepho2 =(0.11*0.11*epho2+0.01*0.01*epho2*epho2)**0.5
		epho3 = threephotons_vec[imind].E()
		sigmaepho3 =(0.11*0.11*epho3+0.01*0.01*epho3*epho3)**0.5
		
		#sigmaalp=MALP*0.5*((sigmaepho1/epho1)*(sigmaepho1/epho1)+(sigmaepho2/epho2)*(sigmaepho2/epho2))**0.5
		#sgimaalp modified to take into account error on DR
		sigmaalp=((MALP*0.5*((sigmaepho1/epho1)*(sigmaepho1/epho1)+(sigmaepho2/epho2)*(sigmaepho2/epho2))**0.5)**2+(0.15*MALP)**2)**0.5
		ephotest = (ecm*ecm-testmass*testmass)/2./ecm

		epho2epho1 = threephotons_vec[ipalp2].E()/threephotons_vec[ipalp1].E()	
		Sg_DeltaR = threephotons_vec[ipalp1].DeltaR(threephotons_vec[ipalp2])
		
		MALPcut = ((MALP-testmass)**2/(sigmaalp**2)+(epho3-ephotest)**2/(sigmaepho3**2))**0.5
		
		histepho3egtest.Fill((MALP-testmass)/sigmaalp)
		histSgEpho3.Fill(epho3)
		histSgEphoMalp.Fill((epho3-ephotest)/sigmaepho3, (MALP-testmass)/sigmaalp )
		histSgepho2epho1.Fill(epho2epho1)
		histSgdeltar.Fill(Sg_DeltaR)
		histSgMALP.Fill(MALP)
		histSgMALPcut.Fill(MALPcut)
		histSgEphoDR.Fill(epho2epho1, Sg_DeltaR)
		if MALPcut<1.5:
			histSgEphoDR_afterMALPcut.Fill(epho2epho1, Sg_DeltaR)

	histTruthSgdeltar.Write()
	histepho3egtest.Write()
	histSgEpho3.Write()
	histSgEphoDR_afterMALPcut.Write()
	histSgEphoMalp.Write()
	histSgepho2epho1.Write()
	histSgdeltar.Write()
	histSgMALP.Write()
	histSgMALPcut.Write()
	histSgEphoDR.Write()

def funcbhisto():
	histBgepho2epho1 = ROOT.TH1F("Bg_Epho2/Epho1", "Bg_Epho2/Epho1", 150, 0.0, 1.1)
	histBgdeltar = ROOT.TH1F("Bg_DeltaR", "Bg_DeltaR", 100, 0., 10.)
	histBgMALP = ROOT.TH1F("Bg_MALP", "Bg", 100, 0., 100.)
	histBgMALPcut = ROOT.TH1F("Bg_MALPcut", "Bg", 100, 0., 50.)
	histBgEphoDR = ROOT.TH2F("Bg", "Bg", 100, 0.00, 1.1, 100, 0., 6.)
	histBgEphoMalp = ROOT.TH2F("Bg_2","Bg_2", 100, -70, 35, 100, -90, 60)
	histBgEphoDR_afterMALPcut = ROOT.TH2F("Bg_afterMALPcut", "Bg_afterMALPcut", 100, 0.0, 1.1, 100, 0.0, 6.0)

	# Loop over signal events
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

		MALP = (threephotons_vec[ipalp1]+threephotons_vec[ipalp2]).M()
		#print MALP
		epho1 = threephotons_vec[ipalp1].E()
		sigmaepho1 =(0.1*0.1*epho1+0.01*0.01*epho1*epho1)**0.5
		epho2 = threephotons_vec[ipalp2].E()
		sigmaepho2 =(0.1*0.1*epho2+0.01*0.01*epho2*epho2)**0.5
		epho3 = threephotons_vec[imind].E()
		sigmaepho3 =(0.1*0.1*epho3+0.01*0.01*epho3*epho3)**0.5
		#sigmaalp = 1.057 #from sigma of reconstructed MALP
		sigmaalp=MALP*0.5*((sigmaepho1/epho1)*(sigmaepho1/epho1)+(sigmaepho2/epho2)*(sigmaepho2/epho2))**0.5
		
		ephotest = (ecm*ecm-testmass*testmass)/2./ecm

		epho2epho1 = threephotons_vec[ipalp2].E()/threephotons_vec[ipalp1].E()	
		Bg_DeltaR = threephotons_vec[ipalp1].DeltaR(threephotons_vec[ipalp2])
		
		MALPcut = ((MALP-testmass)**2/(sigmaalp**2)+(epho3-ephotest)**2/(sigmaepho3)**2)**0.5
		#print (epho3-ephotest)/sigmaepho3, (MALP-testmass)/sigmaalp
		histBgEphoMalp.Fill((epho3-ephotest)/sigmaepho3, (MALP-testmass)/sigmaalp )
		histBgepho2epho1.Fill(epho2epho1)
		histBgdeltar.Fill(Bg_DeltaR)
		histBgMALP.Fill(MALP)
		histBgMALPcut.Fill(MALPcut)
		histBgEphoDR.Fill(epho2epho1, Bg_DeltaR)
		if MALPcut<1.5:
			histBgEphoDR_afterMALPcut.Fill(epho2epho1, Bg_DeltaR)

	histBgEphoDR_afterMALPcut.Write()
	histBgEphoMalp.Write()
	histBgepho2epho1.Write()
	histBgdeltar.Write()
	histBgMALP.Write()
	histBgMALPcut.Write()
	histBgEphoDR.Write()

cutDRarray = array('d',[])
cutEphoarray=array('d',[])
sbarray = array('d',[])
sarray = array('d',[])
barray = array('d',[])

#funcssqrtb(100,100)
'''
for j in range(5):
	cutEpho = 0.8+0.05*j
	for i in range(10):
		cutDR = 3.0+(0.05*i)
		s,b,sb = funcssqrtb(cutDR, cutEpho)
		cutDRarray.append(cutDR)
		cutEphoarray.append(cutEpho)
		sbarray.append(sb)
		sarray.append(s)
		barray.append(b)

n = len(sarray)
SGraph = TGraph2D(n, cutDRarray, cutEphoarray, sarray)
SGraph.SetTitle("Signal")
BGraph = TGraph2D(n, cutDRarray, cutEphoarray, barray)
BGraph.SetTitle("Bg")
SBGraph = TGraph2D(n, cutDRarray, cutEphoarray, sbarray)
SBGraph.SetTitle("S/sqrt(B)")

SBGraph.Write()
BGraph.Write() 
SGraph.Write()
'''
funcshisto()
funcbhisto()

out_root.Close()