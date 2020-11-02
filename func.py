import ROOT
import math

def plotSignificance(s1, b1):
	#significance calculated from https://cds.cern.ch/record/2643488
	# relabel variables to match CDS formula
	observed = s1+b1
	error = 0.05
	nbObs = observed
	nbExp = b1
	nbExpEr = error*nbExp

	#print 'calculating significance from W. Buttinger and M.Lefebvre recommendation'
	factor1 = nbObs*math.log( (nbObs*(nbExp+nbExpEr**2))/(nbExp**2+nbObs*nbExpEr**2) )
	factor2 = (nbExp**2/nbExpEr**2)*math.log( 1 + (nbExpEr**2*(nbObs-nbExp))/(nbExp*(nbExp+nbExpEr**2)) )

	if nbObs < nbExp:
		pull  = -math.sqrt(2*(factor1 - factor2))
	else:
		pull  = math.sqrt(2*(factor1 - factor2))

	#print "pull: "+str(pull)
	return pull


plot = ROOT.TH2F("pull", "pull", 90, 1., 10., 9999, 1., 1000.)

signal = [0.1+0.1*x for x in range(100)]
backgroung = [0.1+0.1*x for x in range(10000)]

for s in signal:
	for b in backgroung:
		pull = plotSignificance(s,b)
		plot.Fill(s,b,pull)

func = ROOT.TF1("f",'(x/2)^2',1.,10.)

ROOT.gStyle.SetPalette(55)
outputfile = ROOT.TFile("output.root","RECREATE")
func.Write()
plot.GetXaxis().SetTitle("Signal")
plot.GetYaxis().SetTitle("Background")
plot.Write()
outputfile.Close()

