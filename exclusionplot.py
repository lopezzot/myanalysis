from ROOT import TGraph, TFile, TCanvas
from array import array

coupligs = array('d')
masses = array('d')
theoretical = array('d')
rootfile = TFile("output.root","RECREATE")
c = [10000, 0.00446, 0.00540, 0.00686, 0.00765, 0.00885, 0.0143, 0.018, 0.046, 10000, 10000]
m = [10, 10, 20, 30, 40, 50, 60, 70, 80, 80, 10]

t = [10000, 0.00010584473200797015, 0.00011195767438563123, 0.00012333718369766136, 0.00014287565957617757, 0.0001768077913388226, 0.00024087176728341438, 0.0003875681881371373, 0.0008997139223959913, 10000, 10000]

for i in c:
	coupligs.append(i)
for x in m:
	masses.append(x)
for y in t:
	theoretical.append(y)

canvas=TCanvas("canvas","canvas",800,800)
canvas.SetLogy()
canvas.SetLogx()
Exclusiongraph = TGraph(len(m), masses, coupligs)
Exclusiongraph.GetHistogram().SetMaximum(2000.)
Exclusiongraph.GetHistogram().SetMinimum(0.1e-6)
XAxis = Exclusiongraph.GetXaxis()
XAxis.SetLimits(1.e-15,2000)


Exclusiongraph2 = TGraph(len(m), masses, theoretical)
Exclusiongraph2.SetTitle("")
Exclusiongraph2.GetHistogram().SetMaximum(2000.)
Exclusiongraph2.GetHistogram().SetMinimum(0.1e-6)
XAxis = Exclusiongraph2.GetXaxis()
XAxis.SetLimits(1.e-15,2000)
XAxis.SetNdivisions(105)
XAxis.SetTitle("M_{a} [GeV]")
XAxis.SetTitleOffset(1.5)
YAxis = Exclusiongraph2.GetYaxis()
YAxis.SetNdivisions(105)
YAxis.SetTitle("c")
Exclusiongraph.SetFillColorAlpha(7, 0.5)
Exclusiongraph.SetLineColor(7)
Exclusiongraph.SetLineWidth(3)
Exclusiongraph2.SetFillColorAlpha(6, 0.5)
Exclusiongraph2.SetLineColor(6)
Exclusiongraph2.SetLineWidth(3)
Exclusiongraph2.Draw("afl")
Exclusiongraph.Draw("same fl")

Exclusiongraph.Write()
Exclusiongraph2.Write()
canvas.SaveAs("canvas.pdf")