from ROOT import TGraph, TFile, TCanvas, TMultiGraph
from array import array

coupligs = array('d')
masses = array('d')
theoretical = array('d')

coupling365 = array('d')
masses365 = array('d')

rootfile = TFile("output.root","RECREATE")
c = [10000, 0.00446, 0.00540, 0.00686, 0.00765, 0.00885, 0.0143, 0.018, 0.046, 10000, 10000]
m = [10, 10, 20, 30, 40, 50, 60, 70, 80, 80, 10]

t = [10000, 0.00010584473200797015, 0.00011195767438563123, 0.00012333718369766136, 0.00014287565957617757, 0.0001768077913388226, 0.00024087176728341438, 0.0003875681881371373, 0.0008997139223959913, 10000, 10000]

c365 = [10000, 0.078215, 0.09351, 0.122627, 0.16467, 0.25945, 0.4986, 23.347, 10000]
m365 = [50, 50, 100, 150, 200, 240, 300, 360, 360]

for i in c:
	coupligs.append(i)
for x in m:
	masses.append(x)
for y in t:
	theoretical.append(y)

for r in c365:
	coupling365.append(r)
for t in m365:
	masses365.append(t)

canvas=TCanvas("canvas","canvas",800,800)
canvas.SetLogy()
canvas.SetLogx()

Exclusiongraph = TGraph(len(m), masses, coupligs)
Exclusiongraph.GetHistogram().SetMaximum(2000.)
Exclusiongraph.GetHistogram().SetMinimum(0.1e-6)
XAxis = Exclusiongraph.GetXaxis()
XAxis.SetLimits(1.e-15,2000)
Exclusiongraph.SetFillColorAlpha(7, 0.5)
Exclusiongraph.SetLineColor(7)
Exclusiongraph.SetLineWidth(3)

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
Exclusiongraph2.SetFillColorAlpha(7, 0.5)
Exclusiongraph2.SetLineColor(7)
Exclusiongraph2.SetLineWidth(3)
Exclusiongraph2.SetLineStyle(3)
Exclusiongraph2.SetFillStyle(3001)

Exclusiongraph365 = TGraph(len(m365), masses365, coupling365)
Exclusiongraph365.SetTitle("")
Exclusiongraph365.GetHistogram().SetMaximum(2000.)
Exclusiongraph365.GetHistogram().SetMinimum(0.1e-6)
XAxis = Exclusiongraph365.GetXaxis()
XAxis.SetLimits(1.e-15,2000)
XAxis.SetNdivisions(105)
XAxis.SetTitle("M_{a} [GeV]")
XAxis.SetTitleOffset(1.5)
YAxis = Exclusiongraph365.GetYaxis()
YAxis.SetNdivisions(105)
YAxis.SetTitle("c")
Exclusiongraph365.SetFillColorAlpha(3, 0.5)
Exclusiongraph365.SetLineColor(3)
Exclusiongraph365.SetLineWidth(3)
'''
Exclusiongraph365.Draw("afl")
Exclusiongraph2.Draw("same afl")

Exclusiongraph.Draw("same fl")
'''

mg = TMultiGraph()
mg.Add(Exclusiongraph)
mg.Add(Exclusiongraph2)
mg.Add(Exclusiongraph365)
mg.SetTitle("")
mg.GetHistogram().SetMaximum(2000.)
mg.GetHistogram().SetMinimum(0.1e-6)
XAxis = mg.GetXaxis()
XAxis.SetLimits(1.e-15,2000)
XAxis.SetNdivisions(105)
XAxis.SetTitle("M_{a} [GeV]")
XAxis.SetTitleOffset(1.5)
YAxis = mg.GetYaxis()
YAxis.SetNdivisions(105)
YAxis.SetTitle("c")
mg.Draw("afl")
canvas.SaveAs("canvas.pdf")

Exclusiongraph.Write()
Exclusiongraph2.Write()
Exclusiongraph365.Write()