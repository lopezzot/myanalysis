from ROOT import TGraph, TFile, TCanvas, TMultiGraph
from array import array

coupligs = array('d')
masses = array('d')
theoretical = array('d')
theoretical365 = array('d')
theoretical240 = array('d')
theoretical160 = array('d')

coupling365 = array('d')
masses365 = array('d')

coupling240 = array('d')
masses240 = array('d')

coupling160 = array('d')
masses160 = array('d')

rootfile = TFile("output.root","RECREATE")
c = [10000, 0.007146, 0.00446, 0.00540, 0.00686, 0.00765, 0.00885, 0.0143, 0.018, 0.046, 10000, 10000]
m = [5, 5, 10, 20, 30, 40, 50, 60, 70, 80, 80, 5]

t = [10000,0.00010531464629643753, 0.00010584473200797015, 0.00011195767438563123, 0.00012333718369766136, 0.00014287565957617757, 0.0001768077913388226, 0.00024087176728341438, 0.0003875681881371373, 0.0008997139223959913, 10000, 10000]

c365 = [10000, 0.078215, 0.09351, 0.122627, 0.16467, 0.25945, 0.4986, 23.347, 10000]
m365 = [50, 50, 100, 150, 200, 240, 300, 360, 360]

m240 = [100, 100, 150, 180, 220, 235, 235]
c240 = [10000, 0.122954, 0.156492, 0.271, 1.5352, 13.147, 10000]
t240 = [10000, 0.008302866109567658, 0.014043267960524476, 0.02400387934037744, 0.1135393636633489, 0.8767545798463641, 10000]

t365 = [10000, 0.011805359356280826, 0.01291984002946909, 0.016237318278912834, 0.022378778938076393, 0.03174363246917037, 0.07593830445852476, 3.1646458992183386, 10000]

m160 = [50, 50, 100, 120, 140, 155, 155]
c160 = [10000, 0.117973, 0.157546, 0.153418175, 0.765888, 8.2, 10000]

t160 = [10000, 0.005868095964649793, 0.00897964053112733, 0.015062350679744235, 0.03967262339072189, 0.30158366478681703, 10000]

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
for p in t365:
	theoretical365.append(p)

for t in c240:
	coupling240.append(t)
for d in m240:
	masses240.append(d)
for i in t240:
	theoretical240.append(i)

for t in c160:
	coupling160.append(t)
for d in m160:
	masses160.append(d)
for i in t160:
	theoretical160.append(i)


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
Exclusiongraph.SetLineWidth(1)

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
Exclusiongraph2.SetLineWidth(1)
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
Exclusiongraph365.SetLineWidth(1)

Exclusiongraph160 = TGraph(len(m160), masses160, coupling160)
Exclusiongraph160.SetTitle("")
Exclusiongraph160.GetHistogram().SetMaximum(2000.)
Exclusiongraph160.GetHistogram().SetMinimum(0.1e-6)
XAxis = Exclusiongraph160.GetXaxis()
XAxis.SetLimits(1.e-15,2000)
XAxis.SetNdivisions(105)
XAxis.SetTitle("M_{a} [GeV]")
XAxis.SetTitleOffset(1.5)
YAxis = Exclusiongraph160.GetYaxis()
YAxis.SetNdivisions(105)
YAxis.SetTitle("c")
Exclusiongraph160.SetFillColorAlpha(2, 0.5)
Exclusiongraph160.SetLineColor(2)
Exclusiongraph160.SetLineWidth(1)

Exclusiongraph240 = TGraph(len(m240), masses240, coupling240)
Exclusiongraph240.SetTitle("")
Exclusiongraph240.GetHistogram().SetMaximum(2000.)
Exclusiongraph240.GetHistogram().SetMinimum(0.1e-6)
XAxis = Exclusiongraph240.GetXaxis()
XAxis.SetLimits(1.e-15,2000)
XAxis.SetNdivisions(105)
XAxis.SetTitle("M_{a} [GeV]")
XAxis.SetTitleOffset(1.5)
YAxis = Exclusiongraph240.GetYaxis()
YAxis.SetNdivisions(105)
YAxis.SetTitle("c")
Exclusiongraph240.SetFillColorAlpha(5, 0.5)
Exclusiongraph240.SetLineColor(5)
Exclusiongraph240.SetLineWidth(1)

Exclusiongrapht365 = TGraph(len(masses365), masses365, theoretical365)
Exclusiongrapht365.SetTitle("")
Exclusiongrapht365.GetHistogram().SetMaximum(2000.)
Exclusiongrapht365.GetHistogram().SetMinimum(0.1e-6)
XAxis = Exclusiongrapht365.GetXaxis()
XAxis.SetLimits(1.e-15,2000)
XAxis.SetNdivisions(105)
XAxis.SetTitle("M_{a} [GeV]")
XAxis.SetTitleOffset(1.5)
YAxis = Exclusiongrapht365.GetYaxis()
YAxis.SetNdivisions(105)
YAxis.SetTitle("c")
Exclusiongrapht365.SetFillColorAlpha(3, 0.5)
Exclusiongrapht365.SetLineColor(3)
Exclusiongrapht365.SetLineWidth(1)
Exclusiongrapht365.SetLineStyle(3)
Exclusiongrapht365.SetFillStyle(3001)

Exclusiongrapht240 = TGraph(len(masses240), masses240, theoretical240)
Exclusiongrapht240.SetTitle("")
Exclusiongrapht240.GetHistogram().SetMaximum(2000.)
Exclusiongrapht240.GetHistogram().SetMinimum(0.1e-6)
XAxis = Exclusiongrapht240.GetXaxis()
XAxis.SetLimits(1.e-15,2000)
XAxis.SetNdivisions(105)
XAxis.SetTitle("M_{a} [GeV]")
XAxis.SetTitleOffset(1.5)
YAxis = Exclusiongrapht240.GetYaxis()
YAxis.SetNdivisions(105)
YAxis.SetTitle("c")
Exclusiongrapht240.SetFillColorAlpha(5, 0.5)
Exclusiongrapht240.SetLineColor(5)
Exclusiongrapht240.SetLineWidth(1)
Exclusiongrapht240.SetLineStyle(5)
Exclusiongrapht240.SetFillStyle(3001)

Exclusiongrapht160 = TGraph(len(masses160), masses160, theoretical160)
Exclusiongrapht160.SetTitle("")
Exclusiongrapht160.GetHistogram().SetMaximum(2000.)
Exclusiongrapht160.GetHistogram().SetMinimum(0.1e-6)
XAxis = Exclusiongrapht160.GetXaxis()
XAxis.SetLimits(1.e-15,2000)
XAxis.SetNdivisions(105)
XAxis.SetTitle("M_{a} [GeV]")
XAxis.SetTitleOffset(1.5)
YAxis = Exclusiongrapht160.GetYaxis()
YAxis.SetNdivisions(105)
YAxis.SetTitle("c")
Exclusiongrapht160.SetFillColorAlpha(2, 0.5)
Exclusiongrapht160.SetLineColor(2)
Exclusiongrapht160.SetLineWidth(1)
Exclusiongrapht160.SetLineStyle(5)
Exclusiongrapht160.SetFillStyle(3001)

'''
Exclusiongraph365.Draw("afl")
Exclusiongraph2.Draw("same afl")

Exclusiongraph.Draw("same fl")
'''

mg = TMultiGraph()
mg.Add(Exclusiongraph365)
mg.Add(Exclusiongrapht365)
mg.Add(Exclusiongraph240)
mg.Add(Exclusiongrapht240)
mg.Add(Exclusiongrapht160)
mg.Add(Exclusiongraph160)
mg.Add(Exclusiongraph)
mg.Add(Exclusiongraph2)
mg.SetTitle("")
mg.GetHistogram().SetMaximum(2000.)
mg.GetHistogram().SetMinimum(0.1e-6)
XAxis = mg.GetXaxis()
#XAxis.SetLimits(1.e-15,2000)
XAxis.SetLimits(3,1000)
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
