def ass_3l(vp, ecm, testmass, imind = 0, ipalp1 = 0, ipalp2 = 0, egtest = 0):
	'''' take an array of 3lorentz vector for the three photons
	and based on the test mass value in input and on the ecm of accelerator
	assign two of the three photons to ALP decay. 
	return vector of photons assigned to the alp (ipalp1 and ipalp2) and to 
	the third photon (imind), it also returns the ALP mass as calculated from
	the third photon from energy conservation (egtest)'''

	ecm2 = ecm*ecm
	ephotest = (ecm2-testmass*testmass)/2./ecm
	egtest1 = -999999.
	egtest2 = -999999.
	egtest3 = -999999.

	if vp[0].E() < ecm/2:
		egtest1 = (ecm2-2*vp[0].E()*ecm)**0.5
	if vp[1].E() < ecm/2:
		egtest2 = (ecm2-2*vp[1].E()*ecm)**0.5
	if vp[2].E() < ecm/2:
		egtest3 = (ecm2-2*vp[2].E()*ecm)**0.5

	minmd = abs(ephotest-vp[0].E())
	imind = 0
	ipalp1 = 1
	ipalp2 = 2
	egtest = egtest1

	if abs(ephotest-vp[1].E()) < minmd:
		minmd = abs(ephotest - vp[1].E())
		imind = 1
		ipalp1 = 0
		ipalp2 = 2
		egtest = egtest2

	if abs(ephotest-vp[2].E()) < minmd: #fixed there was an error in abs()
		minmd = abs(ephotest-vp[2].E())
		imind = 2
		ipalp1 = 0
		ipalp2 = 1
		egtest = egtest3

	return ipalp1, ipalp2, imind, egtest