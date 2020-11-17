def sabetta(Mw,Re,Re1,isito,isig,dt,nacc,tot_dur,scale):
	import numpy as np
	nacc = int(nacc)
	'''======================================================================= 
	function for the generation of non-stationary accelerograms 
	adapted from original fortran program of Sabetta F. and Pugliese A. (1995)
	======================================================================= 
	input parameters
	Mw = moment magnitude
	Re = epicentral distance (km)
	isito = site conditions (0=rock; 1=shallow all.; 2=deep alluvium)
	isig = st. dev. of GMPE (0=median value,1=84th percentile)
	dt = time step of output accelerogram
	nacc = number of accelerograms to be generated
	tot_dur = total duration of accelerogram
	scale = scale factor (1 for cm/s/s)
	======================================================================= 
	output
	t = time vector (s)
	acc = acceleration time history (in cm/s/s)
	vel = velocity time history (in cm/s)
	disp = displacement time history (in cm)
	======================================================================='''

	if isito == 0:
		S1=0
		S2=0
	elif isito == 1:
		S1=1
		S2=0
	elif isito == 2:
		S1=0
		S2=1

	''' empirical relationship for strong ground motion duration (DV) and 
	Arias intensity (Ia) developed by Sabetta & Pugliese (1996) 
	Definition of duration according to the definition given by Vanmarcke & Lai (1980) '''
	DV = 10**(-0.783 + 0.193*Mw + 0.208*np.log10(np.sqrt(Re1**2 + 5.1**2)) - 0.133*S1 + 0.138*S2 + 0.247*isig)
	Ia = 10**(0.729 + 0.911*Mw - 1.818*np.log10(np.sqrt(Re1**2 + 5.3**2)) + 0.244*S1 + 0.139*S2 + 0.397*isig)
	# time delay in seconds between S ans P waves (VP*VS/(VP-VS) = 7 km/s)
	T1 = Re/7

	''' values of times and sigma for the definition of Pa = lognormal time 
	envelope desribing the amplitude variation of ground motion 
	(see S&P96, p. 346) '''
	T2 = T1 + 0.5*DV
	T3 = T1 + 2.5*DV
	TFc = T2 - 3.5 - Re/50
	T_cost = T2 + 0.5*DV
	T4 = tot_dur - T1
	T_fond = T4/3
	fo = 1/T_fond

	t = np.arange(0,tot_dur,dt)
	nt = len(t)
	# Nyquist frequency
	fmax = 1/(2*dt)

	# lognormal time envelope function desribing the amplitude variation of
	# ground motion 
	# Pa [cm^2/s^4]
	# stand. dev. (NB: /2.5 nell'articolo ma /3 nel .for!) 
	sqm_Pa = np.log(T3/T2)/3 
	# mean value 
	med_Pa = np.log(T2) + sqm_Pa**2 
	Pa = np.zeros(nt) 
	Pa = Ia* (np.exp(-(np.log(t) - med_Pa)**2 / (2 * sqm_Pa**2)) / (t * sqm_Pa * np.sqrt(2 * np.pi)))
	t_val = t - TFc
	for i in range(len(t)):
		if t_val[i] < 1:
			t_val[i] = 1
		if t[i] > T_cost:
			t_val[i] = t_val[i-1]

	# empirical regression for Fc, i.e., central frequency 
	# Fc [Hz]
	Fc = np.exp(3.4 - 0.35*np.log(t_val) - 0.218*Mw - 0.15*S2)

	# empirical regression for the ratio Fb/Fc between the frequency bandwidth 
	# and the central frequency
	Fb_Fc = 0.44 + 0.07*Mw - 0.08*S1 + 0.03*S2

	delta = np.sqrt(np.log(1+Fb_Fc**2))
	ln_beta = np.log(Fc) - 0.5*delta**2
	   
	f = np.arange(fo,fmax,fo)
	nf = len(f)
	ind_f = np.arange(1,nf+1) 

	R = np.zeros(nf) 
	PS = np.zeros(nf)
	Ccos = np.zeros(nf)
	Ccos_vel = np.zeros(nf) 
	Ccos_dis = np.zeros(nf) 
	acc = np.zeros((nt,nacc))
	vel = np.zeros((nt,nacc))
	dis = np.zeros((nt,nacc))

	for k in range(nacc):
		R = np.random.uniform(0,2*np.pi,nf)
		for i in range(nt):
			# PS in cm^2 / s^3
			PS = (Pa[i]/(ind_f*np.sqrt(2*np.pi)*delta))*np.exp(-(np.log(f) - ln_beta[i])**2/(2*delta**2))
			# Ccos in cm / s^2
			Ccos = np.sqrt(2*PS)*np.cos(2.*np.pi*f*t[i] + R)
			# Ccos_vel in cm / s      
			Ccos_vel = Ccos/(2*np.pi*f)
			# Ccos_dis in cm      
			Ccos_dis = Ccos/(2*np.pi*f)**2
			# acc in cm/s/s
			acc[i,k] = np.sum(Ccos)
			# vel in cm/s
			vel[i,k] = np.sum(Ccos_vel)
			# dis in cm
			dis[i,k] = np.sum(Ccos_dis)
	acc = acc*scale 
	vel = vel*scale 
	dis = dis*scale
	return t, acc, vel, dis
