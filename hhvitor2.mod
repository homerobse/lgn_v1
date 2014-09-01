TITLE HH currents Pospischil et al

UNITS {
        (mA) = (milliamp)
        (mV) = (millivolt)
		
}
 
NEURON {
	SUFFIX hhvitor2
	NONSPECIFIC_CURRENT ileak, ina, ikd, im, il, it
	RANGE gleak, gna, gkd, gm, Vt, gl, ekd, gt
	
}

PARAMETER {
    v						(mV)
    gleak = 0.00000000273	(mho/cm2)
	eleak = -70				(mV)
	
	gna = 0.039  			(mho/cm2)
	ena = 50				(mV)
	
	gkd = 0.006				(mho/cm2)
	ekd = -90				(mV)

	gm = 0.00013			(mho/cm2)
	ptaumax = 1123			(ms)
	
	gl = 0.0001				(mho/cm2)
	eca = 120				(mV)
	
	Vt = -60
	
	gt = 0.0008				(mho/cm2)
	
	}
	
STATE {
	m
	h
	n
	p
	q
	r
	s
	u
}

ASSIGNED {
	i 		(mA/cm2)
	ileak   (mA/cm2)
	ina		(mA/cm2)
	ikd		(mA/cm2)
	im		(mA/cm2)
	il		(mA/cm2)
	it		(mA/cm2)
	minf
	mtau    (ms)
	hinf
	htau    (ms)
	ninf
	ntau	(ms)
	pinf
	ptau	(ms)
	rinf
	rtau	(ms)
	qinf
	qtau	(ms)
	sinf    
	stau	(ms)
	uinf	
	utau	(ms)
}

BREAKPOINT {
        SOLVE states METHOD cnexp
		ileak = gleak*(v-eleak)
		ina = gna*m*m*m*h*(v-ena)
		ikd = gkd*n*n*n*n*(v-ekd)	
		im = gm*p*(v-ekd)
		il = gl*q*q*r*(v-eca)
		it = gt*sinf*sinf*u*(v-eca)
		i=ileak+ina+ikd+im+il
}

DERIVATIVE states { 
	rates(v)
	   m' = (minf-m)/mtau
	   h' = (hinf-h)/htau
	   n' = (ninf-n)/ntau
	   p' = (pinf-p)/ptau
	   q' = (qinf-q)/qtau
	   r' = (rinf-r)/rtau
	   u' = (uinf-u)/utau
}

INITIAL { 
	rates(v)
	h = hinf
	m = minf
	n = ninf
	p = pinf
	q = qinf
	r = rinf
	s = sinf
	u = uinf
	}
	
PROCEDURE rates(v (mV)) {
LOCAL  a, b
UNITSOFF

a = -0.32*(v-Vt-13)/(exp(-(v-Vt-13)/4)-1)
b = 0.28*(v-Vt-40)/(exp((v-Vt-40)/5)-1)
minf = a/(a+b)
mtau = 1/(a+b)

a = 0.128*exp(-(v-Vt-17)/18)
b = 4/(1+exp(-(v-Vt-40)/5))
hinf = a/(a+b)
htau = 1/(a+b)

a = -0.032*(v-Vt-15)/(exp(-(v-Vt-15)/5)-1)
b = 0.5*exp(-(v-Vt-10)/40)
ninf = a/(a+b)
ntau = 1/(a+b)

pinf = 1/(1+exp(-(v+35)/10))
ptau = ptaumax/(3.3*exp((v+35)/20)+exp(-(v+35)/20))

a = 0.055*(-27-v)/(exp((-27-v)/3.8)-1)
b = 0.94*exp((-75-v)/17)
qinf = a/(a+b)
qtau = 1/(a+b)

a = 0.000457*exp((-13-v)/50)
b = 0.0065/(exp((-15-v)/28)+1)
rinf = a/(a+b)
rtau = 1/(a+b)

sinf = 1/(1+exp(-(v+2+57)/6.2))
uinf = 1/(1+exp((v+2+81)/4))
utau = 30.8+(211.4+exp((v+2+113.2)/5))/(3.7*(1+exp((v+2+84)/32)))

UNITSON
}
