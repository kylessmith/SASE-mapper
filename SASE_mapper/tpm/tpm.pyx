cdef extern from "stdlib.h" nogil:
	double	drand48()
	void	srand48(long)
cdef extern from "math.h":
	double log(double x)
	double pow(double x, double y)
	double exp(double x)
cdef enum :
	TableSize = 1500
	


cdef double GamLn(double x):
	cdef double ser = 1.000000000190015
	cdef double *Coe = [76.18009172947146, -86.50532032941677, 24.01409824083091,
        						-1.231739572450155, 0.1208650973866179e-2, -0.5395239384953e-5]
	cdef int j			
	for j in xrange(6): ser += Coe[j] / (x + j + 1.0)
	return log(ser/x) - x + (x + 0.5)*log(x + 5.5) - 4.581061466795327


cdef double Lfa(int n):
	cdef double FacTable[TableSize]
	
	if n <= 1: return 0.0
	elif n < TableSize:
		return FacTable[n-2] if FacTable[n-2] else GamLn(n+1.0)
	else: return GamLn(n+1.0)



cdef double DoubleSumEquSimple(double w, int L):
	cdef double innersum = 0.0
	cdef int k=L
	cdef int s
	cdef double A

	for s in xrange(k):
		A = -log(w)
		innersum += pow(A, s) / exp(Lfa(s))
	innersum = w*innersum

	return innersum


cdef double DoubleSumEqu(double w, int L, double tau):
	if tau == 1: return DoubleSumEquSimple(w, L)
	
	cdef double r = 0.0
	cdef int k, s
	cdef double r1, r2, innersum, A
	
	for k in xrange(1,L+1):
		r1 = Lfa(L) - Lfa(L-k) - Lfa(k)
		r2 = (L-k)*log(1-tau);
		innersum = 0.0
		if w <= pow(tau, k):
			for s in xrange(k):
				A = k*log(tau) - log(w)
				innersum += pow(A, s) / exp(Lfa(s))
			innersum = w*innersum
		else:
			innersum = pow(tau, k)
		r += exp(r1+r2) * innersum
		
	return r


cdef double empir_mult(long sd, double t, long n, double ta, long o):
	srand48(sd)
	cdef long j,v,c=0
	cdef double ct,p
	
	for j in xrange(o):
		ct = 0
		for v in xrange(n):
			p = drand48()
			if p <= ta: ct -= log(p)
		if ct != 0 and ct >= t: c+=1
		
	return float(c)/o

def tpm(double ptau, int ploops, int pseed, list pvals):
	cdef int upL = 100
	cdef double tau = ptau
	cdef long loops = ploops
	cdef unsigned long sd = pseed
	
	cdef double p=1.0, t=0
	cdef int L=len(pvals), k=0, i
		
	for i in xrange(L):
		if pvals[i] <= tau:
			k+=1
			p *= pvals[i]
			t -= log(pvals[i])

	if k:
		if L <= upL:
			cout = DoubleSumEqu (p, L, tau)
		else:
			cout = empir_mult(sd, t, L, tau, loops)
			
		return cout
	
	else:
		return 1
	