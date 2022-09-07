import numpy
import math
from gmpy2 import bincoef, log2,log
from gmpy2 import mpz,mpfr

from scipy.optimize import brentq
def H2_G(x,a):
	if(x == 0 or x == 1):
		return -a
	return -x*log2(x) - (1-x)*log2(1-x) - a
	
def H2(x):
	assert x >= 0
	assert x <= 1
	return H2_G(x,0)

def H2_I(x):
	if(x == 0):
		return 0

	if(x == 1):
		return 0.5
	
	return brentq(H2_G,a=0,b=0.5,args=(x))
def kraw(N,W,X):
	res = 0
	for j in range(0,W+1):
		res += mpz(((-1)**j))*bincoef(mpz(X),mpz(j))*bincoef(mpz(N-X),mpz(W-j))
	return res
def eps(N,W,X):
	return kraw(N,W,X)/bincoef(N,W)

def find_GV1(k,n):
	#start = int(H2_I(1-k/n)*n)
	#start = int(0.95*H2_I(1-k/n)*n)
	start = 0
	tot = 0
	dual = 2**(n-k)
	for i in range(start,n):
		tot += bincoef(n,i)
		if(tot > dual):
			return i
def binom(n,k):
	return math.sqrt(1/(2*math.pi*n*(k/n)*(1-k/n)))*(2**(H2(k/n)*n))

def find_GV2(k,n):
	start = 1
	tot = 1
	dual = 2**(n-k)
	for i in range(start,n):
		tot += binom(n,i)
		if(tot > dual):
			return i


def find_GV(k,n):
	if (n <= 2**16):
		return find_GV2(k,n)
	else:
		return find_GV2(k,n)




def find_u(n,k,s,w):
	if (n-s < 1):
		return -1
	nbEq = bincoef(n-s,w)/2**(k-s)
	lenH = 2**s
	lenCode = int(lenH*(1-(1-1/lenH)**nbEq))

	GV_lenCode = find_GV(s,lenCode)
	Walsh_GV = lenCode - 2*GV_lenCode
	for u in range(int((n-s)/2),1,-1):
		epsil = eps(n-s,w,u)
		probFlip = (1-epsil)/2
		Walsh_Theorical = lenCode - 2*probFlip*lenCode
		K = Walsh_Theorical/Walsh_GV

		if K > s:
			return u

def find_n(k,s,w):
	for n in range(k,2**16-1):
		nbEq = bincoef(n-s,w)/2**(k-s)
		lenH = 2**s
		K = nbEq/(2**s)
		if K > 3:
			return n
	return -1

for s in [12,16,19]:
	for k in [s,s+7,s+14,s+21]:
		for w in [2,4,6,8,10]:
			n = find_n(k,s,w)
			if(n != -1):
				u = find_u(n,k,s,w)
				if(u != -1 and u != None):
					sb = math.floor(s/3)
					t = sb+u
					name = "[" + str(w) + "," + str(s) + "," + str(k) + "," + str(n) + "," + str(u) + "," + str(t) +"],"
					print(name)

