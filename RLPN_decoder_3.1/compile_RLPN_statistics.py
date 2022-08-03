
"""
Execute this script : python compile_RLPN_statistics.py (after possibly modifying the input parameters that are given below)

This python script compiles and executes main_RLPN_statistics.cpp


We run here Algorithm 3.1 without taking into account the acceptation criteria. We print in a file the relevant information for each iteration (highest value of the Walsh transform, |e_P| the error weight of e on the set P choosen at this iteration).

--- 
INPUT
---

-----
[n,k,t,w,s,u] below are the problem parameters as in Section 3 of the article.

nbCodes : Number of times with apply Algorithm 3.1 (each time with a different random code)

N_iter : Number of iterations of the RLPN decoding Algorithm 3.1 ran for each nbCodes random codes.

-----

---
OUTPUT : A file "n_k_t_w_s_u.py" 
---

n_k_t_w_s_u.py contains two lists : L_Parity and L_TrueRandom.

L_parity is the list containing the execution trace of successive iterations of the RLPN decoder as described in Algorithm 3.1.
L_TrueRandom is the list containing the execution trace of successive iterations of the RLPN decoder but where we replace <e_N,h_N> by a Bernouilli variable which has the same parameter as the one given in Lemma 3.2 where u is given by t-|e_P|

Both lists are multidimentionnal -> L[0],L[1],...,L[nbCode] represent the execution trace of each random code
								-> L[0][0], L[0][1], ..., L[0][nbIt] contains the trace of each iteration of the RLPN decoder.

More precisly L[.][.] is of the form :
[|e_P| : the weight of the error e on the set P (P is different for each iteration)
,|e_N| 
, N : Number of disctinct parity checks created
, Theorical bias
, Practical bias
, Position of e_P in L
, \hat{f}{e_P} : Walsh transform of e_P
, \hat{f}{L[0]} : Highest coefficient of the Walsh transform
, \hat{f}{L[1]} : Second highest coefficient of the Walsh transform
]

Where L denotes the ordered table of Walsh values at each iteration: L[0] is the first (biggest) value, L[1] the second, etc

"""


# CHANGE
[n,k,t,w,s,u] =[100, 20, 24, 4, 16, 19]
[nbIt,nbCodes] = [100,10]

name = "decode_trace.out"
import re
import subprocess
import sys
from math import floor, ceil
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
def compile(N,K,T,U,V,P,name):
	prec = " -D PARAM_N=" + str(N) + " -D PARAM_K=" + str(K) + " -D PARAM_T=" + str(T) + " -D PARAM_U=" + str(U) + " -D PARAM_P=" + str(P) + " -D PARAM_V=" + str(V)
	optim = " -O3 -march=native -fopenmp"
	fil = " main_RLPN_statistics.cpp -std=c++20"
	compiler = "g++"
	command = compiler + fil + optim + prec + " -o" + name
	comp = subprocess.Popen(command,shell=True)
	comp.wait()



def launch(name,nbIt,nbCodes):
	command = "./" +name + " " + str(nbIt) + " " + str(nbCodes)
	comp = subprocess.Popen(command,shell=True)
	comp.wait()


def main(args):
	#[N,K,T,W,U,V,nbIt,nbCodes,nbItDumer] =[100, 20, 24, 4, 16, 19,100,10,40]
	#[200, 16, 67, 2, 15, 54,1,10000,40,0]
	#[200, 16, 67, 2, 13, 59,1,10000,20,0]
	#[2000, 20, 832, 2, 15, 826,1,10000,20,0]
	# [201, 40, 42, 6, 19, 37,1,100,10,10]
	#[200, 16, 67, 2, 14, 59,1, 500,100,0]
	#[300, 24, 95, 3, 17, 90,1,1000,1,0]
	#[200, 16, 57, 2, 14, 46,1,10000,20]
	#[200, 16, 67, 2, 14, 59,1, 500,100]
	#[301, 24, 100, 3, 17, 90,1,1000,20]
	#[400, 32, 114, 4, 18, 102,1,10000,20,50] #LAUNCH 1
	#[400, 32, 114, 4, 18, 108,1,100,1,50]     
	#[351, 28, 97, 4, 15, 92,1,100,1,50]
	#[401, 32, 114, 4, 17, 103,8888,100,1,1]
	#[401, 32, 114, 4, 17, 104,8888,1000,1,10]
	#[401, 28, 130, 4, 17, 122,950,300,10,4]

	#[401, 28, 118, 4, 17, 108,1400]
	#[350, 24, 101, 2, 20, 89,0.95]
	#[280, 19, 77, 2, 16, 68,0.9]
	#[200, 14, 49, 2, 10, 42,0.95] 
	#[481, 33, 146, 4, 21, 135]
	#[401, 28, 118, 4, 17, 108] # GOOD PARAM BECAUSE T LESS (-20) 
	#[401, 28, 138, 4, 19, 130]
	#[701, 21, 278, 2, 19, 271]
	#[201, 30, 55, 6, 21, 46]
	#[120, 18, 33, 4, 14, 22]
	#[120, 18, 33, 4, 16, 23]
	#[81, 12, 22, 4, 11, 13]
	#[151, 22, 41, 4, 11, 30]
	
	#[201, 20, 63, 4, 15, 53]
		
	#[200, 20, 63, 4, 16, 53]
	
	#[200, 30, 55, 6, 18, 45]
	#[250, 25, 59, 4, 12, 53]

	#[202, 20, 44, 4, 10, 34]


	#t = int(0.25 * n)
	print("PARAMS : ")
	print("N : " + str(n))
	print("K : " + str(k))
	print("T : " + str(t))

	print("U : " + str(s))
	print("V : " + str(u))


	print("W : " + str(w))
	compile(n,k,t,s,u,w,name)
	print("EXEC")

	launch(name,nbIt,nbCodes)

if __name__ == '__main__':
    import sys
    sys.exit(main(sys.argv))
