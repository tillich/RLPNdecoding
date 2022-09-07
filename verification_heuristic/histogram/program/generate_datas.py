# PARAMS_TO_GENERATE contains the list of the simulations to be made. It will create a file "w_s_k_n_u_t.py"
# Parameters given in the following form [w,s,k,n,u,t]

PARAMS_TO_GENERATE=[[2,12,12,170,14,18],
[2,12,19,1787,161,165],
[4,12,19,92,3,7],
[2,12,26,20079,1807,1811],
[4,12,26,278,12,16],
[6,12,26,87,2,6],
[4,12,33,901,41,45],
[6,12,33,178,5,9],
[2,16,16,644,151,156],
[4,16,16,65,6,11],
[2,16,23,7111,1712,1717],
[4,16,23,175,21,26],
[6,16,23,70,4,9],
[8,16,23,52,2,7],
[4,16,30,545,73,78],
[6,16,30,134,11,16],
[8,16,30,78,4,9],
[10,16,30,62,2,7],
[4,16,37,1792,248,253],
[6,16,37,277,25,30],
[8,16,37,127,7,12],
[10,16,37,87,3,8],
[2,19,19,1794,573,579],
[4,19,19,99,15,21],
[6,19,19,54,4,10],
[8,19,19,46,2,8],
[4,19,26,285,53,59],
[6,19,26,94,10,16],
[8,19,26,64,4,10],
[10,19,26,55,2,8],
[4,19,33,908,179,185],
[6,19,33,185,23,29],
[8,19,33,98,8,14],
[10,19,33,74,4,10],
[4,19,40,3004,605,611],
[6,19,40,387,53,59],
[8,19,40,161,15,21],
[10,19,40,105,7,13]
]
[nbIt,nbCodes] = [4,1]


name = "generate_data.out"
import re
import subprocess
import sys
from math import floor, ceil
from gmpy2 import log2, bincoef
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
	optim = " -O3 -march=native"
	fil = " main_RLPN_histograms.cpp -std=c++20"
	compiler = "g++"
	command = compiler + fil + optim + prec + " -o" + name
	comp = subprocess.Popen(command,shell=True,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	ou,err = comp.communicate()
	if err:
		print(err)
		print("ERROR COMPILATION !")
		exit()
	if ou:
		print("ERROR COMPILATION !")
		exit()



def launch(name,nbIt,nbCodes):
	command = "./" +name + " " + str(nbIt) + " " + str(nbCodes)
	comp = subprocess.Popen(command,shell=True)
	comp.wait()

def main(args):
	for param in PARAMS_TO_GENERATE:
		w,s,k,n,u,t = param
		print("PARAMS : ")
		print("n : " + str(n))
		print("k : " + str(k))
		print("t : " + str(t))

		print("s : " + str(s))
		print("u : " + str(u))
		print("w : " + str(w))
		compile(n,k,t,s,u,w,name)
		print("EXEC")
		launch(name,nbIt,nbCodes)

if __name__ == '__main__':
    import sys
    sys.exit(main(sys.argv))
