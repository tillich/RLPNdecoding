#The parameters are given in the following form : [w,s,k,n,u,t]
#This script generates histogram_w_s_k_n_u_t.pdf and histogram_w_s_k_n_u_t_zoom.pdf from the file w_s_k_n_u_t.py.
#The files w_s_k_n_u_t.py should be previously generated with compile_RLPN_histogram.py

paramG = [[[2,12,12,170,14,18],
[2,12,19,1787,161,165],
[4,12,19,92,3,7],
[2,12,26,20079,1807,1811],
[4,12,26,278,12,16],
[6,12,26,87,2,6],
[4,12,33,901,41,45],
[6,12,33,178,5,9]],
[[2,16,16,644,151,156],
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
[10,16,37,87,3,8]],
[[2,19,19,1794,573,579],
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
[10,19,40,105,7,13]]
]

binwidthG = [8,32,128]
import matplotlib.pyplot as plt
import importlib
import numpy
import math

from gmpy2 import bincoef, log2,log,fac
from gmpy2 import mpz
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
		return find_GV1(k,n)
	else:
		return find_GV2(k,n)


nbBin = 0

for params in paramG:
	binwidth = binwidthG[nbBin]
	nbBin += 1
	file_name = []
	for param in params:
		#[_n,_k,_t,_s,_w,_u] = param
		[_w,_s,_k,_n,_u,_t] = param
		st = str(_w)+"_"+str(_s)+"_"+str(_k)+"_"+str(_n)+"_"+str(_u)+"_"+str(_t)
		file_name.append(st)

	for name in file_name:
		data = importlib.import_module(name)


		n = data.n
		k = data.k
		t = data.t
		s = data.s
		w = data.w
		u = data.u

		L_parity = data.L_Parity
		L_TrueRandom = data.L_TrueRandom


		nbDistinctParity_P = L_parity[0][0][0][2]

		nbError = [[],[]]
		lenCode = []
		WalshGV = []
		WalshSolution = [[],[]]
		WalshTheoricSolution = []
		WalshBound = []
		countSuperiorBound = [[],[]]
		Kratio = []


		minValue = math.inf
		maxValue = -math.inf

		minValueZoom = math.inf
		maxValueZoom = -math.inf


		x = [[],[]]
		zoom = [[],[]]
		SecondHighest = [[],[]]
		L = [L_parity[0], L_TrueRandom[0]]
		for z in L_parity[0]:
			_lenCode = z[0][2]
			lenCode.append(_lenCode)
			GV_C = find_GV(s,_lenCode)
			_WalshGV = _lenCode - 2*GV_C
			WalshGV.append(_WalshGV)

			eps = kraw(n-s,w,u)/bincoef(n-s,w)
			theoricError = (1-eps)/2

			_WalshTheoricSolution = _lenCode-2*theoricError*_lenCode
			WalshTheoricSolution.append(_WalshTheoricSolution)
			#Kratio.append(((_lenCode)* (eps**2))/s)
			Kratio.append(_WalshTheoricSolution/_WalshGV)
			WalshBound.append((_WalshGV+_WalshTheoricSolution)/2)

		for l in range(0,len(L)):
			for i in range(0,len(L[l])):

				z = L[l][i]
				x[l].append(z[1])
				mi = min(z[1])
				ma = max(z[1])
				minValue = min(mi,minValue)
				maxValue = max(ma,maxValue)
				
				nbError[l].append(z[0][4])
				WalshSolution[l].append(z[0][6])
				
				j = 0
				while(z[1][j] >= WalshBound[i]):
					j+= 1
				countSuperiorBound[l].append(j)

				j = 0
				zoom[l].append([])
				while(z[1][j] >= 0.6*WalshGV[i]):
					zoom[l][i].append(z[1][j])
					j+=1
				
				mi = min(zoom[l][i])
				ma = max(zoom[l][i])
				minValueZoom = min(mi,minValueZoom)
				maxValueZoom = max(ma,maxValueZoom)
				SecondHighest[l].append(z[1][1])
				

				





		#bins=numpy.arange(minValue, maxValue + binwidth, binwidth)
		bins = numpy.arange(minValue, 2*WalshGV[0] + binwidth, binwidth)
		fig, axs = plt.subplots(2, 2)

		for i in range(0,2):
			for j in range(0,2):
				txt = 'local min'
				idx = i + 2*j
				#axs[i,j].set_title(r'$\mathcal{{H}} = {0}, \; \frac{{ \#\mathcal{{H}} \epsilon^2}}{{s}}={1}, \; \mbox{{Experimental bias}}={2}, \mbox{{Theoric walsh of word at distance GV}}={3}, \; \mbox{{Theoric walsh value of e_P}}={4} ,\; \mbox{{Walsh value of e_P}}={5} \; ; \; {6}$'.format(lenCode[idx],Kratio[idx],round(nbError[0][idx] / lenCode[idx],5),WalshGV[idx],WalshTheoricSolution[idx],WalshSolution[0][idx],WalshSolution[1][idx] ))
				#axs[i,j].set_title((r'$\#\mathcal{{H}} = {0}, \; \frac{{ \#\mathcal{{H}} \epsilon^2}}{{s}}={1} $' + ' Theoric Walsh value of E_P={2}' + ' Walsh value of E_P={3}').format(lenCode[idx],round(Kratio[idx],1),round(WalshTheoricSolution[idx],0),WalshSolution[0][idx]))
				#txt = (r'$\#\mathcal{{H}} = {}, \; \frac{{ \#\mathcal{{H}} \epsilon^2}}{{s}}={} $' + ', Theorical Walsh value of solution:{:.1f}' + ' \n Experimental Walsh value of solution: {} (Parity Checks) ; {} (BSC)' + '\n Second highest walsh coefficient: {} (Parity Checks) ; {} (BSC)').format(lenCode[idx],int(Kratio[idx]),round(WalshTheoricSolution[idx],0),WalshSolution[0][idx],WalshSolution[1][idx],SecondHighest[0][idx],SecondHighest[1][idx])
				txt = (r'$\#\mathcal{{H}} = {}, \;$' + 'Theoretical values : ' +r'$\frac{{\mathcal{{F}}(\epsilon)}}{{\mathcal{{F}}(GV)}}={}, \; \mathcal{{F}}(\epsilon) = {}, \; \mathcal{{F}}(GV) = {}$' + ' \n Experimental values :' + r'$\mathcal{{F}}(e_{{P}})$' +' :  {} (Parity Checks) ; {} (BSC)' + '\n Second highest walsh coefficient: {} (Parity Checks) ; {} (BSC) \n Number of Walsh coefficient greater than ' + r'$\frac{{\mathcal{{F}}(GV) + \mathcal{{F}}(\epsilon)  }}{{2}}$' + ': {} (Parity Checks) ; {} (BSC)').format(lenCode[idx],int(Kratio[idx]),int(WalshTheoricSolution[idx]),int(WalshGV[idx]),WalshSolution[0][idx],WalshSolution[1][idx],SecondHighest[0][idx],SecondHighest[1][idx],countSuperiorBound[0][idx],countSuperiorBound[1][idx])
				#txt = r'$\#\mathcal{{H}} = {}, \; \frac{{ \#\mathcal{{H}} \epsilon^2}}{{s}}={} \\mbox{{IHI}}$'.format(1,2)
				#axs[i,j].text(minValue,0,txt)
				axs[i,j].set_title(txt)
				axs[i,j].hist(x[0][idx], bins, alpha=0.5, label='Parity Checks')
				axs[i,j].hist(x[1][idx], bins, alpha=0.5, label='BSC')
				if(i == 0 and j == 0):
					axs[i,j].legend(loc='upper right')
		#plt.hist(walshId, bins, alpha=0.5, label='IdealWalsh')


		#plt.legend(loc='upper right')

		plt.rcParams['text.usetex'] = True
		#plt.suptitle(r'$n={0}, \; k={1}, \; s={2}, \; w={3}, \; u={4},\; t={5}, \; \frac{{ \#\mathcal{{H}} \epsilon^2}}{{s}}={6} \\ \\ \frac{{1-\epsilon}}{{2}}={7}, \; \mbox{{Experimental bias}}={8}, \; \#\mathcal{{H}}={9}, \; GV(s,\#\mathcal{{H}})={10} \\ \\  \mbox{{Theoric walsh of word at distance GV}}={11}, \; \mbox{{Theoric best walsh}}={12} ,\; \mbox{{Best walsh value Parity Checks}}={13}, \; \mbox{{Best walsh value BSC}}={14},\; \mbox{{Average 2...11 highest walsh values of Parity Checks}}={15}, \; \mbox{{Average of 2...11 highest walsh values of BSC}}={16} $'.format(n,k,s,w,u,t,Kratio,theoricError,moyError,nbDistinctParity_P,GV_H,theoricWalsh,theoricWalshBest,bestWalshValue_x1,bestWalshValue_x2,best_10_x1,best_10_x2))


		
		plt.suptitle(r'$w = {0}, \; s = {1} \; k = {2}, \; n={3}, \: |e_P| = {4}, \; |e_N| = {5}, \; \; \; \frac{{1-\epsilon}}{{2}}={6:.6f}$'.format(w,s,k,n,t-u,u,round(theoricError,6)))
		#plt.suptitle(r'$n={0}, \; k={1}, \; s={2}, \; w={3}, \; u={4},\; t={5}, \; \; \; \frac{{1-\epsilon}}{{2}}={6:.6f}$'.format(n,k,s,w,u,t,round(theoricError,6)))


		#plt.suptitle(r'n= is Number $\displaystyle\sum_{n=1}^\infty' r'\frac{-e^{i\pi}}{2^n}$!', fontsize=16, color='r')


		#plt.show()
		[n,k,t,s,w,u]
		name = "histogram_"+str(w)+"_"+str(s)+"_"+str(k)+"_"+str(n)+"_"+str(u)+"_"+str(t)+".pdf"
		#print(plt.rcParams["backend"])
		plt.tight_layout()
		plt.gcf().set_size_inches(20, 9)
		fig.savefig(name,dpi=fig.dpi,bbox_inches='tight')
		plt.close()


		fig2, axs2 = plt.subplots(2, 2)


		#bins=numpy.arange(minValueZoom, maxValueZoom + binwidth, binwidth)
		bins = numpy.arange(minValueZoom, 2*WalshGV[0] + binwidth, binwidth)

		for i in range(0,2):
			for j in range(0,2):
				idx = i + 2*j
				#axs[i,j].set_title(r'$\mathcal{{H}} = {0}, \; \frac{{ \#\mathcal{{H}} \epsilon^2}}{{s}}={1}, \; \mbox{{Experimental bias}}={2}, \mbox{{Theoric walsh of word at distance GV}}={3}, \; \mbox{{Theoric walsh value of e_P}}={4} ,\; \mbox{{Walsh value of e_P}}={5} \; ; \; {6}$'.format(lenCode[idx],Kratio[idx],round(nbError[0][idx] / lenCode[idx],5),WalshGV[idx],WalshTheoricSolution[idx],WalshSolution[0][idx],WalshSolution[1][idx] ))
				#axs[i,j].set_title((r'$\#\mathcal{{H}} = {0}, \; \frac{{ \#\mathcal{{H}} \epsilon^2}}{{s}}={1} $' + ' Theoric Walsh value of E_P={2}' + ' Walsh value of E_P={3}').format(lenCode[idx],round(Kratio[idx],1),round(WalshTheoricSolution[idx],0),WalshSolution[0][idx]))
				txt = ('Walsh transform of a word at distance GV: ' + r'$\mathcal{{F}}(GV):$' + ' {:.1f} \n Number Walsh coefficient greater than ' + r'$\frac{{\mathcal{{F}}(GV) + \mathcal{{F}}(\epsilon)  }}{{2}}:$' + ' {} (Parity Checks) ; {} (BSC)').format(WalshGV[idx],countSuperiorBound[0][idx],countSuperiorBound[1][idx])
				#txt = r'$\#\mathcal{{H}} = {}, \; \frac{{ \#\mathcal{{H}} \epsilon^2}}{{s}}={} \\mbox{{IHI}}$'.format(1,2)
				#axs[i,j].text(minValue,0,txt)
				axs2[i,j].set_title(txt)
				axs2[i,j].hist(zoom[0][idx], bins, label="Parity Checks" , alpha=0.5)
				axs2[i,j].hist(zoom[1][idx], bins, label="BSC", alpha=0.5)
				if(i == 0 and j == 0):
					axs2[i,j].legend(loc='upper right')

		
		
		plt.suptitle((r'$w = {0}, \; s = {1} \; k = {2}, \; n={3}, \: |e_P| = {4}, \; |e_N| = {5}, \; \; \; \frac{{1-\epsilon}}{{2}}={6:.6f},\;\;$' + '   Tail distribution ' + r'$0.6*\mathcal{{F}}(GV) $' ).format(w,s,k,n,t-u,u,round(theoricError,6)))
		#plt.suptitle((r'$n={0}, \; k={1}, \; s={2}, \; w={3}, \; u={4},\; t={5}, \; \; \; \frac{{1-\epsilon}}{{2}}={6:.6f}$' + '  Tail distribution (60 Percent GV)').format(n,k,s,w,u,t,round(theoricError,6)))
		
		#name = "histogram_"+str(n)+"_"+str(k)+"_"+str(t)+"_"+str(s)+"_"+str(w)+"_"+str(u)+"_zoom.pdf"
		name = "histogram_"+str(w)+"_"+str(s)+"_"+str(k)+"_"+str(n)+"_"+str(u)+"_"+str(t)+"_zoom.pdf"
		#print(plt.rcParams["backend"])
		plt.tight_layout()
		plt.gcf().set_size_inches(20, 9)
		fig2.savefig(name,dpi=fig.dpi,bbox_inches='tight')
		plt.close()