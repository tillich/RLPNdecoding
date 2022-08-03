#include <cassert>
#include <cstring>
#include <iostream>
#include "linear_algebra.hpp"
template<typename T>
static constexpr inline T pow(T x, unsigned p) {
    T result = 1;

    while (p) {
        if (p & 0x1) {
            result *= x;
        }
        x *= x;
        p >>= 1;
    }

    return result;
}

using namespace std;
template <typename eval,typename bit>
void walsh(bit *f, eval* e, size_t n){
	assert(n < 64);
	//assert(sizeof(eval)*4 - 2 >= n);
	eval *ep = new eval[(size_t) pow(2,n)];
	for(size_t i = 0; i < pow(2,n); ++i){
		/*if( f[i] == 0){
			e[i] = 1;
		}
		if(f[i]==1){
			e[i] == -1;
		}
		if(f[i] == 2){
			e[i] == 0;
		}*/
		if(f[i] == 0){
			e[i] = 0;
		}
		else if (f[i] == 1){
			e[i] = -1;
		}
		else if(f[i] == 2){
			e[i] = 0;
		}
		//e[i] = (f[i] == 0) ? 1 : -1;
	}
	for(size_t k = 1; k <= n; ++k){
		size_t kpow = pow(2,k);
		for(size_t i = 0; i < pow(2,n-k) ; ++i){
			size_t kpowi = kpow*i;
			
			for(size_t j = 0; j < pow(2,k-1); ++j){

				ep[kpowi + j] = e[kpowi+j] + e[kpowi+j + kpow/2];
				ep[kpowi + kpow/2 + j] = e[kpowi+j] - e[kpowi+j + kpow/2];
			}
		}
		//copy ep into e
		//copy(ep,e);
		memcpy(e,ep,sizeof(eval)*((size_t)  pow(2,n)));
	}

	delete[] ep;
}

template <size_t T_U>
long int walsh_naive_coeff(uint8_t *f, t_vector<uint64_t,T_U>& a){
	t_vector<uint64_t,T_U> x;
	long int w = 0;
	for(size_t i = 0; i < pow(2,T_U); ++i){
		bool res = a.dot_product(x) ^ ((bool) f[i]);
		if(res == 0){
			w += 1;
		}else{
			w -= 1;
		}
		x.increment();
	}
	return w;
}