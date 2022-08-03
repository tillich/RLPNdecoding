#ifndef TOOLS_HEADER
#define TOOLS_HEADER

#include "linear_algebra.hpp"
#include "permutation.hpp"
#include <iostream>
#include <cmath>
using namespace std;

inline int popcount(uint8_t x){
  return __builtin_popcount(x);
}
inline int popcount(uint16_t x){
  return __builtin_popcount(x);
}
inline int popcount(uint64_t x){
  return __builtin_popcountll(x);
}
inline int popcount(uint32_t x){
  return __builtin_popcountl(x);
}

inline size_t binom(size_t n, size_t k) noexcept
{
    return
      (        k> n  )? 0 :          // out of range
      (k==0 || k==n  )? 1 :          // edge
      (k==1 || k==n-1)? n :          // first
      (     k+k < n  )?              // recursive:
      (binom(n-1,k-1) * n)/k :       //  path to k=1   is faster
      (binom(n-1,k) * n)/(n-k);      //  path to k=n-1 is faster
}


double kraw(size_t N,size_t W,size_t X){
  double res = 0;
  for(size_t i = 0; i <= W; ++i){
    res += std::pow(-1,(int) i)*binom(X,i)*binom(N-X,W-i);
  }
  return res;
}

double varepsilon(size_t LEN,size_t W,size_t V){


  double k = ((double)kraw(LEN,W,V) )/ ((double)binom(LEN,W));
  return k;
}



#endif