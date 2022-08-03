#ifndef __INCLUDE__MATH
#define __INCLUDE__MATH

using namespace std;
constexpr uint64_t power(uint64_t a, uint64_t b)
{
    uint64_t res = 1;
    for(int i = 0; i < b ; ++i){
        res *= a;
    }
    return res;
}
// PROBLEM WITH BOOL WHEN EXTRACTING LAST./...
/*
template <size_t n>
struct type_cell_vector{
    using type = typename std::conditional<
        n <= 1,bool,typename std::conditional<
             (n <= 8),uint8_t,typename std::conditional< 
                (n <= 16),uint16_t,typename std::conditional< 
                    (n <= 32),uint32_t,uint64_t>::type>::type>::type>::type;
};*/


template <size_t n>
struct type_cell_vector{
    using type = typename std::conditional<
             (n <= 8),uint8_t,typename std::conditional< 
                (n <= 16),uint16_t,typename std::conditional< 
                    (n <= 32),uint32_t,uint64_t>::type>::type>::type;
};


template <size_t n>
struct type_cell_walsh{
    using type = typename std::conditional<
             (n <= 6),signed char,typename std::conditional< 
                (n <= 14),short,typename std::conditional< 
                    (n <= 30),signed int,long long>::type>::type>::type;
};

/*
template <size_t n>
struct type_cell_vector{
    using type = typename std::conditional<
        n <= 1,bool,typename std::conditional<
             (n <= 8,uint8_t,typename std::conditional< 
                (n <= 16,uint16_t,typename std::conditional< 
                    (n <= 32,uint32_t,uint64_t>::type>::type>::type>::type;
};*/


constexpr inline size_t binomial(size_t n, size_t k) noexcept
{
    return
      (        k> n  )? 0 :          // out of range
      (k==0 || k==n  )? 1 :          // edge
      (k==1 || k==n-1)? n :          // first
      (     k+k < n  )?              // recursive:
      (binomial(n-1,k-1) * n)/k :       //  path to k=1   is faster
      (binomial(n-1,k) * n)/(n-k);      //  path to k=n-1 is faster
}

constexpr inline size_t error_per_syndrome(size_t L1,size_t l) noexcept{
    size_t temp = L1/power(2,l);
    return temp + 1;
}

constexpr inline size_t size_merge(size_t L1, size_t L2,size_t l) noexcept{
    size_t temp = L1*L2/(power(2,l));
    if(temp < 1){
        cout << "ERROR MERGE TOO LITTLE" << endl;
    }
	return temp;
}
/*
constexpr inline auto adaptedType(size_t value){
    if(value <= power(2,0)){
        return (bool) 0;
    }
    if(value <= power(2,8)){
        return (uint8_t) 0;
    }
    if(value <= power(2,16)){
        return (uint16_t) 0;
    }
    if(value <= power(2,32)){
        return (uint32_t) 0;
    }else{
        return (uint64_t) 0;
    }
}*/



#endif
