#ifndef ENUMERATE_HEADER
#define ENUMERATE_HEADER
#include <iostream>
template<size_t depth, size_t depth_max, class Callable, typename T_TYPE_BLOCK, size_t T_ROWS, size_t T_COLS,size_t T_P>
constexpr void _enumerate(t_matrix_col<T_TYPE_BLOCK,T_ROWS,T_COLS>& H, t_vector<T_TYPE_BLOCK,T_ROWS>& tmp,permutation<T_P>& perm,size_t begin, size_t end, Callable&& c){
    if constexpr (depth == 0 ){

        c(tmp,perm);

    } else {
        for (size_t i =begin; i <= end-depth; ++i){
            perm[depth_max-depth] = i;
           
            t_vector<T_TYPE_BLOCK,T_ROWS> tmp1 = tmp ^ H[i];

            _enumerate<depth- 1,depth_max>(H,tmp1,perm,i+1,end,c);
        }
    }
}

template<size_t depth, class Callable, typename T_TYPE_BLOCK, size_t T_ROWS, size_t T_COLS,size_t T_P>
constexpr void enumerate(t_matrix_col<T_TYPE_BLOCK,T_ROWS,T_COLS>& H, t_vector<T_TYPE_BLOCK,T_ROWS>& tmp,permutation<T_P>& perm,size_t begin, size_t end, Callable&& c){
    _enumerate<depth,depth>(H,tmp,perm,begin,end,c);
}


template<size_t depth, class Callable, typename T_TYPE_BLOCK, size_t T_ROWS, size_t T_COLS>
constexpr void enumerate2(t_matrix_col<T_TYPE_BLOCK,T_ROWS,T_COLS>& H, t_vector<T_TYPE_BLOCK,T_ROWS>& tmp,size_t begin, size_t end, Callable&& c){
    permutation<depth> perm;
    _enumerate<depth,depth>(H,tmp,perm,begin,end,c);
}


#endif
