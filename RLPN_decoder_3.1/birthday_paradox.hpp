#ifndef BIRTHDAY_PARADOX_HEADER
#define BIRTHDAY_PARADOX_HEADER
template <size_t T_K, size_t T_N, size_t T_W>
struct birthday_paradox{


    static constexpr size_t N = T_N;
    static constexpr size_t K = T_K;

    static constexpr size_t N1 = N/2;
    static constexpr size_t N2 = N-N1;
    static constexpr size_t P1 = T_W/2;
    static constexpr size_t P2 = T_W-P1;

    static constexpr size_t LEN_L1 = binomial(N1,P1);
    static constexpr size_t LEN_L2 = binomial(N2,P2);

    
    using TYPE_CELL_VECTOR = typename type_cell_vector<T_K>::type;
    using TYPE_VECTOR = t_vector<TYPE_CELL_VECTOR,T_K>;
    using TYPE_MATRIX = t_matrix_col<TYPE_CELL_VECTOR,T_K,T_N>;

    using TYPE_VALUE_1 = permutation<P1>;
    using TYPE_VALUE_2 = permutation<P2>;
    using TYPE_LIST_1 = list_SOA_fixed_length<TYPE_VECTOR,TYPE_VALUE_1,LEN_L1>;
    using TYPE_LIST_2 = list_SOA_fixed_length<TYPE_VECTOR,TYPE_VALUE_2,LEN_L2>;

    birthday_paradox(TYPE_MATRIX& _H, TYPE_VECTOR& _S) : H(_H), S(_S) {

    }
    
    void solve(auto&& callback){
        
        TYPE_LIST_1& var_L1 = L1;
        TYPE_LIST_2& var_L2 = L2;



        L1.init();
        L2.init();

        auto add_L1 = [&var_L1](TYPE_VECTOR& key, TYPE_VALUE_1 & value){
            var_L1.insert(key,value);
        };
        auto add_L2 = [&var_L2](TYPE_VECTOR& key, TYPE_VALUE_2 & value){
            var_L2.insert(key,value);
        };

        permutation<P1> perm1;
        TYPE_VECTOR S0;
        enumerate<P1>(H,S0,perm1,0,N1,add_L1);
        L1.sort();

        permutation<P2> perm2;
        enumerate<P2>(H,S,perm2,N1,N,add_L2);

        L2.sort();

        auto callback1 = [&callback](TYPE_VALUE_1& v1, TYPE_VALUE_2& v2){
           auto perm_concat = v1+v2;
            callback(perm_concat);  
        };
        naive_merge_list(L1,L2,callback1);
    }

    TYPE_MATRIX& H;
    TYPE_VECTOR& S;

    TYPE_LIST_1 L1;
    TYPE_LIST_2 L2;


};

#endif