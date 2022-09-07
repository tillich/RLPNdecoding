#ifndef LINEAR_HEADER
#define LINEAR_HEADER

#include <iostream>
#include <cassert>
#include <random>
#include "permutation.hpp"
#include "tools.hpp"
using namespace std;


// Compute the number of required blocks of type T_TYPE_BLOCK for a vector of size n
template<typename T_TYPE_BLOCK>
constexpr size_t nb_Blocks(size_t n){
    if (n % (sizeof(T_TYPE_BLOCK) * 8) == 0){
        return n / (sizeof(T_TYPE_BLOCK) * 8);
    }
    else{
        return 1 + n / (sizeof(T_TYPE_BLOCK) * 8);
    }
}

template<typename T_TYPE_BLOCK>
constexpr size_t offset(size_t n){
    return ((n-1) % (sizeof(T_TYPE_BLOCK)*8));
}


template <typename T_TYPE_BLOCK, size_t T_LEN>
struct t_vector{
    using TYPE_BLOCK = T_TYPE_BLOCK;
    static constexpr size_t nb_block = nb_Blocks<T_TYPE_BLOCK>(T_LEN);


    protected:
    



    public:

    size_t get_k(size_t i) const {
        const int off = sizeof(T_TYPE_BLOCK)*8 - (T_LEN % (sizeof(T_TYPE_BLOCK)*8));
        if (0 == (T_LEN % (sizeof(T_TYPE_BLOCK)*8))){
            return i/(sizeof(T_TYPE_BLOCK)*8);
        }else{
            return (off + i)/(sizeof(T_TYPE_BLOCK)*8);
        }
    }

    size_t get_mod(size_t i) const {
        const int off = (sizeof(T_TYPE_BLOCK)*8) - (T_LEN % (sizeof(T_TYPE_BLOCK)*8));
        return ((off+i) % (sizeof(T_TYPE_BLOCK)*8));
    }


        void init(){
            for(size_t i = 0; i < nb_block; ++i){
                vect[i] = 0;
            }
        }
        t_vector() : vect(){
           //init();
        }
        t_vector(const t_vector& v) : vect() {
            for(size_t i = 0; i < nb_block; ++i){
                vect[i] = 0;
                vect[i] ^= v.vect[i];
            }
            }


        template <size_t T_P>
        t_vector(const permutation<T_P>& p) : vect(){
            for(size_t i = 0; i < T_P; ++i){
                invert(p[i]);
                //set(p[i],true);
            }
        }
        size_t len(){
            return T_LEN;
        }
        
    
        t_vector& operator^=(const t_vector& right){
            for(size_t i = 0; i < nb_block; ++i){
                this->vect[i] ^= right.vect[i];
            }
            return *this;
        }

        t_vector operator^(const t_vector& right){
            t_vector v;
            for(size_t i = 0; i < nb_block; ++i){
                v.vect[i] = this->vect[i] ^ right.vect[i];
            }
            return v;
        }

        t_vector operator&(const t_vector& right){
            t_vector v;
            for(size_t i = 0; i < nb_block; ++i){
                v.vect[i] = this->vect[i] & right.vect[i];
            }
            return v;
        }


        t_vector& operator=(const t_vector& right){
            for(size_t i = 0; i < nb_block; ++i){
                this->vect[i] = right.vect[i];
            }
            return *this;
        }

        template <size_t T_P>
        t_vector& operator=(const permutation<T_P>& p) {
            for(size_t i = 0; i < T_P; ++i){
                this->invert(p[i]);
                //set(p[i],true);
            }
            return (*this);
        }
        T_TYPE_BLOCK get_block (size_t i) const{
            return vect[i];
        }

    void set(size_t i, bool value){
        int k = get_k(i);
        int pos = get_mod(i);

        vect[k] &= ~((T_TYPE_BLOCK)1 << (sizeof(T_TYPE_BLOCK)*8 - 1 - pos));
        if (value){
            vect[k] |= ((T_TYPE_BLOCK)1 << (sizeof(T_TYPE_BLOCK)*8 - 1 - pos));
        }
    }

    void invert(size_t i){
        int k = get_k(i);
        int pos = get_mod(i);

        vect[k] ^= ((T_TYPE_BLOCK)1 << (sizeof(T_TYPE_BLOCK)*8 - 1 - pos));
    }
    
        size_t hamming(){
            size_t res = 0;
            for(size_t i = 0; i < nb_block; ++i){
                res += popcount(vect[i]);
            }
            return res;
        }

        bool dot_product(const t_vector& right){
            size_t res = ((*this) & right).hamming();
            return res % 2;
        }
        bool operator[](size_t i) const {
            size_t k = get_k(i);
            size_t m = get_mod(i);
            T_TYPE_BLOCK mask = ((T_TYPE_BLOCK) 1 << (sizeof(T_TYPE_BLOCK)*8 - 1 - m));
            T_TYPE_BLOCK c = get_block(k);
            return (( c & mask) >> (sizeof(T_TYPE_BLOCK)*8 - 1 - m));
        }

    

        bool operator==(const t_vector<T_TYPE_BLOCK,T_LEN> &other) const {
        for(int i = 0; i < nb_block; ++i){
            if(vect[i] != other.vect[i]){
                return false;
            }
        }
            return true;
        }
        
        bool operator!=(const t_vector<T_TYPE_BLOCK,T_LEN>  &other) const {
        return !(*this == other);
        }
        
        bool operator<(const t_vector<T_TYPE_BLOCK,T_LEN>  &other) const {
        for(int i = 0; i < nb_block; ++i){
            if(vect[i] < other.vect[i]){
                return true;
            }else if (vect[i] > other.vect[i]){
                return false;
            }
        }
            return false;
        }
        
        bool operator>(const t_vector<T_TYPE_BLOCK,T_LEN> &other) const {
        for(int i = 0; i < nb_block; ++i){
            if(vect[i] > other.vect[i]){
                return true;
            }else if (vect[i] < other.vect[i]){
                return false;
            }
        }
            return false;
        }
        bool operator<=(const t_vector<T_TYPE_BLOCK,T_LEN>  &other) const {
        return !(*this > other);
        }
        bool operator>=(const t_vector<T_TYPE_BLOCK,T_LEN>  &other) const {
        return !(*this < other);
        }
        

        uint8_t get_byte(uint16_t i) const {
            if constexpr (nb_block != 1){
                std::cout << "2 ERREUR GET_BYTE NOT IMPLEMENTED : " << std::endl;
            }
            return (vect[0] >> (sizeof(T_TYPE_BLOCK)-i-1)*8);
        }

    
            // SLOW TO -> to be bitsliced
            // Extract last L bits of this into v
        template <size_t T_PADDING, typename T_TYPE_BLOCK_1,size_t L>
        void extract_last(t_vector<T_TYPE_BLOCK_1,L>& v){
            v.init();
            size_t start = T_LEN - L - T_PADDING;
            for(size_t i = 0; i < L; ++i){
                v.set(i,(*this)[start+i]);
            }
        }

        

        void print(){
            for (int i = 0; i < T_LEN ; ++i){
                cout << (uint64_t) (*this) [i]; //<< endl;
            }
            cout << endl;
        }

        // ONLU FOR < 64
        void rand(){
            size_t j = this->get_mod(0);
            /*
            uint64_t r = 0;
            for(size_t i = 0; i < nb_block; ++i){
                r = (((uint64_t) std::rand()) << 33);
                r |= std::rand();
                vect[i] = r;
            }*/
            for(size_t i = 0; i < T_LEN; ++i){
                set(i,std::rand()%2);
            }

            vect[0] = vect[0] << j;
            vect[0] = vect[0] >> j;

        }
        void permute(t_vector& v, permutation_matrix<T_LEN>& p){
       // init();
       // cout << "INITIAL " << endl;
       // v.print();
        for(size_t i= 0; i < T_LEN; ++i){
            set(p[i],v[i]);
            //cout << p[i] << "    ";
            //set(i,v[p[i]]);
        }
        //cout << endl;
        }

        TYPE_BLOCK get_value0(){
            static_assert(nb_block == 1);
            return vect[0];
        }

        void increment(){
            vect[0] +=1;
        }
        
        T_TYPE_BLOCK vect[nb_Blocks<T_TYPE_BLOCK>(T_LEN)];


};

template <typename T_TYPE_BLOCK, size_t T_ROWS, size_t T_COLS, bool T_COL_REPR>
struct t_matrix
{

    using TYPE_VECTOR = t_vector<T_TYPE_BLOCK,T_ROWS>;
     static constexpr size_t COLS = T_COLS;
     static constexpr size_t ROWS = T_ROWS;

     TYPE_VECTOR& operator[](size_t c){
         return mat[c];
     }

    t_matrix() : mat(){

    }
    template <size_t T_START, size_t T_V_PADDING, typename T_TYPE_BLOCK_1,size_t T_ROWS_1, size_t T_COLS_1>
    void extract_last(t_matrix<T_TYPE_BLOCK_1,T_ROWS_1,T_COLS_1,T_COL_REPR>& m){
        for (int i = T_START; i < T_COLS_1 + T_START ; ++i){
            //cout << "COL" << endl;
            //mat[i].print();
            mat[i].template extract_last<T_V_PADDING>(m[i-T_START]);
        }
    }

    void print(){
        bool res;
        for(int i= 0; i < T_ROWS; ++i){
            for(int j = 0; j < T_COLS; ++j){
                res = mat[j][i];
                cout << (uint16_t) res;
            }
            cout << endl;
        }
    }

    t_matrix<T_TYPE_BLOCK,T_COLS,T_ROWS,!T_COL_REPR> inverse_representation(){
        t_matrix<T_TYPE_BLOCK,T_COLS,T_ROWS,!T_COL_REPR> M;
        for(size_t i = 0; i < T_ROWS; ++i){
            for(size_t j = 0; j < T_COLS; ++j){
                M[i].set(j,mat[j][i]);
            }
        }
        return M;
    }

    size_t rows(){
        if constexpr (T_COL_REPR){
            return T_ROWS;
        }else{
            return T_COLS;
        }
    }

    size_t cols(){
        if constexpr (!T_COL_REPR){
            return T_ROWS;
        }else{
            return T_COLS;
        }
    }

    void rand(){
        for(size_t i = 0; i < T_COLS; ++i){
            mat[i].rand();
            
        }
    }
    void permute(t_matrix& A, permutation_matrix<T_COLS>& p){
        for(size_t i= 0; i < T_COLS; ++i){
            mat[p[i]] = A[i];
        }
    }
    template <typename TYPE_BLOCK_V_MUL>
    TYPE_VECTOR operator*(t_vector<TYPE_BLOCK_V_MUL,T_COLS>& v){
        TYPE_VECTOR mul;
        for(size_t i = 0; i < T_COLS; ++i){
            if(v[i]){
                mul ^= mat[i];
            }
        }
        return mul;
    }

    template <size_t T_P>
    TYPE_VECTOR operator*(permutation<T_P>& p){
        TYPE_VECTOR mul;
        for(size_t i = 0; i < T_P; ++i){
                mul ^= mat[p[i]];
        }
        return mul;
    }
    TYPE_VECTOR mat[T_COLS];


};


template <typename T_TYPE_BLOCK, size_t T_ROWS, size_t T_COLS>
using t_matrix_col = t_matrix<T_TYPE_BLOCK,T_ROWS,T_COLS,true>;

template <typename T_TYPE_BLOCK, size_t T_ROWS, size_t T_COLS>
using t_matrix_row = t_matrix<T_TYPE_BLOCK,T_COLS,T_ROWS,false>;




// Parameters 

// In/Out : Matrix  A
// In/Out : Vector  S
// In     : Integer nb_cols_matrix& A

// Output

// A <- U*A with U the matrix such that U*A has its nb_cols first colunms echelonnized
// S <- U*S
//t_vector<T_TYPE_BLOCK,T_ROWS>& S, 
template <typename T_TYPE_BLOCK, size_t T_ROWS, size_t T_COLS>
bool partial_echelon_form(t_matrix_row<T_TYPE_BLOCK,T_ROWS,T_COLS>& A,size_t nb_cols, size_t& stop){
    //assert(A.rows() == S.len());
    assert(nb_cols <= A.rows() && nb_cols <= A.cols());
    t_vector<T_TYPE_BLOCK,T_COLS> tmp;
    size_t current_col = 0;
    while (current_col < nb_cols){
        size_t j = current_col;
        stop = current_col;
        while((j < A.rows()) && !A[j][current_col]){
            j += 1;
        }
        if (j == A.rows()){
            return false;
        }
        if(j != current_col){
            //tmp = A[current_col];

            //A[current_col] = A[j];
            //A[j] = tmp;
            A[current_col] ^= A[j];
            //S[niveauCourant] ^= S[j];
        }
        for(size_t z = 0; z < A.rows(); ++z){
            if((z != current_col) && A[z][current_col]){
                    A[z] ^= A[current_col];
                    //S[z] ^= S[niveauCourant];
            }
        }
        current_col +=1;
        stop = current_col;
    }
    return true;
}


template <typename TYPE_CELL,size_t T_K, size_t T_N, size_t T_P>
t_vector<TYPE_CELL,T_K> mat_mult(t_matrix_col<TYPE_CELL,T_K,T_N>& A,permutation<T_P>& p){
	t_vector<TYPE_CELL,T_K> res;
	for(size_t i = 0; i < T_P; ++i){
		res ^= A[p[i]];
	}
	return res;
}	



//using matrix_col = t_matrix_col<uint64_t>;
#endif






