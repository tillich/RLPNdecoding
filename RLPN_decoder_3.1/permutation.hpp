#ifndef PERMUTATION_MATRIX
#define PERMUTATION_MATRIX
#include <iostream>

using namespace std;

template <int T_N>
struct permutation_matrix{
    
  uint16_t permu[T_N];
  static constexpr int P = T_N;

    void init(){
        for(int i = 0; i < T_N; ++i){
            permu[i] = i;
        }
    }
  permutation_matrix(){
      init();
  };

    uint16_t& operator[](size_t i) {
    return permu[i];
  }

  uint16_t operator[](size_t i) const {
    return permu[i];
  }
  void rand(){
    init();
    int index;
    int temp;
    for(int i = 0; i< T_N; ++i){
        int tmp = std::rand();
        index = tmp % (T_N-i);
        temp = permu[T_N-i-1];
        permu[T_N-i-1] = permu[index];
        permu[index] = temp;
    }
  }

  permutation_matrix inverse(){
    permutation_matrix p_inv;
    for(size_t i = 0; i < T_N; ++i){
      p_inv[permu[i]] = i;
    }
    return p_inv;
  }

    void print(){
    for(size_t i = 0; i < T_N; ++i){
      cout << permu[i] << "   ";
    }
    cout << endl;
  }
  

      template <size_t T_LEN_E>
  void rand_E(size_t *E){
      permutation_matrix<T_N - T_LEN_E> tmp;
      tmp.rand();

      
      for(size_t i = 0; i < T_N; ++i){
        permu[i] = T_N;
      }

      for(size_t i = 0; i < T_LEN_E; ++i){
        permu[E[i]] = i; 
      }

      size_t offset = 0;
      for(size_t i = 0; i < T_N; ++i){
        if(permu[i] == T_N){
          permu[i] = tmp[i-offset]+T_LEN_E;
        }else{
          offset +=1;
        }
      }
  }



};


template <size_t T_P>
struct permutation;

template <size_t T_P>
struct permutation{

  void init(){
    for(int i = 0; i < T_P; ++i){
      permu[i] = 0;
    }
  }

  permutation(){
    init();
  }

  template <size_t T_P1>
  permutation<T_P + T_P1> operator+(permutation<T_P1>& p){
    permutation<T_P + T_P1> new_p;
    for(size_t i = 0; i < T_P; ++i){
      new_p[i] = permu[i];
    }
    for(size_t i = 0; i < T_P1; ++i){
      new_p[i+T_P] = p.permu[i];
    }
    return new_p;
  }

  void print(){
    for(size_t i = 0; i < T_P; ++i){
      cout << permu[i] << "   ";
    }
    cout << endl;
  }

  uint16_t& operator[](size_t i) {
    return permu[i];
  }

  uint16_t operator[](size_t i) const {
    return permu[i];
  }


  uint16_t permu[T_P];

  static constexpr size_t P = T_P;
};











#endif