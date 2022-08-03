#ifndef __INCLUDE__ListTriAverage
#define __INCLUDE__ListTriAverage

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

#include "math.h"

#include <string.h> 

#include "linear_algebra.hpp"

#include <iostream>
using namespace std;


template <typename T_TYPE_KEY, typename T_TYPE_VALUE, size_t T_LEN>
struct list_SOA_fixed_length{


    using TYPE_THIS = list_SOA_fixed_length<T_TYPE_KEY,T_TYPE_VALUE,T_LEN>;
    using TYPE_KEY = T_TYPE_KEY;
    using TYPE_VALUE = T_TYPE_VALUE;

    T_TYPE_KEY *key;
    T_TYPE_VALUE *value;

    size_t count;


    static constexpr bool KEY_VALUE_SEPARATED = true;
    static constexpr size_t MAX_LEN = T_LEN;

    list_SOA_fixed_length(){
      key = new T_TYPE_KEY[T_LEN];
      value = new T_TYPE_VALUE[T_LEN];
      count = 0;
    }

    ~list_SOA_fixed_length(){
     delete[] key;
     delete[] value;
      count = 0;
    }

    size_t len(){
		return count;
	}
    size_t max_len(){
		return T_LEN;
	}

    void init(){
        count = 0;
    }
    TYPE_THIS& operator=(TYPE_THIS& L){
        count = L.len();
        for(size_t i = 0; i < count; ++i){
            value[i] = L.get_value(i);
            key[i] = L.get_key(i);
        }
        return *this;
    }
    


    T_TYPE_VALUE& get_value(size_t i){
        return value[i];
    }
    T_TYPE_KEY& get_key(size_t i){
        return key[i];
    }
    void insert(T_TYPE_KEY& k, T_TYPE_VALUE& v){
        key[count] = k;
        value[count] = v;
        count+=1;
    }
    void sort(){
        //cout << "GLO "<< endl;
        const int RADIX = 8;
        const int nbCount = power(2,RADIX);

        using TYPE_COUNT = size_t;
        TYPE_COUNT radix_count[nbCount];
        TYPE_THIS list_aux;
        list_aux.count = count;
        uint8_t rad;

        bool b_current = true;

        const int numberRadix = TYPE_KEY::nb_block * sizeof(typename TYPE_KEY::TYPE_BLOCK);
        for(int pass = numberRadix - 1; pass >= 0; --pass){

            TYPE_THIS& current = b_current ? *this : list_aux;
            TYPE_THIS& aux = b_current ? list_aux : *this;
            b_current = !b_current;


            //memset(radix_count,0,sizeof(TYPE_COUNT)*nbCount);
            for(int i = 0; i < nbCount; ++i){
                radix_count[i] = 0;
            }
            for (size_t i = 0; i < len(); ++i){
                rad = current.get_key(i).get_byte(pass);
                radix_count[rad] += 1;
            }

            for(size_t i = 1; i < nbCount; ++i){
                radix_count[i] += radix_count[i-1];
            }

            for(size_t i = len(); i-->0 ;){
                rad = current.get_key(i).get_byte(pass);
                size_t place = --radix_count[rad];

                aux.get_value(place) = current.get_value(i);
                aux.get_key(place) = current.get_key(i);
            }
        }
        if(!b_current){
            *this = list_aux;
        }
    }

};


// ATTENTION BON POUR DENSITE <= 1
template <typename T_TYPE_KEY1, typename T_TYPE_VALUE1, size_t T_LEN1,typename T_TYPE_KEY2, typename T_TYPE_VALUE2, size_t T_LEN2, typename T_FUNCTION>
void merge_oneElementValue_list(list_SOA_fixed_length<T_TYPE_KEY1,T_TYPE_VALUE1,T_LEN1>& __restrict list1,list_SOA_fixed_length<T_TYPE_KEY2,T_TYPE_VALUE2,T_LEN2>&  __restrict list2, auto&& __restrict fun){
    size_t i1 = 0;
    size_t i2 = 0;
    size_t cursor1 = 0;
    size_t cursor2 = 0;
    size_t len_l1 = list1.len();
    size_t len_l2 = list2.len();
    
    while(i1 < len_l1 && i2 < len_l2){
        if(list1.get_key(i1) == list2.get_key(i2)){
            cursor2 = 1;
            cursor1 = 0;
            while((i2 + cursor2 < len_l2) && (list2.get_key(i2) == list2.get_key(i2+cursor2))){
                cursor2 += 1;
            }
            while((i1 + cursor1 < len_l1) && (list1.get_key(i1) == list1.get_key(i1+cursor1))){
                fun(list1.get_value(i1+cursor1),list2.get_value(i2),cursor2);
                cursor1 += 1;
            }
            i1 += cursor1;
            i2 += cursor2;
        }else if(list1.get_key(i1)<list2.get_key(i2)){
            i1 += 1;
        }else{
            i2 += 1;
        }
    }
}

// ATTENTION BON POUR DENSITE <= 1
template <typename T_TYPE_KEY1, typename T_TYPE_VALUE1, size_t T_LEN1,typename T_TYPE_KEY2, typename T_TYPE_VALUE2, size_t T_LEN2>
void naive_merge_list(list_SOA_fixed_length<T_TYPE_KEY1,T_TYPE_VALUE1,T_LEN1>& __restrict list1,list_SOA_fixed_length<T_TYPE_KEY2,T_TYPE_VALUE2,T_LEN2>&  __restrict list2, auto&& __restrict fun){
    size_t i1 = 0;
    size_t i2 = 0;
    size_t cursor1 = 0;
    size_t cursor2 = 0;
    size_t len_l1 = list1.len();
    size_t len_l2 = list2.len();
    
    while(i1 < len_l1 && i2 < len_l2){
        
        if(list1.get_key(i1) == list2.get_key(i2)){

            cursor2 = 1;
            cursor1 = 1;
            while((i2 + cursor2 < len_l2) && (list2.get_key(i2) == list2.get_key(i2+cursor2))){
                cursor2 += 1;
            }
            while((i1 + cursor1 < len_l1) && (list1.get_key(i1+cursor1) == list2.get_key(i2))){
                cursor1 += 1;
            }

            for(size_t c1 = 0; c1 < cursor1; ++c1){
                for(size_t c2 = 0; c2 < cursor2; ++c2){
                   // std::cout << "HH" << std::endl;
                /*    list1.get_key(i1+c1).print();
                    list2.get_key(i2+c2).print();
                    list1.get_value(i1+c1).print();
                    list2.get_value(i2+c2).print();*/
                    fun(list1.get_value(i1+c1),list2.get_value(i2+c2));
                }
            }
            i1 += cursor1;
            i2 += cursor2;
        }else if(list1.get_key(i1)<list2.get_key(i2)){
            i1 += 1;
        }else{
            i2 += 1;
        }
    }
}
/*
template <typename T_TYPE_KEY1, typename T_TYPE_VALUE1, size_t T_LEN1,typename T_TYPE_KEY2, typename T_TYPE_VALUE2, size_t T_LEN2, typename T_FUNCTION>
void merge_oneElementSimple_list(list_tri_SOA_average<T_TYPE_KEY1,T_TYPE_VALUE1,T_LEN1> * __restrict list1,list_tri_SOA_average<T_TYPE_KEY2,T_TYPE_VALUE2,T_LEN2> *  __restrict list2, T_FUNCTION&& __restrict fun){
    size_t i1 = 0;
    size_t i2 = 0;


    while(i1 < list1->get_length() && i2 < list2->get_length()){
        if(*get_key(list1,i1) == *get_key(list2,i2)){
            fun(get_value(list1,i1),get_value(list2,i2),1);
            i1 += 1;
            i2 += 1;
        }else if(*get_key(list1,i1) < *get_key(list2,i2)){
            i1 += 1;
        }else{
            i2 += 1;
        }
    }
}
*/
#endif
