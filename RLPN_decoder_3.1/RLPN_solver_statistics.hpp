
#ifndef ISD_STAT_HEADER
#define ISD_STAT_HEADER

#include "math.hpp"
#include "linear_algebra.hpp"
#include "enumerate.hpp"
#include "permutation.hpp"
#include "sort.hpp"
#include "walsh.h"
#include <omp.h>

#include <algorithm>
#include <functional>
#include <array>
#include <random>
#include <vector>
#include <fstream>
#include <iterator>
#include <map>
#include "informationDecoding.hpp"
#include "birthday_paradox.hpp"

#define MAX_TRHEAD 10


template <size_t T_N,size_t T_K,size_t T_U,size_t T_P>
struct isd_statistical_decoding{
    static constexpr size_t K = T_K;
    static constexpr size_t N = T_N;
    
    static constexpr size_t P = T_P;
     static constexpr size_t U = T_U;
    static constexpr size_t D_U = T_U;
    static constexpr size_t L = T_K-T_U;
    static constexpr size_t T_L = T_K-T_U;
    using TYPE_CELL_VECTOR = uint64_t;
    using TYPE_CELL_VECTOR_BIRTHDAY = typename type_cell_vector<T_L>::type;


    using TYPE_VECTOR = t_vector<TYPE_CELL_VECTOR,T_K>;

    using TYPE_VECTOR_CODEWORD = t_vector<uint64_t,T_N>;
    using TYPE_VECTOR_BIRTHDAY = t_vector<TYPE_CELL_VECTOR_BIRTHDAY,L>;
    using TYPE_VECTOR_PRIME_1 = t_vector<uint64_t,N-K+L>;
    using TYPE_VECTOR_PRIME_2 = t_vector<uint64_t,T_K-T_L>;

    using TYPE_VECTOR_CODEWORD_PRIME_1 = t_vector<uint64_t,N-K+L>;
    using TYPE_VECTOR_CODEWORD_PRIME_2 = t_vector<uint64_t,K-L>;

    using TYPE_MATRIX_G = t_matrix_col<TYPE_CELL_VECTOR,T_K,T_N>;

    using TYPE_MATRIX_G_ROW = t_matrix_row<TYPE_CELL_VECTOR,T_K,T_N>;

    using TYPE_MATRIX_G_PRIME_1 =  t_matrix_col<TYPE_CELL_VECTOR_BIRTHDAY,L,N-K+L>;
    using TYPE_MATRIX_G_PRIME_2 =  t_matrix_col<TYPE_CELL_VECTOR,K-L,N-K+L>;
    
    using TYPE_VECTOR_E = t_vector< typename type_cell_vector<D_U>::type,D_U>;

    using TYPE_VECTOR_PRIME_2_EBAR = t_vector<uint64_t,T_K-T_L-T_U>;


    using TYPE_VECTOR_EBAR = t_vector<uint64_t,T_N-T_U>;

    using TYPE_VECTOR_EBAR_1 = t_vector<uint64_t,(T_N-T_U)/2>;
    using TYPE_VECTOR_EBAR_2 = t_vector<uint64_t,T_N-T_U - (T_N-T_U)/2>;
    size_t D_nbEquation;


    permutation_matrix<T_N> perm;
    size_t set_E[D_U];

    // GENERATOR MATRIX
    TYPE_MATRIX_G G_mat_col;
    TYPE_VECTOR_CODEWORD y;

    TYPE_VECTOR_CODEWORD e;
    TYPE_VECTOR_CODEWORD c;

    TYPE_VECTOR_CODEWORD_PRIME_1 y_permuted_prime_1;
    TYPE_VECTOR_CODEWORD_PRIME_2 y_permuted_prime_2;

    
    
    TYPE_MATRIX_G G_mat_col_permuted;
    

    TYPE_MATRIX_G_PRIME_1  G_mat_prime_1;

    TYPE_MATRIX_G_PRIME_2  G_mat_prime_2;

    isd_statistical_decoding(TYPE_MATRIX_G& _generator_matrix, TYPE_VECTOR_CODEWORD& _y,t_vector<uint64_t,PARAM_N>& static_e,t_vector<uint64_t,PARAM_N>& static_c){
        G_mat_col = _generator_matrix;
        y = _y;
        e = static_e;
        c = static_c;
    }

    // PARAM THE SET E
    //

    


    bool solve_naive(size_t E[D_U], auto&& callback, informationIteration& info){
        
        

            TYPE_MATRIX_G_ROW G_mat_row;
            TYPE_VECTOR_CODEWORD v_y(y);
            //size_t nb_eq_per_iteration = binomial((size_t) ((N-K+L)/2),(size_t) P/2) * binomial((size_t) ((N-K+L)/2),(size_t) P/2);
            //nb_eq_per_iteration /= pow(2,L);
            //size_t nb_IT = (nb_equation / nb_eq_per_iteration) + 1;
            //std::cout << "Equation per dumer iteration : " <<  nb_eq_per_iteration << std::endl;
            //std::cout << "Number of Dumer iteration : " <<  nb_IT << std::endl;
            TYPE_VECTOR_CODEWORD y_permuted;
            //for(size_t i = 0; i < nb_IT; ++i){
                //cout << "IT" << endl;
                //cout << nbEq_produced << endl;
                size_t nbSupress = 0;
                bool sucess = false;
                TYPE_MATRIX_G G_mat_col_permuted_before_gauss;
                perm.template rand_E<D_U>(E);
                while(!sucess){
                    //perm.template rand_E<D_U>(E);
                    //perm.print();
                    G_mat_col_permuted.permute(G_mat_col,perm);
                    G_mat_col_permuted_before_gauss.permute(G_mat_col,perm);
                    G_mat_row = G_mat_col_permuted.inverse_representation();
                    sucess = partial_echelon_form(G_mat_row, I_K - I_L,nbSupress);
                    if(nbSupress < D_U){
                        return false;
                    }
                    if(!sucess){
                        //cout << "NO SUCCESS" << endl;
                        //int h;
                        //cin >> h;
                    }
                }
            
            //  perm.print();
        

                y_permuted.permute(y,perm);
                TYPE_VECTOR_CODEWORD c_permuted;
                c_permuted.permute(c,perm);

            // cout << "PERMUTED Y - C : " << (y_permuted^c_permuted).hamming() << endl;
            //  y.print();
            //  perm.print();
            // y_permuted.print();

    
                TYPE_VECTOR_CODEWORD e_permuted;
                e_permuted.permute(e,perm);

                TYPE_VECTOR_EBAR e_permuted_EBAR;
                e_permuted.template extract_last<0>(e_permuted_EBAR);

             

                TYPE_VECTOR_E e_permuted_E;
                e_permuted.template extract_last<(T_N-T_U)>(e_permuted_E);

                size_t weight_U = e_permuted_E.hamming();

               
                size_t v_real = e_permuted_EBAR.hamming();
              
           
                size_t T_T = e_permuted.hamming();




                info.weight_on_U = weight_U;
                info.weight_on_N = T_T-weight_U;
                double vareps = varepsilon( T_N - T_U ,T_P,info.weight_on_N);
                double proba = (1-vareps)/2.0;
                info.theoric_error = proba;



                
            
                
                //info2 = info1;
                
                //y_permuted.print();
                //int b;
                //cin >> b;

                y_permuted.template extract_last<0>(y_permuted_prime_1);


                y_permuted.template extract_last<N-K+L>(y_permuted_prime_2);


                //THOMAS
                TYPE_VECTOR_PRIME_2_EBAR y_permuted_prime_2_EBAR;
                y_permuted.template extract_last<N-K-L>(y_permuted_prime_2_EBAR);

    

                G_mat_col_permuted = G_mat_row.inverse_representation();
                

                G_mat_col_permuted.template extract_last<I_K - I_L,0>(G_mat_prime_1);
                G_mat_col_permuted.template extract_last<I_K - I_L,I_L>(G_mat_prime_2);

                TYPE_MATRIX_G v_G_mat_col(G_mat_col);

                //G_mat_col_permuted.print();
                //cout << "G PRIME 1   : " << I_K - I_L << endl;
                //G_mat_prime_1.print();
                TYPE_VECTOR_BIRTHDAY S01;
                
                
                
                multimap<uint32_t, permutation<T_P>> good_elems; // empty multimap container
               
                // insert elements in random order
               
                TYPE_MATRIX_G_PRIME_2& v_G_mat_prime_2(G_mat_prime_2);
                TYPE_VECTOR_CODEWORD_PRIME_2& v_y_permuted_prime_2 = y_permuted_prime_2;
                TYPE_VECTOR_CODEWORD_PRIME_1& v_y_permuted_prime_1 = y_permuted_prime_1;
                permutation_matrix<T_N> perm_matrix(perm);

                auto store_parity = [&S01,&good_elems,&v_G_mat_prime_2](TYPE_VECTOR_BIRTHDAY& k,permutation<T_P>& permu){
                     if(k == S01){
                        TYPE_VECTOR_PRIME_2 sol;
                        for(size_t i = 0; i < T_P; ++i){
                            sol ^= v_G_mat_prime_2[permu[i]];
                        }
                        good_elems.insert(pair<uint32_t, permutation<T_P>>(sol.vect[0], permu));
                     }
                };
                TYPE_VECTOR_BIRTHDAY S02;
                enumerate2<T_P>(G_mat_prime_1,S02,0,T_N-T_U,store_parity);


                
                size_t *v_E = E;

                bool stop = false;
                auto itr = good_elems.begin();

                std::random_device dev;
                std::mt19937 rng(dev());

                auto sub_solution = [&v_G_mat_prime_2,&v_y_permuted_prime_2,&v_y_permuted_prime_1,&e_permuted_EBAR,&callback](permutation<T_P>& permu){
                        TYPE_VECTOR_PRIME_2 sol;
                        for(size_t i = 0; i < T_P; ++i){
                            sol ^= v_G_mat_prime_2[permu[i]];
                        }
                        bool part1 = sol.dot_product(v_y_permuted_prime_2);


                        TYPE_VECTOR_CODEWORD_PRIME_1 h_EBAR;
                        for(size_t i = 0; i < T_P; ++i){
                            h_EBAR.set(permu[i],1);
                        }
                        bool part2 = h_EBAR.dot_product(v_y_permuted_prime_1);

                        bool error = h_EBAR.dot_product(e_permuted_EBAR);

                        size_t value_E = sol.vect[0];

                        bool trueValue = error^part1^part2;
                        callback(value_E, part1^part2,trueValue,error);
                };

                for(auto itr = good_elems.begin(); itr != good_elems.end(); itr = good_elems.upper_bound(itr->first)) {
                    auto itbis = good_elems.equal_range(itr->first);
                    size_t ct = good_elems.count(itr->first);
                    std::uniform_int_distribution<std::mt19937::result_type> dist(0,ct-1);
                    size_t nb = dist(rng);

                    size_t counter = 0;
                    for(auto  kv = itbis.first; kv != itbis.second; ++kv){
                        if(counter == nb){
                            sub_solution(kv->second);
                            /*cout << "TOTAL : " << ct << endl;
                            cout << "CHOSEN : " << nb + 1 << endl;
                            cout << "---" << endl;*/
                        }
                        counter++;
                    }
                }

            return true;
    }


    bool solve_dumer(size_t E[D_U], auto&& callback, informationIteration& info){
        // TO BE DYNAMICALLY MODIFIED.....
        size_t nbIt = 40;
        
        // ATTENTION !!! SI SUPPRIMER CONTRAINTE ALORS METTRE VARPESILON PLUS COMPLIQUE PLUS TARD
        //assert((N-U) / 2 == N-U - (N-U) / 2 );
        //assert((T_P) / 2 == T_P - T_P/2);

        
        for(size_t it = 0; it < nbIt; ++it){

            TYPE_MATRIX_G_ROW G_mat_row;
            TYPE_VECTOR_CODEWORD v_y(y);

            TYPE_VECTOR_CODEWORD y_permuted;

                size_t nbSupress = 0;
                bool sucess = false;
                TYPE_MATRIX_G G_mat_col_permuted_before_gauss;
                perm.template rand_E<D_U>(E);
                while(!sucess){

                    G_mat_col_permuted.permute(G_mat_col,perm);
                    G_mat_col_permuted_before_gauss.permute(G_mat_col,perm);
                    G_mat_row = G_mat_col_permuted.inverse_representation();
                    sucess = partial_echelon_form(G_mat_row, I_K - I_L,nbSupress);
                    if(nbSupress < D_U){
                        return false;
                    }
                    if(!sucess){
                        //cout << "NO SUCCESS" << endl;
                        //int h;
                        //cin >> h;
                    }
                }
            
            //  perm.print();
        

                y_permuted.permute(y,perm);
                TYPE_VECTOR_CODEWORD c_permuted;
                c_permuted.permute(c,perm);

            // cout << "PERMUTED Y - C : " << (y_permuted^c_permuted).hamming() << endl;
            //  y.print();
            //  perm.print();
            // y_permuted.print();

    
                TYPE_VECTOR_CODEWORD e_permuted;
                e_permuted.permute(e,perm);

                TYPE_VECTOR_EBAR e_permuted_EBAR;
                e_permuted.template extract_last<0>(e_permuted_EBAR);

                TYPE_VECTOR_EBAR_2 e_permuted_EBAR_2;
                e_permuted.template extract_last<0>(e_permuted_EBAR_2);

                TYPE_VECTOR_EBAR_1 e_permuted_EBAR_1;
                e_permuted.template extract_last<(T_N-T_U)-(T_N-T_U)/2>(e_permuted_EBAR_1);
                //cout << (T_N-T_U)/2 << endl;
                //cout << (T_N-T_U)-(T_N-T_U)/2 << endl;
                //cout << "---" << endl;
                // TO VERIFY
                //e_permuted.print();
                //e_permuted_EBAR_1.print();
                TYPE_VECTOR_E e_permuted_E;
                e_permuted.template extract_last<(T_N-T_U)>(e_permuted_E);

                size_t weight_U = e_permuted_E.hamming();

                size_t T_V1 = e_permuted_EBAR_1.hamming();
                size_t v_real = e_permuted_EBAR.hamming();
                //assert((T_N - T_U) % 2 == 0);
                //assert(T_P %2 == 0);

                size_t T_T = e_permuted.hamming();
                size_t T_V2 = e_permuted_EBAR_2.hamming();
                assert(T_T = T_V1 + T_V2 + weight_U);


                if(it == 0){
                    info.weight_on_U = weight_U;
                    info.weight_on_N = T_T-weight_U;
                    double vareps = varepsilon( T_N - T_U ,T_P,info.weight_on_N);
                    double proba = (1-vareps)/2.0;
                    info.theoric_error = proba;
                   /* if(info.weight_on_U == 10){
                        cout << info.weight_on_N << endl;
                        cout << vareps << endl;
                        cout << info.theoric_error << endl;
                       cout << "H" << endl;
                        cout << kraw(T_N-T_U,T_P,info.weight_on_N)  << endl;
                        cout << binom(T_N-T_U,T_P) << endl;
                        cout << "H" << endl;
                    }*/
                    
                }

                
            
         

                y_permuted.template extract_last<0>(y_permuted_prime_1);


                y_permuted.template extract_last<N-K+L>(y_permuted_prime_2);


                //THOMAS
                TYPE_VECTOR_PRIME_2_EBAR y_permuted_prime_2_EBAR;
                y_permuted.template extract_last<N-K-L>(y_permuted_prime_2_EBAR);

    

                G_mat_col_permuted = G_mat_row.inverse_representation();
                

                G_mat_col_permuted.template extract_last<I_K - I_L,0>(G_mat_prime_1);
                G_mat_col_permuted.template extract_last<I_K - I_L,I_L>(G_mat_prime_2);

                TYPE_MATRIX_G v_G_mat_col(G_mat_col);

                //G_mat_col_permuted.print();
                //cout << "G PRIME 1   : " << I_K - I_L << endl;
                //G_mat_prime_1.print();


                TYPE_VECTOR_BIRTHDAY S0;
                birthday_paradox<T_L,N-K+L,T_P> sub_problem(G_mat_prime_1,S0);

                TYPE_MATRIX_G_PRIME_2& v_G_mat_prime_2(G_mat_prime_2);
                TYPE_VECTOR_CODEWORD_PRIME_2& v_y_permuted_prime_2 = y_permuted_prime_2;
                TYPE_VECTOR_CODEWORD_PRIME_1& v_y_permuted_prime_1 = y_permuted_prime_1;
                permutation_matrix<T_N> perm_matrix(perm);

                
                size_t *v_E = E;
                //G_mat_prime_1.print();
                
                auto sub_solution = [&v_G_mat_prime_2,&v_y_permuted_prime_2,&v_y_permuted_prime_1,&e_permuted_EBAR,&e_permuted_E,&callback,it](permutation<T_P>& permu){

                        TYPE_VECTOR_PRIME_2 sol;
                        for(size_t i = 0; i < T_P; ++i){
                            sol ^= v_G_mat_prime_2[permu[i]];
                        }
                        bool part1 = sol.dot_product(v_y_permuted_prime_2);


                        TYPE_VECTOR_CODEWORD_PRIME_1 h_EBAR;
                        for(size_t i = 0; i < T_P; ++i){
                            h_EBAR.set(permu[i],1);
                        }
                        bool part2 = h_EBAR.dot_product(v_y_permuted_prime_1);

                        bool error = h_EBAR.dot_product(e_permuted_EBAR);

                        size_t value_E = sol.vect[0];

                        bool trueValue = error^part1^part2;

                       
                        callback(value_E, part1^part2,trueValue,error);
                    

                };
                
            
            
                sub_problem.solve(sub_solution);
                
        }
        return true;
        
    }


    static constexpr size_t I_N = T_N;
    static constexpr size_t I_K = T_K;
    static constexpr size_t I_L = T_L;
    static constexpr size_t I_P = T_P;
    static constexpr size_t I_U = T_U;


};
struct debug{
    debug() : iteration(0), falseNegative(0), falsePositive(0), found(0) {}
    void print(){
        cout << "NB ITERATION : " << iteration << endl;
        cout << "FOUND : " << (size_t) found << endl;
        cout << "falseNegative : " << falseNegative << endl;
        cout << "falsePositive : " << falsePositive << endl;
    }
    size_t iteration;
    size_t falseNegative;
    size_t falsePositive;
    bool found;
};


enum stop_RLPN {SUCCESS, FAIL_FALSE_POSITIVE, FAIL_FALSE_NEGATIVE, FAIL_BAD_BET};
template <size_t T_N,size_t T_K,size_t T_U,size_t T_P>
struct stat_decoding{


    using TYPE_VECTOR_CODEWORD = t_vector<uint64_t,T_N>;


    using TYPE_MATRIX_G = t_matrix_col<uint64_t,T_K,T_N>;
       
    stat_decoding(TYPE_MATRIX_G& _generator_matrix, TYPE_VECTOR_CODEWORD& _y,t_vector<uint64_t,PARAM_N>& _static_e,t_vector<uint64_t,PARAM_N>& _static_c,size_t _t) : G_mat_col(_generator_matrix), y(_y), static_e(_static_e),static_c(_static_c), t(_t) {}

    
    ~stat_decoding(){
        
    }

    void gather_info_walsh(long int *f_walsh, size_t E_SOL, informationIteration& info){
        size_t first_idx;
        long int first_value = -1;
        
        long int walsh_value_solution = f_walsh[E_SOL];
        size_t position_solution = 0;
        size_t nbEqualWalshSolution = 0;
        for(size_t i = 0; i < pow(2,T_U); ++i){
            if(f_walsh[i] > first_value){
                first_idx = i;
                first_value = f_walsh[i];
            }
            if(f_walsh[i] > walsh_value_solution){
                position_solution += 1;
            }
            if(f_walsh[i] == walsh_value_solution){
                nbEqualWalshSolution += 1;
            }
        }
        assert(nbEqualWalshSolution >= 1);
        std::random_device dev;
        std::mt19937 rng(dev());
        std::uniform_int_distribution<std::mt19937::result_type> dist(1,nbEqualWalshSolution);
        size_t offset = dist(rng);
        position_solution += offset - 1;

        size_t second_idx;
        long int second_value = -1;
        

        for(size_t i = 0; i < pow(2,T_U); ++i){
            if(f_walsh[i] > second_value && i != first_idx){
                second_idx = i;
                second_value = f_walsh[i];
            }
        }

        info.position_solution = position_solution;
        info.walsh_value_solution = walsh_value_solution;
        info.first_value_walsh = first_value;
        info.second_value_walsh = second_value;

    }
    bool solve(size_t E[T_U], informationIteration& info1, informationIteration& info2){

        isd_statistical_decoding<T_N,T_K,T_U,T_P> isd(G_mat_col,y,static_e,static_c);


        uint8_t *f = new uint8_t[power(2,T_U)]();
		uint8_t *f_trueValue = new uint8_t[power(2,T_U)]();


        for(size_t i = 0; i < power(2,T_U); ++i){
            f[i] = 2;
            //f[i] = rand()%2;
        }

        size_t nbEq = 0;

        uint8_t *filling = new uint8_t[power(2,T_U)]();
		double moyErr = 0;



        auto callback = [f,&nbEq,filling,&f_trueValue,&moyErr](auto&& idx,bool res,bool trueValue, bool error){
            assert(idx < power(2,T_U));
            if(filling[idx] == 0){
                f[idx] = res;
                //moy += res;
				f_trueValue[idx] = trueValue;
				moyErr += error;
                
                //nbEq_produced +=1;
				nbEq +=1;
                //cout <<"H" << endl;
                filling[idx] = 1;
            }
        };


        //moy /= nbEq;


       
        
        bool ok = false;
        //ok = isd.solve2(E,callback,info1,nbItPerDumer);
       /* if(nbItPerDumer == 0){
            ok = isd.solve_naive(E,callback,info1);
        }else{
            ok = isd.solve_dumer(E,callback,info1,nbItPerDumer);
        }*/

         ok = isd.solve_dumer(E,callback,info1);
       
        
        if(!ok){

            return false;
            //break;
        }


		info1.number_distinct_parityChecks = nbEq;

		//info1.perfect_weight_distribution = ((info1.weight_on_N1 == PARAM_V/2) && (info1.weight_on_N2 == PARAM_V - PARAM_V/2));
		info2 = info1;

		info1.average_proba_error = ((double) moyErr) / ((double) nbEq);


         //using TYPE_VALUE = typename type_cell_vector<T_U>::type;
        using TYPE_F_WALSH = long int;
        using TYPE_VALUE = long int;
        TYPE_F_WALSH *f_walsh = new TYPE_F_WALSH[power(2,T_U)];
        walsh(f,f_walsh,T_U);
        


    
		
        size_t E_SOL = 0;
        for(size_t i = 0; i < T_U; ++i){
            if(static_e[E[i]]){
                E_SOL += power(2,T_U-i-1);
            }
        }

        gather_info_walsh(f_walsh,E_SOL,info1);
		

		std::random_device rd{}; // use to seed the rng 
		std::mt19937 rng{rd()}; // rng
		std::bernoulli_distribution bernou(info1.theoric_error);
		double total_Err = 0;
		 for(size_t i = 0; i < power(2,T_U); ++i){
            if(f[i] != 2){
				bool err = bernou(rng);
				f[i] = f_trueValue[i]^err;
				total_Err += err;
			}
            //f[i] = rand()%2;
        }

		info2.average_proba_error = total_Err / ((double) nbEq);

		walsh(f,f_walsh,T_U);
        
        
        gather_info_walsh(f_walsh,E_SOL,info2);

		
		delete[] f_trueValue;

        delete[] f;

        delete[] f_walsh;
  
        delete[] filling;

        return true;
    }
	void make_iteration(informationIteration& info1,informationIteration& info2){
		bool goodE = false;

		while(!goodE){
			size_t E[T_U];
			permutation_matrix<T_N> p;

			p.rand();
			for(size_t i = 0; i < T_U; ++i){
				E[i] = p[i];
			}

			goodE = solve(E,info1,info2);
		}

	}
	void decode_code_stat(size_t nbIt,informationDecode& dec1, informationDecode& dec2){
        dec1.iterations.resize(nbIt);
        dec2.iterations.resize(nbIt);
        
		for(int i = 0; i < nbIt; ++i){
			//dec1.iterations.push_back({});
			//dec2.iterations.push_back({});
            make_iteration(dec1.iterations[i],dec2.iterations[i]);
            //cout << i << endl;
            //int tid = omp_get_thread_num();
           // cout << tid << endl;
		}
		
	}
    size_t decode_code(size_t nbIt,size_t u,long int treshold){

        informationIteration info1;
        informationIteration info2;
        size_t i = 0;
        bool stop = false;
        //stop_RLPN stopR = FAIL_BAD_BET;
        size_t stopR = 3;
		while(i < nbIt && !stop){
			//dec1.iterations.push_back({});
			//dec2.iterations.push_back({});
            make_iteration(info1,info2);
            if(info1.weight_on_N == u){
                if(info1.position_solution == 0 && info1.walsh_value_solution >= treshold){
                    //stopR = SUCCESS;
                    stopR = 0;
                }else if(info1.position_solution != 0 && info1.walsh_value_solution >= treshold){
                    //stopR = FAIL_FALSE_POSITIVE;
                   stopR = 1;
                }
                else{
                   //stopR = FAIL_FALSE_NEGATIVE;
                    stopR = 2;
                }
            }else{
                if(info1.position_solution != 0 && info1.walsh_value_solution >= treshold){
                    //stopR = FAIL_FALSE_POSITIVE;
                    stopR = 1;
                }
            }
            //if(stop_RLPN == SUCCESS || stop_RLPN == FAIL_FALSE_POSITIVE){
            if(stopR == 0 || stopR== 1){
                stop = true;
            }

            //cout << i << endl;
            ++i;
		}
		return stopR;
	}


    // ONLY CALLABLE WHEN ERROR IS IN STANDARD FORM 111100 on E, 1111000000001111 on N
    // SLIDE -> 1 LESS FALSE POSITIVE
    double find_treshhold(double T_T, double T_V,size_t nbIt,double slide){
            double it = 0;
            double meanGoodmaxvalue = 0;
            for(size_t i = 0; i < nbIt; ++i){

                size_t E[T_U];
                for(size_t i = 0; i < T_U;++i){
                    E[i] = i;
                }

                size_t E_SOL = 0;
                size_t wt_e_E = 0;
                for(size_t i = 0; i < T_U; ++i){
                    if(static_e[E[i]]){
                        wt_e_E += 1;
                        E_SOL += power(2,T_U-i-1);
                    }
                }

   

                double maxvalue = 0;
                size_t index_maxvalue = 0;
                bool res = solve(E,maxvalue,index_maxvalue);

               
                if(wt_e_E == T_T - T_V){
                    //cout << maxvalue << endl;
                    if(index_maxvalue == E_SOL){
                        meanGoodmaxvalue += maxvalue;
                        it += 1;
                    }
                }
            }
            if(it == 0){
                exit(0);
            }

            double averageGoodmaxvalue = meanGoodmaxvalue / ((double) it);
            

            it = 0;

            double meanBadmaxvalue = 0;
            for(size_t i = 0; i < nbIt; ++i){
                permutation_matrix<T_N> p;
                p.rand();

                size_t E[T_U];

                for(size_t i = 0; i < T_U; ++i){
                    E[i] = p[i];
                }


                size_t E_SOL = 0;
                size_t wt_e_E = 0;
                for(size_t i = 0; i < T_U; ++i){
                    if(static_e[E[i]]){
                        wt_e_E += 1;
                        E_SOL += power(2,T_U-i-1);
                    }
                }

   
                double maxvalue = 0;
                size_t index_maxvalue = 0;
                bool res = solve(E,maxvalue,index_maxvalue);
               /* cout << index_maxvalue << endl;
                cout << E_SOL<< endl;
                cout << "F" << endl;*/
                if(index_maxvalue != E_SOL){
                    meanBadmaxvalue += maxvalue;
                    it += 1;
                }
                

            }

            if(it == 0){
                exit(0);
            }

             double averageBadmaxvalue = meanBadmaxvalue / ((double) it);

            cout << "Average good max : " << averageGoodmaxvalue << endl;
            cout << "Average bad max :"<< averageBadmaxvalue << endl;
            return averageBadmaxvalue + (averageGoodmaxvalue - averageBadmaxvalue)*slide;

           //}
    }


    TYPE_MATRIX_G G_mat_col;
    TYPE_VECTOR_CODEWORD y;
    t_vector<uint64_t,PARAM_N> static_e;
    t_vector<uint64_t,PARAM_N> static_c;
    size_t t;
};


	auto make_problem(size_t t){
		static constexpr int T_N = PARAM_N;
		static constexpr int T_K = PARAM_K;
		static constexpr int T_U = PARAM_U;

		static constexpr int T_P = PARAM_P;


		
		using TYPE_CELL_VECTOR = typename type_cell_vector<T_K>::type;
		srand(time(NULL));
		// WARNING PROBLEM TYPE_CELL_VECTOR !!!!!!!!!!!!
		t_matrix_col<uint64_t,T_K,T_N> M;
		t_matrix_row<uint64_t,T_K,T_N> M_row;

		t_vector<uint64_t,T_K> message;
		t_vector<uint64_t,T_N> codeword;
		t_vector<uint64_t,T_N> error_codeword;
		t_vector<uint64_t,T_N> error_base;
		t_vector<uint64_t,T_N> error;
		
		M.rand();
		M_row = M.inverse_representation();



		for(size_t i = 0;i < T_K; ++i){
			if(std::rand()%2){
				codeword ^= M_row[i];
			}
		}

		/*for(size_t i = 0; i < T_T; ++i){
			error_base.set(i,1);
		}*/
		for(size_t i = 0; i < t; ++i){
			error_base.set(i,1);
		}

		error_codeword = codeword ^ error_base;
		t_vector<uint64_t,PARAM_N> static_e = error_base;
		t_vector<uint64_t,PARAM_N> static_c = codeword;

		stat_decoding<T_N,T_K,T_U,T_P> stat_dec(M,error_codeword,static_e,static_c,t);
		return stat_dec;
	}

    

	void make_stat(size_t nbIt, size_t nbCode, informationAverage& avg1, informationAverage& avg2){
		size_t T_T = PARAM_T;
		//size_t T_V = PARAM_V;

		informationAverage av1;
		informationAverage av2;
        avg1.codes.resize(nbCode);
        avg2.codes.resize(nbCode);
        #pragma omp parallel for num_threads(MAX_TRHEAD) shared(avg1,avg2)
		for(int i = 0; i < nbCode;++i){
            
			auto problem = make_problem(T_T);
			problem.decode_code_stat(nbIt,avg1.codes[i],avg2.codes[i]);
		}
		string st = avg1.value_print_python();
		//cout << st << endl;
        string F1 = "L_Parity=" +st;
        string st2 = avg2.value_print_python();
        string F2 = "L_TrueRandom="+st2;
        ofstream myfile;
        
        string header = to_string(PARAM_N) + "_" + to_string(PARAM_K) + "_" +to_string(PARAM_T) + "_" +to_string(PARAM_U) + "_" +to_string(PARAM_P)+".py";
        /*if(nbItDumer == 0){
            header += "_naiveEnumeration";
        }else{
            header += "_Dumer";
        }*/
        myfile.open(header);
        myfile << "[n,k,t,s,w] = " + header << endl;
        myfile << "nbIterationPerCode=" + to_string(nbIt) << endl;
        myfile << "nbRandomCodes=" + to_string(nbCode) << endl;

        myfile << F1 << endl;
        myfile << F2 << endl;
        myfile.close();
		//cout << st2 << endl;
	}

#endif