// -*- C++ -*-
//==============================================================================================
//
//      This file is part of LiDIA --- a library for computational number theory
//
//      Copyright (c) 1994--2003 the LiDIA Group.  All rights reserved.
//
//      See http://www.informatik.tu-darmstadt.de/TI/LiDIA/
//
//----------------------------------------------------------------------------------------------
//
//
//  File    : prime_proof_factorize.cc
//  Author  : Jochen Hechler (JH)
//      Changes : See CVS log
//
//==============================================================================================

#include        <iostream>
#include        <cstdio>
#include        <fstream>
#include        <cstdlib>
#include        <cassert>
#include        "LiDIA/prime_proof.h"
#include        "LiDIA/bigint.h"
#include        "LiDIA/LiDIA.h"
#include        "LiDIA/base_vector.h"
#include	"LiDIA/factorization.h"
#include	"LiDIA/single_factor.h"

#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif


// this function is used to test whether the n-1 or the n+1 test can be applied to n

    bool prime_proof::factorize_n_pm_one()
    {
        static char const env_var_name[] = "LIDIA_PRIME_DIFFS";

        char const* env_lidia_primes = getenv(env_var_name);
        std::string in_file;
	if(env_lidia_primes) {
	    in_file = env_lidia_primes;
	}
	if(in_file.empty()) {
	    in_file = LIDIA_PRIME_DIFFS;
	}

        std::ifstream in;
	in.open(in_file.c_str());
	if(!in.good()) {
	    std::string msg("Can't open file ");
	    msg += in_file;
	    lidia_error_handler("prime_proof::factorize_n_pm_one()",
				msg.c_str());
	}
	    
	bool plus_test = false;
	bool minus_test = false;
	bigint n_plus_one_tmp = n+1;
	bigint n_plus_one = n+1;
	bigint n_minus_one = n-1;
	while(remainder(n_minus_one,2) == 0) {
	    n_minus_one = n_minus_one/2;
	}
	bigint biggest_prime_factor_plus_one = 2;
	bigint biggest_prime_factor_minus_one = 2;

	int prime=0;
	int i = 0;
	bigint ii;
	int append;
	while(in.good())
	{
	    in >> append;
	    if(in.fail() && in.eof()) {
		break;
	    }
	    else if(in.fail()) {
		std::string msg("Error while reading from file ");
		msg += in_file;
		lidia_error_handler("prime_proof::factorize_n_pm_one()",
				    msg.c_str());
	    }

	    prime += append;
	    assert(prime > 1);
	    if(prime > (n-1)) {
		break;
	    }
	    ii = remainder(n_plus_one_tmp,prime);
	    if(ii==0)
	    {
		while(remainder(n_plus_one, prime) == 0) {
		    n_plus_one = n_plus_one/prime;
		}
		biggest_prime_factor_plus_one = prime;
	    }
	    else if(ii == 2)
	    {
		while(remainder(n_minus_one,prime)==0) {
		    n_minus_one = n_minus_one/prime;
		}
		biggest_prime_factor_minus_one = prime;
	    }
	}
	in.clear();
	in.close();
	if(!in.good()) {
	    std::string msg("Error while closing file ");
	    msg += in_file;
	    lidia_error_handler("prime_proof::factorize_n_pm_one()",
				msg.c_str());
	}

	// try to factor n-1
	single_factor <bigint> prime_compo;
	bigint large_prime;
	int expo = 0;
	int no_prime_compo = 0;
	i = 0;
	ii = 0;
	int a = pollard_roh_start;
	int b = pollard_roh_rounds + a;

	if(n_plus_one.is_prime()) {
	    biggest_prime_factor_plus_one=n_plus_one;
	}
	if(biggest_prime_factor_plus_one >
	   (n+1)/biggest_prime_factor_plus_one) {
	    // make the n+1 test
	    if(verbose) {
		std::cout <<" Make the N+1 Test"<<std::endl;
	    }
	    plus_test = true;
	    if(plus_one_test(biggest_prime_factor_plus_one)) {
		return true;
	    }
	}

	if(n_minus_one.is_prime()) {
	    biggest_prime_factor_minus_one=n_minus_one;
	}
        if(biggest_prime_factor_minus_one >
	   (n-1)/biggest_prime_factor_minus_one) {
	    // make the n-1 test
	    if(verbose) {
		std::cout <<" Make the N-1 Test"<<std::endl;
	    }
	    minus_test = true;
	    if(minus_one_test(biggest_prime_factor_minus_one)) {
		return true;
	    }
        }

	factorization <bigint> f;
	while(!n_plus_one.is_prime() &&
	      (n+1)/n_plus_one < n_plus_one &&
	      a < b &&
	      a*2 < n_plus_one.decimal_length() &&
	      !plus_test) {
	    f = PollardRho(n_plus_one,a);
	    no_prime_compo = f.no_of_prime_components();	
	    while(no_prime_compo > i) {
		expo=f.prime_exponent(i);
		prime_compo = f.prime_base(i);	
		large_prime = prime_compo.base(); 
		biggest_prime_factor_plus_one = large_prime;
		while( expo > ii ) {
		    n_plus_one = n_plus_one / large_prime;
		    ii++;
		}
		ii = 0;
		i++;
	    }
	    i = 0;
	    a++;
	}
        if(biggest_prime_factor_plus_one >
	   (n+1)/biggest_prime_factor_plus_one && !plus_test) {
	    // make the n+1 test
	    if(verbose) {
		std::cout <<" Make the N+1 Test"<<std::endl;
	    }
	    if(plus_one_test(biggest_prime_factor_plus_one)) {
		return true;
	    }
        }

	i=0;
	a=pollard_roh_start;
	ii=0;
        while(!n_minus_one.is_prime() &&
	      (n-1)/n_minus_one < n_minus_one &&
	      a < b &&
	      a*2 < n_minus_one.decimal_length() &&
	      !minus_test) {
	    f = PollardRho(n_minus_one,a);
	    no_prime_compo = f.no_of_prime_components();
	    while(no_prime_compo > i) {
		expo=f.prime_exponent(i);
		prime_compo = f.prime_base(i);
		large_prime = prime_compo.base();
		biggest_prime_factor_minus_one = large_prime;
		while( expo > ii ) {
		    n_minus_one = n_minus_one / large_prime;
		    ii++;
		}
		ii = 0;
		i++;
	    }
	    i = 0;
	    a++;
        }
        if(biggest_prime_factor_minus_one >
	   (n-1)/biggest_prime_factor_minus_one && !minus_test) {
	    // make the n-1 test
	    if(verbose) {
		std::cout <<" Make the N-1 Test"<<std::endl;
	    }
	    if(minus_one_test(biggest_prime_factor_minus_one)) {
		return true;
	    }
        }
	return false;
    }



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif
