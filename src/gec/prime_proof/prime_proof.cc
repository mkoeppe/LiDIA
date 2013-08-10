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
//  File    : prime_proof.cc
//  Author  : Jochen Hechler (JH)
//      Changes : See CVS log
//
//==============================================================================================


#include        <iostream>
#include 	<cstdio>
#include	<fstream>
#include	<cstdlib>
#include        "LiDIA/prime_proof.h"
#include        "LiDIA/bigint.h"
#include        "LiDIA/LiDIA.h"
#include        "LiDIA/base_vector.h"
#include	"LiDIA/certificate.h"
#include	"LiDIA/factorization.h"
#include	"LiDIA/single_factor.h"

#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif


    const bigfloat prime_proof::alpha = 0.015;
    const bigint prime_proof::default_delta = -21311;
    const bigint prime_proof::default_delta_4 = -21311;
    const bigint prime_proof::default_delta_1 = -125579;
    const bigint prime_proof::default_delta_2 = -53444;
    // as default_delta is odd, the corresponding case is set to 1
    const short prime_proof::default_delta_case = 1;
    // the default number of subintervals is 2^40
    const long prime_proof::bitlength_no_of_intervals = 40;
    // some constant static variables for the OEF case
    const lidia_size_t prime_proof::default_word_length = 32;
    const lidia_size_t prime_proof::default_degree_oef = 5;


//
// CONSTRUCTORS
//
    prime_proof::prime_proof (const bigint & p)
    {
	verbose = false;
	set_prime(p);
	top_level = true;
	ecpp_order_mode = 0;
	cl_min = 1;
	cl_max = 50;
	pollard_roh_rounds = 4;
	pollard_roh_start = 6;
	pl_length = 100;
    }

    prime_proof::prime_proof ()
    {
        verbose = false;
        top_level = true;
        ecpp_order_mode = 0;
        cl_min = 1;
        cl_max = 50;
        pollard_roh_rounds = 4;
        pollard_roh_start = 6;
        pl_length = 100;
    }

//
// DESTRUCTOR
//

    prime_proof ::~prime_proof()
    {
    }


//
// FUNCTIONS
//
	
    void prime_proof::set_prime(const bigint & p)
    {
	n = p;
    }



    bool prime_proof::prove_prime(const bigint & new_p)
    {
	set_prime(new_p);
	return prove_prime();
    }

    bool prime_proof::prove_prime()
    {
	if(verbose)std::cout<<" Primelength:"<<n.decimal_length()<<std::endl;
	if(!n.is_prime()) {
	    return(false);
	}
	if(n.decimal_length() < 15)
	{
	    if(verbose) {
		std::cout <<" Make the SPP Test"<<std::endl;
	    }
	    if(spp(n))
	    {
		down_cert = spp_get_certificate();
		return true;
	    }
	    return false;
	}
	if(verbose) {
	    std::cout <<" Test for the N-1 and N+1 test"<<std::endl;
	}
	if(factorize_n_pm_one()) {
	    return true;
	}
	if(verbose) {
	    std::cout <<" N-1 and N+1 can't be used"<<"\n";
	    std::cout <<" Make the ECPP Test"<<std::endl;
	}
	return ecpp_test();
    }

    // verifies a given certificate
    // returns if the certificate is valid

    bool prime_proof::verify_certificate(certificate & cert)
    {
	int nr = cert.get_number_of_tests();
	int i = 0;
	base_vector <bigint> bv;
	while(i<nr)
	{
	    bv = cert.get_cert_vector(i);
	    i++;
	    if(bv[0]==1)
	    {
		if(verbose)std::cout<<"VERIFY: n - plus - one"<<std::endl;
		if(npo_verify_proof_prime(bv))continue;
		else return false;
	    }
	    if(bv[0]==0)
	    {
		if(verbose)std::cout<<"VERIFY: spp"<<std::endl;
		if(spp_verify_proof_prime(bv))continue;
		else return false;
	    }
	    if(bv[0]==-1)
	    {
		if(verbose)std::cout<<"VERIFY: n - minus - one"<<std::endl;
		if(nmo_verify_proof_prime(bv))continue;
		else return false;
	    }
	
	    if(verbose)std::cout<<"VERIFY: ecpp"<<std::endl;
	    if(ecpp_verify_proof(bv))continue;
	    else return false;
	}
	return true;
    }
	
    // returns the certificate 
    certificate prime_proof::get_certificate()
    {
	return down_cert;
    }

    // if verbose is set true a lot of useful information are displayed during the proof

    void prime_proof::set_verbose(bool verbose_)
    {
	verbose = verbose_;
    }

    // calls the n-1 test as discriped in diploma thesis.
    // if the test is successfull a downstep is executed

    bool prime_proof::minus_one_test(bigint & biggest_prime_factor)
    {
	certificate c; 
	bool test = n_minus_one(  (n-1)/ biggest_prime_factor ,  biggest_prime_factor);
	if(test)down_cert = nmo_get_certificate();
	if(test)
	{
	    prime_proof *next_test;
	    next_test = new prime_proof();
	    (*next_test).set_verbose(verbose);	
	    (*next_test).set_top(false);
	    (*next_test).set_ecpp_mode(ecpp_order_mode);
	    (*next_test).set_prime(biggest_prime_factor);
	    test = (*next_test).prove_prime();
	    if(!test)
	    {
		delete next_test;
		return false;
	    }
	    c = (*next_test).get_certificate();
	    down_cert.add(c); 
	    delete next_test;
	    return true;
	}
	return false;
    }


    // calls the n+1 test as discriped in diploma thesis.
    // if the test is successfull a downstep is executed

    bool prime_proof::plus_one_test(bigint & biggest_prime_factor)
    {
	bool test = n_plus_one( (n+1)/ biggest_prime_factor , biggest_prime_factor);
	certificate c;
	if(test)down_cert = npo_get_certificate();	
	if(test)
	{
	    prime_proof *next_test;
	    next_test = new prime_proof;
	    (*next_test).set_verbose(verbose);
	    (*next_test).set_ecpp_mode(ecpp_order_mode);
	    (*next_test).set_top(false);
	    (*next_test).set_prime(biggest_prime_factor);
	    test = (*next_test).prove_prime();
	    if(!test)
	    {
		delete next_test;
		return false;
	    }
	    c = (*next_test).get_certificate();
	    down_cert.add(c);
	    delete next_test;
	    return true;
	}
	return false;
    }

    // if this object is the first in the donstep list, is_top should be true
    // if is_top is true, and a test will speak against n to be prime, the proof will stop

    void prime_proof::set_top(bool is_top)
    {
	top_level = is_top;
    }

    // calls the ecpp test
    // if the ecpp test has found a divisor of n, and n is the prime given in the 
    // top of the downstep list, the test stops and the certificate holds the divisor

    bool prime_proof::ecpp_test()
    {
	bool test = ecpp();
	if(test)down_cert  = ecpp_get_certificate();
	if(!test)
	{
	    if(top_level)
	    {
		//return the gcd as a divisor of no	
		if(verbose)std::cout<<"Divisor of n found:"<<found_divisor()<<std::endl;
		divisior_of_n = found_divisor();
		down_cert.delete_all();
		base_vector <bigint> bv;
		bv.set_capacity(1);
		bv[0]=divisior_of_n;
		down_cert.set_certificate(bv);
	    }
	    return false;
	}
	return true;

    }

    // sets the mode of the order choosing algorithm. For further infomartion see
    // diploma thesis chapter 6

    void prime_proof::set_ecpp_mode(int set_mode)
    {
        ecpp_order_mode = set_mode;
        if(ecpp_order_mode>2 || ecpp_order_mode<0)ecpp_order_mode=0;

    }

#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif

