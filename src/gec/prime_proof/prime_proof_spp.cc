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
//  File    : prime_proof_spp.cc
//  Author  : Jochen Hechler (JH)
//      Changes : See CVS log
//
//==============================================================================================

#include	<iostream>
#include 	<cstdio>
#include        <fstream>
#include        <cstdlib>
#include	"LiDIA/bigint.h"
#include	"LiDIA/LiDIA.h"
#include	"LiDIA/base_vector.h"
#include	"LiDIA/certificate.h"
#include	"LiDIA/prime_proof.h"


#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif


//
// CONSTRUCTORS
//


//this is an implementation of the spp test as described in diploma thesis

    bool prime_proof::spp(const bigint & n)
    {
	if(n<17)
	{
	    if(verbose && n.is_prime())std::cout <<"SPP: Test was successful "<<std::endl;	
	    return n.is_prime();
	}
	if(!n.is_prime())
	{
	    if(verbose && top_level)std::cout <<"SPP: Test wasnot successful, n is not prime "<<std::endl;
	    return false;
	}
	bigint limit = 341;
	limit = limit * 100000;
	limit = limit + 55007;
	limit = limit * 100000;
	limit = limit + 17283;
	limit = limit * 100;
	limit = limit + 21;
	
	bigint exp=0;
	bigint b;
	bigint a = 0;
	int tries = 0;
	int test = 0;
	r_i.set_capacity(12);
	r_i[0]=0; // type of test
	r_i[1]=n; // number to be tested
	if(n>=limit)return false;
  	// n-1 = 2^l * q
	bigint n_minus_one = n-1;
	while(n_minus_one.is_even())
	{
	    n_minus_one.divide_by_2();
	    exp++;
	}	 
	bigint q = n_minus_one;
	r_i[2]=exp;
	r_i[3]=q;
	if(verbose)std::cout <<"SPP: n-1 = 2^l *q, l: "<<exp<<" q: "<<q<<std::endl;
	bool tested =false;
	bool base_test = false;
	while(!tested)
	{
	    // b = a^q mod n
	    base_test = false;
	    if(a==13)a=17;
	    if(a==11)a=13;
	    if(a==7)a=11;
	    if(a==5)a=7;
	    if(a==3)a=5;
	    if(a==2)a=3;
	    if(a==0)a=2;
	    power_mod(b,a,q,n);
	    if(b==1 || b==-1 || b==(n-1))base_test=true;
	    r_i[test+4]=0;
	    tries = 0;
	    while( (exp > tries) && !base_test )		
	    {
		power_mod(b,b,2,n);
		tries++;
		if(b==-1 || b==(n-1))base_test=true;
		if(b==-1 || b==(n-1))r_i[test+4]=tries;
	    }
	    if(!base_test)
	    {
		if(verbose)std::cout <<"SPP: Test was not successful "<<std::endl;
		return false;	
	    }
	    if(a==17 & base_test)tested=true;
	    test++;
	}
	if(verbose)std::cout <<"SPP: Test was successful "<<std::endl;
	return true;
    }

    bool prime_proof::spp_verify_proof_prime(const base_vector <bigint> cert_vector)
    {
	
	n = cert_vector[1];
        bigint exp=cert_vector[2];
        bigint q = cert_vector[3];
	bigint b,a,tries;
	int test = 0;
	bigint xx;
	bool base_test = false;
	bool tested = false;
        while(!tested)
        {
	    // b = a^q mod n
	    base_test = false;
	    if(a==13)a=17;
	    if(a==11)a=13;
	    if(a==7)a=11;
	    if(a==5)a=7;
	    if(a==3)a=5;
	    if(a==2)a=3;
	    if(a==0)a=2;
	    power_mod(b,a,q,n,0);
	    if(b==1 || b==-1 || b==(n-1))base_test=true;
	    tries = 0;
	    xx = cert_vector[test+4];
	    while( (xx> tries) && !base_test )
	    {
		power_mod(b,b,2,n);
		tries++;
		if(b==-1 || b==(n-1))base_test=true;
	    }
	    if(!base_test)return false;
	    if(a==17 && base_test)tested=true;
	    test++;
        }
        return true;


    }


    certificate prime_proof::spp_get_certificate()
    {
	certificate c;
	c.set_certificate(r_i);
	return c;
    }


#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif
