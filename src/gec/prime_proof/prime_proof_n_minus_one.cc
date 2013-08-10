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
//  File    : prime_proof_n_minus_one.cc
//  Author  : Jochen Hechler (JH)
//      Changes : See CVS log
//
//==============================================================================================

#include	<iostream>
#include	"LiDIA/prime_proof.h"
#include	"LiDIA/certificate.h"
#include	"LiDIA/base_vector.h"
#include	"LiDIA/bigint.h"
#include	"LiDIA/LiDIA.h"


#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif


// This an implementation of the n-1 test as described in diploa thesis


//
// FUNCTIONS
//

    bool prime_proof::n_minus_one( const bigint & new_q, const bigint & new_p)
    {
	bigint u,k,l;
	u = 8;
	k = 9;
	l = 7;
	l = dgcd(u,k);
	q = new_q;
	p = new_p;
	if(verbose)std::cout<<"N-1 TEST n:"<<n<<" q:"<<q<<" p:"<<p<<std::endl;
        if(!n.is_prime())
        {
	    if(verbose)std::cout<<"n is not prime by bigint.is_prime()"<<std::endl;
	    return false;
        }
	if(!p.is_prime())
	{
	    if(verbose)std::cout<<"p is not prime by bigint.is_prime()"<<std::endl;
	    return false;
	}
	if(q>p)
	{
	    if(verbose)std::cout<<"p < q. Test does not work for this."<<std::endl;
	    return false;
	}
	bigint c; 
	bigint d;
	a = 2;
	int i = 0; // this is the counter for the tries.
	while(i < 100){
	    if(verbose)std::cout<<"a: "<<a<<std::endl;
	    power_mod(c,a,(n-1),n);
	    if(verbose)std::cout<<"a^n-1 mod n = "<<c<<std::endl;
	    if(c==1)
	    {
		power_mod(d,a,q,n);
		if(verbose)std::cout<<"a^q="<<d<<std::endl;
		d = dgcd( d-1 ,n);
		if(verbose)std::cout<<"gcd(a^q -1,n)="<<d<<std::endl;
	    }
	    if(c==1 & d==1)
	    {
		if(verbose)std::cout<<"N-1 TEST was successfull"<<std::endl;
		return true;
	    }
	    i++;
	    a = randomize(n-3);
	    a = a+3;	
	}
	if(verbose)std::cout<<"N-1 TEST was not successfull"<<std::endl;
	return false;
    }

    bool prime_proof::nmo_verify_proof_prime(base_vector <bigint> cert_vector)
    {
	n = cert_vector[1];
	p = cert_vector[2];
	q = cert_vector[3];
	a = cert_vector[4];
	if(verbose)std::cout<<"N-1 TEST VERIFY n:"<<n<<" q:"<<q<<" p:"<<p<<" a:"<<a<<std::endl;	
	bigint c;
	// c = a^(n-1) should be one modulo n;
	power_mod(c,a,(n-1),n);
	if(!(c==1))return false;
	// c = gcd(a^q -1,n) should be one;
	power_mod(c,a,q,n);
	c = dgcd( c-1 ,n); 
	if(!(c==1))return false;
	return true;
    }


    certificate  prime_proof::nmo_get_certificate()
    {
	base_vector<bigint> bv;
        bv.set_capacity(5);
        bv[0]=-1; // type of test
        bv[1]=n; // prime to test
        bv[2]=p; // p
        bv[3]=q; // q
	bv[4]=a; // a	 base of test
	certificate cert;
	cert.set_certificate(bv);
	return cert;
    }


#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif
