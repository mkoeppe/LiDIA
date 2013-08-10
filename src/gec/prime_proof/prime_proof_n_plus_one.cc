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
//  File    : prime_proof_n_plus_one.cc
//  Author  : Jochen Hechler (JH)
//      Changes : See CVS log
//
//==============================================================================================

#include	<iostream>
#include        "LiDIA/certificate.h"
#include	"LiDIA/prime_proof.h"
#include	"LiDIA/bigint.h"
#include	"LiDIA/LiDIA.h"
#include	"LiDIA/bigmod.h"
#include	"LiDIA/base_vector.h"


#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif


// This is an implementation of the n+1 test as described in diploma thesis

  
// makes the test for $n=q*p$ 
// after 50 tests it will return false
    bool prime_proof::n_plus_one( const bigint & new_q, const bigint & new_p)
    {
	q = new_q;
	p = new_p;
	if(verbose)std::cout<<"N+1 TEST n:"<<n<<" q:"<<q<<" p:"<<p<<std::endl;
	q.divide_by_2();
	bigint g;
	g=n+1;
	g.divide_by_2();
	npo_r;
	npo_s=3;
	int j;
	int i = 0; // this is the counter for the tries.
	while(i < 50){
	    i++;
	    npo_s = randomize(n);
	    npo_s++;
	    npo_r = 1;
	    if(npo_s.is_odd())npo_r=2;
	    j = jacobi(npo_r-4*npo_s,n);
	    if(!(j==-1))continue;
	    if(lucas_sequence(npo_r,npo_s,g)!=0)continue;
	    if(lucas_sequence(npo_r,npo_s,q)!=0)continue;			
	    if(verbose)std::cout<<"N+1 TEST was successful "<<std::endl;
	    return true;
	}
	if(verbose)std::cout<<"N+1 TEST was not successful "<<std::endl;
	return false;
    }
 
  

    bool prime_proof::npo_verify_proof_prime(base_vector <bigint> cert_vector)
    {
	n = cert_vector[1];
	p = cert_vector[2];
	q = cert_vector[3];
	npo_s = cert_vector[4];	
	if(verbose)std::cout<<"N+1 TEST VERIFY n:"<<n<<" q:"<<q<<" p:"<<p<<" s:"<<npo_s<<std::endl;	
	bigint r;
	q.divide_by_2();
	r=1;
	if(npo_s.is_odd())r=2;
	if(jacobi((r-4*npo_s),n)!=-1)return false;
	bigint g;
	g = n+1;
	g.divide_by_2();
	bigint rem;
	if(lucas_sequence(r,npo_s,g)!=0)return false;
	if(lucas_sequence(r,npo_s,q)!=0)return false;
	if(verbose)std::cout<<"N+1 TEST was successful "<<std::endl;
	return true;
    }


    certificate prime_proof::npo_get_certificate()
    {
	base_vector<bigint> bv;
	bv.set_capacity(5);
	bv[0]=1; // type of test
	bv[1]=n; // prime to test
	bv[2]=p; // p
	bv[3]=2*q; // q
	bv[4]=npo_s; // s
	certificate cert;
	cert.set_certificate(bv);
	return cert;
    }


    bigmod prime_proof::lucas_sequence(const bigint & r, const bigint & s, const bigint & d)
    {
	// s even => r = 1;
	// s odd => r = 2;
	
	bigmod::set_modulus(n);

        bigint P=r;
        bigint Q=s;
	bigmod power_Q;
	bigmod U_m=1;
	bigmod U_h;
	bigmod V_m=P;
	bigmod V_h;
	bigmod D = P*P - 4*Q;
	bigmod t;
	bigint s_j = 1;

	
	int length = d.bit_length();	
	int i = length - 2;
	while(i>=0)
	{
	    square(V_h,V_m);
		
	    power(power_Q,Q,s_j);
	    multiply(power_Q,2,power_Q);
	    multiply(U_m,U_m,V_m);
	    multiply(s_j,2,s_j);
	    subtract(V_h,V_h,power_Q);

	    if(d.bit(i))
	    {
		multiply(V_m,P,V_h);
		multiply(t,D,U_m);
		add(V_m,V_m,t);
		V_m.divide_by_2();

		multiply(U_m,P,U_m);
		add(U_m,V_h,U_m);
		U_m.divide_by_2();
		V_h = V_m;
		add(s_j,s_j,1);
	    }
	    V_m = V_h;
	    i--;
	}
	return U_m; 

    }



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif
