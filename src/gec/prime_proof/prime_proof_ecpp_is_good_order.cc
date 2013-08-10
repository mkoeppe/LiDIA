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
//  File    : prime_proof_ecpp_is_good_order.cc
//  Author  : Jochen Hechler (JH)
//      Changes : See CVS log
//
//==============================================================================================

#include	<iostream>
#include        "LiDIA/path.h"
#include	"LiDIA/prime_proof.h"
#include	"LiDIA/bigint.h"
#include	"LiDIA/bigfloat.h"
#include	"LiDIA/base_vector.h"
#include	"LiDIA/factorization.h"
#include	"LiDIA/LiDIA.h"
#include	"LiDIA/error.h"



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif


// different modi are as follow
// ecpp_order_mode 0: first try prime_list, then try factorizing (default)
// ecpp_order_mode 1: try only prime_list, if there is no good_discriminant compute the next primelist and start ove
// ecpp_order_mode 2: try only factorization 
// the downstep will be activated

    bool prime_proof::is_good_order(bigint & t, bigint & y)
    {
	// m_{\pm} = n + 1 \pm t;
	// try for both m if there is a divisor in the primelist
	bigint m_plus_t = n + 1 + t;
	bigint m_minus_t = n + 1 - t;

	// the following orders are only for D=-3 and D=-4
	bigint m_3,m_4,m_5,m_6;
	// D = -4
	if(D==-4){
	    m_3 = n + 1 + 2*y;
	    m_4 = n + 1 - 2*y;
	}
	if(D==-3){
	    m_3 = n + 1 + t + 3*y;
	    m_4 = n + 1 - t - 3*y;
	    m_5 = n + 1 + t - 3*y;
	    m_6 = n + 1 - t + 3*y;
	}
	int i = 0;
	int jj;
	if(ecpp_order_mode==0 || ecpp_order_mode==1){
	    prime_factor = start_list;
	    bigint rest_1;
	    bigint rest_2;
	    bigint rest_3;
	    bigint rest_4;
	    bigint rest_5;
	    bigint rest_6;
	    while(i < pl_length)
	    {
		remainder(rest_1,m_plus_t,prime_factor);
		remainder(rest_2,m_minus_t,prime_factor);
		if(D==-4){
		    remainder(rest_3,m_3,prime_factor);
		    remainder(rest_4,m_4,prime_factor);
		    if(rest_3 == 0 || rest_4 == 0)
		    {
			if(rest_3 == 0)order = m_3;
			if(rest_4 == 0)order = m_4;
			if(prove_downstep(prime_factor))return prove_n_with_curve();
		    }
				
		}
		if(D==-3){
		    remainder(rest_3,m_3,prime_factor);
		    remainder(rest_4,m_4,prime_factor);
		    remainder(rest_5,m_5,prime_factor);
		    remainder(rest_6,m_6,prime_factor); 
		    if(rest_3 == 0 || rest_4 == 0 || rest_5 == 0 || rest_6 == 0)
		    {
			if(rest_3 == 0)order = m_3;
			if(rest_4 == 0)order = m_4;
			if(rest_5 == 0)order = m_5;
			if(rest_6 == 0)order = m_6;
			if(prove_downstep(prime_factor))return prove_n_with_curve();
		    }
		}
		if(rest_1 == 0 || rest_2 == 0)
		{
		    if(rest_1==0)order = m_plus_t;	
		    if(rest_2==0)order = m_minus_t;
		    if(prove_downstep(prime_factor))return prove_n_with_curve();
		}
		prime_factor = prime_factor + prime_list[i];
		i++;		
		
	    }
	}
	if(ecpp_order_mode==0 || ecpp_order_mode==2 ||  mode_one_tried){
	    if(verbose) {
		std::cout<<"ecpp: Begin factorization"<<std::endl;
	    }

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
		lidia_error_handler("prime_proof::is_good_order()",
				    msg.c_str());
	    }
	    
	    bigint n_plus_one = n + 1;
	    bigint rem_t,rem;

	    bigint prime = 0;
	    int append;
	    // first divide the possible orders by the stored primes
	    while(in.good())
	    {
		in >> append;
		if(in.fail() && in.eof()) {
		    break;
		}
		else if(in.fail()) {
		    std::string msg("Error while reading from file ");
		    msg += in_file;
		    lidia_error_handler("prime_proof::is_good_order()",
					msg.c_str());
		}

		add(prime,prime,append);
		if(prime > (m_minus_t)) {
		    break;
		}
		remainder(rem,n_plus_one,prime);
		remainder(rem_t,t,prime);
		if(rem==rem_t) {
		    while((m_minus_t%prime)==0) {
			m_minus_t=m_minus_t/prime;
		    }
		}
		rem.negate();
		add(rem,rem,prime);
		if(rem==rem_t) {
		    while((m_plus_t%prime)==0) {
			m_plus_t=m_plus_t/prime;
		    }
		}
		if(D==-4 || D==-3)
		{
		    while((m_3%prime)==0) {
			m_3=m_3/prime;
		    }
		    while((m_4%prime)==0) {
			m_4=m_4/prime;
		    }
		    if(D==-3)
		    {
			while((m_5%prime)==0) {
			    m_5=m_5/prime;
			}
			while((m_6%prime)==0) {
			    m_6=m_6/prime;
			}
		    }
		}
	    }
	    in.clear();
	    in.close();
	    if(!in.good()) {
		std::string msg("Error while closing file ");
		msg += in_file;
		lidia_error_handler("prime_proof::is_good_order()",
				    msg.c_str());
	    }

	    if(!m_minus_t.is_prime()) {
		m_minus_t = pollard_rho(m_minus_t);
	    }
	    if(!m_plus_t.is_prime()) {
		m_plus_t = pollard_rho(m_plus_t);
	    }
	    if(D==-4 || D==-3)
	    {
		if(!m_3.is_prime()) {
		    m_3 = pollard_rho(m_3);
		}
		if(!m_4.is_prime()) {
		    m_4 = pollard_rho(m_4);
		}
		if(D==-3)
		{
		    if(!m_5.is_prime()) {
			m_5 = pollard_rho(m_5);
		    }
		    if(!m_6.is_prime()) {
			m_6 = pollard_rho(m_6);
		    }
		}
	    }
	    if(m_plus_t.is_prime() && m_plus_t > s)
	    {
		order = n + 1 + t;
		prime_factor = m_plus_t;
		if(order!=prime_factor && prove_downstep(prime_factor)) 
		{
		    return prove_n_with_curve();
		}
			
			
	    }
	    if(m_minus_t.is_prime() && m_minus_t > s)
	    {
		order = n + 1 - t;
		prime_factor = m_minus_t;
		if(order!=prime_factor && prove_downstep(prime_factor))
		{
		    return prove_n_with_curve();
		}
	    }
	    if(D==-3)
	    {
		if(m_3.is_prime() && m_3  > s)
		{
		    order = n + 1 + t + 3*y;
		    prime_factor = m_3;
		    if(order!=prime_factor && prove_downstep(prime_factor))
		    {
			return prove_n_with_curve();
		    }
		}
		if(m_4.is_prime() && m_4  > s)
		{
		    order = n + 1 -t -3*y;
		    prime_factor = m_4;
		    if(order!=prime_factor && prove_downstep(prime_factor))
		    {
			return prove_n_with_curve();
		    }
		} 
		if(m_5.is_prime() && m_5  > s)
		{
		    order = n + 1 + t - 3*y;
		    prime_factor = m_5;
		    if(order!=prime_factor && prove_downstep(prime_factor)) 
		    {
			return prove_n_with_curve();
		    }
		}
		if(m_6.is_prime() && m_6  > s)
		{
		    order = n + 1 -t + 3*y;
		    prime_factor = m_6;
		    if(order!=prime_factor && prove_downstep(prime_factor))
		    {
			return prove_n_with_curve();
		    }
		}
	    }
	    if(D==-4)
	    {
		if(m_3.is_prime() && m_3  > s)
		{
		    order = n + 1 + 2*y;
		    prime_factor = m_3;
		    if(order!=prime_factor && prove_downstep(prime_factor))
		    {
			return prove_n_with_curve();
		    }
		}
		if(m_4.is_prime() && m_4  > s)
		{
		    order = n + 1 - 2*y;
		    prime_factor = m_4;
		    if(order!=prime_factor && prove_downstep(prime_factor))
		    {
			return prove_n_with_curve();
		    }
		}
	    }
	    return false;
	}	
	return false;
    }

    bigint prime_proof::pollard_rho(bigint & t)
    {
	single_factor <bigint> prime_compo;
	int expo = 0;
	int no_prime_compo = 0;
	int ii = 0;
	int i = 0;
	int l = t.decimal_length();
	l=l/2;
	if(l>9) {
	    l=16;
	}
	factorization <bigint> f;
	f = PollardRho(t,9);
	no_prime_compo = f.no_of_prime_components();
	while(no_prime_compo > i)
	{
	    expo=f.prime_exponent(i);
	    prime_compo = f.prime_base(i);
	    prime_factor = prime_compo.base();
	    while( expo > ii )
	    {
		t = t / prime_factor;
		ii++;
	    }
	    ii=0;
	    i++;
	}
	return t; 
    }

#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif
