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
//  File    : prime_proof_ecpp.cc
//  Author  : Jochen Hechler (JH)
//      Changes : See CVS log
//
//==============================================================================================

#include	<iostream>
#include	"LiDIA/prime_proof.h"
#include	"LiDIA/bigint.h"
#include	"LiDIA/bigfloat.h"
#include	"LiDIA/base_vector.h"
#include	"LiDIA/factorization.h"
#include	"LiDIA/LiDIA.h"
#include	"LiDIA/error.h"
#include	"LiDIA/quadratic_order.h"
#include	"LiDIA/certificate.h"


#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif



//
// FUNCTIONS
//

    bool prime_proof::ecpp()
    {
	cert_vector.set_capacity(8);
	if(verbose)std::cout<<"ecpp: Test for n:"<<n<<std::endl;
	mode_one_tried= false;
        calculate_lower_bound_s();
        if(ecpp_order_mode!=2)calculate_prime_list();
	if(verbose)std::cout<<"ecpp: Searching for Discriminant"<<std::endl;
	return find_next_discriminant(cl_min,cl_max);
    }
  
    void prime_proof::set_classnumbers(const int min,const int max)
    {
	cl_min = min;
	cl_max = max;
    }

    void prime_proof::set_pollard_rho_parameter(const int start_length,const int rounds)
    {
	pollard_roh_start = start_length;
	pollard_roh_rounds = rounds;
    }


	
    // if during the point operation on a curve over F_n a divisor of n is found
    // the following function will return the divisor

    bigint prime_proof::found_divisor()
    {
	if(gcd>1)return gcd;
	return false;
    }
	
    // this function starts the generating algorithm for the elliptic curves.
    // there are two cases where we alreadzy know the j-invariant of the curve
    // and we know the possible curves. so for these two cases we skip the following steps 
    // and go straight to the curves

    bool prime_proof::prove_n_with_curve ()
    {
	delta_case = 1;
	bigint abs_D = -D;
	if( (abs_D.bit(0) == 0 && abs_D.bit(1) == 0) )delta_case = 0;
	F_q    = galois_field( abs_D );	
	return generate_curve();
    }

    // calls the downstep for the ecpp test
    // note: the downstep in the ecpp test will be executed before the ecpp test is completed
    // moreover the downstep is executed as soon as a useable order of an elliptic curve is
    // found. therefor the curve has not been generated.
    // for more details see diploma thesis

    bool prime_proof::prove_downstep(const bigint & r)
    {
	prime_proof *next_test;
	next_test = new prime_proof();
	(*next_test).set_verbose(verbose);
	(*next_test).set_top(false);
	(*next_test).set_ecpp_mode(ecpp_order_mode);
	(*next_test).set_prime(r);
	bool test = (*next_test).prove_prime();
	if(!test)
	{
	    delete next_test;
	    return false;
	}
	down_cert = (*next_test).get_certificate();
	return true;
    }


    certificate prime_proof::ecpp_get_certificate()
    {
	certificate cert;
	cert.set_certificate(cert_vector);
	cert.add(down_cert);
	return cert;
    }



    // verifies a given certificate 
    bool prime_proof::ecpp_verify_proof(base_vector <bigint> cert_vector)
    {
	bool b1,b2,b3;
	D = cert_vector[0];
	p = cert_vector[1];
	order = cert_vector[2];
	prime_factor = cert_vector[3];
	F_q = galois_field(-D);
	gf_element a4(F_q), a6(F_q), Px(F_q), Py(F_q);  
	a4.assign(cert_vector[4]);
	a6.assign(cert_vector[5]);
	Px.assign(cert_vector[6]);
	Py.assign(cert_vector[7]);
	elliptic_curve<gf_element> e_0( a4, a6 );
	point<gf_element> P_0(Px,Py, e_0);
	point<gf_element> Q_0(e_0);
        bigint div_order;
        divide(div_order,order,prime_factor);
	multiply(Q_0, order, P_0);
	b1= Q_0.is_zero();
	multiply(Q_0, div_order, P_0);
	b2= Q_0.is_zero();
	b3=Q_0.on_curve();
	if(verbose && b1 && !b2 && b3){
	    std::cout<<"ecpp:Test success"<<std::endl;
	    return true;
	}
	if(verbose)std::cout<<"ecpp: Test was not successfull"<<std::endl;
	return false;
    }

#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif
