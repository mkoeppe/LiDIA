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
//  File    : prime_proof_ecpp_primelist.cc
//  Author  : Jochen Hechler (JH)
//      Changes : See CVS log
//
//==============================================================================================

#include	<iostream>
#include	"LiDIA/prime_proof.h"
#include	"LiDIA/bigint.h"
#include	"LiDIA/bigfloat.h"
#include	"LiDIA/base_vector.h"
#include	"LiDIA/LiDIA.h"
#include	"LiDIA/error.h"
#include	<time.h>


#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif


    void prime_proof::calculate_prime_list()
    {
	if(verbose)std::cout<<"ecpp: calculating primelist"<<std::endl;
	prime_list.set_capacity(pl_length);
	start_list = s.next_prime();
	bigint tmp = start_list;
	bigint tmp2;
	int difference;
	for(int i = 0; i < pl_length ; i++)
	{
	    tmp2 = tmp.next_prime();
	    tmp = tmp2 - tmp; 	
	    tmp.intify(difference);
	    prime_list[i] = difference;
	    tmp = tmp2;
		
	}
	if(verbose)std::cout<<"ecpp: Primelist calculated"<<std::endl;
    }


    void prime_proof::calculate_lower_bound_s ()
    {
	bigfloat t = n;
	bigfloat tmp = sqrt(t);
	bigfloat tmp2 = sqrt(tmp);
	tmp = tmp + (2*tmp2) + 2;
	tmp.bigintify(s);	
	if(verbose)std::cout<<"ecpp: lower bound s: "<<s<<std::endl;
    }

#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif
