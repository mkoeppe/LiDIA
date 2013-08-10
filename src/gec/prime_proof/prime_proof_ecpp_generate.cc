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
//  File    : prime_proof_ecpp_generate.cc
//  Author  : Harald Baier (HB),
//	      small changes by Jochen Hechler (JH)
//      Changes : See CVS log
//
//==============================================================================================



#include        "LiDIA/prime_proof.h"


#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif


//
// generate() implements the virtual function of the base class
// in the scope of the complex multiplication approach
//

    bool prime_proof::generate_curve()
    {
	// the reduced representatives are stored in the base_vector class_group
	class_group.set_mode(EXPAND);
	class_group = compute_class_group( D ); 
	set_generation_mode();
	set_complex_precision();
	compute_class_polynomial();

	//
	// Compute some irreducible factor of degree of class polynomial
	//
	// initialize the class polynomial mod p

	Fp_polynomial class_polynomial_mod_p( class_polynomial, n );
	bigint root_of_class_polynomial_mod_p = find_root( class_polynomial_mod_p ); 

	gf_element class_invariant_mod_p( F_q );
	class_invariant_mod_p.assign( root_of_class_polynomial_mod_p );

	
	   
	// Determine the curve parameters and a base point:
	// This is accomplished by the function 
	// assign_class_invariant_to_curve.
	if( ! assign_class_invariant_to_curve( class_invariant_mod_p ) )
	{
	    std::cout<<"Could not assign class_invariant_mod_p to curve. Returning FALSE"<<std::endl;
	    return false;
	}

	return true;
    }

#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif
