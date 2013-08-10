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
//  File    : prime_proof_ecpp_set_generation_mode.cc
//  Author  : Harald Baier (HB),
//	      small changes by Jochen Hechler (JH)
//      Changes : See CVS log
//
//==============================================================================================

#include       "LiDIA/prime_proof.h"

// The Code was taken from the gec package
// Originally by Harald Baier (HB)
// Some Changes due to variables have been apllied

#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif


// Determine the generation mode which requires the lowest floating point precision 

    void prime_proof::set_generation_mode()
    {
	// compute the remainders of D modulo 3 and 8, respectively
	long remainder_delta_mod_3 = remainder( D, (long) 3 );
	long remainder_delta_mod_8 = remainder( D, (long) 8 );
	// if no mode is explicitely set, set as follows:
	// delta % 3 != 0 && delta % 8 = 1: generation_mode = 4 (polynomial W)
	// delta % 3 != 0 && delta % 8 != 8: generation_mode = 3 (polynomial G)
	// delta % 3 = 0: generation_mode = 1 (polynomial H)
	if( remainder_delta_mod_3 )
	{  
	    if( remainder_delta_mod_8 == -7 ) 
		generation_mode = 4;    
	    else
		generation_mode = 3;
	}
	else generation_mode = 1;
    }

#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif
