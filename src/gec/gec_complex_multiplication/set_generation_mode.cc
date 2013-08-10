// -*- C++ -*-
//==============================================================================================
//
//	This file is part of LiDIA --- a library for computational number theory
//
//	Copyright (c) 1994--2002 the LiDIA Group.  All rights reserved.
//
//	See http://www.informatik.tu-darmstadt.de/TI/LiDIA/
//
//----------------------------------------------------------------------------------------------
//
//
//  File    : set_generation_mode.cc
//  Author  : Harald Baier (HB)
//	Changes	: See CVS log
//
//==============================================================================================

#ifdef HAVE_CONFIG_H
# include       "config.h"
#endif
# include        "LiDIA/gec_complex_multiplication.h"
#ifdef VERBOSE
# include        <cassert>
#endif


#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif

//
// mode = 0: Determine the generation mode which requires the 
// lowest floating point precision
// mode > 0: Set generation_mode <- mode, if possible
//
void
gec_complex_multiplication::set_generation_mode( unsigned int mode )
{
	// if delta is not initialized, invoke the lidia_error_handler
	if( delta == 0 )
		lidia_error_handler("gec_complex_multiplication",
							"set_generation_mode: delta is not initialized");

	// compute the remainders of delta modulo 3 and 8, respectively
	long remainder_delta_mod_3 = remainder( delta, (long) 3 );
	long remainder_delta_mod_8 = remainder( delta, (long) 8 );

	if( VERBOSE )
		std::cout << "Setting the generation mode....";
	
	//
	// set the generation mode
	// 

	// if no mode is explicitely set, set as follows:
	// delta % 3 != 0 && delta % 8 = 1: generation_mode = 4 (polynomial W)
	// delta % 3 != 0 && delta % 8 != 8: generation_mode = 3 (polynomial G)
	// delta % 3 = 0: generation_mode = 1 (polynomial H)
	if( mode == 0 )
	{
		if( remainder_delta_mod_3 )
		{  
			if( remainder_delta_mod_8 == -7 ) 
				generation_mode = 4;    
			else
				generation_mode = 3;
		}
		else
			generation_mode = 1;     
	}
	// otherwise try to assign generation_mode <- mode
	else
	{
		if( mode == 3 )
		{
			if( remainder_delta_mod_3 == 0 )
			{
				std::cout << "Cannot assign generation_mode <- 3 " <<
					"as delta is divisible by 3." << std::endl;
				std::cout << "We set generation_mode <- 1." << std::endl; 
			  
				generation_mode = 1;
			}
			else
				generation_mode = mode;
		}
		else if( mode == 4 )
		{
			if( remainder_delta_mod_3 == 0 || remainder_delta_mod_8 != -7 )
			{
				std::cout << "Cannot assign generation_mode <- 4 " <<
					"as delta is either divisible by 3" << std::endl;
				std::cout << "or delta != 1 mod 8. " << 
					"We set generation_mode <- 1." << std::endl; 
				
				generation_mode = 1;
			}
			else
				generation_mode = mode;
		}
		else if( mode == 1 )
			generation_mode = 1;

		// if mode \notin {1,3,4} output an error message
		else
		{
			std::cout << "We do not support generation_mode = " 
					  << generation_mode << std::endl;
			std::cout << "We set generation_mode <- 1." << std::endl; 
			generation_mode = 1;
		}
	}
	
	if( VERBOSE )
		std::cout << "generation_mode = " 
				  << generation_mode << std::endl;
}


#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif
