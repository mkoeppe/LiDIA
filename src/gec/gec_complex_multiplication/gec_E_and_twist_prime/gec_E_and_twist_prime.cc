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
//  File    : gec_E_and_twist_prime.cc
//  Author  : Harald Baier (HB)
//	Changes	: See CVS log
//
//==============================================================================================

#ifdef HAVE_CONFIG_H
# include       "config.h"
#endif
# include        "LiDIA/gec_E_and_twist_prime.h"
#ifdef VERBOSE
# include        <cassert>
#endif


#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif


//
// Initialize static variables
//
bigint gec_E_and_twist_prime::default_delta_E_and_twist_prime = -356131;
short gec_E_and_twist_prime::default_delta_case_E_and_twist_prime = 1; 

//
// Constructor / Destructor
//
gec_E_and_twist_prime::
gec_E_and_twist_prime()
{
	r_tw.assign_zero();

	k.assign_one();
	k_tw.assign_one();

	a2_tw.assign_zero();
	a4_tw.assign_zero();
	a6_tw.assign_zero();
}

gec_E_and_twist_prime::
gec_E_and_twist_prime( std::ostream & out ) : gec_complex_multiplication( out )
{
	r_tw.assign_zero();

	k.assign_one();
	k_tw.assign_one();

	a2_tw.assign_zero();
	a4_tw.assign_zero();
	a6_tw.assign_zero();
}


gec_E_and_twist_prime::
~gec_E_and_twist_prime()
{
}


//
// Accessors: The lidia_error_handler is called if
// is_initialized is false.
//
const bigint &
gec_E_and_twist_prime::get_r_tw() const
{
	return( r_tw );
}

const bigint & 
gec_E_and_twist_prime::get_k_tw() const
{
	return( k_tw );
}

const point<gf_element> &
gec_E_and_twist_prime::get_G_tw() const
{
	if( ! is_initialized )
		lidia_error_handler("gec_E_and_twist_prime", 
							"get_G_tw(): twist not initialized");

	return( G_tw );
}

//
// Accessors for the twist parameters
//

const gf_element &
gec_E_and_twist_prime::get_a2_tw() const
{
	return( a2_tw );
}

const gf_element & 
gec_E_and_twist_prime::get_a4_tw() const
{
	return( a4_tw );
}

const gf_element & 
gec_E_and_twist_prime::get_a6_tw() const
{
	return( a6_tw );
}


//
// is_cryptographically_strong() - 
// decide whether or not curve and twist are strong
//
bool
gec_E_and_twist_prime::is_cryptographically_strong()
{
	if( VERBOSE )
		std::cout << "Testing for cryptographical strength..." << std::endl;

	// test whether r is large enough
	if( VERBOSE )
		std::cout << "Testing whether r is large enough....\t";

	if( r.bit_length() < lower_bound_bitlength_r )
	{
		if( VERBOSE )
		{
			std::cout << "r is not large enough!";
			std::cout << "We test a new curve." << std::endl;
		}
		
		return( false );
	}

	if( r_tw.bit_length() < lower_bound_bitlength_r )
	{
		if( VERBOSE )
		{
			std::cout << "r_tw is not large enough!";
			std::cout << "We test a new curve." << std::endl;
		}

		return( false );
	}

	if( VERBOSE )
		std::cout << "Yes! r and r_tw are large enough." << std::endl;

	// test for anomalous curve
	if( VERBOSE )
		std::cout << "Testing for anomalous curve....\t"; 
	if( r == q )
	{
		if( VERBOSE )
			std::cout << "Curve is anomalous!\t We test a new curve."
					  << std::endl;

		return( false ); 
	}

	if( r_tw == q )
	{
		if( VERBOSE )
			std::cout << "Twist is anomalous!\t We test a new curve."
					  << std::endl;

		return( false ); 
	}

	if( VERBOSE )
		std::cout << "Curve and twist are not anomalous!" << std::endl;

	// test for MOV condition
	if( VERBOSE )
		std::cout << "Testing for MOV-condition....\n";
  
	lidia_size_t lower_bound_extension_degree;

	if( according_to_BSI )
		lower_bound_extension_degree = BSI_lower_bound_extension_degree;
	else
		lower_bound_extension_degree = 
			default_lower_bound_extension_bitlength / ( q.bit_length() - 1);

	if( VERBOSE )
		std::cout << "Extension degree should be at least " <<
			lower_bound_extension_degree << std::endl;
 
	bigint tmp_bigint( 1 ), tmp_bigint_tw( 1 );

	for(int i = 1; i <= lower_bound_extension_degree; i++)
	{
		multiply( tmp_bigint, tmp_bigint, q );
		multiply( tmp_bigint_tw, tmp_bigint_tw, q );

		remainder( tmp_bigint, tmp_bigint, r );
		remainder( tmp_bigint_tw, tmp_bigint_tw, r_tw );

		if( tmp_bigint == 1 )
		{
			if( VERBOSE )
			{
				std::cout << "Extension degree of curve is just " 
					 << i << " ." << std::endl;
				std::cout << "We test a new curve!" << std::endl;
			  
			}

			return( false );
		}

		if( tmp_bigint_tw == 1 )
		{
			if( VERBOSE )
			{
				std::cout << "Extension degree of twist is just " 
					 << i << " ." << std::endl;
				std::cout << "We test a new curve!" << std::endl;
			  
			}

			return( false );
		}

	}

	if( VERBOSE )
	{
		std::cout << "Extension degree of curve and twist are at least " 
			 << lower_bound_extension_degree;
		std::cout << " . MOV-attack is not feasible." << std::endl;
	}

	//
	// test whether the class number is large enough
	// i.e. ensure h_field >= BSI_lower_bound_h_field
	//
	if( according_to_BSI )
	{
		if( VERBOSE )
			std::cout << "Testing h_field is at least " 
				 << BSI_lower_bound_h_field << "...\n";

		//
		// We do not have to handle the CM-case, as delta_field and h_field
		// are already computed in this case
		//
		if( delta_field == 0 )
			compute_delta_field();

		if( lower_bound_h_field == 0 )
			if( ! compute_lower_bound_h_field() )
				return( false );
	}

	return( true );
}


#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif
