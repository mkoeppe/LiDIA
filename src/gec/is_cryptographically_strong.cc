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
//  File    : is_cryptographically_strong.cc
//  Author  : Harald Baier (HB)
//	Changes	: See CVS log
//
//==============================================================================================


#ifdef HAVE_CONFIG_H
# include       "config.h"
#endif
# include        "LiDIA/gec.h"
#ifdef VERBOSE
# include        <cassert>
#endif


#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif


//
// is_cryptographically_strong():
// the method decides whether the curve to test is
// suitable for use in cryptography or not
//

// Currently, the following tests are performed:
// * r > 2^(lower_bound_bitlength_r - 1)
// * k < upper_bound_k
// * curve is not anomalous
// * MOV/FR-reduction: dependent on according_to_BSI
// * class number of maximal order containing End(E) is large enough
// * degree is either 1 or a prime (in case of p = 2)
//
bool
gec::is_cryptographically_strong( const bigint & order )
{
	if( VERBOSE )
	{
		std::cout << std::endl;
		std::cout << "Testing for cryptographical strength..." << std::endl;
		std::cout << "order = " << order << std::endl;
		std::cout << "delta = " << delta << std::endl;
	}

	//
	// temporary variable for r
	//
	bigint tmp_r( order );

	// initialize cofactor with 1
	k.assign_one();

	if( VERBOSE )
	{
		std::cout << "Computing r and k." << std::endl;
		std::cout << "Testing whether k <= " << upper_bound_k << "....\t";
	}
  
	//
	// factor out the even part of tmp_r = order
	//

	// tmp_int is the exponent of 2-power dividing order
	unsigned int tmp_int = 0;
	while( tmp_r.bit( tmp_int ) == 0 )
		tmp_int++;

	// divide order by 2^tmp_int and assign it again to tmp_r
	shift_right( tmp_r, tmp_r, tmp_int );
	// multiply k by 2^tmp_int
	shift_left( k, k, tmp_int );
	// if k > upper_bound_k return false
	if( k > upper_bound_k )
	{
		if( VERBOSE )
			std::cout << "k = " << k << " is too big." << std::endl;
      
		return( false );
	}

	//
	// factor out the 3-part of tmp_r
	//
	while( remainder( tmp_r, 3 ) == 0 )
	{
		divide( tmp_r, tmp_r, 3 );
		multiply( k, k, 3 );
		if( k > upper_bound_k )
		{
			if( VERBOSE )
				std::cout << "k = " << k << " is too big." << std::endl;
      
			return( false );
		}
	}

	//
	// if tmp_r is not prime and upper_bound_k < 5, the curve
	// is not cryptographically strong. Hence return( false )
	// Otherwise factor it. We use trial division
	// in the range [5, upper_bound_k]. 
	//  
	if( ! is_prime( tmp_r, nr_of_prob_prime_tests ) )
	{
		if( upper_bound_k < 5 )
			return( false );
		else
		{
			// assign upper_bound_k to an integer: if this is not possible
			// upper_bound_k is too large and we return false.
			int upper_bound_k_as_int;
			if( upper_bound_k.intify( upper_bound_k_as_int ) )
			{
				if( VERBOSE )
					std::cout << "Upper bound of k is too big." << std::endl;

				return( false );
			}     

			// use tmp_r_factored to store factorization of tmp_r
			rational_factorization tmp_r_factored = 
				trialdiv( tmp_r, upper_bound_k_as_int, 5 );

		// if we cannot factor tmp_r by trialdivision, return false
			if( ! tmp_r_factored.is_prime_factorization() )
			{
				if( VERBOSE )
					std::cout << "Cannot factor order by trial division. "
							  << std::endl;
	  
				return( false );
			}      
		
			// assign the larges prime factor of order to r and compute k
			r = tmp_r_factored.base( tmp_r_factored.no_of_comp() - 1 );
			divide( tmp_r, tmp_r, r);
			multiply( k, k, tmp_r );

			// if k > upper_bound_k return false
			if( k > upper_bound_k )
			{
				if( VERBOSE )
					std::cout << "k = " << k << " is too big." << std::endl;
      
				return( false );
			}
		}
	}
	else
		r.assign( tmp_r );

	if( VERBOSE )
		std::cout << "k = " << k << " is not too large." << std::endl;
	
	//
	// test whether r is large enough
	//
	if( VERBOSE )
		std::cout << "Testing whether r is large enough....\t";

	if( r.bit_length() < lower_bound_bitlength_r )
	{
		if( VERBOSE )
			std::cout << "r is not large enough!" << std::endl;

		return( false );
	}
  
	if( VERBOSE )
		std::cout << "Yes! r is large enough." << std::endl;

	//
	// test for anomalous curve
	//
	if( VERBOSE )
		std::cout << "Testing for anomalous curve....\t"; 
	
	if( r == p )
	{
		if( VERBOSE )
			std::cout << "Curve is anomalous!\t" << std::endl;

		return( false ); 
	}
	
	if( VERBOSE )
		std::cout << "Curve is not anomalous!" << std::endl;

	//
	// test if MOV/FR-reduction is feasible
	//  
	if( VERBOSE )
		std::cout << "Testing for MOV/FR-reduction....\n";

	// lower_bound_extension_degree stores a lower bound of the
	// field degree over F_q where the MOV/FR-reduction maps to
	lidia_size_t lower_bound_extension_degree;

	// if BSI requirements are relevant, take
	// BSI_lower_bound_extension_degree as lower bound
	if( according_to_BSI )
		lower_bound_extension_degree = BSI_lower_bound_extension_degree;
	// otherwise the DLP in fields of bitlength
	// default_lower_bound_extension_bitlength
	// is considered to be intractable
	else
		lower_bound_extension_degree = 
			default_lower_bound_extension_bitlength / ( q.bit_length() - 1);

	if( VERBOSE )
		std::cout << "Extension degree should be at least " <<
			lower_bound_extension_degree << std::endl;

	// tmp_bigint is a temporary variable storing the current power q^i mod r
	bigint tmp_bigint;
	tmp_bigint.assign_one();

	// test if q^i = 1 mod r
	for(int i = 1; i <= lower_bound_extension_degree; i++)
	{
		multiply( tmp_bigint, tmp_bigint, q );
		remainder( tmp_bigint, tmp_bigint, r );

		if( tmp_bigint.is_one() )
		{
			if( VERBOSE )
			{
				std::cout << "Extension degree of MOV/FR-reduction "
						  << " is only " << i << " ." << std::endl;
			}

			return( false );
		}
	}

	if( VERBOSE )
	{
		std::cout << "Extension degree is at least "
				  << lower_bound_extension_degree;
		std::cout << " . MOV/FR-attack is not feasible." << std::endl;
	}

	//
	// test whether the class number is large enough
	// i.e. ensure h_field >= BSI_lower_bound_h_field
	//

	//
	// We do not have to handle the CM-case at this point
	// as delta_field and h_field are already computed:
	// see gec_complex_multiplication/generate.cc
	//
	if( PC_FLAG && according_to_BSI )
	{
		if( VERBOSE )
			std::cout << "Testing whether h_field is at least " 
					  << BSI_lower_bound_h_field << "...\n";
		
		compute_delta_field();
		lower_bound_h_field = compute_lower_bound_h_field();
		
		if( lower_bound_h_field < BSI_lower_bound_h_field )
			return( false );
	}
	
	// if p = 2, test for Weil-descent:
	// if char(F_q) = 2, then the degree of F_q over F_2 has to be prime
	if( q.is_even() && degree > 1 &&
	    ! is_prime( bigint( degree ), nr_of_prob_prime_tests ) )
	    return( false );
	
	return( true );
}


#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif
