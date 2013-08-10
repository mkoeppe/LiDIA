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
//  File    : gec.cc
//  Author  : Harald Baier (HB)
//	Changes	: See CVS log
//
//==============================================================================================

#ifdef HAVE_CONFIG_H
# include       "config.h"
#endif
# include        "LiDIA/gec.h"
# include        <ctime>
#ifdef VERBOSE
# include        <cassert>
#endif


#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif


//
// initialize the static variables: their meaning is explained in gec.h
//
lidia_size_t gec::
default_lower_bound_bitlength_r         = 159;

bigint gec::
default_upper_bound_k                   = 4;

lidia_size_t gec::
BSI_lower_bound_h_field                 = 200;

lidia_size_t gec::
BSI_lower_bound_extension_degree        = 10000;

lidia_size_t gec::
default_lower_bound_extension_bitlength = 2000;

lidia_size_t gec::
upper_bound_of_degree                   = 600;

int gec::
nr_of_prob_prime_tests                  = 50; // see X9.62, Annex A.2.1

//
// constructor and destructor
//
gec::gec() : out_( std::cout )
{
    p = q                   = 0;
	degree                  = 0;
    k                       = 0;
    r                       = 0;
    lower_bound_bitlength_r = 0;
    upper_bound_k           = 0;
    lower_bound_h_field = h = 0;
	h_field                 = 0;
    delta_field = delta     = 0;
	abs_delta               = 0;
    according_to_BSI        = false;
    CM_FLAG                 = false;
    PC_FLAG                 = false;
	timings                 = false;

	VERBOSE                 = true;
}

gec::gec( std::ostream & out ) : out_( out )
{
    p = q                   = 0;
	degree                  = 0;
    k                       = 0;
    r                       = 0;
    lower_bound_bitlength_r = 0;
    upper_bound_k           = 0;
    lower_bound_h_field = h = 0;
	h_field                 = 0;
    delta_field = delta     = 0;
	abs_delta               = 0;
    according_to_BSI        = false;
    CM_FLAG                 = false;
    PC_FLAG                 = false;
	timings                 = true;

	VERBOSE                 = true;
}


gec::~gec()
{
}

//
// Accessors:
// The methods return 0 if the corresponding attribute is not set.
//
const bigint & 
gec::get_p() const
{
	return( p );
}


const bigint & 
gec::get_q() const
{
	return( q );
}

lidia_size_t
gec::get_degree() const
{
	return( degree );
}

const bigint & 
gec::get_r() const
{
	return( r );
}

const bigint & 
gec::get_k() const
{
	return( k );
}

const point<gf_element> &
gec::get_G() const
{
	return( G );
}

//
// Accessors for the curve parameters:
// The lidia_error_handler is invoked if curve is not initialized.
//
const gf_element &
gec::get_a2() const
{
	if( ! is_initialized )
		lidia_error_handler("gec", 
							"get_a2(): curve not initialized");

	return( a2 );
}

const gf_element & 
gec::get_a4() const
{
	if( ! is_initialized )
		lidia_error_handler("gec", 
							"get_a4(): curve not initialized");
	
	return( a4 );
}

const gf_element & 
gec::get_a6() const
{
	if( ! is_initialized )
		lidia_error_handler("gec", 
							"get_a6(): curve not initialized");
	
	return( a6 );
}

//
// Accessors for default values of bounds of bitlength_r and upper_bound_k
//
lidia_size_t 
gec::get_default_lower_bound_bitlength_r() const
{
	return( default_lower_bound_bitlength_r );
}

const bigint & 
gec::get_default_upper_bound_k() const
{
	return( default_upper_bound_k );
}

//
// Accessors for bounds of bitlength_r and upper_bound_k:
// If they are not set, the method returns 0.
//
lidia_size_t
gec::get_lower_bound_bitlength_r() const
{
	return( lower_bound_bitlength_r );
}


const bigint & 
gec::get_upper_bound_k() const
{
	return( upper_bound_k );
}

//
// Accessors for BSI settings
//
bool 
gec::get_according_to_BSI() const
{
	return( according_to_BSI );
}

lidia_size_t
gec::get_BSI_lower_bound_h_field() const
{
	return( BSI_lower_bound_h_field );
}

lidia_size_t
gec::get_BSI_lower_bound_extension_degree() const
{
	return( BSI_lower_bound_extension_degree );
}

//
// Accessor for default_lower_bound_extension_bitlength
//
lidia_size_t 
gec::get_default_lower_bound_extension_bitlength() const
{
	return( default_lower_bound_extension_bitlength );
}

//
// Accessors for discriminants, class numbers and lower bound of class number 
// of Quot(End(E)):
// They return 0 if the corresponding attribute is not set.
//
const bigint & 
gec::get_delta_field() const
{
	return( delta_field );
}

lidia_size_t
gec::get_h_field() const
{
	return( h_field );
}

lidia_size_t
gec::get_lower_bound_h_field() const
{
	return( lower_bound_h_field );
}

//
// Accessors for discriminant and class number of End(E).
// If they are not computed (e.g. because of PC-method), 0 is returned.
//
const bigint & 
gec::get_delta() const
{
	return( delta );
}

lidia_size_t
gec::get_h() const
{
	return( h );
}

//
// Accessors of verbosity and efficient curve parameters
//
bool
gec::get_verbose_level()
{
	return( VERBOSE );
}
bool
gec::get_efficient_curve_parameters()
{
	return( efficient_curve_parameters );
}

//
// Mutators of field and degree
//
void
gec::set_field( const bigint & field_size )
{
	if( is_prime( field_size, nr_of_prob_prime_tests ) && field_size > 3 )
	{
		q      = field_size;
		degree = 1;
		F_q    = galois_field( field_size );
	}
	else
		lidia_error_handler("gec", "set_field( const bigint & q ): q is not a prime >= 5");	
}

void
gec::set_degree( lidia_size_t d )
{
	if( degree > 0 )
		degree = d;
	else
		lidia_error_handler("gec", "set_degree( lidia_size_t d ): d is not positive");
}

//
// Mutators of security level
//
void 
gec::set_lower_bound_bitlength_r( lidia_size_t bitlength_r )
{
	if( bitlength_r > 0 )
		lower_bound_bitlength_r = bitlength_r;
	else
		lidia_error_handler("gec", "set_lower_bound_bitlength_r( lidia_size_t bitlength_r ): bitlength_r is not positive");
}

void
gec::set_upper_bound_k( const bigint & bound_k )
{
	if( bound_k > 0 )
		upper_bound_k = bound_k;
	else
		lidia_error_handler("gec", "set_upper_bound_k( const bigint & bound_k ): bound_k is not positive");
}

//
// Mutators of default security bounds
//
void
gec::set_default_lower_bound_bitlength_r( lidia_size_t bitlength_r )
{
	if( bitlength_r > 0 )
		default_lower_bound_bitlength_r = bitlength_r;
	else
		lidia_error_handler("gec", "set_default_lower_bound_bitlength_r( lidia_size_t bitlength_r ): bitlength_r is not positive");
}

void 
gec::set_default_upper_bound_k( const bigint & bound_k )
{
	if( bound_k > 0 )
		default_upper_bound_k = bound_k;
	else
		lidia_error_handler("gec", "set_default_upper_bound_k( const bigint & bound_k ): bound_k is not positive");
}

void 
gec::set_default_lower_bound_extension_bitlength( lidia_size_t default_extension_bitlength )
{
	if( default_extension_bitlength > 0 )
		default_lower_bound_extension_bitlength = default_extension_bitlength;
	else
		lidia_error_handler("gec", "set_default_lower_bound_extension_bitlength( lidia_size_t default_extension_bitlength ): default_extension_bitlength is not positive");
}

//
// Mutators of BSI-attributes
//
void
gec::set_according_to_BSI( bool bsi )
{
	according_to_BSI = bsi;
}

void 
gec::set_BSI_lower_bound_h_field( lidia_size_t bsi_h_field )
{
	if( bsi_h_field > 0 )
		BSI_lower_bound_h_field = bsi_h_field;
	else
		lidia_error_handler("gec", "set_BSI_lower_bound_h_field( lidia_size_t bsi_h_field ): bsi_h_field is not positive");
}

void 
gec::set_BSI_lower_bound_extension_degree( lidia_size_t bsi_degree )
{
	if( bsi_degree > 0 )
		BSI_lower_bound_extension_degree = bsi_degree;
	else
		lidia_error_handler("gec", "set_BSI_lower_bound_extension_degree( lidia_size_t bsi_degree ): bsi_degree is not positive");
}

	//
	// Mutators of verbosity and efficient curve parameters
	//
void
gec::set_verbose_level( bool verbose_level )
{
	VERBOSE = verbose_level;
}

void
gec::set_efficient_curve_parameters( bool efficiency )
{
	efficient_curve_parameters = efficiency;
}

//
// High level functions: 
// * compute_delta_field
// * get_lower_bound_h_field
// * check_the_settings

//
// delta_field is the maximal imaginary quadratic discriminant dividing either
// delta (CM-approach) or
// (q+1-k*r)^2 - 4*q (PC-approach)
//
bool
gec::compute_delta_field()
{
	// ensure that q, k and r are known in case of random approach
	if( PC_FLAG && ( q == 0 || k == 0 || r == 0 ) )
		lidia_error_handler("gec", 
							"compute_delta_field(): q, k, or r not initialized");

	// ensure that delta is known in case of CM-approach
	if( CM_FLAG && delta.is_zero() )
	{
		lidia_error_handler( "gec_complex_multiplication",
							 "compute_delta_field(): delta not initialized.");
	}

	if( VERBOSE )
		std::cout << "Computing delta_field...." << std::endl;
	
	// 
	// temporary variables
	//
	bigint tmp_delta, tmp_bigint;

	// if CM-approach: initialize tmp_delta with delta
	if( CM_FLAG )
		tmp_delta.assign( delta );
	// if PC-approach: initialize tmp_delta with (q+1-k*r)^2 - 4*q
	else if( PC_FLAG )
	{
		multiply( tmp_delta, k, r );
		tmp_delta.negate();
		add( tmp_delta, tmp_delta, q + bigint(1) );
		square( tmp_delta, tmp_delta );
		shift_left( tmp_bigint, q, 2 );
		subtract( tmp_delta, tmp_delta, tmp_bigint );
	}
	else
		lidia_error_handler("gec", 
							"compute_delta_field(): PC_FLAG and CM_FLAG not set");

	//
	// find factorization of tmp_delta
	//
	rational_factorization delta_factored = 
		factor( tmp_delta );

    // if current factorization is not prime factorization, refine it
	if( ! delta_factored.is_prime_factorization() )
	{
		delta_factored.trialdiv();
		delta_factored.ecm();
	   
		int counter = 0;
		
		// We refine it at most 5 times. If we are not successful
		// after 5 trials, we terminate with an error message.
		for( int counter = 0; counter < 5; counter++ )
		{
			delta_factored.refine();
		   
			if( delta_factored.is_prime_factorization() ) 
				break;
		   
			delta_factored.trialdiv();
			delta_factored.ecm();
		}

		if( counter == 5 )
		{
			std::cout << "gec::" <<
				"compute_delta_field(): " <<
				"cannot find prime factorization" <<
				"of discriminant.";
			delta_field = 0;
		   
			return( false );
		}
	}

	// assign the field discriminant:
	// it is the maximal imaginary quadratic discriminant dividing
	// delta_factored (hence delta_field is minimal with respect
	// to the absolute value

	// tmp_delta will store the accumulated product of the prime
	// factors of delta_field
	tmp_delta = -1;

	// find the prime divisors of delta_field by choosing the
	// primes dividing delta_factored with an odd power
	for( lidia_size_t i = 0; i < delta_factored.no_of_comp(); i++ )
	{
		if( delta_factored.exponent(i)%2 )
			multiply( tmp_delta, tmp_delta, delta_factored.base(i) );
	}	      

	// we next have to find the correct power of 2 dividing delta_field
	abs_delta.assign( tmp_delta );
	abs_delta.negate();

	// If tmp_delta is even we currently have tmp_delta = 2 mod 4
	// As tmp_delta and delta_field differ by a factor f^2 (f the conductor)
	// we have delta_field = 4*tmp_delta
	if( ! abs_delta.bit( 0 ) )
		shift_left( tmp_delta, tmp_delta, 2 );
	// if tmp_delta is odd we have to ensure delta_field = 0 or 1 mod 4
	else
	{
		// abs_delta = 1 mod 4 ==> tmp_delta = -abs_delta = 3 mod 4
		if( ! ( abs_delta.bit( 1 ) ) )
			shift_left( tmp_delta, tmp_delta, 2 );
	}
	
	delta_field.assign( tmp_delta );

	if( VERBOSE )
		std::cout << "delta_field = " << delta_field << std::endl;
	
	return( true );
}

//
// Depending on the order of magnitude of |delta_field|
// we proceed as follows:
// * if |delta_field| < 10^9 we make use of the function
//   class_number() of the class quadratic_order
// * otherwise we choose a random element x of the class group
//   and show that x^k != 1 for all 1 <= k <= 200
lidia_size_t
gec::compute_lower_bound_h_field()
{
	// if |delta_field| is small we compute h_field exactly
	if( abs( delta_field ) < bigint( 1000000000 ) )
	{
		// to test whether class number fits in a long
		long tmp_l;

		quadratic_order QO;
		QO.assign( delta_field );

		// If h_field fits in a long assign it to lower_bound_h_field.
		// Otherwise set it to the BSI default value.
		if( QO.class_number().longify( tmp_l ) )
			lower_bound_h_field = BSI_lower_bound_h_field;
		else
			lower_bound_h_field = tmp_l;

		return( lower_bound_h_field );
	}

	// If |delta_field| is large we make use of a random element
	// of the class group and try to show its order to be at least
	// the BSI value.
	//
	// The basic idea is as follows: For each rational prime p <= 1000
	// we try to find a quadratic form of discriminant delta_field
	// of the form (p,*.*). We test whether its order is at least
	// the BSI bound.
	//
	// We proceed as follows: If delta_field is even, we set b_0 = 2.
	// Otherwise, we set b_0 = 3 (hence b_0 is minimal with corresponding
	// form different from unit form).
	// Next, we compute the fraction fraction = (b_0*b_0 - delta_field)/4.
	// We test for each prime divisor p of fraction (starting with
	// p = 2 or p = 3 in case of odd or even delta_field) whether the
	// order of (p, b_0, fraction/p) is at least the BSI bound.
	//
	// If p is not a prime divisor of fraction, we try to find a form
	// (p,*,*) be testing whether delta_field mod p is a square mod p.

	// assign unit_form to the identity element in class group
	quadratic_form unit_form;
	unit_form.assign_one( delta_field );

	// quadratic forms for computations
	quadratic_form qf_current, qf_initial;

	// the parameters of the current quadratic form (a,b,c)
	bigint a, b, c;
	// b_0 is equal to 2 if delta_field is even, and 3 otherwise
	bigint b_0;
	// delta_field mod p
	bigint delta_field_mod_p;
	// r will be a root of delta_field_mod_p
	bigint r;

	// the parameter a will be the prime prime
	bigint prime;
	// fraction stores (b_0*b_0 - delta_field)/4.
	bigint fraction;
	// fraction_prime is fraction mod prime
	bigint fraction_mod_prime;

	// current multiple stores the current scalar factor of the form
	int current_multiple = 0;
	// boolean to store whether current form will be tested
	bool test_this_form;

	// set b_0 = 2 if delta_field is even, and 1 otherwise
	// the first tested parameter a will be set to next_prime( b_0 + 1 )
	if( delta_field.bit( 0 ) == 0 )
	{
		b_0 = 2;
		prime = 3;
	}
	else
	{
		b_0 = 1;
		prime = 2;
	}

	// compute the fraction; a will be set to prime factors of factor
	fraction = ( b_0 * b_0 - delta_field ) / 4 ;

	if( VERBOSE )
	{
		std::cout << std::endl << "***************************************"
				  << std::endl;
		std::cout << "delta_field = " << delta_field << std::endl;
		std::cout << "b_0 = " << b_0 << std::endl;
		std::cout << "(b_0*b_0-delta_field)/4 = " << fraction << std::endl;
	}

	//
	// The loop to find a quadratic form of order
	// at least BSI_lower_bound_h_field.
	//
	while( prime < 1000 && current_multiple < BSI_lower_bound_h_field )
	{
		test_this_form = false;

		if( VERBOSE )
			std::cout << std::endl << "prime = " << prime << std::endl;

		// compute whether current prime divides fraction
		remainder( fraction_mod_prime, fraction, prime );
		
		// if yes we will make use of the form (prime, b_0, fraction/prime)
		if( fraction_mod_prime == 0 )
		{
			a.assign( prime );
			b.assign( b_0 );

			test_this_form = true;
		}
		// if not we test whether there is a form (prime,*,*)
		else if( prime >= 3 )
		{
			remainder( delta_field_mod_p, delta_field, prime );
			// delta_field = b^2 - 4 * prime * c = b^2 mod prime
			// hence delta_field has to be square mod prime
			if( jacobi( delta_field_mod_p, prime ) != -1 )
			{
				// r = sqrt( delta_field ) mod p
				ressol( r, delta_field_mod_p, prime );
				
				// we have to ensure that 4 | r^2 - delta_field
				if( r.bit( 0 ) != delta_field.bit( 0 ) )
					subtract( r, prime, r );
				
				a.assign( prime );
				b.assign( r );
						   
				test_this_form = true;
			   
				if( VERBOSE )
					std::cout << "r = " << r << std::endl;
			}
		}

		// if we have found a quadratic form of discriminant delta_field,
		// we test whether multiples up to BSI bound of it are the unit form
		if( test_this_form == true )
		{
			// assign c
			c = ( b * b - delta_field ) / (4 * a );

			// initialize the initial and current form with (a,b,c)
			qf_initial.assign( a, b, c );
			current_multiple = 1;
			qf_initial.reduce();
			qf_current.assign( qf_initial );

			if( VERBOSE )
			{
				std::cout << "Current initial form qf = " << 
					qf_current << std::endl;   
				std::cout << "delta(qf) = " << 
					qf_current.discriminant() << std::endl;
			}
			if( qf_current.discriminant() != delta_field )
				lidia_error_handler( "gec", "compute_lower_bound_h_field(): wrong discrimant of quadratic form." );

			// test for the multiples of qf_initial to be different
			// from unit form
			while( ! qf_current.is_equal( unit_form ) )
			{
				if( current_multiple >= BSI_lower_bound_h_field )
					break;

				compose_reduce( qf_current, qf_current, qf_initial );
				current_multiple++;
			}
		}

		prime = next_prime( prime );
	}

	if( current_multiple < BSI_lower_bound_h_field )
	{
		if( VERBOSE )
			std::cout << "We cannot show that h_field is at least " 
					  << BSI_lower_bound_h_field << " ." << std::endl;

	   
		return( 0 );
	}

	if( VERBOSE )
		std::cout << "h_field is at least " << current_multiple << " ." << std::endl;

	return( current_multiple );
}

//
// check_the_settings:
// this method tests whether the necessary security are set.
// If no, the method uses default values.
// If yes, it tests whether they fit.
//
void
gec::check_the_settings()
{
	// if no degree was set, generate a curve over a prime field
	if( degree <= 0 )
		degree = 1;

	// test whether the field is already set
	if( q.is_zero() )
	{
		// test whether the lower bound of the bitlength of r is set
		if( lower_bound_bitlength_r == 0 )
			lower_bound_bitlength_r = default_lower_bound_bitlength_r;

		// test whether the upper bound of of k is set
		if( upper_bound_k == 0 )
			upper_bound_k = default_upper_bound_k;
	}
	// the field has already been set
	else
	{
		// test whether bitlength of bound of r and q fit
		if( lower_bound_bitlength_r > q.bit_length() )
			lidia_error_handler("gec", "check_the_settings(): lower_bound_bitlength_r > q.bit_length()");

		// if bound of k is not set use the default value
		if( upper_bound_k == 0 )
			upper_bound_k = default_upper_bound_k;

		// if bound of bitlength of r is not set use the default value of k
		if( lower_bound_bitlength_r == 0 )
			lower_bound_bitlength_r = q.bit_length() - 
				( upper_bound_k.bit_length() - 1 );
	}

	//
	// check if the set bounds match the BSI requirements
	//
	if( ( lower_bound_bitlength_r < default_lower_bound_bitlength_r || 
		  upper_bound_k > default_upper_bound_k ) && 
		according_to_BSI == true )
		lidia_error_handler( "gec", "check_the_settings()::the set values do not satisfy the BSI conditions" );

	//
	// if we use the CM approach check the following:
	// * delta = 0,1 mod 4
	// * delta and upper_bound_k fit
	//
	if( CM_FLAG && ! delta.is_zero() )
	{
		// check whether delta is negative
		if( delta.is_ge_zero() )
			lidia_error_handler("gec",
								"check_the_settings(): delta is not negative");
		
		// check whether delta = 0,1 mod 4
		abs_delta.assign( delta );
		abs_delta.negate();
	
		if( ( abs_delta.bit(0) == 1 && abs_delta.bit(1) == 0 ) ||
			( abs_delta.bit(0) == 0 && abs_delta.bit(1) == 1 ) )
		{
			lidia_error_handler("gec",
								"check_the_settings(): delta is not = 0,1 mod 4");
		}
		
		//
		// check whether delta and upper_bound_k fit
		//

		// upper_bound_k == 1  ==>  delta = 5 mod 8  ==>  abs_delta = 3 mod 8
		if( upper_bound_k.is_one() )
		{
			if(	! ( abs_delta.bit(0) == 1 && abs_delta.bit(1) == 1 &&
					abs_delta.bit(2) == 0 ) )
				
				lidia_error_handler("gec",
									"check_the_settings(): upper_bound_k = 1 AND delta is not = 5 mod 8");
		}
		// delta = 1 mod 8  or  delta = 0,4 mod 16 ==>  upper_bound_k >= 4
		else if( ( abs_delta.bit(0) == 1 && abs_delta.bit(1) == 1 &&
				   abs_delta.bit(2) == 1 ) || 
				 ( abs_delta.bit(0) == 0 && abs_delta.bit(1) == 0 &&
				   abs_delta.bit(2) == 0 && abs_delta.bit(3) == 0 ) || 
				 ( abs_delta.bit(0) == 0 && abs_delta.bit(1) == 0 &&
				   abs_delta.bit(2) == 1 && abs_delta.bit(3) == 1 ) )
		{
			if( upper_bound_k < 4 )
				lidia_error_handler("gec",
									"check_the_settings(): ( delta = 1 mod 8 OR delta = 0,4 mod 16 AND  upper_bound_k < 4 ) ");
			
		}
	}
}

#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif
