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
//  File    : modular_functions.cc
//  Author  : Harald Baier (HB)
//	Changes	: See CVS log
//
//==============================================================================================

#ifdef HAVE_CONFIG_H
# include       "config.h"
#endif

# include        "LiDIA/modular_functions.h"

#ifdef VERBOSE
# include        <cassert>
#endif

#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif


/********************************************************************
 *
 * Dedekind's eta-function computed via the Euler sum.
 * See PhD for details.
 * The current implementation corresponds to the algorithm
 * computeEtaViaEulerSum.
 *
 *******************************************************************/
//
// dedekind_eta bases on the algorithm of PhD
//
bigcomplex 
dedekind_eta( const bigcomplex & tau, bool is_fundamental )
{
	// complex variables: 
	// * sum: stores 1 + \sum_{n=1}^oo (-1)^n*(q^{(3n^2-n)/2}+q^{(3n^2+n)/2})
	// * q: stores e^(2\pi i\tau)
	// * q_to_n: stores current power of q
	// * q_to_diff: stores q^diff, where diff stores the difference 
	//              between two non-trivial terms
	bigcomplex sum( 1 ), local_tau( tau ), q, q_to_n, q_to_diff;
	bigcomplex factor( 1 ), result, tau_inverse, square_root;
	
	// if is_fundamental = false translate tau to fundamental area
	// and compute appropriate factor
	if( is_fundamental == false )
	{
		//
		// ensure that |Re( tau )| <= 1/2 and |tau| >= 1
		// otherwise shift tau to this stripe
		//
		bigint real_part_as_bigint;

		while( true )
		{
			// assign to real_part_as_bigint the nearest integer
			// of Re(local_tau)
			local_tau.real().bigintify( real_part_as_bigint );
			
			// if it is not zero, shift it in stripe |Re( tau )| <= 1/2
			// eta(tau+1) = zeta_24 * eta(tau)
			if( ! real_part_as_bigint.is_zero() )
			{
				subtract( local_tau, local_tau, 
						  bigcomplex( bigfloat( real_part_as_bigint ) ) );
				multiply( factor, factor, 
						  exp( 2 * bigcomplex( 0, Pi() ) * 
							   bigfloat( real_part_as_bigint, 24 ) ) );
			}
			
			// if |local_tau| < 1 invert it to put it out of unit circle
			// eta(-1/tau) = \sqrt(-i*tau) * eta(tau)
			if( abs( local_tau ) <= 0.9999 )
			{
				invert( tau_inverse, local_tau );
				multiply( square_root, bigcomplex( 0, 1 ), tau_inverse );
				sqrt( square_root, square_root );
				
				if( square_root.real() < 0 )
					square_root.negate();
				
				multiply( factor, factor, square_root );
				
				local_tau.assign( tau_inverse );
				local_tau.negate();
			}
			// otherwise we are done
			else
				break;
		}

	}

	// factor stores the prefactor of Euler sum
	multiply( factor, factor, 
			  exp( 2 * bigcomplex( 0, Pi() ) * local_tau * 
				   bigfloat( 1, 24 ) ) );

	exp( q, 2 * bigcomplex( 0, Pi() ) * local_tau );

	// to store the powers of q in use:
	// * N_minus = ( 3 * n_minus^2 - n_minus ) / 2
	// * N_plus = ( 3 * n_plus^2 - n_plus ) / 2
	// * previous_power stores the previous used power of q
	bigint N_minus( 5 ), N_plus( 7 ); 
	int n_minus, n_plus;
	n_minus = n_plus = 2;
	bigint previous_power( 2 ), diff;

	// assign the correct values: sum = 1 - q - q^2
	square( q_to_n, q );
	subtract( sum, sum, q );
	subtract( sum, sum, q_to_n );

	while( true )
	{
		if( N_minus < N_plus )
		{
	  
			subtract( diff, N_minus, previous_power );
			power( q_to_diff, q, diff );
			multiply( q_to_n, q_to_n, q_to_diff );

			if( n_minus & 1 )
				subtract( sum, sum, q_to_n );
			else
				add( sum, sum, q_to_n );

			if( q_to_n.is_approx_zero() )
			{
				multiply( result, factor, sum );

				return( result );
			}

			previous_power.assign( N_minus );
			N_minus += 3 * n_minus + 1;
			n_minus++;
		}

		else
		{
			subtract( diff, N_plus, previous_power );
			power( q_to_diff, q, diff );
			multiply( q_to_n, q_to_n, q_to_diff );
	  
			if( n_plus & 1 )
				subtract( sum, sum, q_to_n );
			else
				add( sum, sum, q_to_n );

			if( q_to_n.is_approx_zero() )
			{
				multiply( result, factor, sum );
		
				return( result );
			}

			previous_power.assign( N_plus );
			N_plus += 3 * n_plus + 2;
			n_plus++;
		}
	}
}


/********************************************************************
 *
 * Weber functions computed via the eta-function.
 * See PhD for details.
 *
 *******************************************************************/
//
// weber_f bases on the algorithm described in PhD
//
bigcomplex
weber_f( const bigcomplex & tau, bool is_fundamental )
{
	// if tau is not in the fundamental domain, we make use of the
	// relation weber_f = eta( (tau+1)/2 ) / ( zeta_48 * eta( tau ) )
	if( is_fundamental == false )
	{
		// declare the necessary variables
		bigcomplex current_tau, numerator, denominator, zeta_48, result;
		
		// compute zeta_48 = e^{2\pi i / 48}
		exp( zeta_48, bigcomplex( 0, Pi() ) * bigfloat( 1, 24 ) );

		// compute argument of numerator, numerator, and denominator
		current_tau = ( tau + 1 ) / 2;
		numerator = dedekind_eta( current_tau, false );
		denominator = dedekind_eta( tau, false );

		// compute the result and return it
		divide( result, numerator, zeta_48 );
		divide( result, result, denominator );

		return( result );
	}

	// if is_fundamental is true, use the algorithm of PhD

	// complex variables for computation of eta((tau+1)/2): 
	// * eta_1_2: stores 1 + \sum_{n=1}^oo (-1)^n*(q^{(3n^2-n)/4}+q^{(3n^2+n)/4})
	//          hence the current approximation of \eta(\tau/2)/q^{1/48}
	// * q_1_2: to store q^{1/2} = e^(\pi i\tau)
	// * q_1_2_to_n: to store q^{(3n^2+n)/4}
	// * q_to_diff: stores q^diff, where diff stores the difference 
	//              between two non-trivial terms
	bigcomplex eta_1_2( 1 ), q_1_2, q_1_2_to_n, q_to_diff, result;
	exp( q_1_2, bigcomplex( 0, Pi() ) * tau );

	// complex variables for computation of eta(tau): 
	// * eta: stores 1 + \sum_{n=1}^oo (-1)^n*(q^{(3n^2-n)/2}+q^{(3n^2+n)/2}),
	//        hence the current approximation of \eta(\tau)/q^{1/24}
	// * q_to_n: stores current power of q
	bigcomplex eta( 1 ), q_to_n;

	// to store the powers of q_1_2 in use:
	// * N_minus = ( 3 * n_minus^2 - n_minus ) / 2
	// * N_plus = ( 3 * n_plus^2 - n_plus ) / 2
	// * previous_power stores the previous used power of q
	bigint N_minus( 5 ), N_plus( 7 ); 
	int n_minus, n_plus;
	n_minus = n_plus = 2;
	bigint previous_power( 2 ), diff;

	// compute the prefactor q^{-1/48} of the result
	bigcomplex prefactor;
	exp( prefactor, bigcomplex( 0, Pi() ) * tau * bigfloat( -1, 24 ) );

	// assign the correct values: eta_1_2 = 1 + q_1_2 - q_1_2^2
	square( q_1_2_to_n, q_1_2 );
	add( eta_1_2, eta_1_2, q_1_2 );
	subtract( eta_1_2, eta_1_2, q_1_2_to_n );

	// assign the correct values: eta = 1 - q - q^2
	subtract( eta, eta, q_1_2_to_n );
	subtract( eta, eta, square( q_1_2_to_n ) );

	// eta_false is true iff q_to_n is not zero within precision
	bool eta_false = true;

	// only evaluate this loop iff q_to_n is not zero within precision
	while( eta_false )
	{
		if( N_minus < N_plus )
		{
			// compute the difference of the exponents
			subtract( diff, N_minus, previous_power );
			power( q_to_diff, q_1_2, diff );
			multiply( q_1_2_to_n, q_1_2_to_n, q_to_diff );

			square( q_to_n, q_1_2_to_n );

			// take care of the sign and update eta_1_2 and eta
			if( n_minus & 1 )
			{
				if( N_minus.bit( 0 ) )
					add( eta_1_2, eta_1_2, q_1_2_to_n );
				else
					subtract( eta_1_2, eta_1_2, q_1_2_to_n );

				subtract( eta, eta, q_to_n );
			}
			else
			{
				if( N_minus.bit( 0 ) )
					subtract( eta_1_2, eta_1_2, q_1_2_to_n );
				else
					add( eta_1_2, eta_1_2, q_1_2_to_n );

				add( eta, eta, q_to_n );
			}

			// check if q_to_n is zero within precision
			if( q_to_n.is_approx_zero() )
				eta_false = false;

			// adapt the exponent N_minus and the counter n_minus
			previous_power.assign( N_minus );
			N_minus += 3 * n_minus + 1;
			n_minus++;
		}
		else
		{
			// compute the difference of the exponents
			subtract( diff, N_plus, previous_power );
			power( q_to_diff, q_1_2, diff );
			multiply( q_1_2_to_n, q_1_2_to_n, q_to_diff );

			square( q_to_n, q_1_2_to_n );

			// take care of the sign and update eta_1_2 and eta
			if( n_plus & 1 )
			{
				if( N_plus.bit( 0 ) )
					add( eta_1_2, eta_1_2, q_1_2_to_n );
				else
					subtract( eta_1_2, eta_1_2, q_1_2_to_n );

				subtract( eta, eta, q_to_n );
			}
			else
			{
				if( N_plus.bit( 0 ) )
					subtract( eta_1_2, eta_1_2, q_1_2_to_n );
				else
					add( eta_1_2, eta_1_2, q_1_2_to_n );

				add( eta, eta, q_to_n );
			}

			// check if q_to_n is zero within precision
			if( q_to_n.is_approx_zero() )
				eta_false = false;

			// adapt the exponent N_plus and the counter n_plus
			previous_power.assign( N_plus );
			N_plus += 3 * n_plus + 2;
			n_plus++;
		}
	}

	// evaluate this loop until q_to_n is zero within precision
	while( true )
	{
		if( N_minus < N_plus )
		{
			// compute the difference of the exponents
			subtract( diff, N_minus, previous_power );
			power( q_to_diff, q_1_2, diff );
			multiply( q_1_2_to_n, q_1_2_to_n, q_to_diff );

			square( q_to_n, q_1_2_to_n );

			// take care of the sign and update eta_1_2
			if( n_minus & 1 )
			{
				if( N_minus.bit( 0 ) )
					add( eta_1_2, eta_1_2, q_1_2_to_n );
				else
					subtract( eta_1_2, eta_1_2, q_1_2_to_n );
			}
			else
			{
				if( N_minus.bit( 0 ) )
					subtract( eta_1_2, eta_1_2, q_1_2_to_n );
				else
					add( eta_1_2, eta_1_2, q_1_2_to_n );
			}

			// if q_1_2_to_n is zero within precision return the result:
			// q^{-1/48} * eta_1_2 / eta = prefactor * eta_1_2 / eta
			if( q_1_2_to_n.is_approx_zero() )
			{
				multiply( result, prefactor, eta_1_2 );
				divide( result, result, eta );
		
				return( result );
			}

			previous_power.assign( N_minus );
			N_minus += 3 * n_minus + 1;
			n_minus++;
		}
		else
		{
			// compute the difference of the exponents
			subtract( diff, N_plus, previous_power );
			power( q_to_diff, q_1_2, diff );
			multiply( q_1_2_to_n, q_1_2_to_n, q_to_diff );

			square( q_to_n, q_1_2_to_n );

			// take care of the sign and update eta_1_2
			if( n_plus & 1 )
			{
				if( N_plus.bit( 0 ) )
					add( eta_1_2, eta_1_2, q_1_2_to_n );
				else
					subtract( eta_1_2, eta_1_2, q_1_2_to_n );
			}
			else
			{
				if( N_plus.bit( 0 ) )
					subtract( eta_1_2, eta_1_2, q_1_2_to_n );
				else
					add( eta_1_2, eta_1_2, q_1_2_to_n );
			}

			// if q_1_2_to_n is zero within precision return the result:
			// q^{-1/48} * eta_1_2 / eta = prefactor * eta_1_2 / eta
			if( q_1_2_to_n.is_approx_zero() )
			{
				multiply( result, prefactor, eta_1_2 );
				divide( result, result, eta );
		
				return( result );
			}

			previous_power.assign( N_plus );
			N_plus += 3 * n_plus + 2;
			n_plus++;
		}
	}
}

//
// weber_f1 bases on the algorithm in PhD
//
bigcomplex
weber_f1( const bigcomplex & tau, bool is_fundamental )
{
	// if tau is not in the fundamental domain, we make use of the
	// relation weber_f1 = eta( tau/2 ) / eta( tau )
	if( is_fundamental == false )
	{
		// declare the necessary variables
		bigcomplex current_tau, numerator, denominator, result;
		
		// compute argument of numerator, numerator, and denominator
		divide( current_tau, tau, 2 );
		numerator = dedekind_eta( current_tau, false );
		denominator = dedekind_eta( tau, false );

		// compute the result and return it
		divide( result, numerator, denominator );

		return( result );
	}

	// if is_fundamental is true, use the algorithm of PhD

	// complex variables for computation of eta(tau/2): 
	// * eta_1_2: stores 1 + \sum_{n=1}^oo (-1)^n*(q^{(3n^2-n)/4}+q^{(3n^2+n)/4})
	//          hence the current approximation of \eta(\tau/2)/q^{1/48}
	// * q_1_2: to store q^{1/2} = e^(\pi i\tau)
	// * q_1_2_to_n: to store q^{(3n^2+n)/4}
	// * q_to_diff: stores q^diff, where diff stores the difference 
	//              between two non-trivial terms
	bigcomplex eta_1_2( 1 ), q_1_2, q_1_2_to_n, q_to_diff, result;
	exp( q_1_2, bigcomplex( 0, Pi() ) * tau );

	// complex variables for computation of eta(tau): 
	// * eta: stores 1 + \sum_{n=1}^oo (-1)^n*(q^{(3n^2-n)/2}+q^{(3n^2+n)/2}),
	//        hence the current approximation of \eta(\tau)/q^{1/24}
	// * q_to_n: stores current power of q
	bigcomplex eta( 1 ), q_to_n;

	// to store the powers of q_1_2 in use:
	// * N_minus = ( 3 * n_minus^2 - n_minus ) / 2
	// * N_plus = ( 3 * n_plus^2 - n_plus ) / 2
	// * previous_power stores the previous used power of q
	bigint N_minus( 5 ), N_plus( 7 ); 
	int n_minus, n_plus;
	n_minus = n_plus = 2;
	bigint previous_power( 2 ), diff;

	// compute the prefactor q^{-1/48} of the result
	bigcomplex prefactor;
	exp( prefactor, bigcomplex( 0, Pi() ) * tau * bigfloat( -1, 24 ) );

	// assign the correct values: eta_1_2 = 1 - q_1_2 - q_1_2^2
	square( q_1_2_to_n, q_1_2 );
	subtract( eta_1_2, eta_1_2, q_1_2 );
	subtract( eta_1_2, eta_1_2, q_1_2_to_n );

	// assign the correct values: eta = 1 - q - q^2
	subtract( eta, eta, q_1_2_to_n );
	subtract( eta, eta, square( q_1_2_to_n ) );

	// eta_false is true iff q_to_n is not zero within precision
	bool eta_false = true;

	// only evaluate this loop iff q_to_n is not zero within precision
	while( eta_false )
	{
		if( N_minus < N_plus )
		{
			// compute the difference of the exponents
			subtract( diff, N_minus, previous_power );
			power( q_to_diff, q_1_2, diff );
			multiply( q_1_2_to_n, q_1_2_to_n, q_to_diff );

			square( q_to_n, q_1_2_to_n );

			// take care of the sign and update eta_1_2 and eta
			if( n_minus & 1 )
			{
				subtract( eta_1_2, eta_1_2, q_1_2_to_n );
				subtract( eta, eta, q_to_n );
			}
			else
			{
				add( eta_1_2, eta_1_2, q_1_2_to_n );
				add( eta, eta, q_to_n );
			}

			// check if q_to_n is zero within precision
			if( q_to_n.is_approx_zero() )
				eta_false = false;

			// adapt the exponent N_minus and the counter n_minus
			previous_power.assign( N_minus );
			N_minus += 3 * n_minus + 1;
			n_minus++;
		}
		else
		{
			// compute the difference of the exponents
			subtract( diff, N_plus, previous_power );
			power( q_to_diff, q_1_2, diff );
			multiply( q_1_2_to_n, q_1_2_to_n, q_to_diff );

			square( q_to_n, q_1_2_to_n );

			// take care of the sign and update eta_1_2 and eta
			if( n_plus & 1 )
			{
				subtract( eta_1_2, eta_1_2, q_1_2_to_n );
				subtract( eta, eta, q_to_n );
			}
			else
			{
				add( eta_1_2, eta_1_2, q_1_2_to_n );
				add( eta, eta, q_to_n );
			}

			// check if q_to_n is zero within precision
			if( q_to_n.is_approx_zero() )
				eta_false = false;

			// adapt the exponent N_plus and the counter n_plus
			previous_power.assign( N_plus );
			N_plus += 3 * n_plus + 2;
			n_plus++;
		}
	}

	// evaluate this loop until q_to_n is zero within precision
	while( true )
	{
		if( N_minus < N_plus )
		{
			// compute the difference of the exponents
			subtract( diff, N_minus, previous_power );
			power( q_to_diff, q_1_2, diff );
			multiply( q_1_2_to_n, q_1_2_to_n, q_to_diff );

			square( q_to_n, q_1_2_to_n );

			// take care of the sign and update eta_1_2
			if( n_minus & 1 )
				subtract( eta_1_2, eta_1_2, q_1_2_to_n );
			else
				add( eta_1_2, eta_1_2, q_1_2_to_n );

			// if q_1_2_to_n is zero within precision return the result:
			// q^{-1/48} * eta_1_2 / eta = prefactor * eta_1_2 / eta
			if( q_1_2_to_n.is_approx_zero() )
			{
				multiply( result, prefactor, eta_1_2 );
				divide( result, result, eta );
		
				return( result );
			}

			previous_power.assign( N_minus );
			N_minus += 3 * n_minus + 1;
			n_minus++;
		}
		else
		{
			// compute the difference of the exponents
			subtract( diff, N_plus, previous_power );
			power( q_to_diff, q_1_2, diff );
			multiply( q_1_2_to_n, q_1_2_to_n, q_to_diff );

			square( q_to_n, q_1_2_to_n );

			// take care of the sign and update eta_1_2
			if( n_plus & 1 )
				subtract( eta_1_2, eta_1_2, q_1_2_to_n );
			else
				add( eta_1_2, eta_1_2, q_1_2_to_n );

			// if q_1_2_to_n is zero within precision return the result:
			// q^{-1/48} * eta_1_2 / eta = prefactor * eta_1_2 / eta
			if( q_1_2_to_n.is_approx_zero() )
			{
				multiply( result, prefactor, eta_1_2 );
				divide( result, result, eta );
		
				return( result );
			}

			previous_power.assign( N_plus );
			N_plus += 3 * n_plus + 2;
			n_plus++;
		}
	}
}

//
// weber_f2 is described in PhD
//
bigcomplex 
weber_f2( const bigcomplex & tau, bool is_fundamental )
{
	// if tau is not in the fundamental domain, we make use of the
	// relation weber_f2 = sqrt( 2 ) * eta( 2*tau ) / eta( tau )
	if( is_fundamental == false )
	{
		// declare the necessary variables
		bigcomplex current_tau, numerator, denominator, result;
		bigfloat sqrt_2;

		// compute sqrt_2 = sqrt( 2 )
		sqrt_2 = sqrt( bigfloat( 2.0 ) );

		// compute argument of numerator, numerator, and denominator
		multiply( current_tau, tau, 2 );
		numerator = dedekind_eta( current_tau, false );
		denominator = dedekind_eta( tau, false );

		// compute the result and return it
		multiply( result, sqrt_2, numerator );
		divide( result, result, denominator );

		return( result );
	}

	// otherwise use algorithm of PhD

	// complex variables for computation of eta(tau): 
	// * eta: stores 1 + \sum_{n=1}^oo (-1)^n*(q^{(3n^2-n)/2}+q^{(3n^2+n)/2}),
	//        hence the current approximation of \eta(\tau)/q^{1/24}
	// * q: stores e^(2\pi i\tau)
	// * q_to_n: stores current power of q
	// * q_to_diff: stores q^diff, where diff stores the difference 
	//              between two non-trivial terms
	bigcomplex eta( 1 ), q, q_to_n, q_to_diff, result;
	exp( q, 2 * bigcomplex( 0, Pi() ) * tau );

	// to store the powers of q in use:
	// * N_minus = ( 3 * n_minus^2 - n_minus ) / 2
	// * N_plus = ( 3 * n_plus^2 - n_plus ) / 2
	// * previous_power stores the previous used power of q
	bigint N_minus( 5 ), N_plus( 7 ); 
	int n_minus, n_plus;
	n_minus = n_plus = 2;
	bigint previous_power( 2 ), diff;

	// compute the prefactor \sqrt{2} * q^{1/24} of the result
	bigcomplex prefactor;
	exp( prefactor, 2 * bigcomplex( 0, Pi() ) * tau * bigfloat( 1, 24 ) );
	multiply( prefactor, prefactor, sqrt( bigfloat( 2 ) ) );

	// assign the correct values: eta = 1 - q - q^2
	square( q_to_n, q );
	subtract( eta, eta, q );
	subtract( eta, eta, q_to_n );

	// complex variables for computation of eta(2tau): 
	// * eta_2: stores 1 + \sum_{n=1}^oo (-1)^n*(q^{(3n^2-n)}+q^{(3n^2+n)})
	//          hence the current approximation of \eta(2\tau)/q^{1/12}
	// * q_to_n_square: to store q_to_n^2
	bigcomplex eta_2( 1 ), q_to_n_square;

	// assign the correct values: eta_2 = 1 - q^2 - q^4
	subtract( eta_2, eta_2, q_to_n );
	subtract( eta_2, eta_2, square( q_to_n ) );

	// eta_2_false is true iff q_to_n_square is not zero within precision
	bool eta_2_false = true;

	// only evaluate this loop iff q_to_n_square is not zero within precision
	while( eta_2_false )
	{
		if( N_minus < N_plus )
		{
			// compute the difference of the exponents
			subtract( diff, N_minus, previous_power );
			power( q_to_diff, q, diff );
			multiply( q_to_n, q_to_n, q_to_diff );

			square( q_to_n_square, q_to_n );

			// take care of the sign and update eta and eta_2
			if( n_minus & 1 )
			{
				subtract( eta, eta, q_to_n );
				subtract( eta_2, eta_2, q_to_n_square );
			}
			else
			{
				add( eta, eta, q_to_n );
				add( eta_2, eta_2, q_to_n_square );
			}

			// check if q_to_n_square is zero within precision
			if( q_to_n_square.is_approx_zero() )
				eta_2_false = false;

			// adapt the exponent N_minus and the counter n_minus
			previous_power.assign( N_minus );
			N_minus += 3 * n_minus + 1;
			n_minus++;
		}
		else
		{
			// compute the difference of the exponents
			subtract( diff, N_plus, previous_power );
			power( q_to_diff, q, diff );
			multiply( q_to_n, q_to_n, q_to_diff );
	  
			square( q_to_n_square, q_to_n );

			// take care of the sign and update eta and eta_2
			if( n_plus & 1 )
			{
				subtract( eta, eta, q_to_n );
				subtract( eta_2, eta_2, q_to_n_square );
			}
			else
			{
				add( eta, eta, q_to_n );
				add( eta_2, eta_2, q_to_n_square );
			}

			// check if q_to_n_square is zero within precision
			if( q_to_n_square.is_approx_zero() )
				eta_2_false = false;

			// adapt the exponent N_plus and the counter n_plus
			previous_power.assign( N_plus );
			N_plus += 3 * n_plus + 2;
			n_plus++;
		}
	}

	// evaluate this loop until q_to_n is zero within precision
	while( true )
	{
		if( N_minus < N_plus )
		{
			subtract( diff, N_minus, previous_power );
			power( q_to_diff, q, diff );
			multiply( q_to_n, q_to_n, q_to_diff );

			if( n_minus & 1 )
				subtract( eta, eta, q_to_n );
			else
				add( eta, eta, q_to_n );

			// if q_to_n is zero within precision return the result:
			// \sqrt(2) * q^{1/24} * eta_2 / eta = prefactor * eta_2 / eta
			if( q_to_n.is_approx_zero() )
			{
				multiply( result, prefactor, eta_2 );
				divide( result, result, eta );
		
				return( result );
			}

			previous_power.assign( N_minus );
			N_minus += 3 * n_minus + 1;
			n_minus++;
		}

		else
		{
			subtract( diff, N_plus, previous_power );
			power( q_to_diff, q, diff );
			multiply( q_to_n, q_to_n, q_to_diff );
	  
			if( n_plus & 1 )
				subtract( eta, eta, q_to_n );
			else
				add( eta, eta, q_to_n );

			// if q_to_n is zero within precision return the result:
			// \sqrt(2) * q^{1/24} * eta_2 / eta = prefactor * eta_2 / eta
			if( q_to_n.is_approx_zero() )
			{
				multiply( result, prefactor, eta_2 );
				divide( result, result, eta );
		
				return( result );
			}

			previous_power.assign( N_plus );
			N_plus += 3 * n_plus + 2;
			n_plus++;
		}
	}
}


/********************************************************************
 *
 * Compute \gamma_2 and the modular function j using weber_f2.
 *
 *******************************************************************/
//
// The function gamma2 bases on the direct computation of 
// weber_f2 via the eta-function.
//
bigcomplex 
gamma_2( const bigcomplex & tau, bool is_fundamental )
{
	// assign:
	// f2_to_8 = weber_f2( tau ) ^ 8
	// assign gamma2 = (f2^24 + 16 ) / f2^8 to result
	bigcomplex f2_to_8, result, local_tau( tau ), tau_inverse, factor( 1 );

	// if is_fundamental = false translate tau to fundamental area
	// and compute appropriate factor
	if( is_fundamental == false )
	{
		//
		// ensure that |Re( tau )| <= 1/2 and |tau| >= 1
		// otherwise shift tau to this stripe
		//
		bigint real_part_as_bigint;
		
 		while( true )
		{
			// assign to real_part_as_bigint the nearest integer
			// of Re(local_tau)
			local_tau.real().bigintify( real_part_as_bigint );
			
			// if it is not zero, shift it in stripe |Re( tau )| <= 1/2
			// gamma_2(tau+1) = 1/zeta_3 * gamma(tau)
			if( ! real_part_as_bigint.is_zero() )
			{
				subtract( local_tau, local_tau, 
						  bigcomplex( bigfloat( real_part_as_bigint ) ) );
				multiply( factor, factor, 
						  exp( -2 * bigcomplex( 0, Pi() ) * 
							   bigfloat( real_part_as_bigint, 3 ) ) );
			}
			
			// if |local_tau| < 1 invert it to put it out of unit circle
			// gamma_2(-1/tau) = gamma(tau)
			if( abs( local_tau ) <= 0.9999 )
			{
				invert( tau_inverse, local_tau );
				local_tau.assign( tau_inverse );
				local_tau.negate();
			}
			// otherwise we are done
			else
				break;
		}
	}

	power( f2_to_8, weber_f2( local_tau, true ), 8 );
	power( result, f2_to_8, 3 );
	add( result, result, bigint( 16 ) );
	divide( result, result, f2_to_8 );

	multiply( result, result, factor );

	return( result );
}


bigcomplex 
modular_j( const bigcomplex & tau, bool is_fundamental )
{
	// assign:
	// f2_to_24 = weber_f2( tau ) ^ 24
	// j = (f2^24 + 16 )^3 / f2^24 to result
	//
	bigcomplex f2_to_24, result, local_tau( tau ), tau_inverse;

	if( is_fundamental == false )
	{
		//
		// ensure that |Re( tau )| <= 1/2 and |tau| >= 1
		// otherwise shift tau to this stripe
		//
		bigint real_part_as_bigint;

		while( true )
		{
			// assign to real_part_as_bigint the nearest integer
			// of Re(local_tau)
			local_tau.real().bigintify( real_part_as_bigint );
			
			// if it is not zero, shift it in stripe |Re( tau )| <= 1/2
			// j(tau+1) = j(tau)
			if( ! real_part_as_bigint.is_zero() )
			{
				subtract( local_tau, local_tau, 
						  bigcomplex( bigfloat( real_part_as_bigint ) ) );
			}
			
			// if |local_tau| < 1 invert it to put it out of unit circle
			// j(-1/tau) = j(tau)
			if( abs( local_tau ) <= 0.9999 )
			{
				invert( tau_inverse, local_tau );
				local_tau.assign( tau_inverse );
				local_tau.negate();
			}
			// otherwise we are done
			else
				break;
		}

	}
	
	power( f2_to_24, weber_f2( local_tau, true ), 24 );
	add( result, f2_to_24, bigint( 16 ) );
	power( result, result, 3 );
	divide( result, result, f2_to_24 );

	return( result );
}


#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif
