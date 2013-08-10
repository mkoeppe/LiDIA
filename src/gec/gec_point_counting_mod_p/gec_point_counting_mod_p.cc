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
//  File    : gec_point_counting_mod_p.cc
//  Author  : Harald Baier (HB)
//	Changes	: See CVS log
//
//==============================================================================================

#ifdef HAVE_CONFIG_H
# include       "config.h"
#endif
# include        "LiDIA/gec_point_counting_mod_p.h"
#ifdef VERBOSE
# include        <cassert>
#endif


#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif


//
// constructor and destructor
//
gec_point_counting_mod_p::gec_point_counting_mod_p() : gec()
{
	PC_FLAG        = true;
	tries          = 0;
	both_abort     = 0;
	no_abort       = 0;
	is_initialized = false;
}

gec_point_counting_mod_p::~gec_point_counting_mod_p()
{
}

//
// Accessors
//
lidia_size_t
gec_point_counting_mod_p::get_number_of_initializations() const
{
	return( tries );
}

lidia_size_t
gec_point_counting_mod_p::get_number_of_both_aborts() const
{
	return( both_abort );
}

lidia_size_t
gec_point_counting_mod_p::get_number_of_no_aborts() const
{
	return( no_abort );
}


//
// High level functions
//
void 
gec_point_counting_mod_p::
get_twist_coeff( gf_element & new_a4, gf_element & new_a6 )
{
	gf_element d( a4 );
	
	do
	{
		d.randomize();
	}
	while (jacobi(d.polynomial_rep().const_term(), d.characteristic()) != -1);
	
	new_a4 = a4*d*d;
	new_a6 = a6*d*d*d;
}

void
gec_point_counting_mod_p::
generate()
{
	bigint tmp(1);

	// Set field cardinality q to a random prime
	if( q == 0 )
	{
		shift_left( tmp, tmp, lower_bound_bitlength_r - 2 );
		q.randomize( tmp );
      
		shift_left( tmp, tmp, 1 );
		multiply( tmp, tmp, upper_bound_k );

		add( q, q, tmp );
		q = next_prime( q - 1 );
	}

	// check whether the bounds are set appropriate
	check_the_settings();

	// some variables for early abort strategy
	bigint found_co_factor, found_co_factor_tw, order;

	//
	// Initialize counters
	//
	tries      = 1;
	both_abort = 0;
	no_abort   = 0;

	//
	// computation control
	//
	bool good_curve_Etw;
	bool good_curve_E;
	bool stop_computation;

	// assign the prime field to theField of type galois_field
	galois_field theField( q );
	bigmod::set_modulus( q );

	// coordinates of twisted curve
	gf_element a4tw( theField ), a6tw( theField );

	// initialize computation
	elliptic_curve<gf_element> e;
	trace_list tl;
	trace_mod tm;
	udigit ll;
	eco_prime ep, eptw;
	bigint g;

	if( VERBOSE )
	{
		ep.set_info_mode(1);
		trace_list::set_info_mode(1);

		std::cout << "\nCharacteristic : q = " << q;
		std::cout << " (" << q.bit_length() << " bits)" << std::endl;
		std::cout << "Upper bound for k is " << upper_bound_k;
		std::cout << " (" << upper_bound_k.bit_length() << " bits)"
				  << std::endl;
		std::cout << std::endl;
	}
	else
	{	  
		ep.set_info_mode(0);
		trace_list::set_info_mode(0);
	}

	//
	// Choose curves until strong one is found.
	// Use early abort strategy on curve and its twist.
	//
	do 
	{
		// Choose next curve
		a4.assign_zero( theField );
		a6.assign_zero( theField );
		a4.randomize(); 
		a6.randomize(); 

		// Determine coefficients of twist.
		get_twist_coeff( a4tw, a6tw );

		if( VERBOSE )
		{
			std::cout << std::endl << "Current E : a4 = " << a4 << std::endl;
			std::cout << "a6 = " << a6 << std::endl << std::endl;
			std::cout << "Current E_twist : a4tw = " << a4tw << std::endl;
			std::cout << "a6tw = "
					  << a6tw << std::endl << std::endl;
		}

		// initialize for new curve
		tries++; 
		found_co_factor = 1;
		found_co_factor_tw = 1;

		good_curve_E   = true;
		good_curve_Etw = true;
		stop_computation = false;

		e.set_coefficients( a4, a6 );
		ep.set_curve( a4, a6 );
		eptw.set_curve( a4tw, a6tw );

      // Test supersingularity: Remark that for a large prime
      // E is supersingular iff Etw is so
      //
		if( VERBOSE )
			std::cout << "\nTesting E for supersingularity.... " << std::flush;

		if( ep.is_supersingular( order ) )
		{
			if( VERBOSE )
				std::cout << "YES" << std::endl;

			continue;
		}
		else if( VERBOSE )
			std::cout << "NO" << std::endl;


      // Test for isogeny to curve with j-invariant 0 or 1728
      //
		if( VERBOSE )
			std::cout << "Testing whether E is Fq-isogenous to curve with "
				"j-invariant 0 or 1728 ... " << std::flush;

		if( ep.check_j_0_1728( order ) )
		{
			if( VERBOSE )
				std::cout << "YES" << std::endl;
			continue;
		}
		else if( VERBOSE )
			std::cout << "NO" << std::endl;


      // Test twist for isogeny to curve with j-invariant 0 or 1728
      //
		if( VERBOSE )
			std::cout << "Testing whether twist(E) is Fq-isogenous to curve "
				"with j-invariant 0 or 1728 ... " << std::flush;

		if( eptw.check_j_0_1728( order ) )
		{
			if( VERBOSE )
				std::cout << "YES" << std::endl;
			continue;
		}
		else if( VERBOSE )
			std::cout << "NO" << std::endl;


      // Compute trace modulo 2
      //
		tl.clear();
		tl.set_curve( e );

      // trace mod 2
      //
		if( VERBOSE )
		{
			std::cout << std::endl;
			std::cout << "--------------------------------------" << std::endl;
			std::cout << "Working on prime l = 2" << std::endl;
		}

		ep.compute_mod_2_power();
		tm.set_vector( ep.get_prime(), ep.get_relation() );
		tl.append( tm );

	// There is no need to use qq as the found cofactors
	// of E and its twist are computed in the lines that follow		
//      qq = remainder( q, ep.get_prime() );
		g = gcd( q + 1 - (ep.get_relation())[0], ep.get_prime() );
		if (! g.is_one())
			multiply(found_co_factor, found_co_factor, g);
		g = gcd(q + 1 + (ep.get_relation())[0], ep.get_prime());
		if (! g.is_one())
			multiply(found_co_factor_tw, found_co_factor_tw, g);


      // for each prime ll compute trace mod ll
      // and check for divisor of group order
      //
		ll = 3;

		do 
		{
			if( found_co_factor > upper_bound_k && 
				found_co_factor_tw > upper_bound_k )
			{
				both_abort ++;
				good_curve_E = good_curve_Etw = false;
				stop_computation = true;
			}
			else 
			{
				if( VERBOSE )
				{
					std::cout << std::endl;
					std::cout << "--------------------------------------" << std::endl;
					std::cout << "Working on prime l = " << ll <<std::endl;
				}

				if( ep.set_prime(ll) )
				{
					ep.compute_splitting_type();
		   
					if ( ep.is_elkies() )
					{
						ep.compute_trace_elkies();
	
						// see remark above
//		     qq = remainder(q, ep.get_prime());
						g = gcd(q + 1 - (long) (ep.get_relation())[0], ep.get_prime());
						if (! g.is_one())
							multiply(found_co_factor, found_co_factor, g);
						g = gcd(q + 1 + (long) (ep.get_relation())[0], ep.get_prime());
						if (! g.is_one())
							multiply(found_co_factor_tw, found_co_factor_tw, g);
					}
					else
						ep.compute_trace_atkin();
		
					tm.set_vector( ep.get_prime(), ep.get_relation() );
					stop_computation = tl.append( tm );
				}
				else 
				{
					if( ll >= ep.MAX_MEQ_PRIME )
						lidia_error_handler("gec_point_counting_mod_p",
											"Not enough modular equations");
					else
						lidia_warning_handler("gec_point_counting_mod_p",
											  "Modular equation could not be"
											  " found, ignored ...");
				}
		
				ll = next_prime( ll + 1 );
			}
		}
		while( ! stop_computation );

      // Early abort
      //
		if( ! good_curve_E && ! good_curve_Etw ) 
			continue;
      
      // Determine group order of E
      //     
		if( VERBOSE )
			std::cout<< std::endl;
      
		order = tl.bg_search_for_order();
      
		no_abort++;

      // Test for strong curve E
      //
		if( found_co_factor <= upper_bound_k )
			if( is_cryptographically_strong( order ) )
				good_curve_E = true;
			else
				good_curve_E = false;
		else 
		{
		    if( VERBOSE )
			std::cout << "found_co_factor too big." << std::endl;
		    good_curve_E = false;
		}

		// Test for strong twist Etw
		//
		if( ! good_curve_E ) 
		{
			if (found_co_factor_tw <= upper_bound_k) 
			{ 
				order = q + 1 + ( q + 1 - order );
				if( is_cryptographically_strong( order ) )
				{
					good_curve_Etw = true;
					a4 = a4tw;
					a6 = a6tw;
				}
				else
					good_curve_Etw = false;
			}
			else
			{
				if( VERBOSE )
				    std::cout << "found_co_factor_tw too big." << std::endl;

				good_curve_Etw = false;
			}
		} 
	}
	while( ! good_curve_E && ! good_curve_Etw );
  
  
  // Verification of the group order
  //
	e.set_coefficients( a4, a6 );

	if ( e.probabilistic_test_of_group_order( order ) )
	{
		// search point G of order r
		//      
		point<gf_element> P;
		do
		{
			P = e.random_point();
			multiply(G, k, P);	
		}
		while(G.is_zero());

		if( VERBOSE )
		{
			std::cout << std::endl <<"FOUND GOOD CURVE : " << std::endl;
			std::cout << "Characteristic : " << q << " (" << q.bit_length();
			std::cout << " bits )" << std::endl;
			std::cout << "r = " << r << std::endl;
			std::cout << "k = " << k << std::endl;
			std::cout << "k*r = " << order << std::endl;
			std::cout << "Coefficient a4 : " << a4 << std::endl;
			std::cout << "Coefficient a6 : " << a6 << std::endl;
			std::cout << "Point G of order r : " << G << std::endl;
		}

		is_initialized = true; 
	}
	else
	{
		std::cout<<"\n\nERROR: Probabilistic correctness test rejects candidate!!\n";
		exit(1);
	}
}


#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif
