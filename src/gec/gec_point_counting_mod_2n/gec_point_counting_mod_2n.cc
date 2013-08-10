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
//  File    : gec_point_counting_mod_2n.cc
//  Author  : Harald Baier (HB)
//	Changes	: See CVS log
//
//==============================================================================================

#ifdef HAVE_CONFIG_H
# include       "config.h"
#endif
# include        "LiDIA/gec_point_counting_mod_2n.h"
#ifdef VERBOSE
# include        <cassert>
#endif


#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif


//
// constructor and destructor
//
gec_point_counting_mod_2n::
gec_point_counting_mod_2n() :
gec()
{
  PC_FLAG        = true;
  tries          = 0;
  both_abort     = 0;
  no_abort       = 0;
  degree         = 0;
  
  is_initialized = false;
}

gec_point_counting_mod_2n::
~gec_point_counting_mod_2n()
{
}


//
// Mutators
//
void
gec_point_counting_mod_2n::
set_degree( const unsigned int n )
{
   if( n > 0 && n < upper_bound_of_degree )
   {
       //
       // Check whether bounds for r are appropriate set.
       // If no bound is set, take an upper bound of 4 for k.
       if( lower_bound_bitlength_r == 0 )
		   lower_bound_bitlength_r = n - 3;

       if( lower_bound_bitlength_r > n - 2 )
		   lidia_error_handler("gec_point_counting_mod_2n", "set_degree( const unsigned int n ): lower_bound_bitlength_r > n - 2");

       if( according_to_BSI && 
		   default_lower_bound_bitlength_r > n - 2 )
	   lidia_error_handler("gec_point_counting_mod_2n", "set_degree( const unsigned int n ): default_lower_bound_bitlength_r > n - 2");    

       degree = n;

       q = 1;
       shift_left( q, q, degree);

       F_q.assign( galois_field( bigint( 2 ), degree ) );
   }
   else
       lidia_error_handler("gec_point_counting_mod_2n", "set_degree( const unsigned int ): degree is either <= 0 or too large");
  
}


//
// High level functions
//
void
gec_point_counting_mod_2n::
generate()
{
  check_the_settings();

  bigint found_co_factor, found_co_factor_tw, order;

  //
  // Set field cardinality and bounds if not already done
  //
  if( q == 0 )
  {
      degree = lower_bound_bitlength_r + upper_bound_k.bit_length() + 1;

      q = 1;
      shift_left( q, q, degree);
      F_q.assign( galois_field( bigint( 2 ), degree ) );
  }
  else if( q != 0 && degree == 0 )
  {
     lidia_size_t bit_length_q = q.bit_length();

     for( int i = 0; i < bit_length_q; i++ )
	  if( q.bit( i ) == 1 )
	      lidia_error_handler("gec_point_counting_mod_2n","generate(): q is not a 2 power");

     degree = bit_length_q - 1;
  }
  
  if( VERBOSE )
  {
      std::cout << "lower_bound_bitlength_r = " << lower_bound_bitlength_r << std::endl;
      std::cout << "upper_bound_k = " << upper_bound_k << std::endl;
      std::cout << "degree = " << degree << std::endl;
  }

  //
  // Initialize counters
  //
  tries      = 1;
  both_abort = 0; 
  no_abort   = 0;

  //
  // computation control
  //
  bool good_curve, stop_computation;
  bool no_mod_equations;

  //
  // Initialize the GF2N arithmetic with the sparse polynomial
  // from the LiDIA GF2N.database
  //
  gf2n_init( degree );

  if( degree > 300 )
    {
      std::cout<<"\n\nInitializing table for solving quadratic equations. This";
      std::cout<<"\nmay take some time ... "<< std::flush;
      gf2n::initialize_table_for_solve_quadratic();
      std::cout<<" done.\n"<< std::flush;
    }
  
  //
  // initialize computation
  //
  gf2n a2_gf2n, a6_gf2n;

  a2_gf2n.assign_zero();
  
  elliptic_curve<gf_element> e;
  trace_list tl;
  trace_mod tm;
  udigit ll;
  eco_gf2n ep;
  bigint g;

  ep.set_info_mode( VERBOSE );
  ep.set_schoof_bound(11);
  ep.set_atkin_bound(70);
  ep.set_power_bound(100);
  ep.set_strategy(eco_gf2n::COMPUTE_SP_DEGREE_IF_ATKIN);

  //
  // Choose curves until good one is found.
  // Use early abort strategy on curve and twist.
  //
  do 
  {
      //
      // Choose next curve
      // 
      std::cout << " \n-------------------------------------------------";
      a6_gf2n.randomize(); 
      std::cout << "\nChoose random E :   a6 = " << a6_gf2n << std::endl;

      //
      // initialize for new curve
      //
      tries ++; 
      found_co_factor = 1;
      found_co_factor_tw = 1;

      good_curve       = false;
      stop_computation = false;
      no_mod_equations = false;

      e = ep.generateCurve(a6_gf2n);

      ep.set_curve(a6_gf2n);
	  ep.compute_group_order();

      tl.clear();
      tl.set_curve(e);

/*    // trace mod 2: the case ll = 2 has to be debugged
      //
	ll = 2;

      if ( VERBOSE )
      {
	  std::cout << "\n\n-------------------------------------------------------\n";
	  std::cout << "Working on prime l = " << ll << std::endl;
      }
      ep.compute_mod_2_power();
      tm.set_vector(ep.get_prime(), ep.get_relation());
      tl.append(tm);

      g = gcd(q + 1 - (ep.get_relation())[0], ep.get_prime());
      if (! g.is_one())
         multiply(found_co_factor, found_co_factor, g);
      g = gcd(q + 1 + (ep.get_relation())[0], ep.get_prime());
      if (! g.is_one())
         multiply(found_co_factor_tw, found_co_factor_tw, g);
*/


      // for each prime ll compute trace mod ll
      // and check for divisor of group order
      //

      bool enough_relations = false;
	ll = 3;
       do 
       {
	  if (found_co_factor > upper_bound_k && found_co_factor_tw > upper_bound_k)
	  {
	       found_co_factor    = 1; 
	       found_co_factor_tw = 1;
	       both_abort++;
               stop_computation = true;
	  }
          else 
	  {
	       if( VERBOSE )
	       {
		  std::cout << "\n--------------------------------------\n";
		  std::cout << "Working on prime l = " << ll;
	       }

	       if( ep.set_prime(ll) )
	       {
		   no_mod_equations = false;
		   ep.compute_splitting_type();

		   if( ep.is_elkies() )
		   {
		       ep.compute_trace_elkies();

		       g = gcd(q + 1 - (long) (ep.get_relation())[0], ep.get_prime());
		       if (! g.is_one())
			 multiply(found_co_factor, found_co_factor, g);
		       g = gcd(q + 1 + (long) (ep.get_relation())[0], ep.get_prime());
		       if (! g.is_one())
			 multiply(found_co_factor_tw, found_co_factor_tw, g);
		     }
		   else
		     ep.compute_trace_atkin();
      
		   tm.set_vector(ep.get_prime(), ep.get_relation());
		   enough_relations = tl.append(tm);
	       }
	       else if( ll >= 1020 )
	       {
		   std::cout << "\nNot enough modular equations. Aborting." << std::endl;
		   stop_computation = true;
	       }
	       
	       ll = next_prime( ll+1 );
	  }
       }
       while( ! enough_relations && ! stop_computation );

       // Exit, if there are no modular equations
       if (no_mod_equations)
       {
	   std::cout << "\nNot enough modular equations." << std::endl;
	   std::cout << "\nExiting program." << std::endl;
	   exit(1);
       }

       //
       // Early abort
       //
       if( stop_computation )   // found factor
           continue;

       //
       // Determine group order
       //
       order = tl.bg_search_for_order();

       // Test for good curve
       //
       if( found_co_factor <= upper_bound_k )
       {
	   if( is_cryptographically_strong( order ) )
	       good_curve = true;
       }
       else 
	   std::cout<<"\nearly abort for E"<< std::flush;
       
       // Test for good twist
       //
       if( ! good_curve ) 
       {
	   if( found_co_factor_tw <= upper_bound_k ) 
	   {
	       order = q + 1 + ( q + 1 - order );
	   
	       if( is_cryptographically_strong( order ) )
	       {
		   good_curve = true;
	
		   // as a2_gf2n = 0, a representative of the twist class
		   // is given by a2_gf2n = 1 in case of odd degree
		   if( degree % 2 )
			   a2_gf2n.assign_one();
		   else
		   {
			   do 
			   {
				   a2_gf2n.randomize();
			   }
			   while( a2_gf2n.trace() != 1 );
		   }
	       }
	   }
	   else 
	       std::cout<<"\nearly abort for E^tw"<< std::flush;
       }
      
      no_abort++;
  }
  while( ! good_curve );

  gf_element One( F_q );
  gf_element Zero( F_q );
  One.assign_one();
  Zero.assign_zero();
  a2 = eco_gf2n::convertToGFElement( a2_gf2n );
  a6 = eco_gf2n::convertToGFElement( a6_gf2n );
  
  e.set_coefficients(One, a2, Zero, Zero, a6);

  //
  // Probabilistic correctness proof: Choose five times a random
  // point P on e and check order*P = O.
  //
  point<gf_element> P(e), H(e);

  int i;
  for (i = 0; i < 5; i++)
  {  
      P = e.random_point();
      multiply(H, order, P);
      
      if( ! H.is_zero() )
      {
	  std::cout<<"\nProbabilistic correctness rejects candidate !!"<< std::flush;
	  std::cout<<"\n\n"; exit(1);
      }
  }
  
  do
  {
      P = e.random_point();
      
      multiply( H, k, P);
  }
  while( H.is_zero() );

  if( VERBOSE )
  {
      std::cout << "q = 2^"<< degree << " = " << q << std::endl;
      std::cout << "Coefficient a2 : " << a2 << std::endl;
      std::cout << "Coefficient a6 : " << a6 << std::endl;
      std::cout << "E: " << e << std::endl;
      std::cout << "Group order is " << order << std::endl;
      std::cout << "             = " << k << " * " << r << std::endl;
  }
}


#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif
