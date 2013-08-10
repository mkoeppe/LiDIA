//==============================================================================================
//
//	This file is part of LiDIA --- a library for computational number theory
//
//	Copyright (c) 1994--2001 the LiDIA Group.  All rights reserved.
//
//	See http://www.informatik.tu-darmstadt.de/TI/LiDIA/
//
//----------------------------------------------------------------------------------------------
//
//	$Id$
//
//	Author	: Volker Mueller (VM), Markus Maurer (MM)
//	Changes	: See CVS log
//
//==============================================================================================

#ifdef HAVE_CONFIG_H
# include	"config.h"
#endif
#include	"LiDIA/eco_prime.h"
#include	"LiDIA/wep_rat_function.h"



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



// these functions search for the eigenvalue of Frobenius after
// Elkies polynomial is known. Several search algorithms are supported.
// returns true <=> something has been found

//-----------------------------------------------------------------
//  alpha    : eigenvalue to be computed
//  fC       : Elkies polynomial 
//  lpower   : compute eigenvalue mod l^lpower
//  prev_ev  : eigenvalue mod l^(lpower-1)
//
//  test_all : controls the search for alpha in 
//             Phi (X, Y) == alpha * (X, Y)
//
//             0 --> stop searching after first match
//             1 --> test all possible values for alpha
//                   or stop after second match
//
//

bool eco_prime::find_eigenvalue (lidia_size_t & alpha,
                                 const ff_pol & fC,
                                 lidia_size_t prev_ev,
                                 lidia_size_t lpower,
                                 bool test_all)
{
  lidia_size_t *klist; // list of possible eigenvalues
  lidia_size_t lpow = 1;
  bool  rc =  false;
  
  
  //****** computation for l == 3 ******
  
  if (l == 3 && lpower == 1)
    // deg(fC) == (l-1)/2 == 1 --> computation in Fq
    return eigenv_mod3 (alpha, fC);
  
  
  //****** computation for l > 3 or lpower > 1 ******
  
  if (sp_type == ONE_SUBGRP_INV || sp_type == ALL_SUBGRPS_INV) {
    if (lpower == 1) {
      // eigenvalue^2 = p^n = q mod l
      klist = new lidia_size_t [2];
      klist[0] = 1;
      
      if (jacobi (static_cast<udigit>(eco_prime::q), static_cast<udigit>(this->l))) {
	klist[1] = (lidia_size_t)
	  sqrt_mod(static_cast<udigit>(eco_prime::q),
		   static_cast<udigit>(this->l));
	if (klist[1] > static_cast<int>((this->l-1)/2))
	  klist[1] = this->l - klist[1];
	
	if (test_x_coordinate()) {
	  alpha = klist[1];
	  if (test_x_coordinate_uniquely())
	    sign_dewaghe(alpha, fC);
	  rc = true;
	}
	else {
	  char ev = ev_strategy & 0x0f;
	  
	  if (ev == eco_prime::EV_RATIONAL_FUNCTION ||
	      ev == eco_prime::EV_RATIONAL_FUNCTION_TABLE)
	    rc = schoofpart_rat_function(alpha, fC, klist, test_all);
	  else
	    rc = schoofpart(alpha, fC, klist, test_all);
	  
	  if (alpha < 0)
	    alpha += l;
	}
	
	delete[] klist;
	return rc;
      }
      else {
	prev_ev ++;
	lidia_error_handler ("eco_prime::find_eigenvalue",
			     "p^n is not a square mod l");
	return false;
      }
    }
    else {
      lidia_error_handler ("eco_prime::find_eigenvalue()",
			   "Sorry, eigenvalue modulo powers of l not"
			   " yet implemented");
      return  false;
    }
  }
  
  //------------------------------------------------------------
  // now we have Elkies type (11ddd)
  // first we determine list of candidates for ev.
  
  
  if (lpower == 1) 
    {
      lidia_size_t * d_list;
      
      if (sp_degree_known) 
	{
	  d_list = new lidia_size_t[2];
	  d_list[0] = 1;
	  d_list[1] = sp_degree;
	}
      else
	compute_d_list (d_list);

      // use d for computing possible candidates for ev

      compute_ev_list (d_list, klist);
      refine_ev_list(fC, klist);
      delete[] d_list;
    }
  else 
    {
      lidia_error_handler ("eco_prime::find_eigenvalue()",
			   "Sorry, eigenvalue modulo prime powers not"
			   " yet implemented.");
      return false;
    }
  
  //---------------------------------------------------------
  // now we search for the ev.
  
  if (lpower == 1) {
    char ev = ev_strategy & 0x0f;
    
#ifdef TIMING
    timer t;
    t.set_print_mode();
    t.start_timer();
#endif
    
    if (ev == eco_prime::EV_RATIONAL_FUNCTION ||
	ev == eco_prime::EV_RATIONAL_FUNCTION_TABLE)
      rc = schoofpart_rat_function(alpha, fC, klist, test_all);
    else
      if (ev == eco_prime::EV_DIVISION_POLYNOMIAL ||
	  (ev == eco_prime::EV_OPTIMAL &&
	   ((l % 4 == 1) || l <= lower_bound_for_FBG)))
	rc = schoofpart(alpha, fC, klist, test_all);
      else
	if (ev == eco_prime::EV_BABYSTEP_GIANTSTEP)
	  rc = schoofpart_BG(alpha, fC, klist);
	else
	  if (ev == eco_prime::EV_FUNNY_BABYSTEP_GIANTSTEP ||
	      (ev == eco_prime::EV_OPTIMAL && l % 4 == 3))
	    rc = schoofpart_FBG(alpha, fC, klist);
	  else
	    lidia_error_handler("eco_prime",
				"wrong ev strategy chosen");
    
#ifdef TIMING
    t.stop_timer();
    if (info)
      std::cout << "\nTotal time needed to find eigenvalue : " << t << std::flush;
#endif
    
    if (test_x_coordinate_uniquely()) {
      sign_dewaghe(alpha, fC);
      
      if (info) {
	std::cout << "\nDewaghe's method: Eigenvalue of Frobenius is ";
	std::cout << alpha << std::flush;
      }
    }
    
    delete[] klist;
    
    if (rc) {
      if (alpha < 0) {
	switch (lpower) {
	case 1  : lpow = l;
	  break;
	  
	case 2  : lpow = l * l;
	  break;
	  
	default :
	  lidia_error_handler("eco_prime", "find_eigenvalue::Can't handle prime");
	}
	alpha += lpow;
      }
    }
    return rc;
  }
  return rc;
}
  
  

//
// this is the special version for finding the eigenvalue of Frobenius
// for l==3. In this case, all computations can be done in Fq, since
// degree of Elkies polynomial is 1.
//

//-------------------------------------------------------------------
// fC = (monic) Elkies polynomial for l == 3
//

bool eco_prime::eigenv_mod3 (lidia_size_t & alpha, const ff_pol &fC)
{
	ff_element rootfC, x3axb_pm1h;
	bigint pm1h; // == (q-1)/2

	shift_right (pm1h, pn, 1);
	negate (rootfC, ff_element(fC.const_term()));

	// x3axb_pm1h = (rootfC^3 + a*rootfC + b)^((q-1)/2)

	square (x3axb_pm1h, rootfC);
	add (x3axb_pm1h, eco_prime::A, x3axb_pm1h);
	multiply (x3axb_pm1h, x3axb_pm1h, rootfC);
	add      (x3axb_pm1h, x3axb_pm1h, eco_prime::B);
	power  (x3axb_pm1h, x3axb_pm1h, pm1h);

	// test whether Y^(q-1) == +-1

	if (x3axb_pm1h.is_one())
		alpha = 1;
	else
		alpha = 2;

	return true;
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
