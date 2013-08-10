//==============================================================================================
//
//      This file is part of LiDIA --- a library for computational number theory
//
//      Copyright (c) 1994--2001 the LiDIA Group.  All rights reserved.
//
//      See http://www.informatik.tu-darmstadt.de/TI/LiDIA/
//
//----------------------------------------------------------------------------------------------
//
//      $Id$
//
//      Author  : Volker Mueller (VM)
//      Changes : See CVS log
//
//==============================================================================================



#ifdef HAVE_CONFIG_H
# include	"config.h"
#endif
#include	"LiDIA/eco_gf2n.h"
#include	"LiDIA/timer.h"



#ifdef LIDIA_NAMESPACE
namespace LiDIA
{
#endif



  void eco_gf2n::compute_powers_of_l (const ff_pol & fC_in, lidia_size_t ev,
				      lidia_size_t ev1_order,
				      lidia_size_t ev2_order,
				      lidia_size_t & ev_lpower,
				      lidia_size_t & lpower)
  {
    ff_element fth;
    ff_pol f, fC;

    lpower = l;
    ev_lpower = ev;

    if (ev2_order < ev1_order)
      {
	fth = ftau[0];
	ftau[0] = ftau[1];
	ftau[1] = fth;
	ff1::set_characteristic (l);	// compute new eigenvalue
	ff1 h (ev);
	  invert (h, h);
	  multiply (h, h, q);
	  ev = static_cast < lidia_size_t > (h.as_udigit ());
	  divisor_of_divpol_via_lercier (fC, eco_gf2n::WITH_VELU);
	  ev1_order = ev2_order;
      }
    else
        fC.assign (fC_in);

#ifdef DEBUG
    gf2n_poly_modulus fff (fC);
    gf2n_polynomial div;

    compute_psi (div, l, fff);
    assert (div.is_zero ());
#endif

#ifdef TIMING
    timer t;
    t.set_print_mode ();
    t.start_timer ();
#endif

    if (ev1_order < static_cast < int >((l - 1) / 2))
      {
	gf2n_poly_modulus fCpm (fC);
	gf2n_polynomial xq;

	Xq (xq, fCpm);
	EDF_one_factor (fC, xq, fCpm, ev1_order);

#ifdef DEBUG
	assert ((fCpm.modulus () % fC).is_zero ());
#endif
      }

    ff_element old_A6 (A6), new_A6 (inverse (ftau[1]));

    prime_square (f, fC, new_A6);

#ifdef DEBUG
    set_curve (new_A6);
    fff.build (f);
    compute_psi (div, l, fff);
    assert (!div.is_zero ());
    compute_psi (div, l * l, fff);
    assert (div.is_zero ());
#endif


    if (info)
      {
	std::cout << "\nDivisor of l^2-division polynomial with degree ";
	std::cout << f.degree () << " computed. " << std::flush;
      }

#ifdef TIMING
    t.stop_timer ();
    if (info)
      {
	std::cout << "\nComputation of divisor of l^2-division"
	  "polynomial takes time ";
	std::cout << t << std::flush;
      }
#endif

    set_curve (new_A6);
    char ev_strategy_old = ev_strategy;
    unsigned int FBG_bound_old = lower_bound_for_FBG;

    lower_bound_for_FBG = 10;
    ev_strategy = eco_gf2n::EV_FUNNY_BABYSTEP_GIANTSTEP;

    if (ev < 0)
      ev = l - ev;

    find_eigenvalue (ev_lpower, f, ev, 2, false);

    lpower = l * l;
    ev_strategy = ev_strategy_old;
    lower_bound_for_FBG = FBG_bound_old;
    set_curve (old_A6);

    if ((c.size () == 1 && c[0] != 0) || test_y_coordinate ())
      {
	if (static_cast < int >(ev_lpower % l)
	    != ev && (static_cast < int >((l * l - ev_lpower) % l) == ev))
	  ev_lpower = l * l - ev_lpower;

	if (info)
	  std::cout << "\nCombining result mod l and mod l^2:"
	    " Eigenvalue is " << ev_lpower;
      }
  }



#ifdef LIDIA_NAMESPACE
}				// end of namespace LiDIA
#endif
