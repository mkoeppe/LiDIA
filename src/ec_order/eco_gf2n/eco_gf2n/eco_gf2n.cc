//==============================================================================================
//
// This file is part of LiDIA --- a library for computational number theory
//
// Copyright (c) 1994--2001 the LiDIA Group.  All rights reserved.
//
// See http://www.informatik.tu-darmstadt.de/TI/LiDIA/
//
//----------------------------------------------------------------------------------------------
//
// $Id$
//
// Author : Volker Mueller
// Changes : See CVS log
//
//==============================================================================================



#ifdef HAVE_CONFIG_H
# include	"config.h"
#endif
#include	"LiDIA/path.h"
#include	"LiDIA/eco_gf2n.h"
#include	"LiDIA/trace_list.h"
#include	"LiDIA/rational_factorization.h"
#include	"LiDIA/timer.h"
#include	"LiDIA/gf2n.h"
#include	"LiDIA/gf_element.h"
#include	"LiDIA/galois_field.h"
#include	<fstream>
#include	<cstring>
#include	<cstdlib>

#include "LiDIA/osstream.h"
#include "LiDIA/isstream.h"


#ifdef LIDIA_NAMESPACE
namespace LiDIA
{
#endif



// WARNING: The following conversion routine only works if
// - the defining polynomial of degree d for gf2n was taken from GF2n.database
//   and
// - either no galois_field of degree d and characteristic 2 is in use or,
//   if such a galois_field is in use, its defining polynomial was taken from
//   GF2n.database, too.
// If several galois_fields with 2^d elements (i.e. with different defining
// polynomials) are in use, the behaviour of the conversion routine is
// undefined.

  gf_element eco_gf2n::convertToGFElement (const gf2n & a)
  {
    galois_field GF2n (2, gf2n::get_absolute_degree ());
    gf_element x (GF2n);

    static char prefix[10];
    char *oldprefix = gf2nIO::showprefix ();
    if (oldprefix != NULL)
        strcpy (prefix, oldprefix);
      gf2nIO::setprefix (gf2nIO::Dec);

      gf2nIO::base oldbase = gf2nIO::showbase ();

      osstream os;
      os << "(" << a << ")" << '\0';
      isstream is (extractString(os));
      is >> x;

      gf2nIO::setbase (oldbase);
      if (oldprefix != NULL)
        gf2nIO::setprefix (prefix);

      return x;
  }



  gf2n eco_gf2n::convertToGF2n (const gf_element & a)
  {
    osstream os;
    os << a << '\0';

    isstream is (extractString(os));
    char c;

    is >> c >> c >> c >> c >> c;

    gf2n x;

    is >> x;

    return x;
  }



  elliptic_curve < gf_element > eco_gf2n::generateCurve (const gf2n & a6)
  {
    galois_field theField (2, gf2n::get_absolute_degree ());
    gf_element One (theField);
    gf_element Zero (theField);

    One.assign_one ();
    Zero.assign_zero ();
    gf_element gfElementA6 = eco_gf2n::convertToGFElement (a6);

    elliptic_curve < gf_element > e (One, Zero, Zero, Zero, gfElementA6);
    return e;
  }



  gf2n eco_gf2n::j_invariant (const gf2n & a1, const gf2n & a2,
			      const gf2n & a3, const gf2n & a4,
			      const gf2n & a6)
  {
    debug_handler ("eco_gf2n", "j_invariant()");

    gf2n b2, b4, b6, b8, c4, c6, delta;

    square (b2, a1);
    multiply (b4, a1, a3);
    square (b6, a3);
    b8 = a1 * (a1 * a6 + a3 * a4) + a2 * a3 * a3 + a4 * a4;
    square (c4, b2);
    multiply (c6, c4, b2);
    delta = b2 * (b2 * b8 + b4 * b6) + b6 * b6;

    return (c4 * c4 * c4) / delta;
  }



  const char eco_gf2n::UNKNOWN = 5;
  const char eco_gf2n::ALL_SUBGRPS_INV = 1;
  const char eco_gf2n::ONE_SUBGRP_INV = 2;
  const char eco_gf2n::TWO_SUBGRPS_INV = 3;
  const char eco_gf2n::NO_SUBGRP_INV = 4;

  const char eco_gf2n::COMPUTE_SP_DEGREE = 5;
  const char eco_gf2n::DONT_COMPUTE_SP_DEGREE = 6;
  const char eco_gf2n::COMPUTE_SP_DEGREE_IF_ATKIN = 7;
  const char eco_gf2n::COMPUTE_SP_DEGREE_IF_ELKIES = 8;
  const char eco_gf2n::COMPUTE_SP_DEGREE_AND_ROOTS = 9;

  const char eco_gf2n::EV_RATIONAL_FUNCTION = 0;
  const char eco_gf2n::EV_DIVISION_POLYNOMIAL = 1;
  const char eco_gf2n::EV_RATIONAL_FUNCTION_TABLE = 2;
  const char eco_gf2n::EV_BABYSTEP_GIANTSTEP = 3;
  const char eco_gf2n::EV_FUNNY_BABYSTEP_GIANTSTEP = 4;
  const char eco_gf2n::EV_FUNNY_BABYSTEP_GIANTSTEP_SIGN = 5;

  const char eco_gf2n::TEST_X_COORDINATE = 0x80;
  const char eco_gf2n::TEST_Y_COORDINATE = 0x00;



//
// constructors /destructor
//

  eco_gf2n::eco_gf2n ()
  {
    c.set_mode (EXPAND);
    sp_type = UNKNOWN;
    degree_mode = DONT_COMPUTE_SP_DEGREE;
    ev_strategy = EV_DIVISION_POLYNOMIAL;
    schoof_bound = 0;
    atkin_bound = 50;
    power_bound = 25;
    degree_max_2_power = 16;
    use_corr_test = true;
    info = 0;
  }



  eco_gf2n::eco_gf2n (const ff_element & a6)
  {
    c.set_mode (EXPAND);
    sp_type = UNKNOWN;
    degree_mode = DONT_COMPUTE_SP_DEGREE;
    ev_strategy = EV_DIVISION_POLYNOMIAL;
    schoof_bound = 0;
    atkin_bound = 50;
    power_bound = 25;
    degree_max_2_power = 16;
    use_corr_test = true;
    set_curve (a6);
    info = 0;
  }



  eco_gf2n::~eco_gf2n ()
  {
  }



//
// assignments and access
//

  bool eco_gf2n::set_prime (udigit ll)
  {
    q = static_cast < udigit > (remainder (pn, static_cast < long >(ll)));

    l = ll;
    sp_type = UNKNOWN;
    return true;
  }



  eco_gf2n & eco_gf2n::operator = (const eco_gf2n & p)
  {
    if (this != &p)
      {
	pn = p.pn;
	A6 = p.A6;
	jinv = p.jinv;
	meq_pol = p.meq_pol;
	ftau[0].assign (p.ftau[0]);
	ftau[1].assign (p.ftau[1]);
	sp_type = p.sp_type;
	sp_degree = p.sp_degree;
	sp_degree_known = p.sp_degree_known;
	degree_mode = p.degree_mode;
	ev_strategy = p.ev_strategy;
	schoof_bound = p.schoof_bound;
	atkin_bound = p.atkin_bound;
	power_bound = p.power_bound;
	degree_max_2_power = p.degree_max_2_power;
	use_corr_test = p.use_corr_test;
	info = p.info;
	q = p.q;
	c = p.c;
	l = p.l;
      }
    return *this;
  }



  void eco_gf2n::set_curve (const ff_element & a6)
  {
    shift_left (pn, bigint (1), a6.relative_degree ());
    A6.assign (a6);
    compute_jinv (jinv, a6);
  }



  void eco_gf2n::set_strategy (char m)
  {
    if (m != COMPUTE_SP_DEGREE &&
	m != DONT_COMPUTE_SP_DEGREE &&
	m != COMPUTE_SP_DEGREE_IF_ATKIN &&
	m != COMPUTE_SP_DEGREE_IF_ELKIES && m != COMPUTE_SP_DEGREE_AND_ROOTS)
      lidia_error_handler ("eco_gf2n", "set_strategy::wrong mode");
    degree_mode = m;
  }


  void eco_gf2n::reset ()
  {
    sp_type = UNKNOWN;
    degree_mode = DONT_COMPUTE_SP_DEGREE;
    ev_strategy = EV_DIVISION_POLYNOMIAL;
    schoof_bound = 0;
    atkin_bound = 50;
    power_bound = 25;
    lower_bound_for_BG = 40;
    lower_bound_for_FBG = 40;
    degree_max_2_power = 16;
    use_corr_test = true;
    info = 0;
  }



  void eco_gf2n::set_ev_search_strategy (char m)
  {
    if (m != (EV_RATIONAL_FUNCTION | TEST_Y_COORDINATE) &&
	m != (EV_DIVISION_POLYNOMIAL | TEST_Y_COORDINATE) &&
	m != (EV_RATIONAL_FUNCTION | TEST_X_COORDINATE) &&
	m != (EV_DIVISION_POLYNOMIAL | TEST_X_COORDINATE) &&
	m != (EV_RATIONAL_FUNCTION_TABLE | TEST_Y_COORDINATE) &&
	m != (EV_RATIONAL_FUNCTION_TABLE | TEST_X_COORDINATE) &&
	m != EV_BABYSTEP_GIANTSTEP && m != EV_FUNNY_BABYSTEP_GIANTSTEP &&
	m != EV_FUNNY_BABYSTEP_GIANTSTEP_SIGN)
      lidia_error_handler ("eco_gf2n", "set_ev_search_strategy::wrong mode");

    ev_strategy = m;
  }



//****************************************************************
// read modular equation from File.
//
// File format :
//
//  prime l
//  X^0: < all exponents where coeff == 1
//  X^1: < all exponents where coeff == 1
//  ...
//  X^{l+1}: 0
//  checksum: sum over all exponents
//

  void eco_gf2n::read_meq (ff_pol & meq, int l)
  {
    power_tab < ff_element > jinv_power;
    jinv_power.initialize (jinv, l + 1);

    meq = ff_pol (l + 1);

    osstream oss;
    char const* name_env = getenv("LIDIA_ECO2");
    if (name_env != NULL && name_env[0] != '\0') {
	oss << name_env;
    }
    else {
	oss << LIDIA_ECO2_DB;
    }
    oss << MEQ2_PREFIX << l;
    std::string name = extractString(oss);

    std::ifstream ifs (name.c_str(), std::ios::in);
    if (!ifs)
      {
	lidia_error_handler ("read_meq", "Meq file can't be opened.");
	return;
      }

    int ll, check = 0;
    ff_element var;

    ifs >> ll;
    if (ll != l || ifs.fail())
      {
	lidia_error_handler ("read_meq", "Meq-file corrupted.");
	return;
      }

    std::string input;
    ifs >> input;

    for (int i = 0; i <= l + 1; i++)
      {
	var.assign_zero ();
	ifs >> input;

	while (!ifs.fail() && !input.empty() &&
	       input[0] >= '0' && input[0] <= '9') {
	    ll = strtol (input.c_str(), NULL, 10);
	    check += ll;
	    if (ll > l + 1)
	      {
		lidia_error_handler ("read_meq", "Meq file corrupted.");
		return;
	      }
	    add (var, var, jinv_power.get_power (ll));
	    ifs >> input;
	  }
	meq.set_coefficient (var, i);
      }

    ifs >> ll;

    if (ifs.fail() || check + l != ll)
      {
	lidia_error_handler ("read_meq", "Meq file corrupted.");
	return;
      }

    ifs.clear();
    ifs.close ();
  }



//***********************************************************************
// splitting type computation
//
// Comment: the same structure as in eco_prime.

  void eco_gf2n::compute_splitting_type ()
  {
    gf2n_polynomial f_xp, f_gcd, f_x;
    gf2n_poly_modulus meq_polmod;
    int def_degree = A6.relative_degree ();

    if (sp_type != UNKNOWN)
      return;

    read_meq (meq_pol, l);
    meq_polmod.build (meq_pol);

#ifdef TIMING
    timer t;

    t.set_print_mode ();
    t.start_timer ();
#endif

    Xq (f_xp, meq_polmod, def_degree);

    f_x.assign_x ();

    gcd (f_gcd, f_xp + f_x, meq_pol);


#ifdef TIMING
    t.stop_timer ();
    if (info)
      std::
	cout << "\nComputing X^q mod " << l << "-th MEQ takes time : " << t
	<< std::flush;
#endif

#ifdef DEBUG
    std::cout << "\nsplitting type: GCD = " << f_gcd << std::flush;
#endif

    lidia_size_t gcd_degree = f_gcd.degree ();

    if (gcd_degree == meq_pol.degree ())
      {
	ftau[0] = find_root (meq_pol, def_degree);
	sp_type = ALL_SUBGRPS_INV;

#ifdef DEBUG
	assert (meq_pol (ftau[0]).is_zero ());
	assert (def_degree % ftau[0].relative_degree () == 0);
#endif

	if (info)
	  {
	    std::
	      cout << "\nSplitting type: all " << l << "-groups invariant ";
	    std::cout << "==> Type (1 ... 1)" << std::flush;
	  }
	sp_degree = 1;
	sp_degree_known = true;
      }
    else if (gcd_degree == 1)
      {
	ftau[0] = find_root (f_gcd, def_degree);

#ifdef DEBUG
	assert (meq_pol (ftau[0]).is_zero ());
	assert (def_degree % ftau[0].relative_degree () == 0);
#endif
	sp_type = ONE_SUBGRP_INV;
	if (info)
	  {
	    std::
	      cout << "\nSplitting type: one " << l <<
	      "-group invariant ==> Type (1 ";
	    std::cout << l << ")" << std::flush;
	  }
	sp_degree_known = true;
	sp_degree = l;
      }
    else if (gcd_degree == 2)
      {
	ftau[0] = find_root (f_gcd, def_degree);
	divide (ftau[1], ff_element (f_gcd.const_term ()), ftau[0]);

#ifdef DEBUG
	assert (meq_pol (ftau[0]).is_zero ());
	assert (meq_pol (ftau[1]).is_zero ());
	assert (def_degree % ftau[0].relative_degree () == 0);
	assert (def_degree % ftau[1].relative_degree () == 0);
#endif
	sp_type = TWO_SUBGRPS_INV;
	if (info)
	  {
	    std::cout << "\nSplitting type: two " << l << "-groups invariant";
	    std::cout << std::flush;
	  }

	if (degree_mode == COMPUTE_SP_DEGREE ||
	    degree_mode == COMPUTE_SP_DEGREE_AND_ROOTS ||
	    degree_mode == COMPUTE_SP_DEGREE_IF_ELKIES)
	  {

	    lidia_size_t *d_list;

	    compute_d_list (d_list);

	    if (d_list[0] == 1)
	      sp_degree = static_cast < udigit > (d_list[1]);
	    else
	      {
		divide (meq_pol, meq_pol, f_gcd);
		meq_pol.make_monic ();
		meq_polmod.build (meq_pol);
		gf2n_polynomial h;

		remainder (h, f_xp, meq_pol);
		f_xp.assign (h);
		udigit lcm_value, hlcm;

		lcm_value = LiDIA::gcd (static_cast < udigit > (d_list[1]),
					static_cast < udigit > (d_list[2]));
		lcm_value = (static_cast < udigit > (d_list[1]) / lcm_value) *
		  static_cast < udigit > (d_list[2]);

		for (int i = 3; i <= d_list[0]; i++)
		  {
		    hlcm =
		      LiDIA::gcd (lcm_value,
				  static_cast < udigit > (d_list[i]));
		    lcm_value =
		      (lcm_value / hlcm) * static_cast < udigit > (d_list[i]);
		  }


#ifdef TIMING
		timer t;

		t.set_print_mode ();
		t.start_timer ();
#endif
		sp_degree =
		  static_cast < udigit >
		  (compute_degree
		   (f_xp, meq_polmod,
		    static_cast < lidia_size_t > (lcm_value)));
#ifdef TIMING
		t.stop_timer ();
		if (info)
		  std::
		    cout << "\nComputing splitting type takes time : " << t <<
		    std::flush;
#endif
	      }
	    delete[]d_list;
	    sp_degree_known = true;
	    if (info)
	      if (sp_degree < static_cast < unsigned int >(l - 1))
		std::
		  cout << " ==> Type (1 1 " << sp_degree << " ... " <<
		  sp_degree << ")" << std::flush;
	      else
		std::cout << " ==> Type (1 1 " << sp_degree << ")" << std::
		  flush;
	  }
	else
	  {
	    sp_degree_known = false;
	    if (info)
	      std::cout << " ==> Type (1 1 ? ... ?)" << std::flush;
	  }
      }
    else
      {
	sp_type = NO_SUBGRP_INV;
	if (info)
	  {
	    std::cout << "\nSplitting type: no " << l << "-group invariant";
	    std::cout << std::flush;
	  }

	if ((degree_mode == eco_gf2n::COMPUTE_SP_DEGREE ||
	     degree_mode == eco_gf2n::COMPUTE_SP_DEGREE_AND_ROOTS ||
	     degree_mode == eco_gf2n::COMPUTE_SP_DEGREE_IF_ATKIN)
	    && l <= atkin_bound)
	  {
	    lidia_size_t *d_list;

	    compute_d_list (d_list);

	    if (d_list[0] == 1)
	      sp_degree = static_cast < udigit > (d_list[1]);
	    else
	      {
		udigit lcm;

		lcm = LiDIA::gcd (static_cast < udigit > (d_list[1]),
				  static_cast < udigit > (d_list[2]));
		lcm =
		  (static_cast < udigit > (d_list[1]) / lcm) * static_cast <
		  udigit > (d_list[2]);

		for (int i = 3; i <= d_list[0]; i++)
		  {
		    udigit h;

		    h = LiDIA::gcd (lcm, static_cast < udigit > (d_list[i]));
		    lcm = (lcm / h) * static_cast < udigit > (d_list[i]);
		  }
#ifdef TIMING
		timer t;

		t.set_print_mode ();
		t.start_timer ();
#endif
		sp_degree =
		  static_cast < udigit >
		  (compute_degree
		   (f_xp, meq_polmod, static_cast < lidia_size_t > (lcm)));
#ifdef TIMING
		t.stop_timer ();
		if (info)
		  std::
		    cout << "\nComputing splitting type takes time : " << t <<
		    std::flush;
#endif
	      }
	    delete[]d_list;
	    sp_degree_known = true;
	    if (info)
	      if (sp_degree < static_cast < unsigned int >(meq_pol.degree ()))
		std::
		  cout << " ==> Type (" << sp_degree << " ... " << sp_degree
		  << ")" << std::flush;
	      else
		std::cout << " ==> Type (" << sp_degree << ")" << std::flush;
	  }
	else
	  {
	    sp_degree_known = false;
	    if (info)
	      std::cout << " ==> Type (? ... ?)" << std::flush;
	  }
      }
  }



  bool eco_gf2n::is_elkies (void)
  {
    return (sp_type == ONE_SUBGRP_INV ||
	    sp_type == TWO_SUBGRPS_INV || sp_type == ALL_SUBGRPS_INV);
  }



//-----------------------------------------------------------
// Assumes, that l is an elkies-prime and
// that ftau[0] contains a root of the
// modular polynomial.
//

  void eco_gf2n::compute_trace_elkies ()
  {
    ff_pol fC;			// the Elkies - polynomial
    ff_polmod fCpm;		// and as gf2n_poly_modulus
    lidia_size_t ev;		// eigenvalue of Frobenius-endomorphism
    udigit u_ev;		// ev as udigit

    //* * *  compute the Elkies polynomial * * * *

#ifdef TIMING
    timer t;

    t.set_print_mode ();
    t.start_timer ();
#endif

    divisor_of_divpol_via_lercier (fC, eco_gf2n::WITH_VELU);

#ifdef TIMING
    t.stop_timer ();
    if (info)
      {
	std::
	  cout <<
	  "\nComputing Elkies polynomial with Lercier's method takes ";
	std::cout << "time : " << t << std::flush;
      }
#endif


#ifdef DEBUG
    if (fC.degree () >= 2)
      {
	fCpm.build (fC);
	ff_pol div;

	compute_psi (div, static_cast < lidia_size_t > (l), fCpm);
	assert (div.is_zero ());
      }
#endif

    //* * * compute the eigenvalue  * * * *

    if (!find_eigenvalue (ev, fC, 0, 1, false))
      {
	lidia_error_handler ("eco_gf2n::compute_trace_elkies(...)",
			     "No eigenvalue found.");
      }
    else
      {
	udigit ev1_order, ev2_order;

	u_ev = static_cast < udigit > (ev);

	ff1::set_characteristic (l);
	ff1 ff1_ev (u_ev);

	ev1_order = ff1_ev.multiplicative_order ();
	if (!(ev1_order & 1))
	  ev1_order >>= 1;

	invert (ff1_ev, ff1_ev);
	multiply (ff1_ev, ff1_ev, q);

	ev2_order = ff1_ev.multiplicative_order ();
	if (!(ev2_order & 1))
	  ev2_order >>= 1;


	if (test_y_coordinate () || test_x_coordinate_uniquely ())
	  {
	    // only one possibility
	    c.set_capacity (1);
	    c.set_size (1);

	    c[0] = divide_mod (eco_gf2n::q, u_ev, this->l);
	    c[0] = add_mod (c[0], u_ev, this->l);
	    if (info)
	      std::cout << "\n==> c = " << c[0] << " mod " << l << std::flush;
	  }
	else
	  {
	    // test on x-coordinate ==> two possibilities
	    c.set_capacity (2);
	    c.set_size (2);

	    c[0] = divide_mod (eco_gf2n::q, u_ev, this->l);
	    c[0] = add_mod (c[0], u_ev, this->l);
	    if (c[0] == 0)
	      {
		c.set_capacity (1);
		c.set_size (1);
	      }
	    else
	      c[1] = this->l - c[0];

	    if (info)
	      if (c[0] != 0)
		std::cout << "\n==> c = +/- " << c[0] << " mod " << l << std::
		  flush;
	      else
		std::cout << "\n==> c = 0 mod " << l << std::flush;
	  }

	if (l * comparator < lidia_size_t >::min (ev1_order, ev2_order)
	    <= power_bound && sp_type == TWO_SUBGRPS_INV)
	  {
	    int lpower = l * l, ev_lpower = 0;

	    if (info)
	      std::cout << "\n\nWorking on " << l << "^2 : " << std::flush;

	    compute_powers_of_l (fC, ev, ev1_order, ev2_order, ev_lpower,
				 lpower);
	    eco_gf2n::q =
	      static_cast < udigit >
	      (remainder (pn, static_cast < long >(lpower)));
	    eco_gf2n::l = lpower;

	    if (test_y_coordinate () || (c.size () == 1 && c[0] != 0))
	      {
		// only one possibility

		c[0] = divide_mod (eco_gf2n::q, ev_lpower, lpower);
		c[0] = add_mod (c[0], ev_lpower, lpower);
		if (info)
		  std::cout << "\n==> c = " << c[0] << " mod " << eco_gf2n::
		    l << std::flush;
	      }
	    else
	      {
		// test on x-coordinate ==> two possibilities

		std::cout<<"\nhere ";

		c.set_size(2);
		c[0] = divide_mod (eco_gf2n::q, ev_lpower, lpower);
		c[0] = add_mod (c[0], ev_lpower, lpower);
		if (c[0] == 0)
		  c.set_size (1);
		else
		  c[1] = lpower - c[0];

		if (info)
		  if (c[0] != 0)
		    std::
		      cout << "\n==> c = +/- " << c[0] << " mod " <<
		      eco_gf2n::l << std::flush;
		  else
		    std::cout << "\n==> c = 0 mod " << l << std::flush;
	      }
	  }
      }
  }



//---------------------------------------------------------------------
// Atkin's ideas are used to determine values for c mod l.
//

//-----------------------------------------------------------
// computes for the given value of d all possibilities for the
// Frobenius trace and adds these to the internal vector c.
//

  void eco_gf2n::compute_c_for_d (lidia_size_t d)
  {
    if (sp_type == NO_SUBGRP_INV)
      {
	base_vector < ff2_element > adl2;
	ff2_element yl2;
	ff2_element Q;
	long i, sz;
	int pos;
	udigit h;

	ff2_element::init_field (l);
	nearly_all_of_order (d, adl2);

	sz = adl2.size ();
	Q = ff1 (q);

	for (i = 0; i < sz; i++)
	  {
	    add (yl2, adl2[i], ff2 (ff1 (1)));
	    multiply (yl2, yl2, yl2);
	    multiply (yl2, yl2, Q);
	    divide (yl2, yl2, adl2[i]);
	    yl2.as_udigit (h);
	    h = sqrt_mod (h, l);
	    if (!c.bin_search (h, pos))
	      c.insert (h);

	    if (h != 0)
	      {
		h = negate_mod (h, l);
		if (!c.bin_search (h, pos))
		  c.insert (h);
	      }
	  }
      }				// end of sp_TYPE = (d...d)
    else
      {
	// sp_TYPE = (1 1 d...d)

	// c = +-  (x+1) * sqrt [ q / x ], x in GF(l), ord(x) = d

	base_vector < ff1_element > ad;
	ff1_element y;
	ff1_element Q;
	long i, sz;
	int pos;
	udigit h;


	ff1_element::set_characteristic (l);
	nearly_all_of_order (d, ad);

	sz = ad.size ();
	Q = ff1 (q);

	for (i = 0; i < sz; i++)
	  {
	    invert (y, ad[i]);
	    multiply (y, Q, y);
	    sqrt (y, y);
	    add (ad[i], ad[i], ff1_element (1));
	    multiply (y, y, ad[i]);
	    h = y.as_udigit ();
	    if (!c.bin_search (h, pos))
	      c.insert (h);

	    if (h != 0)
	      {
		negate (y, y);
		h = y.as_udigit ();
		if (!c.bin_search (h, pos))
		  c.insert (h);
	      }
	  }
      }				// end of split type (11d...d)
  }



  void eco_gf2n::compute_trace_atkin ()
  {
    if (!sp_degree_known)
      {
	lidia_size_t *d_list = NULL;

	c.set_size (0);
	compute_d_list (d_list);

	if (d_list[0] == 1)
	  compute_c_for_d (d_list[1]);
	else
	  {
	    for (int i = 1; i <= d_list[0]; i++)
	      {
		compute_c_for_d (d_list[i]);
	      }
	  }

	delete[]d_list;
      }
    else
      {
	// splitting degree known
	// Consider Splitting type is (1...1) or (1 l)

	if (sp_type == ALL_SUBGRPS_INV || sp_type == ONE_SUBGRP_INV)
	  {
	    //  c = +- 2 * sqrt (q) in GF(l)
	    udigit h;

	    c.set_capacity (2);
	    h = sqrt_mod (q, l);
	    h = add_mod (h, h, l);	// h = 2 * sqrt (p^n mod l)
	    c[0] = h;
	    c[1] = negate_mod (h, l);
	  }
	else if (sp_type == NO_SUBGRP_INV)
	  {
	    // Splitting type is (d...d), d = sp_degree
	    if (sp_degree == 2)
	      {
		c.set_capacity (1);
		c[0] = 0;
	      }
	    else
	      {
		c.set_size (0);
		compute_c_for_d (sp_degree);
	      }
	  }
	else if (sp_type == TWO_SUBGRPS_INV)
	  {
	    // c = +-  (x+1) * sqrt [ q / x ], x in GF(l), ord(x) = d

	    if (sp_degree == 2)
	      {
		c.set_capacity (1);
		c[0] = 0;
	      }
	    else
	      {
		c.set_size (0);
		compute_c_for_d (sp_degree);
	      }			// end of split type (11d...d)
	  }
	else
	  {
	    lidia_error_handler ("eco_gf2n::compute_trace_atkin()",
				 "Unknown splitting type for Atkin prime.");
	  }
      }

    if (info)
      {
	c.sort ();
	std::cout << "\n==> c mod " << l << " in " << c << std::flush;
      }


    if (l <= schoof_bound)
      if (sp_degree & 1)
	schoof_algorithm (sp_degree);
      else
	schoof_algorithm (sp_degree / 2);
  }



  bool eco_gf2n::
    probabilistic_correctness_proof (const bigint & res, lidia_size_t tests)
  {
    elliptic_curve < gf_element > e = eco_gf2n::generateCurve (A6);
    return e.probabilistic_test_of_group_order (res, tests);
  }



  bigint eco_gf2n::compute_group_order ()
  {
    elliptic_curve < gf_element > e = eco_gf2n::generateCurve (A6);

    trace_list::set_info_mode (info);
    trace_list tl;
    trace_mod tm;
    udigit ll;
    bigint ec_order;

    //
    // initialize timer
    //

    timer t;

    t.set_print_mode (HMS_MODE);

    if (A6.relative_degree () < 15)
      {
	if (info)
	  {
	    std::
	      cout << "\nElliptic curve is already defined over the subfield";
	    std::cout << " GF(2^" << A6.relative_degree () << ") of GF(2^";
	    std::cout << gf2n::get_absolute_degree () << "). ";
	    std::
	      cout <<
	      "\n==> Using the standard Babystep Giantstep algorithm in"
	      " GF(2^" << A6.relative_degree () << ") and lifting ... \n";
	  }

	t.start_timer ();
	ec_order = e.group_order ();
	shift_left (pn, bigint (1), gf2n::extension_degree ());
	t.stop_timer ();

	if (info)
	  {
	    std::cout << "\n\nGroup Order (over GF(2^" << gf2n::
	      get_absolute_degree ();
	    std::cout << ")) is " << ec_order;
	    std::cout << "\nComputing Time : " << t << std::flush;
	  }

	if (use_corr_test)
	  if (probabilistic_correctness_proof (ec_order))
	    {
	      if (info)
		{
		  std::
		    cout << "\n\nProbabilistic correctness proof agrees !!";
		  std::cout << "\n\nField Size : ";
		  std::cout << pn;
		  std::cout << "\nCoefficient a6 : ";
		  std::cout << A6;
		  std::cout << "\nGroup order is " << ec_order << "." << std::
		    endl;
		}
	    }
	  else
	    {
	      std::cout << "\n\nField Size : ";
	      std::cout << pn;
	      std::cout << "\nCoefficient a6 : ";
	      std::cout << A6;
	      std::
		cout << "\n\nOrder candidate : " << ec_order << "." << std::
		endl;
	      lidia_error_handler ("eco_gf2n::compute_group_order",
				   "Probabilistic correctness rejects "
				   "candidate !!");
	    }
	return ec_order;
      }

    if (A6.relative_degree () < gf2n::get_absolute_degree ())
      {
	if (info)
	  {
	    std::
	      cout <<
	      "\nElliptic curve is already defined over the subfield GF(2^";
	    std::cout << A6.relative_degree () << ")";
	    std::cout << "\n==> Computing group order over subfield, then ";
	    std::cout << "lift group order ...\n\n";
	  }
      }

    //
    // compute trace mod power of 2
    //

    if (info)
      {
	std::
	  cout <<
	  "\n\n-------------------------------------------------------\n";
	std::cout << "Working on prime l = 2" << std::endl;
      }

    tl.set_curve (e);
    compute_mod_2_power ();
    tm.set_vector (get_prime (), get_relation ());
    tl.append (tm);

    //
    // compute trace mod primes and store the residues
    // in the trace_list
    //

    bool end_loop;

    ll = 3;
    end_loop = false;

    do
      {
	if (info)
	  {
	    std::
	      cout <<
	      "\n\n-------------------------------------------------------\n";
	    std::cout << "Working on prime l = " << ll << std::endl;
	  }

	if (set_prime (ll))
	  {
	    if (info)
	      std::cout << std::endl;

	    compute_splitting_type ();

	    if (is_elkies ())
	      compute_trace_elkies ();
	    else
	      compute_trace_atkin ();

	    tm.set_vector (get_prime (), get_relation ());
	    end_loop = tl.append (tm);
	  }
	else if (info)
	  std::cout << " ... not available!" << std::endl;

	ll = next_prime (ll + 1);
      }
    while ((!end_loop));


    //
    // compute the group order with babystep/giantstep algorithm
    //

    if (info)
      {
	std::
	  cout <<
	  "\n-------------------------------------------------------\n";
	std::cout << "\nStarting Babystep Giantstep Phase";
	std::cout << "\n\nFinal TraceList : " << tl << std::endl;
      }
    ec_order = tl.bg_search_for_order ();

    if (info)
      {
	t.stop_timer ();
	std::cout << "\n\nOutput of BG algorithm = " << ec_order << std::
	  flush;
	std::cout << "\nComputing Time : " << t << std::flush;
      }

    if (A6.relative_degree () < gf2n::get_absolute_degree ())
      {
	bigint c1, ci, cim1 (2), cip1;
	int loops = gf2n::get_absolute_degree () / A6.relative_degree ();

	subtract (c1, pn + 1, ec_order);
	ci.assign (c1);

	for (int i = 1; i < loops; i++)
	  {
	    cip1.assign (ci * c1 - pn * cim1);
	    cim1.assign (ci);
	    ci.assign (cip1);
	  }

	shift_left (pn, bigint (1), gf2n::get_absolute_degree ());
	ec_order.assign (pn + 1 - ci);

	if (info)
	  {
	    std::cout << "\nLifted group order (over GF(2^" << gf2n::
	      extension_degree ();
	    std::cout << ")) = " << ec_order;
	  }
      }

    //
    // Verification of the group order
    //

    if (use_corr_test)
      if (probabilistic_correctness_proof (ec_order))
	{
	  if (info)
	    {
	      std::cout << "\n\nProbabilistic correctness proof agrees !!";
	      std::cout << "\n\nField Size : ";
	      std::cout << pn;
	      std::cout << "\nCoefficient a6 : ";
	      std::cout << A6;
	      std::cout << "\nGroup order is " << ec_order << "." << std::
		endl;
	    }
	}
      else
	{
	  std::cout << "\n\nField Size : ";
	  std::cout << pn;
	  std::cout << "\nCoefficient a6 : ";
	  std::cout << A6;
	  std::cout << "\n\nOrder candidate : " << ec_order << "." << std::
	    endl;
	  lidia_error_handler ("eco_gf2n::compute_group_order",
			       "Probabilistic correctness rejects candidate !!");
	}

    return ec_order;
  }



//--------------------------------------------------------------------
// finally a general function which can be used generally.

  bigint compute_group_order (const gf2n & a6, int info)
  {
    eco_gf2n ep (a6);

    ep.set_info_mode (info);
    return ep.compute_group_order ();
  }



  bigint compute_group_order (const gf2n & a2, const gf2n & a6, int info)
  {
    eco_gf2n ep (a6);
    bigint eord, pn;
    gf2n A2 (a2);

    if (info && A2.trace () == 1)
      std::cout << "\nComputing the group order of the twist ....\n\n";

    ep.set_info_mode (info);
    eord = ep.compute_group_order ();

    if (A2.trace () == 1)
      {
	shift_left (pn, bigint (1), gf2n::extension_degree ());
	eord = 2 * (pn + 1) - eord;
      }
    return eord;
  }



  bigint compute_group_order (const gf2n & aa1, const gf2n & aa2,
			      const gf2n & aa3, const gf2n & aa4,
			      const gf2n & aa6, int info)
  {
    const int number_of_tests = 5;

    gf_element gfElaa1 = eco_gf2n::convertToGFElement (aa1);
    gf_element gfElaa2 = eco_gf2n::convertToGFElement (aa2);
    gf_element gfElaa3 = eco_gf2n::convertToGFElement (aa3);
    gf_element gfElaa4 = eco_gf2n::convertToGFElement (aa4);
    gf_element gfElaa6 = eco_gf2n::convertToGFElement (aa6);

    elliptic_curve < gf_element > e (gfElaa1, gfElaa2, gfElaa3, gfElaa4,
				     gfElaa6);

    point < gf_element > P (e);
    int n (gf2n::extension_degree ());
    bigint q, eord, q_sqrt;
    int i;
    bool success;

    if (!aa1.is_zero ())
      {
	// transform to (a2, a6)
	gf2n a2, a6;

	power (a6, aa1, 3);
	a2 = (aa1 * aa2 + aa3) / a6;
	invert (a6, eco_gf2n::j_invariant (aa1, aa2, aa3, aa4, aa6));
	if (info)
	  {
	    std::
	      cout << "\nTransformation to elliptic curve Y^2 + X Y = X^3 + ";
	    std::cout << a2 << " X^2 + " << a6;
	    std::cout << "\nComputing group order of that curve ...";
	  }
	return compute_group_order (a2, a6, info);
      }

    add (eord, q, 1);
    success = true;

    for (i = 0; i < number_of_tests; i++)
      {
	// possible order q+1
	P = e.random_point ();
	multiply (P, eord, P);
	if (!P.is_zero ())
	  {
	    success = false;
	    break;
	  }
      }

    if (success)
      return eord;

    if (n & 1)
      {
	// possible orders q +1 + / - sqrt(2*q)
	shift_left (q_sqrt, bigint (1), (n + 1) / 2);
	add (eord, q, 1);
	add (eord, eord, q_sqrt);
	success = true;

	for (i = 0; i < number_of_tests; i++)
	  {
	    P = e.random_point ();
	    multiply (P, eord, P);
	    if (!P.is_zero ())
	      {
		success = false;
		break;
	      }
	  }

	if (success)
	  return eord;

	add (eord, q, 1);
	subtract (eord, eord, q_sqrt);
	success = true;

	for (i = 0; i < number_of_tests; i++)
	  {
	    P = e.random_point ();
	    multiply (P, eord, P);
	    if (!P.is_zero ())
	      {
		success = false;
		break;
	      }
	  }

	if (success)
	  return eord;
      }
    else
      {
	shift_left (q_sqrt, bigint (1), n / 2);
	add (eord, q, 1);
	add (eord, eord, q_sqrt);
	success = true;

	for (i = 0; i < number_of_tests; i++)
	  {
	    P = e.random_point ();
	    multiply (P, eord, P);
	    if (!P.is_zero ())
	      {
		success = false;
		break;
	      }
	  }

	if (success)
	  return eord;

	add (eord, q, 1);
	subtract (eord, eord, q_sqrt);
	success = true;

	for (i = 0; i < number_of_tests; i++)
	  {
	    P = e.random_point ();
	    multiply (P, eord, P);
	    if (!P.is_zero ())
	      {
		success = false;
		break;
	      }
	  }

	if (success)
	  return eord;

	shift_left (q_sqrt, bigint (1), n / 2 + 1);
	add (eord, q, 1);
	add (eord, eord, q_sqrt);
	success = true;

	for (i = 0; i < number_of_tests; i++)
	  {
	    P = e.random_point ();
	    multiply (P, eord, P);
	    if (!P.is_zero ())
	      {
		success = false;
		break;
	      }
	  }

	if (success)
	  return eord;

	add (eord, q, 1);
	subtract (eord, eord, q_sqrt);
	success = true;

	for (i = 0; i < number_of_tests; i++)
	  {
	    P = e.random_point ();
	    multiply (P, eord, P);
	    if (!P.is_zero ())
	      {
		success = false;
		break;
	      }
	  }

	if (success)
	  return eord;
      }
    lidia_error_handler ("eco_gf2n",
			 "compute_group_order::no solution found");
    return 0;
  }



#ifdef LIDIA_NAMESPACE
}				// end of namespace LiDIA
#endif
