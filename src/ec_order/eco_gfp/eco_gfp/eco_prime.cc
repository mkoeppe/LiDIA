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
//	Author	: Markus Maurer, Volker Mueller
//	Changes	: See CVS log
//
//==============================================================================================


#ifdef HAVE_CONFIG_H
# include	"config.h"
#endif
#include	"LiDIA/eco_prime.h"
#include	"LiDIA/trace_list.h"
#include	"LiDIA/rational_factorization.h"
#include	"LiDIA/timer.h"
#ifdef DEBUG
#include	<cassert>
#endif


#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif

//
// constants
//

const char eco_prime::UNKNOWN          = 5;
const char eco_prime::ALL_SUBGRPS_INV  = 1;
const char eco_prime::ONE_SUBGRP_INV   = 2;
const char eco_prime::TWO_SUBGRPS_INV  = 3;
const char eco_prime::NO_SUBGRP_INV    = 4;

const char eco_prime::COMPUTE_SP_DEGREE	   = 5;
const char eco_prime::DONT_COMPUTE_SP_DEGREE      = 6;
const char eco_prime::COMPUTE_SP_DEGREE_IF_ATKIN  = 7;
const char eco_prime::COMPUTE_SP_DEGREE_IF_ELKIES = 8;
const char eco_prime::COMPUTE_SP_DEGREE_AND_ROOTS = 9;

const char eco_prime::EV_RATIONAL_FUNCTION = 0;
const char eco_prime::EV_DIVISION_POLYNOMIAL = 1;
const char eco_prime::EV_RATIONAL_FUNCTION_TABLE = 2;
const char eco_prime::EV_BABYSTEP_GIANTSTEP = 3;
const char eco_prime::EV_FUNNY_BABYSTEP_GIANTSTEP = 4;
const char eco_prime::EV_OPTIMAL = 5;

const char eco_prime::TEST_X_COORDINATE = 0x80;
const char eco_prime::TEST_Y_COORDINATE = 0x00;



//
// constructors /destructor
//

eco_prime::eco_prime ()
{
	c.set_mode(EXPAND);
	sp_type = UNKNOWN;
	degree_mode = DONT_COMPUTE_SP_DEGREE;
	ev_strategy = EV_OPTIMAL;
	schoof_bound = 0;
	atkin_bound = 50;
	lower_bound_for_BG = 40;
	lower_bound_for_FBG = 40;
	info = 0;
	use_corr_test = true;
}



eco_prime::eco_prime (const gf_element & a4, const gf_element & a6)
{
	c.set_mode(EXPAND);
	sp_type = UNKNOWN;
	degree_mode = DONT_COMPUTE_SP_DEGREE;
	ev_strategy = EV_OPTIMAL;
	schoof_bound = 0;
  	atkin_bound = 50;
	lower_bound_for_BG = 40;
	lower_bound_for_FBG = 40;
	info = 0;
	use_corr_test = true;

	this->set_curve(a4, a6);
}



eco_prime::~eco_prime () 
{
}



//
// assignments and access
//

bool eco_prime
::set_prime(udigit ll)
{
	if (meq_prime::set_prime(ll)) {
		q = static_cast<udigit>(remainder(pn, static_cast<long>(ll)));
		l = ll;
		sp_type = UNKNOWN;
		return true;
	}
	else
		return false;
}



eco_prime & eco_prime::operator = (const eco_prime & p)
{
	if (this != &p) {
		meq_pol = p.meq_pol;
		meq_pol_mod = p.meq_pol_mod;
		ftau[0].assign(p.ftau[0]);
		ftau[1].assign(p.ftau[1]);
		sp_type = p.sp_type;
		sp_degree = p.sp_degree;
		degree_mode = p.degree_mode;
		sp_degree_known = p.sp_degree_known;
		atkin_bound = p.atkin_bound;
		schoof_bound = p.schoof_bound;
		lower_bound_for_BG = p.lower_bound_for_BG;
		lower_bound_for_FBG = p.lower_bound_for_FBG;
		info = p.info;
		use_corr_test = p.use_corr_test;
		q = p.q;
		c = p.c;
		l = p.l;
	}
	return *this;
}



void eco_prime::reset()
{
	sp_type = UNKNOWN;
	degree_mode = DONT_COMPUTE_SP_DEGREE;
	ev_strategy = EV_OPTIMAL;
  	atkin_bound = 50;
	schoof_bound = 0;
	lower_bound_for_BG = 40;
	lower_bound_for_FBG = 40;
	info = 0;
	use_corr_test = true;
}



void eco_prime::set_curve (const gf_element & a, const gf_element & b)
{
  p.assign(a.characteristic());
  pn.assign(a.get_field().number_of_elements());
  if (p != pn) {
    lidia_error_handler("eco_prime", "finite fields of extension degree > 1"
			" not yet implemented");
  }

  bigmod::set_modulus(p);

  A.assign(a.polynomial_rep().const_term());
  B.assign(b.polynomial_rep().const_term());
  compute_jinv(jinv, A, B);
  
  E_.set_coefficients(a, b);
}



void eco_prime::set_strategy (char m)
{
	if (m != COMPUTE_SP_DEGREE &&
	    m != DONT_COMPUTE_SP_DEGREE &&
	    m != COMPUTE_SP_DEGREE_IF_ATKIN &&
	    m != COMPUTE_SP_DEGREE_IF_ELKIES &&
	    m != COMPUTE_SP_DEGREE_AND_ROOTS)
		lidia_error_handler("eco_prime", "set_strategy::wrong mode");
	degree_mode = m;
}



void eco_prime::set_ev_search_strategy(char m)
{
	if (m != (EV_RATIONAL_FUNCTION   | TEST_Y_COORDINATE) &&
	    m != (EV_DIVISION_POLYNOMIAL | TEST_Y_COORDINATE) &&
	    m != (EV_RATIONAL_FUNCTION   | TEST_X_COORDINATE) &&
	    m != (EV_DIVISION_POLYNOMIAL | TEST_X_COORDINATE) &&
	    m != (EV_RATIONAL_FUNCTION_TABLE   | TEST_Y_COORDINATE) &&
	    m != (EV_RATIONAL_FUNCTION_TABLE   | TEST_X_COORDINATE) &&
	    m != EV_BABYSTEP_GIANTSTEP && m != EV_FUNNY_BABYSTEP_GIANTSTEP &&
	    m != EV_OPTIMAL)
		lidia_error_handler("eco_prime", "set_ev_search_strategy::wrong mode");
	ev_strategy = m;
}



//--------------------------------------------------------------------
//
// splitting type computation
//

void eco_prime::compute_splitting_type()
{
  ff_pol f_xp, f_gcd, f_x;
  
  if (sp_type != UNKNOWN)
    return;
  
  meq_pol.set_modulus(p);
  f_x.set_modulus(p);
  f_xp.set_modulus(p);
  f_gcd.set_modulus(p);
  
  f_x.assign_x();
  this->build_poly_in_X(meq_pol, jinv);
  meq_pol_mod.build(meq_pol);
  
#ifdef TIMING
  timer t;
  t.set_print_mode();
  t.start_timer();
#endif
  
  power_x(f_xp, pn, meq_pol_mod);
  f_gcd = gcd(f_xp - f_x, meq_pol);
  
#ifdef TIMING
  t.stop_timer();
  if (info)
    std::cout << "\nComputing X^q mod " << l << "-th MEQ takes "
      "time : " << t << std::flush;
#endif
  
#ifdef DEBUG
  std::cout << "\nsplitting type: pn = " << pn << std::flush;
  std::cout << "\nsplitting type: f_xp = " << f_xp << std::flush;
  std::cout << "\nsplitting type: GCD = " << f_gcd << std::flush;
#endif
  
  lidia_size_t gcd_degree = f_gcd.degree();
  
  if (gcd_degree == meq_pol.degree()) {
    ftau[0] = find_root(f_gcd);
    sp_type = ALL_SUBGRPS_INV;
    
#ifdef DEBUG
    assert(meq_pol(ftau[0].mantissa()).is_zero());
#endif
    
    if (info) 
      {
	std::cout << "\nSplitting type: all " << l << "-groups invariant ";
	std::cout << "==> Type (1 ... 1)" << std::flush;
      }
    sp_degree = 1;
    sp_degree_known = true;
  }
  else
    if (gcd_degree == 1) 
      {
	ftau[0] = find_root(f_gcd);
	
#ifdef DEBUG
	assert(meq_pol(ftau[0].mantissa()).is_zero());
#endif
	sp_type = ONE_SUBGRP_INV;
	if (info) {
	  std::cout << "\nSplitting type: one " << l << "-group invariant ==> Type (1 ";
	  std::cout << l << ")" << std::flush;
	}
	sp_degree_known = true;
	sp_degree = l;
      }
    else
      if (gcd_degree == 2) {
	ftau[0] = find_root(f_gcd);
	divide(ftau[1], ff_element(f_gcd.const_term()), ftau[0]);
#ifdef DEBUG
	assert(meq_pol(ftau[0].mantissa()).is_zero());
	assert(meq_pol(ftau[1].mantissa()).is_zero());
#endif
	sp_type = TWO_SUBGRPS_INV;
	if (info) {
	  std::cout << "\nSplitting type: two " << l << "-groups invariant";
	  std::cout << std::flush;
				}
	if (degree_mode == COMPUTE_SP_DEGREE ||
	    degree_mode == COMPUTE_SP_DEGREE_AND_ROOTS ||
	    degree_mode == COMPUTE_SP_DEGREE_IF_ELKIES) 
	  {
	    
	    lidia_size_t *d_list;
	    
	    compute_d_list(d_list);
	    
	    if (d_list[0] == 1)
	      sp_degree = static_cast<udigit>(d_list[1]);
	    else {
	      divide(meq_pol, meq_pol, f_gcd);
	      meq_pol_mod.build(meq_pol);
	      remainder(f_xp, f_xp, meq_pol);
	      
	      udigit lcm, hlcm;
	      
	      lcm = LiDIA::gcd(static_cast<udigit>(d_list[1]),
			       static_cast<udigit>(d_list[2]));
	      lcm = (static_cast<udigit>(d_list[1]) / lcm) * 
		static_cast<udigit>(d_list[2]);
	      
	      for (int i = 3; i <= d_list[0]; i++)
		{
		  hlcm = LiDIA::gcd(lcm, static_cast<udigit>(d_list[i]));
		  lcm = (lcm / hlcm) * static_cast<udigit>(d_list[i]);
		}
	      
#ifdef TIMING
	      timer t;
	      t.set_print_mode();
	      t.start_timer();
#endif
	      sp_degree = static_cast<udigit>(compute_degree 
					      (f_xp, meq_pol_mod, lcm));
	      
#ifdef TIMING
	  t.stop_timer();
	  if (info)
	    std::cout << "\nComputing splitting type takes time : " << t << std::flush;
#endif
	}
	delete[] d_list;
	sp_degree_known = true;
	if (info)
	  if (sp_degree < static_cast<unsigned int>(l)-1)
	    std::cout << " ==> Type (1 1 " << sp_degree << " ... " << sp_degree << ")" << std::flush;
	  else std::cout << " ==> Type (1 1 " << sp_degree << ")" << std::flush;
				}
				else {
				  sp_degree_known = false;
				  if (info)
				    std::cout << " ==> Type (1 1 ? ... ?)" << std::flush;
				}
			}
			else {
			  sp_type = NO_SUBGRP_INV;
			  if (info) {
			    std::cout << "\nSplitting type: no " << l << "-group invariant";
			    std::cout << std::flush;
			  }
			  
			  if ((degree_mode == eco_prime::COMPUTE_SP_DEGREE ||
			       degree_mode == eco_prime::COMPUTE_SP_DEGREE_AND_ROOTS ||
			       degree_mode == eco_prime::COMPUTE_SP_DEGREE_IF_ATKIN) &&
			      l <= atkin_bound) {
			    lidia_size_t *d_list;
			    compute_d_list(d_list);
			    
			    if (d_list[0] == 1)
			      sp_degree = static_cast<udigit>(d_list[1]);
			    else {
			      udigit lcm;
			      
			      lcm = LiDIA::gcd(static_cast<udigit>(d_list[1]),
					       static_cast<udigit>(d_list[2]));
			      lcm = (static_cast<udigit>(d_list[1]) / lcm) * 
				static_cast<udigit>(d_list[2]);
			      for (int i = 3; i <= d_list[0]; i++)
				{
				  udigit hlcm;
				  hlcm = LiDIA::gcd(lcm, static_cast<udigit>(d_list[i]));
				  lcm = (lcm / hlcm) * static_cast<udigit>(d_list[i]);
				}
#ifdef TIMING
			      timer t;
			      t.set_print_mode();
			      t.start_timer();
#endif
			      
			      sp_degree = static_cast<udigit>(compute_degree (f_xp, meq_pol_mod,
									      lcm));
#ifdef TIMING
			      t.stop_timer();
			      if (info)
				std::cout << "\nComputing splitting type takes time : " << t << std::flush;
#endif
			    }
			    delete[] d_list;
			    sp_degree_known = true;
			    if (info)
			      if (sp_degree < static_cast<unsigned int>(meq_pol.degree()))
				std::cout << " ==> Type (" << sp_degree << " ... " << sp_degree << ")" << std::flush;
			      else std::cout << " ==> Type (" << sp_degree << ")" << std::flush;
			  }
			  else {
			    sp_degree_known = false;
			    if (info)
			      std::cout << " ==> Type (? ... ?)" << std::flush;
			  }
			}
}



bool eco_prime::is_elkies(void)
{
	return (sp_type == ONE_SUBGRP_INV  ||
		sp_type == TWO_SUBGRPS_INV ||
		sp_type == ALL_SUBGRPS_INV);
}



//-----------------------------------------------------------
// Assumes, that l is an elkies-prime and
// that ftau[0] contains a root of the
// modular polynomial.
//

void eco_prime::compute_trace_elkies ()
{
  debug_handler ("eco_prime", "compute_trace_elkies()");
  
  ff_pol     fC; 	// the Elkies - polynomial
  ff_element E2a, E2b; // the curve E/C corresponding to l-isogeny
  
  lidia_size_t ev; // eigenvalue of Frobenius-endomorphism
  udigit       u_ev; // ev as udigit
  
  //* * *  compute the Elkies polynomial * * * *
  
#ifdef TIMING
  timer t;
  t.set_print_mode();
  t.start_timer();
#endif
  
  compute_divisor_of_division_polynomial (fC, E2a, E2b);
  
#ifdef TIMING
  t.stop_timer();
  if (info)
    std::cout << "\nComputing Elkies polynomial takes time : " << t << std::flush;
#endif
  
#ifdef DEBUG
  ff_polmod fCpm;
  ff_pol div;
  
  fCpm.build (fC);
  compute_psi (div, static_cast<lidia_size_t>(l), fCpm);
  assert(div.is_zero());
#endif
  
  //* * * compute the eigenvalue  * * * *
  
  if (!find_eigenvalue (ev, fC, 0, 1, false)) 
    {
      lidia_error_handler ("eco_prime::compute_trace_elkies(...)",
			   "No eigenvalue found.");
    }
  else 
    {
      u_ev = static_cast<udigit>(ev);
      
      if (test_y_coordinate() || test_x_coordinate_uniquely()) {
	c.set_capacity (1);
	
	c[0] = divide_mod (eco_prime::q, u_ev, this->l);
	c[0] = add_mod    (c[0], u_ev, this->l);
	if (info)
	  std::cout << "\n==> c = " << c[0] << " mod " << l << std::flush;
      }
      else {
	// test on x-coordinate ==> two possibilities
	c.set_capacity (2);
	
	c[0] = divide_mod (eco_prime::q, u_ev, this->l);
	c[0] = add_mod    (c[0], u_ev, this->l);
	if (c[0] == 0)
	  c.set_capacity(1);
	else
	  c[1] = this->l - c[0];
	
	if (info)
	  if (c[0] != 0)
	    std::cout << "\n==> c = +/- " << c[0] << " mod " << l << std::flush;
	  else
	    std::cout << "\n==> c = 0 mod " << l << std::flush;
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

void eco_prime::compute_c_for_d(lidia_size_t d)
{
	if (sp_type == NO_SUBGRP_INV) {
		base_vector< ff2_element > adl2;
		ff2_element              yl2;
		ff2_element              Q;
		long                     i, sz;
		int                      pos;
		udigit                   h;

		ff2_element::init_field (l);
		nearly_all_of_order (d, adl2);

		sz = adl2.size ();
		Q = ff1(q);

		for (i = 0; i < sz; i++) {
			add(yl2, adl2[i], ff2(ff1(1)));
			multiply(yl2, yl2, yl2);
			multiply(yl2, yl2, Q);
			divide(yl2, yl2, adl2[i]);
			yl2.as_udigit(h);
			h = sqrt_mod(h, l);
			if (!c.bin_search(h, pos))
				c.insert(h);

			if (h != 0) {
				h = negate_mod(h, l);
				if (!c.bin_search(h, pos))
					c.insert(h);
			}
		}
	} // end of sp_TYPE = (d...d)
	else {
		// sp_TYPE = (1 1 d...d)

		// c = +-  (x+1) * sqrt [ q / x ], x in GF(l), ord(x) = d

		base_vector< ff1_element > ad;
		ff1_element              y;
		ff1_element              Q;
		long                  i, sz;
		int pos;
		udigit h;


		ff1_element::set_characteristic (l);
		nearly_all_of_order (d, ad);

		sz = ad.size ();
		Q = ff1(q);

		for (i = 0; i < sz; i++) {
			invert (y, ad[i]);
			multiply (y, Q, y);
			sqrt (y, y);
			add  (ad[i], ad[i], ff1_element(1));
			multiply (y, y, ad[i]);
			h = y.as_udigit();
			if (!c.bin_search(h, pos))
				c.insert(h);

			if (h != 0) {
				negate (y, y);
				h = y.as_udigit();
				if (!c.bin_search(h, pos))
					c.insert(h);
			}
		}
	} // end of split type (11d...d)
}



void eco_prime::compute_trace_atkin ()
{
	if (!sp_degree_known) {
		lidia_size_t *d_list;

		c.set_size(0);
		compute_d_list(d_list);

		if (d_list[0] == 1)
			compute_c_for_d(d_list[1]);
		else {
			for (int i = 1; i <= d_list[0]; i++) {
				compute_c_for_d(d_list[i]);
			}
		}
	}
	else {
		// splitting degree known

		// Consider Splitting type is (1...1) or (1 l)

		if (sp_type == ALL_SUBGRPS_INV || sp_type == ONE_SUBGRP_INV) {
			//  c = +- 2 * sqrt (q) in GF(l)
			udigit h;
			c.set_capacity (2);
			h = sqrt_mod(q, l);
			h = add_mod(h, h, l); // h = 2 * sqrt (p^n mod l)
			c[0] = h;
			c[1] = negate_mod(h, l);
		}
		else
			if (sp_type == NO_SUBGRP_INV) {
				// Splitting type is (d...d), d = sp_degree
				// c = +-  (x+1) * sqrt [ q / x ], x in GF(l^2), ord(x) = d
				if (sp_degree == 2) {
					c.set_capacity (1);
					c[0] = 0;
				}
				else {
					c.set_size(0);
					compute_c_for_d (sp_degree);
				}
			}
			else
				if (sp_type == TWO_SUBGRPS_INV) {    // c = +-  (x+1) * sqrt [ q / x ], x in GF(l), ord(x) = d

					if (sp_degree == 2) {
						c.set_capacity (1);
						c[0] = 0;
					}
					else {
						c.set_size(0);
						compute_c_for_d (sp_degree);
					} // end of split type (11d...d)
				}
				else {
					lidia_error_handler ("eco_prime::compute_trace_atkin()",
							     "Unknown splitting type for Atkin prime.");
				}
	}
	if (info) {
		c.sort();
		std::cout << "\n==> c mod " << l << " in " << c << std::flush;
	}

	if (l <= schoof_bound)
		schoof_algorithm();
}



bool eco_prime::
probabilistic_correctness_proof(const bigint & res,
				lidia_size_t tests)
{
	return E_.probabilistic_test_of_group_order(res, tests);
}



bool eco_prime
::is_supersingular(bigint & res)
{
	unsigned int d = get_absolute_degree(A);
	bigint m;
	bigint p  = characteristic(A);
	bigint pn = number_of_elements(A);

	if ((d & 1) || (!(d&1) && ((p.least_significant_digit() & 3) == 3)))
		if (probabilistic_correctness_proof(pn+1)) {
			res.assign(pn+1);
			return true;
		}

	if (!(d & 1)) {
		power(m, p, d/2);
		if (probabilistic_correctness_proof(pn+1-m)) {
			res.assign(pn+1-m);
			return true;
		}
		if (probabilistic_correctness_proof(pn+1+m)) {
			res.assign(pn+1+m);
			return true;
		}

		m.multiply_by_2();
		if (probabilistic_correctness_proof(pn+1-m)) {
			res.assign(pn+1-m);
			return true;
		}
		if (probabilistic_correctness_proof(pn+1+m)) {
			res.assign(pn+1+m);
			return true;
		}
	}
	return false;
}



//
// check whether E is Fq-isogen to a curve of j-invariant 0 or
// 1728. If so, then res is set to the (probable) group order and
// true is returned, otherwise false is returned.
//

bool eco_prime
::check_j_0_1728 (bigint & res)

{
	bigint pn = number_of_elements(A);

	bigint qp1(pn+1);
	bigint x, y, cc;
	point< gf_element > P(E_), Q(E_);

	P = E_.random_point();

	if (cornacchia(x, y, -3, pn)) {
		multiply(Q, qp1+x, P);
		if (Q.is_zero())
			if (probabilistic_correctness_proof(qp1 + x)) {
				add(res, qp1, x);
				return true;
			}
		multiply(Q, qp1-x, P);
		if (Q.is_zero())
			if (probabilistic_correctness_proof(qp1 - x)) {
				subtract(res, qp1, x);
				return true;
			}

		multiply(cc, y, 3);
		subtract(cc, cc, x);
		cc.divide_by_2();

		multiply(Q, qp1+cc, P);
		if (Q.is_zero())
			if (probabilistic_correctness_proof(qp1 + cc)) {
				add(res, qp1, cc);
				return true;
			}
		multiply(Q, qp1-cc, P);
		if (Q.is_zero())
			if (probabilistic_correctness_proof(qp1 - cc)) {
				subtract(res, qp1, cc);
				return true;
			}

		multiply(cc, y, 3);
		add(cc, cc, x);
		cc.divide_by_2();

		multiply(Q, qp1+cc, P);
		if (Q.is_zero())
			if (probabilistic_correctness_proof(qp1 + cc)) {
				add(res, qp1, cc);
				return true;
			}
		multiply(Q, qp1-cc, P);
		if (Q.is_zero())
			if (probabilistic_correctness_proof(qp1 - cc)) {
				subtract(res, qp1, cc);
				return true;
			}
	}

	if (cornacchia(x, y, -4, pn)) {
		multiply(Q, qp1+x, P);
		if (Q.is_zero())
			if (probabilistic_correctness_proof(qp1 + x)) {
				add(res, qp1, x);
				return true;
			}
		multiply(Q, qp1-x, P);
		if (Q.is_zero())
			if (probabilistic_correctness_proof(qp1 - x)) {
				subtract(res, qp1, x);
				return true;
			}

		shift_left(cc, y, 1);

		multiply(Q, qp1+cc, P);
		if (Q.is_zero())
			if (probabilistic_correctness_proof(qp1 + cc)) {
				add(res, qp1, cc);
				return true;
			}
		multiply(Q, qp1-cc, P);
		if (Q.is_zero())
			if (probabilistic_correctness_proof(qp1 - cc)) {
				subtract(res, qp1, cc);
				return true;
			}
	}
	return false;
}



//
// determine the j-invariant of short weierstrass curve, should
// be replaced by function in elliptic_curve class.
//

void eco_prime
::compute_jinv (ff_element & jinv,
		const ff_element & A,
		const ff_element & B)
{
	ff_element t1, t2, t3;

	// t1 = 4 * A^3
	square  (t1, A);
	multiply (t1, t1, A);
	multiply (t1, t1, ff_element(4));

	// t2 = 27 * B^2
	square(t3, B);
	multiply (t2, t3, ff_element(27));

	// t2 = 4 * A^3 + 27 * B^2
	add (t2, t2, t1);

	// t1 = 6912 * A^3

	multiply (t1, t1, ff_element(1728));

	// jinv = t1 / t2

	if (t2.is_zero())
		lidia_error_handler ("eco_prime", "compute_jinv::singular curve");
	else
		divide (jinv, t1, t2);
}



//
//  This function determines Y^(q-1) mod (f, Y^2-X^3-aX-b) and sets
//  res to this polynomial.
//

void eco_prime::Ytop_f(ff_pol & res, const ff_polmod & f)

{
  bigint pm1h;
  ff_pol    h;
  ff_polmult H;
  lidia_size_t i, n;
  
  shift_right (pm1h, eco_prime::pn, 1); // pm1h = (q-1)/2

  // H = h = X^3 + AX + B
  
  CurveEqn(h, eco_prime::p, A, B, f);
  H.build (h, f); // poly_multiplier for faster arithmetic

  // compute res = H^((q-1)/2), left to right exponentiation
  
  n = pm1h.bit_length();
  res.assign_one();

  for (i = n - 1; i >= 0; i--) {
    square (res, res, f);
    if (pm1h.bit(i))
      multiply (res, res, H, f);
  }
}



//---------------------------------------------------------------
// this function returns the equation of the curve modulo
// input modulus (used in Elkies or Schoof algorithm) as
// polynomial.

void eco_prime::CurveEqn (ff_pol & pol, const bigint & p,
			  const ff_element & a, const ff_element & b,
			  const ff_polmod & f)
{
	ff_element tmp;

	pol.set_modulus(p);
	pol.assign_zero();
	pol.set_coefficient(3);
	tmp.assign(a);
	pol.set_coefficient(tmp.mantissa(), 1);
	tmp.assign(b);
	pol.set_coefficient(tmp.mantissa(), 0);

	if (f.modulus().degree() <= 3)
		remainder (pol, pol, f);
}



//----------------------------------------------------------
// This function determines all points Q with 2*Q = P and
// returns the number of such points. The result is returned
// in the array res of suitable size.
//

int eco_prime::extract_sqrt(point< gf_element > * & res,
                            const point< gf_element > & P)
{
	base_vector< bigint > root_vector;
	point< gf_element > Q(P), H(P);
	ff_element h;
	int i, j;

	meq_pol.set_modulus(p);
	meq_pol.assign_zero();
	meq_pol.set_coefficient(4);

	multiply(h, ff_element(-4), P.get_x().polynomial_rep().const_term());
	meq_pol.set_coefficient(h.mantissa(), 3);

	multiply(h, ff_element(-2), A);
	meq_pol.set_coefficient(h.mantissa(), 2);

	multiply(h, A, P.get_x().polynomial_rep().const_term());
	multiply(h, h, ff_element(4));
	add(h, h, ff_element(8) * B);
	negate(h, h);
	meq_pol.set_coefficient(h.mantissa(), 1);

	square(h, A);
	subtract(h, h, ff_element(4)* B * P.get_x().polynomial_rep().const_term());
	meq_pol.set_coefficient(h.mantissa(), 0);

	root_vector = find_roots(meq_pol);

	if (root_vector.size() == 0)
		return 0;


	gf_element hx(E_.get_a4().get_field());
	gf_element hy(E_.get_a4().get_field());

	j = 0;

	for (i = 0; i < root_vector.size(); i++) {
		square(h, root_vector[i]);
		add(h, h, A);
		multiply(h, h, root_vector[i]);
		add(h, h, B);
		if (is_square(h)) {
			h = sqrt(h);

			hx.assign(root_vector[i]);
			hy.assign(h.mantissa());
			H.assign(hx, hy);

			multiply_by_2(Q, H);
			if (Q == P)
				res[j++] = H;
			negate(Q, Q);
			if (Q == P)
				res[j++] = -H;
		}
	}
	return j;
}



int eco_prime::max_two_power (const point< gf_element > & P)
{
	point< gf_element > * T;
	T = new point< gf_element > [4];

	if (extract_sqrt(T, P) ==  0) {
		delete[] T;
		return 1;
	}
	else {
		int l0, l1, l2;

		l0 = max_two_power(T[0]);
		l1 = max_two_power(T[1]);
		l2 = max_two_power(T[2]);

		delete[] T;

		if (l0 < l1)
			l0 = l1;
		return (l0 < l2 ? 2*l2 : 2*l0);
	}
}



void eco_prime::compute_mod_2_power()
{
	ff_pol res, f_x;

	// f_x = x mod p
	f_x.set_modulus(p);
	f_x.assign_zero();
	f_x.set_coefficient(1);

	// meq_pol = X^3 + aX + b mod p
	meq_pol.set_modulus(p);
	meq_pol.assign_zero();
	meq_pol.set_coefficient(3);
	meq_pol.set_coefficient(A.mantissa(), 1);
	meq_pol.set_coefficient(B.mantissa(), 0);

	// poly_modulus for meq_pol

	meq_pol_mod.build(meq_pol);

	// res = gcd(x^p - x, x^3 + ax + b)
	power_x(res, pn, meq_pol_mod);
	res = gcd(res-f_x, meq_pol);

	// no root in Fp, no 2-torsion point in E(Fp)
	if (res.is_one()) {
		l = 2;
		c[0] = 1;
		if (info)
			std::cout << "\n==> c = 1 mod 2 " << std::flush;
	}

	// at least one 2-torsion point exists in E(Fp)

	else if (res.degree() == 1) {
		point< gf_element > P(E_);
		point< gf_element > *T;

		l = 2;
		T = new point< gf_element > [4];

		// P.assign(find_root(res), ff_element(0));
		//
		gf_element hx(E_.get_a4().get_field());
		gf_element hy(E_.get_a4().get_field());
		hx.assign(find_root(res));
		hy.assign_zero();
		P.assign(hx, hy);

		while (extract_sqrt(T, P) > 0) {
			l <<= 1;
			P = T[0];
		}
		delete [] T;
		c[0] = (pn.least_significant_digit() - l + 1) % (2*l);
		l <<= 1;

		if (info)
			std::cout << "\n==> c = " << c[0] << " mod " << l << std::flush;
	}
	else {
		// three 2-TP exist, E[2] = <P1> x <P2>
		udigit l1 = 2, l2 = 2, l3 = 2;
		//elliptic_curve< ff_element > e(A, B);
		point< gf_element > P1(E_), P2(E_), P3(E_);

		base_vector< bigint > roots = find_roots(res);

		gf_element hx(E_.get_a4().get_field());
		gf_element hy(E_.get_a4().get_field());

		// P1.assign(roots[0], ff_element(0));
		//
		hx.assign(roots[0]);
		hy.assign_zero();
		P1.assign(hx, hy);

		// P2.assign(roots[1], ff_element(0));
		//
		hx.assign(roots[1]);
		hy.assign_zero();
		P2.assign(hx, hy);

		// P3 = P1 + P2
		//
		add(P3, P1, P2);

		l1 = 2*max_two_power(P1);
		l2 = 2*max_two_power(P2);
		l3 = 2*max_two_power(P3);

		if (l1 == l3)
			l = l1 * l2;
		else
			l = l1 * l3;
		c[0] = (pn.least_significant_digit() - l + 1) % (2*l);
		l <<= 1;

		if (info)
			std::cout << "\n==> c = " << c[0] << " mod " << l << std::flush;
	}
}



bigint eco_prime::compute_group_order(int info)
{
	//
	// initialize timer
	//

	timer t;
	t.set_print_mode (HMS_MODE);

	//
	// Initialize curve, static variables of eco_prime,
	// trace_mod, and trace_list< ff_element >
	//
	set_info_mode(info);
	trace_list::set_info_mode(info);
	trace_list tl;
	trace_mod tm;
	udigit ll;
	bigint ec_order;

	tl.set_curve(E_);

	//
	// Test supersingularity and isogeny to curves with
	// j-invariant 0 or 1728.
	//

	if (info) {
		t.start_timer();
		std::cout << "\n\n-------------------------------------------------------\n";
		std::cout << "Checking supersingularity ... " << std::flush;
	}
	if (is_supersingular(ec_order)) {
		if (info) {
			std::cout << " ==> Test successful \n" << std::flush;
			std::cout << "\n\nField Size : "; std::cout << pn;
			std::cout << "\nCoefficient a4 : "; std::cout << A;
			std::cout << "\nCoefficient a6 : "; std::cout << B;
			std::cout << "\n\nGroup order is " << ec_order << "." << std::endl;
		}
		return ec_order;
	}
	else if (info)
		std::cout << " ==> Test failed \n" << std::flush;

	if (info)
		std::cout << "Checking isogeny to j=0 or j=1728 curves ... " << std::flush;

	if (check_j_0_1728(ec_order)) {
		if (info) {
			std::cout << " ==> Test successful \n" << std::flush;
			std::cout << "\n\nField Size : "; std::cout << pn;
			std::cout << "\nCoefficient a4 : "; std::cout << A;
			std::cout << "\nCoefficient a6 : "; std::cout << B;
			std::cout << "\n\nGroup order is " << ec_order << "." << std::endl;
		}
		return ec_order;
	}
	else if (info)
		std::cout << " ==> Test failed \n" << std::flush;

	if (info)
		std::cout << "\n\n-------------------------------------------------------\n";


	//
	// compute trace mod power of 2
	//

	if (info)
		std::cout << "Working on prime l = 2\n" << std::flush;
	compute_mod_2_power();
	tm.set_vector(get_prime(), get_relation());
	tl.append(tm);


	//
	// compute trace mod primes and store the residues
	// in the trace_list
	//

	bool end_loop;

	ll = 3;
	end_loop = false;

 	 while ((!end_loop) && (ll <= meq_prime::MAX_MEQ_PRIME))
	 {
		if (info) {
			std::cout << "\n\n-------------------------------------------------------\n";
			std::cout << "Working on prime l = " << ll;
		}

                if (p % ll == 0)
                  {
		  if (info) 
			std::cout << " ... skipped since equal to field characteristic.\n";
		  ll = next_prime(ll+1);
                  continue;
                  }  

		if (set_prime(ll)) {
			if (info)
				std::cout << "\n" << std::flush;

			compute_splitting_type();

			if (is_elkies())
				compute_trace_elkies();
			else
				compute_trace_atkin();

			tm.set_vector(get_prime(), get_relation());
			end_loop = tl.append(tm);
		}
		else if (info)
			std::cout << " ... not available!" << std::endl;

		ll = next_prime(ll+1);
	}

	if (ll > meq_prime::MAX_MEQ_PRIME) {
		lidia_error_handler ("eco_prime::compute_group_order",
				     "Not enough primes available.");
	}


	//
	// compute the group order with babystep/giantstep algorithm
	//

	if (info) {
		std::cout << "\n-------------------------------------------------------\n";
		std::cout << "\nStarting Babystep Giantstep Phase";
		std::cout << "\n\nFinal TraceList : "<< tl << "\n" << std::flush;
	}

	ec_order = tl.bg_search_for_order();

	if (info) {
		t.stop_timer();
		std::cout << "\n\nOutput of BG algorithm = " << ec_order << std::flush;
		std::cout << "\nComputing Time : " << t << std::flush;
	}


	//
	// Verification of the group order, if desired
	//

	if (use_corr_test)
		if (probabilistic_correctness_proof(ec_order)) {
			if (info) {
				std::cout << "\n\nProbabilistic correctness proof agrees !!";
				std::cout << "\n\nField Size : "; std::cout << pn;
				std::cout << "\nCoefficient a4 : "; std::cout << A;
				std::cout << "\nCoefficient a6 : "; std::cout << B;
				std::cout << "\n\nGroup order is " << ec_order << "." << std::endl;
			}
		}
		else {
			std::cout << "\n\nField Size : "; std::cout << pn;
			std::cout << "\nCoefficient a4 : "; std::cout << A;
			std::cout << "\nCoefficient a6 : "; std::cout << B;
			std::cout << "\n\nOrder candidate : " << ec_order << "." << std::endl;
			lidia_error_handler("eco_prime",
					    "Probabilistic correctness rejects candidate !!");
		}
	return ec_order;
}



//--------------------------------------------------------------------
// finally a general function which can be used generally.

bigint compute_group_order(const gf_element & a4, const gf_element & a6,
			   int info)
{
	eco_prime ep(a4, a6);
	return ep.compute_group_order(info);
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
