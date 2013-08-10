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
//	Author	: Volker Mueller (VM), Markus Maurer (MM), Andrea Rau (AR)
//	Changes	: See CVS log
//
//==============================================================================================



#ifdef HAVE_CONFIG_H
# include	"config.h"
#endif
#include	"LiDIA/rational_factorization.h"
#include	"LiDIA/eco_prime.h"
#include	"LiDIA/trace_list.h"
#include	"LiDIA/EC_domain_parameters_P1363.h"
#include	"LiDIA/timer.h"

#include        <cstdlib>  // declares exit()


#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



int EC_domain_parameters_P1363::GFP = 1;
int EC_domain_parameters_P1363::GF2N = 1;


lidia_size_t EC_domain_parameters_P1363::defaultBitsize_ = 160;
lidia_size_t EC_domain_parameters_P1363::defaultPercentage_ = 10;



//
// Constructor / destructor
//

EC_domain_parameters_P1363::EC_domain_parameters_P1363()
{
	initialized = false;
}



EC_domain_parameters_P1363::~EC_domain_parameters_P1363()
{
}



//
// assignments and access
//

lidia_size_t EC_domain_parameters_P1363::default_bitsize()
{
	return EC_domain_parameters_P1363::defaultBitsize_;
}



lidia_size_t EC_domain_parameters_P1363::default_percentage()
{
	return EC_domain_parameters_P1363::defaultPercentage_;
}



void EC_domain_parameters_P1363::assign(const EC_domain_parameters_P1363 & I)
{

	if (this != &I) {
		initialized = I.initialized;

		if (initialized) {
			a = I.a;
			b = I.b;
			q = I.q;
			r = I.r;
			k = I.k;
			G = I.G;
		}
	}
}



EC_domain_parameters_P1363 &
EC_domain_parameters_P1363::operator = (const EC_domain_parameters_P1363 & I)
{
	this->assign(I);
	return *this;
}



const bigint & EC_domain_parameters_P1363
::get_q () const
{
	if (!initialized)
		lidia_error_handler("EC_domain_parameters_P1363", "get_q::not initialized");

	return q;
}



const gf_element & EC_domain_parameters_P1363
::get_a () const
{
	if (!initialized)
		lidia_error_handler("EC_domain_parameters_P1363", "get_a::not initialized");

	return a;
}



const gf_element & EC_domain_parameters_P1363
::get_b () const
{
	if (!initialized)
		lidia_error_handler("EC_domain_parameters_P1363", "get_b::not initialized");

	return b;
}



const bigint & EC_domain_parameters_P1363
::get_k () const
{
	if (!initialized)
		lidia_error_handler("EC_domain_parameters_P1363", "get_k::not initialized");

	return k;
}



const bigint & EC_domain_parameters_P1363
::get_r () const
{
	if (!initialized)
		lidia_error_handler("EC_domain_parameters_P1363", "get_r::not initialized");

	return r;
}



const point< gf_element > & EC_domain_parameters_P1363
::get_G () const
{
	if (!initialized)
		lidia_error_handler("EC_domain_parameters_P1363", "get_G::not initialized");

	return G;
}



void EC_domain_parameters_P1363
::get_twist_coeff(gf_element & new_a, gf_element & new_b)
{
	gf_element d(a);
	do {
		d.randomize();
	} while (jacobi(d.polynomial_rep().const_term(), d.characteristic()) != -1);

	new_a = a*d*d;
	new_b = b*d*d*d;
}



//
// Curve generation
//

void EC_domain_parameters_P1363
::generate_parameters(int field, int info)
{
	this->generate_parameters(field, this->default_bitsize(),
				  this->default_percentage(), info);
}



void EC_domain_parameters_P1363
::generate_parameters(int field, int bitsize_factor, int info)
{
	this->generate_parameters(field, bitsize_factor,
				  this->default_percentage(), info);
}



//--------------------------------------------------------------------
//


void EC_domain_parameters_P1363::generate_parameters(int field,
						     int bitsize_factor,
						     int percentage,
						     int info)
{
	bigint p, co_factor, found_co_factor, found_co_factor_tw, order;
	int tries = 1;
	rational_factorization rf_order;

	// here we need something like
	// if(field != GFP)
	// lidia_error_handler("not yet implemented");

	// Initialize timer
	//
	timer t; t.set_print_mode (HMS_MODE);
	timer t2; t2.set_print_mode(HMS_MODE);

	// Initialize counters
	//
	int both_abort = 0;
	int no_abort   = 0;

	// computation control
	//
	bool good_curve_Etw;
	bool good_curve_E;
	bool stop_computation;

	shift_left(p, bigint(1), bitsize_factor - 1);

	// determine cofactor = 2^{n}, n = "percentage" percent of bitsize_factor.
	//
	if (percentage >= 100)
		co_factor = p;
	else {
		co_factor.assign_one();
		if (percentage != 0 && (static_cast<int>((static_cast<double>(bitsize_factor) / 100.0)
							* static_cast<double>(percentage)) - 1 > 0))
			shift_left(co_factor, co_factor,
				   static_cast<int>((static_cast<double>(bitsize_factor)/100.0) *
						    static_cast<double>(percentage)) - 1);

		if (co_factor > p)
			co_factor = p;
	}

	// determine characteristic p of prime field
	//
	multiply(p, p, co_factor);
	p = next_prime(p-1);

	// The field
	//
	galois_field theField(p);
	bigmod::set_modulus(p);

	// Twist coordinates
	//
	gf_element a4tw(theField), a6tw(theField);


	// initialize computation
	//
	elliptic_curve< gf_element > e;
	trace_list tl;
	trace_mod tm;
	udigit ll;
	long qq;
	eco_prime ep;
	bigint g;

	if (info >= MUCH_INFO) {
		ep.set_info_mode(1);
		trace_list::set_info_mode(1);
	}
	else {
		ep.set_info_mode(0);
		trace_list::set_info_mode(0);
	}
	rf_order.verbose(0);
	t.start_timer();


	//
	// Choose curves until strong one is found.
	// Use early abort strategy on curve and its twist.
	//

	std::cout << "\nCharacteristic : q = " << p;
	std::cout << " (" << p.bit_length() << " bits)" << std::endl;
	std::cout << "Upper bound for the cofactor is " << co_factor;
	std::cout << " (" << co_factor.bit_length() << " bits)" << std::endl;
	std::cout << std::endl;

	do {
		// Choose next curve
		//
		a.assign_zero(theField);
		b.assign_zero(theField);
		a.randomize();
		b.randomize();

		if (info >= LITTLE_INFO) {
			std::cout << "\n###############################################################";
			std::cout << "\nChoosing random E : \n\na = " << a << "\nb = " << b << std::endl;
		}

		// initialize for new curve
		//
		tries ++;
		found_co_factor = 1;
		found_co_factor_tw = 1;

		good_curve_E   = true;
		good_curve_Etw = true;
		stop_computation = false;

		e.set_coefficients(a, b);
		ep.set_curve(a, b);
		t2.start_timer();

		// Test supersingularity
		//
		if (info >= MUCH_INFO)
			std::cout << "\nTesting E for supersingularity ... " << std::flush;

		if (ep.is_supersingular(order)) {
			if (info >= MUCH_INFO)
				std::cout << "YES" << std::endl;
			continue;
		}
		else
			if (info >= MUCH_INFO)
				std::cout << "NO" << std::endl;

		// Test for isogeny to curve with j-invariant 0 or 1728
		//

		if (info >= MUCH_INFO)
			std::cout << "Testing whether E is Fq-isogenous to curve with "
				"j-invariant 0 or 1728 ... " << std::flush;

		if (ep.check_j_0_1728 (order)) {
			if (info >= MUCH_INFO)
				std::cout << "YES" << std::endl;
			continue;
		}
		else
			if (info >= MUCH_INFO)
				std::cout << "NO" << std::endl;


		// Determine coefficients of twist.
		//
		get_twist_coeff(a4tw, a6tw);

		// Test twist for supersingularity
		//
		if (info >= MUCH_INFO)
			std::cout << "Testing twist(E) for supersingularity ... " << std::flush;

		if (ep.is_supersingular(order)) {
			if (info >= MUCH_INFO)
				std::cout << "YES" << std::endl;
			continue;
		}
		else
			if (info >= MUCH_INFO)
				std::cout << "NO" << std::endl;


		// Test twist for isogeny to curve with j-invariant 0 or 1728
		//
		if (info >= MUCH_INFO)
			std::cout << "Testing whether twist(E) is Fq-isogenous to curve "
				"with j-invariant 0 or 1728 ... " << std::flush;

		if (ep.check_j_0_1728 (order)) {
			if (info >= MUCH_INFO)
				std::cout << "YES" << std::endl;
			continue;
		}
		else
			if (info >= MUCH_INFO)
				std::cout << "NO" << std::endl;


		// Compute trace modulo 2
		//
		tl.clear();
		tl.set_curve(e);

		// trace mod 2
		//
		if (info >= MUCH_INFO) {
			std::cout << std::endl;
			std::cout << "-------------------------------------------------------" << std::endl;
			std::cout << "Working on prime l = 2" << std::endl;
		}

		ep.compute_mod_2_power();
		tm.set_vector(ep.get_prime(), ep.get_relation());
		tl.append(tm);

		qq = remainder(p, ep.get_prime());
		g = gcd(p + 1 - (ep.get_relation())[0], ep.get_prime());
		if (! g.is_one())
			multiply(found_co_factor, found_co_factor, g);
		g = gcd(p + 1 + (ep.get_relation())[0], ep.get_prime());
		if (! g.is_one())
			multiply(found_co_factor_tw, found_co_factor_tw, g);


		// for each prime ll compute trace mod ll
		// and check for divisor of group order
		//
		ll = 3;

		do {
			if (found_co_factor > co_factor && found_co_factor_tw > co_factor) {
				t2.stop_timer();
				if (info >= LITTLE_INFO) {
					if (info >= MUCH_INFO)
						std::cout << "\n";
					std::cout << "\nTest No. " << tries-1;
					std::cout <<" (early abort for E and twist(E)) took time " << t2 << std::endl;
				}
				both_abort ++;
				good_curve_E = good_curve_Etw = false;
				stop_computation = true;
			}
			else {
				if (info >= MUCH_INFO) {
					std::cout << std::endl;
					std::cout << "-------------------------------------------------------" << std::endl;
					std::cout << "Working on prime l = " << ll <<std::endl;
				}

				if (ep.set_prime(ll)) {
					ep.compute_splitting_type();

					if (ep.is_elkies()) {
						ep.compute_trace_elkies();

						qq = remainder(p, ep.get_prime());
						g = gcd(p + 1 - static_cast<long>(ep.get_relation()[0]), ep.get_prime());
						if (! g.is_one())
							multiply(found_co_factor, found_co_factor, g);
						g = gcd(p + 1 + static_cast<long>(ep.get_relation()[0]), ep.get_prime());
						if (! g.is_one())
							multiply(found_co_factor_tw, found_co_factor_tw, g);
					}
					else
						ep.compute_trace_atkin();

					tm.set_vector(ep.get_prime(), ep.get_relation());
					stop_computation = tl.append(tm);
				}
				else
					if (ll >= 1020)
						lidia_error_handler("EC_domain_parameters_P1363", "Not "
								    "enough modular equations");
					else
						lidia_warning_handler("EC_domain_parameters_P1363",
								      "Modular equation could not be"
								      " found, ignored ...");

				ll = next_prime(ll+1);
			}
		} while (!stop_computation);

		// Early abort
		//
		if (!good_curve_E && !good_curve_Etw)
			continue;


		// Determine group order of E
		//

		if (info >= MUCH_INFO)
			std::cout << "\n";
		order = tl.bg_search_for_order();
		no_abort ++;
		t2.stop_timer();
		if (info >= LITTLE_INFO) {
			std::cout << "\nTest No. " << tries-1 << " took time " << t2 <<"\n"<< std::endl;
		}

		// Test for strong curve E
		//
		if (found_co_factor <= co_factor) {
		        t2.start_timer();
			rf_order.assign(order);
			good_curve_E = this->is_strong_curve(rf_order, order, co_factor,
							     p, info);
			t2.stop_timer();
			if (info >= LITTLE_INFO) {
			    std::cout << "\nTest for strong curve took time "
				      << t2 << std::endl;
			}
		}
		else
			good_curve_E = false;


		// Test for strong twist Etw
		//
		if (!good_curve_E) {
			if (found_co_factor_tw <= co_factor) {
   			        t2.start_timer();
				order = p + 1 + (p+1-order);
				rf_order.assign(order);

				good_curve_Etw = this->is_strong_curve(rf_order, order,
								       co_factor, p, info);
			        t2.stop_timer();
				if (info >= LITTLE_INFO) {
				    std::cout << "\nTest for strong twist "
					      << "took time "
					      << t2 << std::endl;
				}

				if (good_curve_Etw) {
					a = a4tw;
					b = a6tw;
				}
			}
			else
				good_curve_Etw = false;
		}
	} while (!good_curve_E && !good_curve_Etw);
	t.stop_timer();


	// Verification of the group order
	//
	e.set_coefficients(a, b);

	if (e.probabilistic_test_of_group_order(order)) {
		if (info >= LITTLE_INFO) {
			std::cout << "\n\nFOUND GOOD CURVE : ";
			std::cout << "\nCharacteristic : " << p << " (" << p.bit_length();
			std::cout << " bits)" << std::flush;
			std::cout << "\nCoefficient a : "; std::cout << a;
			std::cout << "\nCoefficient b : "; std::cout << b;
			std::cout << "\n\nGroup order is " << order << std::endl;
			std::cout << "             = " << rf_order << std::endl;
		}
		// Initialize P1363 parameters
		//
		q = p;
		r = rf_order.base(rf_order.no_of_comp()-1);
		k = order / r;

		// search point G of order r
		//

		point< gf_element > P;
		do {
			P = e.random_point();
			multiply(G, k, P);
		} while (G.is_zero());

		if (info >= LITTLE_INFO) {
			std::cout << "\nPrime Factor r of group order is : " << r << std::endl;
			std::cout << "Point G of order r : " << G << std::endl;
			std::cout << "Cofactor k is : " << k << std::endl;
			std::cout << "\n\nComplete Test took " << tries-1 << " tries and "
				"time " << t << "\n";
			std::cout << "\nearly abort for both curves was used in ";
			std::cout << both_abort << " tests ";
			std::cout << "\nfull group order comp. was used in "<< no_abort << " tests \n\n";
		}

		initialized = true;
	}
	else {
		std::cout << "\n\nERROR: Probabilistic correctness test rejects candidate!!\n";
		std::exit(1);
	}
}



bool EC_domain_parameters_P1363::is_strong_curve(rational_factorization &
						 rf_order,
                                                 const bigint & order,
						 const bigint & co_factor,
                                                 const bigint & p,
                                                 int info) const
{
	bool good_curve = true;


	if (info >= MUCH_INFO)
		std::cout << "Checking whether small cofactor is too large ... " << std::flush;

	if (!rf_order.is_prime_factorization()) {
		rf_order.trialdiv();
		rf_order.ecm();
	}

	if (!rf_order.is_prime_factorization()) {
		good_curve = false;
		if (info >= MUCH_INFO)
			std::cout << " --> can not factor order." << std::endl;
	}
	else
		if (order/rf_order.base(rf_order.no_of_comp()-1) <= co_factor) {
			if (info >= MUCH_INFO)
				std::cout << "NO" << std::endl;
			good_curve = true;
		}
		else {
			if (info >= MUCH_INFO) {
				std::cout << "YES (";
				std::cout << order/rf_order.base(rf_order.no_of_comp()-1) << ")" << std::endl;
			}
			good_curve = false;
		}

	if (good_curve) {
		if (info >= MUCH_INFO)
			std::cout << "\nTesting for anomalous curve ... " << std::flush;

		if (order == p) {
			if (info >= MUCH_INFO)
				std::cout << "YES" << std::endl;
			good_curve = false;
		}
		else
			if (info >= MUCH_INFO)
				std::cout << "NO" << std::endl;
	}

	// Test that order does not divide (p^B)-1 for B small
	// (Reduction to finite fields attack).
	//

	if (good_curve) {
		int i;
		int log2 = p.bit_length() - 1;
		int stop = static_cast<int>(std::floor(2000.0 / log2));

		if (info >= MUCH_INFO) {
			std::cout << "Test whether order divides q^k-1 for some 1 <= k <= ";
			std::cout << stop << " ... " << std::flush;
		}

		multi_bigmod p_mod(p, order);
		multi_bigmod ppower_mod(1, order);

		for (i = 1; i <= stop; i++) {
			multiply(ppower_mod, ppower_mod, p_mod);
			if (ppower_mod.is_one()) {
				if (info >= MUCH_INFO)
					std::cout << "YES (k = " << i << ")" << std::endl;
				good_curve = false;
				i = stop+1;
			}
		}

		if (good_curve && info >= MUCH_INFO)
			std::cout << "NO" << std::endl;
	}

	return good_curve;
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
