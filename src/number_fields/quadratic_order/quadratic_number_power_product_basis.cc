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
//	Author	: Markus Maurer (MM)
//	Changes	: See CVS log
//
//==============================================================================================


//#define PPBQN_WATCH_RECOMP

#ifdef HAVE_CONFIG_H
# include	"config.h"
#endif
#include	"LiDIA/quadratic_number_power_product_basis.h"
#include "LiDIA/precondition_error.h"



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



xbigfloat quadratic_number_power_product_basis::dummy;



//
//  constructor
//

quadratic_number_power_product_basis::quadratic_number_power_product_basis ()
{
	debug_handler ("quadratic_number_power_product_basis",
		       "quadratic_number_power_product_basis()");

	basis.set_mode(vector_flags::expand);
	ln_basis.set_mode(vector_flags::expand);
	prec_ln_basis.set_mode(vector_flags::expand);
	init.set_mode(vector_flags::expand);
}



//
// destructor
//

quadratic_number_power_product_basis::~quadratic_number_power_product_basis ()
{
	debug_handler ("quadratic_number_power_product_basis",
		       "~quadratic_number_power_product_basis()");
}



//
// Assignment
//

quadratic_number_power_product_basis &
quadratic_number_power_product_basis::operator = (const quadratic_number_power_product_basis & b)
{
	debug_handler ("quadratic_number_power_product_basis",
		       "operator = (const qnppb&)");

	this->assign(b);
	return *this;
}



void
quadratic_number_power_product_basis::assign (const quadratic_number_power_product_basis & b)
{
	debug_handler ("quadratic_number_power_product_basis",
		       "assign(const qnppb&)");

	if (this != &b) {
		basis.assign(b.basis);
		ln_basis.assign(b.ln_basis);
		prec_ln_basis.assign(b.prec_ln_basis);
		init.assign(b.init);
	}
}



void
quadratic_number_power_product_basis::set_basis (const base_vector< quadratic_number_standard > & b)
{
	debug_handler ("quadratic_number_power_product_basis",
		       "set_basis(const base_vector< quadratic_number_standard > &)");

	basis.assign(b);
	basis.set_mode(vector_flags::expand);
	ln_basis.set_capacity(0);
	prec_ln_basis.set_capacity(0);
	init.set_capacity(0);
}



void
quadratic_number_power_product_basis::set_basis (const quadratic_number_standard & b)
{
	debug_handler ("quadratic_number_power_product_basis",
		       "set_basis(const quadratic_number_standard &)");

	basis.set_capacity(1);
	basis[0].assign(b);
	ln_basis.set_capacity(0);
	prec_ln_basis.set_capacity(0);
	init.set_capacity(0);
}



void
quadratic_number_power_product_basis::concat (const quadratic_number_power_product_basis & b,
					      const quadratic_number_standard & q)
{
	debug_handler ("quadratic_number_power_product_basis",
		       "concat(const qnppb &, const quadratic_number_standard &)");

	lidia_size_t sz = b.basis.get_size()+1;

	if (this != &b) {
		// avoid double allocation for *this
		if (basis.get_capacity() < sz) {
			basis.set_size(sz);
			ln_basis.set_size(sz);
			prec_ln_basis.set_size(sz);
			init.set_size(sz);
		}

		basis.assign(b.basis);
		ln_basis.assign(b.ln_basis);
		prec_ln_basis.assign(b.prec_ln_basis);
		init.assign(b.init);
	}

	basis[b.basis.get_size()] = q;
	if (init.get_size() > 0)
		init[b.basis.get_size()] = 0;
}



void
quadratic_number_power_product_basis::concat (const quadratic_number_power_product_basis & b,
					      const quadratic_number_power_product_basis & c)
{
	debug_handler ("quadratic_number_power_product_basis",
		       "concat (const qnppb &, const qqnpp&)");

	basis.concat(b.basis, c.basis);
	ln_basis.concat(b.ln_basis, c.ln_basis);
	prec_ln_basis.concat(b.prec_ln_basis, c.prec_ln_basis);
	init.concat(b.init, c.init);
}



//
// Access
//

const base_vector< quadratic_number_standard > &
quadratic_number_power_product_basis::get_basis () const
{
	debug_handler ("quadratic_number_power_product_basis",
		       "get_basis() const");
	return basis;
}



const quadratic_number_standard &
quadratic_number_power_product_basis::operator [] (lidia_size_t i) const
{
	debug_handler ("quadratic_number_power_product_basis",
		       "operator[](lidia_size_t) const");
	return basis[i];
}



quadratic_number_standard &
quadratic_number_power_product_basis::operator [] (lidia_size_t i)
{
	debug_handler ("quadratic_number_power_product_basis",
		       "operator[](lidia_size_t)");

	if (i < init.get_size())
		init[i] = 0;
	return basis[i];
}



lidia_size_t
quadratic_number_power_product_basis::get_size () const
{
	debug_handler ("quadratic_number_power_product_basis",
		       "get_size() const");
	return basis.get_size();
}



const quadratic_order &
quadratic_number_power_product_basis::get_order () const
{
	debug_handler ("quadratic_number_power_product_basis",
		       "get_order() const");

	if (basis.get_size() == 0) {
		lidia_error_handler("quadratic_number_power_product_basis::get_order()",
				    "power product_basis not initialized.");
	}
	return (basis[0]).get_order();
}



//
// Comparison
//

bool
operator == (const quadratic_number_power_product_basis & a,
	     const quadratic_number_power_product_basis & b)
{
	debug_handler ("quadratic_number_power_product_basis",
		       "operator == (const qnppb&, const qnppb&)");
	if (&a == &b)
		return true;
	else
		return false;
}



//
// Ln approximations
//

const xbigfloat &
quadratic_number_power_product_basis::get_absolute_Ln_approximation (long k, lidia_size_t i) const
{
	debug_handler ("quadratic_number_power_product_basis",
		       "get_absolute_Ln_approximation(long, lidia_size_t)const");

#ifdef PPBQN_WATCH_RECOMP
	std::cout << "quadratic_number_power_product_basis:: ";
	std::cout << "(index:" << i << ", prec:" << k << ") " << std::flush;
#endif

	// verify index range
	//
	if (i< 0 || i >= basis.get_size()) {
		precondition_error_handler(i, "i", "0 <= i < size",
				    basis.get_size(), "size", "",
				    "quadratic_number_power_product_basis::"
				    "get_absolute_Ln_approximation",
				    "quadratic_number_power_product_basis",
				    "Index out of range.");
#ifdef PPBQN_WATCH_RECOMP
		std::cout << std::endl;
#endif
		return quadratic_number_power_product_basis::dummy;
	}

	// check for trivial case
	//
	if (basis[i].is_one()) {
		quadratic_number_power_product_basis::dummy.assign_zero();
#ifdef PPBQN_WATCH_RECOMP
		std::cout << "is one" << std::endl;
#endif
		return quadratic_number_power_product_basis::dummy;
	}

	// recompute, if uninitialized or not accurate enough
	//
	if (init.get_size() != basis.get_size()) {
		lidia_size_t j, n;

		n = basis.get_size();

		ln_basis.set_capacity(n);
		prec_ln_basis.set_capacity(n);
		init.set_capacity(n);

		for (j = 0; j < n; j++)
			init[j] = 0;
	}

	if ((!init[i]) || (k > prec_ln_basis[i])) {
		if (k < -1) k = -1;
		basis[i].get_absolute_Ln_approximation(ln_basis[i], k);
		prec_ln_basis[i] = k;
		init[i] = 1;

#ifdef PPBQN_WATCH_RECOMP
		std::cout << "recomputation" << std::endl;
#endif
		return ln_basis[i];
	}

	// truncate, if more accurate than needed
	//
	else if (k < prec_ln_basis[i]) {
		// k += 1+b_value(log_q[i])

		if (check_overflow(k, k, 1)) {
			lidia_error_handler ("quadratic_number_power_product_basis",
					     "::get_absolute_Ln_approximation",
					     "precision overflow.");

			return quadratic_number_power_product_basis::dummy;
		}

		if (check_overflow(k, k, ln_basis[i].b_value())) {
			lidia_error_handler ("quadratic_number_power_product_basis",
					     "::get_absolute_Ln_approximation",
					     "precision overflow.");

			return quadratic_number_power_product_basis::dummy;
		}

		truncate (quadratic_number_power_product_basis::dummy, ln_basis[i], k);

#ifdef PPBQN_WATCH_RECOMP
		std::cout << std::endl;
#endif
		return quadratic_number_power_product_basis::dummy;
	}

	// reuse, if required precision matches exact.
	//
	else {
#ifdef PPBQN_WATCH_RECOMP
		std::cout << std::endl;
#endif
		return ln_basis[i];
	}
}


//
//
//

void
quadratic_number_power_product_basis::absolute_Ln_approximations (const matrix< bigint > & M, long k)

	//
	// For each i, determine an absolute prec_ln_basis[i] approximation
	// ln_basis[i] to Ln(basis[i]), such that sum_i e_ij ln_basis[i] is an
	// absolute k approximation to sum_i e_ij Ln(basis[i]) for each j, i.e.,
	//
	//  prec_ln_basis[i] = k + 2 + b(nofr) + max(j) b(e_ij).
	//
	// M = (e_i,j))
	// nofr = M.get_no_of_rows
	// nofc = M.get_no_of_columns
	//
	// 0 <= i < nofr,  0 <= j < nofc,
	//

{
	debug_handler ("quadratic_number_power_product_basis",
		       "absolute_Ln_approximations (const matrix< bigint > &, long)");

	lidia_size_t nofr, nofc, i, j;
	long t, kbn2;
	bigint max_row_entry;

	// initialize
	//
	nofr = M.get_no_of_rows();
	nofc = M.get_no_of_columns();

	if (nofr != basis.get_size()) {
		lidia_error_handler ("quadratic_number_power_product_basis::"
				     "absolute_Ln_approximations(matrix< bigint >, long)",
				     "Number of exponents does not match size of "
				     "basis.");
		return;
	}

	// If no precomputations yet, initialize arrays.
	//
	if (init.get_size() != basis.get_size()) {
		ln_basis.set_capacity(nofr);
		prec_ln_basis.set_capacity(nofr);
		init.set_capacity(nofr);

		for (j = 0; j < nofr; j++)
			init[j] = 0;
	}


	// Determine the precision kbn2
	//
	// kbn2 = k+b_value(nofr)+2
	//
	if (check_overflow(kbn2, k, 2)) {
		lidia_error_handler (
			"quadratic_number_power_product_basis::"
			"absolute_Ln_approximations(matrix< bigint >, long)",
			"Precision overflow. (1)");
		return;
	}

	if (check_overflow(kbn2, kbn2, b_value(nofr))) {
		lidia_error_handler (
			"quadratic_number_power_product_basis::"
			"absolute_Ln_approximations(matrix< bigint >, long)",
			"Precision overflow. (2)");
		return;
	}


	// Determine the approximation for each quadratic number basis[i].
	//
	for (i = 0; i < nofr; i++) {
		// find abs max row entry
		max_row_entry.assign_zero();

		for (j = 0; j < nofc; j++)
			if (max_row_entry.abs_compare(M.member(i, j)) < 0)
				max_row_entry.assign(M.member(i, j));

		// determine precision
		// t = kbn2 + b_value(max_row_entry);
		//
		if (check_overflow(t, kbn2, b_value(max_row_entry))) {
			lidia_error_handler (
				"quadratic_number_power_product_basis::"
				"absolute_Ln_approximations(matrix< bigint >, long)",
				"Precision overflow. (3)");
			return;
		}

		// Approximate Ln(basis[i]), if no precomputation exists or if
		// precomputation is not accurate enough.
		//
		if (t < -1)
			t = -1;

		if (!init[i] || (init[i] && (t > prec_ln_basis[i]))) {
			basis[i].get_absolute_Ln_approximation(ln_basis[i], t);
			prec_ln_basis[i] = t;
			init[i] = 1;
		}
	}
}



//
//  check_Ln_correctness
//

bool
quadratic_number_power_product_basis::check_Ln_correctness (lidia_size_t trials) const
{
	debug_handler ("quadratic_number_power_product_basis",
		       "check_Ln_correctness(lidia_size_t)const");

	lidia_size_t i, n;
	n = basis.get_size();

	for (i = 0; i < n; i++) {
		std::cout << " " << i << std::flush;
		if (!basis[i].check_Ln_correctness(trials))
			return false;
	}
	return true;
}



//
// conjugate
//

void
quadratic_number_power_product_basis::conjugate ()
{
	debug_handler ("quadratic_number_power_product_basis",
		       "conjugate()");

	// Compute conjugates.
	//
	lidia_size_t i, n;
	n = basis.get_size();

	for (i = 0; i < n; i++)
		basis[i].conjugate();

	// Correct Ln approximations also.
	//
	if (init.get_size() == n) {
		for (i = 0; i < n; i++)
			ln_basis[i].negate();
	}
}

//
//  input / output
//

std::istream &
operator >> (std::istream & in, quadratic_number_power_product_basis & b)
{
	debug_handler ("quadratic_number_power_product_basis",
		       "operator >> (std::istream, qnppb&)");

	b.basis.set_capacity(0);
	in >> b.basis;

	b.ln_basis.set_capacity(0);
	b.prec_ln_basis.set_capacity(0);
	b.init.set_capacity(0);

	return in;
}



std::ostream &
operator << (std::ostream & out, const quadratic_number_power_product_basis & b)
{
	debug_handler ("quadratic_number_power_product_basis",
		       "operator << (std::ostream, const qnbbp&)");
	out << b.basis;
	return out;
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
