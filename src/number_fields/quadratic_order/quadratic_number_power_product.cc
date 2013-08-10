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


//#define PPQN_ASSIGN_COMPACT_DEBUG
//#define QNPP_SPIG_DEBUG

#ifdef HAVE_CONFIG_H
# include	"config.h"
#endif
#include	"LiDIA/quadratic_number_power_product.h"
#include	"LiDIA/quadratic_number_standard.h"
#include	"LiDIA/quadratic_ideal.h"
#include	"LiDIA/quadratic_ideal_power_product.h"
#include	"LiDIA/quadratic_order.h"
#include	"LiDIA/prime_list.h"

#ifdef PPQN_ASSIGN_COMPACT_DEBUG
#include	<cassert>
#endif



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



xbigfloat quadratic_number_power_product::xbigfloat_dummy;
quadratic_number_power_product_basis quadratic_number_power_product::basis_dummy;

residue_class_list< quadratic_number_power_product_basis > quadratic_number_power_product::basis_list;

int quadratic_number_power_product::info = 0;
int quadratic_number_power_product::default_debug_verification_value = 1;


quadratic_number_power_product_basis*
quadratic_number_power_product::verify_preconditions (char * s) const
{
	debug_handler ("quadratic_number_power_product",
		       "verify_preconditions(char*)");

	quadratic_number_power_product_basis* b;
	lidia_size_t n;

	// verify preconditions
	//
	if (basis == NULL) {
		lidia_error_handler (s, "Basis not initialized (1).");
		return NULL;
	}

	b = static_cast<quadratic_number_power_product_basis*>(&(basis->get_mod()));
	n = b->get_size();

	if (b == NULL) {
		lidia_error_handler (s, "b == NULL.");
		return NULL;
	}

	if (n == 0) {
		lidia_error_handler (s, "Basis not initialized (2).");
		return NULL;
	}

	if (n != exp.get_size()) {
		lidia_error_handler (s,
				     "Different number of exponents "
				     "and base elements.");
		return NULL;
	}
	return b;
}



//
//  constructor
//
quadratic_number_power_product::quadratic_number_power_product ()
{
	debug_handler ("quadratic_number_power_product",
		       "quadratic_number_power_product()");
	basis = NULL;
	exp.set_mode(vector_flags::expand);
	do_debug_verification = quadratic_number_power_product::default_debug_verification_value;
}



//
//  destructor
//
quadratic_number_power_product::~quadratic_number_power_product()
{
	debug_handler ("quadratic_number_power_product",
		       "~quadratic_number_power_product");

	quadratic_number_power_product::basis_list.clear(basis);
}



//
//  Assignment
//

void
quadratic_number_power_product::set_basis (const quadratic_number_power_product_basis & b)
{
	debug_handler ("quadratic_number_power_product",
		       "set_basis(const quadratic_number_power_product_basis &)");

	quadratic_number_power_product::basis_list.clear(basis);
	basis = quadratic_number_power_product::basis_list.insert(b);
}



void
quadratic_number_power_product::set_basis (const quadratic_number_power_product & x)
{
	debug_handler ("quadratic_number_power_product",
		       "set_basis(const qnpp &)");

	x.verify_preconditions("quadratic_number_power_product::"
			       "set_basis(quadratic_number_power_product)");

	quadratic_number_power_product::basis_list.clear(basis);
	basis = quadratic_number_power_product::basis_list.set_to(x.basis);
}



void
quadratic_number_power_product::set_exponents (const base_vector< exp_type > & e)
{
	debug_handler ("quadratic_number_power_product",
		       "set_exponents(const base_vector< exp_type > &)");

	exp.assign(e);
	exp.set_mode(vector_flags::expand);
}



void
quadratic_number_power_product::reset ()
{
	debug_handler ("quadratic_number_power_product",
		       "reset()");

	exp.reset();
	quadratic_number_power_product::basis_list.clear(basis);
	basis = NULL;
}



quadratic_number_power_product &
quadratic_number_power_product::operator = (const quadratic_number_power_product & x)
{
	debug_handler ("quadratic_number_power_product",
		       "operator = (const qnpp&)");

	if (this != &x)
		assign(x);
	return *this;
}



quadratic_number_power_product &
quadratic_number_power_product::operator = (const quadratic_number_standard & x)
{
	debug_handler ("quadratic_number_power_product",
		       "operator = (const quadratic_number_standard&)");
	assign(x);
	return *this;
}



void
quadratic_number_power_product::assign (const quadratic_number_power_product & x)
{
	debug_handler ("quadratic_number_power_product",
		       "assign(const qnpp&)");

	if (this != &x) {
		exp.assign(x.exp);
		quadratic_number_power_product::basis_list.clear(basis);

		if (x.basis == NULL)
			basis = NULL;
		else
			basis = quadratic_number_power_product::basis_list.set_to(x.basis);
	}
}



void
quadratic_number_power_product::assign (const quadratic_number_standard & g)
{
	debug_handler ("quadratic_number_power_product",
		       "assign(const quadratic_number_standard&)");

	quadratic_number_power_product::basis_list.clear(basis);
	quadratic_number_power_product_basis b;
	b.set_basis(g);
	basis = quadratic_number_power_product::basis_list.insert(b);
	exp.set_size(1);
	exp[0] = 1;
}



void
quadratic_number_power_product::assign_one (const quadratic_order & O)
{
	debug_handler ("quadratic_number_power_product",
		       "assign_one(const quadratic_order &)");

	quadratic_number_standard alpha;
	alpha.assign_one(O);
	this->assign(alpha);
}



void
quadratic_number_power_product::assign (const quadratic_ideal & I,
					const xbigfloat & a,
					int s)
{
	debug_handler ("quadratic_number_power_product",
		       "assign(const quadratic_ideal&, const xbigfloat&, int)");

	// Let O = I.get_order(), I principal ideal
	// with alpha * O = I, |a - Ln alpha| < 2^{-4}, sign(alpha) = s.
	//
	// This function initializes the element by a compact
	// representation of alpha.
	//

	xbigfloat b, c, d, m, h;
	quadratic_ideal L, J, K;
	quadratic_number_standard gamma, delta;
	quadratic_number_power_product beta;
	long k;

	// Determine the accuracy, 2^{-k+2} < ln2.
	k = 3;

	// Reduce I
	K = I;
	K.reduce(gamma);
	gamma.get_absolute_Ln_approximation(c, k+1);

	// Find gamma / alpha.
	//
	// Determine minimum in O close to l
	subtract(m, c, a);
	J.assign(K.get_order());
	J.order_close(beta, b, m, k+1);

	// Found K again ?
	if (K != J) {
#ifdef PPQN_ASSIGN_COMPACT_DEBUG
		quadratic_ideal A1, A2;
		A1 = J; A1.rho();
		A2 = J; A2.inverse_rho();
		assert(K == A1 || K == A2);
#endif
		// If not, move one step to the left
		L = J;
		L.rho(delta);
		delta.get_absolute_Ln_approximation(d, k+1);

		// Check whether rho^{-2}(K) was found.
		//
		h.assign_one();
		shift_left(h, h, -k+1);

		if (b + d < m - h)
			J.inverse_rho(delta);
		else
			J.rho(delta);
		beta.multiply(beta, delta);
	}

	// Invert the base elements.
	//
	beta.invert_base_elements();
	this->multiply (beta, gamma);

	// Adjust sign
	if (s == -1)
		this->negate();
}



void
swap (quadratic_number_power_product & x,
      quadratic_number_power_product & y)
{
	debug_handler ("quadratic_number_power_product",
		       "swap(qnpp &, qnpp &)");

	swap (x.exp, y.exp);
	residue_class< quadratic_number_power_product_basis > *h = x.basis;
	x.basis = y.basis;
	y.basis = h;
}



//
//  Access
//
const quadratic_number_power_product_basis &
quadratic_number_power_product::get_basis () const
{
	debug_handler ("quadratic_number_power_product",
		       "get_basis() const");

	if (basis == NULL) {
		lidia_error_handler ("quadratic_number_power_product::get_basis",
				     "Basis not initialized.");

		return quadratic_number_power_product::basis_dummy;
	}
	else
		return basis->get_mod();
}



const base_vector< bigint > &
quadratic_number_power_product::get_exponents () const
{
	debug_handler ("quadratic_number_power_product",
		       "get_exponents() const");
	return exp;
}



bool
quadratic_number_power_product::is_initialized () const
{
	debug_handler ("quadratic_number_power_product",
		       "is_initialized() const");

	quadratic_number_power_product_basis* b;
	lidia_size_t n;

	if (basis == NULL)
		return false;

	b = static_cast<quadratic_number_power_product_basis*>(&(basis->get_mod()));
	n = b->get_size();

	if (b == NULL)
		return false;

	if (n == 0)
		return false;

	if (n != exp.get_size())
		return false;

	return true;
}



bool
quadratic_number_power_product::is_one_simple (const quadratic_number_power_product_basis* b) const
{
	debug_handler ("quadratic_number_power_product",
		       "is_one_simple(const qnppb *)");

	if (exp.get_size() != 1)
		return false;
	else if (((b->get_basis())[0]).is_one())
		return true;
	else
		return false;
}



const quadratic_order &
quadratic_number_power_product::get_order () const
{
	debug_handler ("quadratic_number_power_product",
		       "get_order() const");

	if (!is_initialized()) {
		lidia_error_handler("quadratic_number_power_product::get_order()",
				    "power product not initialized.");
	}
	return (basis->get_mod()).get_order();
}



//
// get_sign
//
// Returns 1, if the power product is greater or equal to zero.
// Otherwise, it returns -1.
//

int
quadratic_number_power_product::get_sign () const
{
	debug_handler ("quadratic_number_power_product",
		       "get_sign() const");

	quadratic_number_power_product_basis *b;
	lidia_size_t i, n;
	int sign;

	b = verify_preconditions("quadratic_number_power_product::"
				 "get_sign() const");
	if (exp[0].is_even())
		sign = 1;
	else
		sign = (b->get_basis())[0].get_sign();

	n = b->get_size();
	for (i = 1; i < n; i++) {
		if (exp[i].is_odd())
			sign *= (b->get_basis())[i].get_sign();
	}
	return sign;
}



//
//  negate
//

void
quadratic_number_power_product::negate ()
{
	debug_handler ("quadratic_number_power_product",
		       "negate()");

	quadratic_number_power_product_basis *b;

	b = verify_preconditions("quadratic_number_power_product::"
				 "negate()");
	// Multiply by -1.
	//
	quadratic_number_standard q;
	q.assign_one((b->get_basis())[0].get_order());
	q.negate();
	this->multiply (*this, q);
}



//
//  multiply
//

void
quadratic_number_power_product::multiply (const quadratic_number_standard & q,
					  const quadratic_number_power_product & x)
{
	debug_handler ("quadratic_number_power_product",
		       "multiply(const qns &, const qnpp &)");

	this->multiply(x, q);
}



void
quadratic_number_power_product::multiply (const quadratic_number_power_product & x,
					  const quadratic_number_standard & q)
{
	debug_handler ("quadratic_number_power_product",
		       "multiply(const qnpp &, const qns &)");

	quadratic_number_power_product_basis *bz, *bx;

	bx = x.verify_preconditions("quadratic_number_power_product::"
				    "multiply(pp, quadratic_number)(1)");
	if (q.is_one()) {
		this->assign(x);
		return;
	}

	if (x.is_one_simple(bx)) {
		this->assign(q);
		return;
	}

	if (basis != NULL) {
		// Object initialized. Reuse in case of just one reference.
		//
		if (basis->get_ref() == 1) {
			if (&x == this)
				bz = bx;
			else
				bz = this->verify_preconditions(
					"quadratic_number_power_product::"
					"multiply(pp, quadratic_number)(2)");
			bz->concat(*bx, q);
		}
		else // There are more references. Create a new basis, assign,
			// and clear old one. Remember, &z could be equal to &x.
		{
			bz = new quadratic_number_power_product_basis;
			bz->concat(*bx, q);
			quadratic_number_power_product::basis_list.clear(basis);
			basis = quadratic_number_power_product::basis_list.insert(*bz);
		}
	}
	else // z.basis == NULL
		// Create a new basis.
	{
		bz = new quadratic_number_power_product_basis;
		bz->concat(*bx, q);
		basis = quadratic_number_power_product::basis_list.insert(*bz);
	}

	if (this != &x) {
		exp.set_size(x.exp.get_size()+1);
		exp.assign(x.exp);
	}
	exp.insert_at(bigint(1), x.exp.get_size());
}



void
quadratic_number_power_product::multiply (const quadratic_number_power_product & x,
					  const quadratic_number_power_product & y)
{
	debug_handler ("quadratic_number_power_product",
		       "const qnpp &, const qnpp &)");

	quadratic_number_power_product_basis *bz, *bx, *by;

	bx = x.verify_preconditions("quadratic_number_power_product::"
				    "multiply(pp, pp) (1)");
	by = y.verify_preconditions("quadratic_number_power_product::"
				    "multiply(pp, pp) (2)");

	// If the basis of x and y are the same object
	// then just handle the exponents.
	//
	if (bx == by) {
		if (this != &x)
			quadratic_number_power_product::basis_list.clear(basis);

		basis = quadratic_number_power_product::basis_list.set_to(x.basis);
		add(exp, x.exp, y.exp);
	}

	// If the basis of x is a different object than the basis of y
	// then invert and concatenate.
	//
	else {
		if (basis != NULL) {
			if (basis->get_ref() == 1) {
				bz = verify_preconditions("quadratic_number_power_product::"
							  "multiply(pp, pp) (3)");
				bz->concat(*bx, *by);
			}
			else {
				bz = new quadratic_number_power_product_basis;
				bz->concat(*bx, *by);
				quadratic_number_power_product::basis_list.clear(basis);
				basis = quadratic_number_power_product::basis_list.insert(*bz);
			}
		}
		else {
			bz = new quadratic_number_power_product_basis;
			bz->concat(*bx, *by);
			basis = quadratic_number_power_product::basis_list.insert(*bz);
		}

		exp.concat(x.exp, y.exp);
	}
}



void
quadratic_number_power_product::invert (const quadratic_number_power_product & x)
{
	debug_handler ("quadratic_number_power_product",
		       "invert(const qnpp &)");

	this->assign(x);

	lidia_size_t i, n;

	n = exp.get_size();
	for (i = 0; i < n; i++)
		exp[i].negate();
}



void
quadratic_number_power_product::invert_base_elements ()
{
	debug_handler ("quadratic_number_power_product",
		       "invert_base_elements ()");

	lidia_size_t i, n;
	quadratic_number_power_product_basis *b;

	b = verify_preconditions("quadratic_number_power_product::"
				 "invert_base_elements()");

	n = exp.get_size();
	for (i = 0; i < n; i++) {
		(*b)[i].invert();
	}
}



void
quadratic_number_power_product::divide (const quadratic_number_standard & q,
					const quadratic_number_power_product & x)
{
	debug_handler ("quadratic_number_power_product",
		       "divide(const qns &, const qnpp &)");

	// Note that q = 1 is ignored in multiply.
	//
	lidia_size_t i, n;
	n = x.exp.get_size();

	this->multiply(x, q);

	for (i = 0; i < n; i++)
		exp[i].negate();

}



void
quadratic_number_power_product::divide (const quadratic_number_power_product & x,
					const quadratic_number_standard & q)
{
	debug_handler ("quadratic_number_power_product",
		       "divide(const qnpp &, const qns &)");

	// Note that q = 1 is ignored in multiplication.
	//
	lidia_size_t n;

	n = x.exp.get_size();
	this->multiply(x, q);

	if (n < exp.get_size())
		exp[n+1].negate();
}



void
quadratic_number_power_product::divide (const quadratic_number_power_product & x,
					const quadratic_number_power_product & y)
{
	debug_handler ("quadratic_number_power_product",
		       "divide(const qnpp &, const qnpp &)");

	quadratic_number_power_product_basis *bz, *bx, *by;

	bx = x.verify_preconditions("quadratic_number_power_product::"
				    "divide(pp, pp) (1)");
	by = y.verify_preconditions("quadratic_number_power_product::"
				    "divide(pp, pp) (2)");

	// If the basis of x and y are the same object
	// then just handle the exponents.
	//
	if (bx == by) {
		if (this != &x)
			quadratic_number_power_product::basis_list.clear(basis);

		basis = quadratic_number_power_product::basis_list.set_to(x.basis);
		subtract(exp, x.exp, y.exp);
	}

	// If the basis of x is a different object than the basis of y
	// then concatenate.
	//
	else {
		if (basis != NULL) {
			if (basis->get_ref() == 1) {
				bz = verify_preconditions("quadratic_number_power_product::"
							  "divide(pp, pp) (3)");
				bz->concat(*bx, *by);
			}
			else {
				bz = new quadratic_number_power_product_basis;
				bz->concat(*bx, *by);
				quadratic_number_power_product::basis_list.clear(basis);
				basis = quadratic_number_power_product::basis_list.insert(*bz);
			}
		}
		else {
			bz = new quadratic_number_power_product_basis;
			bz->concat(*bx, *by);
			basis = quadratic_number_power_product::basis_list.insert(*bz);
		}

		exp.concat(x.exp, y.exp);

		lidia_size_t i, n;
		n = exp.get_size();

		for (i = x.exp.get_size(); i < n; i++)
			exp[i].negate();
	}
}



void
quadratic_number_power_product::square (const quadratic_number_power_product & x)
{
	debug_handler ("quadratic_number_power_product",
		       "square(const qnpp&)");

	lidia_size_t i, n;

	if (this != &x)
		this->assign(x);

	n = exp.get_size();
	for (i = 0; i < n; i++)
		exp[i].multiply_by_2();
}



void
quadratic_number_power_product::power (const quadratic_number_power_product & x,
				       const bigint & e)
{
	debug_handler ("quadratic_number_power_product",
		       "power(qnpp &, const bigint &)");

	if (e.is_zero())
		this->assign_one(x.get_order());
	else if (e.is_one())
		this->assign(x);
	else {
		if (this != &x) {
			quadratic_number_power_product::basis_list.clear(basis);
			basis = quadratic_number_power_product::basis_list.set_to(x.basis);
		}

		lidia_size_t i, n;
		n = x.exp.get_size();
		for (i = 0; i < n; i++)
			LiDIA::multiply(exp[i], x.exp[i], e);

		exp.set_size(n);
	}
}



quadratic_number_power_product &
quadratic_number_power_product::operator *= (const quadratic_number_standard & x)
{
	debug_handler ("quadratic_number_power_product",
		       "operator *= (const qns &)");

	this->multiply(*this, x);
	return *this;
}



//
// conjugate
//

void
quadratic_number_power_product::conjugate ()
{
	debug_handler ("quadratic_number_power_product",
		       "conjugate()");

	quadratic_number_power_product_basis *b;

	b = this->verify_preconditions("quadratic_number_power_product::"
				       "conjugate()");

	if (basis->get_ref() > 1) {
		// There are more than one reference.
		// Copy basis and clear old one.
		//
		quadratic_number_power_product_basis *bz;

		bz = new quadratic_number_power_product_basis;
		bz->assign(*b);
		quadratic_number_power_product::basis_list.clear(basis);
		basis = quadratic_number_power_product::basis_list.insert(*bz);
		b = bz;
	}

	// Now, compute conjugates of base elements.
	//
	b->conjugate();
}



//
// norm
//

void
quadratic_number_power_product::norm_modulo (bigint & num,
					     bigint & den,
					     const bigint & m) const
{
	debug_handler ("quadratic_number_power_product",
		       "norm_modulo(bigint&, bigint&, const bigint&)");

	lidia_size_t i, n;
	quadratic_number_power_product_basis* b;

	bigmod Num, Den, NumTmp, DenTmp;
	bigrational BaseNorm;
	bigint e;

	if (m.is_zero()) {
		lidia_error_handler("quadratic_number_power_product::"
				    "norm_modulo(num, den, m)",
				    "Zero modulus m.");
		return;
	}


	// verify preconditions
	//
	b = verify_preconditions("quadratic_number_power_product::"
				 "norm_modulo(bigint&, bigint&, const bigint&)");
	if (b != NULL) {
		// Set bigmod modulus
		//
		bigint old_modulus = bigmod::modulus();
		bigmod::set_modulus(m);

		// Compute norm modulo m
		//
		n = b->get_size();
		Num.assign_one();
		Den.assign_one();

		for (i = 0; i < n; i++) {
			((b->get_basis())[i]).norm(BaseNorm);

			if (exp[i].is_positive()) {
				LiDIA::power(NumTmp, bigmod(BaseNorm.numerator()), exp[i]);
				LiDIA::power(DenTmp, bigmod(BaseNorm.denominator()), exp[i]);
			}
			else {
				e.assign(exp[i]);
				e.negate();
				LiDIA::power(DenTmp, bigmod(BaseNorm.numerator()), e);
				LiDIA::power(NumTmp, bigmod(BaseNorm.denominator()), e);
			}
			Num *= NumTmp;
			Den *= DenTmp;
		}

		num.assign(Num.mantissa());
		den.assign(Den.mantissa());

		// Restore bigmod modulus
		//
		if (!old_modulus.is_zero())
			bigmod::set_modulus(old_modulus);
	}
}



//
// Evaluation
//

quadratic_number_standard
quadratic_number_power_product::evaluate () const
{
	debug_handler ("quadratic_number_power_product",
		       "evaluate() const");

	lidia_size_t i, n;
	quadratic_number_power_product_basis* b;

	quadratic_number_standard res, h;

	// verify preconditions
	b = verify_preconditions("quadratic_number_power_product::evaluate()");

	if (b != NULL) {
		// evaluate
		n = b->get_size();

		LiDIA::power(res, (b->get_basis())[0], exp[0]);
		for (i = 1; i < n; i++) {
			LiDIA::power(h, (b->get_basis())[i], exp[i]);
			LiDIA::multiply(res, res, h);
		}
	}
	return res;
}



#if 0
quadratic_number_standard
quadratic_number_power_product::evaluate_modulo (const bigint & m) const
{
	debug_handler ("quadratic_number_power_product",
		       "evaluate_modulo(const bigint &)");

	lidia_size_t i, n;
	quadratic_number_power_product_basis* b;

	quadratic_number_standard res, h;

	// verify preconditions
	b = verify_preconditions("quadratic_number_power_product::evaluate()");

	if (b != NULL) {
		// evaluate
		n = b->get_size();

		LiDIA::power(res, (b->get_basis())[0], exp[0]);
		for (i = 1; i < n; i++) {
			LiDIA::power(h, (b->get_basis())[i], exp[i]);
			LiDIA::multiply(res, res, h);
		}
	}
	return res;
}



#endif




//
//  Comparison
//

bool
operator == (const quadratic_number_power_product & x,
	     const quadratic_number_power_product & y)
{
	debug_handler ("quadratic_number_power_product",
		       "operator == (const qnpp &, const qnpp &)");

	quadratic_number_standard alpha, beta;

	alpha = x.evaluate();
	beta = y.evaluate();
	return (alpha == beta);
}



bool
operator == (const quadratic_number_power_product & x,
	     const quadratic_number_standard & q)
{
	debug_handler ("quadratic_number_power_product",
		       "operator == (const qnpp &, const qns &)");

	quadratic_number_standard alpha;

	alpha = x.evaluate();
	return (q == alpha);

}



bool
operator == (const quadratic_number_standard & q,
	     const quadratic_number_power_product & x)
{
	debug_handler ("quadratic_number_power_product",
		       "operator == (const qns &, const qnpp &)");

	return (x == q);
}



bool
operator != (const quadratic_number_power_product & x,
	     const quadratic_number_power_product & y)
{
	debug_handler ("quadratic_number_power_product",
		       "operator != (const qnpp &, const qnpp &)");

	if (x == y)
		return false;
	else
		return true;
}



bool
operator != (const quadratic_number_power_product & x,
	     const quadratic_number_standard & q)
{
	debug_handler ("quadratic_number_power_product",
		       "operator != (const qnpp &, const qns &)");

	if (x == q)
		return false;
	else
		return true;
}



bool
operator != (const quadratic_number_standard & q,
	     const quadratic_number_power_product & x)
{
	debug_handler ("quadratic_number_power_product",
		       "operator != (const qns &, const qnpp &)");

	if (x == q)
		return false;
	else
		return true;
}



//
//  absolute k-approximation to Ln of power product
//

xbigfloat
quadratic_number_power_product::get_absolute_Ln_approximation (long k) const
{
	debug_handler ("quadratic_number_power_product",
		       "get_absolute_Ln_approximation(long)");

	xbigfloat l;
	this->get_absolute_Ln_approximation(l, k);
	return l;
}



xbigfloat
quadratic_number_power_product::absolute_Ln_approximation (long k) const
{
	debug_handler ("quadratic_number_power_product",
		       "absolute_Ln_approximation(long)");

	xbigfloat l;
	this->get_absolute_Ln_approximation(l, k);
	return l;
}



void
quadratic_number_power_product::get_absolute_Ln_approximation (xbigfloat & l, long k) const
{
	debug_handler ("quadratic_number_power_product",
		       "absolute_Ln_approximation(xbigfloat&, long)");

	quadratic_number_power_product_basis* b;
	lidia_size_t i, n;
	long t, kbn2;
	xbigfloat h;

	b = verify_preconditions("quadratic_number_power_product::"
				 "get_absolute_Ln_approximation"
				 "(xbigfloat, long)");
	if (b == NULL)
		return;

	n = b->get_size();

	//
	// Determine the precision kbn2
	//
	// kbn2 = k+b_value(n)+2
	//
	//
	if (check_overflow(kbn2, k, 2)) {
		lidia_error_handler ("quadratic_number_power_product::"
				     "get_absolute_Ln_approximation",
				     "precision overflow.");
		return;
	}

	if (check_overflow(kbn2, kbn2, b_value(n))) {
		lidia_error_handler ("quadratic_number_power_product::"
				     "get_absolute_Ln_approximation",
				     "precision overflow.");
		return;
	}

	//
	// Determine an absolute k approximation to Ln of power product
	//

	l.assign_zero();

	for (i = 0; i < n; i++) {
		if (!(exp[i].is_zero())) {
			// absolute precision t for Ln(basis[i])
			// t = kbn2 + b_value(exp[i])
			//
			if (check_overflow(t, kbn2, b_value(exp[i]))) {
				lidia_error_handler ("quadratic_number_power_product::"
						     "get_absolute_Ln_approximation",
						     "precision overflow.");
				return;
			}
			LiDIA::multiply (h, b->get_absolute_Ln_approximation(t, i), exp[i]);

			// truncate and add h to the sum
			// t = kbn2 + b_value(h)

			if (check_overflow(t, kbn2, b_value(h))) {
				lidia_error_handler ("quadratic_number_power_product::"
						     "get_absolute_Ln_approximation",
						     "precision overflow.");
				return;
			}

			truncate (h, h, t);
			add (l, l, h);
		}
	}

	// truncate Ln approximation
	// t = k + 1 + b_value(l)

	if (check_overflow(t, k+1, b_value(l))) {
		lidia_error_handler ("quadratic_number_power_product::"
				     "get_absolute_Ln_approximation",
				     "precision overflow.");
		return;
	}

	truncate(l, l, t);
}



void
quadratic_number_power_product::absolute_Ln_approximation (xbigfloat & l, long k) const
{
	debug_handler ("quadratic_number_power_product",
		       "absolute_Ln_approximation(xbigfloat&, long)");
	this->get_absolute_Ln_approximation(l, k);
}



//
//  relative k-approximation to Ln of power product
//

xbigfloat
quadratic_number_power_product::get_relative_Ln_approximation (long k, long m) const
{
	debug_handler ("quadratic_number_power_product",
		       "get_relative_Ln_approximation(long, long)");
	xbigfloat l;
	this->get_relative_Ln_approximation(l, k, m);
	return l;
}



xbigfloat
quadratic_number_power_product::relative_Ln_approximation (long k, long m) const
{
	debug_handler ("quadratic_number_power_product",
		       "relative_Ln_approximation(long, long)");
	xbigfloat l;
	this->get_relative_Ln_approximation(l, k, m);
	return l;
}



void
quadratic_number_power_product::get_relative_Ln_approximation (xbigfloat & l, long k, long m) const
{
	debug_handler ("quadratic_number_power_product",
		       "get_relative_Ln_approximation(xbigfloat&, long, long)");
	long t;

	// t = k+1-m
	if (check_overflow(t, k, 1)) {
		lidia_error_handler ("quadratic_number_power_product::relative_Ln_approximation(xbig, long, long)",
				     "precision overflow.");
		return;
	}
	if (check_overflow(t, t, -m)) {
		lidia_error_handler ("quadratic_number_power_product::relative_Ln_approximation(xbig, long, long)",
				     "precision overflow.");
		return;
	}

	this->get_absolute_Ln_approximation (l, t);

	// t = k+1+b_value(l)
	k++;
	if (check_overflow(t, k, b_value(l))) {
		lidia_error_handler ("quadratic_number_power_product::relative_Ln_approximation(xbig, long, long)",
				     "precision overflow.");
		return;
	}

	l.truncate (t);
}



//
//  compact_representation
//

void
quadratic_number_power_product::compact_representation (const quadratic_ideal & I)
{
	debug_handler ("quadratic_number_power_product",
		       "compact_representation(const quadratic_ideal &)");

	// Let a = *this, O = this->get_order().
	// We assume that a * O = I.
	//

	// Check for same order.
	//
	if (this->get_order() != I.get_order()) {
		lidia_error_handler("quadratic_number_power_product::"
				    "compact_representation(const quadratic_ideal&)",
				    "Order of element and order of ideal are not the same!");
		return;
	}

	// Compute absolute 4 approx. l to Ln a and sign of a.
	//
	xbigfloat l;
	this->get_absolute_Ln_approximation(l, 4);
	this->assign(I, l, this->get_sign());
}



//
//  Short principal ideal generator.
//

void
quadratic_number_power_product::get_short_principal_ideal_generator (xbigfloat & l,
								     bigint    & z,
								     const quadratic_number_power_product & rho,
								     long b,
								     long k) const
{
	debug_handler ("quadratic_number_power_product",
		       "short_principal_ideal_generator("
		       "const xbigfloat&, const bigint&, "
		       "const quadratic_number_power_product&, long, long)const");

	xbigfloat w, r, h1, h2;

	// approximate Ln(omega)
	//  w = this->absolute_Ln_approximation(-b+4+max(b-2+k,4));
	//
	if (b+k-2 >= 4)
		w = this->get_absolute_Ln_approximation(k+2);
	else
		w = this->get_absolute_Ln_approximation(-b+8);

	// Check for |w| < |r|/2
	//
	// if (|w| - 2^(b-5) < 2^(b-4))
	//  return truncate(w,w.exponent()+k+1);
	//
	h1.assign(2);
	shift_left(h1, h1, b-5);
	h2.assign(2);
	shift_left(h2, h2, b-4);
	h2 += h1;

	h1.assign(w);
	h1.absolute_value();

	if (h1 < h2) {
		w.truncate(w.exponent()+5);
		l = w;
		z = 0;
		return;
	}

	// approximate regulator Ln(rho)
	//
	// r = rho.absolute_Ln_approximation(-2b+b(w)+8+max(8,b+3+k));
	//
	if (b+k+3 >= 8)
		r = rho.get_absolute_Ln_approximation(-b + w.b_value() + 11 +k);
	else
		r = rho.get_absolute_Ln_approximation(-2*b + w.b_value() + 16);

	r.absolute_value();

	// | floor(Ln(omega)/|Ln(rho)|) - z | <= 1
	//
	LiDIA::divide (z, w, r);

#ifdef QNPP_SPIG_DEBUG
	std::cout << "Found z = " << z << std::endl;
#endif

	// Try to find -1 <= i <= 2 such that
	//
	//   -|Ln(rho)|/2 < Ln(omega) - (z+i) |Ln(rho)| <= |Ln(rho)|/2
	//
	// Set l = w - (z+i) * r.
	//
	long s, t;
	xbigfloat TwoToS, TwoToT, TwoToTM1, TwoToTP1, r2;

	// s = min(-k-1, b-6);
	//
	if (-k-1 <= b-6)
		s = -k-1;
	else
		s = b-6;

	// t = 2b-b(w)-8-max(8,b+k+3);
	if (b+k+3 >= 8)
		t = -(-b + w.b_value() + 11 +k);
	else
		t = -(-2*b + w.b_value() + 16);

	l = w - z * r;

	// TwoToS = 2^s
	TwoToS.assign(2);
	shift_left(TwoToS, TwoToS, s);

	// TwoToT = 2^t
	TwoToT.assign(2);
	shift_left(TwoToT, TwoToT, t);

	// TwoToTM1 = 2^(t-1)
	TwoToTM1.assign(2);
	shift_left(TwoToTM1, TwoToTM1, t-1);

	// TwoToTP1 = 2^(t+1)
	TwoToTP1.assign(2);
	shift_left(TwoToTP1, TwoToTP1, t+1);

	// r2 = r/2
	r2.assign(r);
	shift_right(r2, r2, 1);

	//if (l + 2^s <= -r/2 - 2^(t-1))
	// i = -1, l += r;
	//
	if (l + TwoToS <= -r2 - TwoToTM1) {
		l += r;
		z -= 1;
#ifdef QNPP_SPIG_DEBUG
		std::cout << "l <= -r/2, z = " << z << std::endl;
#endif
	}

	// else if (-r/2 + 2^(t-1) <= l - 2^s &&
	//    l + 2^s <= r/2 - 2^(t-1))
	//    i = 0
	//
	else if (-r2 + TwoToTM1 <= l - TwoToS &&
		 l + TwoToS <= r2 - TwoToTM1) {
#ifdef QNPP_SPIG_DEBUG
		std::cout << "-r/2 < l <= r/2, z = " << z << std::endl;
#endif
	}
	//else if (r/2 + 2^(t-1) <= l - 2^s &&
	//	    l + 2^s <= 3/2 * r - 2^(t+1))
	//    i = 1, l = l - r;
	//
	else if (r2 + TwoToTM1 <= l - TwoToS &&
		 l + TwoToS <= 3 * r2 - TwoToTP1) {
		l -= r;
		z += 1;
#ifdef QNPP_SPIG_DEBUG
		std::cout << "r/2 < l < 3/2 r, z = " << z << std::endl;
#endif
	}
	//else if (3/2 * r + 2^(t+1) <= l - 2^s)
	//   i = 2, l = l - 2*r;
	//
	else if (3 * r2 + TwoToTP1 <= l - TwoToS) {
		l -= 2*r;
		z += 2;
#ifdef QNPP_SPIG_DEBUG
		std::cout << "3/2 r < l, z = " << z << std::endl;
#endif
	}

	// Ln(alpha) is too close to |Ln(rho)|/2. Find -1 <= i <= 1
	// such that
	//
	//  0 <= Ln(omega) - (z+i) |Ln(rho)| < |Ln(rho)|
	//
	// Set l = w - (z+i) * r.
	//

	// else if (r + 2^t <= l - 2^s)
	//  i = 1, l = l - r;
	//
	else if (r + TwoToT <= l - TwoToS) {
		l -= r;
		z += 1;
#ifdef QNPP_SPIG_DEBUG
		std::cout << "r < l, z = " << z << std::endl;
#endif
	}
	// else if (l + 2^s <= 0)
	//   i = -1, l = l + r;
	//
	else if (l + TwoToS <= 0) {
		l += r;
		z -= 1;
#ifdef QNPP_SPIG_DEBUG
		std::cout << "l < 0, z = " << z << std::endl;
#endif
	}

	l.truncate(l.exponent() + k+1);
}



//
//   generating unit
//

void
quadratic_number_power_product::generating_unit (const char * units_input_file)
{
	debug_handler ("quadratic_number_power_product",
		       "generating_unit(const char*)");

	std::ifstream in(units_input_file);

	if (!in)
		lidia_error_handler("quadratic_number_power_product::generating_unit(char*)",
				    "Can't open input file");
	else {
		matrix< bigint > M;
		base_vector< quadratic_number_standard > q;
		long m;

		in >> M >> q >> m;

		std::cout << "m = " << m << std::endl;

		this->generating_unit(M, q, m);
	}
}



void
quadratic_number_power_product::generating_unit (const matrix< bigint > & M,
						 const base_vector< quadratic_number_standard > & q,
						 long m,
						 int strategy)
{
	debug_handler ("quadratic_number_power_product",
		       "generating_unit(const matrix< bigint > &, "
		       "const base_vector< quadratic_number_standard > &, "
		       "long, int)");

	//std::ofstream out ("dd56.units");
	//out << M;
	//out << q;
	//out << m;

	lidia_size_t i, c, r;
	quadratic_number_power_product_basis b;
	base_vector< quadratic_number_power_product > p;
	base_vector< bigint > column;

	// initialize
	//
	c = M.get_no_of_columns();
	r = M.get_no_of_rows();
	b.set_basis(q);

	if (strategy == 0) {
		// Precompute approximations to logarithms of quadratic numbers.
		//
		// Note: Based on an upper bound on the absolute values
		//       of the Ln of the quadratic numbers, a heuristic upper
		//       bound on the largest required absolute precision is
		//       computed and used to compute approximations to the Ln's
		//       of the quadratic numbers.
		//
		//       Nevertheless, during the computation all sufficient
		//       accuracies are determined and if necessary, the
		//       approximations to the Ln's are recomputed to a higher
		//       precision. Hence, the algorithm really yields the
		//       generating unit.
		//
		xbigfloat L, U, Umax;

		// Determine upper bound Umax on max(i) |Ln q[i]|
		//
		q[0].get_Ln_estimate(L, U);

		for (i = 1; i < r; i++) {
			q[i].get_Ln_estimate(L, U);
			if (U > Umax)
				U.assign(Umax);
		}

		b.absolute_Ln_approximations (M, quadratic_number_power_product::estimate_pp_accuracy (M, m, Umax));
	}

	// Transform matrix and numbers into a vector of numbers.
	//
	p.set_capacity(c);

	M.get_column_vector(column, 0);
	p[0].set_exponents(column);
	p[0].set_basis(b);

	for (i = 1; i < c; i++) {
		M.get_column_vector(column, i);
		p[i].set_exponents(column);
		p[i].set_basis(p[0]);
	}

	// Compute generating unit
	//
	if (strategy == 0)
		this->generating_unit(p, m, 1);
	else
		this->generating_unit(p, m, 0);
}



void
quadratic_number_power_product::generating_unit (base_vector< quadratic_number_power_product > & q,
						 long m,
						 int strategy)
{
	debug_handler ("quadratic_number_power_product",
		       "generating_unit(base_vector< qu > &, long, int)");

#ifdef DEBUG
	if (info < 3) info = 3;
#endif

	lidia_size_t i, j, n;
	lidia_size_t step_size, neighbour, iterations;
	quadratic_ideal O;
	bigint x, y, M1, M2;

	// check preconditions
	//
	n = q.get_size();

	if (n == 0) {
		lidia_error_handler("quadratic_number_power_product::"
				    "generating_unit",
				    "No units given.");
		return;
	}
	else if (info > 3) {
		std::cout << "quadratic_number_power_product::generating_unit::";
		std::cout << std::endl;
		std::cout << q.get_size() << " units given." << std::endl;
	}

	if (do_debug_verification) {
		xbigfloat l;
		quadratic_number_standard u_qns;
		quadratic_number_power_product a, b;

#if 0
		  std::cout << ": Ln test " << std::flush;
		  if (q[0].get_basis().check_Ln_correctness())
			  std::cout << "OK" << std::endl;
		  else
			  lidia_error_handler("qu::generating_unit()",
					      "FAILED");
#endif

		for (i = 0; i < n; i++) {
#if 0
			  std::cout << "q[" << i << "]" << std::flush;

			  std::cout << ": close test " << std::flush;
			  if (q[i].could_be_quasi_unit_close_test())
				  std::cout << "OK" << std::flush;
			  else
				  std::cout << "FAILED" << std::flush;
#endif
#if 0
			  std::cout << ": norm test " << std::flush;
			  if (q[i].could_be_unit_norm_test())
				  std::cout << "OK" << std::flush;
			  else
				  std::cout << "FAILED" << std::flush;
#endif

#if 0
			  std::cout << ": refinement test " << std::flush;
			  b.assign(q[i]);
			  b.conjugate();
			  a.divide(q[i], b);

			  if (a.could_be_unit_refinement_test())
				  std::cout << "OK" << std::flush;
			  else
				  std::cout << "FAILED" << std::flush;
#endif
#if 0
			  std::cout << std::endl;
			  l = q[i].get_absolute_Ln_approximation(3);
			  std::cout << "abs. 3-approx. to Ln is l = " << l << std::endl;
#endif
		}
	}

	// check strategy
	//
	if (strategy == 1) {
		if (info > 3) {
			std::cout << "quadratic_number_power_product::generating_unit";
			std::cout << "(base_vector< qu >, long, int)::" << std::endl;
			std::cout << "strategy 1, no compact representations." << std::endl;
		}
	}
	else {
		if (info > 3) {
			std::cout << "quadratic_number_power_product::generating_unit";
			std::cout << "(base_vector< qu >, long, int)::" << std::endl;
			std::cout << "default strategy, compact representations." << std::endl;
		}

		// Compute compact representations
		//
		timer t;

		O.assign_one(q[0].get_order());
		n = q.get_size();
		t.start_timer();

		for (i = 0; i < n; i++) {
			if (info > 3) {
				std::cout << i << "-th power product, size before reduction = ";
				std::cout << q[i].get_basis().get_size() << std::endl;
			}

			q[i].compact_representation_of_unit();

			if (info > 3) {
				std::cout << i << "-th power product, size after reduction = ";
				std::cout << q[i].get_basis().get_size() << std::endl;
			}
		}

		t.stop_timer();
		if (info > 3)
			std::cout << "Reduction time is "; t.print(); std::cout << std::endl;
	}


	// remove +-1
	//
	quadratic_number_power_product::remove_rationals_from_quasi_units (q, m);
	if (info > 3) {
		std::cout << "quadratic_number_power_product::generating_unit::";
		std::cout << std::endl;
		std::cout << q.get_size() << " units remain after removing +-1." << std::endl;
	}

	if (q.get_size() == 0) {
		this->assign_one(O.get_order());
		return;
	}
	else if (q.get_size() == 1) {
		this->assign(q[0]);
		return;
	}


	// binary tree method with application of rgcd
	//
	n = q.get_size();
	iterations = b_value(static_cast<long>(n));
	step_size = 2;
	neighbour = 1;

	for (i = 0; i < iterations; i++) {
		if (info > 3) {
			std::cout << "quadratic_number_power_product::generating_unit::" << std::endl;
			std::cout << "STEP " << i << std::endl;
		}

		for (j = 0; j+neighbour < n; j += step_size) {
			if (info > 3) {
				std::cout << "combining ";
				std::cout << "(" << j << ", " << j+neighbour << ")" << std::endl;
			}

			// representation of real gcd
			// resulting unit is q[j]^x * q[j+neighbour]^y
			//
			q[j].generating_unit(q[j], q[j+neighbour], m);
		}

		step_size <<= 1;
		neighbour <<= 1;
	}

	this->assign(q[0]);
}



void
quadratic_number_power_product::generating_unit (const quadratic_number_power_product & q1,
						 const quadratic_number_power_product & q2,
						 long m)
{
	debug_handler ("quadratic_number_power_product",
		       "const qu&, const qu&, long)");
	bigint M1, M2, x, y;
	this->generating_unit(x, y, M1, M2, q1, q2, m);
}



void
quadratic_number_power_product::generating_unit (bigint & x,
						 bigint & y,
						 bigint & M1,
						 bigint & M2,
						 const quadratic_number_power_product & q1,
						 const quadratic_number_power_product & q2,
						 long m)
{
	debug_handler ("quadratic_number_power_product",
		       "x, y, M1, M2, q1, q2, m");

	if (&q1 == &q2) {
		x.assign_one();
		y.assign_zero();
		M1.assign_one();
		M2.assign_one();
		this->assign(q1);
	}
	else {
		// Compute x,y such that
		//
		// Ln q1 = M1 * M * R,
		// Ln q2 = M2 * M * R, gcd (M1,M2) = 1,
		//
		// x * Ln q1 + y * Ln q2 = M * R,
		//
		// where R is the regulator, and q1, q2 are in
		// Q * O^{*}.
		//
		quadratic_number_power_product::rgcd(x, y, M1, M2, q1, q2, m);

		// Compute u = q1^{x} * q2^{y} in Q * O^{*}.
		//
		// Then Ln u = M * R.
		//

		// no verification for debugging purposes.
		// optimize.
		//
		if (!do_debug_verification) {
			if (x.is_zero())
				this->power(q2, y);
			else if (y.is_zero())
				this->power(q1, x);
			else {
				quadratic_number_power_product p1, p2;
				p1.power(q1, x);
				p2.power(q2, y);
				this->multiply(p1, p2);
			}
		}

		// verification for debugging.
		// Compute result into a new object
		// and test the results.
		//
		else {
			quadratic_number_power_product p1, p2;
			quadratic_number_power_product u, w1, w2;

			if (x.is_zero())
				u.power(q2, y);
			else if (y.is_zero())
				u.power(q1, x);
			else {
				p1.power(q1, x);
				p2.power(q2, y);
				u.multiply(p1, p2);
			}

			// We have Ln u^M1 = Ln q1
			// and     Ln u^M2 = Ln q2
			//
			// Hence, u^M1 / q1, u^M2 / q2
			// must be rational numbers.
			//
			w1.power(u, M1);
			w2.power(u, M2);
			p1.divide(w1, q1);
			p2.divide(w2, q2);

			if (!p1.is_rational(m)) {
				// std::cout << "q1 = " << q1 << std::endl;
				// std::cout << "q2 = " << q2 << std::endl;
				std::cout << "m = " << m << std::endl;
				std::cout << "x = " << x << std::endl;
				std::cout << "y = " << y << std::endl;
				std::cout << "M1 = " << M1 << std::endl;
				std::cout << "M2 = " << M2 << std::endl;
				lidia_error_handler("quadratic_number_power_product"
						    "::generating_unit"
						    "(quadratic_number_power_product, "
						    "quadratic_number_power_product, long)",
						    "wrong generator (1)");
			}

			if (!p2.is_rational(m)) {
				//std::cout << "q1 = " << q1 << std::endl;
				//std::cout << "q2 = " << q2 << std::endl;
				std::cout << "m = " << m << std::endl;
				std::cout << "x = " << x << std::endl;
				std::cout << "y = " << y << std::endl;
				std::cout << "M1 = " << M1 << std::endl;
				std::cout << "M2 = " << M2 << std::endl;
				lidia_error_handler("quadratic_number_power_product"
						    "::generating_unit"
						    "(quadratic_number_power_product, "
						    "quadratic_number_power_product, long)",
						    "wrong generator (2)");
			}

			this->assign(u);
		}
	}
}



//
//  is_rational
//
//  Returns true, if this in Q and false otherwise.
//  The function assumes that this in Q * O^{*}, i.e.,
//  is a rational number times a unit.
//
//  The result is derived by computing an absolute -m+1 approximation l
//  to Ln(this) and verifying |l| < 2^{m-1}. 2^{m} must be a lower bound
//  to the regulator.
//
//  The caller must provide m.
//

bool
quadratic_number_power_product::is_rational (long m) const
{
	debug_handler ("quadratic_number_power_product",
		       "is_rational(long) const");

	xbigfloat l;
	long mmp1; // minus m plus 1
	long mm1; // m minus 1

	if (check_overflow(mmp1, -m, 1)) {
		lidia_error_handler ("quadratic_number_power_product::is_rational::",
				     "precision overflow.");
		return false;
	}

	mm1 = -mmp1;

	// absolute -m+1 approx. to Ln
	this->get_absolute_Ln_approximation(l, mmp1);

	// this == +-1 iff |l| < 2^{m-1} iff (l == 0 || b(l) <= m-1)
	if (l.is_zero() || b_value(l) <= mm1)
		return true;
	else
		return false;
}



//
//  remove_rationals_from_quasi_units
//

void
quadratic_number_power_product::remove_rationals_from_quasi_units (
	base_vector< quadratic_number_power_product > & q,
	long m)
{
	debug_handler ("quadratic_number_power_product",
		       "remove_rationals_from_quasi_units"
		       "(base_vector< qnpp > &, long)");

	lidia_size_t c, i;
	lidia_size_t nzi; // non-zero index

	c = q.get_size();
	nzi = 0;

	for (i = 0; i < c; i++) {
		if (q[i].is_rational(m)) {
			if (quadratic_number_power_product::info > 3) {
				std::cout << "quadratic_number_power_product::";
				std::cout << "remove_rationals_from_quasi_units::" << std::endl;
				std::cout << "unit " << i << " is +-1." << std::endl;
			}
		}
		else {
			// store element i at index nzi
			if (i > nzi)
				swap(q[nzi], q[i]);

			nzi++;
		}
	}

	q.set_size(nzi);
}



//
// ::rgcd
//
// Determines a pair (x, y) with
//
// x Ln q1 + y Ln q2 = rgcd(Ln q1, Ln q2) = r.
//
// Ln q1 = M1 r, Ln q2 = M2 r, 1 = x M1 + y M2.
//
// Regulator > 2^m.
//
// Condition: q1, q2 not rational.
//
//

void
quadratic_number_power_product::rgcd (bigint & x,
				      bigint & y,
				      bigint & M1,
				      bigint & M2,
				      const quadratic_number_power_product & q1,
				      const quadratic_number_power_product & q2,
				      long m)
{
	debug_handler ("quadratic_number_power_product",
		       "rgcd(x, y, M1, M2, q1, q2, m)");

	xbigfloat s, t, h;
	bigint    S, num, den, gcdM1M2;
	bigint    big_k, big_ks, big_kt;
	bigint    bs, bt;
	long      k, ks, kt;
	bool      s_is_negative, t_is_negative;

	// s = relative 1-approximation to Ln(q1)
	// t = relative 1-approximation to Ln(q2)
	// |b(s) - b(Ln(q1))| <= 1
	// |b(t) - b(Ln(q2))| <= 1
	//
	big_k = -m;
	++big_k;

	if (big_k.longify(k)) {
		lidia_error_handler ("quadratic_number_power_product::rgcd",
				     "precision overflow.");
		return;
	}
	if (quadratic_number_power_product::info > 3) {
		std::cout << "quadratic_number_power_product::rgcd::" << std::endl;
		std::cout << "COMPUTING absolute " << k << " approximations";
		std::cout << std::endl;
	}

	q1.get_absolute_Ln_approximation (s, k);
	q2.get_absolute_Ln_approximation (t, k);
	bs = b_value(s);
	bt = b_value(t);

	if (quadratic_number_power_product::info > 3) {
		std::cout << "quadratic_number_power_product::rgcd::" << std::endl;
		std::cout << "FOUND s = " << s << std::endl;
		std::cout << "FOUND t = " << t << std::endl;
		std::cout << "FOUND b_value(s) = " << bs << std::endl;
		std::cout << "FOUND b_value(t) = " << bt << std::endl;
	}

	// upper bound S on |M2|
	// S = ceil( 2*|t| / 2^m );
	//
	h.assign(t);
	h.absolute_value();
	if (k >= 0)
		shift_left (h, h, k);
	else
		shift_right (h, h, -k);

	ceil(S, h);

	if (quadratic_number_power_product::info > 3) {
		std::cout << "quadratic_number_power_product::rgcd::" << std::endl;
		std::cout << "FOUND upper bound h on |M2|, h = " << h << std::endl;
		std::cout << "S = ceil(h), S = " << S << std::endl;
	}

	// precisions to guarantee |R1/R2 - s/t| < 1/(2S^2)
	// ks = 2*b_value(S)
	// kt = b_value(t)-b_value(s)
	// k = max {ks, kt} +1
	//
	big_ks = 2 * b_value(S);
	big_kt = bt-bs;

	if (big_ks > big_kt)
		big_k = big_ks;
	else
		big_k = big_kt;

	++big_k;

	// s = relative ks approx. to Ln(q1)
	// t = relative kt approx. to Ln(q2)
	//
	// ks = b_value(s)-b_value(t)+4+k;
	// kt = b_value(s)-b_value(t)+6+k;
	//
	big_ks = bs - bt + 4 + big_k;
	big_kt = bs - bt + 6 + big_k;

	if (quadratic_number_power_product::info > 3) {
		std::cout << "quadratic_number_power_product::rgcd::" << std::endl;
		std::cout << "COMPUTING relative " << big_ks << " approximation s" << std::endl;
		std::cout << "COMPUTING relative " << big_kt << " approximation t" << std::endl;
	}

	// ks += -b_value(s)+2
	// kt += -b_value(t)+2
	//
	big_ks -= (bs - 2);
	big_kt -= (bt - 2);

	if (big_ks.longify(ks)) {
		lidia_error_handler ("quadratic_number_power_product::rgcd",
				     "precision overflow.");
		return;
	}
	if (big_kt.longify(kt)) {
		lidia_error_handler ("quadratic_number_power_product::rgcd",
				     "precision overflow.");
		return;
	}

	if (quadratic_number_power_product::info > 3) {
		std::cout << "quadratic_number_power_product::rgcd::" << std::endl;
		std::cout << "COMPUTING absolute " << ks << " approximation s" << std::endl;
		std::cout << "COMPUTING absolute " << kt << " approximation t" << std::endl;
	}

	q1.get_absolute_Ln_approximation (s, ks);
	q2.get_absolute_Ln_approximation (t, kt);

	if (quadratic_number_power_product::info > 3) {
		std::cout << "quadratic_number_power_product::rgcd::" << std::endl;
		std::cout << "FOUND s = " << s << std::endl;
		std::cout << "FOUND t = " << t << std::endl;
	}

	// Find the convergent
	//
	// k =   s.get_exponent() - b_value(s.get_mantissa());
	// k -= (t.get_exponent() - b_value(t.get_mantissa()));
	//
	big_k = s.get_exponent();
	big_k -= b_value(s.get_mantissa());
	big_k -= t.get_exponent();
	big_k += b_value(t.get_mantissa());

	if (big_k.longify(k)) {
		lidia_error_handler ("quadratic_number_power_product::rgcd",
				     "precision overflow.");
		return;
	}

	num = s.get_mantissa();
	den = t.get_mantissa();

	if (k > 0)
		shift_left (num, num, k);
	else
		shift_left (den, den, -k);

	if (num.is_negative()) {
		s_is_negative = true;
		num.negate();
	}
	else
		s_is_negative = false;

	if (den.is_negative()) {
		t_is_negative = true;
		den.negate();
	}
	else
		t_is_negative = false;

	if (quadratic_number_power_product::info > 3) {
		std::cout << "quadratic_number_power_product::rgcd::" << std::endl;
		std::cout << "numerator = " << num << std::endl;
		std::cout << "denominator = " << den << std::endl;
		std::cout << "denominator bound = " << S << std::endl;
	}

	quadratic_number_power_product::cfrac_expansion_bounded_den (M1, M2, num, den, S);

	if (s_is_negative)
		M1.negate();

	if (t_is_negative)
		M2.negate();

	if (quadratic_number_power_product::info > 3) {
		std::cout << "quadratic_number_power_product::rgcd::" << std::endl;
		std::cout << "Found convergent " << M1 << " / " << M2 << std::endl;
	}

	// compute the x,y with 1 = x M1 + y M2
	//
	gcdM1M2 = xgcd (x, y, M1, M2);

	if (!gcdM1M2.is_one()) {
		lidia_error_handler ("quadratic_number_power_product::rgcd",
				     "gcd(M1, M2) != 1");
	}

	if (!(x*M1+y*M2 == 1)) {
		lidia_error_handler ("quadratic_number_power_product::rgcd",
				     "x*M1 + y*M2) != 1");
	}

	if (quadratic_number_power_product::info > 3) {
		std::cout << "quadratic_number_power_product::rgcd::" << std::endl;
		std::cout << "FOUND x = " << x << std::endl;
		std::cout << "FOUND y = " << y << std::endl;
	}
}



//
// ::cfrac_expansion_bounded_den
//
// Computes the last convergent conv_num/conv_den
// of the continued fraction expansion of num/den,
// of which absolute value of the denominator is
// less or equal to S, where S is positive.
//

void
quadratic_number_power_product::cfrac_expansion_bounded_den (bigint & conv_num,
							     bigint & conv_den,
							     bigint num,
							     bigint den,
							     const bigint & S)

{

	debug_handler ("quadratic_number_power_product",
		       "cfrac_expansion_bounded_den"
		       "(conv_num, conv_den, num, den, S)");

	bigint 	 r; // contains the rest
	// in the euclidian step

	bigint 	 q; // contains the quotient
	// in the euclidean step
	int 	 sign;

	bigint 	prev1_conv_num, prev1_conv_den;
	bigint 	prev2_conv_num, prev2_conv_den;

	//
	// sign * (num/den) = n/d,
	// num, den are positive
	//

	sign = 1;

	if (num.is_negative()) {
		num.negate();
		sign *= -1;
	}

	if (den.is_negative()) {
		den.negate();
		sign *= -1;
	}


	// compute the partial quotients
	// and the convergents

	// first convergent: q[0]/1

	if (S.is_negative())
		lidia_error_handler("quadratic_number_power_product"
				    "::cfrac_expansion_bounded_den",
				    "Negative upper bound.");

	if (S.is_zero())
		lidia_error_handler("quadratic_number_power_product"
				    "::cfrac_expansion_bounded_den",
				    "Upper bound is zero.");


	// first quotient q
	// first convergent: conv_num / d_n
	// first convergent: q / 1

	div_rem (q, r, num, den);
	num = den;
	den = r;

	conv_num = q;
	conv_den.assign_one();

	if (quadratic_number_power_product::info > 3) {
		std::cout << "quadratic_number_power_product::";
		std::cout << "cfrac_expansion_bounded_den::" << std::endl;
		std::cout << "first convergent: " << conv_num << " / ";
		std::cout << conv_den << std::endl;
	}

	if (!den.is_zero()) {
		// store previous convergent

		prev1_conv_num = conv_num;
		prev1_conv_den = conv_den;

		conv_num = q;

		// second quotient q

		div_rem (q, r, num, den);
		num = den;
		den = r;

		// second convergent conv_num / conv_den
		// second convergent: (q[0]q[1]+1)/q[1]

		LiDIA::multiply (conv_num, conv_num, q);
		inc (conv_num);

		conv_den = q;

		if (quadratic_number_power_product::info > 3) {
			std::cout << "quadratic_number_power_product::";
			std::cout << "cfrac_expansion_bounded_den::" << std::endl;
			std::cout << "second convergent: " << conv_num << " / ";
			std::cout << conv_den << std::endl;
		}

		// stop ?

		if (conv_den > S) {
			conv_num = prev1_conv_num;
			conv_den = prev1_conv_den;
		}
		else if (!den.is_zero()) {
			do {
				prev2_conv_num = prev1_conv_num;
				prev2_conv_den = prev1_conv_den;

				prev1_conv_num = conv_num;
				prev1_conv_den = conv_den;

				// next quotient q

				div_rem (q, r, num, den);
				num = den;
				den = r;

				// next convergent conv_num / conv_den
				// conv_num(k) = q(k) * conv_num(k-1) + conv_num(k-2)
				// conv_den(k) = q(k) * conv_den(k-1) + conv_den(k-2)

				LiDIA::multiply (conv_num, q, prev1_conv_num);
				LiDIA::multiply (conv_den, q, prev1_conv_den);

				add (conv_num, conv_num, prev2_conv_num);
				add (conv_den, conv_den, prev2_conv_den);

				if (quadratic_number_power_product::info > 3) {
					std::cout << "quadratic_number_power_product::";
					std::cout << "cfrac_expansion_bounded_den::" << std::endl;
					std::cout << "convergent: " << conv_num << " / ";
					std::cout << conv_den << std::endl;
				}
			} while (!den.is_zero() && conv_den <= S);

			if (conv_den > S) {
				conv_num = prev1_conv_num;
				conv_den = prev1_conv_den;
			}
		}
	}

	if (sign == -1)
		conv_num.negate();
}



//
//  estimate_pp_accuracy
//

long
quadratic_number_power_product::estimate_pp_accuracy (const matrix< bigint > & M,
						      long m,
						      xbigfloat c
	)

	//
	// estimate accuracy of power product
	//
	// Returns k such that k heuristically is an upper bound on the
	// largest absolute accuracy needed in the rgcd algorithm.
	// 2^m should be a lower bound on the regulator and
	// c an upper bound on the absolute values of Ln(q), where q
	// runs over all principal ideal generators, which form the
	// bases in the power products of the units.
	//
	//  k = b(c) + max(j) b(sum |e_ij|) - 2m + 9
	//

{
	debug_handler ("quadratic_number_power_product",
		       "estimate_pp_accuracy"
		       "(const matrix< bigint > &, long, xbigfloat)");

	bigint col_sum, max_col_sum;
	bigint big_k;
	long   k;

	lidia_size_t i, j;
	lidia_size_t nofr, nofc;


	// initialize
	//
	nofr = M.get_no_of_rows();
	nofc = M.get_no_of_columns();
	max_col_sum.assign_zero();


	// determine column sum norm
	//
	for (j = 0; j < nofc; j++) {
		col_sum.assign_zero();
		for (i = 0; i < nofr; i++) {
			// add absolute value
			if (M.member(i, j).is_negative())
				subtract(col_sum, col_sum, M.member(i, j));
			else
				add(col_sum, col_sum, M.member(i, j));
		}

		// compare max
		if (col_sum > max_col_sum)
			max_col_sum.assign(col_sum);
	}


	// determine precision
	//
	big_k = bigint(c.b_value()) + bigint(b_value(max_col_sum));
	big_k += 9 - 2 * bigint(m);

	if (big_k.longify(k)) {
		lidia_error_handler ("quadratic_number_power_product::"
				     "estimate_pp_accuracy",
				     "Precision overflow.");
	}

	return k;
}



//
//  compact_representation_of_unit
//

void
quadratic_number_power_product::compact_representation_of_unit ()
{
	debug_handler ("quadratic_number_power_product",
		       "compact_representation_of_unit()");

	// I = order of the unit.
	//
	quadratic_ideal I;
	I.assign_one(this->get_order());

	this->compact_representation(I);
}



//
//  Unit tests
//

bool
quadratic_number_power_product::could_be_unit_norm_test (int trials)
{
	debug_handler ("quadratic_number_power_product",
		       "could_be_unit_norm_test(int)");

	int i;
	bigint m, num, den;

	m.randomize(100);
	for (i = 0; i < trials; i++) {
		m = next_prime(m);
		this->norm_modulo(num, den, m);

		num.absolute_value();
		den.absolute_value();

		if (num != den)
			return false;
	}
	return true;
}



bool
quadratic_number_power_product::could_be_quasi_unit_close_test ()
{
	debug_handler ("quadratic_number_power_product",
		       "could_be_quasi_unit_close_test()");

	// Compute absolute 3 approx. l to Ln a
	//
	xbigfloat l;
	this->get_absolute_Ln_approximation(l, 3);

	quadratic_order O = this->get_order();
	return O.could_be_regulator_multiple(l);
}



int
quadratic_number_power_product::could_be_unit_refinement_test ()
{
	debug_handler ("quadratic_number_power_product",
		       "could_be_unit_refinement_test()");

	quadratic_ideal_power_product p;
	p.assign((this->get_basis()).get_basis(), this->get_exponents());
	int overorder = p.factor_refinement();

	if (!p.is_gcd_free()) {
		lidia_error_handler("quadratic_number_power_product::"
				    "could_be_unit_refinement_test()",
				    "p not gcd free.");
		return 0;
	}

	lidia_size_t i, n;

	n = p.get_no_of_components();
	for (i = 0; i < n; i++) {
		if (!p.get_exponent(i).is_zero())
			return 0;
	}

	if (overorder)
		return -1;
	else
		return 1;
}



//
//  Fundamental unit computation.
//

bool
quadratic_number_power_product::fundamental_unit_if_regulator_leq (const xbigfloat & x,
								   const quadratic_order & O)
{
	debug_handler ("quadratic_number_power_product",
		       "fundamental_unit_if_regulator_leq("
		       "const xbigfloat&, const quadratic_order&)");

	quadratic_number_standard mu;
	quadratic_number_power_product beta;
	quadratic_ideal I;
	xbigfloat b;
	long k;

	if (x.is_negative())
		return false;

	I.assign(O);
	beta.assign_one(O);
	b.assign_zero();
	k = 4;

	do {
		I.inverse_rho(mu);
		beta *= mu;
		b = beta.absolute_Ln_approximation(k);
	} while (I != O && b < x);

	if (I == O) {
		this->assign(O, b, 1);
		return true;
	}
	else {
		I.inverse_rho(mu);
		if (I == O) {
			beta *= mu;
			b = beta.absolute_Ln_approximation(k);
			this->assign(O, b, 1);
			return true;
		}
		else
			return false;
	}
}



void
quadratic_number_power_product::fundamental_unit (const quadratic_number_power_product & beta,
						  const bigint & h)
{
	debug_handler ("quadratic_number_power_product",
		       "fundamental_unit(const qu&, const bigint&)");

	quadratic_ideal I;
	quadratic_number_power_product gamma;
	quadratic_number_standard mu;
	xbigfloat a, c, d, l, x, tmp, t;
	bigint B, Delta;
	long k, LRB;
	bool new_unit;
	lidia_size_t i;
	unsigned long p;

	// Note: We need a reference here, because a
	// local order is deleted from the qo_list at
	// the end of the scope. So "this" below, initialized
	// with the local order has no valid order at the end
	// of the scope anymore.
	//
	const quadratic_order & O = beta.get_order();
	Delta = O.discriminant();

	// Test for regulator < 3/2 ln Delta + 2^{-2}
	//
	// x = (1+1/16) * b(Delta) + 2^{-2}
	//
	tmp.assign_one();
	shift_right(tmp, tmp, 4);
	inc(tmp);
	LiDIA::multiply(x, tmp, xbigfloat(b_value(Delta)));

	tmp.assign_one();
	shift_right(tmp, tmp, 2);
	x += tmp;

#ifdef QU_FU_DEBUG
	std::cout << "Testing for regulator smaller than ";
	x.print_as_bigfloat();
	std::cout << " ... " << std::flush;
#endif

	if (this->fundamental_unit_if_regulator_leq(x, O))
		return;

#ifdef QU_FU_DEBUG
	std::cout << "DONE. Regulator large enough. Going on." << std::endl;
#endif
	// Regulator > 2^{LRB}. Needed for test for 1 below.
	//
	LRB = x.b_value()-1;
#ifdef QU_FU_DEBUG
	std::cout << "LRB = " << LRB << std::endl;
#endif

	// If regulator is large enough, start searching
	// for order in distance (Ln beta)/p.
	//
	l = O.relative_L1chi_approximation(1);
	LiDIA::sqrt(d, Delta, 4);

	// l = l * d/2
	//
	shift_right(d, d, 1);
	l *= d;

	gamma = beta;

	do {
		new_unit = false;
		*this = gamma;

#ifdef QU_FU_DEBUG
		std::cout << "current unit is " << *this << std::endl;
		std::cout << "Starting search for new unit." << std::endl;
#endif
		// Let alpha = *this.
		//
		// Upper bound B >= m
		//
		a = this->absolute_Ln_approximation(3);

		// x = (|a|+2^(-3)) * h;
		//
		x = a;
		x.absolute_value();
		tmp.assign_one();
		shift_right(tmp, tmp, 3);
		tmp += x;
		LiDIA::multiply(x, tmp, xbigfloat(h));

		// divide (x,x,l*d/2,4);
		//
		LiDIA::divide(x, x, l, 4);

		// B = floor(x*1.75);
		//
		tmp.assign_one();
		shift_right(tmp, tmp, 1);
		inc(tmp);
		shift_right(tmp, tmp, 1);
		inc(tmp);
		x *= tmp;
		floor(B, x);
		++B;

#ifdef QU_FU_DEBUG
		std::cout << "Upper bound for primes is B = " << B << std::endl;
#endif
		// Initialize prime array and accuracy
		// InitPrimes(P,B);
		//
		if (B.sizetify(i))
			lidia_error_handler("quadratic_number_power_product::fundamental_unit(qu, bigint)",
					    "Number of primes too large for lidia_size_t.");
		prime_list P(i);

		// k = 5+b(a)
		//
		k = 5 + a.b_value();

		// Test (Ln alpha) / p for each prime p
		// and stop if a new unit is found.
		//
		for (i = P.get_number_of_primes()-1; i >= 0 && !new_unit; i--) {
			p = P[i];
#ifdef QU_FU_DEBUG
			std::cout << "Testing prime p = " << p << std::endl;
#endif
			// absolute 3-approx. to (Ln alpha)/p
			//
			LiDIA::divide (t, a, p, k-b_value(p));

			// search for order
			//
			I.assign(O);
			I.order_close(gamma, c, t, 4);
#ifdef QU_FU_DEBUG
			std::cout << "t = "; t.print_as_bigfloat(); std::cout << std::endl;
#endif
			I.rho(mu);
			if (I == O) {
				gamma *= mu; new_unit = true;
#ifdef QU_FU_DEBUG
				std::cout << "Found with I.rho(); I = " << I << std::endl;
#endif
			}
			else {
				I.inverse_rho();
				if (I == O) {
					new_unit = true;
#ifdef QU_FU_DEBUG
					std::cout << "Found with I; I = " << I << std::endl;
#endif
				}
				else {
					I.inverse_rho(mu);
					if (I == O) {
						gamma *= mu; new_unit = true;
#ifdef QU_FU_DEBUG
						std::cout << "Found with I.inverse_rho(); I = " << I << std::endl;
#endif
					}
				}
			}

			if (new_unit) {
				if (gamma.is_rational(LRB)) {
					new_unit = false;
#ifdef QU_FU_DEBUG
					std::cout << "|Unit| is one. Ignoring." << std::endl;
#endif
				}
			}
		}
	} while (new_unit);

	this->compact_representation_of_unit();
}



//
//   input operator
//

std::istream & operator >> (std::istream & in,
			    quadratic_number_power_product & x)
{
	debug_handler("quadratic_number_power_product",
		      "operator >> (std::istream&, quadratic_number_power_product&)");

	x.read(in);
	return in;
}



//
//  input function
//

void
quadratic_number_power_product::read (std::istream & in)
{
	debug_handler ("quadratic_number_power_product",
		       "read(std::istream&)");

	quadratic_number_power_product_basis b;
	char c;

	// read white spaces
	in >> c;
	while (c == ' ') in >> c;

	// read '('
	if (c != '(') {
		lidia_error_handler("quadratic_number_power_product::read",
				    "'(' expected");
	}
	else {
		// read basis
		in >> b;

		// read white spaces
		in >> c;
		while (c == ' ') in >> c;

		// read ','
		if (c != ',') {
			lidia_error_handler("quadratic_number_power_product::read",
					    "(b read; ',' expected.");
		}
		else {
			// read exponents
			in >> exp;

			// read white spaces
			in >> c;
			while (c == ' ') in >> c;

			// read ')'
			if (c != ')') {
				lidia_error_handler ("quadratic_number_power_product::read",
						     "(b, e read; ')' expected.");
			}
		}
		// end if ( c != ',' )
	}
	// end if ( c != '(' )

	basis = quadratic_number_power_product::basis_list.insert(b);
	verify_preconditions("quadratic_number_power_product::read(std::istream)");
}



//
//   output operator
//

std::ostream &
operator << (std::ostream & out,
	     const quadratic_number_power_product & b)
{
	debug_handler ("quadratic_number_power_product",
		       "operator << (std::ostream &, const qnpp &)");
	b.write(out);
	return out;
}



//
//  output function
//

void
quadratic_number_power_product::write (std::ostream & out) const
{
	debug_handler ("quadratic_number_power_product",
		       "write(std::ostream&)");
	out << "(" << basis->get_mod() << ", " << exp << ")";
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
