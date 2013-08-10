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
//	Author	: Stefan Neis (SN)
//	Changes	: See CVS log
//
//==============================================================================================


#ifdef HAVE_CONFIG_H
# include	"config.h"
#endif
#include	"LiDIA/alg_number.h"
#include	"LiDIA/debug.h"
#include	"LiDIA/base_vector.h"
#include	<cctype>



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



#ifdef LIDIA_DEBUG
int alg_number::count = 0;
#endif

// Constructors & destructor:
alg_number::alg_number(const nf_base * O1)
	: den(1),
	  coeff(O1->degree(), O1->degree()),
	  O(const_cast<nf_base *>(O1)) // Initialised with 0
{
	debug_handler_c("alg_number", "in alg_number(const nf_base *)", 1,
			count++;
			std::cout << "\nNow we have " << count << " alg_numbers!\n");
	O->inc_ref();
}



alg_number::alg_number(const bigint & i, const nf_base * O1)
	: den(1),
	  coeff((const_cast<nf_base *>(O1))->get_one() * i),
	  O(const_cast<nf_base *>(O1))
{
	debug_handler_c("alg_number", "in alg_number(const bigint &, nf_base *)",
			1, count++;
			std::cout << "\nNow we have " << count << " alg_numbers!\n");
	O->inc_ref();
}



alg_number::alg_number(const base_vector< bigint > & a,
		       const bigint &i, const nf_base * O1)
	: den(i),
	  coeff(a),
	  O(const_cast<nf_base *>(O1))

{
	if (degree() != a.size()) {
		lidia_error_handler("alg_number", "in alg_number "
				    "(const base_vector< bigint > &, const bigint &, "
				    "const nf_base *)::base_vector of wrong size");
		assign(alg_number());
	}
	debug_handler_c("alg_number",
			"in alg_number"
			"(const base_vector &, const bigint &, const nf_base *)", 1,
			count++;
			std::cout << "\nNow we have " << count << " alg_numbers!\n");
	normalize();
	O->inc_ref();
}



alg_number::alg_number(const bigint * a, const bigint &i, const nf_base * O1)
	: den(i),
	  coeff(a, O1->degree()),
	  O(const_cast<nf_base *>(O1))

{
	debug_handler_c("alg_number",
			"in alg_number(bigint *, const bigint &, const nf_base *)",
			1, count++;
			std::cout << "\nNow we have " << count << " alg_numbers!\n");
	normalize();
	O->inc_ref();
}



alg_number::alg_number(const alg_number &a)
	: den(a.den),
	  coeff(a.coeff),
	  O(a.O)

{
	debug_handler_c("alg_number", "in alg_number(const a_n &)", 1,
			count++;
			std::cout << "\nNow we have " << count << " alg_numbers!\n");
	O->inc_ref();
}



alg_number::~alg_number()
{
	debug_handler_c("alg_number", "in ~alg_number()", 1,
			count--;
			std::cout << "\nNow we have only " << count << " alg_numbers!\n");
	O->dec_ref();
}



// member-functions
bigfloat alg_number::get_conjugate(lidia_size_t j) const
{
	if (j< 1 || j > degree()) {
		lidia_error_handler("field", "get_conjugate(lidia_size_t i, lidia_size_t j)"
				    "::j out of range (must be between 1 and n)");
		return bigfloat();
	}
	math_vector< bigfloat > tmp1(degree(), degree()), tmp(degree(), degree());
	(static_cast<nf_base *>(O))->get_conjugates().get_row_vector(tmp1, --j);
	for (register lidia_size_t i = 0; i < degree(); i++)
		tmp[i].assign(coeff[i]);

	bigfloat res;
	multiply(res, tmp1, tmp);
	divide(res, res, den);
	return res;
}



math_vector< bigfloat > alg_number::get_conjugates() const
{
	math_vector< bigfloat > tmp(degree(), degree());
	for (register lidia_size_t i = 0; i < degree(); i++)
		tmp[i].assign(coeff[i]);
	multiply(tmp, (static_cast<nf_base *>(O))->get_conjugates(), tmp);
	divide(tmp, tmp, bigfloat(den));
	return tmp;
}



bigint_matrix rep_matrix(const alg_number & a)
	// RepMatrix of a * denominator
{
	debug_handler("alg_number", "in function rep_matrix(const alg_number &)");
	bigint_matrix A(a.degree(), a.degree()); // initialised with 0 !
	// construct the matrix using either MT or the polynomial
	nf_base * O = a.O;
	if (O->using_necessary()) {
		// If necessary construct MT
		if (!(O->table_computed())) {
			O->compute_table();
		}
		debug_handler("rep_matrix", "Constructing Matrix");

		base_vector< bigint > MT_vec(a.degree() * a.degree(),
					     a.degree() * a.degree());
		for (register lidia_size_t k = 0; k < a.degree(); k++) {
			O->table.get_column_vector(MT_vec, k);
			lidia_size_t l = 0;
			for (register lidia_size_t i = 0; i < a.degree(); i++) {
				for (register lidia_size_t j = 0; j < i; j++, l++) {
					A.sto(k, i, A.member(k, i) + MT_vec[l] * a.coeff.member(j));
					A.sto(k, j, A.member(k, j) + MT_vec[l] * a.coeff.member(i));
				}
				A.sto(k, i, A.member(k, i) + MT_vec[l++] *a.coeff.member(i));
			}
		}
		debug_handler("rep_matrix", "Constructed Matrix");
	}
	else {
		debug_handler("rep_matrix", "Constructing Matrix by Polynomials");
		bigint * tmp = new bigint[2];
		tmp[1] = 1;

		polynomial< bigint > x(tmp, 1);
		delete[] tmp;

		polynomial< bigint > p(a.coeff.get_data_address(), a.degree()-1);

		for (register lidia_size_t i = 0; i < a.degree();
		     i++, p = p*x%(O->f)) {
			register lidia_size_t l = 0;
			for (; l <= p.degree(); l++)
				A.sto(l, i, p[l]);
			for (; l < a.degree(); l++)
				A.sto(l, i, 0);
		}
		debug_handler("rep_matrix", "Constructed Matrix by Polynomials");
	}
	return A;
}



bool alg_number::is_zero() const
{
	for (register lidia_size_t i = 0; i < degree(); i++)
		if (!(coeff.member(i).is_zero())) return false;
	return true;
}



bool alg_number::is_one() const
{
	if (!den.is_one()) return false;
	if (coeff != O->get_one()) return false;
	return true;
}



void alg_number::normalize()
{
	bigint d(den);
	register lidia_size_t i;
	for (i = 0; i < degree(); i++)
		d = gcd (coeff.member(i), d);
	if (den.is_negative()) d.negate();
	if (!d.is_one()) {
		divide(den, den, d);
		for (i = 0; i < degree(); i++)
			coeff[i] /= d;
	}
}



void alg_number::negate()
{
	debug_handler("alg_number", "in member - function negate()");
	for (lidia_size_t i = 0; i < degree(); i++)
		coeff[i].negate();
}



void alg_number::multiply_by_2()
{
	debug_handler("alg_number", "in member - function multiply_by_2()");
	if (den.is_even())
		den.divide_by_2();
	else
		for (lidia_size_t i = 0; i < degree(); i++)
			coeff[i].multiply_by_2();
}



void alg_number::divide_by_2()
{
	debug_handler("alg_number", "in member - function divide_by_2()");
	den.multiply_by_2();
	normalize();
}



void alg_number::invert()
{
	debug_handler("alg_number", "in member - function invert()");
	divide(*this, bigint(1), *this);
}



void alg_number::assign_zero()
{
	for (register lidia_size_t i = 0; i < degree(); i++)
		coeff[i].assign_zero();
	den = 1;
}



void alg_number::assign_one()
{
	debug_handler("alg_number", "in member - function assign_one()");
	coeff.assign(O->get_one());
	den = 1;
}



void alg_number::assign(const bigint & a)
{
	multiply(coeff, O->get_one(), a);
	den = 1;
}



void alg_number::assign(const alg_number & a)
{
	debug_handler("alg_number", "in member - function assign(a_n &)");
	if (O != a.O) {
		O->dec_ref();
		O = a.O;
		O->inc_ref();
	}
	coeff = a.coeff;
	den = a.den;
}



// Procedural versions:
void add(alg_number &c, const alg_number &a, const alg_number &b)
{
	if (a.O != b.O) {
		lidia_error_handler("alg_number", "add(...)::addition of numbers from "
				    "different orders is not yet implemented");
		return;
	}
	if (c.O != a.O) {
		c.O->dec_ref();
		c.O = a.O;
		c.O->inc_ref();
	}
	bigint d = gcd(a.den, b.den);
	bigint e = a.den / d;

	add(c.coeff, (a.coeff) * (b.den / d), (b.coeff) * e);
	multiply(c.den, e, b.den);
	c.normalize();
}



void subtract(alg_number &c, const alg_number &a, const alg_number &b)
{
	if (a.O != b.O) {
		lidia_error_handler("alg_number", "subtract(...)::subtraction of numbers "
				    "from different orders is not yet implemented");
		return;
	}
	if (c.O != a.O) {
		c.O->dec_ref();
		c.O = a.O;
		c.O->inc_ref();
	}
	bigint d = gcd(a.den, b.den);
	bigint e = a.den / d;

	subtract(c.coeff, a.coeff * (b.den / d), b.coeff * e);
	multiply(c.den, e, b.den);
	c.normalize();
}



void multiply(alg_number &c, const alg_number &a, const alg_number &b)
{
	debug_handler("alg_number", "in function multiply(a_n &, const a_n &, const a_n &)");
	if (&a == &b) {
		square(c, a);
		return;
	}
	if (a.O != b.O) {
		lidia_error_handler("alg_number", "multiply(...)::"
				    "multiplication of numbers from "
				    "different orders is not yet implemented");
		return;
	}
	if (c.O != a.O) {
		c.O->dec_ref();
		c.O = a.O;
		c.O->inc_ref();
	}
	debug_handler_c("alg_number",
			"in function multiply(a_n &, const a_n &, const a_n &)", 0,
			std::cout << "multiply " << a << " by " << b << std::endl);
	nf_base * O = a.O;
	multiply(c.den, a.den, b.den);
	if (O->using_necessary()) {
		// If necessary construct MT
		if (!(O->table_computed())) {
			O->compute_table();
		}
		bigint tmp1, tmp2, tmp3;
		base_vector< bigint > MT_vec(a.degree() * a.degree(),
					     a.degree() * a.degree());
		if ((&c != &a) && (&c != &b))
			for (register lidia_size_t k = 0; k < a.degree(); k++) {
				c.coeff[k].assign_zero();
				O->table.get_column_vector(MT_vec, k);
				lidia_size_t l = 0;
				for (register lidia_size_t i = 0; i < a.degree(); i++) {
					tmp1.assign_zero();
					tmp2.assign_zero();
					for (register lidia_size_t j = 0; j < i; j++, l++) {
						// tmp1 += MT_vec[l] * b.coeff.member(j);
						multiply(tmp3, MT_vec[l], b.coeff.member(j));
						add(tmp1, tmp1, tmp3);

						// tmp2 += MT_vec[l] * a.coeff.member(j);
						multiply(tmp3, MT_vec[l], a.coeff.member(j));
						add(tmp2, tmp2, tmp3);
					}
					// tmp1 *= a.coeff.member(i);
					multiply(tmp1, tmp1, a.coeff.member(i));

					// tmp2 += MT_vec[l] * a.coeff.member(i);
					multiply(tmp3, MT_vec[l++], a.coeff.member(i));
					add(tmp2, tmp2, tmp3);

					// tmp2 *= b.coeff.member(i);
					multiply(tmp2, tmp2, b.coeff.member(i));

					// c.coeff[k] += tmp1+tmp2;
					add(tmp1, tmp1, tmp2);
					add(c.coeff[k], c.coeff[k], tmp1);
				}
			}
		else {
			debug_handler_l("alg_number",
					"in function multiply(a_n &, const a_n &, const a_n &):"
					" use MT and overwrite argument", 0);
			math_vector< bigint > tmp(a.degree(), a.degree());
			for (register lidia_size_t k = 0; k < a.degree(); k++) {
				tmp[k].assign_zero();
				O->table.get_column_vector(MT_vec, k);
				lidia_size_t l = 0;
				for (register lidia_size_t i = 0; i < a.degree(); i++) {
					tmp2.assign_zero();
					tmp1.assign_zero();
					for (register lidia_size_t j = 0; j < i; j++, l++) {
						// tmp1 += MT_vec[l] * b.coeff.member(j);
						multiply(tmp3, MT_vec[l], b.coeff.member(j));
						add(tmp1, tmp1, tmp3);

						// tmp2 += MT_vec[l] * a.coeff.member(j);
						multiply(tmp3, MT_vec[l], a.coeff.member(j));
						add(tmp2, tmp2, tmp3);
					}
					// tmp1 *= a.coeff.member(i);
					multiply(tmp1, tmp1, a.coeff.member(i));

					// tmp2 += MT_vec[l] * a.coeff.member(i);
					multiply(tmp3, MT_vec[l++], a.coeff.member(i));
					add(tmp2, tmp2, tmp3);

					// tmp2 *= b.coeff.member(i);
					multiply(tmp2, tmp2, b.coeff.member(i));

					// tmp[k] += + tmp2+tmp1;
					add(tmp1, tmp1, tmp2);
					add(tmp[k], tmp[k], tmp1);
				}
			}
			swap(c.coeff, tmp);
		}
	}
	else {
		// Polynomarithmetik verwenden!
		bigint * tmp;

		polynomial< bigint > pa(a.coeff.get_data_address(),
					a.degree()-1);

		polynomial< bigint > pb(b.coeff.get_data_address(),
					a.degree()-1);

		polynomial< bigint > f = pa*pb%(O->f);
		register lidia_size_t l = 0;
		for (; l <= f.degree(); l++)
			c.coeff[l] = f[l];
		for (; l < a.degree(); l++)
			c.coeff[l] = 0;
	}
	c.normalize();
	debug_handler_c("alg_number",
			"in function multiply(a_n &, const a_n &, const a_n &)", 0,
			std::cout << "result is " << c << std::endl);
}



void divide(alg_number &c, const alg_number &a, const alg_number &b)
{
	debug_handler("alg_number",
		      "in function divide(a_n &, const a_n &, const a_n &)");
	if (a.O != b.O) {
		lidia_error_handler("alg_number", "divide(...)::division of numbers from "
				    "different orders is not yet implemented");
		return;
	}

	bigint_matrix A = rep_matrix(b);

	// solve linear equation system
	bigint_matrix tmp(a.degree(), 1);
	bigint_matrix tmp2(a.degree()+1, 1);
	bigint * temp;
	tmp.sto_column_vector(a.coeff, a.degree(), 0);
	tmp2.reginvimage(A, tmp);

	// initialise result with solution
	if (c.O != a.O) {
		c.O->dec_ref();
		c.O = a.O;
		c.O->inc_ref();
	}
	c.coeff.set_data(temp = tmp2.get_column(0), a.degree());

	// get denominator
	bigint d1 = gcd(a.den, b.den);
	bigint e = b.den / d1;
	bigint d2 = gcd(temp[a.degree()], e);

	divide(c.den, a.den, d1);
	c.normalize();
	//These factors won't need normalization:
	multiply(c.den, c.den, temp[a.degree()] / d2);
	multiply(c.coeff, c.coeff, e / d2);
	delete[] temp;
}



void add(alg_number &c, const alg_number &a, const bigint & b)
{
	if (c.O != a.O) {
		c.O->dec_ref();
		c.O = a.O;
		c.O->inc_ref();
	}
	alg_number to_add((c.O->get_one()) * b * a.den, 1, c.O);
	add(c.coeff, a.coeff, to_add.coeff);
	c.den = a.den;
	c.normalize();
}



void subtract(alg_number &c, const alg_number &a, const bigint & b)
{
	if (c.O != a.O) {
		c.O->dec_ref();
		c.O = a.O;
		c.O->inc_ref();
	}
	alg_number to_sub((c.O->get_one()) * b * a.den, 1, c.O);
	subtract(c.coeff, a.coeff, to_sub.coeff);
	c.den = a.den;
	c.normalize();
}



void multiply(alg_number &c, const alg_number &a, const bigint & b)
{
	bigint d = gcd(a.den, b);

	if (c.O != a.O) {
		c.O->dec_ref();
		c.O = a.O;
		c.O->inc_ref();
	}
	divide(c.den, a.den, d);
	multiply(c.coeff, a.coeff, (b / d));
}



void divide(alg_number &c, const alg_number &a, const bigint & b)
{
	debug_handler("alg_number",
		      "in function divide(a_n &, const a_n &, const bigint &)");

	if (c.O != a.O) {
		c.O->dec_ref();
		c.O = a.O;
		c.O->inc_ref();
	}
	c.coeff = a.coeff;

	if (&c != &a) {
		c.den = b;
		c.normalize();
		multiply(c.den, c.den, a.den);
	}
	else {
		bigint tmp = a.den;
		c.den = b;
		c.normalize();
		multiply(c.den, c.den, tmp);
	}
}



void add(alg_number &c, const bigint & b, const alg_number &a)
{
	if (c.O != a.O) {
		c.O->dec_ref();
		c.O = a.O;
		c.O->inc_ref();
	}
	alg_number to_add((c.O->get_one()) * b * a.den, 1, c.O);
	add(c.coeff, a.coeff, to_add.coeff);
	c.den = a.den;
	c.normalize();
}



void subtract(alg_number &c, const bigint & b, const alg_number &a)
{
	if (c.O != a.O) {
		c.O->dec_ref();
		c.O = a.O;
		c.O->inc_ref();
	}
	alg_number sub((c.O->get_one()) * b * a.den, 1, c.O);
	subtract(c.coeff, sub.coeff, a.coeff);
	c.den = a.den;
	c.normalize();
}



void multiply(alg_number &c, const bigint &b, const alg_number &a)
{
	bigint d = gcd(a.den, b);

	if (c.O != a.O) {
		c.O->dec_ref();
		c.O = a.O;
		c.O->inc_ref();
	}
	divide(c.den, a.den, d);
	multiply(c.coeff, a.coeff, b / d);
}



void divide(alg_number &c, const bigint &b, const alg_number &a)
{
	debug_handler("alg_number",
		      "in function divide(a_n &, const bigint &, const a_n &)");
	nf_base * O = a.O;
	alg_number dividend(O->get_one() * b, 1, O);

	divide(c, dividend, a);
}



void power(alg_number &c, const alg_number &a, const bigint & b)
{
	debug_handler("alg_number",
		      "in function power(a_n &, const a_n &, const bigint &)");
	bigint expo;
	if (c.O != a.O) {
		c.O->dec_ref();
		c.O = a.O;
		c.O->inc_ref();
	}
	alg_number multiplier(c.O);

	if (b.is_negative())
		power(c, inverse(a), -b);
	else if (b.is_zero() || a.is_one())
		c.assign_one();
	else {
		expo.assign(b);
		multiplier.assign(a);
		c.assign_one();
		while (expo.is_gt_zero()) {
			if (!expo.is_even()) {
				debug_handler_c("alg_number", "in function power"
						"(a_n &, const a_n &, const bigint &)", 4,
						std::cout << "multiply " << c << " by " << multiplier;
						std::cout << std::endl << std::flush);
				multiply(c, c, multiplier);
			}
			debug_handler_c("alg_number", "in function power"
					"(a_n &, const a_n &, const bigint &)", 4,
					std::cout << "square " << multiplier << std::endl << std::flush);
			square(multiplier, multiplier);
			expo.divide_by_2();
			debug_handler_c("alg_number", "in function power"
					"(a_n &, const a_n &, const bigint &)", 4,
					std::cout << "exponent is now " << expo << std::endl << std::flush);

		}
	}
}



void power_mod_p(alg_number &c, const alg_number &a, const bigint & b,
		 const bigint &p)
{
	debug_handler("alg_number",
		      "in function power(a_n &, const a_n &, const bigint &, "
		      "const bigint &)");
	bigint expo;
	if (c.O != a.O) {
		c.O->dec_ref();
		c.O = a.O;
		c.O->inc_ref();
	}
	alg_number multiplier(c.O);

	if (b.is_negative())
		power_mod_p(c, inverse(a), -b, p);
	else if (b.is_zero() || a.is_one())
		c.assign_one();
	else {
		expo.assign(b);
		multiplier.assign(a);
		for (register lidia_size_t i = 0; i < multiplier.coeff.size(); i++)
			multiplier.coeff[i] %= p;
		c.assign_one();
		while (expo.is_gt_zero()) {
			if (!expo.is_even()) {
				debug_handler_c("alg_number", "in function power"
						"(a_n &, const a_n &, const bigint &)", 4,
						std::cout << "multiply " << c << " by " << multiplier;
						std::cout << std::endl << std::flush);
				multiply(c, c, multiplier);
				for (register lidia_size_t i = 0; i < c.coeff.size(); i++)
					c.coeff[i] %= p;
			}
			debug_handler_c("alg_number", "in function power"
					"(a_n &, const a_n &, const bigint &)", 4,
					std::cout << "square " << multiplier << std::endl << std::flush);
			square(multiplier, multiplier);
			for (register lidia_size_t i = 0; i < multiplier.coeff.size(); i++)
				multiplier.coeff[i] %= p;
			expo.divide_by_2();
			debug_handler_c("alg_number", "in function power"
					"(a_n &, const a_n &, const bigint &)", 4,
					std::cout << "exponent is now " << expo << std::endl << std::flush);

		}
	}
}



// Comparision:
// By now, only comparision of numbers over the same order is implemented.
bool operator == (const alg_number & a, const alg_number & b)
{
	if (a.O != b.O) {
		lidia_error_handler("alg_number", "operator == :: "
				    "You tried to compare numbers over "
				    "different orders!");
		return false;
	}
	if (a.den != b.den) return false;
	if (a.coeff != b.coeff) return false;
	return true;
}



// Some number-theoretic functions:
bigrational norm(const alg_number & a)
{
	bigint normdenominator;
	power(normdenominator, a.den, bigint(a.degree()));
	return bigrational(rep_matrix(a).det(), normdenominator);
}



bigrational trace(const alg_number & a)
{
	debug_handler("alg_number", "in function trace(const alg_number &)");
	bigint_matrix A(rep_matrix(a));
	bigint tmp;
	A.trace(tmp);
	return bigrational(tmp, a.den);
}



polynomial< bigint > charpoly(const alg_number & a)
{
	debug_handler_l("alg_number",
			"in friend function `charpoly(const alg_number &)'", 5);
	bigint * tmp = rep_matrix(a).charpoly();
	debug_handler_l("alg_number",
			"charpoly(const alg_number &)::computed poly of matrix", 5);
	polynomial< bigint > p(tmp, a.degree());
	delete[] tmp;
	bigint multiplier = a.den;
	for (register lidia_size_t i = 1; i <= a.degree(); i++)
		multiply(p[i], p[i], multiplier);
	multiply(multiplier, multiplier, a.den);
	return p;
}



// Other functions:
void negate(alg_number & b, const alg_number & a)
{
	b.assign(alg_number(-a.coeff, a.den, a.O));
}



void multiply_by_2(alg_number & c, const alg_number & a)
{
	c.assign(a);
	c.multiply_by_2();
}



void divide_by_2(alg_number & c, const alg_number & a)
{
	c.assign(a);
	c.divide_by_2();
}



void invert(alg_number &c, const alg_number & a)
{
	c.assign(a);
	c.invert();
}



alg_number inverse(const alg_number & a)
{
	alg_number c = a;
	c.invert();
	return c;
}



void square(alg_number & c, const alg_number & a)
{
	debug_handler("alg_number",
		      "in function square(a_n &, const a_n &)");
	if (c.O != a.O) {
		c.O->dec_ref();
		c.O = a.O;
		c.O->inc_ref();
	}
	debug_handler_c("alg_number",
			"in function square(a_n &, const a_n &)", 0,
			std::cout << "squaring " << a << std::endl);
	nf_base * O = a.O;
	square(c.den, a.den);
	if (O->using_necessary()) {
		// If necessary construct MT
		if (!(O->table_computed())) {
			O->compute_table();
		}
		bigint tmp1, tmp2;
		base_vector< bigint > MT_vec(a.degree() * a.degree(),
					     a.degree() * a.degree());
		if (&c != &a)
			for (register lidia_size_t k = 0; k < a.degree(); k++) {
				c.coeff[k].assign_zero();
				O->table.get_column_vector(MT_vec, k);
				lidia_size_t l = 0;
				for (register lidia_size_t i = 0; i < a.degree(); i++) {
					tmp1.assign_zero();
					for (register lidia_size_t j = 0; j < i; j++) {
						// tmp1 += MT_vec[l] * a.coeff.member(j);
						multiply(tmp2, MT_vec[l++], a.coeff.member(j));
						add(tmp1, tmp1, tmp2);
					}
					// tmp1 *= 2;
					tmp1.multiply_by_2();
					// tmp1 += MT_vec[l] * a.coeff.member(i);
					multiply(tmp2, MT_vec[l++], a.coeff.member(i));
					add(tmp1, tmp1, tmp2);
					// tmp1 *= a.coeff.member(i);
					multiply(tmp1, tmp1, a.coeff.member(i));

					// c.coeff[k] += tmp1;
					add(c.coeff[k], c.coeff[k], tmp1);
				}
			}
		else {
			debug_handler_l("alg_number",
					"in function multiply(a_n &, const a_n &, const a_n &):"
					" use MT and overwrite argument", 0);
			math_vector< bigint > tmp(a.degree(), a.degree());
			for (register lidia_size_t k = 0; k < a.degree(); k++) {
				tmp[k].assign_zero();
				O->table.get_column_vector(MT_vec, k);
				lidia_size_t l = 0;
				for (register lidia_size_t i = 0; i < a.degree(); i++) {
					tmp1.assign_zero();
					for (register lidia_size_t j = 0; j < i; j++) {
						// tmp1 += MT_vec[l] * b.coeff.member(j);
						multiply(tmp2, MT_vec[l++], a.coeff.member(j));
						add(tmp1, tmp1, tmp2);
					}
					// tmp1 *= 2;
					tmp1.multiply_by_2();
					// tmp1 += MT_vec[l] * a.coeff.member(i);
					multiply(tmp2, MT_vec[l++], a.coeff.member(i));
					add(tmp1, tmp1, tmp2);
					// tmp1 *= a.coeff.member(i);
					multiply(tmp1, tmp1, a.coeff.member(i));

					// tmp[k] += tmp1;
					add(tmp[k], tmp[k], tmp1);
				}
			}
			swap(c.coeff, tmp);
		}
	}
	else {
		// Polynomarithmetik verwenden!
		bigint * tmp;

		polynomial< bigint > pa(a.coeff.get_data_address(),
					a.degree()-1);

		polynomial< bigint > f = (pa*pa)%(O->f);
		register lidia_size_t l = 0;
		for (; l <= f.degree(); l++)
			c.coeff[l] = f[l];
		for (; l < a.degree(); l++)
			c.coeff[l] = 0;
	}
	c.normalize();
	debug_handler_c("alg_number",
			"in function multiply(a_n &, const a_n &, const a_n &)", 0,
			std::cout << "result is " << c << std::endl);
}



void swap(alg_number & a, alg_number & b)
{
	swap (a.den, b.den);
	//     math_vector <bigint> help = a.coeff;
	//     a.coeff = b.coeff;
	//     b.coeff = help;

	swap(a.coeff, b.coeff);

	nf_base * O = a.O;
	a.O = b.O;
	b.O = O;
}



// random numbers
void alg_number::randomize(const bigint & b)
{
	for (lidia_size_t i = 0; i < degree(); i++)
		coeff[i] = LiDIA::randomize(b);
	den = LiDIA::randomize(b)+1;
}



// In-/Output:
std::ostream& operator << (std::ostream & s, const alg_number & a)
{
	s << a.coeff;
	if (!(a.is_zero() || a.den.is_one()))
		s << " / " << a.den;
	return s;
}



std::istream& operator >> (std::istream & s, alg_number & a)
{
	if (nf_base::current_base == nf_base::dummy_base) {
		lidia_error_handler ("alg_number",
				     "operator >>::No number_field or order!");
		return s;
	}
	if (a.O != nf_base::current_base) {
		a.O->dec_ref();
		a.O = nf_base::current_base;
		a.O->inc_ref();
	}
	s >> a.coeff;
	char c;
	do {
		s.get(c);
	} while (isspace(c) && c != '\n');
	if (c == '/') {
		s >> a.den;
		a.normalize();
	}
	else {
		a.den = 1;
		if (c != '\n' && c != '\r')
			s.putback(c);
	}
	return(s);
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
