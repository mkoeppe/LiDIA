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
//	Author	: Volker Mueller (VM)
//	Changes	: See CVS log
//
//==============================================================================================


#ifdef HAVE_CONFIG_H
# include	"config.h"
#endif
#include	"LiDIA/Fp_rational_function.h"
#include	<cctype>




#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



// ---------- CONSTRUCTORS / DESTRUCTORS ----------

Fp_rational_function::Fp_rational_function()
{
	num = new Fp_polynomial;
	den = new Fp_polynomial;
	set_modulus(3);
	den->assign_one();
	num->assign_zero();
}



Fp_rational_function::Fp_rational_function(const bigint & p)
{
	num = new Fp_polynomial;
	den = new Fp_polynomial;
	if (!is_prime(p) || p.is_even())
		lidia_error_handler("Fp_rational_function", "ct::no odd prime number");
	set_modulus(p);
	den->assign_one();
	num->assign_zero();
}



Fp_rational_function::Fp_rational_function(const Fp_polynomial & f)
{
	num = new Fp_polynomial;
	den = new Fp_polynomial;
	set_modulus(f.modulus());
	num->assign(f);
	den->assign_one();
}



Fp_rational_function::Fp_rational_function(const Fp_polynomial & f,
					   const Fp_polynomial & g)
{
	if (g.modulus() != f.modulus())
		lidia_error_handler("Fp_rational_function", "ct::different moduli");
	num = new Fp_polynomial;
	den = new Fp_polynomial;
	set_modulus(f.modulus());
	num->assign(f);
	den->assign(g);
}



Fp_rational_function::Fp_rational_function(const Fp_rational_function & f)
{
	num = new Fp_polynomial;
	den = new Fp_polynomial;
	set_modulus(f.modulus());
	num->assign(*f.num);
	den->assign(*f.den);
}



Fp_rational_function::~Fp_rational_function()
{
	if (num != NULL)
		delete num;

	if (den != NULL)
		delete den;
}



//******************************************************************
// some simple assignment functions   etc.
//******************************************************************

void
Fp_rational_function::kill()
{
	num->kill();
	den->kill();
}



void
Fp_rational_function::set_modulus(const bigint & p)
{
	num->set_modulus(p);
	den->set_modulus(p);
}



const bigint &
Fp_rational_function::modulus() const
{
	return num->modulus();
}



lidia_size_t
Fp_rational_function::degree_numerator() const
{
	return num->degree();
}



lidia_size_t
Fp_rational_function::degree_denominator() const
{
	if (num->is_zero())
		return -1;
	else
		return den->degree();
}



void
Fp_rational_function::get_coefficient_numerator(bigint & a, lidia_size_t i) const
{
	num->get_coefficient(a, i);
}



void
Fp_rational_function::get_coefficient_denominator(bigint & a, lidia_size_t i) const
{
	den->get_coefficient(a, i);
}



void
Fp_rational_function::set_coefficient_numerator(const bigint& a, lidia_size_t i)
{
	num->set_coefficient(a, i);
}



void
Fp_rational_function::set_coefficient_denominator(const bigint& a, lidia_size_t i)
{
	den->set_coefficient(a, i);
}



void
Fp_rational_function::set_coefficient_numerator(lidia_size_t i)
{
	num->set_coefficient(i);
}



void
Fp_rational_function::set_coefficient_denominator(lidia_size_t i)
{
	den->set_coefficient(i);
}



const bigint&
Fp_rational_function::lead_coefficient_numerator() const
{
	return num->lead_coeff();
}



const bigint&
Fp_rational_function::lead_coefficient_denominator() const
{
	if (num->is_zero()) {
		lidia_error_handler("Fp_rational_function",
				    "lead_coeff_denominator::zero rational function");
		return num->modulus();
	}
	else
		return den->lead_coeff();
}



const bigint&
Fp_rational_function::const_term_numerator() const
{
	return num->const_term();
}



const bigint&
Fp_rational_function::const_term_denominator() const
{
	if (num->is_zero()) {
		lidia_error_handler("Fp_rational_function",
				    "const_term_denominator::zero rational function");
		return num->modulus(); }
	else
		return den->const_term();
}



//********** Assignments **************************************

Fp_rational_function &
Fp_rational_function::operator = (const Fp_rational_function & f)
{
	if (this != &f)
		this->assign(*f.num, *f.den);
	return *this;
}



Fp_rational_function &
Fp_rational_function::operator = (const Fp_polynomial & f)
{
	num->assign(f);
	den->assign_one();
	return *this;
}



void
Fp_rational_function::assign(const Fp_rational_function & f)
{
	num->assign(*f.num);
	den->assign(*f.den);
}



void
Fp_rational_function::assign(const Fp_polynomial & f, const Fp_polynomial & g)
{
	if (f.modulus() != g.modulus())
		lidia_error_handler("Fp_rational_function", "assign::different moduli");

	num->assign(f);
	den->assign(g);
}



void
Fp_rational_function::assign(const Fp_polynomial & f)
{
	num->assign(f);
	den->assign_one();
}



void
Fp_rational_function::assign_numerator(const Fp_polynomial &a)
{
	if (a.modulus() != num->modulus())
		lidia_error_handler("Fp_rational_function",
				    "assign_numerator::different moduli");
	num->assign(a);
}



void
Fp_rational_function::assign_denominator (const Fp_polynomial &a)
{
	if (a.modulus() != num->modulus())
		lidia_error_handler("Fp_rational_function",
				    "assign_denominator::different moduli");
	den->assign(a);
}



void
Fp_rational_function::assign_zero()
{
	num->assign_zero();
	den->assign_one();
}



void
Fp_rational_function::assign_one()
{
	num->assign_one();
	den->assign_one();
}



void
Fp_rational_function::assign_x()
{
	num->assign_x();
	den->assign_one();
}



void
Fp_rational_function::randomize(lidia_size_t deg_num, lidia_size_t deg_denom)
{
	num->randomize(deg_num);
	den->randomize(deg_denom);
}



Fp_polynomial &
Fp_rational_function::numerator()
{
	return static_cast<Fp_polynomial &>(*num);
}



Fp_polynomial & Fp_rational_function::denominator()
{
	return static_cast<Fp_polynomial &>(*den);
}



const Fp_polynomial &
Fp_rational_function::numerator() const
{
	return static_cast<const Fp_polynomial &>(*num);
}



const  Fp_polynomial & Fp_rational_function::denominator() const
{
	return static_cast<const Fp_polynomial &>(*den);
}



//************ Comparisons *************************************

bool
operator == (const Fp_rational_function & a, const Fp_rational_function & b)
{
	Fp_polynomial h1, h2;

	multiply(h1, *a.num , *b.den);
	multiply (h2 , *a.den , *b.num);
	return (h1 == h2);
}



bool
operator != (const Fp_rational_function & a,
		  const Fp_rational_function & b)
{
	return (!(a == b));
}



bool
Fp_rational_function::is_zero() const
{
	return num->is_zero();
}



bool
Fp_rational_function::is_one() const
{
	return ((*num) == (*den));
}



bool
equal_mod (const Fp_rational_function & a, const Fp_rational_function & b,
	   const Fp_polynomial & f)
{
	Fp_polynomial h1, h2;

	if ((b.is_zero() && !a.is_zero()) || (a.is_zero() && !b.is_zero()))
		return false;

	if (a.is_zero() && b.is_zero())
		return true;

	multiply_mod(h1, *a.num, *b.den, f);
	multiply_mod(h2, *a.den, *b.num, f);
	return (h1 == h2);
}



bool
equal (const Fp_rational_function & a, const Fp_rational_function & b,
       const Fp_poly_modulus & f)
{
	Fp_polynomial h1, h2;

	if ((b.is_zero() && !a.is_zero()) || (a.is_zero() && !b.is_zero()))
		return false;

	if (a.is_zero() && b.is_zero())
		return true;

	multiply(h1, *a.num, *b.den, f);
	multiply(h2, *a.den, *b.num, f);
	return (h1 == h2);
}



bigint
Fp_rational_function::operator() (const bigint & a) const
{
	bigint p(modulus());
	bigint res, h, aa(a % p);

	res = (*den)(aa);
	h = xgcd_left(res, res, p);

	if (!h.is_one()) {
		lidia_error_handler("Fp_rational_function", "operator(const bigint &a)::denominator not invertible");
		return 0;
	}

	multiply(res, res, (*num)(aa));
	remainder(res, res, p);
	return res;
}



//******** procedural versions for arithmetic *****************

void
add (Fp_rational_function & x, const Fp_rational_function & a,
          const Fp_rational_function & b)
{
	if (&a == &b) {
		add(*x.num, *a.num, *a.num);
		x.den->assign(*a.den);
	}
	else
		if ((*a.den) == (*b.den)) {
			x.den->assign(*a.den);
			add(*x.num, *a.num, *b.num);
		}
		else {
			Fp_polynomial h1, h2;

			multiply(h1, *a.num, *b.den);
			multiply(h2, *a.den, *b.num);

			multiply(*x.den, *a.den, *b.den);
			add(*x.num, h1, h2);
		}
}



void
add_mod (Fp_rational_function & x, const Fp_rational_function & a,
              const Fp_rational_function & b, const Fp_polynomial & f)
{
	if (&a == &b) {
		add(*x.num, *a.num, *a.num);
		x.den->assign(*a.den);
	}
	else
		if ((*a.den) == (*b.den)) {
			x.den->assign(*a.den);
			add(*x.num, *a.num, *b.num);
		}
		else {
			Fp_polynomial h1, h2;
			multiply_mod(h1, *a.num, *b.den, f);
			multiply_mod(h2, *a.den, *b.num, f);
			multiply_mod(*x.den, *a.den, *b.den, f);
			add(*x.num, h1, h2);
		}
}



void
add (Fp_rational_function & x, const Fp_rational_function & a,
          const Fp_rational_function & b, const Fp_poly_modulus & F)
{
	if (&a == &b) {
		add(*x.num, *a.num, *a.num);
		x.den->assign(*a.den);
	}
	else
		if ((*a.den) == (*b.den)) {
			x.den->assign(*a.den);
			add(*x.num, *a.num, *b.num);
		}
		else {
			Fp_polynomial h1, h2;
			multiply(h1, *a.num, *b.den, F);
			multiply(h2, *a.den, *b.num, F);
			multiply(*x.den, *a.den, *b.den, F);
			add(*x.num, h1, h2);
		}
}



void
subtract (Fp_rational_function & x, const Fp_rational_function & a,
               const Fp_rational_function & b)
{
	if (&a == &b) {
		x.num->assign_zero();
		x.den->assign_one();
	}
	else
		if ((*a.den) == (*b.den)) {
			x.den->assign(*a.den);
			subtract(*x.num, *a.num, *b.num);
		}
		else {
			Fp_polynomial h1, h2;

			multiply(h1, *a.num, *b.den);
			multiply(h2, *a.den, *b.num);
			multiply(*x.den, *a.den, *b.den);
			subtract(*x.num, h1, h2);
		}
}



void
subtract_mod (Fp_rational_function & x, const Fp_rational_function & a,
                   const Fp_rational_function & b, const Fp_polynomial & f)
{
	if (&a == &b) {
		x.num->assign_zero();
		x.den->assign_one();
	}
	else
		if ((*a.den) == (*b.den)) {
			x.den->assign(*a.den);
			subtract(*x.num, *a.num, *b.num);
		}
		else {
			Fp_polynomial h1, h2;
			multiply_mod(h1, *a.num, *b.den, f);
			multiply_mod(h2, *a.den, *b.num, f);
			multiply_mod(*x.den, *a.den, *b.den, f);
			subtract(*x.num, h1, h2);
		}
}



void
subtract (Fp_rational_function & x, const Fp_rational_function & a,
               const Fp_rational_function & b, const Fp_poly_modulus & F)
{
	if (&a == &b) {
		x.num->assign_zero();
		x.den->assign_one();
	}
	else
		if ((*a.den) == (*b.den)) {
			x.den->assign(*a.den);
			subtract(*x.num, *a.num, *b.num);
		}
		else {
			Fp_polynomial h1, h2;
			multiply(h1, *a.num, *b.den, F);
			multiply(h2, *a.den, *b.num, F);
			multiply(*x.den, *a.den, *b.den, F);
			subtract(*x.num, h1, h2);
		}
}



void
divide(Fp_rational_function & q, const Fp_rational_function & a,
            const Fp_rational_function & b)
{
	if (&a == &b) {
		q.num->assign_one();
		q.den->assign_one();
	}
	else
		if (&q == & b) {
			Fp_rational_function r(*b.den, *b.num);
			multiply(q, a, r);
		}
		else {
			multiply(*q.num, *a.num, *b.den);
			multiply(*q.den, *a.den, *b.num);
		}
}



void
divide(Fp_rational_function & q, const Fp_polynomial & a,
            const Fp_rational_function & b)
{
	if (&q == &b) {
		Fp_rational_function r(*b.den, *b.num);
		multiply(q, r, a);
	}
	else {
		q.den->assign(*b.num); multiply(*q.num, a, *b.den);
	}
}



void
invert(Fp_rational_function & c, const Fp_rational_function & a)
{
	if (&a == &c) {
		Fp_rational_function r(*a.den, *a.num);
		c.assign(r);
	}
	else {
		c.num->assign(*a.den);
		c.den->assign(*a.num);
	}
}



void
Fp_rational_function::invert()
{
	Fp_rational_function c(*this);

	num->assign(c.denominator());
	den->assign(c.numerator());
}



//************ Miscellaneous functions ************************

void
derivative(Fp_rational_function &f, const Fp_rational_function & g)
{
	square(*f.den, *g.den);
	subtract(*f.num, derivative(*g.num) * (*g.den),
		 derivative(*g.den) * (*g.num));
}



//**************************** IO *********************************

void
Fp_rational_function::read(std::istream & s)
{
	bigint p;
	char c, cm, co, cd;

	s >> std::ws >> c;

	if (c != '[') {
		s.putback(c);
		pretty_read(s);
		return;
	}

	s.putback(c);
	base_vector< bigint > help_coeff1, help_coeff2;

	s >> help_coeff1;

	s >> cm;

	if (cm != '/' && cm != 'm') {
		lidia_error_handler("Fp_rational_function", "read(std::istream&)::'/' expected");
		return;
	}

	if (cm == '/')  // denominator not one
	{
		s >> help_coeff2;
		s >> cm >> co >> cd;

		if (cm != 'm' || co != 'o' || cd != 'd') {
			lidia_error_handler("Fp_rational_function", "read(std::istream&)::'mod' expected");
			return;
		}

		s >> p;
		set_modulus(p);

		help_coeff1.reverse(); help_coeff2.reverse();
		num->assign(help_coeff1, p);
		den->assign(help_coeff2, p);
		if (den->is_zero())
			lidia_error_handler("Fp_rational_function", "read()::denominator is zero");
		if (num->is_zero())
			den->assign_one();
	}
	else {
		s >> co >> cd;
		if (cm != 'm' || co != 'o' || cd != 'd') {
			lidia_error_handler("Fp_rational_function", "read(std::istream&)::'mod' expected");
			return;
		}
		s >> p;
		set_modulus(p);
		help_coeff1.reverse();
		num->assign(help_coeff1, p);
		den->assign_one();
	}
}



void
Fp_rational_function::print(std::ostream & s) const
{
	s << "[";
	if (is_zero()) {
		s << "] mod " << modulus();
		return;
	}
	lidia_size_t i;
	for (i = num->degree(); i >= 0; i--) {
		s << (*num)[i];
		if (i != 0)
			s << " ";
	}
	s << "] ";
	if (den->is_one()) {
		s << "mod " << modulus();
		return;
	}
	else s << "/ [";
	for (i = den->degree(); i >= 0; i--) {
		s << (*den)[i];
		if (i != 0)
			s << " ";
	}
	s << "] mod " << modulus();
}



void
Fp_rational_function::pretty_print(std::ostream &os) const
{
	if (is_zero()) {
		os << "0 mod " << modulus();
		return;
	}

	lidia_size_t j, deg = num->degree();
	bigint coeff;

	for (j = deg; j >= 0; j--) {
		num->get_coefficient(coeff, j);
		if (!coeff.is_zero()) {
			if (coeff > 0 && j < deg)
				os << " + ";
			if (j == 0 || !coeff.is_one()) {
				os << coeff;
				if (j > 0)
					os << "*";
			}
			if (j > 0)
				os << "x";
			if (j > 1)
				os << "^" << j;
		}
	}
	if (den->is_one())
		os << " mod " << modulus();
	else {
		os << " / ";
		deg = den->degree();
		for (j = deg; j >= 0; j--) {
			den->get_coefficient(coeff, j);
			if (!coeff.is_zero()) {
				if (j < deg)
					os << " + ";
				if (j == 0 || !coeff.is_one()) {
					os << coeff;
					if (j > 0)
						os << "*";
				}
				if (j > 0)
					os << "x";
				if (j > 1)
					os << "^" << j;
			}
		}
		os << " mod " << modulus();
	}
}



void
Fp_rational_function::pretty_read(std::istream & s)
{
	// This function reads a univariate rational function in any variable.
	// input format : a_n*x^n+ ... + a_1*x + a_0 / b_m*x^m + ... b_0 mod p
	// Monomials need not be sorted, and powers of x may even appear
	// repeated, '*' may be omitted and coefficients may follow the variable:
	//        -4 + 8x^5 + 2 - x^2 3 + x^5 + x^5*17 mod 5
	// Note however, that the routine will work faster, if the leading monomial
	// of the numerator, denominator is read first.
	// All coefficients will be reduced mod p.

	lidia_size_t sz;
	char c;
	base_vector< bigint > help_coeff(8, EXPAND);
	base_vector< bigint > help_coeff2(8, EXPAND);
	bool mod_read = false, jump_to_den = false;

	char variable = 0;
	bigint coeff_tmp(1);
	bigint tmp;

	// Read a monomial, i.e. "x^k" or "- x^k"
	// or "a*x^k" or "a x^k" or "x^k*a" or "x^k a"

	do {
		c = s.get();
	} while (isspace(c) && c != '\n');

	while (c != '\n' && c != EOF && s.good() && !jump_to_den) {
		sz = 0; // Assume we read coeffizient of x^0;
		if (c == '+') {
			coeff_tmp.assign_one();
			do {
				c = s.get();
			} while (isspace(c) && c != '\n');
		}
		if (c == '-') {
			coeff_tmp = -1;
			do {
				c = s.get();
			} while (isspace(c) && c != '\n');
		}
		if (c >= '0' && c <= '9' || c == '(') {
			s.putback(c);
			s >> tmp;
			multiply(coeff_tmp, coeff_tmp, tmp);
			do {
				c = s.get();
			} while (isspace(c) && c != '\n');
			if (c == '*') {
				do {
					c = s.get();
				} while (isspace(c) && c != '\n');
				if (c == 'm') {
					c = s.get();
					if (c == 'o') {
						lidia_error_handler("Fp_rational_function",
								    "pretty_read(...)::wrong input format");
						return;
					}
					else s.putback(c);
				}
			}
		}
		if (c == 'm') {
			c = s.get();
			if (c == 'o') {
				c = s.get();
				if (c != 'd') {
					lidia_error_handler("Fp_rational_function",
							    "pretty_read (...):: 'mod' expected");
					return;
				}
				c = '\n';
				mod_read = true;
			}
			else s.putback(c);
		}
		if ((c >= 'A' && c <= 'Z') || (c >= 'a' && c <= 'z')) {
			if (variable == 0)
				variable = c;
			else if (variable != c) {
				lidia_error_handler_c("Fp_rational_function", "pretty_read (...)::"
						      "The given string is not recognized to be"
						      " a univariate polynomial",
						      std::cout << "Variable name seemed to be " << variable;
						      std::cout << " and now you used " << c << "." << std::endl);
				return;
			}
			do {
				c = s.get();
			} while (isspace(c) && c != '\n');

			if (c != '^') sz = 1;
			else {
				s >> sz;
				do {
					c = s.get();
				} while (isspace(c) && c != '\n');
			}
			if (c == '*') {
				s >> tmp;
				multiply(coeff_tmp, coeff_tmp, tmp);
				do {
					c = s.get();
				} while (isspace(c) && c != '\n');
			}

			if (c >= '0' && c <= '9' || c == '(') {
				s.putback(c);
				s >> tmp;
				multiply(coeff_tmp, coeff_tmp, tmp);
				do {
					c = s.get();
				} while (isspace(c) && c != '\n');
			}
		}
		if (c == 'm') {
			c = s.get();
			if (c != 'o') {
				lidia_error_handler("Fp_rational_function",
						    "pretty_read (...):: 'mod' expected");
				return;
			}
			c = s.get();
			if (c != 'd') {
				lidia_error_handler("Fp_rational_function",
						    "pretty_read (...):: 'mod' expected");
				return;
			}

			c = '\n';
			mod_read = true;
		}

		if (c == '/') {
			c = '\n';
			jump_to_den = true;
		}
		if (c != '+' && c != '-' && c != '\n') {
			// No next monomial, so assume end of input is reached
			s.putback(c);
			c = '\n'; // set c to end--marker
		}
		add(help_coeff[sz], help_coeff[sz], coeff_tmp);
	}

	bigint p;

	if (mod_read)  // denominator is one
	{
		s >> p;
		set_modulus(p);

		num->assign(help_coeff, p);
		den->assign_one();
		return;
	}

	do {
		c = s.get();
	} while (isspace(c) && c != '\n');

	while (c != '\n' && c != EOF && s.good()) {
		sz = 0; // Assume we read coeffizient of x^0;
		coeff_tmp.assign_one();
		if (c == '+')
			do {
				c = s.get();
			} while (isspace(c) && c != '\n');
		if (c == '-') {
			coeff_tmp = -1;
			do {
				c = s.get();
			} while (isspace(c) && c != '\n');
		}
		if (c >= '0' && c <= '9' || c == '(') {
			s.putback(c);
			s >> tmp;
			multiply(coeff_tmp, coeff_tmp, tmp);
			do {
				c = s.get();
			} while (isspace(c) && c != '\n');
			if (c == '*') {
				do {
					c = s.get();
				} while (isspace(c) && c != '\n');
				if (c == 'm') {
					c = s.get();
					if (c == 'o') {
						lidia_error_handler("Fp_rational_function",
								    "pretty_read(...)::wrong input format");
						return;
					}
					else s.putback(c);
				}
			}
		}
		if (c == 'm') {
			c = s.get();
			if (c == 'o') {
				c = s.get();
				if (c != 'd') {
					lidia_error_handler("Fp_rational_function",
							    "pretty_read (...):: 'mod' expected");
					return;
				}
				c = '\n';
			}
			else s.putback(c);
		}
		if ((c >= 'A' && c <= 'Z') || (c >= 'a' && c <= 'z')) {
			if (variable == 0)
				variable = c;
			else if (variable != c) {
				lidia_error_handler_c("Fp_rational_function", "pretty_read (...)::"
						      "The given string is not recognized to be"
						      " a univariate polynomial",
						      std::cout << "Variable name seemed to be " << variable;
						      std::cout << " and now you used " << c << "." << std::endl);
				return;
			}
			do {
				c = s.get();
			} while (isspace(c) && c != '\n');

			if (c != '^') sz = 1;
			else {
				s >> sz;
				do {
					c = s.get();
				} while (isspace(c) && c != '\n');
			}
			if (c == '*') {
				s >> tmp;
				multiply(coeff_tmp, coeff_tmp, tmp);
				do {
					c = s.get();
				} while (isspace(c) && c != '\n');
			}

			if (c >= '0' && c <= '9' || c == '(') {
				s.putback(c);
				s >> tmp;
				multiply(coeff_tmp, coeff_tmp, tmp);
				do {
					c = s.get();
				} while (isspace(c) && c != '\n');
			}
		}
		if (c == 'm') {
			c = s.get();
			if (c != 'o') {
				lidia_error_handler("Fp_rational_function",
						    "pretty_read (...):: 'mod' expected");
				return;
			}
			c = s.get();
			if (c != 'd') {
				lidia_error_handler("Fp_rational_function",
						    "pretty_read (...):: 'mod' expected");
				return;
			}

			c = '\n';
		}
		if (c != '+' && c != '-' && c != '\n') {
			// No next monomial, so assume end of input is reached
			s.putback(c);
			c = '\n'; // set c to end--marker
		}
		add(help_coeff2[sz], help_coeff2[sz], coeff_tmp);
	}

	s >> p;
	set_modulus(p);

	num->assign(help_coeff, p);
	den->assign(help_coeff2, p);
	if (den->is_zero())
		lidia_error_handler("Fp_rational_function", "pretty_read::denominator is zero");
	if (num->is_zero())
		den->assign_one();
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
