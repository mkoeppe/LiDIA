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
//	Author	: Volker Mueller
//	Changes	: See CVS log
//
//==============================================================================================


#ifdef HAVE_CONFIG_H
# include	"config.h"
#endif
#include	"LiDIA/gf2n_rational_function.h"



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



// ---------- CONSTRUCTORS / DESTRUCTORS ----------

gf2n_rational_function::gf2n_rational_function()
{
	num = new gf2n_polynomial;
	den = new gf2n_polynomial;
	den->assign_one();
	num->assign_zero();
}



gf2n_rational_function::gf2n_rational_function(const gf2n_polynomial & f)
{
	num = new gf2n_polynomial;
	den = new gf2n_polynomial;
	num->assign(f);
	den->assign_one();
}



gf2n_rational_function::gf2n_rational_function(const gf2n_polynomial & f,
					       const gf2n_polynomial & g)
{
	num = new gf2n_polynomial;
	den = new gf2n_polynomial;
	num->assign(f);
	den->assign(g);
}



gf2n_rational_function::
gf2n_rational_function(const gf2n_rational_function & f)
{
	num = new gf2n_polynomial;
	den = new gf2n_polynomial;
	num->assign(*f.num);
	den->assign(*f.den);
}



gf2n_rational_function::~gf2n_rational_function()
{
	if (num != NULL)
		delete num;

	if (den != NULL)
		delete den;
}



//******************************************************************
//  some simple assignment functions   etc.
//******************************************************************

lidia_size_t gf2n_rational_function::degree_numerator() const
{
	return num->degree();
}



lidia_size_t gf2n_rational_function::degree_denominator() const
{
	if (num->is_zero())
		return -1;
	else
		return den->degree();
}



void gf2n_rational_function::get_coefficient_numerator(gf2n & a, lidia_size_t i) const
{
	num->get_coefficient(a, i);
}



void gf2n_rational_function::get_coefficient_denominator(gf2n & a, lidia_size_t i) const
{
	den->get_coefficient(a, i);
}



void gf2n_rational_function::set_coefficient_numerator(const gf2n& a, lidia_size_t i)
{
	num->set_coefficient(a, i);
}



void gf2n_rational_function::set_coefficient_denominator(const gf2n& a, lidia_size_t i)
{
	den->set_coefficient(a, i);
}



void gf2n_rational_function::set_coefficient_numerator(lidia_size_t i)
{
	num->set_coefficient(i);
}



void gf2n_rational_function::set_coefficient_denominator(lidia_size_t i)
{
	den->set_coefficient(i);
}



gf2n gf2n_rational_function::lead_coefficient_numerator() const
{
	return num->lead_coeff();
}



gf2n  gf2n_rational_function::lead_coefficient_denominator() const
{
	if (num->is_zero()) {
		lidia_error_handler("gf2n_rational_function",
				    "lead_coeff_denominator::zero rational function");
		return gf2n();
	}
	else
		return den->lead_coeff();
}



gf2n gf2n_rational_function::const_term_numerator() const
{
	return num->const_term();
}



gf2n gf2n_rational_function::const_term_denominator() const
{
	if (num->is_zero()) {
		lidia_error_handler("gf2n_rational_function",
				    "const_term_denominator::zero rational function");
		return gf2n();
	}
	else
		return den->const_term();
}



//********** Assignments **************************************

gf2n_rational_function &
gf2n_rational_function::operator = (const gf2n_rational_function & f)
{
	if (this != &f)
		this->assign(*f.num, *f.den);
	return *this;
}



gf2n_rational_function &
gf2n_rational_function::operator = (const gf2n_polynomial & f)
{
	num->assign(f);
	den->assign_one();
	return *this;
}



void gf2n_rational_function::assign(const gf2n_rational_function & f)
{
	num->assign(*f.num);
	den->assign(*f.den);
}



void gf2n_rational_function::assign(const gf2n_polynomial & f,
				    const gf2n_polynomial & g)
{
	num->assign(f);
	den->assign(g);
}



void gf2n_rational_function::assign(const gf2n_polynomial & f)
{
	num->assign(f);
	den->assign_one();
}



void gf2n_rational_function::assign_numerator(const gf2n_polynomial &a)
{
	num->assign(a);
}



void gf2n_rational_function::assign_denominator (const gf2n_polynomial &a)
{
	den->assign(a);
}



void gf2n_rational_function::assign_zero()
{
	num->assign_zero();
	den->assign_one();
}



void gf2n_rational_function::assign_one()
{
	num->assign_one();
	den->assign_one();
}



void gf2n_rational_function::assign_x()
{
	num->assign_x();
	den->assign_one();
}



void gf2n_rational_function::randomize(lidia_size_t deg_num,
				       lidia_size_t deg_denom)
{
	num->randomize(deg_num);
	den->randomize(deg_denom);
}



gf2n_polynomial & gf2n_rational_function::numerator()
{
	return static_cast<gf2n_polynomial &>(*num);
}



gf2n_polynomial & gf2n_rational_function::denominator()
{
	return static_cast<gf2n_polynomial &>(*den);
}



const  gf2n_polynomial & gf2n_rational_function::numerator() const
{
	return static_cast<const gf2n_polynomial &>(*num);
}



const  gf2n_polynomial & gf2n_rational_function::denominator() const
{
	return static_cast<gf2n_polynomial &>(*den);
}



//************ Comparisons *************************************

bool operator == (const gf2n_rational_function & a,
		  const gf2n_rational_function & b)
{
	gf2n_polynomial h1, h2;

	multiply(h1, *a.num , *b.den);
	multiply (h2 , *a.den , *b.num);
	return (h1 == h2);
}



bool operator != (const gf2n_rational_function & a,
		  const gf2n_rational_function & b)
{
	return (!(a == b));
}



bool gf2n_rational_function::is_zero() const
{
	return num->is_zero();
}



bool gf2n_rational_function::is_one() const
{
	return ((*num) == (*den));
}



bool equal_mod (const gf2n_rational_function & a,
		const gf2n_rational_function & b,
                const gf2n_polynomial & f)
{
	gf2n_polynomial h1, h2;

	if ((a.is_zero() && !b.is_zero()) || (!a.is_zero() && b.is_zero()))
		return false;

	multiply_mod(h1, *a.num, *b.den, f);
	multiply_mod(h2, *a.den, *b.num, f);
	return (h1 == h2);
}



bool equal (const gf2n_rational_function & a, const gf2n_rational_function & b,
            gf2n_poly_modulus & f)
{
	gf2n_polynomial h1, h2;

	if ((a.is_zero() && !b.is_zero()) || (!a.is_zero() && b.is_zero()))
		return false;

	multiply(h1, *a.num, *b.den, f);
	multiply(h2, *a.den, *b.num, f);
	return (h1 == h2);
}



gf2n gf2n_rational_function::operator() (const gf2n & a) const
{
	gf2n res;

	res = (*den)(a);
	res.invert();
	multiply(res, res, (*num)(a));
	return res;
}



//******** procedural versions for arithmetic *****************

void add (gf2n_rational_function & x, const gf2n_rational_function & a,
          const gf2n_rational_function & b)
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
			gf2n_polynomial h1, h2;

			multiply(h1, *a.num, *b.den);
			multiply(h2, *a.den, *b.num);

			multiply(*x.den, *a.den, *b.den);
			add(*x.num, h1, h2);
		}
}



void add_mod (gf2n_rational_function & x, const gf2n_rational_function & a,
              const gf2n_rational_function & b, const gf2n_polynomial & f)
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
			gf2n_polynomial h1, h2;
			multiply_mod(h1, *a.num, *b.den, f);
			multiply_mod(h2, *a.den, *b.num, f);
			multiply_mod(*x.den, *a.den, *b.den, f);
			add(*x.num, h1, h2);
		}
}



void add (gf2n_rational_function & x, const gf2n_rational_function & a,
          const gf2n_rational_function & b, gf2n_poly_modulus & F)
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
			gf2n_polynomial h1, h2;
			multiply(h1, *a.num, *b.den, F);
			multiply(h2, *a.den, *b.num, F);
			multiply(*x.den, *a.den, *b.den, F);
			add(*x.num, h1, h2);
		}
}



void divide(gf2n_rational_function & q, const gf2n_rational_function & a,
            const gf2n_rational_function & b)
{
	if (&a == &b) {
		q.num->assign_one();
		q.den->assign_one();
	}
	else
		if (&q == & b) {
			gf2n_rational_function r(*b.den, *b.num);
			multiply(q, a, r);
		}
		else {
			multiply(*q.num, *a.num, *b.den);
			multiply(*q.den, *a.den, *b.num);
		}
}



void divide(gf2n_rational_function & q, const gf2n_polynomial & a,
            const gf2n_rational_function & b)
{
	if (&q == &b) {
		gf2n_rational_function r(*b.den, *b.num);
		multiply(q, r, a);
	}
	else {
		q.den->assign(*b.num); multiply(*q.num, a, *b.den);
	}
}



void invert(gf2n_rational_function & c, const gf2n_rational_function & a)
{
	if (&a == &c) {
		gf2n_rational_function r(*a.den, *a.num);
		c.assign(r);
	}
	else {
		c.num->assign(*a.den);
		c.den->assign(*a.num);
	}
}



void gf2n_rational_function::invert()
{
	gf2n_rational_function c(*this);

	num->assign(c.denominator());
	den->assign(c.numerator());
}



//**************************** IO *********************************

void gf2n_rational_function::read(std::istream & s)
{
	char c, cm;

	s >> std::ws >> c;

	if (c != '[') {
		s.putback(c);
		//pretty_read(s);
		lidia_error_handler("gf2n_rational_function", "pretty_read::sorry, not"
				    " yet implemented");
		return;
	}

	s.putback(c);
	base_vector< gf2n > help_coeff1, help_coeff2;

	s >> help_coeff1;

	s >> cm;

	if (cm == '/') {
		// denominator not one
		s >> help_coeff2;

		help_coeff1.reverse(); help_coeff2.reverse();
		num->assign(help_coeff1);
		den->assign(help_coeff2);
		if (den->is_zero())
			lidia_error_handler("gf2n_rational_function",
					    "read()::denominator is zero");
		if (num->is_zero())
			den->assign_one();
	}
	else {
		help_coeff1.reverse();
		num->assign(help_coeff1);
		den->assign_one();
	}
}



void gf2n_rational_function::print(std::ostream & s) const
{
	s << "[";
	if (is_zero()) {
		s << "]";
		return;
	}
	lidia_size_t i;
	for (i = num->degree(); i >= 0; i--) {
		s << (*num)[i];
		if (i != 0)
			s << " ";
	}
	s << "] ";
	if (den->is_one())
		return;
	else
		s << "/ [";

	for (i = den->degree(); i >= 0; i--) {
		s << (*den)[i];
		if (i != 0)
			s << " ";
	}
	s << "]";
}



void gf2n_rational_function::pretty_print(std::ostream &os) const
{
	if (is_zero()) {
		os << "0"; return;
	}

	lidia_size_t j, deg = num->degree();
	gf2n coeff;

	std::cout << "(";

	for (j = deg; j >= 0; j--) {
		num->get_coefficient(coeff, j);
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
	std::cout << ")";
	if (!den->is_one()) {
		os << " / (";
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
		std::cout << ")";
	}
}



#if 0
void gf2n_rational_function::pretty_read(std::istream & s)
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
	base_vector< gf2n > help_coeff(8, EXPAND);
	base_vector< gf2n > help_coeff2(8, EXPAND);
	bool mod_read = false, jump_to_den = false;

	char variable = 0;
	gf2n coeff_tmp(1);
	gf2n tmp;

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
						lidia_error_handler("gf2n_rational_function",
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
					lidia_error_handler("gf2n_rational_function",
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
				lidia_error_handler_c("gf2n_rational_function", "pretty_read (...)::"
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
				lidia_error_handler("gf2n_rational_function",
						    "pretty_read (...):: 'mod' expected");
				return;
			}
			c = s.get();
			if (c != 'd') {
				lidia_error_handler("gf2n_rational_function",
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

	gf2n p;

	if (mod_read) {
		// denominator is one
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
						lidia_error_handler("gf2n_rational_function",
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
					lidia_error_handler("gf2n_rational_function",
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
				lidia_error_handler_c("gf2n_rational_function", "pretty_read (...)::"
						      "The given string is not recognized to be"
						      " a univariate polynomial",
						      std::cout << "Variable name seemed to be " << variable;
						      std::cout << " and now you used " << c << "." << std::endl);
				return;
			}
			do {
				c = s.get();
			} while (isspace(c) && c != '\n');

			if (c != '^')
				sz = 1;
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
				lidia_error_handler("gf2n_rational_function",
						    "pretty_read (...):: 'mod' expected");
				return;
			}
			c = s.get();
			if (c != 'd') {
				lidia_error_handler("gf2n_rational_function",
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
		lidia_error_handler("gf2n_rational_function", "pretty_read::denominator is zero");
	if (num->is_zero())
		den->assign_one();
}
#endif



void compose(gf2n_rational_function & f,
	     const gf2n_rational_function & g,
	     const gf2n_polynomial & p)
{
	if (p.is_zero())
		f.assign_zero();
	else {
		if (g.is_zero())
			f.assign(p.const_term());
		else {
			// We compute a powertable g^0, ... , g^(deg(p))
			// (squaring is cheaper than multiplication)

			gf2n_rational_function* gpowers;

			gpowers = new gf2n_rational_function [p.degree()];
			gpowers[0].assign_one();

			if (p.degree() >= 1)
				gpowers[1].assign(g);
			if (p.degree() >= 2) {
				gpowers[2].assign(g);
				square(gpowers[2], gpowers[2]);
			}

			int i, j;

			i = 2; j = 3; // compute table with squaring trick
			while (j <= p.degree()) {
				while (i <= p.degree()) {
					gpowers[i] = gpowers[i/2];
					square(gpowers[i], gpowers[i]);
					i <<= 1;
				}
				i = j; j += 2;
				multiply(gpowers[i], gpowers[i-1], gpowers[1]);
				i <<= 1;
				std::cout << "Squaringtrick: " << j << std::endl;
			}

			f.assign_zero();
			add(f, f, p.const_term());

			std::cout << "f = " << f << std::endl;

			for (i = 1; i <= p.degree(); i++) {


				multiply(gpowers[i], gpowers[i], p[i]);
				add(f, f, gpowers[i]);

				std::cout << "f = " << f << std::endl;


			}
			delete [] gpowers;

		}
	}
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
