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
//	Author	: Victor Shoup (VS), Thomas Pfahler (TPf)
//	Changes	: See CVS log
//
//==============================================================================================


#ifdef HAVE_CONFIG_H
# include	"config.h"
#endif
#include	"LiDIA/Fp_polynomial.h"
#include	<cctype>



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



residue_class_list< bigint > *Fp_polynomial::L = 0;

bigint Fp_polynomial::ZERO = 0;




//************************************************************
//                  class Fp_pol_ref
//  tries to differentiate between the use of
//  Fp_polynomial::operator[] as an l-value and
//  as an r-value
//
//  see: Coplien, Advanced C++ programming styles and idioms
//
//  const bigint & Fp_polynomial::operator[] (lidia_size_t i) const
//      always operates on R-values and therefore returns a const reference to
//      the coefficient of X^i
//  Fp_pol_ref Fp_polynomial::operator[] (lidia_size_t i)
//      always returns an instance of type Fp_pol_ref, which itself consists
//      only of a reference to a polynomial (Fp_polynomial &) and the number i,
//      which are both initialized by the only constructor.
//      In case the expression is used where an R-value of type bigint
//      is expected, a direct type conversion using
//          Fp_pol_ref::operator bigint()
//      is executed. But if we use the construction in the following contexts
//          a[i] = c;
//          a[i].assign(c);
//          a[i].assign_one();
//          a[i].assign_zero();
//      (a[i] is an L-value expression !), then the corresponding functions
//      within the class Fp_pol_ref are called, which in turn call the
//      appropriate functions for the coefficients.
//      ADVANTAGE   : we can always control the values the user assignes to
//                    coefficients
//                    this would not be possible if we just used
//                    bigint & Fp_polynomial::operator[] (lidia_size_t)
//      DISADVANTAGE: in order to make this technique completely transparent,
//                    we would have to implement all the functions which modify
//                    bigints again for the class Fp_pol_ref; we only
//                    implemented the most simple ones.
//                    It is FORBIDDEN to use Fp_polynomial::operator[]
//                    when a non-const reference to a bigint is expected,
//                    e.g. in expressions such as
//                        square(a[i], a[i]);
//                        multiply(a[i], b, c);
//
//************************************************************

Fp_pol_ref & Fp_pol_ref::operator = (const bigint & a)	// L-VALUE
{
	p.set_coefficient(a, ix);
	return *this;
}



void Fp_pol_ref::assign_one()
{
	p.set_coefficient(ix);
}



void Fp_pol_ref::assign_zero()
{
	p.set_coefficient(bigint(0), ix);
}



void Fp_pol_ref::assign(const bigint & a)
{
	p.set_coefficient(a, ix);
}



Fp_pol_ref::operator bigint()				// R-VALUE
{
	if (p.Mp == 0) {
		lidia_error_handler("Fp_polynomial", "operator[]::modulus is zero");
		return bigint();
	}

	if (ix < 0) {
		lidia_error_handler_c ("Fp_polynomial",
				       "operator[]::index out of range", std::cout << "index = " << ix << std::endl;);
		return bigint();
	}

	if (ix > p.degree())
		return bigint(0);

	return p.coeff[ix];
}



//***************************************************************
//
//	Constructors, Destructors, and Basic Functions
//
//***************************************************************

Fp_polynomial::Fp_polynomial(const Fp_polynomial & a) :
    	coeff(0),
    	c_length(0),
    	c_alloc(0),
	Mp(0)
{
	debug_handler_l("Fp_polynomial", "Fp_polynomial(const Fp_polynomial&)", 1);

	assign(a);
}



Fp_polynomial::Fp_polynomial(const polynomial< bigint > & f, const bigint & p) :
	coeff(0),
	c_length(0),
	c_alloc(0),
	Mp(0)
{
	debug_handler_l("Fp_polynomial", "Fp_polynomial(const polynomial< bigint > &, const bigint &)", 1);
	lidia_size_t i, deg = f.degree();
	set_modulus(p);
	set_degree(deg);
	for (i = 0; i <= deg; i++)
		Remainder(coeff[i], f[i], p);
	remove_leading_zeros();
}



Fp_polynomial::~Fp_polynomial()
{
	debug_handler_l("Fp_polynomial", "~Fp_polynomial()", 1);
	kill();
}



void Fp_polynomial::kill()
{
	debug_handler_l("Fp_polynomial", "kill (void)", 1);
	delete[] coeff;
	coeff = 0;
	c_alloc = 0;
	c_length = 0;
	if (Mp != 0)			//which means that Fp_polynomial::L
	{					//cannot be Nil
		Fp_polynomial::L->clear(Mp);
		Mp = 0;
	}
}



const bigint & Fp_polynomial::modulus() const
{
	if (Mp != 0)
		return Mp->get_mod();
	else
		return ZERO;
}



void Fp_polynomial::set_modulus(const Fp_polynomial & f)
{
	debug_handler_l("Fp_polynomial", "set_modulus(const Fp_polynomial&)", 2);
	if (f.Mp == 0) {
		lidia_error_handler("Fp_polynomial", "set_modulus(Fp_polynomial&)::"
				    "modulus = 0");
		return;
	}

	if (Mp != f.Mp) {
		if (Fp_polynomial::L == 0) {
			lidia_error_handler("Fp_polynomial", "set_modulus(Fp_polynomial&)"
					    "::internal error");
			//because this means that at least one of Mp and f.Mp is not Nil,
			//and thus Fp_polynomial::L cannot be Nil
			return;
		}

		if (Mp != 0)
			Fp_polynomial::L->clear(Mp);
		Mp = Fp_polynomial::L->set_to(f.Mp);
	}
	assign_zero();
	//we preserve the memory allocated for the coefficient array
}



void Fp_polynomial::set_modulus(const bigint & p)
{
	debug_handler_l("Fp_polynomial", "set_modulus(const bigint&)", 2);

	if (p < 2) {
		lidia_error_handler("Fp_polynomial", "set_modulus(bigint&)::"
				    "modulus must be >= 2");
		return;
	}

	if (modulus() != p) {
	//if (!is_prime(p,8))		// speed !
	//	lidia_error_handler( "Fp_polynomial", "set_modulus( bigint& )::modulus must be prime" );

		if (Fp_polynomial::L == 0) {
			Fp_polynomial::L = new residue_class_list< bigint >;
			memory_handler(Fp_polynomial::L, "Fp_polynomial",
				       "set_modulus(bigint &)::Error in memory allocation");
		}
		else
			if (Mp != 0)
				Fp_polynomial::L->clear(Mp);
		Mp = Fp_polynomial::L->insert(p);
	}
	assign_zero();
	//we preserve the memory allocated for the coefficient array
}



void Fp_polynomial::comp_modulus(const Fp_polynomial & x, const char* fctname)
	const
{
	//if (modulus() == 0)
	//	lidia_error_handler_c("Fp_polynomial", "comp_modulus:: modulus = 0",
	//		std::cout<<"in function "<<fctname<<std::endl;);

	if (Mp == x.Mp || modulus().is_zero())
		return;
	static char msg[60];
	strcpy(msg, fctname);
	strcat(msg, "::different moduli");
	lidia_error_handler("Fp_polynomial", msg);
}


void Fp_polynomial::set_max_degree(lidia_size_t n)
{
	// nothing happens if degree() > n
	debug_handler_l("Fp_polynomial", "set_max_degree (lidia_size_t)", 2);

	n++;
	if (c_alloc < n) {
		bigint* tmp = new bigint[n];
		memory_handler(tmp, "Fp_polynomial",
			       "set_max_degree :: Error in memory allocation");

		bigint *xp = tmp, *ap = coeff;
		lidia_size_t i;
		for (i = 0; i < c_length; i++, xp++, ap++)
			xp->assign(*ap);
		delete[] coeff;
		coeff = tmp;
		c_alloc = n;
	}
}



//***************************************************************
//
//	routines for manipulating coefficients
//
//***************************************************************

// a = rep[i], or zero if i not in range
void Fp_polynomial::get_coefficient(bigint & a, lidia_size_t i) const
{
	debug_handler_l ("Fp_polynomial", "get_coefficient (bigint&, lidia_size_t)", 3);
	if (i< 0 || i > degree())
		a.assign_zero();
	else
		a.assign(coeff[i]);
}



// returns coeff[ degree() ]
// zero if a == 0
const bigint & Fp_polynomial::lead_coeff() const
{
	debug_handler_l ("Fp_polynomial", "lead_coeff (void)", 3);
	if (is_zero())
		return ZERO;
	else
		return coeff[degree()];
}



//returns coeff[0]
// zero if a == 0
const bigint & Fp_polynomial::const_term() const
{
	debug_handler_l ("Fp_polynomial", "const_term (void)", 3);
	if (is_zero())
		return ZERO;
	else
		return coeff[0];
}



// rep[i] = a, error is raised if i < 0
void Fp_polynomial::set_coefficient(const bigint & a, lidia_size_t i)
{
	debug_handler_l ("Fp_polynomial", "set_coefficient (bigint&, lidia_size_t)", 2);

	if (Mp == 0) {
		lidia_error_handler("Fp_polynomial",
				    "set_coefficient(bigint&, lidia_size_t)::modulus = 0");
		return;
	}

	if (i < 0) {
		lidia_error_handler_c("Fp_polynomial",
				      "set_coefficient(bigint&, lidia_size_t)::negative index",
				      std::cout << "index = " << i << std::endl;);
		return;
	}

	lidia_size_t j, m = degree();
	bigint tmp;
	Remainder(tmp, a, modulus());

	if (i > m) {
		if (!tmp.is_zero()) {
			debug_handler_l ("Fp_polynomial", "set_coefficient(bigint&, lidia_size_t)::enlarging coefficient vector", 1);
			set_degree(i);
			for (j = m+1; j < i; j++)
				coeff[j].assign_zero();
			coeff[i].assign(tmp);
		}
	}
	else {
		coeff[i].assign(tmp);
		if (i == m)
			remove_leading_zeros();
	}
}



// coeff[i] = 1, error is raised if i < 0
void Fp_polynomial::set_coefficient(lidia_size_t i)
{
	debug_handler_l ("Fp_polynomial", "set_coefficient (lidia_size_t)", 2);

	if (Mp == 0) {
		lidia_error_handler("Fp_polynomial",
				    "set_coefficient(lidia_size_t)::modulus = 0");
		return;
	}

	if (i < 0) {
		lidia_error_handler_c ("Fp_polynomial",
				       "set_coefficient(lidia_size_t)::negative index",
				       std::cout << "index = " << i << std::endl;);
		return;
	}

	lidia_size_t j, m = degree();

	if (i > m) {
		debug_handler_l ("Fp_polynomial", "set_coefficient(lidia_size_t)::enlarging coefficient vector", 1);
		set_degree(i);
		for (j = m+1; j < i; j++)
			coeff[j].assign_zero();
	}

	coeff[i].assign_one();
}



void Fp_polynomial::remove_leading_zeros()
{
	debug_handler_l ("Fp_polynomial", "remove_leading_zeros (void)", 1);

	if (Mp == 0) {
		lidia_error_handler("Fp_polynomial", "remove_leading_zeros(void)::"
				    "modulus = 0");
		return;
	}

	lidia_size_t n = c_length;
	bigint *v = coeff - 1 + n;

	while (n > 0 && (*v).is_zero()) {
		v--;
		n--;
	}
	c_length = n;
}



#if 0
bigint & Fp_polynomial::operator[] (lidia_size_t i)
{
	debug_handler_l ("Fp_polynomial", "operator [] (lidia_size_t)", 3);
	if (modulus().is_zero()) {
		lidia_error_handler("Fp_polynomial", "operator[]::modulus is zero");
		return bigint();
	}

	if (i < 0) {
		lidia_error_handler_c ("Fp_polynomial", "operator[]::negative index",
				       std::cout << "index = " << i << std::endl;);
		return bigint();
	}

	if (i >= c_length) {
		lidia_error_handler_c ("Fp_polynomial", "operator[]::"
				       "index out of range",
				       std::cout << "index = " << i << "     c_length = " << c_length << std::endl;);
		return bigint();
	}

	//never: return ZERO;
	//it would discard "const" state of ZERO

	return coeff[i];
}
#endif



const bigint & Fp_polynomial::operator[] (lidia_size_t i) const
{
	debug_handler_l ("Fp_polynomial", "(const T&) operator [] (lidia_size_t)", 3);

	if (i < 0) {
		lidia_error_handler_c ("Fp_polynomial",
				       "(const bigint&) operator[]::negative index",
				       std::cout << "index = " << i << std::endl;);
		return ZERO;
	}

	if (i < c_length)
		return coeff[i];
	else
		return ZERO;
}



//***************************************************************
//
//		assignments
//
//***************************************************************


void Fp_polynomial::assign_x()
{
	debug_handler_l ("Fp_polynomial", "assign_x (void)", 2);
	if (Mp == 0) {
		lidia_error_handler("Fp_polynomial", "assign_x(void)::modulus = 0");
		return;
	}

	set_degree(1);
	coeff[0].assign_zero();
	coeff[1].assign_one();
}



void Fp_polynomial::assign_zero(const bigint & p)
{
	debug_handler_l ("Fp_polynomial", "assign_zero (bigint&)", 2);
	set_modulus(p); //assigns zero
}



void Fp_polynomial::assign_one(const bigint & p)
{
	debug_handler_l ("Fp_polynomial", "assign_one (bigint&)", 2);
	set_modulus(p);
	assign_one();
}



void Fp_polynomial::assign_x(const bigint & p)
{
	debug_handler_l ("Fp_polynomial", "assign_x (bigint&)", 2);
	set_modulus(p);
	assign_x();
}



void Fp_polynomial::assign(const Fp_polynomial & a)
{
	debug_handler_l ("Fp_polynomial", "assign (Fp_polynomial&)", 2);

	if (a.Mp == 0) {
		lidia_error_handler("Fp_polynomial", "assign(Fp_polynomial&)::"
				    "modulus = 0");
		return;
	}

	if (&a == this)
		return;

	set_modulus(a); //assigns zero

	if (!a.is_zero()) {
		set_degree(a.degree());
		lidia_size_t i;
		for (i = 0; i <= a.degree(); i++)
			coeff[i].assign(a.coeff[i]);
	}
	//this function tries to preserve the memory allocated
	//by this->coeff
}



void Fp_polynomial::assign(const base_vector< bigint > & a, const bigint & p)
{
	debug_handler_l ("Fp_polynomial", "assign (base_vector< bigint > &, const bigint&)", 2);

	lidia_size_t j, deg = a.size() - 1;
	set_modulus(p);
	set_degree(deg);

	for (j = 0; j <= deg; j++)
		Remainder(coeff[j], a[j], p);
	remove_leading_zeros();
}



void Fp_polynomial::assign(const bigint & c)
{
	debug_handler_l ("Fp_polynomial", "assign (bigint&)", 2);
	set_coefficient(c, 0);
}



Fp_polynomial & Fp_polynomial::operator = (const Fp_polynomial & p)
{
	debug_handler_l ("Fp_polynomial", "operator = (Fp_polynomial&)", 2);
	if (this != &p)
		assign(p);
	return *this;
}



Fp_polynomial & Fp_polynomial::operator = (const bigint & c)
{
	debug_handler_l ("Fp_polynomial", "operator = (Fp_polynomial&)", 2);
	set_coefficient(c, 0);
	return *this;
}



void Fp_polynomial::randomize(lidia_size_t n)
	//random polynomial of degree n
{
	debug_handler_l ("Fp_polynomial", "randomize (lidia_size_t)", 2);
	if (Mp == 0 || n < 0) {
		lidia_error_handler("Fp_polynomial",
				    "randomize (lidia_size_t)::modulus = 0  or  degree < 0");
		return;
	}

	const bigint & p = modulus();
	set_degree(n);

	lidia_size_t i;
	for (i = 0; i <= n; i++)
		coeff[i].assign(LiDIA::randomize(p));

	if (coeff[n].is_zero())
		coeff[n].assign_one();
}



void randomize(Fp_polynomial & x, const bigint & p, lidia_size_t n)
{
	debug_handler_l("Fp_polynomial", "randomize (Fp_polynomial&, const bigint&, lidia_size_t)", 3);

	x.set_modulus(p);
	x.randomize(n);
}



// swap x & y (only pointers are swapped)
void swap(Fp_polynomial & x, Fp_polynomial & y)
{
	debug_handler_l ("Fp_polynomial", "swap (Fp_polynomial&, Fp_polynomial&)", 2);

	if (&x == &y)
		return;

	bigint* b_tmp;
	lidia_size_t l_tmp;
	residue_class< bigint > *r_tmp;

	b_tmp = x.coeff;
	x.coeff = y.coeff;
	y.coeff = b_tmp;

	l_tmp = x.c_length;
	x.c_length = y.c_length;
	y.c_length = l_tmp;

	l_tmp = x.c_alloc;
	x.c_alloc = y.c_alloc;
	y.c_alloc = l_tmp;

	r_tmp = x.Mp;
	x.Mp = y.Mp;
	y.Mp = r_tmp;
}



//***************************************************************
//
//		comparisons
//
//***************************************************************

bool Fp_polynomial::is_zero() const
{
	debug_handler_l ("Fp_polynomial", "is_zero (void)", 4);
	return (c_length == 0);
}



bool Fp_polynomial::is_one() const
{
	debug_handler_l ("Fp_polynomial", "is _one (void)", 4);
	return (c_length == 1 && coeff[0].is_one());
}



bool Fp_polynomial::is_x() const
{
	debug_handler_l ("Fp_polynomial", "is_x (void)", 4);
	return (c_length == 2 && coeff[0].is_zero() && coeff[1].is_one());
}



bool Fp_polynomial::is_binomial() const
{
	lidia_size_t i, non_zeros = 1, n = c_length - 1;;
	//first non-zero coefficient: f(n)
	if (n < 1) return false;
	for (i = 0; i < n; i++)
		if (!coeff[i].is_zero()) {
			non_zeros++;
			if (non_zeros > 2) return false;
		}
	return (non_zeros == 2);
}



bool operator == (const Fp_polynomial & a,
		  const Fp_polynomial & b)
{
	debug_handler_l ("Fp_polynomial", "operator == (Fp_polynomial, Fp_polynomial)", 4);
// this code only works well if both polynomials are normalized
	lidia_size_t i, n = a.c_length;
	if (n != b.c_length || a.modulus() != b.modulus())
		return false;

	const bigint *ap = a.coeff, *bp = b.coeff;
	for (i = 0; i < n; i++)
		if (ap[i] != bp[i]) return false;

	return true;
#if 0
// the following code will work even if polynomials are not normalized
// for operator[] will return zero if index is out of range
	lidia_size_t max_deg = comparator< lidia_size_t >::max(a.c_length, b.c_length);
	lidia_size_t i;
	for (i = 0; i <= max_deg; i++)
		if (a[i] != b[i]) return 0;
	return 1;
#endif
}



bool operator != (const Fp_polynomial & a,
		  const Fp_polynomial & b)
{
	debug_handler_l ("Fp_polynomial", "operator != (Fp_polynomial, Fp_polynomial)", 4);
	return !(a == b);
}



//for class factorization< Fp_polynomial > :
bool operator < (const Fp_polynomial & a, const Fp_polynomial & b)
{
	lidia_size_t i, deg_a = a.degree(), deg_b = b.degree();
	if (deg_a != deg_b)
		return deg_a < deg_b;

	for (i = deg_a; i >= 0; i--) {
		if (a[i] < b[i]) return true;
		if (a[i] > b[i]) return false;
	}
	return false;
}



bool operator <= (const Fp_polynomial & a, const Fp_polynomial & b)
{
	lidia_size_t i, deg_a = a.degree(), deg_b = b.degree();
	if (deg_a != deg_b)
		return deg_a < deg_b;

	for (i = deg_a; i >= 0; i--) {
		if (a[i] < b[i]) return true;
		if (a[i] > b[i]) return false;
	}
	return true;
}



//********************************************************************
//
//		input and output
//
//*********************************************************************/

//
// I/O format:  "[a_0 a_1 ... a_n] mod p"
//   represents the polynomial a_0 + a_1*x + ... + a_n*x^n over Z/pZ.
//   for operator << and >>, write, read.
// pretty_print prints "x^2-x+3 mod 5" instead of "[3 -1 1] mod 5"
//

std::istream & operator >> (std::istream & s, Fp_polynomial & a)
{
	debug_handler_l ("Fp_polynomial", "operator >> (std::istream&, Fp_polynomial&)", 3);
	a.read(s);
	return s;
}



std::ostream & operator << (std::ostream & s, const Fp_polynomial & a)
{
	debug_handler_l ("Fp_polynomial", "operator << (std::ostream&, const Fp_polynomial&)", 3);
	a.pretty_print(s);
	return s;
}



void Fp_polynomial::read(std::istream & s)
{
	debug_handler_l ("Fp_polynomial", "read(std::istream&)", 2);
	bigint p;
	char c, cm, co, cd;

	s >> std::ws >> c;

	if (c != '[') {
		s.putback(c);
		pretty_read(s);
		return;
	}

	base_vector< bigint > help_coeff;
	s >> help_coeff;
	s >> cm >> co >> cd;
	if (cm != 'm' || co != 'o' || cd != 'd') {
		lidia_error_handler("Fp_polynomial", "read(std::istream&)::"
				    "'mod' expected");
		return;
	}

	s >> p;
	set_modulus(p);

	lidia_size_t i, n = help_coeff.size();
	set_degree(n-1);
	for (i = 0; i < n; i++)
		Remainder(coeff[i], help_coeff[n-i-1], p);

	remove_leading_zeros();
}



void Fp_polynomial::write(std::ostream & s) const
{
	debug_handler_l ("Fp_polynomial", "write(std::ostream&)", 2);

	s << "[";
	lidia_size_t i;
	for (i = degree(); i >= 0; i--) {
		s << coeff[i];
		if (i != 0)
			s << " ";
	}
	s << "] mod " << modulus();
}


void Fp_polynomial::pretty_print(std::ostream & os) const
{
	//written by Damian Weber
	//assumes lead_coeff != 0, all coeffs must be > 0

	debug_handler_l ("Fp_polynomial", "pretty_print(std::ostream&)", 2);

	if (is_zero()) {
		os << "0 mod " << modulus();
		return;
	}

	lidia_size_t j, deg = c_length-1;

	for (j = deg; j >= 0; j--) {
		if (coeff[j] != 0) {
			if (coeff[j] > 0 && j < deg)
				os << " + ";
			if (j == 0 || coeff[j] != 1) {
				os << coeff[j];
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



void Fp_polynomial::pretty_read(std::istream & s)
{
	// This function reads a univariate polynomial in any variable.
	// input format : a_n*x^n+ ... + a_1*x + a_0 mod p
	// e.g. :   3*x^2 + 4*x - 1 mod 17
	// Monomials need not be sorted, and powers of x may even appear
	// repeated, '*' may be omitted and coefficients may follow the variable:
	//        -4 + 8x^5 + 2 - x^2 3 + x^5 + x^5*17 mod 5
	// Note however, that the routine will work faster, if the leading monomial
	// is read first.
	// All coefficients will be reduced mod p.

	debug_handler_l("Fp_polynomial", "in member-function "
			"std::istream & read_verbose(std::istream &)", 6);

	lidia_size_t sz;
	char c;

	base_vector< bigint > help_coeff(8, vector_flags(vector_flags::expand));

	char variable = 0;
	bigint coeff_tmp(1);
	bigint tmp;

	// Read a monomial, i.e. "x^k" or "- x^k"
	// or "a*x^k" or "a x^k" or "x^k*a" or "x^k a"

	do {
		c = s.get();
	} while (isspace(c) && c != '\n');
	while (c != '\n' && c != EOF && s.good()) {
		sz = 0; // Assume we read coeffizient of x^0;
		if (c == '+') {
			coeff_tmp = 1;
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
			coeff_tmp *= tmp;
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
						lidia_error_handler("Fp_polynomial",
								    "pretty_read (...)::wrong input format");
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
					lidia_error_handler("Fp_polynomial",
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
				lidia_error_handler_c("Fp_polynomial", "pretty_read (...)::"
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
				coeff_tmp *= tmp;
				do {
					c = s.get();
				} while (isspace(c) && c != '\n');
			}

			if (c >= '0' && c <= '9' || c == '(') {
				s.putback(c);
				s >> tmp;
				coeff_tmp *= tmp;
				do {
					c = s.get();
				} while (isspace(c) && c != '\n');
			}
		}
		if (c == 'm') {
			c = s.get();
			if (c != 'o') {
				lidia_error_handler("Fp_polynomial",
						    "pretty_read (...):: 'mod' expected");
				return;
			}
			c = s.get();
			if (c != 'd') {
				lidia_error_handler("Fp_polynomial",
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
		help_coeff[sz] += coeff_tmp;
	}

	bigint p;
	s >> p;
	set_modulus(p);

	lidia_size_t n = help_coeff.size() - 1;
	set_degree(n);

	for (; n >= 0; n--)
		Remainder(coeff[n], help_coeff[n], p);
	remove_leading_zeros();
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
