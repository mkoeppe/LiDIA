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
//	Author	: Detlef Anton (DA), Thomas Pfahler (TPf)
//	Changes	: See CVS log
//
//==============================================================================================


#ifdef HAVE_CONFIG_H
# include	"config.h"
#endif
#include	"LiDIA/gf_element.h"
#include	"LiDIA/finite_fields/galois_field_rep.h"
#include	"LiDIA/Fp_polynomial.h"
#include	"LiDIA/timer.h"
#include	"LiDIA/power_functions.h"
#include	"LiDIA/gf_polynomial.h"
#include	<cctype>



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



// Since we often have to call member functions (of the class galois_field_rep)
// via function pointers, we use the following #define to increase readibility.
// If you don't understand the syntax, go and learn C++ ...
// Just kidding :-) If A::foo is a pointer to a member function of A (with
// argument x), then for an object "a" of type A the function is called via
// "a.*(a.foo)(x)"; if "p" is a pointer to objects of type A this becomes
// "p->*(p->foo)(x)".
#define CALL(R, fct) (R->*(R->fct))


//***********************************************************************
//*			class gf_element				*
//***********************************************************************


galois_field gf_element::uninitialized_field;
unsigned int gf_element::output_format = 0;


//
// constructors and destructor
//

gf_element::gf_element ()
	: ff(uninitialized_field),
	  rep(0)
{
	galois_field_rep const* R = get_ff_rep();
	CALL(R, construct) (*this);

	if (gf_polynomial::FIELD != 0)
		set_field(*gf_polynomial::FIELD);
}



gf_element::gf_element (const gf_element& a)
	: ff(a.ff),
	  rep(0)
{
	galois_field_rep const* R = get_ff_rep();
	CALL(R, copy) (*this, a);
}



gf_element::gf_element (const galois_field& K)
	: ff(K),
	  rep(0)
{
	galois_field_rep const* R = get_ff_rep();
	CALL(R, construct) (*this);
}



gf_element::~gf_element()
{
	galois_field_rep const* R = get_ff_rep();
	CALL(R, destruct) (*this);
}



//
// access functions
//

const galois_field &
gf_element::get_field () const
{
	return ff;
}



bigint
gf_element::characteristic () const
{
	return ff.characteristic();
}



bigint
gf_element::lift_to_Z () const
{
//	galois_field_rep *R = get_ff_rep();
//	return CALL(R, lift) (*this);
	Fp_polynomial f = polynomial_rep();
	if (f.degree() > 0)
		lidia_error_handler("gf_element",
				    "lift_to_Z: element not in prime field");
	return f.const_term();
}



const Fp_polynomial &
gf_element::polynomial_rep2 () const
{
	galois_field_rep const* R = get_ff_rep();
	return CALL(R, get_pol_rep) (*this);
	// we have immediate access to the representation, but the return value is
	// only valid until the next operation involving (*this), so be careful!
}



Fp_polynomial
gf_element::polynomial_rep () const
{
	return polynomial_rep2(); // return a copy
}



void
gf_element::set_polynomial_rep (const Fp_polynomial& a)
{
	if (a.modulus() != ff.characteristic())
		lidia_error_handler("gf_element",
				    "set_polynomial_rep: wrong characteristic");

	galois_field_rep const* R = get_ff_rep();
        CALL(R, set_pol_rep) (*this, a);
}



//
// assignments
//

void
gf_element::set_field (const galois_field& K)
{
	if (ff == K) {
		assign_zero();
		return;
	}

	galois_field_rep const* R = get_ff_rep();
	CALL(R, destruct) (*this);
	ff.assign(K);
	R = get_ff_rep();
	CALL(R, construct) (*this);
}



void
gf_element::assign_zero ()
{
	galois_field_rep const* R = get_ff_rep();
	CALL(R, as0) (*this);
}



void
gf_element::assign_one ()
{
	galois_field_rep const* R = get_ff_rep();
	CALL(R, as1) (*this);
}



void
gf_element::assign_zero (const galois_field& K)
{
	set_field(K);
	//assign_zero();
}



void
gf_element::assign_one (const galois_field& K)
{
	set_field(K);
	assign_one();
}



void
gf_element::assign (const gf_element& a)
{
	galois_field_rep const* R = get_ff_rep();
	if (ff == a.ff)
		CALL(R, copy) (*this, a);
	else {
//	set_field(a.ff);		simple solution
//	assign(a);
//	return;

		CALL(R, destruct) (*this); // 'fast' solution
		ff.assign(a.ff);
		R = get_ff_rep();
		CALL(R, copy) (*this, a);
	}
}



void
gf_element::assign (const bigint& a)
{
	assign(promote(a));
}



gf_element&
gf_element::operator = (const gf_element& a)
{
	assign(a);
	return *this;
}



gf_element&
gf_element::operator = (const bigint& a)
{
	assign(a);
	return *this;
}



//
// operators
//

gf_element
operator - (const gf_element & a)
{
	gf_element h(a);

	h.negate();
	return h;
}



gf_element
operator + (const gf_element & a, const gf_element & b)
{
	gf_element h;

	add(h, a, b);
	return h;
}



gf_element
operator - (const gf_element & a, const gf_element & b)
{
	gf_element h;

	subtract(h, a, b);
	return h;
}



gf_element
operator * (const gf_element & a, const gf_element & b)
{
	gf_element h;

	multiply(h, a, b);
	return h;
}



gf_element
operator / (const gf_element & a, const gf_element & b)
{
	gf_element h;

	divide(h, a, b);
	return h;
}



//
// operators for type gf_element
//

gf_element &
gf_element::operator += (const gf_element & a)
{
	add(*this, *this, a);
	return *this;
}



gf_element &
gf_element::operator -= (const gf_element & a)
{
	subtract(*this, *this, a);
	return *this;
}



gf_element &
gf_element::operator *= (const gf_element & a)
{
	multiply(*this, *this, a);
	return *this;
}



gf_element &
gf_element::operator /= (const gf_element & a)
{
	divide(*this, *this, a);
	return *this;
}



//
// operators for type bigint
//

gf_element
operator + (const gf_element & a, const bigint & b)
{
	gf_element h;

	add(h, a, b);
	return h;
}



gf_element
operator - (const gf_element & a, const bigint & b)
{
	gf_element h;

	subtract(h, a, b);
	return h;
}



gf_element
operator * (const gf_element & a, const bigint & b)
{
	gf_element h;

	multiply(h, a, b);
	return h;
}



gf_element
operator / (const gf_element & a, const bigint & b)
{
	gf_element h;

	divide(h, a, b);
	return h;
}



gf_element
operator + (const bigint & a, const gf_element & b)
{
	gf_element h;

	add(h, a, b);
	return h;
}



gf_element
operator - (const bigint & a, const gf_element & b)
{
	gf_element h;

	subtract(h, a, b);
	return h;
}



gf_element
operator * (const bigint & a, const gf_element & b)
{
	gf_element h;

	multiply(h, a, b);
	return h;
}



gf_element
operator / (const bigint & a, const gf_element & b)
{
	gf_element h;

	divide(h, a, b);
	return h;
}



gf_element &
gf_element::operator += (const bigint & a)
{
	add(*this, *this, a);
	return *this;
}



gf_element &
gf_element::operator -= (const bigint & a)
{
	subtract(*this, *this, a);
	return *this;
}



gf_element &
gf_element::operator *= (const bigint & a)
{
	multiply(*this, *this, a);
	return *this;
}



gf_element &
gf_element::operator /= (const bigint & a)
{
	divide(*this, *this, a);
	return *this;
}



//
// procedural versions
//

gf_element
gf_element::promote (const bigint& a) const
{
	gf_element res(ff);
	Fp_polynomial pol;
	galois_field_rep const* R = get_ff_rep();

	pol.set_modulus(R->characteristic());
	pol[0] = a;
	CALL(R, set_pol_rep) (res, pol);
	return res;
}



gf_element
gf_element::promote (const gf_element& a) const
{
	lidia_error_handler("gf_element", "promote(gf_element&): not implemented");

	gf_element r;

	return r;
}



void
add (gf_element & c, const gf_element & a, const gf_element & b)
{
	galois_field_rep const* R = a.get_ff_rep();
	const galois_field& Ka = a.ff;
	const galois_field& Kb = b.ff;

	if (Ka == Kb) {
		if (Ka != c.ff)
			c.set_field(Ka);
		CALL(R, add) (c, a, b);
	}
	else if (Kb < Ka) {
		c.set_field(Ka);
		CALL(R, add) (c, a, a.promote(b));
	}
	else if (Ka < Kb) {
		c.set_field(Kb);
		R = b.get_ff_rep();
		CALL(R, add) (c, b.promote(a), b);
	}
	else
		lidia_error_handler("gf_element", "add: different fields");
}



void
subtract (gf_element & c, const gf_element & a, const gf_element & b)
{
	galois_field_rep const* R = a.get_ff_rep();
	const galois_field& Ka = a.ff;
	const galois_field& Kb = b.ff;

	if (Ka == Kb) {
		if (Ka != c.ff)
			c.set_field(a.ff);
		CALL(R, sub) (c, a, b);
	}
	else if (Kb < Ka) {
		c.set_field(Ka);
		CALL(R, sub) (c, a, a.promote(b));
	}
	else if (Ka < Kb) {
		c.set_field(Kb);
		R = b.get_ff_rep();
		CALL(R, sub) (c, b.promote(a), b);
	}
	else
		lidia_error_handler("gf_element", "subtract: different fields");
}



void
multiply (gf_element & c, const gf_element & a, const gf_element & b)
{
	galois_field_rep const* R = a.get_ff_rep();
	const galois_field& Ka = a.ff;
	const galois_field& Kb = b.ff;

	if (Ka == Kb) {
		if (Ka != c.ff)
			c.set_field(a.ff);
		CALL(R, mul) (c, a, b);
	}
	else if (Kb < Ka) {
		c.set_field(Ka);
		CALL(R, mul) (c, a, a.promote(b));
	}
	else if (Ka < Kb) {
		c.set_field(Kb);
		R = b.get_ff_rep();
		CALL(R, mul) (c, b.promote(a), b);
	}
	else
		lidia_error_handler("gf_element", "multiply: different fields");
}



void
divide (gf_element & c, const gf_element & a, const gf_element & b)
{
	if (b.is_zero())
		lidia_error_handler("gf_element", "divide: division by zero");

	if (a.ff != b.ff)
		lidia_error_handler("gf_element", "divide: different fields");

	gf_element inv_b(b);

	inv_b.invert();
	multiply(c, a, inv_b);
}



void
add (gf_element & c, const bigint & a, const gf_element & b)
{
	add(c, b.promote(a), b);
}



void
add (gf_element & c, const gf_element & a, const bigint & b)
{
	add(c, a, a.promote(b));
}



void
subtract (gf_element & c, const bigint & a, const gf_element & b)
{
	subtract(c, b.promote(a), b);
}



void
subtract (gf_element & c, const gf_element & a, const bigint & b)
{
	subtract(c, a, a.promote(b));
}



void
multiply (gf_element & c, const bigint & a, const gf_element & b)
{
	multiply(c, b.promote(a), b);
}



void
multiply (gf_element & c, const gf_element & a, const bigint & b)
{
	multiply(c, a, a.promote(b));
}



void
divide (gf_element & c, const bigint & a, const gf_element & b)
{
	divide(c, b.promote(a), b);
}



void
divide (gf_element & c, const gf_element & a, const bigint & b)
{
	divide(c, a, a.promote(b));
}



void
gf_element::multiply_by_2 ()
{
	galois_field_rep const* R = this->get_ff_rep();

	CALL(R, add) (*this, *this, *this);
}



void
gf_element::divide_by_2 ()
{
	galois_field_rep const* R = this->get_ff_rep();

	if (R->characteristic() == 2)
		lidia_error_handler("gf_element", "divide_by_2"
				    "division by 2 not possible in fields of "
				    "char. 2");
	gf_element h(*this);
	h.assign((R->characteristic() + 1) >> 1);
	CALL(R, mul) (*this, *this, h);
}



void
gf_element::negate ()
{
	galois_field_rep const* R = get_ff_rep();

	CALL(R, neg) (*this);
}



void
negate (gf_element & a, const gf_element & b)
{
	a.assign(b);
	a.negate();
}



void
gf_element::invert ()
{
	if (is_zero())
		lidia_error_handler("gf_element", "invert: element is zero");
	if (is_one())
		return;

	galois_field_rep const* R = get_ff_rep();
	CALL(R, inv) (*this);
}



void
invert (gf_element & a, const gf_element & b)
{
	a.assign(b);
	a.invert();
}



gf_element
inverse (const gf_element & a)
{
	gf_element h(a);

	h.invert();
	return h;
}



void
square (gf_element & a, const gf_element & b)
{
	if (a.ff != b.ff)
		a.set_field(b.ff);
	galois_field_rep const* R = b.get_ff_rep();
	CALL(R, sqr) (a, b);
}



void
gf_element::power (const bigint & e)
{
	lidia_power(*this, *this, e);
}



void
power (gf_element & c, const gf_element & a, const bigint & e)
{
	c.assign(a);
	c.power(e);
}



void
pth_power (gf_element & c, const gf_element & a, lidia_size_t e)
{
	bigint q;
	const bigint &p = a.get_field().characteristic();

	power(q, p, e);
	power(c, a, q);
}



gf_element
sqrt (const gf_element &a)
	// algorithm implemented by Volker Mueller
	// see Bach/Shallit, Alg. of Cipolla
{
	const galois_field& ff = a.ff;
	unsigned int i;
	gf_element t(ff);

	if (a.ff.characteristic() == 2) {
		t.assign(a);
		for (i = 1; i < ff.degree(); i++)
			square (t, t);
		return t;
	}

	gf_element res0(ff), res1(ff), h(ff), h2(ff);
	bigint exp = ff.number_of_elements();

	if (a.is_zero())
		return a;

	multiply(h, a, 4);
	do {
		t.randomize();
	} while ((t*t-h).is_square()); //while is_square(t^2-4a)

	inc(exp); shift_right(exp, exp, 1);

	res1.assign_one(); res0.assign_zero();

	for (i = exp.bit_length()-2; i > 0; i--) {
		// fast exponentiation left-right
		h.assign(res0);
		square(h2, h);
		square(res0, res1);
		multiply(res0, res0, a);
		subtract(res0, h2, res0); // res0 = h*h - a*res1*res1;

		add(h, h, h);
		multiply(h2, res1, t);
		add(h, h, h2);
		multiply(res1, res1, h); // res1 = res1*(res1*t + 2*h);

		if (exp.bit(i)) {
			h.assign(res0);
			multiply(res0, res1, a);
			negate(res0, res0); // res0 = -res1*a;

			multiply(res1, res1, t);
			add(res1, res1, h); // res1 = res1*t + h;
		}
	}

	h.assign(res0);
	square(h2, h);
	square(res0, res1);
	multiply(res0, res0, a);
	subtract(res0, h2, res0); // res0 = h*h - a*res1*res1;

	add(h, h, h);
	multiply(h2, res1, t);
	add(h, h, h2);
	multiply(res1, res1, h); // res1 = res1*(res1*t + 2*h);

	if (exp.is_odd()) {
		multiply(res0, res1, a);
		negate(res0, res0); // res0 = -res1*a;
	}

	return (res0);
}



//
// comparisons
//

bool
gf_element::operator == (const gf_element & a) const
{
	if (ff == a.ff) {
		galois_field_rep const* R = get_ff_rep();
		return CALL(R, iseq) (*this, a);
	}
	return false;
}



bool
gf_element::operator != (const gf_element & a) const
{
	return !(*this == a);
}



bool
gf_element::operator == (const bigint & a) const
{
	if (a.is_zero())
		return this->is_zero();
	if (a.is_one())
		return this->is_one();

	gf_element tmp = promote(a);
	return *this == tmp;
}



bool
gf_element::operator != (const bigint & a) const
{
	return !(*this == a);
}



bool
gf_element::is_zero () const
{
	galois_field_rep const* R = get_ff_rep();

	return CALL(R, is0) (*this);
}



bool
gf_element::is_one () const
{
	galois_field_rep const* R = get_ff_rep();

	return CALL(R, is1) (*this);
}



//
// basic functions
//

void
gf_element::randomize (unsigned int deg)
{
	unsigned int n = ff.degree();

	if (n % deg != 0)
		lidia_error_handler("gf_element", "randomize(unsigned int d): "
				    "d must be a positive divisor of field degree");

	Fp_polynomial pol;
	galois_field_rep const* R = get_ff_rep();

	pol.set_modulus(R->characteristic());
	pol.randomize(n); //pol[n] = 1, pol[n-1], ..., pol[0] = random
	pol[n] = 0; //yield random polynomial of degree < n
	CALL(R, set_pol_rep) (*this, pol);

	if (deg < n) {
		bigint e1, e2;

		e1 = R->number_of_elements();
		dec(e1);

		e2 = static_cast<int>(deg);
		lidia_power(e2, R->characteristic(), e2);
		dec(e2);
		divide(e1, e1, e2);
		power(e1);
	}
}



void
swap (gf_element & a, gf_element & b)
{
	swap(a.ff, b.ff);

	void * p = a.rep;

	a.rep = b.rep;
	b.rep = p;
}



//
// high level functions
//


bigint
gf_element::order () const
{
	if (is_zero()) {
		lidia_error_handler("gf_element", "order(): zero is not an element "
				    "of the multiplicative group");
		return bigint(0);
	}
	if (is_one())
		return bigint(1);


	// now compute order by dividing (p^n-1) by all prime factors of group order

	lidia_size_t i, j;
	bigint q = ff.number_of_elements() - 1;
	bigint e_ord(1), h;
	gf_element res(ff);

	const rational_factorization &F = ff.factorization_of_mult_order();
	for (i = 0; i < F.no_of_comp(); i++) {
		LiDIA::power(h, F.base(i), F.exponent(i));
		divide(h, q, h);

		LiDIA::power(res, *this, h);
		j = 0;

		while (j < F.exponent(i) && !res.is_one()) {
			j++;
			multiply(e_ord, e_ord, F.base(i));
			LiDIA::power(res, res, F.base(i));
		}
		if (!res.is_one()) {
			lidia_error_handler("gf_element",
					    "order()::element order does not divide (q-1)");
			return bigint(0);
		}
	}

	return e_ord;
}



multi_bigmod
gf_element::trace () const
{
	unsigned int d = ff.degree();
	const bigint &p = ff.characteristic();
	bool char2 = (p == 2);

	gf_element a(*this), sum(*this);
	unsigned int i;
	for (i = 1; i < d; i++) {
		if (char2)
			LiDIA::square(a, a);
		else
			LiDIA::power(a, a, p);
		LiDIA::add(sum, sum, a);
	}

//    if (sum.relative_degree()) != 1)
//    	lidia_error_handler("gf_element", "trace(): trace not in prime	field");

	return multi_bigmod(sum.lift_to_Z(), p);
}



multi_bigmod
gf_element::norm () const
{
	const bigint &p = ff.characteristic();
	multi_bigmod m;

	m.set_modulus(p);

	if (p == 2) {
		//N(elem)=0 <==> (elem)=0  for char. 2
		if (!is_zero())
			m.set_mantissa(1);
	}
	else {
		// N(elem) = prod_{i=0}^{n-1} elem^(p^i)
		//         = elem^(sum_{i=0}^{n-1} p^i)
		//         = elem^((p^n-1)/(p-1))
		gf_element N(*this);
		if (ff.degree() > 1)
			N.power((ff.number_of_elements()-1) / (p-1));
		m.set_mantissa(N.lift_to_Z());
	}
	return m;
}



bool
gf_element::is_primitive_element () const
{
	return this->order() == (ff.number_of_elements()-1);
}



bool
gf_element::is_free_element () const
{
	if (ff.degree() == 1)
		return !is_zero(); // all nonzero elts are free elts in Fp

	// Gao[1993] (acc. to D.Anton's thesis):
	// elem is free element <=> gcd(f,g)==1 (in F_p[X](!))
	// where f(X) = X^n - 1
	// and   h(X) = \sum_{i=0}^{n-1}  trace(elem^(p^i + 1))*X^i
	const bigint &p = ff.characteristic();
	lidia_size_t n = ff.degree();
	Fp_polynomial f, h;
	f.set_modulus(p);
	f[n] = 1;
	f[0] = -1; // f(X) = X^n - 1 mod p
	h.set_modulus(f);
	h.set_max_degree(n-1);

	lidia_size_t i;
	gf_element h_i_elem(ff), pp(ff), temp(ff);

	square(h_i_elem, *this);
	h.set_coefficient(h_i_elem.trace().mantissa(), 0);
	for (i = 1; i < n; i++) {
		pth_power(pp, temp, 1);
		multiply(h_i_elem, *this, pp);
		h.set_coefficient(h_i_elem.trace().mantissa(), i);
		temp.assign(pp);
	}		// h(X) = \sum_{i = 0}^{n-1}  trace(elem^(p^i + 1))*X^i

	Fp_polynomial x;
	gcd(x, f, h);
	return (x.is_one());
}



bool
gf_element::is_square () const
{
	if (is_zero() || ff.characteristic() == 2)
		return true;

	if (ff.degree() == 1)
		return (jacobi(lift_to_Z(), ff.characteristic()) != -1);

	bigint e = (ff.number_of_elements()-1) >> 1; // e = (q-1)/2

	gf_element b(*this);
	b.power(e);
	return (b.is_one());
}



static bigint
greatest_divisor_of_1_coprime_to_2 (const bigint & a, const bigint & b)
{
	bigint al(a);
	bigint bl(b);

	while (!(bl.is_one())) {
		bl.assign(gcd(al, bl));
		divide(al, al, bl);
	}
	return al;
}



void
gf_element::assign_primitive_element (const galois_field& K)
{
	set_field(K);

	bigint q = K.number_of_elements() - 1;
	bool char2 = (K.characteristic() == 2);

	if (char2 && K.degree() == 1) {
		// trivial case: GF(2)
		assign_one();
		return;
	}
	if (char2 && is_prime(q))	// if the group order is prime, then any
	{				// nontrivial element is a generator
		do randomize();
		while (is_one() || is_zero());
		return;
	}
	const rational_factorization & F = K.factorization_of_mult_order();

	timer t;
	t.start_timer();

	bigint a, b, f, g, A, B, Bl;

	gf_element x   (ff);
	gf_element y   (ff);
	gf_element z   (ff);
	gf_element h   (ff);
	gf_element pow (ff);
	do
		x.randomize();
	while (x.is_zero());
	bigint e = x.order();
	if (e == q) {
		this->assign(x);
		t.stop_timer(); goto weiter;
		//	return;
	}

	for (;;) {
		do {
			y.randomize();
			LiDIA::power(pow, y, e);
		} while (pow.is_one() || pow.is_zero());
		f = y.order();
		if (f == q) {
			this->assign(y);
			t.stop_timer(); goto weiter;
			//	    return;
		}

		g = lcm(e, f);
		LiDIA::divide(a, f, gcd(e, f));
		A = greatest_divisor_of_1_coprime_to_2(e, a);
		LiDIA::divide(a, e, gcd(e, f));
		Bl = greatest_divisor_of_1_coprime_to_2(f, a);
		LiDIA::divide(B, Bl, gcd(A, Bl));
		LiDIA::divide(a, e, A);
		LiDIA::divide(b, f, B);
		LiDIA::power(h, y, b);
		LiDIA::power(z, x, a);
		LiDIA::multiply(z, z, h);
		x = z;
		e = g;
		if (e == q) {
			this->assign(x);
			t.stop_timer(); goto weiter;
			//	    return;
		}
	}
 weiter:
	return;
	std::cout << "assign_primitive_element, method 1: " << t << std::endl;

	gf_element apow(ff), ahelp(ff);
	t.start_timer();
	int length = F.no_of_comp();
	bigint hh, h1, *v = new bigint[length];

	// compute in v the difference between ord/p_i for prime factors p_i of F
	divide(v[0], q, F.base(length-1));
	hh.assign(v[0]);
	int i = 1;
	while (i < length) {
		divide(h1, q, F.base(length-1-i));
		subtract(v[i++], h1, hh);
		hh.assign(h1);
	}

	bool found = false;
	while (found == false) {
		// check for generator 
		do
			randomize();
		while (is_zero() || is_one());

		LiDIA::power(apow, *this, v[0]);
		if (apow.is_one())
			continue;

		i = 1;
		while (i < length) {
			LiDIA::power(ahelp, *this, v[i++]);
			LiDIA::multiply(apow, apow, ahelp);
			if (apow.is_one())
				break;
		}
		if ((i == length) && (apow.is_one() == false))
			found = true;
	}
	delete[] v;

	t.stop_timer();
	std::cout << "assign_primitive_element, method 2: " << t << std::endl;

	return;

	// this point should never be reached...
	lidia_error_handler("gf_element",
			    "assign_primitive_element(): no primitive element found");
	return;
}



udigit
hash (const gf_element &a)
{
	galois_field_rep const* R = a.get_ff_rep();
	return CALL(R, hash) (a);

//    udigit s = 0;
//    lidia_size_t i;
//    for (i = 0; i < a.pol_rep.degree(); i++)
//	s ^= static_cast<udigit>(a.pol_rep[i].least_significant_digit());
//    return s % max_udigit();
}


unsigned int
gf_element::absolute_degree () const
{
	return ff.degree();
}



unsigned int
gf_element::relative_degree () const
{
	unsigned int i, n = ff.degree();

	if (n == 1 || is_one() || is_zero())
		return 1;

	const bigint &p = ff.characteristic();
	bool char2 = (p == 2);
	gf_element a(ff);
	if (char2)
		LiDIA::square(a, *this);
	else
		LiDIA::power(a, *this, p);

	for (i = 1; i <= n; i++) {
		if ((n % i == 0) && (a == *this))
			return i;
		if (char2)
			LiDIA::square(a, a);
		else
			LiDIA::power(a, a, p);
	}
	lidia_error_handler("galois_field",
			    "relative_degree(): no subfield found");
	return 0;
}



bool
gf_element::solve_quadratic (const gf_element& a1, const gf_element& a0)
// Input : polynomial X^2 + a1*X + a0
// Output: true, and *this is assigned a root, if polynomial splits,
//         false otherwise
{
	if (a1.ff != a0.ff)
		lidia_error_handler("gf_element",
				    "solve_quadratic(...): different fields");
	if (ff != a1.ff)
		set_field(a1.ff);

	// check trivial case

	if (a0.is_zero()) {
		assign_zero();
		return true;
	}

	if (ff.characteristic() == 2) {
		// (copied from an implementation for GF(2^n) by V. Mueller)
		if (a1.is_zero()) {
			assign(LiDIA::sqrt(a0)); // all elements are squares in char. 2
			return true;
		}

		gf_element c(ff), h(ff);
		gf_element c0(ff), c1(ff); // coeffs of c1*z + c0 mod (z^2+z+c)
		gf_element res0(ff), res1(ff); // coeffs of res1*z + res0 mod (z^2+z+c)
		unsigned int i, d, degree = ff.degree(), counter = 0;

		// variable transformation: solve z^2 + z + c = 0 for c= a0/a1^2
		// --> root = z*a1 is solution for original equation
		square(c, a1);
		c.invert();
		multiply(c, c, a0);

		// char = 2 ==> deg 2 - polynomials splits iff trace = 0

		if (c.trace().mantissa() != 0)
			return false;

		d = 1;
		c1.randomize(degree);
		res1.assign(c1);
		c0.randomize(degree);
		res0.assign(c0);


		do {
			// start trace computation for random polynomial
			// c1*z + 0 modulo (z^2+z+c) 
			counter ++;

			for (i = 1; i < degree; i++) {
				// (c1*z+c0) < -- (c1*z+c0)^2
				square(c1, c1);
				square(c0, c0);
				multiply(h, c1, c);
				add(c0, c0, h);
				add(res0, res0, c0);
				add(res1, res1, c1);
			}
			square(c1, c1);
			square(c0, c0);
			multiply(h, c1, c);
			add(c0, c0, h);

			if (!res1.is_zero()) {
				// degree one->check whether root found
				if (!res1.is_one())      // make monic
					divide(res0, res0, res1);

				// z + res0 should be factor of z^2+z+c, check this
				square(h, res0);
				add(h, h, res0);
				add(h, h, c);
				if (h.is_zero())  // root 'res0' found
					break;
			}

			if (counter >= 8) {
				// try smaller subfields
				counter = 0;
				if (d == degree)
					d = 1;
				else
					d++;

				while (degree % d != 0 && d < degree)
					d++;
			}

			c1.randomize(degree / d);
			res1.assign(c1);
			c0.randomize(degree);
			res0.assign(c0);
		} while (true);

		multiply(*this, res0, a1);
		return (true);
	}
	else {
		// odd characteristic
		gf_element delta(ff), c(ff);
		square(delta, a1);
		subtract(delta, delta, 4*a0);

		if (!delta.is_square())
			return false;

		assign(sqrt(delta));

		assign((-a1 + sqrt(delta))/2);
		return true;
	}
}



//
// I/O	   case (1): prime field
//      	(2): extension field of odd characteristic
//	 	(3): GF(2^n), n > 1
//
// I/O format (verbose):
//     (1, 2, 3): (g(X), f(X)), where
//		g(X) is the polynomial representation of the element and
//		f(X) = get_field().irred_polynomial()
//	additionally:
//	(1) (m, p)
//	examples: "(4, X mod 17)",
//		  "(3*X^2-2, X^2-3 mod 5)",
//		  "(X^10+x^7+1, X^11+X^2+1 mod 2)",
//		  "(4, 17)"
// I/O format (short):
//	(1) (m), where m is an integer
//	(2) (g(X)), where g(X) is a polynomial with integer coefficients
//	(3a) (Dec:d), where d = sum_{i} a_i 2^i represents the polynomial
//		sum_{i} a_i X^i mod 2
//	(3b) (Hex:d) (analogously)
//	examples: "(4)", "(3*X^2-2)", "(Dec:1153)", "(Hex:481)"
//
// Default output format is short, and (3a) for GF(2^n).
// Short input is only allowed if the field of the element is already set(!).
//

void gf_element::set_output_format(unsigned int i)
{
	if (i == 0)
		output_format = SHORT_OUTPUT;
	else
		output_format = VERBOSE_OUTPUT;
}



std::istream &
operator >> (std::istream & in, gf_element & a)
{
	a.input(in);
	return in;
}



std::ostream &
operator << (std::ostream & out, const gf_element & a)
{
	a.output(out);
	return out;
}

bigint
gf_element::return_as_bigint() const
{
    if( ff.degree() != 1 )
	{
		std::cerr << "The degree has to be 1 to invoke this method."
				  << std::endl;
		return( 0 );
	}	

    if (ff == uninitialized_field)
    {
		std::cerr << "The field is not initialized."
				  << std::endl;

	    return( 0 );
    }

    Fp_polynomial tmp = polynomial_rep();
    polynomial< bigint > pol;
    pol.set_degree(tmp.degree());
    for (int i = tmp.degree(); i >= 0; i--)
	    pol[i] = tmp[i];
    
    return( pol[0] );
}



void
gf_element::output (std::ostream & out) const
{
	out << "(";

	if (ff == uninitialized_field) {
		if (output_format == SHORT_OUTPUT)
			out << "0";
		else
			out << "0, 0";
	}
	else if (output_format == SHORT_OUTPUT && ff.characteristic() == 2) {
		out << "Dec: ";
		Fp_polynomial tmp = polynomial_rep();
		bigint h;
		for (int i = tmp.degree(); i >= 0; i--) {
			h <<= 1;
			if (tmp[i] != 0)
				h |= 1;
		}
		out << h;
	}
	else {
//	out << polynomial_rep().lift()
		Fp_polynomial tmp = polynomial_rep();
		polynomial< bigint > pol;
		pol.set_degree(tmp.degree());
		for (int i = tmp.degree(); i >= 0; i--)
			pol[i] = tmp[i];
		out << pol;

		if (output_format == VERBOSE_OUTPUT)
			out << ", " << ff.irred_polynomial();
	}
	out << ")";
}



static void
input_error (char c)
{
	static char msg[] = "input: 'X' expected";
	msg[8] = c; //overwrite character 'X'
	lidia_error_handler("gf_element", msg);
}



static void
convert_dec (polynomial< bigint > &pol, const bigint h)
	// h = sum_{i} a_i 2^i  -->  pol = sum_{i} a_i X^i
{
	int i, l = h.bit_length();
	pol.set_degree(l);
	for (i = l; i >= 0; i--) {
		if (h.bit(i) == 1)
			pol[i] = 1;
		else
			pol[i] = 0;
	}
}



static void
read_dec (std::istream &in, polynomial< bigint > &pol)
	// read decimal repr. for elements in GF(2^n)
{
	char c;

	in >> c;
	if (toupper(c) != 'D') input_error('D');
	in >> c;
	if (toupper(c) != 'E') input_error('e');
	in >> c;
	if (toupper(c) != 'C') input_error('c');
	in >> c;
	if (c != ':') input_error(':');
	bigint h;
	in >> h;
	convert_dec(pol, h);
}



static void
read_hex (std::istream &in, polynomial< bigint > &pol)
	// read hex repr. for elements in GF(2^n)
{
	char c;

	in >> c;
	if (toupper(c) != 'H') input_error('H');
	in >> c;
	if (toupper(c) != 'E') input_error('e');
	in >> c;
	if (toupper(c) != 'X') input_error('x');
	in >> c;
	if (c != ':') input_error(':');

	bigint h;
	char s[1000];
	int i = 0, j;

	do {
		in >> c;
		s[i++] = toupper(c);
		if (i == 1000)
			lidia_error_handler("gf_element", "input: input too large while "
					    "reading hex representation");
	} while (c != ')');
	in.putback(c);
	s[--i] = '\0';

	for (i = 0; i < static_cast<int>(strlen(s)); i++) {
		c = s[i];
		multiply(h, h, 16);
		j = c - '0';
		if (j >= 0 && j < 10) {
			add(h, h, j);
			continue;
		}
		j = c - 'A';
		if (j >= 0 && j < 6) {
			add(h, h, j+10);
			continue;
		}
		else {
			lidia_error_handler("gf_element", "input: wrong input character "
					    "while reading hex representation");
			return;
		}
	}
	convert_dec(pol, h);
}



void
gf_element::input (std::istream & in)
{
	if (gf_polynomial::FIELD != 0)
		set_field(*gf_polynomial::FIELD);

	Fp_polynomial        pol_Fp;
	polynomial< bigint > pol_Z;
	bool short_input = false;
	bigint h;
	int i, l;
	char c, *s;

	in >> c;
        if (c != '(')  // Input as a bigint, p-adic representation
         {
         if (! (c == '-' || ( c >= '0' && c <= '9') || c == '+'))
                   input_error('(');
         else
                  in.putback(c);

         bigint pol_as_bigint;  // read the bigint first

	 in >> pol_as_bigint;

	 if (pol_as_bigint.is_negative())
	   lidia_error_handler("gf_element", "p-adic representation of"
			       " negative number not supported");

	 int i = pol_as_bigint.bit_length() / 
	   ff.characteristic().bit_length() + 10; // approx. for pol-degree

	 base_vector < bigint > v(i, EXPAND);
	 
	 i = 0;
	 while(!pol_as_bigint.is_zero())
	   {
	     v.insert_at(pol_as_bigint % ff.characteristic(), i++);
	     pol_as_bigint /= ff.characteristic();
	   }
	 pol_Z = polynomial < bigint > (v);

         short_input = true;
        } 
        else 
        {
	in >> c;
	in.putback(c);
	c = toupper(c);
	switch (c) {
	case 'D':
		//lidia_debug_handler("gf_element", "input: call read_dec");
		short_input = true;
		read_dec(in, pol_Z);
		break;
	case 'H':
		//lidia_debug_handler("gf_element", "input: call read_hex");
		short_input = true;
		read_hex(in, pol_Z);
		break;
	default:
		in >> pol_Z;
		in >> c;
		if (c == ')') { short_input = true; in.putback(c); break; }
		short_input = false;
		if (c != ',') input_error(',');
		in >> c;
		in.putback(c);
		if (isdigit(c))	{
			// expect either a) (..., p) or
			//               b) (..., m*x^i+...)
			//           (we are here ^)
			in >> h;
			in >> c;
			in.putback(c);
			if (c == ')') {
				// yep, case a)
				pol_Fp.assign_x(h);
				break;
			}
			else {
				// no, push integer back->case b)
				s = new char[h.bit_length()/3 + 10];
				l = bigint_to_string(h, s);
				for (i = l-1; i >= 0; i--)
					in.putback(s[i]);
				delete[] s;
			}
		}
		in >> pol_Fp; // read polynomial (also for case b))
	}
	in >> c;
	if (c != ')') input_error(')');
	}
	if (short_input) {
		// make sure that ff is initialized
		if (ff.degree() == 0)
			lidia_error_handler("gf_element", "input: no field specified "
					    "while reading short input format");
	}
	else {
		if ((ff == uninitialized_field) ? 1 : (pol_Fp != ff.irred_polynomial())) {
			// need to assign a new field
			galois_field K(pol_Fp);
			set_field(K);
		}
	}
	Fp_polynomial rep(pol_Z, ff.characteristic()); // reduce pol_Z mod p
	set_polynomial_rep(rep); // assign
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
