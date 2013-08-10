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
//	Author	: Andreas Mueller (AM), Volker Mueller (VM)
//	Changes	: See CVS log
//
//==============================================================================================


#ifdef HAVE_CONFIG_H
# include	"config.h"
#endif
#include	"LiDIA/rational_factorization.h"



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



const unsigned int rf_single_factor::dont_know = 2;
const unsigned int rf_single_factor::prime = 1;
const unsigned int rf_single_factor::not_prime = 0;

int rational_factorization::info = 0;



rf_single_factor &
rf_single_factor::operator = (const rf_single_factor & x)
{
	single_base.assign(x.single_base);
	single_exponent = x.single_exponent;
	factor_state = x.factor_state;
	return *this;
}



bool
operator < (const rf_single_factor & a, const rf_single_factor & b)
{
	return (a.single_base < b.single_base);
}



bool
operator <= (const rf_single_factor & a, const rf_single_factor & b)
{
	return (a.single_base <= b.single_base);
}



bool
operator == (const rf_single_factor & a, const rf_single_factor & b)
{
	return (a.single_base == b.single_base);
}



void
swap (rf_single_factor & a, rf_single_factor & b)
{
	rf_single_factor c;

	c = a;
	a = b;
	b = c;
}



std::istream &
operator >> (std::istream &in, rf_single_factor &f)
{
	bigint hb;
	int he;
	char c;

	in >> c;

	if (c != '(') {
		lidia_error_handler("rational_factorization", "operator >>::(expected");
		return (in);
	}
	else {
		in >> hb >> c;
		if (c != ',') {
			lidia_error_handler("rational_factorization", "operator >>::, expected");
			return (in);
		}
		in >> he >> c;
		while (c != ')')        in >> c;
	}

	if (!he) {
		f.single_exponent = 1;
		f.single_base = 1;
		f.factor_state = rf_single_factor::prime;
	}
	else {
		f.factor_state = rf_single_factor::dont_know;
		f.single_exponent = he;

		if (hb.is_zero()) {
			lidia_error_handler("rational_factorization", "factor 0 not allowed");
			return (in);
		}

		f.single_base = hb;
	}

	return in;
}



void
rational_factorization::pretty_print (std::ostream &out)
{
	int i;

	int length = no_of_comp()-1;
	factors.sort();

	for (i = 0; i < length; i++) {
		out << factors[i].single_base;

		if (factors[i].single_exponent != 1) {
			out << "^" << factors[i].single_exponent;
		}

		out << " * ";
	}

	out << factors[length].single_base;

	if (factors[length].single_exponent != 1) {
		out << "^" << factors[length].single_exponent;
	}

}



std::ostream &
operator << (std::ostream &out, const rf_single_factor &f)
{
	out << "(";
	out << f.single_base;
	out << ",";
	out << f.single_exponent;
	out << ")";

	return out;
}



rational_factorization::rational_factorization ()
{
	factors.set_mode(EXPAND);
	isign = 0;
	decomp_state = rf_single_factor::dont_know;
}



rational_factorization::rational_factorization (int n, int exp)
{
	factors.set_mode(EXPAND);

	if (n == 0) {
		lidia_error_handler("rational_factorization", "factorization of 0 not allowed");
		return;
	}

	if (n > 0) {
		isign = 1;
		factors[0].single_base = bigint(n);
	}
	else {
		isign = -1;
		factors[0].single_base = bigint(-n);
	}
	if (!(exp & 1))
		isign = 1;

	if (factors[0].single_base.is_one() || exp == 0) {
		factors[0].single_base.assign_one();
		exp = 1;
		factors[0].factor_state = rf_single_factor::prime;
		decomp_state = rf_single_factor::prime;
	}
	else {
		decomp_state = rf_single_factor::dont_know;
		factors[0].factor_state = rf_single_factor::dont_know;
	}

	factors[0].single_exponent = exp;
}



rational_factorization::rational_factorization (unsigned int n, int exp)
{
	factors.set_mode(EXPAND);

	if (n == 0) {
		lidia_error_handler("rational_factorization", "factorization of 0 not allowed");
		return;
	}

	isign = 1;

	if (n == 1 || exp == 0) {
		factors[0]. single_base.assign_one();
		exp = 1;
		factors[0].factor_state = rf_single_factor::prime;
		decomp_state = rf_single_factor::prime;
	}
	else {
		factors[0].factor_state = rf_single_factor::dont_know;
		decomp_state = rf_single_factor::dont_know;
	}
	factors[0].single_base = bigint(static_cast<unsigned long>(n));
	factors[0].single_exponent = exp;
}



rational_factorization::rational_factorization (long n, int exp)
{
	factors.set_mode(EXPAND);

	if (n == 0) {
		lidia_error_handler("rational_factorization", "factorization of 0 not allowed");
		return;
	}

	if (n > 0) {
		isign = 1;
		factors[0].single_base = bigint(n);
	}
	else {
		isign = -1;
		factors[0].single_base = bigint(-n);
	}
	if (!(exp & 1))
		isign = 1;

	if (factors[0].single_base.is_one() || exp == 0) {
		factors[0].single_base.assign_one();
		exp = 1;
		factors[0].factor_state = rf_single_factor::prime;
		decomp_state = rf_single_factor::prime;
	}
	else {
		decomp_state = rf_single_factor::dont_know;
		factors[0].factor_state = rf_single_factor::dont_know;
	}

	factors[0].single_exponent = exp;
}



rational_factorization::rational_factorization (unsigned long n, int exp)
{
	factors.set_mode(EXPAND);

	if (n == 0) {
		lidia_error_handler("rational_factorization", "factorization of 0 not allowed");
		return;
	}
	isign = 1;
	if (n == 1 || exp == 0) {
		factors[0].single_base.assign_one();
		exp = 1;
		decomp_state = rf_single_factor::prime;
		factors[0].factor_state = rf_single_factor::prime;
	}
	else {
		decomp_state = rf_single_factor::dont_know;
		factors[0].factor_state = rf_single_factor::dont_know;
	}
	factors[0].single_base = bigint(n);
	factors[0].single_exponent = exp;
}



rational_factorization::rational_factorization (const bigint & n, int expo)
{
	factors.set_mode(EXPAND);

	if (n.is_zero()) {
		lidia_error_handler("rational_factorization", "factorization of 0 not allowed");
		return;
	}

	if (expo == 0) {
		assign(1);
		return;
	}

	if (n.is_gt_zero()) {
		isign = 1;
		factors[0].single_base = n;
	}
	else {
		isign = -1;
		factors[0].single_base = -n;
	}

	if (!(expo & 1))  // exponent is even
		isign = 1;

	if (factors[0].single_base.is_one()) {
		factors[0].factor_state = rf_single_factor::prime;
		factors[0].single_exponent = 1;
		decomp_state = rf_single_factor::prime;
	}
	else {
		decomp_state = rf_single_factor::dont_know;
		factors[0].factor_state = rf_single_factor::dont_know;
		factors[0].single_exponent = expo;
	}
}



rational_factorization::rational_factorization (const bigrational & n, int expo)
{
	factors.set_mode(EXPAND);

	if (!(expo & 1))
		isign = 1;
	else
		isign = n.sign();
	decomp_state = rf_single_factor::dont_know;

	if (expo == 0) {
		assign(1);
		return;
	}

	if (!n.denominator().is_one()) {
		factors[0].single_base = n.denominator();
		factors[0].single_exponent = -expo;
		factors[0].factor_state = rf_single_factor::dont_know;
	}

	if (!abs(n.numerator()).is_one()) {
		factors[1].single_base = abs(n.numerator());
		factors[1].single_exponent = expo;
		factors[1].factor_state = rf_single_factor::dont_know;
	}

	compose();
}



rational_factorization::rational_factorization (const rational_factorization & f,
						int exp)
{
	int i, l = f.no_of_comp();

	if (exp == 0) {
		assign(1);
		return;
	}

	factors.set_mode(EXPAND);
	factors = f.factors;
	for (i = 0; i < l; i++)
		factors[i].single_exponent *= exp;
	isign = f.isign;
	decomp_state = f.decomp_state;
}



rational_factorization::~rational_factorization ()
{
}


rational_factorization &
rational_factorization::operator = (const rational_factorization & f)
{
	factors = f.factors;
	isign = f.isign;
	decomp_state = f.decomp_state;

	return *this;
}



void
rational_factorization::assign (long n, int exp)
{
	if (n == 0) {
		lidia_error_handler("rational_factorization", "factorization of 0 not allowed");
		return;
	}

	factors.set_size(1);

	if (n > 0) {
		isign = 1;
		factors[0].single_base = bigint(n);
	}
	else {
		isign = -1;
		factors[0].single_base = bigint(-n);
	}

	if (!(exp & 1))
		isign = 1;

	if (exp == 0)
		factors[0].single_base = bigint(1);

	if (factors[0].single_base.is_one()) {
		factors[0].factor_state = rf_single_factor::prime;
		exp = 1;
		decomp_state = rf_single_factor::prime;
	}
	else {
		decomp_state = rf_single_factor::dont_know;
		factors[0].factor_state = rf_single_factor::dont_know;
	}
	factors[0].single_exponent = exp;
}



void
rational_factorization::assign (const bigint & n, int expo)
{
	if (n.is_zero()) {
		lidia_error_handler("rational_factorization", "factorization of 0 not allowed");
		return;
	}

	if (expo == 0) {
		assign(1L);
		return;
	}

	factors.set_size(1);

	if (n.is_gt_zero()) {
		isign = 1;
		factors[0].single_base = n;
	}
	else {
		isign = -1;
		factors[0].single_base = -n;
	}

	if (!expo & 1)  // exponent is even
		isign = 1;


	if (factors[0].single_base.is_one()) {
		factors[0].factor_state = rf_single_factor::prime;
		decomp_state = rf_single_factor::prime;
		factors[0].single_exponent = 1;
	}
	else {
		decomp_state = rf_single_factor::dont_know;
		factors[0].factor_state = rf_single_factor::dont_know;
		factors[0].single_exponent = expo;
	}
}



void
rational_factorization::convert(base_vector< bigint > & basis, base_vector< int > & exponent)
{
	lidia_size_t i, len;

	if (basis.size() != exponent.size()) {
		lidia_error_handler("rational_factorization", "different sizes of basis and exponent vectors");
		return;
	}
	else {
		len = basis.size();
		factors.set_size(len);
		for (i = 0; i < len; i++) {
			factors[i].single_base = basis[i];
			factors[i].single_exponent = exponent[i];
			factors[i].factor_state = rf_single_factor::dont_know;
		}
		compose();
	}
}



void
rational_factorization::assign (const bigrational & n, int expo)
{
	if (n.is_zero()) {
		lidia_error_handler("rational_factorization", "factorization of 0 not allowed");
		return;
	}

	if (expo == 0) {
		assign(1L);
		return;
	}

	if (n.denominator().is_one()) {
		assign(n.numerator(), expo);
		return;
	}

	if (!(expo & 1))
		isign = 1;
	else
		isign = n.sign();

	decomp_state = rf_single_factor::dont_know;

	factors.set_size (2);

	factors[0].single_base = n.denominator();
	factors[0].single_exponent = -expo;
	factors[0].factor_state = rf_single_factor::dont_know;

	factors[1].single_base = abs(n.numerator());
	factors[1].single_exponent = expo;
	factors[1].factor_state = rf_single_factor::dont_know;

	compose();
}



void
rational_factorization::assign (const rational_factorization &f, int exp)
{
	int i, l = f.no_of_comp();

	if (exp == 0) {
		assign(1L);
		return;
	}

	factors = f.factors;

	for (i = 0; i < l; i++)
		factors[i].single_exponent *= exp;
	isign = f.isign;
	decomp_state = f.decomp_state;
}



bigint
rational_factorization::base (lidia_size_t index) const
{
	if ((index< 0) || (index >= no_of_comp())) {
		lidia_error_handler("rational_factorization", "base::index out of range");
		return (0);
	}

	return bigint(factors[index].single_base);
}



int
rational_factorization::exponent (lidia_size_t index) const
{
	if ((index< 0) || (index >= no_of_comp())) {
		lidia_error_handler("rational_factorization", "exponent::index out of range");
		return (0);
	}

	return factors[index].single_exponent;
}



void
rational_factorization::set_exponent (lidia_size_t index, int expo)
{
	if ((index< 0) || (index >= no_of_comp())) {
		lidia_error_handler("rational_factorization", "set_exponent::index out of range");
		return;
	}

	factors[index].single_exponent = expo;

	if (expo == 0) {
		compose();
	}
}



bool
rational_factorization::is_prime_factor (lidia_size_t index)
{
	if (index< 0 || index >= no_of_comp()) {
		lidia_error_handler("rational_factorization", "is_prime_factor::index out of range");
		return false;
	}

	if (rf_single_factor::state(factors[index].factor_state) == rf_single_factor::prime)
		return true;

	if (rf_single_factor::state(factors[index].factor_state) == rf_single_factor::not_prime)
		return false;

	if (is_prime(factors[index].single_base, 8)) {
		factors[index].factor_state = rf_single_factor::prime;
		return true;
	}
	else {
		factors[index].factor_state = rf_single_factor::not_prime;
		return false;
	}
}



bool
rational_factorization::is_prime_factorization ()
{

	if (decomp_state == rf_single_factor::prime)
		return true;

	for (lidia_size_t i = no_of_comp()-1; i >= 0; i--) {
		if (rf_single_factor::state(factors[i].factor_state) == rf_single_factor::not_prime) {
			decomp_state = rf_single_factor::not_prime;
			return false;
		}

		if (rf_single_factor::state(factors[i].factor_state) == rf_single_factor::dont_know) {
			if (is_prime(factors[i].single_base, 8))
				factors[i].factor_state = rf_single_factor::prime;
			else {
				decomp_state = rf_single_factor::not_prime;
				factors[i].factor_state = rf_single_factor::not_prime;
				return false;
			}
		}
	}

	decomp_state = rf_single_factor::prime;
	return true;
}



void
rational_factorization::compose ()
{
	lidia_size_t i, j, len = no_of_comp() - 1;

	sort();

	for (i = 0; i < len; i++) {
		j = i+1;
		while ((j <= len) && (factors[i].single_base == factors[j].single_base)) {
			factors[i].single_exponent += factors[j].single_exponent;
			j++;
		}

		if ((j = j - i - 1) > 0) {
			factors.remove_from(i+1, j);
			len = no_of_comp() - 1;
		}
		if (factors[i].single_exponent == 0 || factors[i].single_base.is_one()) {
			factors.remove_from(i);
			len --;
			i --;
		}
	}

        if (factors[len].single_exponent == 0 || factors[len].single_base.is_one()) {
          factors.remove_from(len);
        }
        
	if (no_of_comp() == 0) {
		factors[0].single_base.assign(1);
		factors[0].single_exponent = 1;
		factors[0].factor_state = rf_single_factor::prime;
		decomp_state = rf_single_factor::prime;
	}
}



void
rational_factorization::compose (sort_vector< rf_single_factor > & v)
{
	lidia_size_t i, j, len = v.size() - 1;

	v.sort();

	for (i = 0; i < len; i++) {
		j = i+1;
		while ((j <= len) && (v[i].single_base == v[j].single_base)) {
			v[i].single_exponent += v[j].single_exponent;
			j++;
		}

		if ((j = j - i - 1) > 0) {
			v.remove_from(i+1, j);
			len = v.size() - 1;
		}
		if (v[i].single_exponent == 0 || v[i].single_base.is_one()) {
			v.remove_from(i);
			len --;
		}
	}

        if (v[len].single_exponent == 0 || v[len].single_base.is_one()) {
          v.remove_from(len);
        }

	if (v.size() == 0) {
		v[0].single_base.assign(1);
		v[0].single_exponent = 1;
		v[0].factor_state = rf_single_factor::prime;
	}
	v.sort(); // Testing
}



std::istream &
operator >> (std::istream & in, rational_factorization & f)
{
	lidia_size_t n = 0;
	int sg = 1;
	char c;

	f.decomp_state = rf_single_factor::dont_know;

	f.factors.set_size (0);

	in >> c;

	if (c != '[') {
		lidia_error_handler("rational_factorization", "operator >>::[ expected");
		return(in);
	}

	in >> c;

	while (c != ']') {
		in.putback(c);

		in >> f.factors[n];

		n++;

		in >> c;
	}

	for (lidia_size_t i = 0; i < n; i++) {
		if (f.factors[i].single_base.is_lt_zero()) {
			if (f.factors[i].single_exponent & 1)
				sg = -sg;
			f.factors[i].single_base.negate();
		}

		if ((f.factors[i].single_base).is_one()) {
			f.factors.remove_from(i);
			n--;
			i--;
		}
	}

	f.isign = sg;

	f.compose();

	return in;
}



std::ostream &
operator << (std::ostream & out, const rational_factorization & f)
{
	lidia_size_t i, sz = f.no_of_comp();
	rational_factorization g(f);

	g.factors.sort();

	out << "  [";

	if (g.isign == -1)
		out << "(-1, 1)";

	for (i = 0; i < sz; i++)
		out << g.factors[i];

	out << "]";

	return out;
}



void
rational_factorization::invert ()
{
	lidia_size_t i, len = no_of_comp();

	for (i = 0; i < len; i++)
		factors[i].single_exponent = - factors[i].single_exponent;
}



void
rational_factorization::square ()
{
	lidia_size_t i, len = no_of_comp();

	for (i = 0; i < len; i++)
		factors[i].single_exponent <<= 1;
}



void
multiply (rational_factorization & f,
	  const rational_factorization & a,
	  const rational_factorization & b)
{
	f.factors.concat(a.factors, b.factors);

	if (a.decomp_state == rf_single_factor::prime && b.decomp_state == rf_single_factor::prime)
		f.decomp_state = rf_single_factor::prime;
	else
		f.decomp_state = rf_single_factor::dont_know;

	f.isign = a.isign * b.isign;
	f.compose();
}



void
divide (rational_factorization & f,
	const rational_factorization & a,
	const rational_factorization & b)
{
	rational_factorization h(b);
	h.invert();
	f.factors.concat(a.factors, h.factors);

	if (a.decomp_state == rf_single_factor::prime && b.decomp_state == rf_single_factor::prime)
		f.decomp_state = rf_single_factor::prime;
	else
		f.decomp_state = rf_single_factor::dont_know;

	f.isign = a.isign * b.isign;
	f.compose();
}


#if 0
void gcd(rational_factorization & f,
	 const rational_factorization & a,
	 const rational_factorization & b)
{
	rational_factorization ha(a), hb(b);
	lidia_size_t i, j;

	f.assign(1L);

	ha.refine();
	hb.refine();

	for (i = 0; i < ha.no_of_comp(); i++)
		if (ha.exponent(i) < 0)
			lidia_error_handler("rational_factorization", "gcd::negative exponent");

	for (j = 0; j < hb.no_of_comp(); j++)
		if (hb.exponent(j) < 0)
			lidia_error_handler("rational_factorization", "gcd::negative exponent");

	for (i = 0; i < ha.no_of_comp(); i++)
		for (j = 0; j < hb.no_of_comp(); j++)
			if (ha.base(i) == hb.base(j))
				multiply(f, f, rational_factorization(ha.base(i), min(ha.exponent(i), hb.exponent(j))));
}


void lcm (rational_factorization & f,
	  const rational_factorization & a,
	  const rational_factorization & b)
{
	rational_factorization h;

	multiply(h, a, b);
	gcd(f, a, b);
	divide(f, h, f);
}
#endif



void
rational_factorization::refine2 (sort_vector< rf_single_factor > & v , rf_single_factor & SF ,
				 const rf_single_factor & a, const rf_single_factor & b)

	// we assume vector v to be of expanding size !!!

	// compute refinement of {a,b} and store the
	// last rf_single_factor of this factorization individually
	// in SF
	// thus {v[0],...,v[n],SF} represents the refinement of {a,b}
{
	debug_handler ("rational_factorization" , "refine2()");

	lidia_size_t i = 0, last_index = 2;
	int j;
	rf_single_factor sf;
	bigint q, r;

	v.set_size (2);
	v[0] = a;
	v[1] = b;

	while (i < last_index-1) {
		sf.single_base = gcd(v[i].single_base, v[i+1].single_base);

		if (sf.single_base.is_one())
			i++;
		else {
			div_rem(v[i].single_base, r, v[i].single_base, sf.single_base);
			j = 0;
			while (r.is_zero()) {
				div_rem(q, r, v[i].single_base, sf.single_base);
				if (r.is_zero())
					v[i].single_base.assign(q);
				j++;
			}
			sf.single_exponent = j * v[i].single_exponent;

			div_rem(v[i+1].single_base, r, v[i+1].single_base, sf.single_base);
			j = 0;
			while (r.is_zero()) {
				div_rem(q, r, v[i+1].single_base, sf.single_base);
				if (r.is_zero())
					v[i+1].single_base.assign(q);
				j++;
			}
			sf.single_exponent += j * v[i+1].single_exponent;

			if (is_prime(v[i].single_base, 8))
				v[i].factor_state = rf_single_factor::prime;
			else v[i].factor_state = rf_single_factor::not_prime;

			if (is_prime(v[i+1].single_base, 8))
				v[i+1].factor_state = rf_single_factor::prime;
			else v[i+1].factor_state = rf_single_factor::not_prime;

			if (sf.single_exponent != 0) {
				if (is_prime(sf.single_base, 8))
					sf.factor_state = rf_single_factor::prime;
				else sf.factor_state = rf_single_factor::not_prime;

				if (v[i+1].single_base.is_one())
					v[i+1] = sf;
				else {
					v.insert_at (sf , i+1);
					last_index ++;
				}
			}
		}
	}

	// extract last element and return it as an individual

	i = v.size()-1;
	SF = v[i];
	v.set_size(i);

	compose(v);
}



void
rational_factorization::refine ()
{
	debug_handler ("rational_factorization" , "refine()");

	sort_vector< rf_single_factor > v(0, EXPAND), w(0, EXPAND);
	rf_single_factor b, c, sf;

	lidia_size_t i, j, k, l = factors.size();

	if (l >= 2) {
		refine2 (v , sf , factors[0] , factors[1]);

		k = v.size();
		v[k] = sf;

		for (i = 2; i < l; i++) {
			j = 0;
			b = v[0];
			c = factors[i];

			do
			{
				// compute refinement of {b,c}
				refine2 (w , sf , b , c);


				// insert w into v and proceed with {..,sf}
				k = w.size();

				v.shift_right (j+1 , k-1);

				v.assign      (j , w , 0 , k-1);

				j += k;

				if (j < v.size()) {
					b = v[j];
					c = sf;
				}
				else  if (! sf.single_base.is_one()) {
					v[j] = sf;
					j++;
				}
			}
			while (j < v.size() && ! sf.single_base.is_one());
		}

		swap (v , factors);
		compose ();
	}
}



bool
rational_factorization::refine (const bigint & x)
{
	lidia_size_t i, j, len = no_of_comp();
	int e_akt;

	bool rc = false;

	bigint old, n, res, d;

	rf_single_factor sf;

	j = len;

	for (i = 0; i < len; i++) {
		e_akt = 0;
		old = base(i);
		d = gcd (old , x);

		if (! d.is_one() && d != old) {
			rc = true;

			div_rem(n, res, old, d);

			while (res.is_zero()) {
				e_akt++;
				old = n;

				div_rem(n, res, old, d);
			}

			sf.single_base = d;
			sf.single_exponent = e_akt * exponent(i);

			if (is_prime(d, 8))
				sf.factor_state = rf_single_factor::prime;
			else
				sf.factor_state = rf_single_factor::not_prime;


			if (! old.is_one()) {
				factors[i].single_base = old;

				if (is_prime(old, 8))
					factors[i].factor_state = rf_single_factor::prime;
				else
					factors[i].factor_state = rf_single_factor::not_prime;

				// * * * append sf to the end of factors * * *

				factors[j] = sf;
				j ++;
			}
			else {
				factors[i] = sf;
			}
		}
	}

	compose();

	return rc;
}



bool
rational_factorization::refine_comp (lidia_size_t index, const bigint & x)
{
	int e = 0;
	bigint old, n, res, d;

	rf_single_factor sf;

	bool rc = false;

	if ((index< 0) || (index >= no_of_comp())) {
		lidia_error_handler("rational_factorization", "refine_comp(int, const bigint &)::index out of range");
		return (false);
	}

	old = base(index);
	d = gcd (old , x);

	if (! d.is_one() && d != old) {
		rc = true;

		div_rem(n, res, old, d);

		while (res.is_zero()) {
			e++;
			old = n;
			div_rem(n, res, old, x);
		}

		sf.single_base = d;
		sf.single_exponent = e * exponent(index);

		if (is_prime(d, 8))
			sf.factor_state = rf_single_factor::prime;
		else
			sf.factor_state = rf_single_factor::not_prime;

		if (! old.is_one()) {
			factors[index].single_base = old;

			if (is_prime(old, 8))
				factors[index].factor_state = rf_single_factor::prime;
			else
				factors[index].factor_state = rf_single_factor::not_prime;

			factors[no_of_comp()] = sf;
		}
		else {
			factors[index] = sf;
		}

		compose();
	}

	return rc;
}



bool
operator == (const rational_factorization & a,
	     const rational_factorization & b)
{
	rational_factorization f;

	divide(f, a, b);
	f.refine();

	if (f.no_of_comp() == 1 && f.base(0).is_one() && f.isign == 1)
		return 1;
	else
		return 0;
}



void
rational_factorization::factor_comp (lidia_size_t index, int upper_bound)
{
	if ((index< 0) || (index >= no_of_comp())) {
		lidia_error_handler("rational_factorization", "factor_comp::index out of range");
		return;
	}

	if (is_prime_factor(index)) {
		if (info)
			std::cout << "prime number " << factors[index].single_base << "\n";
		return;
	}

	if (upper_bound > 34) {
		lidia_warning_handler("rational_factorization", "factor_comp::incorrect parameters");
		upper_bound = 34;
	}

	int D, job_buffer[30][5], strategy[30];
	int jobs_avail, k;
	int n = (factors[index].single_base.bit_length() / 3) + 10;

	// try to factor component using TD
	if (upper_bound == 34) {
		char *n_string;
		n_string = new char[n];
		upper_bound = static_cast<int>(bigint_to_string(factors[index].single_base, n_string) * .28) + 1;
		delete [] n_string;

		if (upper_bound > 34)
			upper_bound = 34;
		if (upper_bound < 6)
			upper_bound = 6;
	}

	if ((D = ecm_read_max(upper_bound)) < 1000000)
		D = 1000000;

	ecm_primes prim(1, D+200, 200000);
	trialdiv(index, 1, 1000000, prim);

	if (is_prime_factor(index)) {
		if (info)
			std::cout << "prime number " << factors[index].single_base << "\n";

		compose();
		return;
	}
	// check whether base[index] is an exact power
	bigint N;

	do {
		k = static_cast<int>(is_power(N, factors[index].single_base));

		if (k > 0) {
			factors[index].single_base.assign(N);
			factors[index].single_exponent *= k;

			if (is_prime(N, 8)) {
				factors[index].factor_state = rf_single_factor::prime;
				compose();
				return;
			}
			else
				factors[index].factor_state = rf_single_factor::not_prime;
		}
	}  while (k > 0);

	// try to factor component using ECM
	if (!(rf_single_factor::state(factors[index].factor_state) == rf_single_factor::prime)) {
		k = 0; n = 6; // lower_bound = 6

		while (n < upper_bound) {
			strategy[k++] = n;
			n += 3;
		}

		strategy[k] = upper_bound;
		strategy[k+1] = 0;

		jobs_avail = ecm_job_planing(strategy, job_buffer);

		N.assign(factors[index].single_base);
		int ecm_restlength = static_cast<int>(N.bit_length()*0.177);
		// *0.30 wegen Dezimalstellen und * 0.59 wegen Restlaenge

		if (ecm_restlength < 20)   // always use ecm completely
			ecm_restlength = 1;

		ecm(index, ecm_restlength, jobs_avail, job_buffer, prim);
	}

	if (!(is_prime_factor(index))) {
		prim.resetprimes(1);
		mpqs(index, prim);
	}

	compose();
}



void
rational_factorization::factor (int upper_bound)
{
	if (is_prime_factorization())
		return;

	if (upper_bound > 34) {
		lidia_warning_handler("rational_factorization", "factor::upper_bound > 34");
		upper_bound = 34;
	}

	lidia_size_t k = no_of_comp(), index;

	int D, job_buffer[30][5], strategy[30];
	int jobs_avail;
	int n = (factors[k-1].single_base.bit_length() / 3) + 10;

	// try to factor component using TD
	if (upper_bound == 34) {
		char *n_string;
		n_string = new char[n];
		upper_bound = static_cast<int>(bigint_to_string(factors[k-1].single_base, n_string) * .28) + 1;
		delete [] n_string;

		if (upper_bound > 34)
			upper_bound = 34;
		if (upper_bound < 6)
			upper_bound = 6;
	}

	if ((D = ecm_read_max(upper_bound)) < 1000000)
		D = 1000000;

	ecm_primes prim(1, D+200, 200000);
	bigint N;

	D = 0; n = 6; // lower_bound = 6

	while (n < upper_bound) {
		strategy[D++] = n;
		n += 3;
	}

	strategy[D] = upper_bound;
	strategy[D+1] = 0;

	jobs_avail = ecm_job_planing(strategy, job_buffer);

	for (index = 0; index < k; index++) {
		if (info) {
			std::cout << "\nworking on index " << index;
			std::cout << "\n==================\n";
			std::cout.flush();
		}

		if (rf_single_factor::state(factors[index].factor_state) == rf_single_factor::dont_know)
			if (is_prime_factor(index) == rf_single_factor::prime) {
				if (info)
					std::cout << "prime number " << factors[index].single_base << "\n";
				continue;
			}

		int ecm_restlength = static_cast<int>(factors[index].single_base.bit_length()*0.177); // *0.30 wegen Dezimalstellen und * 0.59 wegen Restlaenge

		if (ecm_restlength < 20)
			ecm_restlength = 1;

		trialdiv(index, 1, 1000000, prim);

		if (is_prime_factor(index)) {
			if (info)
				std::cout << "prime number " << factors[index].single_base << "\n";
			continue;
		}

		// check whether base[index] is an exact power
		do {
			D = static_cast<int>(is_power(N, factors[index].single_base));

			if (D > 0) {
				factors[index].single_base.assign(N);
				factors[index].single_exponent *= D;

				if (is_prime(N, 8)) {
					factors[index].factor_state = rf_single_factor::prime;
					D = -1;
				}
				else
					factors[index].factor_state = rf_single_factor::not_prime;
			}
		} while (D > 0);

		// try to factor component using ECM
		if (!(rf_single_factor::state(factors[index].factor_state) == rf_single_factor::prime)) {
			if (info)   std::cout << "ECM ";

			ecm(index, ecm_restlength, jobs_avail, job_buffer, prim);
		}

		for (D = 0; D < jobs_avail; D++)
			job_buffer[D][3] = job_buffer[D][4];

		if (!(is_prime_factor(index))) {
			prim.resetprimes(1);
			mpqs(index, prim);
		}
	}

	compose();

	if (is_prime_factorization())
		return;

	int trials = 0;

	k = no_of_comp();

	for (index = 0; index < k; index++) {
		if (rf_single_factor::state(factors[index].factor_state) == rf_single_factor::not_prime) {
			N.assign(factors[index].single_base);

			if (info)
				std::cout << "\nindex " << index << " not prime: to factor " << N;

			for (n = 0; n < jobs_avail; n++)
				job_buffer[n][3] = job_buffer[n][4];

			ecm(index, jobs_avail, job_buffer, prim);
			compose();

			if (N == factors[index].single_base)
				trials ++;
			else trials = 0;
			k = no_of_comp();
			index = -1;

			if (trials >= 2) return;
		}
	}
}



rational_factorization
factor (const bigint &N)
{
	rational_factorization f(N);

	f.factor();

	return f;
}



rational_factorization
factor_completely (const bigint &N)
{
	rational_factorization f(N);

	f.factor_completely();

	return f;
}



//********** NEW *********************************************

bigrational
evaluate (const rational_factorization & r)
{
	bigrational x(1);
	bigint h;
	int i, j;

	for (i = 0; i < r.no_of_comp(); i++) {
		j = r.exponent(i);

		if (j > 0) {
			if (j > 1) {
				power(h, r.base(i), j);
				multiply(x, x, h);
			}
			else
				multiply(x, x, r.base(i));
		}
		else {
			j = -j;
			if (j > 1) {
				power(h, r.base(i), j);
				divide(x, x, h);
			}
			else
				divide(x, x, r.base(i));
		}
	}
	if (r.sign() == -1)
		negate(x, x);

	return x;
}



bigint
evaluate_to_bigint (const rational_factorization & r)
{
	bigrational x;

	x = evaluate(r);

	if (!x.denominator().is_one()) {
		lidia_error_handler("rational_factorization", "evaluate::rational_factorization does not evaluate to bigint");
		return 1;
	}

	return x.numerator();
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
