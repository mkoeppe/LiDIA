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
//	Author	:
//	Changes	: See CVS log
//
//==============================================================================================


#ifndef LIDIA_FACTORIZATION_CC_GUARD_
#define LIDIA_FACTORIZATION_CC_GUARD_



#ifndef LIDIA_FACTORIZATION_H_GUARD_
# include	"LiDIA/factorization.h"
#endif



#ifdef LIDIA_NAMESPACE
# ifndef IN_NAMESPACE_LIDIA
namespace LiDIA {
# endif
#endif



#ifdef HAVE_HP_UX_CC

//
// This function is only used in factorization < T >::sort().
// The ppair function compare cannot be used, because
// HP UX 10.20 CC cannot resolve the template parameters.
// MM
//

template< class T >
inline int
ppair_sf_lst_compare (const T &p, const T &q)
{
	int res = 0;
	if (p.left() > q.left()) res = 1;
	if (p.left() < q.left()) res = -1;

	if (res == 0) {
		if (p.right() > q.right()) res = 1;
		if (p.right() < q.right()) res = -1;
	}
	return res;
}

#endif



//******************************************************************************
// sometimes, we have to assign '1' or to compare with '1'
// in such cases, we create a single_factor< T > named ONE, which must be '1'
// by default (see single_factor.h)
//******************************************************************************


//
// constructors, destructor
//

template< class T >
factorization< T >::factorization () :
	prime_component(0, vector_flags(vector_flags::expand)),
	composite_component(0, vector_flags(vector_flags::expand))
{
	epsilon = new single_factor< T >; // default value is '1'
	if (!epsilon)
		lidia_error_handler("factorization< T >", "factorization() :: out of memory");
	check_for_trivial_case();
}



template< class T >
factorization< T >::factorization (const single_factor< T > & f) :
	prime_component(0, vector_flags(vector_flags::expand)),
	composite_component(0, vector_flags(vector_flags::expand))
{
	epsilon = new single_factor< T >; // default value is '1'
	if (!epsilon)
		lidia_error_handler("factorization< T >", "factorization(single_factor< T > &) :: out of memory");
	check_for_trivial_case();
	assign(f);
}



template< class T >
factorization< T >::factorization (const factorization< T > & f) :
	prime_component(f.prime_component, vector_flags(vector_flags::expand)),
	composite_component(f.composite_component, vector_flags(vector_flags::expand)),
	attributes(f.attributes)
{
	epsilon = new single_factor< T >; // default value is '1'
	if (!epsilon)
		lidia_error_handler("factorization< T >", "factorization(factorization< T > &) :: out of memory");
	epsilon->assign(*f.epsilon);
	set_prime_flag(f.prime_flag());
}



template< class T >
void
factorization< T >::reset ()
{
	delete epsilon;
	epsilon = new single_factor< T >; // default value is '1'
	prime_component.set_size(0);
	composite_component.set_size(0);
	set_prime_flag(prime);
	check_for_trivial_case();
}



template< class T >
void
factorization< T >::kill ()
{
	delete epsilon;
	epsilon = new single_factor< T >; // default value is '1'
	prime_component.kill();
	composite_component.kill();
	set_prime_flag(prime);
	check_for_trivial_case();
}



//
// private functions
//

template< class T >
void
factorization< T >::check_for_trivial_case ()
{
	if (no_of_composite_components() == 0)
		set_prime_flag(prime);

	switch(no_of_components()) {
	case (0) :
	case (1) : attributes = (sorted | normalized | refined);
		if (no_of_composite_components() == 1)
			set_prime_flag(composite_component[0].left().prime_flag());
		break;
	case (2) : if (no_of_composite_components() == 1) {
		attributes = (sorted);
//               if (prime_component(0) != composite_component(0))
//                   attributes |= (normalized);
	}
	default  : break;
	}
}



//
// assignment
//

template< class T >
void
factorization< T >::assign (const single_factor< T > &f)
{
	single_factor< T > tmp(f);
	epsilon->assign(tmp.extract_unit());
	if (!tmp.is_one()) {
		if (tmp.is_prime_factor()) {
			prime_component[0].left() = tmp;
			prime_component[0].right() = 1;
			prime_component.set_size(1);
			composite_component.kill();
		}
		else {
			composite_component[0].left() = tmp;
			composite_component[0].right() = 1;
			composite_component.set_size(1);
			prime_component.kill();
		}
	}
	check_for_trivial_case();
}



template< class T >
void
factorization< T >::assign (const factorization< T > &f)
{
	if (this == &f)
		return;

	*epsilon = *f.epsilon;
	prime_component.assign(f.prime_component);
	composite_component.assign(f.composite_component);
	set_prime_flag(f.prime_flag());
	attributes = f.attributes;
}



//
// access functions
//

template< class T >
T
factorization< T >::value () const
{
	if (no_of_components() == 0)
		return unit();

	T x;
	lidia_size_t i, j;
	lidia_size_t npc = no_of_prime_components();
	lidia_size_t ncc = no_of_composite_components();

// this "strange" way of computing the value is necessary for classes
// such as Fp_polynomial where a modulus (or similar information)
// must be set before any computation
// if we did it the usual way, almost always a runtime error would occur

	if (npc > 0) {
		x = prime_base(npc-1).base();
		for (j = 0; j < prime_exponent(npc-1)-1; j++)
			multiply(x, x, prime_base(npc-1).base());
		npc--;
	}
	else {
		x = composite_base(ncc-1).base();
		for (j = 0; j < composite_exponent(ncc-1)-1; j++)
			multiply(x, x, composite_base(ncc-1).base());
		ncc--;
	}

	for (i = 0; i < npc; i++) {
		for (j = 0; j < prime_exponent(i); j++)
			multiply(x, x, prime_base(i).base());
	}
	for (i = 0; i < ncc; i++) {
		for (j = 0; j < composite_exponent(i); j++)
			multiply(x, x, composite_base(i).base());
	}
	if (!epsilon->is_one())
		multiply(x, x, epsilon->base());
	return x;
}



//
// predicates
//

template< class T >
bool
factorization< T >::is_prime_factorization (int test)
{
	switch (prime_flag()) {
	case (not_prime) : return false;

	case (prime) : return true;

	case (unknown) :;
	}

	if (test == 0)
		return (no_of_composite_components() == 0);

	// test != 0 =  > test all composite elements
	lidia_size_t k;
	for (lidia_size_t i = 0; i < no_of_composite_components(); i++) {
		if (composite_component[i].left().is_prime_factor(test)) {
			//found an element that is definitely prime
			k = no_of_prime_components();
			prime_component[k].left().assign(composite_component[i].left());
			prime_component[k].right() = composite_component[i].right();
			composite_component.remove_from(i);
			i--;
		}
		else {
			//found an element that is definitely composite
			set_prime_flag(not_prime);
		}
	}//end for (i)

	check_for_trivial_case();
	return (prime_flag() == prime);
}



//
// input / output
//

template< class T >
void
factorization< T >::read (std::istream &is)
{
	//
	// I/O - format :
	// [ [unit] ]
	// or
	// [ [unit], [base_1, exp_1], [base_2, exp_2], ..., [base_i, exp_i] ]
	// or
	// [ [base_1, exp_1], [base_2, exp_2], ..., [base_i, exp_i] ]
	//

	lidia_size_t exp;
	single_factor< T > base;
	char c;

	kill();

	// read outermost opening '['
	is >> c;
	if (c != '[')
		lidia_error_handler("factorization< T >",
				    "read(std::istream&)::[ expected");

	do {
		is >> c;
		if (c != '[')
			lidia_error_handler("factorization< T >",
					    "read(std::istream&)::[ expected");

		is >> base; // either a unit or a 'normal' component
		is >> c;
		switch (c) {
		case (']') : // read a unit (must be 1st component)
			if (no_of_components() != 0)
				lidia_error_handler("factorization< T >",
						    "read(std::istream&)::unit must be first component");
			*epsilon = base;
			if (epsilon->base() != base.extract_unit())
				lidia_error_handler("factorization< T >",
						    "read(std::istream&)::unit expected");
			break;
		case (',') : // read a 'normal' component
			is >> exp;
			append(base, exp);
			is >> c;
			if (c != ']')
				lidia_error_handler("factorization< T >",
						    "read(std::istream&)::']' expected");
			break;
		default : lidia_error_handler("factorization< T >",
					      "read(std::istream&)::',' or ']' expected");
		}

		is >> c;
	} while (c == ',');

	// read outermost closing ']'
	if (c != ']')
		lidia_error_handler("factorization< T >",
				    "read(std::istream&)::']' expected");

//    check_for_trivial_case(); NICHT NOTWENDIG
}



template< class T >
void
factorization< T >::write (std::ostream &os) const
{
	//
	// I/O - format :
	// [ [unit] ]
	// or
	// [ [unit], [base_1, exp_1], [base_2, exp_2], ..., [base_i, exp_i] ]
	// or
	// [ [base_1, exp_1], [base_2, exp_2], ..., [base_i, exp_i] ]
	//

	lidia_size_t i;
	bool print_comma = false;

	os << "[ ";

	if (!(epsilon->is_one()) || no_of_components() == 0) {
		// if (epsilon != 1)
		os << "[ " << *epsilon << " ]";
		print_comma = true;
	}

	for (i = 0; i < no_of_prime_components(); i++) {
		if (i == 0 && print_comma == false) {
			os << "[";
			print_comma = true;
		}
		else
			os << ",[";
		os << prime_base(i); // base
		os << ",";
		os << prime_exponent(i); // exponent
		os << "]";
	}

	for (i = 0; i < no_of_composite_components(); i++) {
		if (i == 0 && print_comma == false) {
			os << "[";
			print_comma = true;
		}
		else {
			os << ",[";
		}

		os << composite_base(i); // base
		os << ",";
		os << composite_exponent(i); // exponent
		os << "]";
	}

	os << " ]";
}



//
// functions for computing with factorizations
//

template< class T >
void
factorization< T >::invert ()
{
	lidia_size_t i, plen = no_of_prime_components(),
		clen = no_of_composite_components();

	// negate exponents
	for (i = 0; i < plen; i++)
		prime_component[i].right() = - prime_component[i].right();
	for (i = 0; i < clen; i++)
		composite_component[i].right() = - composite_component[i].right();

	// invert unit
	if (!(epsilon->is_one())) {
		single_factor< T > ONE; // default value is '1'
		divide(*epsilon, ONE, *epsilon); // eps = 1/eps
	}
}



template< class T >
void
factorization< T >::square ()
{
	lidia_size_t i, plen = no_of_prime_components(),
		clen = no_of_composite_components();

	// double exponents
	for (i = 0; i < plen; i++)
		prime_component[i].right() <<= 1;
	for (i = 0; i < clen; i++)
		composite_component[i].right() <<= 1;

	// square unit
	if (!(epsilon->is_one()))
		multiply(*epsilon, *epsilon, *epsilon); // eps = eps^2
}



template< class T >
void
factorization< T >::power (lidia_size_t p)
{
	lidia_size_t i, plen = no_of_prime_components(),
		clen = no_of_composite_components();

	if (p == 0) {
		kill();
		return;
	}

	// multiply exponent
	for (i = 0; i < plen; i++)
		prime_component[i].right() *= p;
	for (i = 0; i < clen; i++)
		composite_component[i].right() *= p;

	// multiply unit
	if (!(epsilon->is_one())) {
		single_factor< T > TMP; // default value is '1'
		if (p < 0) {
			divide(TMP, TMP, *epsilon);
			p = -p;
		}
		else
			TMP = *epsilon;

		*epsilon = TMP;
		for (i = 2; i < p; i++)
			multiply(*epsilon, *epsilon, TMP);
		// we could use binary exponentiation here
	}
}



template< class T >
void
factorization< T >::concat (const factorization< T > & a, const factorization< T > & b)
{
	multiply(*epsilon, *a.epsilon, *b.epsilon);
	prime_component.concat(a.prime_component, b.prime_component);
	composite_component.concat(a.composite_component, b.composite_component);

	concat_state(a, b);
	clear_attributes();
	check_for_trivial_case();
}



//
// comparisons
//

template< class T >
bool
factorization< T >::operator == (const factorization< T > & b) const
{
	// *this and b are considered to be equal iff
        // *this and b are equal after normalization
	// which is equivalent to
	// the normalization of *this divided by b equals 1

	if (*epsilon != *b.epsilon)
		return false;

	factorization< T > tmp(*this);
	divide(tmp, tmp, b); // tmp = a/b
	tmp.normalize();

	// if tmp == 1, then no_of_components() == 0
	if (tmp.no_of_components() == 0)
		return true;
	// else
	return false;
}



#if 0
// a possibly faster method
template< class T >
bool
factorization< T >::operator == (const factorization< T > & b) const
{
	if (*epsilon != *b.epsilon)
		return false;

	factorization< T > copy_of_a, copy_of_b;
	const factorization< T > * ptr_a = this;
	const factorization< T > * ptr_b = &b;

	if (!(is_normalized())) {
		copy_of_a.assign(*this);
		copy_of_a.normalize(); //normalize a copy of a
		ptr_a = &copy_of_a;
	}

	if (!(b.is_normalized())) {
		copy_of_b.assign(b);
		copy_of_b.normalize(); //normalize a copy of b
		ptr_b = &copy_of_b;
	}

	//now, *ptr_a and *ptr_b are normalized (and thus sorted)

	// compare length
	lidia_size_t len = ptr_a->component.size();
	if (len != ptr_b->component.size())
		return false;

	// compare components
	register lidia_size_t i;
	for (i = 0; i < len; i++)
		if (ptr_a->component[i] != ptr_b->component[i])
			return false;

	// now, *ptr_a == *ptr_b

	// we might want to store the information we additionally gained
	// by the refinement process, but this would discard 'const'

	return true;
}
#endif   // #if 0



//
// factoring
//

template< class T >
void
factorization< T >::factor_all_components ()
{
	if (prime_flag() == prime)	// already a prime factorization ?
		return;

	lidia_size_t i;
	factorization< T > F;
	for (i = no_of_composite_components()-1; i >= 0; i--) {
		F = composite_base(i).factor();
		F.power(composite_exponent(i));
		// F need not necessarily be a prime factorization

		replace(i, F);
//##### funktioniert, wenn replace hinten anhaengt und auch, wenn
//##### replace an der Stelle i einfuegt
	}

	clear_attributes();
	normalize();
}



//
// high-level functions
//

template< class T >
void
factorization< T >::replace (lidia_size_t pos, const factorization< T > & g)
{
	// deletes composite_component[i] and replaces it by the composite
	// components of 'g'
	// appends the prime components of 'g' to *this.prime_component

	if (pos< 0 || pos >= no_of_composite_components())
		lidia_error_handler("factorization< T >", "replace::wrong index");

	decomp_state fs = prime_flag();
	decomp_state as = composite_base(pos).prime_flag();
	decomp_state gs = g.prime_flag();

	// check state
	if (fs == as) {
		fs = gs;
		if (gs != not_prime) {
			set_prime_flag(unknown);
			if (! is_prime_factorization())
				fs = prime_flag();
		}
	}

	// remove component at position 'pos' from the factorization
	// and append g's components

	lidia_size_t n = g.no_of_composite_components();
	if (n != 0) {
		if (this != &g) {
			composite_component.shift_right(pos, n-1);
			composite_component.assign(pos, g.composite_component, 0, n-1);
		}
		else {
			sort_vector< ppair < single_factor < T >, lidia_size_t > > tmp(composite_component);
			composite_component.shift_right(pos, n-1);
			composite_component.assign(pos, tmp, 0, n-1);
		}
	}
	else
		composite_component.remove_from(pos);

	prime_component.concat(prime_component, g.prime_component);

	LiDIA::multiply(*epsilon, *epsilon, *g.epsilon);

	clear_attributes();
	set_prime_flag(fs);

	check_for_trivial_case();
}



template< class T >
void
factorization< T >::append (const single_factor< T > & a, lidia_size_t exp)
{
	// append single_factor < T > a and lidia_size_t exp to the factorization

	if (exp == 0)
		return;

	if (a.is_one())
	        return;

	lidia_size_t i, nc = no_of_composite_components(), np = no_of_prime_components();
	single_factor< T > *ptr;
	if (!a.is_prime_factor()) {
		composite_component[nc].left() = a;
		composite_component[nc].right() = exp;
		ptr = &composite_component[nc].left();
	}
	else {
		prime_component[np].left() = a;
		prime_component[np].right() = exp;
		ptr = &prime_component[np].left();
	}

	single_factor< T > tmp = ptr->extract_unit();
	if (!tmp.is_one()) {
		if (exp < 0) {
			single_factor< T > ONE; // default value is '1'
			divide(tmp, ONE, tmp);
			exp = -exp;
		}

		for (i = 0; i < exp; i++)
			LiDIA::multiply(*epsilon, *epsilon, tmp);
		// we could use binary exponentiation here
	}

	if (ptr->is_one()) {
		if (a.is_prime_factor())
  		      prime_component.remove_from(np);
		else
		      composite_component.remove_from(nc);
	}
	else
		if (!a.is_prime_factor())  concat_state(*this, a);

	clear_attributes();
	check_for_trivial_case();
}



template< class T >
void
factorization< T >::swapping_append (single_factor< T > & a, lidia_size_t exp)
{
	// append single_factor < T > a and lidia_size_t exp to the factorization

	if (exp == 0)
		return;

//    if (a.is_one()) return;

	lidia_size_t i, nc = no_of_composite_components(), np = no_of_prime_components();
	single_factor< T > *ptr;
	if (!a.is_prime_factor()) {
		LiDIA::swap(composite_component[nc].left(), a);
		composite_component[nc].right() = exp;
		ptr = &composite_component[nc].left();
	}
	else {
		LiDIA::swap(prime_component[np].left(), a);
		prime_component[np].right() = exp;
		ptr = &prime_component[np].left();
	}

	single_factor< T > tmp = ptr->extract_unit();
	if (!tmp.is_one()) {
		if (exp < 0) {
			single_factor< T > ONE; // default value is '1'
			divide(tmp, ONE, tmp);
			exp = -exp;
		}

		for (i = 0; i < exp; i++)
			LiDIA::multiply(*epsilon, *epsilon, tmp);
		// we could use binary exponentiation here
	}
	else
		multiply(*epsilon, *epsilon, tmp);

	if (ptr->is_one()) {
		if (a.is_prime_factor())
			lidia_error_handler("factorization< T >", "append(single_factor< T >, lidia_size_t)::you tried to append a prime factor which only consisted in a unit");

		composite_component.remove_from(nc);
	}
	else
		if (!a.is_prime_factor())  concat_state(*this, a);

	clear_attributes();
	check_for_trivial_case();
}



template< class T >
void
factorization< T >::append (const T & a, lidia_size_t exp)
{
	single_factor< T > aa(a);
	swapping_append(aa, exp);
}



template< class T >
void
factorization< T >::sort ()
{
	// we sort both lists separately

	if (attributes & sorted)		// already sorted ?
		return;

    // sort according to function compare (see ppair.h)
#ifdef HAVE_HP_UX_CC
	prime_component.sort(ppair_sf_lst_compare);
	composite_component.sort(ppair_sf_lst_compare);
#else
	prime_component.sort(compare);
	composite_component.sort(compare);
#endif

	attributes |= sorted; // set flag
}



template< class T >
void
factorization< T >::normalize (sort_vector< ppair < single_factor < T >, lidia_size_t > > &v)
	// sort and collect equal single_factor < T > s in v
{
	v.sort();

	lidia_size_t i = 0, len = v.size() - 1;
	while (i < len) {
		while (i < len ? v[i].left() == v[i+1].left() : 0) {
			// add up exponents
			v[i].right() += v[i+1].right();

			// sum up information
			v[i].left().share_state(v[i+1].left());

			// remove component[i+1] from the factorization
			v.remove_from(i+1);
			len --;
		}

		if (v[i].right() == 0 || v[i].left().is_one()) {
			v.remove_from(i);
			len --;
		}
		else i++;
	}
}



template< class T >
void
factorization< T >::normalize ()	// sort and collect equal single_factor< T > s
{
	if (attributes & normalized)	// already normalized ?
	{
		if (prime_flag() == prime)
			attributes |= refined;
		return;
	}

	lidia_size_t len, i, j;

	sort();

	// STEP 1: normalize prime components
	normalize(prime_component);

	// STEP 2: normalize composite components
	normalize(composite_component);

	attributes |= sorted;


	//  STEP 3: check if any of the prime components occur in the list for
	//  composite components; assumes both lists to be sorted
	len = no_of_prime_components();
	i = j = 0;
	bool ok = true;
	if (len == 0 || no_of_composite_components() == 0)
		ok = false;

	while (ok) {
		if (prime_base(i) == composite_base(j)) {
			if (composite_base(j).prime_flag() == not_prime)
				lidia_error_handler("factorization< T >",
						    "normalize(void)::one component cannot be prime and"
						    "not_prime at the same time");

			// add exponents
			prime_component[i].right() += composite_exponent(j);

			// remove from list of composite components
			composite_component.remove_from(j);

			// if newly gained exponent == 0 : remove
			if (prime_exponent(i) == 0)
				prime_component.remove_from(i);
			else
				i++; // goto next element
		}
		else {
			if (prime_base(i) < composite_base(j))
				i++;
			else	// (prime_base(i) < composite_base(j))
				j++;
		}

		if (i == len || j == no_of_composite_components())
			ok = false;
	}

	check_for_trivial_case();
	if (prime_flag() == unknown) {
		ok = false;
		//  'ok' will be set to true if we find an element with
		//   prime_flag == unknown

		for (i = 0; i < no_of_composite_components(); i++)
			if (composite_base(i).prime_flag() == unknown)
				ok = true;

		if (!ok) {
			// no elements have prime_flag == unknown
			set_prime_flag(not_prime);
		}
	}

	attributes |= normalized; // set flags

	if (prime_flag() == prime)
		attributes |= refined;
}



//
// REFINE
//

template< class T >
void
factorization< T >::refine2 (sort_vector< ppair< single_factor< T >, lidia_size_t > > & v,
			     ppair< single_factor< T >, lidia_size_t > & SF,
			     const ppair< single_factor< T >, lidia_size_t > & a,
			     const ppair< single_factor< T >, lidia_size_t > & b)

	// we assume vector v to be of expanding size !!!

	// compute refinement of {a, b} and store the
	// last single_factor (with exponent) of this factorization
	// individually in SF
	// thus {v[0], ..., v[n], SF} represents the refinement of {a, b}
{
	debug_handler("factorization< T >" , "refine2()");

	lidia_size_t i = 0, last_index = 2;
	ppair< single_factor < T >, lidia_size_t > sf;

	v.set_size(2);
	v[0] = a;
	v[1] = b;

	while (i < last_index-1) {
		gcd(sf.left(), v[i].left(), v[i+1].left());

		if (sf.left().is_one())
			i++;
		else {
			sf.right() = ord_divide(sf.left(), v[i].left()) *   v[i].right();
			sf.right() += ord_divide(sf.left(), v[i+1].left()) * v[i+1].right();

			if (sf.right() != 0) {
				if (v[i+1].left().is_one())
					v[i+1] = sf;
				else {
					v.insert_at(sf, i+1);
					last_index++;
				}
			}
		}
	}

	// extract last element and return it as an individual
	i = v.size()-1;
	SF = v[i];
	v.set_size(i);

	// normalize(v) -- only remove single_factor < T > s which equal '1'
	for (lidia_size_t j = 1; j < i; j++)
		if (v[j].left().is_one()) {
			v.remove_from(j);
			i--; j--;
		}

	//if sorting has advantages
	//v.sort(compare)
}



template< class T >
void
factorization< T >::refine ()			// factor refinement
{
	debug_handler("factorization< T >", "refine()");

	sort_vector< ppair < single_factor < T >, lidia_size_t > >
		v(0, vector_flags(vector_flags::expand)), w(0, vector_flags(vector_flags::expand));
	ppair< single_factor < T >, lidia_size_t > b, c, sf;

	lidia_size_t i, j, k, l = composite_component.size();

	// STEP 1: compute refinement with prime_component
	for (i = 1; i < prime_component.size(); i++)
		for (j = 1; j < composite_component.size(); j++) {
			prime_component[i].right() +=
				ord_divide(prime_component[i].left(), composite_component[j].left())
				* composite_component[j].right();
		}

	// STEP 2: compute refinement on composite_component
	if (l >= 2) {
		factorization< T >::refine2(v, sf, composite_component[0], composite_component[1]); // MM

		k = v.size();
		v[k] = sf;

		for (i = 2; i < l; i++) {
			j = 0;
			b = v[0];
			c = composite_component[i];

			do
			{
				// compute refinement of {b, c}
				factorization< T >::refine2(w, sf, b, c); // MM

				// insert w into v and proceed with {.., sf}
				k = w.size();

				v.shift_right(j+1, k-1);

				v.assign(j, w, 0, k-1);

				j += k;

				if (j < v.size()) {
					b = v[j];
					c = sf;
				}
				else {
					if (!sf.left().is_one()) {
						v[j] = sf;
						j++;
					}
				}
			} while (j < v.size() && !sf.left().is_one());
		}

		LiDIA::swap(v, composite_component);
		attributes &= ~normalized;
	}

	normalize();
	check_for_trivial_case();
	attributes |= refined; // set flag
}



template< class T >
bool
factorization< T >::refine (const single_factor< T > & x)
{
	lidia_size_t i, p_len, c_len, len, e_akt;
	bool rc = false, primefactor = x.is_prime_factor();

	single_factor< T > d;
	ppair< single_factor < T >, lidia_size_t > sf;

	c_len = len = no_of_composite_components();
	p_len = no_of_prime_components();

	for (i = 0; i < len; i++) {
		gcd(d, composite_base(i), x);

		if (!d.is_one() && d != composite_component[i].left()) {
			rc = true;

			e_akt = ord_divide(d, composite_component[i].left());

			sf.left().assign(d);
			sf.left().set_prime_flag(unknown);
			sf.right() = e_akt * composite_exponent(i);

			if (!composite_base(i).is_one())
			{// append sf
				if (primefactor)
					prime_component[p_len] = sf, p_len++;
				else
					composite_component[c_len] = sf, c_len++;
			}
			else {
				if (!primefactor)
					composite_component[i] = sf;
				else {
					prime_component[p_len] = sf;
					p_len++;
					composite_component.remove_from(i);
					len--;
				}
			}
		}//end if d == 1 && ...
	}//end for i

	normalize();

	return rc;
}





// **********************
// NUR ZU TESTZWECKEN !!!
// **********************

template< class T >
void
factorization< T >::output () const
{
	std::cout << "\n--------------------------------------------------------\n";
	std::cout << *this << std::endl;
	lidia_size_t i;

	std::cout << "Einheit :  " << *epsilon << std::endl;

	std::cout << "Primkomponenten :" << std::endl;
	for (i = 0; i < no_of_prime_components(); i++) {
		std::cout << "  " << i << "     (prime) : " << prime_base(i) << " ^ ";
		std::cout << prime_exponent(i) << std::endl;
	}

	std::cout << "Zusammengesetzte Komponenten :" << std::endl;
	for (i = 0; i < no_of_composite_components(); i++) {
		std::cout << "  " << i;
		if (composite_base(i).prime_flag() == not_prime)
			std::cout << " (not_prime) : ";
		else
			std::cout << "   (unknown) : ";
		std::cout << composite_base(i) << " ^ ";
		std::cout << composite_exponent(i) << std::endl;
	}

	std::cout << "Attribute der Faktorisierung :" << std::endl;
	if (prime_flag() == unknown)
		std::cout << "(unknown  ";
	else
		if (prime_flag() == prime)
			std::cout << "(prime    ";
		else
			std::cout << "(not_prime";

	if (is_sorted())
		std::cout << " sorted";
	else
		std::cout << "       ";
	if (is_normalized())
		std::cout << " normalized";
	else
		std::cout << "           ";
	if (is_refined())
		std::cout << " refined)";
	else
		std::cout << ")";
	std::cout << "\n--------------------------------------------------------\n\n";
}



#ifdef LIDIA_NAMESPACE
# ifndef IN_NAMESPACE_LIDIA
}	// end of namespace LiDIA
# endif
#endif



#endif	// LIDIA_FACTORIZATION_CC_GUARD_
