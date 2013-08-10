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


#ifndef LIDIA_SPARSE_POWER_SERIES_CC_GUARD_
#define LIDIA_SPARSE_POWER_SERIES_CC_GUARD_


#ifndef LIDIA_SPARSE_POWER_SERIES_H_GUARD_
# include	"LiDIA/sparse_power_series.h"
#endif



#ifdef LIDIA_NAMESPACE
# ifndef IN_NAMESPACE_LIDIA
namespace LiDIA {
# endif
#endif



//
// ***** arithmetical procedures *****
//

template< class T >
void
sparse_power_series< T >::multiply (const sparse_power_series< T > & a ,
				    const sparse_power_series< T > & b)
{
	debug_handler ("sparse_power_series< T >",
		       "multiply (sp_pow & , const sp_pow & , const sp_pow &)");

	a.init_test ("multiply (3 x sp_pow)");
	b.init_test ("multiply (3 x sp_pow)");

	a.sort_test();
	b.sort_test();

	lidia_size_t cf, cl, i, az, bz, ix, tmp_e;
	lidia_size_t max_size;
	lidia_size_t  all, asz, bsz;
	T zero_T, tmp_c;

	sort_vector< spc< T > > *ac, *bc;
	base_vector< T > cc;

	cf = a.first + b.first;
	cl = comparator< lidia_size_t >::min ((a.last + b.first) , (a.first + b.last));

	zero_T.assign_zero ();

	if (a.is_zero() || b.is_zero()) {
		this->assign_zero (cl);
	}
	else {
		// both series contain non-zero coefficients
		max_size = cl - cf + 1;
		asz = a.coeff->size();
		bsz = b.coeff->size();

		ac = a.coeff;
		bc = b.coeff;

		cc.set_capacity (max_size); // tmp. storage for coeff.
		for (i = 0; i < max_size; i ++)         // init. with 0
			cc[i] = zero_T;

		for (az = 0; az < asz; az ++) {
			for (bz = 0; bz < bsz; bz ++) {
				tmp_e = (*ac)[az].exp + (*bc)[bz].exp;

				if (tmp_e <= cl) {
					// coeff will be valid in result
					LiDIA::multiply (tmp_c , (*ac)[az].coeff , (*bc)[bz].coeff);
					LiDIA::add      (cc[tmp_e-cf], cc[tmp_e-cf], tmp_c);
				}
			}
		}

		// --- counting number of non-zero elts. in product ---

		for (i = 0, all = 0; i <= (cl-cf); i++)
			if (cc[i] != zero_T)  all++;

		if (all > 0) {
			this->coeff->set_capacity (all);
			this->coeff->set_size     (all);

			for (i = 0, ix = 0; i <= (cl-cf); i++) {
				if (cc[i] != zero_T) {
					(*this->coeff)[ix].coeff = cc[i];
					(*this->coeff)[ix].exp = i + cf;

					ix ++;
				}
			}

			this->first = (*this->coeff)[0].exp;
			this->last = cl;
		}
	}
}



template< class T >
void
sparse_power_series< T >::invert (const sparse_power_series< T > & a)
{
	debug_handler ("sparse_power_series< T >",
		       "invert (sp_pow & , const sp_pow &)");

	a.init_test ("invert (2 x sp_pow)");
	a.sort_test();

	lidia_size_t i, j;
	lidia_size_t cf, cl, cx;
	lidia_size_t all, aalloc;

	T tmp, ma0, zero_T, tb;

	sparse_power_series< T > b;
	sort_vector< spc< T > > *ac, *bc;

	zero_T.assign_zero ();

	if (! a.is_zero()) {
		cf = - a.first;
		cl = a.last - 2 * a.first;

		all = (cl - cf + 1);
		aalloc = a.coeff->size();

		b.coeff->set_capacity (all);
		b.first = cf;
		b.last = cl;

		ac = a.coeff;
		bc = b.coeff;

		for (i = 0; i < all; i++) (*bc)[i].coeff = zero_T;

		LiDIA::divide ((*bc)[0].coeff, T(1), (*ac)[0].coeff);
		(*bc)[0].exp = cf;

		LiDIA::negate (ma0, (*bc)[0].coeff); // ma0 == -1/a_0

		tb = (*bc)[0].coeff;

		for (i = cf + 1, all = 1; i <= cl; i++) {
			// a * (b = 1/a) is the one-approximation = 1 + \SUM {i = 0} {cl-cf} {0 * X^i}.
			// The coefficient of this serie
			// that corresponds to X^(n) equals \SUM {a.first <= i <= n+a.first, j = n-i} {a_i * b_j}.
			// This sum is used to compute the b with the greatest index that appears in
			// the sum; this is b with index i+j-a.first = i+j+cf which is stored in
			// bc[i+j+cf-cf] = bc[i+j]. Therefore, the sum a[j].exp + (i-1) of the exponents of a_(a[j].exp)
			// and b_(i-1) gives us the index of the coefficient in bc which can be updated by
			// the product of a_(a[j].exp) and b_(i-1).

			// tb is the coefficient of b = 1/a with exponent 'i-1'.

			if (tb != zero_T)

				for (j = 1; (j < aalloc) && ((cx = i-1+(*ac)[j].exp) <= cl-cf); j++) {
					LiDIA::multiply (tmp , tb , (*ac)[j].coeff);
					LiDIA::add      ((*bc)[ cx ].coeff , (*bc)[ cx ].coeff, tmp);
				}

			// std::cout t << i << " : " << (*bc) << std::endl;

			if ((*bc)[i-cf].coeff != zero_T) {
				LiDIA::multiply ((*bc)[all].coeff, (*bc)[i-cf].coeff, ma0);
				(*bc)[all].exp = i;
				tb = (*bc)[all].coeff;

				all ++;
			}
			else {
				tb = zero_T;
			}
		}

		b.coeff->set_capacity (all);
		this->swap(b);
	}
	else {
		lidia_error_handler ("sparse_power_series< T >" ,
				     "inverting of zero failed");
	}
}



template< class T >
void
sparse_power_series< T >::square (const sparse_power_series< T > & a)
{
	debug_handler ("sparse_power_series< T >",
		       "square (sp_pow & , const sp_pow &)");

	a.init_test ("square (sp_pow & , const sp_pow &)");
	a.sort_test ();

	this->multiply(a, a);
}



template< class T >
void
sparse_power_series< T >::power (const sparse_power_series< T > & a,
				 long n)
{
	debug_handler ("sparse_power_series< T >",
		       "power (sp_pow &, const sp_pow &, long)");

	a.init_test ("power (sp_pow & , const sp_pow & , lidia_size_t)");
	a.sort_test ();

	sparse_power_series< T > z;

	if (n >= 0) {
		z = a;
	}
	else {
		z.invert(a);
		n = -n;
	}

	set (T(1), a.last - a.first); // res = 1 * X^0

	while (n > 1) {
		if (n & 1) {
			this->multiply(*this, z);
		}

		z.square(z);

		n >>= 1;
	}

	if (n == 1)
		this->multiply(*this, z);
}



template< class T >
void
sparse_power_series< T >::divide (const sparse_power_series< T > & a,
				  const sparse_power_series< T > & b)
{
	debug_handler ("sparse_power_series< T >",
		       "divide (sp_pow &, const sp_pow &, const sp_pow &)");

	a.init_test ("divide(3 x sp_pow)");
	b.init_test ("divide(3 x sp_pow)");
	a.sort_test ();
	b.sort_test ();

	sparse_power_series< T > invb;
	invb.invert(b);
	this->multiply(a, invb);
}



template< class T >
void
sparse_power_series< T >::divide (const T & b,
				  const sparse_power_series< T > & a)
{
	debug_handler ("sparse_power_series< T >",
		       "divide (sp_pow & , const T & , const sp_pow &)");

	a.init_test ("divide (sp_pow , cont T & , sp_pow)");
	this->invert(a);
	base_sparse_power_series< T >::multiply(*this, b);
}



#ifdef LIDIA_NAMESPACE
# ifndef IN_NAMESPACE_LIDIA
}	// end of namespace LiDIA
# endif
#endif



#endif	// LIDIA_SPARSE_POWER_SERIES_CC_GUARD_
