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


#ifndef LIDIA_DENSE_POWER_SERIES_CC_GUARD_
#define LIDIA_DENSE_POWER_SERIES_CC_GUARD_


#ifndef LIDIA_DENSE_POWER_SERIES_H_GUARD_
# include	"LiDIA/dense_power_series.h"
#endif



#ifdef LIDIA_NAMESPACE
# ifndef IN_NAMESPACE_LIDIA
namespace LiDIA {
# endif
#endif



//
// ***** arithmetic via functions *****
//

template< class T >
void
dense_power_series< T >::
square (const dense_power_series< T > & a)
{
	debug_handler ("dense_power_series< T >",
                       "square(dense_power_series< T > &, const dense_power_series< T > &)");

	lidia_size_t nc;
	lidia_size_t pc;
	lidia_size_t i, j;
	lidia_size_t ed, n;
	lidia_size_t non_zero_index_a;
	int  ident = 0;

	math_vector< T > *C;
	math_vector< T > *A = a.coeff;

	T x;
	T tmp;
	T zero_elem;

	if (A->size () == 0)
		lidia_error_handler ("dense_power_series< T >::multiply"
				     "(dense_power_series< T > &, dense_power_series< T > &, dense_power_series< T > &",
				     "Argument not initialized.");

	if (a.is_zero (non_zero_index_a))
		this->assign_zero (2 * a.last + 1);

	else {
		// precision and first of c

		pc = A->size() - non_zero_index_a;
		nc = 2 * (a.first + non_zero_index_a);


		// &a == this ?

		if (this != &a) {
			C = this->coeff;
		}
		else {
			ident = 1;
			C = new math_vector< T >;
		}

		C->set_capacity (pc);


		// square a

		zero_elem.assign_zero ();

		for (i = 0, n = 2 * non_zero_index_a; i < pc; i++, n++) {
			if (n & 1)
				ed = (n-1) >> 1;
			else
				ed = (n >> 1) - 1;

			(*C)[i] = zero_elem;

			for (j = non_zero_index_a; j <= ed; j++) {
				LiDIA::multiply (tmp, (*A)[j], (*A)[n-j]);
				LiDIA::add      ((*C)[i], (*C)[i], tmp);
			}

			LiDIA::add ((*C)[i], (*C)[i], (*C)[i]);

			if (!(n&1)) {
				LiDIA::square (x, (*A)[n >> 1]);
				LiDIA::add    ((*C)[i], (*C)[i], x);
			}
		}

		this->first = nc;
		this->last = nc + pc - 1;


		// copy result if necessary

		if (ident) {
			delete this->coeff;
			this->coeff = C;
		}

        } // end else a.is_zero (...)
}



template< class T >
void
dense_power_series< T >::
multiply (const dense_power_series< T > & a  ,
	  const dense_power_series< T > & b)
{
	debug_handler ("dense_power_series< T >",
                       "multiply(dense_power_series< T > &, dense_power_series< T > &, dense_power_series< T > &");

	if (&a == &b) {
		this->square(a);
        }
	else {
		lidia_size_t nc, pc;
		lidia_size_t i, j, n;
		lidia_size_t non_zero_index_a;
		lidia_size_t non_zero_index_b;
		int  zero_a;
		int  zero_b;
		int  ident = 0;

		math_vector< T > *A = a.coeff;
		math_vector< T > *B = b.coeff;
		math_vector< T > *C;

		T zero_elem;
		T tmp;

		zero_elem.assign_zero ();

		if (A->size () == 0 || B->size () == 0)
			lidia_error_handler ("dense_power_series< T >::multiply"
					     "(dense_power_series< T > &, dense_power_series< T > &, dense_power_series< T > &",
					     "Arguments not initialized.");

		zero_a = a.is_zero (non_zero_index_a);
		zero_b = b.is_zero (non_zero_index_b);

		if (zero_a && zero_b)
			this->assign_zero (a.last + b.last + 1);

		else if (zero_a)
			this->assign_zero (a.last + b.first + non_zero_index_b);

		else if (zero_b)
			this->assign_zero (b.last + a.first + non_zero_index_a);

		else {
			// precision and first of c

			pc = A->size() - non_zero_index_a< B->size() - non_zero_index_b ?
				A->size() - non_zero_index_a : B->size() - non_zero_index_b;

			nc = a.first + b.first + non_zero_index_a + non_zero_index_b;


			// c == a or c == b ?

			if ((this != &a) && (this != &b)) {
				C = this->coeff;
			}
			else {
				ident = 1;
				C = new math_vector< T >;
			}

			C->set_capacity (pc);


			// multiply a and b

			for (i = 0, n = non_zero_index_a + non_zero_index_b; i < pc; i++, n++) {
				(*C)[i] = zero_elem;

				for (j = non_zero_index_a; j <= n-non_zero_index_b; j++) {
					LiDIA::multiply (tmp, (*A)[j], (*B)[n-j]);
					LiDIA::add      ((*C)[i], (*C)[i], tmp);
				}
			}

			this->first = nc;
			this->last = nc + pc - 1;


			// copy result if necessary

			if (ident) {
				delete this->coeff;
				this->coeff = C;
			}

		} // end - else if (zero_a && zero_b)

	} // end - else if (&a == &b)
}



template< class T >
void
dense_power_series< T >::
invert (const dense_power_series< T > & a)
{
	debug_handler ("dense_power_series< T >",
                       "invert(dense_power_series< T > &, const dense_power_series< T > &)");

	lidia_size_t non_zero_index_a;
	lidia_size_t nc;
	lidia_size_t pc;
	lidia_size_t i, j, k;
	int ident = 0;

	T x, y;
	T zero_elem;

	math_vector< T > *C;
	math_vector< T > *A = a.coeff;


	// check for invalid argument

	if (A->size() == 0) {
		lidia_error_handler ("dense_power_series< T >::invert(dense_power_series< T > &, dense_power_series< T > &)",
				     "Argument not initialized.");
        }


	// division by zero ?

	if (a.is_zero (non_zero_index_a)) {
		lidia_error_handler ("dense_power_series< T >::invert(dense_power_series< T > &, dense_power_series< T > &)",
				     "Division by zero.");
        }
	else {
		// precision and first of this

		pc = A->size() - non_zero_index_a;
		nc = - (a.first + non_zero_index_a);


		// &a == this ?

		if (this != &a)
			C = this->coeff;
		else {
			ident = 1;
			C = new math_vector< T >;
		}

		C->set_capacity (pc);


		// invert a

		zero_elem.assign_zero ();

		LiDIA::invert ((*C)[0], (*A)[non_zero_index_a]);
		LiDIA::negate (y    , (*C)[0]);

		for (i = 1; i < pc; i++) {
			(*C)[i] = zero_elem;

			for (j = 1, k = 1 + non_zero_index_a; j <= i; j++, k++) {
				LiDIA::multiply (x, (*A)[k], (*C)[i-j]);
				LiDIA::add      ((*C)[i], (*C)[i], x);
			}

			LiDIA::multiply ((*C)[i], (*C)[i], y);
		}

		this->first = nc;
		this->last = nc + pc - 1;


		// copy C if necessary

		if (ident) {
			delete this->coeff;
			this->coeff = C;
		}

        } // end else a.is_zero (...)
}



template< class T >
void
dense_power_series< T >::
power (const dense_power_series< T > & a,
       long n)
{
	debug_handler ("dense_power_series< T >",
                       "power(dense_power_series< T > &, const dense_power_series< T > &, long)");

	dense_power_series< T > z;
	lidia_size_t non_zero_index_a;
	bool zero_a;


	// check for invalid argument

	if ((a.coeff)->size () == 0) {
		lidia_error_handler ("dense_power_series< T >::power"
				     "(dense_power_series< T > &, dense_power_series< T > &, lidia_size_t)",
				     "Argument not initialized.");
        }

	zero_a = a.is_zero (non_zero_index_a);

	if (n == 0) {
		if (zero_a)
			this->assign_one (0);
		else
			this->assign_one (a.last - (a.first + non_zero_index_a));
	}
	else if (zero_a)
		this->assign_zero (static_cast<lidia_size_t>(n * a.last));
	else {
		// initialize z which holds the squares

		if (n < 0) {
			z.invert(a);
			n = -n;
		}
		else
			z = a;

		this->assign_one ((a.coeff)->size() - 1 - non_zero_index_a);


		// repeated squaring

		while (n > 1) {
			// n odd

			if (n&1)
				this->multiply(*this, z);

			z.square(z);

			// divide n by 2

			n = n >> 1;
		}

		if (n == 1)
			this->multiply(*this, z);

	}
}



template< class T >
void
dense_power_series< T >::
divide (const dense_power_series< T > & a ,
	const dense_power_series< T > & b)
{
	debug_handler ("dense_power_series< T >",
                       "divide(dense_power_series< T > &, dense_power_series< T > &, dense_power_series< T > &");

	dense_power_series< T > inv_b;
	inv_b.invert(b);
	multiply(a, inv_b);
}



template< class T >
void
dense_power_series< T >::
divide (const T & b,
	const dense_power_series< T > & a)
{
	debug_handler ("dense_power_series< T >",
                       "divide(dense_power_series< T > &, T&, dense_power_series< T > &)");

	dense_power_series< T > d;
	d.invert(a);
	base_dense_power_series< T >::multiply(b, d);
}



#ifdef LIDIA_NAMESPACE
# ifndef IN_NAMESPACE_LIDIA
}	// end of namespace LiDIA
# endif
#endif



#endif	// LIDIA_DENSE_POWER_SERIES_CC_GUARD_
