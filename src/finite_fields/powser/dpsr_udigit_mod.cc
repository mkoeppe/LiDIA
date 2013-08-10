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
//	Author	: Thorsten Rottschaefer (TR)
//                Markus Maurer (MM)
//	Changes	: See CVS log
//
//==============================================================================================


#ifdef HAVE_CONFIG_H
# include	"config.h"
#endif
#include	"LiDIA/udigit_mod.h"
#include	"LiDIA/math_vector.h"
#include	"LiDIA/dense_power_series.h"



#include	"LiDIA/finite_fields/base_dense_power_series.h"
#include	"LiDIA/finite_fields/base_dense_power_series.cc"



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



template class base_dense_power_series< udigit_mod >;



//
// ******************************************
// **** dense_power_series <udigit_mod> ****
// ******************************************
//

fft_prime dense_power_series< udigit_mod >::q;



dense_power_series< udigit_mod >::
dense_power_series ()
	: base_dense_power_series< udigit_mod > ()
{
	debug_handler ("dense_power_series< udigit_mod >", "dense_power_series()");
}



dense_power_series< udigit_mod >::
dense_power_series (const udigit_mod & a, lidia_size_t l)
	: base_dense_power_series< udigit_mod > (a, l)
{
	debug_handler ("dense_power_series< udigit_mod >",
		       "dense_power_series(const udigit_mod&, lidia_size_t)");
}



dense_power_series< udigit_mod >::
dense_power_series (const base_vector< udigit_mod > & a, lidia_size_t f)
	: base_dense_power_series< udigit_mod > (a, f)
{
	debug_handler ("dense_power_series< udigit_mod >",
		       "dense_power_series(const base_vector< udigit_mod > &, lidia_size_t)");
}



dense_power_series< udigit_mod >::
dense_power_series (const base_dense_power_series< udigit_mod > & a)
	: base_dense_power_series< udigit_mod > (a)
{
	debug_handler ("dense_power_series< udigit_mod >",
		       "dense_power_series(const base_dense_power_series< udigit_mod > &)");
}



dense_power_series< udigit_mod >::
~dense_power_series ()
{
	debug_handler ("dense_power_series< udigit_mod >",
		       "~dense_power_series()");
}



// ************************************************
// ************ assignment - operator *************
// ************************************************


dense_power_series< udigit_mod > &
dense_power_series< udigit_mod >::
operator = (const dense_power_series< udigit_mod > & a)
{
	debug_handler ("dense_power_series< udigit_mod >",
		       "operator = (const dense_power_series< udigit_mod > &)");

	if (&a != this) {
		base_dense_power_series< udigit_mod >::operator = (a);
	}
	return *this;
}



dense_power_series< udigit_mod > &
dense_power_series< udigit_mod >::
operator = (const sparse_power_series< udigit_mod > & a)
{
	debug_handler ("dense_power_series< udigit_mod >",
		       "operator = (const sparse_power_series< udigit_mod > &)");
	base_dense_power_series< udigit_mod >::operator = (a);
	return *this;
}



// ************************************************
// ********** arithmetic via functions ************
// ************************************************


void
dense_power_series< udigit_mod >::
square (const dense_power_series< udigit_mod > & a)
{
	debug_handler ("dense_power_series< udigit_mod >",
                       "square(dense_power_series< udigit_mod > &, const dense_power_series< udigit_mod > &)");
	lidia_size_t nc;
	lidia_size_t pc;
	lidia_size_t i, j;
	lidia_size_t ed, n;
	lidia_size_t non_zero_index_a;
	int  ident = 0;

	math_vector< udigit_mod > *C;
	math_vector< udigit_mod > *A = a.coeff;

	udigit_mod x;
	udigit_mod tmp;
	udigit_mod zero_elem;

	if (A->size () == 0)
		lidia_error_handler ("dense_power_series< udigit_mod >::multiply"
				     "(dense_power_series< udigit_mod > &, dense_power_series< udigit_mod > &, dense_power_series< udigit_mod > &",
				     "Argument not initialized.");

	if (a.is_zero (non_zero_index_a))
		assign_zero (2 * a.last + 1);

	else {
		// precision and first of c

		pc = A->size() - non_zero_index_a;
		nc = 2 * (a.first + non_zero_index_a);


		// &a == this ?

		if (this != &a) {
			C = coeff;
		}
		else {
			ident = 1;
			C = new math_vector< udigit_mod >;
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

		first = nc;
		last = nc + pc - 1;


		// copy result if necessary

		if (ident) {
			delete coeff;
			coeff = C;
		}

        } // end else a.is_zero (...)
}



void
dense_power_series< udigit_mod >::
multiply_plain (const dense_power_series< udigit_mod > & a ,
		const dense_power_series< udigit_mod > & b)
{
	debug_handler ("dense_power_series< udigit_mod >",
                       "multiply(dense_power_series< udigit_mod > &, dense_power_series< udigit_mod > &, dense_power_series< udigit_mod > &");

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

		math_vector< udigit_mod > *A = a.coeff;
		math_vector< udigit_mod > *B = b.coeff;
		math_vector< udigit_mod > *C;

		udigit_mod zero_elem;
		udigit_mod tmp;

		zero_elem.assign_zero ();

		if (A->size () == 0 || B->size () == 0)
			lidia_error_handler ("dense_power_series< udigit_mod >::multiply"
					     "(dense_power_series< udigit_mod > &, dense_power_series< udigit_mod > &, dense_power_series< udigit_mod > &",
					     "Arguments not initialized.");

		zero_a = a.is_zero (non_zero_index_a);
		zero_b = b.is_zero (non_zero_index_b);

		if (zero_a && zero_b)
			assign_zero (a.last + b.last + 1);

		else if (zero_a)
			assign_zero (a.last + b.first + non_zero_index_b);

		else if (zero_b)
			assign_zero (b.last + a.first + non_zero_index_a);

		else {
			// precision and first of c

			pc = A->size() - non_zero_index_a < B->size() - non_zero_index_b ?
				A->size() - non_zero_index_a : B->size() - non_zero_index_b;

			nc = a.first + b.first + non_zero_index_a + non_zero_index_b;


			// c == a or c == b ?

			if ((this != &a) && (this != &b)) {
				C = coeff;
			}
			else {
				ident = 1;
				C = new math_vector< udigit_mod >;
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

			first = nc;
			last = nc + pc - 1;


			// copy result if necessary

			if (ident) {
				delete coeff;
				coeff = C;
			}

		} // end - else if (zero_a && zero_b)

	} // end - else if (&a == &b)
}



void
dense_power_series< udigit_mod >::
invert (const dense_power_series< udigit_mod > & a)
{
	debug_handler ("dense_power_series< udigit_mod >",
                       "invert(dense_power_series< udigit_mod > &, const dense_power_series< udigit_mod > &)");

	lidia_size_t non_zero_index_a;
	lidia_size_t nc;
	lidia_size_t pc;
	lidia_size_t i, j, k;
	int ident = 0;

	udigit_mod x, y;
	udigit_mod zero_elem;

	math_vector< udigit_mod > *C;
	math_vector< udigit_mod > *A = a.coeff;


	// check for invalid argument

	if (A->size() == 0) {
		lidia_error_handler ("dense_power_series< udigit_mod >::invert(dense_power_series< udigit_mod > &, dense_power_series< udigit_mod > &)",
				     "Argument not initialized.");
        }


	// division by zero ?

	if (a.is_zero (non_zero_index_a)) {
		lidia_error_handler ("dense_power_series< udigit_mod >::invert(dense_power_series< udigit_mod > &, dense_power_series< udigit_mod > &)",
				     "Division by zero.");
        }
	else {
		// precision and first of this

		pc = A->size() - non_zero_index_a;
		nc = - (a.first + non_zero_index_a);


		// &a == this ?

		if (this != &a)
			C = coeff;
		else {
			ident = 1;
			C = new math_vector< udigit_mod >;
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

		first = nc;
		last = nc + pc - 1;


		// copy C if necessary

		if (ident) {
			delete coeff;
			coeff = C;
		}

        } // end else a.is_zero (...)
}



void
dense_power_series< udigit_mod >::
power (const dense_power_series< udigit_mod > & a ,
       long n)
{
	debug_handler ("dense_power_series< udigit_mod >",
                       "power(dense_power_series< udigit_mod > &, const dense_power_series< udigit_mod > &, long)");

	dense_power_series< udigit_mod > z;
	lidia_size_t non_zero_index_a;
	bool zero_a;


	// check for invalid argument

	if ((a.coeff)->size () == 0) {
		lidia_error_handler ("dense_power_series< udigit_mod >::power"
				     "(dense_power_series< udigit_mod > &, dense_power_series< udigit_mod > &, lidia_size_t)",
				     "Argument not initialized.");
        }

	zero_a = a.is_zero (non_zero_index_a);

	if (n == 0) {
		if (zero_a)
			assign_one (0);
		else
			assign_one (a.last - (a.first + non_zero_index_a));
	}
	else if (zero_a)
		assign_zero (static_cast<lidia_size_t>(n * a.last));
	else {
		// initialize z which holds the squares

		if (n < 0) {
			z.invert(a);
			n = -n;
		}
		else
			z = a;

		assign_one ((a.coeff)->size() - 1 - non_zero_index_a);


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



void
dense_power_series< udigit_mod >::
divide (const dense_power_series< udigit_mod > & a,
	const dense_power_series< udigit_mod > & b)
{
	debug_handler ("dense_power_series< udigit_mod >",
                       "divide(dense_power_series< udigit_mod > &, dense_power_series< udigit_mod > &, dense_power_series< udigit_mod > &");

	dense_power_series< udigit_mod > inv_b;
	inv_b.invert(b);
	multiply(a, inv_b);
}



void
dense_power_series< udigit_mod >::
divide (const udigit_mod  & b,
	const dense_power_series< udigit_mod > & a)
{
	debug_handler ("dense_power_series< udigit_mod >",
                       "divide(dense_power_series< udigit_mod > &, T&, dense_power_series< udigit_mod > &)");

	dense_power_series< udigit_mod > d;
	d.invert(a);
	base_dense_power_series< udigit_mod >::multiply(b, d);
}



void
dense_power_series< udigit_mod >::
multiply (const dense_power_series< udigit_mod > & a,
	  const dense_power_series< udigit_mod > & b)
{
	lidia_size_t degree_a, degree_b;

	// get the degree of a and b
	degree_a = (a.coeff)->size() -1;
	degree_b = (b.coeff)->size() -1;

	if (degree_a <= fft_cross_over_point || degree_b <= fft_cross_over_point)
		if (&a == &b)
			square(a);
		else
			multiply_plain(a, b);

	else
		multiply_fft(a, b);
}



void
dense_power_series< udigit_mod >::
multiply_fft (const dense_power_series< udigit_mod > & a,
	      const dense_power_series< udigit_mod > & b)
{
	lidia_size_t degree_a, degree_b, k, d, i;

	// initialize q with the modulus of udigit_mod (is the same for a and b)
	q.set_prime(udigit_mod::get_modulus());

	// get the degree of a and b
	degree_a = (a.coeff)->size() -1;
	degree_b = (b.coeff)->size() -1;

	// evaluate degree of a*b
	d = degree_a + degree_b + 1;

	// evaluate for which k 2^k >= d
	// (means the next power of 2 bigger than d)
	i = 1; k = 0;
	while (i < d)
        {
		i = i << 1;
		k++;
        };

	// give this vector new dimensions
	coeff->set_capacity (i);
	coeff->set_size (i);
	//std::cout << i;

	// call the friend function multiply_fft of class fft_prime with a and b.
	if (!LiDIA::multiply_fft(static_cast<void*>(coeff->get_data_address()), 1,
				    static_cast<void*>((a.coeff)->get_data_address()), degree_a,
				    static_cast<void*>((b.coeff)->get_data_address()), degree_b, q))
		lidia_error_handler("dpsr_udigit_mod::multiply_fft", "error in multiplying a and b");

	//if (!::multiply_fft(coeff->get_data_address(),
	//		(a.coeff)->get_data_address(), degree_a,
	//		(b.coeff)->get_data_address(), degree_b, q ) )
	//lidia_error_handler("dpsr_udigit_mod::multiply_fft","error in multiplying a and b");

	// cut the result-vector to the minimum of the length of a and b
	if (degree_a <= degree_b)
		coeff->set_capacity (degree_a + 1);
	else
		coeff->set_capacity (degree_b + 1);

	//set new beginning exponent of result-vector
	first = a.get_first() + b.get_first();
	last = first + coeff->get_capacity ();
}



lidia_size_t dense_power_series< udigit_mod >::fft_cross_over_point = 30;



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
