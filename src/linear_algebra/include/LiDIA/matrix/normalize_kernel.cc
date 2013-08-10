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
//	Author	: Patrick Theobald (PT)
//	Changes	: See CVS log
//
//==============================================================================================


#ifndef LIDIA_NORMALIZE_KERNEL_CC_GUARD_
#define LIDIA_NORMALIZE_KERNEL_CC_GUARD_


#include	<fstream>
#ifndef LIDIA_INFO_H_GUARD_
# include	"LiDIA/info.h"
#endif
#ifndef LIDIA_NORMALIZE_KERNEL_H_GUARD_
# include	"LiDIA/matrix/normalize_kernel.h"
#endif
#include	<cstdlib>


#ifdef LIDIA_NAMESPACE
# ifndef IN_NAMESPACE_LIDIA
namespace LiDIA {
# endif
#endif



#ifdef LIDIA_NAMESPACE
using std::abs;
#endif



#define DMESSAGE "normalize_kernel"


////////////////////////////////////////////////////////////////////////////////////////
// normalization
////////////////////////////////////////////////////////////////////////////////////////

//
// normalize Standard
//

template< class T, class REP, class REP1 >
lidia_size_t *
normalization_kernel< T, REP, REP1 >::normalize_Std (MR< T > &A,
						     lidia_size_t startr,
						     lidia_size_t startc) const
{
	lidia_size_t i, j, k = 1;
	T q, r;

	lidia_size_t *RET = new lidia_size_t[A.rows];

#ifdef LIDIA_COMPUTE_RUNTIME_INFO
	std::ofstream profile("normalize_max.profile");
	T MAX = 0, GMAX = 0;
#endif

	for (i = A.columns - 1, j = A.rows - 1; i >= startc && j >= startr; i--, j--) {
#ifdef LIDIA_COMPUTE_RUNTIME_INFO
		modul.max(A, MAX);
		if (GMAX < MAX)
			GMAX = MAX;
#endif
		if (modul.member(A, j, i) != A.Zero) {
			for (lidia_size_t l = i + 1; l < A.columns; l++)
				if (modul.member(A, j, l) != A.Zero) {
					pos_div_rem(q, r, modul.member(A, j, l), modul.member(A, j, i));
					modul.subtract_multiple_of_column(A, l, q, i, A.rows - 1);
				}
		}
		else {
			RET[k] = j;
			k++;
		}
	}
	if (k == 1) {
		delete [] RET;
		RET = NULL;
	}
	else
		RET[0] = k - 1;

#ifdef LIDIA_COMPUTE_RUNTIME_INFO
	profile << GMAX << std::endl;
	profile.close();
#endif

	return RET;
}



template< class T, class REP, class REP1 >
lidia_size_t *
normalization_kernel< T, REP, REP1 >::normalize_Std (MR< T > &A,
						     matrix< bigint > &TR,
						     lidia_size_t startr,
						     lidia_size_t startc) const
{
	lidia_size_t i, j, k = 1;
	T q, r;

	lidia_size_t *RET = new lidia_size_t[A.rows];

	for (i = A.columns - 1, j = A.rows - 1; i >= startc && j >= startr; i--, j--)
		if (modul.member(A, j, i) != A.Zero) {
			for (lidia_size_t l = i + 1; l < A.columns; l++)
				if (modul.member(A, j, l) != A.Zero) {
					pos_div_rem(q, r, modul.member(A, j, l), modul.member(A, j, i));
					modul.subtract_multiple_of_column(A, l, q, i, A.rows - 1);

					for (lidia_size_t m = 0; m < A.columns; m++)
						TR.sto(m, l, TR.member(m, l) - q * TR.member(m, i));
				}
		}
		else {
			RET[k] = j;
			k++;
		}

	if (k == 1) {
		delete [] RET;
		RET = NULL;
	}
	else
		RET[0] = k - 1;
	return RET;
}



//
// normalize_ChouCollins
//

template< class T, class REP, class REP1 >
lidia_size_t *
normalization_kernel< T, REP, REP1 >::normalize_ChouCollins (MR< T > &A,
							     lidia_size_t startr,
							     lidia_size_t startc) const
{
	lidia_size_t start = startc, i, k, p, l;
	T q, r;

	lidia_size_t *RET = new lidia_size_t[A.rows];

#ifdef LIDIA_COMPUTE_RUNTIME_INFO
	std::ofstream profile("normalize_max.profile");
	T MAX = 0, GMAX = 0;
#endif

	for (i = start + 1, p = startr; i < A.columns; i++, p++) {
//      std::cout << "i = " << i << " p = " << p << std::endl;
		for (k = i - 1, l = p; k >= start; k--, l--) {
//          std::cout << "\t k = " << k << " l = " << l << std::endl;
#ifdef LIDIA_COMPUTE_RUNTIME_INFO
			modul.max(A, MAX);
			if (GMAX < MAX)
				GMAX = MAX;
#endif

			if (modul.member(A, l, i) != A.Zero && modul.member(A, l, k) != A.Zero) {
				pos_div_rem(q, r, modul.member(A, l, i), modul.member(A, l, k));
				modul.subtract_multiple_of_column(A, i, q, k, A.rows - 1);
			}
		}
	}

	k = 1;
	for (i = startr; i < A.rows; i++) {
		if (modul.member(A, i, start + i - startr) == A.Zero) {
			RET[k] = i;
			k++;
		}
	}

	RET[0] = k - 1;

	if (RET[0] == 0) {
		delete[] RET;
		RET = NULL;
	}

#ifdef LIDIA_COMPUTE_RUNTIME_INFO
	profile << GMAX << std::endl;
	profile.close();
#endif

	return RET;
}



template< class T, class REP, class REP1 >
lidia_size_t *
normalization_kernel< T, REP, REP1 >::normalize_ChouCollins (MR< T > &A,
							     matrix< bigint > &TR,
							     lidia_size_t startr,
							     lidia_size_t startc) const
{
	lidia_size_t start = startc, i, k, p, l;
	T q, r;

	lidia_size_t *RET = new lidia_size_t[A.rows];

	for (i = start + 1, p = startr; i < A.columns; i++, p++) {
		for (k = i - 1, l = p; k >= start; k--, l--) {
			if (modul.member(A, l, i) != A.Zero && modul.member(A, l, k) != A.Zero) {
				pos_div_rem(q, r, modul.member(A, l, i), modul.member(A, l, k));
				modul.subtract_multiple_of_column(A, i, q, k, A.rows - 1);

				for (lidia_size_t m = 0; m < A.columns; m++)
					TR.sto(m, i, TR.member(m, i) - q * TR.member(m, k));
			}
		}
	}

	k = 1;
	for (i = startr; i < A.rows; i++)
		if (modul.member(A, i, start + i - startr) == A.Zero) {
			RET[k] = i;
			k++;
		}
	RET[0] = k - 1;

	if (RET[0] == 0) {
		delete[] RET;
		RET = NULL;
	}

	return RET;
}



template< class T, class REP, class REP1 >
lidia_size_t *
normalization_kernel< T, REP, REP1 >::normalize_ChouCollins (MR< T > &A,
							     trans_matrix &TR,
							     matrix< bigint > & tran,
							     lidia_size_t startr,
							     lidia_size_t startc) const
{
	lidia_size_t start = startc, i, k, p, l;
	T q, r;
	matrix< bigint > otran;

	tran.reset();
	otran = tran;
	tran.resize(A.columns, A.columns);
	tran.diag(1, 0);

	lidia_size_t *RET = new lidia_size_t[A.rows];

	for (i = start + 1, p = startr; i < A.columns; i++, p++) {
		for (k = i - 1, l = p; k >= start; k--, l--) {
			if (modul.member(A, l, i) != A.Zero && modul.member(A, l, k) != A.Zero) {
				pos_div_rem(q, r, modul.member(A, l, i), modul.member(A, l, k));
				modul.subtract_multiple_of_column(A, i, q, k, A.rows - 1);

				for (lidia_size_t m = 0; m < A.columns; m++)
					tran.sto(m, i, tran.member(m, i) - q * tran.member(m, k));
			}
		}

		if (i == (A.rows/2)) {
			TR.store_matrix(tran);

			tran = otran;
			tran.resize(A.columns, A.columns);
			tran.diag(1, 0);
		}
	}

	k = 1;
	for (i = startr; i < A.rows; i++)
		if (modul.member(A, i, start + i - startr) == A.Zero) {
			RET[k] = i;
			k++;
		}
	RET[0] = k - 1;

	if (RET[0] == 0) {
		delete[] RET;
		RET = NULL;
	}

	return RET;
}



//
// normalize_ChouCollins_extended
//

template< class T, class REP, class REP1 >
lidia_size_t *
normalization_kernel< T, REP, REP1 >::normalize_ChouCollins_extended (MR< T > &A,
								      lidia_size_t startr,
								      lidia_size_t startc) const
{
	lidia_size_t start = startc, i, k, p, l;
	T q, r;

	lidia_size_t *RET = new lidia_size_t[A.rows];

	for (i = start + 1, p = startr; i < A.columns; i++, p++)
		for (k = i - 1, l = p; k >= start; k--, l--)
			if (modul.member(A, l, i) != A.Zero && modul.member(A, l, k) != A.Zero) {
				pos_div_rem(q, r, modul.member(A, l, i), modul.member(A, l, k));
				modul.subtract_multiple_of_column(A, i, q, k, A.rows - 1);
			}

	// correct diagonal
	for (i = A.columns - 1, l = A.rows - 1; i > start; i--, l--)
		if (modul.member(A, l, i) == A.Zero) {
			k = l;
			while (modul.member(A, k, i) == A.Zero && k > startr)
				k--;
			if (k > startr) {
				if (modul.member(A, k, start + k - startr) == A.Zero) {
					modul.swap_columns(A, i, start + k - startr);
					i++;
					l++;
				}
			}
		}

	k = 1;
	for (i = startr; i < A.rows; i++)
		if (modul.member(A, i, start + i - startr) == A.Zero) {
			RET[k] = i;
			k++;
		}
	RET[0] = k - 1;
	return RET;
}



//
// normalizeHybrid_Std
//

template< class T, class REP, class REP1 >
lidia_size_t *
normalization_kernel< T, REP, REP1 >::normalizeHybrid_Std (MR< T > &A,
							   lidia_size_t startr,
							   lidia_size_t startc) const
{
	lidia_size_t start = startc, i, k;
	T q, r;

#ifdef LIDIA_COMPUTE_RUNTIME_INFO
	std::ofstream profile("normalize_max.profile");
	T MAX = 0, GMAX = 0;
#endif

	// Phase 1
	for (i = startr; i < A.rows; i++) {

#ifdef LIDIA_COMPUTE_RUNTIME_INFO
		modul.max(A, MAX);
		if (GMAX < MAX)
			GMAX = MAX;
#endif

		if (modul.member(A, i, start) == 1)
			for (lidia_size_t l = start + 1; l < A.columns; l++)
				if (modul.member(A, i, l) != A.Zero)
					modul.subtract_multiple_of_column(A, l, modul.member(A, i, l), start, A.rows - 1);
		++start;
	}

	// Phase 2
	lidia_size_t *RET = new lidia_size_t[A.rows];
	k = 1;

	start--;
	for (--i; i >= startr; i--) {

#ifdef LIDIA_COMPUTE_RUNTIME_INFO
		modul.max(A, MAX);
		if (GMAX < MAX)
			GMAX = MAX;
#endif

		if (modul.member(A, i, start) == A.Zero)
			RET[k++] = i;
		else {
			if (modul.member(A, i, start) != 1)
				for (lidia_size_t l = start + 1; l < A.columns; ++l) {
					if (modul.member(A, i, l) != A.Zero) {
						pos_div_rem(q, r, modul.member(A, i, l), modul.member(A, i, start));
						modul.subtract_multiple_of_column(A, l, q, start, A.rows - 1);
					}
				}
		}
		start--;
	}
	if (k == 1) {
		delete [] RET;
		RET = NULL;
	}
	else
		RET[0] = k - 1;

#ifdef LIDIA_COMPUTE_RUNTIME_INFO
	profile << GMAX << std::endl;
	profile.close();
#endif

	return RET;
}



template< class T, class REP, class REP1 >
lidia_size_t *
normalization_kernel< T, REP, REP1 >::normalizeHybrid_Std (MR< T > &A,
							   matrix< bigint > &TR,
							   lidia_size_t startr,
							   lidia_size_t startc) const
{
	lidia_size_t start = startc, i, k;
	T q, r;

	// Phase 1
	for (i = startr; i < A.rows; ++i) {
		if (modul.member(A, i, start) == 1)
			for (lidia_size_t l = start + 1; l < A.columns; ++l)
				if (modul.member(A, i, l) != A.Zero) {
					q = modul.member(A, i, l);
					modul.subtract_multiple_of_column(A, l, q, start, A.rows - 1);

					for (lidia_size_t m = 0; m < A.columns; m++)
						TR.sto(m, l, TR.member(m, l) - q * TR.member(m, start));
				}
		++start;
	}

	lidia_size_t *RET = new lidia_size_t[A.rows];
	k = 1;

	start--;

	// Phase 2
	for (--i; i >= startr; --i) {
		if (modul.member(A, i, start) == A.Zero)
			RET[k++] = i;
		else {
			if (modul.member(A, i, start) != 1)
				for (lidia_size_t l = start + 1; l < A.columns; ++l) {
					if (modul.member(A, i, l) != A.Zero) {
						pos_div_rem(q, r, modul.member(A, i, l), modul.member(A, i, start));
						modul.subtract_multiple_of_column(A, l, q, start, A.rows - 1);

						for (lidia_size_t m = 0; m < A.columns; m++)
							TR.sto(m, l, TR.member(m, l) - q * TR.member(m, start));
					}
				}
		}
		--start;
	}

	if (k == 1) {
		delete [] RET;
		RET = NULL;
	}
	else
		RET[0] = k-1;

	return RET;
}



template< class T, class REP, class REP1 >
lidia_size_t *
normalization_kernel< T, REP, REP1 >::normalizeHybrid_Std (MR< T > &A,
							   trans_matrix &TR,
							   matrix< bigint > & tran,
							   lidia_size_t s)
{
	lidia_size_t start = s, i, k, p, l;
	T q, r;
	matrix< bigint > otran;

	tran.reset();
	otran = tran;
	tran.resize(A.columns, A.columns);
	tran.diag(1, 0);

	lidia_size_t *RET = new lidia_size_t[A.rows];

	for (i = start + 1, p = 0; i < A.columns; i++, p++) {
		for (k = i - 1, l = p; k >= start; k--, l--) {
			if (modul.member(A, l, i) != A.Zero && modul.member(A, l, k) != A.Zero) {
				pos_div_rem(q, r, modul.member(A, l, i), modul.member(A, l, k));
				modul.subtract_multiple_of_column(A, i, q, k, A.rows - 1);

				for (lidia_size_t m = 0; m < A.columns; m++)
					tran.sto(m, i, tran.member(m, i) - q * tran.member(m, k));
			}
		}

		if (i == (A.rows/2)) {
			TR.store_matrix(tran);

			tran = otran;
			tran.resize(A.columns, A.columns);
			tran.diag(1, 0);
		}
	}

	k = 1;
	for (i = 0; i < A.rows; i++)
		if (modul.member(A, i, start + i) == A.Zero) {
			RET[k] = i;
			k++;
		}
	RET[0] = k - 1;

	if (RET[0] == 0) {
		delete[] RET;
		RET = NULL;
	}

	return RET;
}



//
// normalizeHybrid_ChouColllins
//

template< class T, class REP, class REP1 >
lidia_size_t *
normalization_kernel< T, REP, REP1 >::normalizeHybrid_ChouCollins (MR< T > &A,
								   lidia_size_t startr,
								   lidia_size_t startc) const
{
	lidia_size_t start = startc, i, k;
	T q, r;

#ifdef LIDIA_COMPUTE_RUNTIME_INFO
	std::ofstream profile("normalize_max.profile");
	T MAX = 0, GMAX = 0;
#endif

	// Phase 1
	for (i = startr; i < A.rows; i++) {

#ifdef LIDIA_COMPUTE_RUNTIME_INFO
		modul.max(A, MAX);
		if (GMAX < MAX)
			GMAX = MAX;
#endif

		if (modul.member(A, i, start) == 1)
			for (lidia_size_t l = start + 1; l < A.columns; l++)
				if (modul.member(A, i, l) != A.Zero)
					modul.subtract_multiple_of_column(A, l, modul.member(A, i, l), start, A.rows - 1);
		++start;
	}

	lidia_size_t *RET = new lidia_size_t[A.rows];
	k = 1;

	start = startc;

	// Phase 2
	lidia_size_t p, l;

	for (i = start + 1, p = startr; i < A.columns; i++, p++) {
		std::cout << "i = " << i << " p = " << p << std::endl;
		for (k = i - 1, l = p; k >= start; k--, l--) {
			std::cout << "\t k = " << k << " l = " << l << std::endl;
#ifdef LIDIA_COMPUTE_RUNTIME_INFO
			modul.max(A, MAX);
			if (GMAX < MAX)
				GMAX = MAX;
#endif
			if (modul.member(A, l, k) != 1)
				if (modul.member(A, l, i) != A.Zero && modul.member(A, l, k) != A.Zero) {
					pos_div_rem(q, r, modul.member(A, l, i), modul.member(A, l, k));
					modul.subtract_multiple_of_column(A, i, q, k, A.rows - 1);
				}
		}
	}

	k = 1;
	for (i = startr; i < A.rows; i++) {
		if (modul.member(A, i, start + i - startr) == A.Zero) {
			RET[k] = i;
			k++;
		}
	}

	RET[0] = k - 1;

	if (RET[0] == 0) {
		delete[] RET;
		RET = NULL;
	}

#ifdef LIDIA_COMPUTE_RUNTIME_INFO
	profile << GMAX << std::endl;
	profile.close();
#endif

	return RET;
}



template< class T, class REP, class REP1 >
lidia_size_t *
normalization_kernel< T, REP, REP1 >::normalizeHybrid_ChouCollins (MR< T > &A,
								   matrix< bigint > &TR,
								   lidia_size_t startr,
								   lidia_size_t startc) const
{
	lidia_size_t start = startc, i, k;
	T q, r;

	// Phase 1
	for (i = startr; i < A.rows; ++i) {
		if (modul.member(A, i, start) == 1)
			for (lidia_size_t l = start + 1; l < A.columns; ++l)
				if (modul.member(A, i, l) != A.Zero) {
					q = modul.member(A, i, l);
					modul.subtract_multiple_of_column(A, l, q, start, A.rows - 1);

					for (lidia_size_t m = 0; m < A.columns; m++)
						TR.sto(m, l, TR.member(m, l) - q * TR.member(m, start));
				}
		++start;
	}

	lidia_size_t *RET = new lidia_size_t[A.rows];
	k = 1;

	start = startc;

	// Phase 2
	lidia_size_t p, l;

	for (i = start + 1, p = startr; i < A.columns; i++, p++) {
		for (k = i - 1, l = p; k >= start; k--, l--) {
			if (modul.member(A, l, k) != 1)
				if (modul.member(A, l, i) != A.Zero && modul.member(A, l, k) != A.Zero) {
					pos_div_rem(q, r, modul.member(A, l, i), modul.member(A, l, k));
					modul.subtract_multiple_of_column(A, i, q, k, A.rows - 1);

					for (lidia_size_t m = 0; m < A.columns; m++)
						TR.sto(m, i, TR.member(m, i) - q * TR.member(m, k));
				}
		}
	}

	k = 1;
	for (i = startr; i < A.rows; i++)
		if (modul.member(A, i, start + i - startr) == A.Zero) {
			RET[k] = i;
			k++;
		}
	RET[0] = k - 1;

	if (RET[0] == 0) {
		delete[] RET;
		RET = NULL;
	}

	return RET;
}



template< class T, class REP, class REP1 >
lidia_size_t *
normalization_kernel< T, REP, REP1 >::normalizeHybrid_ChouCollins (MR< T > &A,
								   trans_matrix &TR,
								   matrix< bigint > & tran,
								   lidia_size_t s)
{
	lidia_size_t start = s, i, k, p, l;
	T q, r;
	matrix< bigint > otran;

	tran.reset();
	otran = tran;
	tran.resize(A.columns, A.columns);
	tran.diag(1, 0);

	lidia_size_t *RET = new lidia_size_t[A.rows];

	for (i = start + 1, p = 0; i < A.columns; i++, p++) {
		for (k = i - 1, l = p; k >= start; k--, l--) {
			if (modul.member(A, l, i) != A.Zero && modul.member(A, l, k) != A.Zero) {
				pos_div_rem(q, r, modul.member(A, l, i), modul.member(A, l, k));
				modul.subtract_multiple_of_column(A, i, q, k, A.rows - 1);

				for (lidia_size_t m = 0; m < A.columns; m++)
					tran.sto(m, i, tran.member(m, i) - q * tran.member(m, k));
			}
		}

		if (i == (A.rows/2)) {
			TR.store_matrix(tran);

			tran = otran;
			tran.resize(A.columns, A.columns);
			tran.diag(1, 0);
		}
	}

	k = 1;
	for (i = 0; i < A.rows; i++)
		if (modul.member(A, i, start + i) == A.Zero) {
			RET[k] = i;
			k++;
		}
	RET[0] = k - 1;

	if (RET[0] == 0) {
		delete[] RET;
		RET = NULL;
	}

	return RET;
}



//
// normalize modular Standard
//

template< class T, class REP, class REP1 >
lidia_size_t *
normalization_kernel< T, REP, REP1 >::normalizeMod_Std (MR< T > &A,
							lidia_size_t startr,
							lidia_size_t startc) const
{
	lidia_size_t i, j, k = 1, l;
	T q, r;

#ifdef LIDIA_COMPUTE_RUNTIME_INFO
	std::ofstream profile("normalize_max.profile");
	T MAX = 0, GMAX = 0;
#endif

	T mod = 2;
	for (i = startr, j = startc; i < A.rows && j < A.columns; i++, j++)
		LiDIA::multiply(mod, mod, abs(modul.member(A, i, j)));

	lidia_size_t *RET = new lidia_size_t[A.rows];

	for (i = startr, j = startc; i < A.rows && j < A.columns; i++, j++)
		for (l = j + 1; l < A.columns; l++) {
			LiDIA::best_remainder(r, modul.member(A, i, l), mod);
			modul.sto(A, i, l, r);
		}

	for (i = A.columns - 1, j = A.rows - 1; i >= startc && j >= startr; i--, j--) {

#ifdef LIDIA_COMPUTE_RUNTIME_INFO
		modul.max(A, MAX);
		if (GMAX < MAX)
			GMAX = MAX;
#endif

		if (modul.member(A, j, i) != A.Zero) {
			for (lidia_size_t l = i + 1; l < A.columns; l++)
				if (modul.member(A, j, l) != A.Zero) {
					pos_div_rem(q, r, modul.member(A, j, l), modul.member(A, j, i));
					modul.normalize_column_mod(A, l, q, i, A.rows - 1, mod);
				}
		}
		else {
			RET[k] = j;
			k++;
		}
	}

	if (k == 1) {
		delete [] RET;
		RET = NULL;
	}
	else
		RET[0] = k - 1;

#ifdef LIDIA_COMPUTE_RUNTIME_INFO
	profile << GMAX << std::endl;
	profile.close();
#endif

	return RET;
}



//
// normalize modular ChouCollins
//

template< class T, class REP, class REP1 >
lidia_size_t *
normalization_kernel< T, REP, REP1 >::normalizeMod_ChouCollins (MR< T > &A,
								lidia_size_t startr,
								lidia_size_t startc) const
{
	lidia_size_t start = startc, i, j, k, p, l;
	T q, r;

#ifdef LIDIA_COMPUTE_RUNTIME_INFO
	std::ofstream profile("normalize_max.profile");
	T MAX = 0, GMAX = 0;
#endif

	T mod = 2;
	for (i = startr, p = startc; p < A.columns && i < A.rows; i++, p++)
		LiDIA::multiply(mod, mod, abs(modul.member(A, i, p)));

	lidia_size_t *RET = new lidia_size_t[A.rows];

	for (i = startr, j = startc; i < A.rows && j < A.columns; i++, j++)
		for (l = j + 1; l < A.columns; l++) {
			LiDIA::best_remainder(r, modul.member(A, i, l), mod);
			modul.sto(A, i, l, r);
		}

	for (i = start + 1, p = startr; i < A.columns; i++, p++) {
		for (k = i - 1, l = p; k >= start; k--, l--) {

#ifdef LIDIA_COMPUTE_RUNTIME_INFO
			modul.max(A, MAX);
			if (GMAX < MAX)
				GMAX = MAX;
#endif

			if (modul.member(A, l, i) != A.Zero && modul.member(A, l, k) != A.Zero) {
				pos_div_rem(q, r, modul.member(A, l, i), modul.member(A, l, k));
				modul.normalize_column_mod(A, i, q, k, A.rows - 1, mod);
			}
		}
	}

	k = 1;
	for (i = startr; i < A.rows; i++) {
		if (modul.member(A, i, start + i - startr) == A.Zero) {
			RET[k] = i;
			k++;
		}
	}

	RET[0] = k - 1;

	if (RET[0] == 0) {
		delete[] RET;
		RET = NULL;
	}

#ifdef LIDIA_COMPUTE_RUNTIME_INFO
	profile << GMAX << std::endl;
	profile.close();
#endif

	return RET;
}



#undef DMESSAGE



#ifdef LIDIA_NAMESPACE
# ifndef IN_NAMESPACE_LIDIA
}	// end of namespace LiDIA
# endif
#endif



#endif	// LIDIA_NORMALIZE_KERNEL_CC_GUARD_
