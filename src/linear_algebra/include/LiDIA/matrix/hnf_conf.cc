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


#ifndef LIDIA_HNF_CONF_CC_GUARD_
#define LIDIA_HNF_CONF_CC_GUARD_


#ifndef LIDIA_HNF_CONF_H_GUARD_
# include	"LiDIA/matrix/hnf_conf.h"
#endif



#ifdef LIDIA_NAMESPACE
# ifndef IN_NAMESPACE_LIDIA
namespace LiDIA {
# endif
#endif



//
// Hermite normal form configuration: hermite
//

template< class T, class REP, class REP1 >
inline hermite< T, REP, REP1 >::hermite ()
{
}



template< class T, class REP, class REP1 >
inline hermite< T, REP, REP1 >::~hermite ()
{
}



//
// init
//

template< class T, class REP, class REP1 >
inline void
hermite< T, REP, REP1 >::init (MR< T > &, const T &)
{
}



//
// mgcd
//

template< class T, class REP, class REP1 >
inline lidia_size_t
hermite< T, REP, REP1 >::mgcd (MR< T > &A,
			       lidia_size_t startr,
			       lidia_size_t startc) const
{
	if (elim_modul.mgcd_linear(A, startr, startc, startr))
		return 1;
	return 0;
}



template< class T, class REP, class REP1 >
inline lidia_size_t
hermite< T, REP, REP1 >::mgcd (MR< T > &A,
			       matrix< bigint > &TR,
			       lidia_size_t startr,
			       lidia_size_t startc) const
{
	if (elim_modul.mgcd_linear(A, TR, startr, startc, startr))
		return 1;
	return 0;
}



//
// normalization
//

template< class T, class REP, class REP1 >
inline bool
hermite< T, REP, REP1 >::normalize_row (MR< T > &A,
					lidia_size_t startr,
					lidia_size_t startc,
					lidia_size_t index,
					lidia_size_t len) const
{
	return elim_modul.normalize_row(A, startr, startc, index, len);
}



template< class T, class REP, class REP1 >
inline bool
hermite< T, REP, REP1 >::normalize_row (MR< T > &A,
					matrix< bigint > &TR,
					lidia_size_t startr,
					lidia_size_t startc,
					lidia_size_t index,
					lidia_size_t len) const
{
	return elim_modul.normalize_row(A, TR, startr, startc, index, len);
}



//
// Hermite normal form configuration: Bradley
//

template< class T, class REP, class REP1 >
inline
Bradley< T, REP, REP1 >::Bradley ()
{
}



template< class T, class REP, class REP1 >
inline
Bradley< T, REP, REP1 >::~Bradley ()
{
}



//
// init
//

template< class T, class REP, class REP1 >
inline void
Bradley< T, REP, REP1 >::init (MR< T > &, const T &)
{
}



//
// mgcd
//

template< class T, class REP, class REP1 >
inline lidia_size_t
Bradley< T, REP, REP1 >::mgcd (MR< T > &A,
			       lidia_size_t startr,
			       lidia_size_t startc) const
{
	if (elim_modul.mgcd_bradley(A, startr, startc, startr))
		return 1;
	return 0;
}



template< class T, class REP, class REP1 >
inline lidia_size_t
Bradley< T, REP, REP1 >::mgcd (MR< T > &A,
			       matrix< bigint > &TR,
			       lidia_size_t startr,
			       lidia_size_t startc) const
{
	if (elim_modul.mgcd_bradley(A, TR, startr, startc, startr))
		return 1;
	return 0;
}



//
// normalization
//

template< class T, class REP, class REP1 >
inline bool
Bradley< T, REP, REP1 >::normalize_row (MR< T > &A,
					lidia_size_t startr,
					lidia_size_t startc,
					lidia_size_t index,
					lidia_size_t len) const
{
	return elim_modul.normalize_row(A, startr, startc, index, len);
}



template< class T, class REP, class REP1 >
inline bool
Bradley< T, REP, REP1 >::normalize_row (MR< T > &A,
					matrix< bigint > &TR,
					lidia_size_t startr,
					lidia_size_t startc,
					lidia_size_t index,
					lidia_size_t len) const
{
	return elim_modul.normalize_row(A, TR, startr, startc, index, len);
}



//
// Hermite normal form configuration: Iliopoulos
//

template< class T, class REP, class REP1 >
inline
Iliopoulos< T, REP, REP1 >::Iliopoulos ()
{
}



template< class T, class REP, class REP1 >
inline
Iliopoulos< T, REP, REP1 >::~Iliopoulos ()
{
}



//
// init
//

template< class T, class REP, class REP1 >
inline void
Iliopoulos< T, REP, REP1 >::init (MR< T > &, const T &)
{
}



//
// mgcd
//

template< class T, class REP, class REP1 >
inline lidia_size_t
Iliopoulos< T, REP, REP1 >::mgcd (MR< T > &A,
				  lidia_size_t startr,
				  lidia_size_t startc) const
{
	if (elim_modul.mgcd_ilio(A, startr, startc, startr))
		return 1;
	return 0;
}



template< class T, class REP, class REP1 >
inline lidia_size_t
Iliopoulos< T, REP, REP1 >::mgcd (MR< T > &A,
				  matrix< bigint > &TR,
				  lidia_size_t startr,
				  lidia_size_t startc) const
{
	if (elim_modul.mgcd_ilio(A, TR, startr, startc, startr))
		return 1;
	return 0;
}



//
// normalization
//

template< class T, class REP, class REP1 >
inline bool
Iliopoulos< T, REP, REP1 >::normalize_row (MR< T > &A,
					   lidia_size_t startr,
					   lidia_size_t startc,
					   lidia_size_t index,
					   lidia_size_t len) const
{
	return elim_modul.normalize_row(A, startr, startc, index, len);
}



template< class T, class REP, class REP1 >
inline bool
Iliopoulos< T, REP, REP1 >::normalize_row (MR< T > &A,
					   matrix< bigint > &TR,
					   lidia_size_t startr,
					   lidia_size_t startc,
					   lidia_size_t index,
					   lidia_size_t len) const
{
	return elim_modul.normalize_row(A, TR, startr, startc, index, len);
}



//
// Hermite normal form configuration: opt
//

template< class T, class REP, class REP1 >
inline
opt< T, REP, REP1 >::opt ()
{
}



template< class T, class REP, class REP1 >
inline
opt< T, REP, REP1 >::~opt ()
{
}



//
// init
//

template< class T, class REP, class REP1 >
inline void
opt< T, REP, REP1 >::init (MR< T > &, const T &)
{
}



//
// mgcd
//

template< class T, class REP, class REP1 >
inline lidia_size_t
opt< T, REP, REP1 >::mgcd (MR< T > &A,
			   lidia_size_t startr,
			   lidia_size_t startc) const
{
	if (elim_modul.mgcd_opt(A, startr, startc, startr))
		return 1;
	return 0;
}



template< class T, class REP, class REP1 >
inline lidia_size_t
opt< T, REP, REP1 >::mgcd (MR< T > &A,
			   matrix< bigint > &TR,
			   lidia_size_t startr,
			   lidia_size_t startc) const
{
	if (elim_modul.mgcd_opt(A, TR, startr, startc, startr))
		return 1;
	return 0;
}



//
// normalization
//

template< class T, class REP, class REP1 >
inline bool
opt< T, REP, REP1 >::normalize_row (MR< T > &A,
				    lidia_size_t startr,
				    lidia_size_t startc,
				    lidia_size_t index,
				    lidia_size_t len) const
{
	return elim_modul.normalize_row(A, startr, startc, index, len);
}



template< class T, class REP, class REP1 >
inline bool
opt< T, REP, REP1 >::normalize_row (MR< T > &A,
				    matrix< bigint > &TR,
				    lidia_size_t startr,
				    lidia_size_t startc,
				    lidia_size_t index,
				    lidia_size_t len) const
{
	return elim_modul.normalize_row(A, TR, startr, startc, index, len);
}



//
// Hermite normal form configuration: Blankinship
//

template< class T, class REP, class REP1 >
inline
Blankinship< T, REP, REP1 >::Blankinship ()
{
}



template< class T, class REP, class REP1 >
inline
Blankinship< T, REP, REP1 >::~Blankinship ()
{
}



//
// init
//

template< class T, class REP, class REP1 >
inline void
Blankinship< T, REP, REP1 >::init (MR< T > &, const T &)
{
}



//
// mgcd
//

template< class T, class REP, class REP1 >
inline lidia_size_t
Blankinship< T, REP, REP1 >::mgcd (MR< T > &A,
				   lidia_size_t startr,
				   lidia_size_t startc) const
{
	if (elim_modul.mgcd_linear(A, startr, startc, startr))
		return 1;
	return 0;
}



template< class T, class REP, class REP1 >
inline lidia_size_t
Blankinship< T, REP, REP1 >::mgcd (MR< T > &A,
				   matrix< bigint > &TR,
				   lidia_size_t startr,
				   lidia_size_t startc) const
{
	if (elim_modul.mgcd_linear(A, TR, startr, startc, startr))
		return 1;
	return 0;
}



//
// normalization
//

template< class T, class REP, class REP1 >
inline bool
Blankinship< T, REP, REP1 >::normalize_row (MR< T > &A,
					    lidia_size_t startr,
					    lidia_size_t startc,
					    lidia_size_t index,
					    lidia_size_t len) const
{
	return elim_modul.normalize_row(A, startr, startc, index, len);
}



template< class T, class REP, class REP1 >
inline bool
Blankinship< T, REP, REP1 >::normalize_row (MR< T > &A,
					    matrix< bigint > &TR,
					    lidia_size_t startr,
					    lidia_size_t startc,
					    lidia_size_t index,
					    lidia_size_t len) const
{
	return elim_modul.normalize_row(A, TR, startr, startc, index, len);
}



//
// Hermite normal form configuration: Best_remainder
//

template< class T, class REP, class REP1 >
inline
Best_remainder< T, REP, REP1 >::Best_remainder ()
{
}



template< class T, class REP, class REP1 >
inline
Best_remainder< T, REP, REP1 >::~Best_remainder ()
{
}



//
// init
//

template< class T, class REP, class REP1 >
inline void
Best_remainder< T, REP, REP1 >::init (MR< T > &, const T &)
{
}



//
// mgcd
//

template< class T, class REP, class REP1 >
inline lidia_size_t
Best_remainder< T, REP, REP1 >::mgcd (MR< T > &A,
				      lidia_size_t startr,
				      lidia_size_t startc) const
{
	lidia_size_t index, count;
	while ((count = pivot_modul.normal_pivot(A, startr, startc, index)) > 1)
		if (!elim_modul.normalize_row(A, startr, startc, index, startr))
			return -1;
	if (count == 0)
		return 0;

	if (index != startc)
		modul.swap_columns(A, startc, index);

	return 1;
}



template< class T, class REP, class REP1 >
inline lidia_size_t
Best_remainder< T, REP, REP1 >::mgcd (MR< T > &A,
				      matrix< bigint > &TR,
				      lidia_size_t startr,
				      lidia_size_t startc) const
{
	lidia_size_t index, count;
	while ((count = pivot_modul.min_abs_of_row(A, startr, startc, index)) > 1)
		if (!elim_modul.normalize_row(A, TR, startr, startc, index, startr))
			return -1;

	if (count == 0)
		return 0;

	if (index != startc) {
		modul.swap_columns(A, startc, index);
		TR.swap_columns(startc, index);
	}

	return 1;
}



//
// normalization
//

template< class T, class REP, class REP1 >
inline bool
Best_remainder< T, REP, REP1 >::normalize_row (MR< T > &A,
					       lidia_size_t startr,
					       lidia_size_t startc,
					       lidia_size_t index,
					       lidia_size_t len) const
{
	return elim_modul.normalize_row(A, startr, startc, index, len);
}



template< class T, class REP, class REP1 >
inline bool
Best_remainder< T, REP, REP1 >::normalize_row (MR< T > &A,
					       matrix< bigint > &TR,
					       lidia_size_t startr,
					       lidia_size_t startc,
					       lidia_size_t index,
					       lidia_size_t len) const
{
	return elim_modul.normalize_row(A, TR, startr, startc, index, len);
}



//
// Hermite normal form configuration: havas_best_remainder
//

template< class T, class REP, class REP1 >
inline
havas_best_remainder< T, REP, REP1 >::havas_best_remainder ()
{
}



template< class T, class REP, class REP1 >
inline
havas_best_remainder< T, REP, REP1 >::~havas_best_remainder ()
{
}



//
// init
//

template< class T, class REP, class REP1 >
inline void
havas_best_remainder< T, REP, REP1 >::init (MR< T > &, const T &)
{
}



//
// mgcd
//

template< class T, class REP, class REP1 >
inline lidia_size_t
havas_best_remainder< T, REP, REP1 >::mgcd (MR< T > &A,
					    lidia_size_t startr,
					    lidia_size_t startc) const
{
	lidia_size_t index, count;
	while ((count = pivot_modul.min_abs_of_row(A, startr, startc, index)) > 1)
		if (!elim_modul.normalize_row(A, startr, startc, index, startr))
			return -1;

	if (count == 0)
		return 0;

	if (index != startc)
		modul.swap_columns(A, startc, index);

	return 1;
}



template< class T, class REP, class REP1 >
inline lidia_size_t
havas_best_remainder< T, REP, REP1 >::mgcd (MR< T > &A,
					    matrix< bigint > &TR,
					    lidia_size_t startr,
					    lidia_size_t startc) const
{
	lidia_size_t index, count;
	while ((count = pivot_modul.min_abs_of_row(A, startr, startc, index)) > 1)
		if (!elim_modul.normalize_row(A, TR, startr, startc, index, startr))
			return -1;

	if (count == 0)
		return 0;

	if (index != startc) {
		modul.swap_columns(A, startc, index);
		TR.swap_columns(startc, index);
	}

	return 1;
}



//
// normalization
//

template< class T, class REP, class REP1 >
inline bool
havas_best_remainder< T, REP, REP1 >::normalize_row (MR< T > &A,
						     lidia_size_t startr,
						     lidia_size_t startc,
						     lidia_size_t index,
						     lidia_size_t len) const
{
	return elim_modul.normalize_row(A, startr, startc, index, len);
}



template< class T, class REP, class REP1 >
inline bool
havas_best_remainder< T, REP, REP1 >::normalize_row (MR< T > &A,
						     matrix< bigint > &TR,
						     lidia_size_t startr,
						     lidia_size_t startc,
						     lidia_size_t index,
						     lidia_size_t len) const
{
	return elim_modul.normalize_row(A, TR, startr, startc, index, len);
}



//
// Hermite normal form configuration: havas_euclidean_norm
//

template< class T, class REP, class REP1 >
inline
havas_euclidean_norm< T, REP, REP1 >::havas_euclidean_norm ()
{
}



template< class T, class REP, class REP1 >
inline
havas_euclidean_norm< T, REP, REP1 >::~havas_euclidean_norm ()
{
}



//
// init
//

template< class T, class REP, class REP1 >
inline void
havas_euclidean_norm< T, REP, REP1 >::init (MR< T > &, const T &)
{
}



//
// mgcd
//

template< class T, class REP, class REP1 >
inline lidia_size_t
havas_euclidean_norm< T, REP, REP1 >::mgcd (MR< T > &A,
					    lidia_size_t startr,
					    lidia_size_t startc) const
{
	lidia_size_t index, count;
	while ((count = pivot_modul.minimal_norm(A, startr, startc, index, 2)) > 1)
		if (!elim_modul.normalize_row(A, startr, startc, index, startr))
			return -1;

	if (count == 0)
		return 0;

	if (index != startc)
		modul.swap_columns(A, startc, index);
	return 1;

}



template< class T, class REP, class REP1 >
inline lidia_size_t
havas_euclidean_norm< T, REP, REP1 >::mgcd (MR< T > &A,
					    matrix< bigint > &TR,
					    lidia_size_t startr,
					    lidia_size_t startc) const
{
	lidia_size_t index, count;
	while ((count = pivot_modul.minimal_norm(A, startr, startc, index, 2)) > 1)
		if (!elim_modul.normalize_row(A, TR, startr, startc, index, startr))
			return -1;

	if (count == 0)
		return 0;

	if (index != startc) {
		modul.swap_columns(A, startc, index);
		TR.swap_columns(startc, index);
	}
	return 1;
}



//
// normalization
//

template< class T, class REP, class REP1 >
inline bool
havas_euclidean_norm< T, REP, REP1 >::normalize_row (MR< T > &A,
						     lidia_size_t startr,
						     lidia_size_t startc,
						     lidia_size_t index,
						     lidia_size_t len) const
{
	return elim_modul.normalize_row(A, startr, startc, index, len);
}



template< class T, class REP, class REP1 >
inline bool
havas_euclidean_norm< T, REP, REP1 >::normalize_row (MR< T > &A,
						     matrix< bigint > &TR,
						     lidia_size_t startr,
						     lidia_size_t startc,
						     lidia_size_t index,
						     lidia_size_t len) const
{
	return elim_modul.normalize_row(A, TR, startr, startc, index, len);
}



//
// Hermite normal form configuration: havas_minimal_norm
//

template< class T, class REP, class REP1 >
inline
havas_minimal_norm< T, REP, REP1 >::havas_minimal_norm ()
{
}



template< class T, class REP, class REP1 >
inline
havas_minimal_norm< T, REP, REP1 >::~havas_minimal_norm ()
{
}



//
// init
//

template< class T, class REP, class REP1 >
inline void
havas_minimal_norm< T, REP, REP1 >::init (MR< T > &, const T &)
{
}



//
// mgcd
//

template< class T, class REP, class REP1 >
inline lidia_size_t
havas_minimal_norm< T, REP, REP1 >::mgcd (MR< T > &A,
					  lidia_size_t startr,
					  lidia_size_t startc) const
{
	lidia_size_t index, count;
	while ((count = pivot_modul.minimal_norm(A, startr, startc, index, 1)) > 1)
		if (!elim_modul.normalize_row(A, startr, startc, index, startr))
			return -1;

	if (count == 0)
		return 0;

	if (index != startc)
		modul.swap_columns(A, startc, index);
	return 1;

}



template< class T, class REP, class REP1 >
inline lidia_size_t
havas_minimal_norm< T, REP, REP1 >::mgcd (MR< T > &A,
					  matrix< bigint > &TR,
					  lidia_size_t startr,
					  lidia_size_t startc) const
{
	lidia_size_t index, count;
	while ((count = pivot_modul.minimal_norm(A, startr, startc, index, 1)) > 1)
		if (!elim_modul.normalize_row(A, TR, startr, startc, index, startr))
			return -1;

	if (count == 0)
		return 0;

	if (index != startc) {
		modul.swap_columns(A, startc, index);
		TR.swap_columns(startc, index);
	}
	return 1;
}



//
// normalization
//

template< class T, class REP, class REP1 >
inline bool
havas_minimal_norm< T, REP, REP1 >::normalize_row (MR< T > &A,
						   lidia_size_t startr,
						   lidia_size_t startc,
						   lidia_size_t index,
						   lidia_size_t len) const
{
	return elim_modul.normalize_row(A, startr, startc, index, len);
}



template< class T, class REP, class REP1 >
inline bool
havas_minimal_norm< T, REP, REP1 >::normalize_row (MR< T > &A,
						   matrix< bigint > &TR,
						   lidia_size_t startr,
						   lidia_size_t startc,
						   lidia_size_t index,
						   lidia_size_t len) const
{
	return elim_modul.normalize_row(A, TR, startr, startc, index, len);
}



//
// Hermite normal form configuration: havas_sorting_gcd
//

template< class T, class REP, class REP1 >
inline
havas_sorting_gcd< T, REP, REP1 >::havas_sorting_gcd ()
{
}



template< class T, class REP, class REP1 >
inline
havas_sorting_gcd< T, REP, REP1 >::~havas_sorting_gcd ()
{
}



//
// init
//

template< class T, class REP, class REP1 >
inline void
havas_sorting_gcd< T, REP, REP1 >::init (MR< T > &, const T &)
{
}



//
// mgcd
//

template< class T, class REP, class REP1 >
inline lidia_size_t
havas_sorting_gcd< T, REP, REP1 >::mgcd (MR< T > &A,
					 lidia_size_t startr,
					 lidia_size_t startc) const
{
	lidia_size_t index, count;
	while ((count = pivot_modul.pivot_sorting_gcd(A, startr, startc, index)) > 1)
		if (!elim_modul.normalize_row(A, startr, startc, index, startr))
			return -1;

	if (count == 0)
		return 0;

	if (index != startc)
		modul.swap_columns(A, startc, index);
	return 1;

}



template< class T, class REP, class REP1 >
inline lidia_size_t
havas_sorting_gcd< T, REP, REP1 >::mgcd (MR< T > &A,
					 matrix< bigint > &TR,
					 lidia_size_t startr,
					 lidia_size_t startc) const
{
	lidia_size_t index, count;
	while ((count = pivot_modul.pivot_sorting_gcd(A, startr, startc, index)) > 1)
		if (!elim_modul.normalize_row(A, TR, startr, startc, index, startr))
			return -1;

	if (count == 0)
		return 0;

	if (index != startc) {
		modul.swap_columns(A, startc, index);
		TR.swap_columns(startc, index);
	}
	return 1;
}



//
// normalization
//

template< class T, class REP, class REP1 >
inline bool
havas_sorting_gcd< T, REP, REP1 >::normalize_row (MR< T > &A,
						  lidia_size_t startr,
						  lidia_size_t startc,
						  lidia_size_t index,
						  lidia_size_t len) const
{
	return elim_modul.normalize_row(A, startr, startc, index, len);
}



template< class T, class REP, class REP1 >
inline bool
havas_sorting_gcd< T, REP, REP1 >::normalize_row (MR< T > &A,
						  matrix< bigint > &TR,
						  lidia_size_t startr,
						  lidia_size_t startc,
						  lidia_size_t index,
						  lidia_size_t len) const
{
	return elim_modul.normalize_row(A, TR, startr, startc, index, len);
}



//
// Hermite normal form configuration havas_min_no_of_elements
//

template< class T, class REP, class REP1 >
inline
havas_min_no_of_elements< T, REP, REP1 >::havas_min_no_of_elements ()
{
}



template< class T, class REP, class REP1 >
inline
havas_min_no_of_elements< T, REP, REP1 >::~havas_min_no_of_elements ()
{
}



//
// init
//

template< class T, class REP, class REP1 >
inline void
havas_min_no_of_elements< T, REP, REP1 >::init (MR< T > &, const T &)
{
}



//
// mgcd
//

template< class T, class REP, class REP1 >
inline lidia_size_t
havas_min_no_of_elements< T, REP, REP1 >::mgcd (MR< T > &A,
						lidia_size_t startr,
						lidia_size_t startc) const
{
	lidia_size_t index, count;
	while ((count = pivot_modul.minimal_norm(A, startr, startc, index, 0)) > 1)
		if (!elim_modul.normalize_row(A, startr, startc, index, startr))
			return -1;

	if (count == 0)
		return 0;

	if (index != startc)
		modul.swap_columns(A, startc, index);
	return 1;

}



template< class T, class REP, class REP1 >
inline lidia_size_t
havas_min_no_of_elements< T, REP, REP1 >::mgcd (MR< T > &A,
						matrix< bigint > &TR,
						lidia_size_t startr,
						lidia_size_t startc) const
{
	lidia_size_t index, count;
	while ((count = pivot_modul.minimal_norm(A, startr, startc, index, 0)) > 1)
		if (!elim_modul.normalize_row(A, TR, startr, startc, index, startr))
			return -1;

	if (count == 0)
		return 0;

	if (index != startc) {
		modul.swap_columns(A, startc, index);
		TR.swap_columns(startc, index);
	}
	return 1;
}



//
// normalization
//

template< class T, class REP, class REP1 >
inline bool
havas_min_no_of_elements< T, REP, REP1 >::normalize_row (MR< T > &A,
							 lidia_size_t startr,
							 lidia_size_t startc,
							 lidia_size_t index,
							 lidia_size_t len) const
{
	return elim_modul.normalize_row(A, startr, startc, index, len);
}



template< class T, class REP, class REP1 >
inline bool
havas_min_no_of_elements< T, REP, REP1 >::normalize_row (MR< T > &A,
							 matrix< bigint > &TR,
							 lidia_size_t startr,
							 lidia_size_t startc,
							 lidia_size_t index,
							 lidia_size_t len) const
{
	return elim_modul.normalize_row(A, TR, startr, startc, index, len);
}



//
// Hermite normal form configuration: min_abs_of_row + minimal_euclidean_norm
//

template< class T, class REP, class REP1 >
inline
havas_min_abs_of_row_plus_minimal_euclidean_norm< T, REP, REP1 >::havas_min_abs_of_row_plus_minimal_euclidean_norm ()
{
}



template< class T, class REP, class REP1 >
inline
havas_min_abs_of_row_plus_minimal_euclidean_norm< T, REP, REP1 >::~havas_min_abs_of_row_plus_minimal_euclidean_norm ()
{
}



//
// init
//

template< class T, class REP, class REP1 >
inline void
havas_min_abs_of_row_plus_minimal_euclidean_norm< T, REP, REP1 >::init (MR< T > &, const T &)
{
}



//
// mgcd
//

template< class T, class REP, class REP1 >
inline lidia_size_t
havas_min_abs_of_row_plus_minimal_euclidean_norm< T, REP, REP1 >::mgcd (MR< T > &A,
									lidia_size_t startr,
									lidia_size_t startc) const
{
	lidia_size_t index, count;
	while ((count = pivot_modul.min_abs_of_row_plus_minimal_norm(A, startr, startc, index, 2)) > 1)
		if (!elim_modul.normalize_row(A, startr, startc, index, startr))
			return -1;

	if (count == 0)
		return 0;

	if (index != startc)
		modul.swap_columns(A, startc, index);
	return 1;

}



template< class T, class REP, class REP1 >
inline lidia_size_t
havas_min_abs_of_row_plus_minimal_euclidean_norm< T, REP, REP1 >::mgcd (MR< T > &A,
									matrix< bigint > &TR,
									lidia_size_t startr,
									lidia_size_t startc) const
{
	lidia_size_t index, count;
	while ((count = pivot_modul.min_abs_of_row_plus_minimal_norm(A, startr, startc, index, 2)) > 1)
		if (!elim_modul.normalize_row(A, TR, startr, startc, index, startr))
			return -1;

	if (count == 0)
		return 0;

	if (index != startc) {
		modul.swap_columns(A, startc, index);
		TR.swap_columns(startc, index);
	}
	return 1;
}



//
// normalization
//

template< class T, class REP, class REP1 >
inline bool
havas_min_abs_of_row_plus_minimal_euclidean_norm< T, REP, REP1 >::normalize_row (MR< T > &A,
										 lidia_size_t startr,
										 lidia_size_t startc,
										 lidia_size_t index,
										 lidia_size_t len) const
{
	return elim_modul.normalize_row(A, startr, startc, index, len);
}



template< class T, class REP, class REP1 >
inline bool
havas_min_abs_of_row_plus_minimal_euclidean_norm< T, REP, REP1 >::normalize_row (MR< T > &A,
										 matrix< bigint > &TR,
										 lidia_size_t startr,
										 lidia_size_t startc,
										 lidia_size_t index,
										 lidia_size_t len) const
{
	return elim_modul.normalize_row(A, TR, startr, startc, index, len);
}



//
// Hermite normal form configuration: min_abs_of_row + minimal_norm
//

template< class T, class REP, class REP1 >
inline
havas_min_abs_of_row_plus_minimal_norm< T, REP, REP1 >::havas_min_abs_of_row_plus_minimal_norm ()
{
}



template< class T, class REP, class REP1 >
inline
havas_min_abs_of_row_plus_minimal_norm< T, REP, REP1 >::~havas_min_abs_of_row_plus_minimal_norm ()
{
}



//
// init
//

template< class T, class REP, class REP1 >
inline void
havas_min_abs_of_row_plus_minimal_norm< T, REP, REP1 >::init (MR< T > &, const T &)
{
}



//
// mgcd
//

template< class T, class REP, class REP1 >
inline lidia_size_t
havas_min_abs_of_row_plus_minimal_norm< T, REP, REP1 >::mgcd (MR< T > &A,
							      lidia_size_t startr,
							      lidia_size_t startc) const
{
	lidia_size_t index, count;
	while ((count = pivot_modul.min_abs_of_row_plus_minimal_norm(A, startr, startc, index, 1)) > 1)
		if (!elim_modul.normalize_row(A, startr, startc, index, startr))
			return -1;

	if (count == 0)
		return 0;

	if (index != startc)
		modul.swap_columns(A, startc, index);
	return 1;

}



template< class T, class REP, class REP1 >
inline lidia_size_t
havas_min_abs_of_row_plus_minimal_norm< T, REP, REP1 >::mgcd (MR< T > &A,
							      matrix< bigint > &TR,
							      lidia_size_t startr,
							      lidia_size_t startc) const
{
	lidia_size_t index, count;
	while ((count = pivot_modul.min_abs_of_row_plus_minimal_norm(A, startr, startc, index, 1)) > 1)
		if (!elim_modul.normalize_row(A, TR, startr, startc, index, startr))
			return -1;

	if (count == 0)
		return 0;

	if (index != startc) {
		modul.swap_columns(A, startc, index);
		TR.swap_columns(startc, index);
	}
	return 1;
}



//
// normalization
//

template< class T, class REP, class REP1 >
inline bool
havas_min_abs_of_row_plus_minimal_norm< T, REP, REP1 >::normalize_row (MR< T > &A,
								       lidia_size_t startr,
								       lidia_size_t startc,
								       lidia_size_t index,
								       lidia_size_t len) const
{
	return elim_modul.normalize_row(A, startr, startc, index, len);
}



template< class T, class REP, class REP1 >
inline bool
havas_min_abs_of_row_plus_minimal_norm< T, REP, REP1 >::normalize_row (MR< T > &A,
								       matrix< bigint > &TR,
								       lidia_size_t startr,
								       lidia_size_t startc,
								       lidia_size_t index,
								       lidia_size_t len) const
{
	return elim_modul.normalize_row(A, TR, startr, startc, index, len);
}



//
// Hermite normal form configuration: min_abs_of_row + min_no_of_elements
//

template< class T, class REP, class REP1 >
inline
havas_min_abs_of_row_plus_min_no_of_elements< T, REP, REP1 >::havas_min_abs_of_row_plus_min_no_of_elements ()
{
}



template< class T, class REP, class REP1 >
inline
havas_min_abs_of_row_plus_min_no_of_elements< T, REP, REP1 >::~havas_min_abs_of_row_plus_min_no_of_elements ()
{
}



//
// init
//

template< class T, class REP, class REP1 >
inline void
havas_min_abs_of_row_plus_min_no_of_elements< T, REP, REP1 >::init (MR< T > &, const T &)
{
}



//
// mgcd
//

template< class T, class REP, class REP1 >
inline lidia_size_t
havas_min_abs_of_row_plus_min_no_of_elements< T, REP, REP1 >::mgcd (MR< T > &A,
								    lidia_size_t startr,
								    lidia_size_t startc) const
{
	lidia_size_t index, count;
	while ((count = pivot_modul.min_abs_of_row_plus_minimal_norm(A, startr, startc, index, 0)) > 1)
		if (!elim_modul.normalize_row(A, startr, startc, index, startr))
			return -1;

	if (count == 0)
		return 0;

	if (index != startc)
		modul.swap_columns(A, startc, index);
	return 1;

}



template< class T, class REP, class REP1 >
inline lidia_size_t
havas_min_abs_of_row_plus_min_no_of_elements< T, REP, REP1 >::mgcd (MR< T > &A,
								    matrix< bigint > &TR,
								    lidia_size_t startr,
								    lidia_size_t startc) const
{
	lidia_size_t index, count;
	while ((count = pivot_modul.min_abs_of_row_plus_minimal_norm(A, startr, startc, index, 0)) > 1)
		if (!elim_modul.normalize_row(A, TR, startr, startc, index, startr))
			return -1;

	if (count == 0)
		return 0;

	if (index != startc) {
		modul.swap_columns(A, startc, index);
		TR.swap_columns(startc, index);
	}
	return 1;
}



//
// normalization
//

template< class T, class REP, class REP1 >
inline bool
havas_min_abs_of_row_plus_min_no_of_elements< T, REP, REP1 >::normalize_row (MR< T > &A,
									     lidia_size_t startr,
									     lidia_size_t startc,
									     lidia_size_t index,
									     lidia_size_t len) const
{
	return elim_modul.normalize_row(A, startr, startc, index, len);
}



template< class T, class REP, class REP1 >
inline bool
havas_min_abs_of_row_plus_min_no_of_elements< T, REP, REP1 >::normalize_row (MR< T > &A,
									     matrix< bigint > &TR,
									     lidia_size_t startr,
									     lidia_size_t startc,
									     lidia_size_t index,
									     lidia_size_t len) const
{
	return elim_modul.normalize_row(A, TR, startr, startc, index, len);
}



//
// Hermite normal form configuration: minimal_norm + min_abs_of_row
//

template< class T, class REP, class REP1 >
inline
havas_minimal_norm_plus_min_abs_of_row< T, REP, REP1 >::havas_minimal_norm_plus_min_abs_of_row ()
{
}



template< class T, class REP, class REP1 >
inline
havas_minimal_norm_plus_min_abs_of_row< T, REP, REP1 >::~havas_minimal_norm_plus_min_abs_of_row ()
{
}



//
// init
//

template< class T, class REP, class REP1 >
inline void
havas_minimal_norm_plus_min_abs_of_row< T, REP, REP1 >::init (MR< T > &, const T &)
{
}



//
// mgcd
//

template< class T, class REP, class REP1 >
inline lidia_size_t
havas_minimal_norm_plus_min_abs_of_row< T, REP, REP1 >::mgcd (MR< T > &A,
							      lidia_size_t startr,
							      lidia_size_t startc) const
{
	lidia_size_t index, count;
	while ((count = pivot_modul.minimal_norm_plus_min_abs_of_row(A, startr, startc, index, 1)) > 1)
		if (!elim_modul.normalize_row(A, startr, startc, index, startr))
			return -1;

	if (count == 0)
		return 0;

	if (index != startc)
		modul.swap_columns(A, startc, index);
	return 1;

}



template< class T, class REP, class REP1 >
inline lidia_size_t
havas_minimal_norm_plus_min_abs_of_row< T, REP, REP1 >::mgcd (MR< T > &A,
							      matrix< bigint > &TR,
							      lidia_size_t startr,
							      lidia_size_t startc) const
{
	lidia_size_t index, count;
	while ((count = pivot_modul.minimal_norm_plus_min_abs_of_row(A, startr, startc, index, 1)) > 1)
		if (!elim_modul.normalize_row(A, TR, startr, startc, index, startr))
			return -1;

	if (count == 0)
		return 0;

	if (index != startc) {
		modul.swap_columns(A, startc, index);
		TR.swap_columns(startc, index);
	}
	return 1;
}



//
// normalization
//

template< class T, class REP, class REP1 >
inline bool
havas_minimal_norm_plus_min_abs_of_row< T, REP, REP1 >::normalize_row (MR< T > &A,
								       lidia_size_t startr,
								       lidia_size_t startc,
								       lidia_size_t index,
								       lidia_size_t len) const
{
	return elim_modul.normalize_row(A, startr, startc, index, len);
}



template< class T, class REP, class REP1 >
inline bool
havas_minimal_norm_plus_min_abs_of_row< T, REP, REP1 >::normalize_row (MR< T > &A,
								       matrix< bigint > &TR,
								       lidia_size_t startr,
								       lidia_size_t startc,
								       lidia_size_t index,
								       lidia_size_t len) const
{
	return elim_modul.normalize_row(A, TR, startr, startc, index, len);
}



//
// Hermite normal form configuration: minimal_norm + sorting_gcd
//

template< class T, class REP, class REP1 >
inline
havas_minimal_norm_plus_sorting_gcd< T, REP, REP1 >::havas_minimal_norm_plus_sorting_gcd ()
{
}



template< class T, class REP, class REP1 >
inline
havas_minimal_norm_plus_sorting_gcd< T, REP, REP1 >::~havas_minimal_norm_plus_sorting_gcd ()
{
}



//
// init
//

template< class T, class REP, class REP1 >
inline void
havas_minimal_norm_plus_sorting_gcd< T, REP, REP1 >::init (MR< T > &, const T &)
{
}



//
// mgcd
//

template< class T, class REP, class REP1 >
inline lidia_size_t
havas_minimal_norm_plus_sorting_gcd< T, REP, REP1 >::mgcd (MR< T > &A,
							   lidia_size_t startr,
							   lidia_size_t startc) const
{
	lidia_size_t index, count;
	while ((count = pivot_modul.minimal_norm_plus_sorting_gcd(A, startr, startc, index, 1)) > 1)
		if (!elim_modul.normalize_row(A, startr, startc, index, startr))
			return -1;

	if (count == 0)
		return 0;

	if (index != startc)
		modul.swap_columns(A, startc, index);
	return 1;

}



template< class T, class REP, class REP1 >
inline lidia_size_t
havas_minimal_norm_plus_sorting_gcd< T, REP, REP1 >::mgcd (MR< T > &A,
							   matrix< bigint > &TR,
							   lidia_size_t startr,
							   lidia_size_t startc) const
{
	lidia_size_t index, count;
	while ((count = pivot_modul.minimal_norm_plus_sorting_gcd(A, startr, startc, index, 1)) > 1)
		if (!elim_modul.normalize_row(A, TR, startr, startc, index, startr))
			return -1;

	if (count == 0)
		return 0;

	if (index != startc) {
		modul.swap_columns(A, startc, index);
		TR.swap_columns(startc, index);
	}
	return 1;
}



//
// normalization
//

template< class T, class REP, class REP1 >
inline bool
havas_minimal_norm_plus_sorting_gcd< T, REP, REP1 >::normalize_row (MR< T > &A,
								    lidia_size_t startr,
								    lidia_size_t startc,
								    lidia_size_t index,
								    lidia_size_t len) const
{
	return elim_modul.normalize_row(A, startr, startc, index, len);
}



template< class T, class REP, class REP1 >
inline bool
havas_minimal_norm_plus_sorting_gcd< T, REP, REP1 >::normalize_row (MR< T > &A,
								    matrix< bigint > &TR,
								    lidia_size_t startr,
								    lidia_size_t startc,
								    lidia_size_t index,
								    lidia_size_t len) const
{
	return elim_modul.normalize_row(A, TR, startr, startc, index, len);
}



//
// Hermite normal form configuration: minimal_norm + min_no_of_elements
//

template< class T, class REP, class REP1 >
inline
havas_minimal_norm_plus_min_no_of_elements< T, REP, REP1 >::havas_minimal_norm_plus_min_no_of_elements ()
{
}



template< class T, class REP, class REP1 >
inline
havas_minimal_norm_plus_min_no_of_elements< T, REP, REP1 >::~havas_minimal_norm_plus_min_no_of_elements ()
{
}



//
// init
//

template< class T, class REP, class REP1 >
inline void
havas_minimal_norm_plus_min_no_of_elements< T, REP, REP1 >::init (MR< T > &, const T &)
{
}



//
// mgcd
//

template< class T, class REP, class REP1 >
inline lidia_size_t
havas_minimal_norm_plus_min_no_of_elements< T, REP, REP1 >::mgcd (MR< T > &A,
								  lidia_size_t startr,
								  lidia_size_t startc) const
{
	lidia_size_t index, count;
	while ((count = pivot_modul.minimal_norm_plus_min_no_of_elements(A, startr, startc, index, 1)) > 1)
		if (!elim_modul.normalize_row(A, startr, startc, index, startr))
			return -1;

	if (count == 0)
		return 0;

	if (index != startc)
		modul.swap_columns(A, startc, index);
	return 1;

}



template< class T, class REP, class REP1 >
inline lidia_size_t
havas_minimal_norm_plus_min_no_of_elements< T, REP, REP1 >::mgcd (MR< T > &A,
								  matrix< bigint > &TR,
								  lidia_size_t startr,
								  lidia_size_t startc) const
{
	lidia_size_t index, count;
	while ((count = pivot_modul.minimal_norm_plus_min_no_of_elements(A, startr, startc, index, 1)) > 1)
		if (!elim_modul.normalize_row(A, TR, startr, startc, index, startr))
			return -1;

	if (count == 0)
		return 0;

	if (index != startc) {
		modul.swap_columns(A, startc, index);
		TR.swap_columns(startc, index);
	}
	return 1;
}



//
// normalization
//

template< class T, class REP, class REP1 >
inline bool
havas_minimal_norm_plus_min_no_of_elements< T, REP, REP1 >::normalize_row (MR< T > &A,
									   lidia_size_t startr,
									   lidia_size_t startc,
									   lidia_size_t index,
									   lidia_size_t len) const
{
	return elim_modul.normalize_row(A, startr, startc, index, len);
}



template< class T, class REP, class REP1 >
inline bool
havas_minimal_norm_plus_min_no_of_elements< T, REP, REP1 >::normalize_row (MR< T > &A,
									   matrix< bigint > &TR,
									   lidia_size_t startr,
									   lidia_size_t startc,
									   lidia_size_t index,
									   lidia_size_t len) const
{
	return elim_modul.normalize_row(A, TR, startr, startc, index, len);
}



//
// Hermite normal form configuration: sorting_gcd + minimal_euclidean_norm
//

template< class T, class REP, class REP1 >
inline
havas_sorting_gcd_plus_minimal_euclidean_norm< T, REP, REP1 >::havas_sorting_gcd_plus_minimal_euclidean_norm ()
{
}



template< class T, class REP, class REP1 >
inline
havas_sorting_gcd_plus_minimal_euclidean_norm< T, REP, REP1 >::~havas_sorting_gcd_plus_minimal_euclidean_norm ()
{
}



//
// init
//

template< class T, class REP, class REP1 >
inline void
havas_sorting_gcd_plus_minimal_euclidean_norm< T, REP, REP1 >::init (MR< T > &, const T &)
{
}



//
// mgcd
//

template< class T, class REP, class REP1 >
inline lidia_size_t
havas_sorting_gcd_plus_minimal_euclidean_norm< T, REP, REP1 >::mgcd (MR< T > &A,
								     lidia_size_t startr,
								     lidia_size_t startc) const
{
	lidia_size_t index, count;
	while ((count = pivot_modul.sorting_gcd_plus_minimal_norm(A, startr, startc, index, 2)) > 1)
		if (!elim_modul.normalize_row(A, startr, startc, index, startr))
			return -1;

	if (count == 0)
		return 0;

	if (index != startc)
		modul.swap_columns(A, startc, index);
	return 1;

}



template< class T, class REP, class REP1 >
inline lidia_size_t
havas_sorting_gcd_plus_minimal_euclidean_norm< T, REP, REP1 >::mgcd (MR< T > &A,
								     matrix< bigint > &TR,
								     lidia_size_t startr,
								     lidia_size_t startc) const
{
	lidia_size_t index, count;
	while ((count = pivot_modul.sorting_gcd_plus_minimal_norm(A, startr, startc, index, 2)) > 1)
		if (!elim_modul.normalize_row(A, TR, startr, startc, index, startr))
			return -1;

	if (count == 0)
		return 0;

	if (index != startc) {
		modul.swap_columns(A, startc, index);
		TR.swap_columns(startc, index);
	}
	return 1;
}



//
// normalization
//

template< class T, class REP, class REP1 >
inline bool
havas_sorting_gcd_plus_minimal_euclidean_norm< T, REP, REP1 >::normalize_row (MR< T > &A,
									      lidia_size_t startr,
									      lidia_size_t startc,
									      lidia_size_t index,
									      lidia_size_t len) const
{
	return elim_modul.normalize_row(A, startr, startc, index, len);
}



template< class T, class REP, class REP1 >
inline bool
havas_sorting_gcd_plus_minimal_euclidean_norm< T, REP, REP1 >::normalize_row (MR< T > &A,
									      matrix< bigint > &TR,
									      lidia_size_t startr,
									      lidia_size_t startc,
									      lidia_size_t index,
									      lidia_size_t len) const
{
	return elim_modul.normalize_row(A, TR, startr, startc, index, len);
}



//
// Hermite normal form configuration: sorting_gcd + minimal_norm
//

template< class T, class REP, class REP1 >
inline
havas_sorting_gcd_plus_minimal_norm< T, REP, REP1 >::havas_sorting_gcd_plus_minimal_norm ()
{
}



template< class T, class REP, class REP1 >
inline
havas_sorting_gcd_plus_minimal_norm< T, REP, REP1 >::~havas_sorting_gcd_plus_minimal_norm ()
{
}



//
// init
//

template< class T, class REP, class REP1 >
inline void
havas_sorting_gcd_plus_minimal_norm< T, REP, REP1 >::init (MR< T > &, const T &)
{
}



//
// mgcd
//

template< class T, class REP, class REP1 >
inline lidia_size_t
havas_sorting_gcd_plus_minimal_norm< T, REP, REP1 >::mgcd (MR< T > &A,
							   lidia_size_t startr,
							   lidia_size_t startc) const
{
	lidia_size_t index, count;
	while ((count = pivot_modul.sorting_gcd_plus_minimal_norm(A, startr, startc, index, 1)) > 1)
		if (!elim_modul.normalize_row(A, startr, startc, index, startr))
			return -1;

	if (count == 0)
		return 0;

	if (index != startc)
		modul.swap_columns(A, startc, index);
	return 1;

}



template< class T, class REP, class REP1 >
inline lidia_size_t
havas_sorting_gcd_plus_minimal_norm< T, REP, REP1 >::mgcd (MR< T > &A,
							   matrix< bigint > &TR,
							   lidia_size_t startr,
							   lidia_size_t startc) const
{
	lidia_size_t index, count;
	while ((count = pivot_modul.sorting_gcd_plus_minimal_norm(A, startr, startc, index, 1)) > 1)
		if (!elim_modul.normalize_row(A, TR, startr, startc, index, startr))
			return -1;

	if (count == 0)
		return 0;

	if (index != startc) {
		modul.swap_columns(A, startc, index);
		TR.swap_columns(startc, index);
	}
	return 1;
}



//
// normalization
//

template< class T, class REP, class REP1 >
inline bool
havas_sorting_gcd_plus_minimal_norm< T, REP, REP1 >::normalize_row (MR< T > &A,
								    lidia_size_t startr,
								    lidia_size_t startc,
								    lidia_size_t index,
								    lidia_size_t len) const
{
	return elim_modul.normalize_row(A, startr, startc, index, len);
}



template< class T, class REP, class REP1 >
inline bool
havas_sorting_gcd_plus_minimal_norm< T, REP, REP1 >::normalize_row (MR< T > &A,
								    matrix< bigint > &TR,
								    lidia_size_t startr,
								    lidia_size_t startc,
								    lidia_size_t index,
								    lidia_size_t len) const
{
	return elim_modul.normalize_row(A, TR, startr, startc, index, len);
}



//
// Hermite normal form configuration: sorting_gcd + min_no_of_elements
//

template< class T, class REP, class REP1 >
inline
havas_sorting_gcd_plus_min_no_of_elements< T, REP, REP1 >::havas_sorting_gcd_plus_min_no_of_elements ()
{
}



template< class T, class REP, class REP1 >
inline
havas_sorting_gcd_plus_min_no_of_elements< T, REP, REP1 >::~havas_sorting_gcd_plus_min_no_of_elements ()
{
}



//
// init
//

template< class T, class REP, class REP1 >
inline void
havas_sorting_gcd_plus_min_no_of_elements< T, REP, REP1 >::init (MR< T > &, const T &)
{
}



//
// mgcd
//

template< class T, class REP, class REP1 >
inline lidia_size_t
havas_sorting_gcd_plus_min_no_of_elements< T, REP, REP1 >::mgcd (MR< T > &A,
								 lidia_size_t startr,
								 lidia_size_t startc) const
{
	lidia_size_t index, count;
	while ((count = pivot_modul.sorting_gcd_plus_minimal_norm(A, startr, startc, index, 0)) > 1)
		if (!elim_modul.normalize_row(A, startr, startc, index, startr))
			return -1;

	if (count == 0)
		return 0;

	if (index != startc)
		modul.swap_columns(A, startc, index);
	return 1;

}



template< class T, class REP, class REP1 >
inline lidia_size_t
havas_sorting_gcd_plus_min_no_of_elements< T, REP, REP1 >::mgcd (MR< T > &A,
								 matrix< bigint > &TR,
								 lidia_size_t startr,
								 lidia_size_t startc) const
{
	lidia_size_t index, count;
	while ((count = pivot_modul.sorting_gcd_plus_minimal_norm(A, startr, startc, index, 0)) > 1)
		if (!elim_modul.normalize_row(A, TR, startr, startc, index, startr))
			return -1;

	if (count == 0)
		return 0;

	if (index != startc) {
		modul.swap_columns(A, startc, index);
		TR.swap_columns(startc, index);
	}
	return 1;
}



//
// normalization
//

template< class T, class REP, class REP1 >
inline bool
havas_sorting_gcd_plus_min_no_of_elements< T, REP, REP1 >::normalize_row (MR< T > &A,
									  lidia_size_t startr,
									  lidia_size_t startc,
									  lidia_size_t index,
									  lidia_size_t len) const
{
	return elim_modul.normalize_row(A, startr, startc, index, len);
}



template< class T, class REP, class REP1 >
inline bool
havas_sorting_gcd_plus_min_no_of_elements< T, REP, REP1 >::normalize_row (MR< T > &A,
									  matrix< bigint > &TR,
									  lidia_size_t startr,
									  lidia_size_t startc,
									  lidia_size_t index,
									  lidia_size_t len) const
{
	return elim_modul.normalize_row(A, TR, startr, startc, index, len);
}



//
// Hermite normal form configuration: min_no_of_elements + min_abs_of_row
//

template< class T, class REP, class REP1 >
inline
havas_min_no_of_elements_plus_min_abs_of_row< T, REP, REP1 >::havas_min_no_of_elements_plus_min_abs_of_row ()
{
}



template< class T, class REP, class REP1 >
inline
havas_min_no_of_elements_plus_min_abs_of_row< T, REP, REP1 >::~havas_min_no_of_elements_plus_min_abs_of_row ()
{
}



//
// init
//

template< class T, class REP, class REP1 >
inline void
havas_min_no_of_elements_plus_min_abs_of_row< T, REP, REP1 >::init (MR< T > &, const T &)
{
}



//
// mgcd
//

template< class T, class REP, class REP1 >
inline lidia_size_t
havas_min_no_of_elements_plus_min_abs_of_row< T, REP, REP1 >::mgcd (MR< T > &A,
								    lidia_size_t startr,
								    lidia_size_t startc) const
{
	lidia_size_t index, count;
	while ((count = pivot_modul.minimal_norm_plus_min_abs_of_row(A, startr, startc, index, 0)) > 1)
		if (!elim_modul.normalize_row(A, startr, startc, index, startr))
			return -1;

	if (count == 0)
		return 0;

	if (index != startc)
		modul.swap_columns(A, startc, index);
	return 1;

}



template< class T, class REP, class REP1 >
inline lidia_size_t
havas_min_no_of_elements_plus_min_abs_of_row< T, REP, REP1 >::mgcd (MR< T > &A,
								    matrix< bigint > &TR,
								    lidia_size_t startr,
								    lidia_size_t startc) const
{
	lidia_size_t index, count;
	while ((count = pivot_modul.minimal_norm_plus_min_abs_of_row(A, startr, startc, index, 0)) > 1)
		if (!elim_modul.normalize_row(A, TR, startr, startc, index, startr))
			return -1;

	if (count == 0)
		return 0;

	if (index != startc) {
		modul.swap_columns(A, startc, index);
		TR.swap_columns(startc, index);
	}
	return 1;
}



//
// normalization
//

template< class T, class REP, class REP1 >
inline bool
havas_min_no_of_elements_plus_min_abs_of_row< T, REP, REP1 >::normalize_row (MR< T > &A,
									     lidia_size_t startr,
									     lidia_size_t startc,
									     lidia_size_t index,
									     lidia_size_t len) const
{
	return elim_modul.normalize_row(A, startr, startc, index, len);
}



template< class T, class REP, class REP1 >
inline bool
havas_min_no_of_elements_plus_min_abs_of_row< T, REP, REP1 >::normalize_row (MR< T > &A,
									     matrix< bigint > &TR,
									     lidia_size_t startr,
									     lidia_size_t startc,
									     lidia_size_t index,
									     lidia_size_t len) const
{
	return elim_modul.normalize_row(A, TR, startr, startc, index, len);
}



//
// Hermite normal form configuration: min_no_of_elements + minimal_euclidean_norm
//

template< class T, class REP, class REP1 >
inline
havas_min_no_of_elements_plus_minimal_euclidean_norm< T, REP, REP1 >::havas_min_no_of_elements_plus_minimal_euclidean_norm ()
{
}



template< class T, class REP, class REP1 >
inline
havas_min_no_of_elements_plus_minimal_euclidean_norm< T, REP, REP1 >::~havas_min_no_of_elements_plus_minimal_euclidean_norm ()
{
}



//
// init
//

template< class T, class REP, class REP1 >
inline void
havas_min_no_of_elements_plus_minimal_euclidean_norm< T, REP, REP1 >::init (MR< T > &, const T &)
{
}



//
// mgcd
//

template< class T, class REP, class REP1 >
inline lidia_size_t
havas_min_no_of_elements_plus_minimal_euclidean_norm< T, REP, REP1 >::mgcd (MR< T > &A,
									    lidia_size_t startr,
									    lidia_size_t startc) const
{
	lidia_size_t index, count;
	while ((count = pivot_modul.min_no_of_elements_plus_minimal_norm(A, startr, startc, index, 2)) > 1)
		if (!elim_modul.normalize_row(A, startr, startc, index, startr))
			return -1;

	if (count == 0)
		return 0;

	if (index != startc)
		modul.swap_columns(A, startc, index);
	return 1;

}



template< class T, class REP, class REP1 >
inline lidia_size_t
havas_min_no_of_elements_plus_minimal_euclidean_norm< T, REP, REP1 >::mgcd (MR< T > &A,
									    matrix< bigint > &TR,
									    lidia_size_t startr,
									    lidia_size_t startc) const
{
	lidia_size_t index, count;
	while ((count = pivot_modul.min_no_of_elements_plus_minimal_norm(A, startr, startc, index, 2)) > 1)
		if (!elim_modul.normalize_row(A, TR, startr, startc, index, startr))
			return -1;

	if (count == 0)
		return 0;

	if (index != startc) {
		modul.swap_columns(A, startc, index);
		TR.swap_columns(startc, index);
	}
	return 1;
}



//
// normalization
//

template< class T, class REP, class REP1 >
inline bool
havas_min_no_of_elements_plus_minimal_euclidean_norm< T, REP, REP1 >::normalize_row (MR< T > &A,
										     lidia_size_t startr,
										     lidia_size_t startc,
										     lidia_size_t index,
										     lidia_size_t len) const
{
	return elim_modul.normalize_row(A, startr, startc, index, len);
}



template< class T, class REP, class REP1 >
inline bool
havas_min_no_of_elements_plus_minimal_euclidean_norm< T, REP, REP1 >::normalize_row (MR< T > &A,
										     matrix< bigint > &TR,
										     lidia_size_t startr,
										     lidia_size_t startc,
										     lidia_size_t index,
										     lidia_size_t len) const
{
	return elim_modul.normalize_row(A, TR, startr, startc, index, len);
}



//
// Hermite normal form configuration: min_no_of_elements + minimal_norm
//

template< class T, class REP, class REP1 >
inline
havas_min_no_of_elements_plus_minimal_norm< T, REP, REP1 >::havas_min_no_of_elements_plus_minimal_norm ()
{
}



template< class T, class REP, class REP1 >
inline
havas_min_no_of_elements_plus_minimal_norm< T, REP, REP1 >::~havas_min_no_of_elements_plus_minimal_norm ()
{
}



//
// init
//

template< class T, class REP, class REP1 >
inline void
havas_min_no_of_elements_plus_minimal_norm< T, REP, REP1 >::init (MR< T > &, const T &)
{
}



//
// mgcd
//

template< class T, class REP, class REP1 >
inline lidia_size_t
havas_min_no_of_elements_plus_minimal_norm< T, REP, REP1 >::mgcd (MR< T > &A,
								  lidia_size_t startr,
								  lidia_size_t startc) const
{
	lidia_size_t index, count;
	while ((count = pivot_modul.min_no_of_elements_plus_minimal_norm(A, startr, startc, index, 1)) > 1)
		if (!elim_modul.normalize_row(A, startr, startc, index, startr))
			return -1;

	if (count == 0)
		return 0;

	if (index != startc)
		modul.swap_columns(A, startc, index);
	return 1;

}



template< class T, class REP, class REP1 >
inline lidia_size_t
havas_min_no_of_elements_plus_minimal_norm< T, REP, REP1 >::mgcd (MR< T > &A,
								  matrix< bigint > &TR,
								  lidia_size_t startr,
								  lidia_size_t startc) const
{
	lidia_size_t index, count;
	while ((count = pivot_modul.min_no_of_elements_plus_minimal_norm(A, startr, startc, index, 1)) > 1)
		if (!elim_modul.normalize_row(A, TR, startr, startc, index, startr))
			return -1;

	if (count == 0)
		return 0;

	if (index != startc) {
		modul.swap_columns(A, startc, index);
		TR.swap_columns(startc, index);
	}
	return 1;
}



//
// normalization
//

template< class T, class REP, class REP1 >
inline bool
havas_min_no_of_elements_plus_minimal_norm< T, REP, REP1 >::normalize_row (MR< T > &A,
									   lidia_size_t startr,
									   lidia_size_t startc,
									   lidia_size_t index,
									   lidia_size_t len) const
{
	return elim_modul.normalize_row(A, startr, startc, index, len);
}



template< class T, class REP, class REP1 >
inline bool
havas_min_no_of_elements_plus_minimal_norm< T, REP, REP1 >::normalize_row (MR< T > &A,
									   matrix< bigint > &TR,
									   lidia_size_t startr,
									   lidia_size_t startc,
									   lidia_size_t index,
									   lidia_size_t len) const
{
	return elim_modul.normalize_row(A, TR, startr, startc, index, len);
}



//
// Hermite normal form configuration: min_no_of_elements + sorting_gcd
//

template< class T, class REP, class REP1 >
inline
havas_min_no_of_elements_plus_sorting_gcd< T, REP, REP1 >::havas_min_no_of_elements_plus_sorting_gcd ()
{
}



template< class T, class REP, class REP1 >
inline
havas_min_no_of_elements_plus_sorting_gcd< T, REP, REP1 >::~havas_min_no_of_elements_plus_sorting_gcd ()
{
}



//
// init
//

template< class T, class REP, class REP1 >
inline void
havas_min_no_of_elements_plus_sorting_gcd< T, REP, REP1 >::init (MR< T > &, const T &)
{
}



//
// mgcd
//

template< class T, class REP, class REP1 >
inline lidia_size_t
havas_min_no_of_elements_plus_sorting_gcd< T, REP, REP1 >::mgcd (MR< T > &A,
								 lidia_size_t startr,
								 lidia_size_t startc) const
{
	lidia_size_t index = 0, count;
	while ((count = pivot_modul.minimal_norm_plus_sorting_gcd(A, startr, startc, index, 0)) > 1)
		if (!elim_modul.normalize_row(A, startr, startc, index, startr))
			return -1;

	if (count == 0)
		return 0;

	if (index != startc)
		modul.swap_columns(A, startc, index);
	return 1;

}



template< class T, class REP, class REP1 >
inline lidia_size_t
havas_min_no_of_elements_plus_sorting_gcd< T, REP, REP1 >::mgcd (MR< T > &A,
								 matrix< bigint > &TR,
								 lidia_size_t startr,
								 lidia_size_t startc) const
{
	lidia_size_t index = 0, count;
	while ((count = pivot_modul.minimal_norm_plus_sorting_gcd(A, startr, startc, index, 0)) > 1)
		if (!elim_modul.normalize_row(A, TR, startr, startc, index, startr))
			return -1;

	if (count == 0)
		return 0;

	if (index != startc) {
		modul.swap_columns(A, startc, index);
		TR.swap_columns(startc, index);
	}
	return 1;
}



//
// normalization
//

template< class T, class REP, class REP1 >
inline bool
havas_min_no_of_elements_plus_sorting_gcd< T, REP, REP1 >::normalize_row (MR< T > &A,
									  lidia_size_t startr,
									  lidia_size_t startc,
									  lidia_size_t index,
									  lidia_size_t len) const
{
	return elim_modul.normalize_row(A, startr, startc, index, len);
}



template< class T, class REP, class REP1 >
inline bool
havas_min_no_of_elements_plus_sorting_gcd< T, REP, REP1 >::normalize_row (MR< T > &A,
									  matrix< bigint > &TR,
									  lidia_size_t startr,
									  lidia_size_t startc,
									  lidia_size_t index,
									  lidia_size_t len) const
{
	return elim_modul.normalize_row(A, TR, startr, startc, index, len);
}



//
// Hermite normal form configuration: Storjohann
//

template< class T, class REP, class REP1 >
inline
Storjohann< T, REP, REP1 >::Storjohann ()
{
}



template< class T, class REP, class REP1 >
inline
Storjohann< T, REP, REP1 >::~Storjohann ()
{
}



//
// init
//

template< class T, class REP, class REP1 >
inline void
Storjohann< T, REP, REP1 >::init (MR< T > &, const T &)
{
}



//
// mgcd
//

template< class T, class REP, class REP1 >
inline lidia_size_t
Storjohann< T, REP, REP1 >::mgcd (MR< T > &A,
				  lidia_size_t startr,
				  lidia_size_t startc) const
{
	if (elim_modul.mgcd_storjohann(A, startr, startc, startr))
		return 1;
	return 0;
}



template< class T, class REP, class REP1 >
inline lidia_size_t
Storjohann< T, REP, REP1 >::mgcd (MR< T > &A,
				  matrix< bigint > &TR,
				  lidia_size_t startr,
				  lidia_size_t startc) const
{
	if (elim_modul.mgcd_storjohann(A, TR, startr, startc, startr))
		return 1;
	return 0;
}



//
// normalization
//

template< class T, class REP, class REP1 >
inline bool
Storjohann< T, REP, REP1 >::normalize_row (MR< T > &A,
					   lidia_size_t startr,
					   lidia_size_t startc,
					   lidia_size_t index,
					   lidia_size_t len) const
{
	return elim_modul.normalize_row(A, startr, startc, index, len);
}



template< class T, class REP, class REP1 >
inline bool
Storjohann< T, REP, REP1 >::normalize_row (MR< T > &A,
					   matrix< bigint > &TR,
					   lidia_size_t startr,
					   lidia_size_t startc,
					   lidia_size_t index,
					   lidia_size_t len) const
{
	return elim_modul.normalize_row(A, TR, startr, startc, index, len);
}



//
// Hermite normal form configuration: Heuristik
//

template< class T, class REP, class REP1 >
inline
Heuristik< T, REP, REP1 >::Heuristik ()
{
}



template< class T, class REP, class REP1 >
inline
Heuristik< T, REP, REP1 >::~Heuristik ()
{
}



//
// init
//

template< class T, class REP, class REP1 >
inline void
Heuristik< T, REP, REP1 >::init (MR< T > &, const T &)
{
}



//
// mgcd
//

template< class T, class REP, class REP1 >
inline lidia_size_t
Heuristik< T, REP, REP1 >::mgcd (MR< T > &A,
				 lidia_size_t startr,
				 lidia_size_t startc) const
{
	lidia_size_t index, count;
	if ((count = pivot_modul.min_abs_of_row(A, startr, startc, index)) > 1)
		if (!elim_modul.xgcd_elimination(A, startr, startc, index, startr))
			return -1;

	if (count == 0)
		return 0;

	while ((count = pivot_modul.min_abs_of_row(A, startr, startc, index)) > 1)
		if (!elim_modul.normalize_row(A, startr, startc, index, startr))
			return -1;

	if (count == 0)
		return 0;

	if (index != startc)
		modul.swap_columns(A, startc, index);
	return 1;
}



template< class T, class REP, class REP1 >
inline lidia_size_t
Heuristik< T, REP, REP1 >::mgcd (MR< T > &A,
				 matrix< bigint > &TR,
				 lidia_size_t startr,
				 lidia_size_t startc) const
{
	lidia_size_t index, count;
	if ((count = pivot_modul.min_abs_of_row(A, startr, startc, index)) > 1)
		if (!elim_modul.xgcd_elimination(A, TR, startr, startc, index, startr))
			return -1;

	if (count == 0)
		return 0;

	while ((count = pivot_modul.min_abs_of_row(A, startr, startc, index)) > 1)
		if (!elim_modul.normalize_row(A, TR, startr, startc, index, startr))
			return -1;

	if (count == 0)
		return 0;

	if (index != startc) {
		modul.swap_columns(A, startc, index);
		TR.swap_columns(startc, index);
	}
	return 1;
}



//
// normalization
//

template< class T, class REP, class REP1 >
inline bool
Heuristik< T, REP, REP1 >::normalize_row (MR< T > &A,
					  lidia_size_t startr,
					  lidia_size_t startc,
					  lidia_size_t index,
					  lidia_size_t len) const
{
	return elim_modul.normalize_row(A, startr, startc, index, len);
}



template< class T, class REP, class REP1 >
inline bool
Heuristik< T, REP, REP1 >::normalize_row (MR< T > &A,
					  matrix< bigint > &TR,
					  lidia_size_t startr,
					  lidia_size_t startc,
					  lidia_size_t index,
					  lidia_size_t len) const
{
	return elim_modul.normalize_row(A, TR, startr, startc, index, len);
}



///////////////////////////////////////////////////////////////////////////////////////////

//
// Hermite normal form configuration: havas_best_remainder_ext
//

template< class T, class REP, class REP1 >
inline
havas_best_remainder_ext< T, REP, REP1 >::havas_best_remainder_ext ()
{
	max_array = NULL;
}



template< class T, class REP, class REP1 >
inline
havas_best_remainder_ext< T, REP, REP1 >::~havas_best_remainder_ext ()
{
	delete[] max_array;
}



//
// init
//

template< class T, class REP, class REP1 >
inline void
havas_best_remainder_ext< T, REP, REP1 >::init (MR< T > &A, const T &B)
{
	BOUND = B;
	max_array = init_modul.init_max_array(A);
}



//
// mgcd
//

template< class T, class REP, class REP1 >
inline lidia_size_t
havas_best_remainder_ext< T, REP, REP1 >::mgcd (MR< T > &A,
						lidia_size_t startr,
						lidia_size_t startc) const
{
	lidia_size_t index, count;
	while ((count = pivot_modul.min_abs_of_row(A, startr, startc, index)) > 1)
		if (!elim_modul.normalize_row(A, startr, startc, index, startr, max_array, BOUND))
			return -1;

	if (count == 0)
		return 0;

	if (index != startc) {
		modul.swap_columns(A, startc, index);
		LiDIA::swap(max_array[startc], max_array[index]);
	}
	return 1;
}



template< class T, class REP, class REP1 >
inline lidia_size_t
havas_best_remainder_ext< T, REP, REP1 >::mgcd (MR< T > &A,
						matrix< bigint > &TR,
						lidia_size_t startr,
						lidia_size_t startc) const
{
	lidia_size_t index, count;
	while ((count = pivot_modul.min_abs_of_row(A, startr, startc, index)) > 1)
		if (!elim_modul.normalize_row(A, TR, startr, startc, index, startr, max_array, BOUND))
			return -1;

	if (count == 0)
		return 0;

	if (index != startc) {
		modul.swap_columns(A, startc, index);
		LiDIA::swap(max_array[startc], max_array[index]);
		TR.swap_columns(startc, index);
	}
	return 1;
}



//
// normalization
//

template< class T, class REP, class REP1 >
inline bool
havas_best_remainder_ext< T, REP, REP1 >::normalize_row (MR< T > &A,
							 lidia_size_t startr,
							 lidia_size_t startc,
							 lidia_size_t index,
							 lidia_size_t len) const
{
	return elim_modul.normalize_row(A, startr, startc, index, len);
}



template< class T, class REP, class REP1 >
inline bool
havas_best_remainder_ext< T, REP, REP1 >::normalize_row (MR< T > &A,
							 matrix< bigint > &TR,
							 lidia_size_t startr,
							 lidia_size_t startc,
							 lidia_size_t index,
							 lidia_size_t len) const
{
	return elim_modul.normalize_row(A, TR, startr, startc, index, len);
}



//
// Hermite normal form configuration: havas_euclidean_norm_ext
//

template< class T, class REP, class REP1 >
inline
havas_euclidean_norm_ext< T, REP, REP1 >::havas_euclidean_norm_ext ()
{
	max_array = NULL;
}



template< class T, class REP, class REP1 >
inline
havas_euclidean_norm_ext< T, REP, REP1 >::~havas_euclidean_norm_ext ()
{
	delete[] max_array;
}



//
// init
//

template< class T, class REP, class REP1 >
inline void
havas_euclidean_norm_ext< T, REP, REP1 >::init (MR< T > &A, const T &B)
{
	BOUND = B;
	max_array = init_modul.init_max_array(A);
}



//
// mgcd
//

template< class T, class REP, class REP1 >
inline lidia_size_t
havas_euclidean_norm_ext< T, REP, REP1 >::mgcd (MR< T > &A,
						lidia_size_t startr,
						lidia_size_t startc) const
{
	lidia_size_t index, count;
	while ((count = pivot_modul.minimal_norm(A, startr, startc, index, 2)) > 1)
		if (!elim_modul.normalize_row(A, startr, startc, index, startr, max_array, BOUND))
			return -1;

	if (count == 0)
		return 0;

	if (index != startc) {
		modul.swap_columns(A, startc, index);
		LiDIA::swap(max_array[startc], max_array[index]);
	}
	return 1;
}



template< class T, class REP, class REP1 >
inline lidia_size_t
havas_euclidean_norm_ext< T, REP, REP1 >::mgcd (MR< T > &A,
						matrix< bigint > &TR,
						lidia_size_t startr,
						lidia_size_t startc) const
{
	lidia_size_t index, count;
	while ((count = pivot_modul.minimal_norm(A, startr, startc, index, 2)) > 1)
		if (!elim_modul.normalize_row(A, TR, startr, startc, index, startr, max_array, BOUND))
			return -1;

	if (count == 0)
		return 0;

	if (index != startc) {
		modul.swap_columns(A, startc, index);
		LiDIA::swap(max_array[startc], max_array[index]);
		TR.swap_columns(startc, index);
	}
	return 1;
}



//
// normalization
//

template< class T, class REP, class REP1 >
inline bool
havas_euclidean_norm_ext< T, REP, REP1 >::normalize_row (MR< T > &A,
							 lidia_size_t startr,
							 lidia_size_t startc,
							 lidia_size_t index,
							 lidia_size_t len) const
{
	return elim_modul.normalize_row(A, startr, startc, index, len);
}



template< class T, class REP, class REP1 >
inline bool
havas_euclidean_norm_ext< T, REP, REP1 >::normalize_row (MR< T > &A,
							 matrix< bigint > &TR,
							 lidia_size_t startr,
							 lidia_size_t startc,
							 lidia_size_t index,
							 lidia_size_t len) const
{
	return elim_modul.normalize_row(A, TR, startr, startc, index, len);
}



//
// havas_kernel: Hermite normal form configuration 3e
//

template< class T, class REP, class REP1 >
inline
nf_conf3e< T, REP, REP1 >::nf_conf3e ()
{
	max_array = NULL;
}



template< class T, class REP, class REP1 >
inline
nf_conf3e< T, REP, REP1 >::~nf_conf3e ()
{
	delete[] max_array;
}



//
// init
//

template< class T, class REP, class REP1 >
inline void
nf_conf3e< T, REP, REP1 >::init (MR< T > &A, const T &B)
{
	BOUND = B;
	max_array = init_modul.init_max_array(A);
}



//
// mgcd
//

template< class T, class REP, class REP1 >
inline lidia_size_t
nf_conf3e< T, REP, REP1 >::mgcd (MR< T > &A,
				 lidia_size_t startr,
				 lidia_size_t startc) const
{
	lidia_size_t index, count;
	while ((count = pivot_modul.min_abs_of_row_plus_minimal_norm(A, startr, startc, index, 2)) > 1)
		if (!elim_modul.normalize_row(A, startr, startc, index, startr, max_array, BOUND))
			return -1;

	if (count == 0)
		return 0;

	if (index != startc) {
		modul.swap_columns(A, startc, index);
		LiDIA::swap(max_array[startc], max_array[index]);
	}
	return 1;

}



template< class T, class REP, class REP1 >
inline lidia_size_t
nf_conf3e< T, REP, REP1 >::mgcd (MR< T > &A,
				 matrix< bigint > &TR,
				 lidia_size_t startr,
				 lidia_size_t startc) const
{
	lidia_size_t index, count;
	while ((count = pivot_modul.min_abs_of_row_plus_minimal_norm(A, startr, startc, index, 2)) > 1)
		if (!elim_modul.normalize_row(A, TR, startr, startc, index, startr, max_array, BOUND))
			return -1;

	if (count == 0)
		return 0;

	if (index != startc) {
		modul.swap_columns(A, startc, index);
		LiDIA::swap(max_array[startc], max_array[index]);
		TR.swap_columns(startc, index);
	}
	return 1;
}



//
// normalization
//

template< class T, class REP, class REP1 >
inline bool
nf_conf3e< T, REP, REP1 >::normalize_row (MR< T > &A,
					  lidia_size_t startr,
					  lidia_size_t startc,
					  lidia_size_t index,
					  lidia_size_t len) const
{
	return elim_modul.normalize_row(A, startr, startc, index, len);
}



template< class T, class REP, class REP1 >
inline bool
nf_conf3e< T, REP, REP1 >::normalize_row (MR< T > &A,
					  matrix< bigint > &TR,
					  lidia_size_t startr,
					  lidia_size_t startc,
					  lidia_size_t index,
					  lidia_size_t len) const
{
	return elim_modul.normalize_row(A, TR, startr, startc, index, len);
}



//
// kannan_kernel: Hermite normal form configuration 1
//

template< class T, class REP, class REP1 >
inline
Standard_normalization< T, REP, REP1 >::Standard_normalization ()
{
}



template< class T, class REP, class REP1 >
inline
Standard_normalization< T, REP, REP1 >::~Standard_normalization ()
{
}



//
// normalize
//

template< class T, class REP, class REP1 >
inline lidia_size_t *
Standard_normalization< T, REP, REP1 >::normalize (MR< T > &A,
						   lidia_size_t x,
						   lidia_size_t y) const
{
	return modul.normalize_Std(A, x, y);
}



template< class T, class REP, class REP1 >
inline lidia_size_t *
Standard_normalization< T, REP, REP1 >::normalize (MR< T > &A,
						   matrix< bigint > &TR,
						   lidia_size_t x,
						   lidia_size_t y) const
{
	return modul.normalize_Std(A, TR, x, y);
}



//
// kannan_kernel: Hermite normal form configuration 2
//

template< class T, class REP, class REP1 >
inline
ChouCollins_normalization< T, REP, REP1 >::ChouCollins_normalization ()
{
}



template< class T, class REP, class REP1 >
inline
ChouCollins_normalization< T, REP, REP1 >::~ChouCollins_normalization ()
{
}



//
// normalize
//

template< class T, class REP, class REP1 >
inline lidia_size_t *
ChouCollins_normalization< T, REP, REP1 >::normalize (MR< T > &A,
						      lidia_size_t x,
						      lidia_size_t y) const
{
	return modul.normalize_ChouCollins(A, x, y);
}



template< class T, class REP, class REP1 >
inline lidia_size_t *
ChouCollins_normalization< T, REP, REP1 >::normalize (MR< T > &A,
						      matrix< bigint > &TR,
						      lidia_size_t x,
						      lidia_size_t y) const
{
	return modul.normalize_ChouCollins(A, TR, x, y);
}



//
// kannan_kernel: Hermite normal form configuration 3
//

template< class T, class REP, class REP1 >
inline
Jacobson_normalization< T, REP, REP1 >::Jacobson_normalization ()
{
}



template< class T, class REP, class REP1 >
inline
Jacobson_normalization< T, REP, REP1 >::~Jacobson_normalization ()
{
}



//
// normalize
//

template< class T, class REP, class REP1 >
inline lidia_size_t *
Jacobson_normalization< T, REP, REP1 >::normalize (MR< T > &A,
						   lidia_size_t x,
						   lidia_size_t y) const
{
	return modul.normalizeHybrid_Std(A, x, y);
}



template< class T, class REP, class REP1 >
inline lidia_size_t *
Jacobson_normalization< T, REP, REP1 >::normalize (MR< T > &A,
						   matrix< bigint > &TR,
						   lidia_size_t x,
						   lidia_size_t y) const
{
	return modul.normalizeHybrid_Std(A, TR, x, y);
}



#ifdef LIDIA_NAMESPACE
# ifndef IN_NAMESPACE_LIDIA
}	// end of namespace LiDIA
# endif
#endif



#endif	// LIDIA_HNF_CONF_CC_GUARD_
