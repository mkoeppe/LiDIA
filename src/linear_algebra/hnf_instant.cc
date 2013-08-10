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


#ifdef HAVE_CONFIG_H
# include	"config.h"
#endif
#include	"LiDIA/bigint_matrix.h"
#include	"LiDIA/matrix/hnf_conf.h"
#include	"LiDIA/matrix/hnf_conf.cc"
#include	"LiDIA/matrix/hnf_kernel.h"
#include	"LiDIA/matrix/hnf_kernel.cc"
#include	"LiDIA/matrix/normalize_kernel.h"
#include	"LiDIA/matrix/normalize_kernel.cc"
#include	"LiDIA/matrix/dense_bigint_matrix_modules.h"
#include	"LiDIA/matrix/dense_bigint_matrix_modules.cc"
#include	"LiDIA/matrix/sparse_bigint_matrix_modules.h"
#include	"LiDIA/matrix/sparse_bigint_matrix_modules.cc"
#include	"LiDIA/matrix/bigint_matrix_algorithms.h"
#include	"LiDIA/matrix/bigint_matrix_algorithms.cc"
#include	"LiDIA/matrix/crt_and_prime_handling.h"



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif


template class COSMM< bigint >;
template class RODMM< bigint >;

template class hermite< bigint, RODMM< bigint >, RODMM< bigint > >;
template class Bradley< bigint, RODMM< bigint >, RODMM< bigint > >;
template class Iliopoulos< bigint, RODMM< bigint >, RODMM< bigint > >;
template class opt< bigint, RODMM< bigint >, RODMM< bigint > >;
template class Blankinship< bigint, RODMM< bigint >, RODMM< bigint > >;
template class Best_remainder< bigint, RODMM< bigint >, RODMM< bigint > >;
template class havas_best_remainder< bigint, RODMM< bigint >, RODMM< bigint > >;
template class havas_euclidean_norm< bigint, RODMM< bigint >, RODMM< bigint > >;
template class havas_minimal_norm< bigint, RODMM< bigint >, RODMM< bigint > >;
template class havas_sorting_gcd< bigint, RODMM< bigint >, RODMM< bigint > >;
template class havas_min_no_of_elements< bigint, RODMM< bigint >, RODMM< bigint > >;

template class havas_min_abs_of_row_plus_minimal_euclidean_norm< bigint, RODMM< bigint >, RODMM< bigint > >;
template class havas_min_abs_of_row_plus_minimal_norm< bigint, RODMM< bigint >, RODMM< bigint > >;
template class havas_min_abs_of_row_plus_min_no_of_elements< bigint, RODMM< bigint >, RODMM< bigint > >;

template class havas_minimal_norm_plus_min_abs_of_row< bigint, RODMM< bigint >, RODMM< bigint > >;
template class havas_minimal_norm_plus_sorting_gcd< bigint, RODMM< bigint >, RODMM< bigint > >;
template class havas_minimal_norm_plus_min_no_of_elements< bigint, RODMM< bigint >, RODMM< bigint > >;

template class havas_sorting_gcd_plus_minimal_euclidean_norm< bigint, RODMM< bigint >, RODMM< bigint > >;
template class havas_sorting_gcd_plus_minimal_norm< bigint, RODMM< bigint >, RODMM< bigint > >;
template class havas_sorting_gcd_plus_min_no_of_elements< bigint, RODMM< bigint >, RODMM< bigint > >;

template class havas_min_no_of_elements_plus_min_abs_of_row< bigint, RODMM< bigint >, RODMM< bigint > >;
template class havas_min_no_of_elements_plus_minimal_euclidean_norm< bigint, RODMM< bigint >, RODMM< bigint > >;
template class havas_min_no_of_elements_plus_minimal_norm< bigint, RODMM< bigint >, RODMM< bigint > >;
template class havas_min_no_of_elements_plus_sorting_gcd< bigint, RODMM< bigint >, RODMM< bigint > >;

template class Storjohann< bigint, RODMM< bigint >, RODMM< bigint > >;
template class Heuristik< bigint, RODMM< bigint >, RODMM< bigint > >;

template class havas_best_remainder_ext< bigint, RODMM< bigint >, RODMM< bigint > >;
template class havas_euclidean_norm_ext< bigint, RODMM< bigint >, RODMM< bigint > >;


template class hermite< bigint, COSMM< bigint >, COSMM< bigint > >;
template class Bradley< bigint, COSMM< bigint >, COSMM< bigint > >;
template class Iliopoulos< bigint, COSMM< bigint >, COSMM< bigint > >;
template class opt< bigint, COSMM< bigint >, COSMM< bigint > >;
template class Blankinship< bigint, COSMM< bigint >, COSMM< bigint > >;
template class Best_remainder< bigint, COSMM< bigint >, COSMM< bigint > >;
template class havas_best_remainder< bigint, COSMM< bigint >, COSMM< bigint > >;
template class havas_euclidean_norm< bigint, COSMM< bigint >, COSMM< bigint > >;
template class havas_minimal_norm< bigint, COSMM< bigint >, COSMM< bigint > >;
template class havas_sorting_gcd< bigint, COSMM< bigint >, COSMM< bigint > >;
template class havas_min_no_of_elements< bigint, COSMM< bigint >, COSMM< bigint > >;

template class havas_min_abs_of_row_plus_minimal_euclidean_norm< bigint, COSMM< bigint >, COSMM< bigint > >;
template class havas_min_abs_of_row_plus_minimal_norm< bigint, COSMM< bigint >, COSMM< bigint > >;
template class havas_min_abs_of_row_plus_min_no_of_elements< bigint, COSMM< bigint >, COSMM< bigint > >;

template class havas_minimal_norm_plus_min_abs_of_row< bigint, COSMM< bigint >, COSMM< bigint > >;
template class havas_minimal_norm_plus_sorting_gcd< bigint, COSMM< bigint >, COSMM< bigint > >;
template class havas_minimal_norm_plus_min_no_of_elements< bigint, COSMM< bigint >, COSMM< bigint > >;

template class havas_sorting_gcd_plus_minimal_euclidean_norm< bigint, COSMM< bigint >, COSMM< bigint > >;
template class havas_sorting_gcd_plus_minimal_norm< bigint, COSMM< bigint >, COSMM< bigint > >;
template class havas_sorting_gcd_plus_min_no_of_elements< bigint, COSMM< bigint >, COSMM< bigint > >;

template class havas_min_no_of_elements_plus_min_abs_of_row< bigint, COSMM< bigint >, COSMM< bigint > >;
template class havas_min_no_of_elements_plus_minimal_euclidean_norm< bigint, COSMM< bigint >, COSMM< bigint > >;
template class havas_min_no_of_elements_plus_minimal_norm< bigint, COSMM< bigint >, COSMM< bigint > >;
template class havas_min_no_of_elements_plus_sorting_gcd< bigint, COSMM< bigint >, COSMM< bigint > >;

template class Storjohann< bigint, COSMM< bigint >, COSMM< bigint > >;
template class Heuristik< bigint, COSMM< bigint >, COSMM< bigint > >;

template class havas_best_remainder_ext< long, COSMM< long >, COSMM< bigint > >;
template class havas_best_remainder_ext< bigint, COSMM< bigint >, COSMM< bigint > >;
template class havas_euclidean_norm_ext< long, COSMM< long >, COSMM< bigint > >;
template class nf_conf3e< long, COSMM< long >, COSMM< bigint > >;


//
// dense representation
//

template class havas_kernel< bigint, RODMM< bigint >, RODMM< bigint >,
	hermite< bigint, RODMM< bigint >, RODMM< bigint > > >;

template class havas_kernel< bigint, RODMM< bigint >, RODMM< bigint >,
	Bradley< bigint, RODMM< bigint >, RODMM< bigint > > >;

template class havas_kernel< bigint, RODMM< bigint >, RODMM< bigint >,
	Iliopoulos< bigint, RODMM< bigint >, RODMM< bigint > > >;

template class havas_kernel< bigint, RODMM< bigint >, RODMM< bigint >,
	opt< bigint, RODMM< bigint >, RODMM< bigint > > >;

template class havas_kernel< bigint, RODMM< bigint >, RODMM< bigint >,
	Blankinship< bigint, RODMM< bigint >, RODMM< bigint > > >;

template class havas_kernel< bigint, RODMM< bigint >, RODMM< bigint >,
	Best_remainder< bigint, RODMM< bigint >, RODMM< bigint > > >;

template class havas_kernel< bigint, RODMM< bigint >, RODMM< bigint >,
	havas_best_remainder< bigint, RODMM< bigint >, RODMM< bigint > > >;

template class havas_kernel< bigint, RODMM< bigint >, RODMM< bigint >,
	havas_euclidean_norm< bigint, RODMM< bigint >, RODMM< bigint > > >;

template class havas_kernel< bigint, RODMM< bigint >, RODMM< bigint >,
	havas_minimal_norm< bigint, RODMM< bigint >, RODMM< bigint > > >;

template class havas_kernel< bigint, RODMM< bigint >, RODMM< bigint >,
	havas_sorting_gcd< bigint, RODMM< bigint >, RODMM< bigint > > >;

template class havas_kernel< bigint, RODMM< bigint >, RODMM< bigint >,
	havas_min_no_of_elements< bigint, RODMM< bigint >, RODMM< bigint > > >;


template class havas_kernel< bigint, RODMM< bigint >, RODMM< bigint >,
	havas_min_abs_of_row_plus_minimal_euclidean_norm< bigint, RODMM< bigint >, RODMM< bigint > > >;

template class havas_kernel< bigint, RODMM< bigint >, RODMM< bigint >,
	havas_min_abs_of_row_plus_minimal_norm< bigint, RODMM< bigint >, RODMM< bigint > > >;

template class havas_kernel< bigint, RODMM< bigint >, RODMM< bigint >,
	havas_min_abs_of_row_plus_min_no_of_elements< bigint, RODMM< bigint >, RODMM< bigint > > >;


template class havas_kernel< bigint, RODMM< bigint >, RODMM< bigint >,
	havas_minimal_norm_plus_min_abs_of_row< bigint, RODMM< bigint >, RODMM< bigint > > >;

template class havas_kernel< bigint, RODMM< bigint >, RODMM< bigint >,
	havas_minimal_norm_plus_sorting_gcd< bigint, RODMM< bigint >, RODMM< bigint > > >;

template class havas_kernel< bigint, RODMM< bigint >, RODMM< bigint >,
	havas_minimal_norm_plus_min_no_of_elements< bigint, RODMM< bigint >, RODMM< bigint > > >;


template class havas_kernel< bigint, RODMM< bigint >, RODMM< bigint >,
	havas_sorting_gcd_plus_minimal_euclidean_norm< bigint, RODMM< bigint >, RODMM< bigint > > >;

template class havas_kernel< bigint, RODMM< bigint >, RODMM< bigint >,
	havas_sorting_gcd_plus_minimal_norm< bigint, RODMM< bigint >, RODMM< bigint > > >;

template class havas_kernel< bigint, RODMM< bigint >, RODMM< bigint >,
	havas_sorting_gcd_plus_min_no_of_elements< bigint, RODMM< bigint >, RODMM< bigint > > >;


template class havas_kernel< bigint, RODMM< bigint >, RODMM< bigint >,
	havas_min_no_of_elements_plus_min_abs_of_row< bigint, RODMM< bigint >, RODMM< bigint > > >;

template class havas_kernel< bigint, RODMM< bigint >, RODMM< bigint >,
	havas_min_no_of_elements_plus_minimal_euclidean_norm< bigint, RODMM< bigint >, RODMM< bigint > > >;

template class havas_kernel< bigint, RODMM< bigint >, RODMM< bigint >,
	havas_min_no_of_elements_plus_minimal_norm< bigint, RODMM< bigint >, RODMM< bigint > > >;

template class havas_kernel< bigint, RODMM< bigint >, RODMM< bigint >,
	havas_min_no_of_elements_plus_sorting_gcd< bigint, RODMM< bigint >, RODMM< bigint > > >;


template class havas_kernel< bigint, RODMM< bigint >, RODMM< bigint >,
	Storjohann< bigint, RODMM< bigint >, RODMM< bigint > > >;

template class havas_kernel< bigint, RODMM< bigint >, RODMM< bigint >,
	Heuristik< bigint, RODMM< bigint >, RODMM< bigint > > >;


template class havas_kernel< bigint, RODMM< bigint >, RODMM< bigint >,
	havas_best_remainder_ext< bigint, RODMM< bigint >, RODMM< bigint > > >;

template class havas_kernel< bigint, RODMM< bigint >, RODMM< bigint >,
	havas_euclidean_norm_ext< bigint, RODMM< bigint >, RODMM< bigint > > >;

template class havas_kernel< bigint, RODMM< bigint >, RODMM< bigint >,
	nf_conf3e< bigint, RODMM< bigint >, RODMM< bigint > > >;


template class kannan_kernel< bigint, RODMM< bigint >, RODMM< bigint >,
	Standard_normalization< bigint, RODMM< bigint >, RODMM< bigint > > >;

template class kannan_kernel< bigint, RODMM< bigint >, RODMM< bigint >,
	ChouCollins_normalization< bigint, RODMM< bigint >, RODMM< bigint > > >;

template class kannan_kernel< bigint, RODMM< bigint >, RODMM< bigint >,
	Jacobson_normalization< bigint, RODMM< bigint >, RODMM< bigint > > >;

template class normalization_kernel< bigint, RODMM< bigint >, RODMM< bigint > >;

template class Standard_normalization< bigint, RODMM< bigint >, RODMM< bigint > >;
template class ChouCollins_normalization< bigint, RODMM< bigint >, RODMM< bigint > >;
template class Jacobson_normalization< bigint, RODMM< bigint >, RODMM< bigint > >;

//
// sparse representation
//

template class havas_kernel< bigint, COSMM< bigint >, COSMM< bigint >,
	hermite< bigint, COSMM< bigint >, COSMM< bigint > > >;

template class havas_kernel< bigint, COSMM< bigint >, COSMM< bigint >,
	Bradley< bigint, COSMM< bigint >, COSMM< bigint > > >;

template class havas_kernel< bigint, COSMM< bigint >, COSMM< bigint >,
	Iliopoulos< bigint, COSMM< bigint >, COSMM< bigint > > >;

template class havas_kernel< bigint, COSMM< bigint >, COSMM< bigint >,
	opt< bigint, COSMM< bigint >, COSMM< bigint > > >;

template class havas_kernel< bigint, COSMM< bigint >, COSMM< bigint >,
	Blankinship< bigint, COSMM< bigint >, COSMM< bigint > > >;

template class havas_kernel< bigint, COSMM< bigint >, COSMM< bigint >,
	Best_remainder< bigint, COSMM< bigint >, COSMM< bigint > > >;

template class havas_kernel< bigint, COSMM< bigint >, COSMM< bigint >,
	havas_best_remainder< bigint, COSMM< bigint >, COSMM< bigint > > >;

template class havas_kernel< bigint, COSMM< bigint >, COSMM< bigint >,
	havas_euclidean_norm< bigint, COSMM< bigint >, COSMM< bigint > > >;

template class havas_kernel< bigint, COSMM< bigint >, COSMM< bigint >,
	havas_minimal_norm< bigint, COSMM< bigint >, COSMM< bigint > > >;

template class havas_kernel< bigint, COSMM< bigint >, COSMM< bigint >,
	havas_sorting_gcd< bigint, COSMM< bigint >, COSMM< bigint > > >;

template class havas_kernel< bigint, COSMM< bigint >, COSMM< bigint >,
	havas_min_no_of_elements< bigint, COSMM< bigint >, COSMM< bigint > > >;


template class havas_kernel< bigint, COSMM< bigint >, COSMM< bigint >,
	havas_min_abs_of_row_plus_minimal_euclidean_norm< bigint, COSMM< bigint >, COSMM< bigint > > >;

template class havas_kernel< bigint, COSMM< bigint >, COSMM< bigint >,
	havas_min_abs_of_row_plus_minimal_norm< bigint, COSMM< bigint >, COSMM< bigint > > >;

template class havas_kernel< bigint, COSMM< bigint >, COSMM< bigint >,
	havas_min_abs_of_row_plus_min_no_of_elements< bigint, COSMM< bigint >, COSMM< bigint > > >;


template class havas_kernel< bigint, COSMM< bigint >, COSMM< bigint >,
	havas_minimal_norm_plus_min_abs_of_row< bigint, COSMM< bigint >, COSMM< bigint > > >;

template class havas_kernel< bigint, COSMM< bigint >, COSMM< bigint >,
	havas_minimal_norm_plus_sorting_gcd< bigint, COSMM< bigint >, COSMM< bigint > > >;

template class havas_kernel< bigint, COSMM< bigint >, COSMM< bigint >,
	havas_minimal_norm_plus_min_no_of_elements< bigint, COSMM< bigint >, COSMM< bigint > > >;

template class havas_kernel< bigint, COSMM< bigint >, COSMM< bigint >,
	havas_sorting_gcd_plus_minimal_euclidean_norm< bigint, COSMM< bigint >, COSMM< bigint > > >;

template class havas_kernel< bigint, COSMM< bigint >, COSMM< bigint >,
	havas_sorting_gcd_plus_minimal_norm< bigint, COSMM< bigint >, COSMM< bigint > > >;

template class havas_kernel< bigint, COSMM< bigint >, COSMM< bigint >,
	havas_sorting_gcd_plus_min_no_of_elements< bigint, COSMM< bigint >, COSMM< bigint > > >;


template class havas_kernel< bigint, COSMM< bigint >, COSMM< bigint >,
	havas_min_no_of_elements_plus_min_abs_of_row< bigint, COSMM< bigint >, COSMM< bigint > > >;

template class havas_kernel< bigint, COSMM< bigint >, COSMM< bigint >,
	havas_min_no_of_elements_plus_minimal_euclidean_norm< bigint, COSMM< bigint >, COSMM< bigint > > >;

template class havas_kernel< bigint, COSMM< bigint >, COSMM< bigint >,
	havas_min_no_of_elements_plus_minimal_norm< bigint, COSMM< bigint >, COSMM< bigint > > >;

template class havas_kernel< bigint, COSMM< bigint >, COSMM< bigint >,
	havas_min_no_of_elements_plus_sorting_gcd< bigint, COSMM< bigint >, COSMM< bigint > > >;


template class havas_kernel< bigint, COSMM< bigint >, COSMM< bigint >,
	Storjohann< bigint, COSMM< bigint >, COSMM< bigint > > >;

template class havas_kernel< bigint, COSMM< bigint >, COSMM< bigint >,
	Heuristik< bigint, COSMM< bigint >, COSMM< bigint > > >;


template class havas_kernel< long, COSMM< long >, COSMM< bigint >,
	havas_best_remainder_ext< long, COSMM< long >, COSMM< bigint > > >;

template class havas_kernel< bigint, COSMM< bigint >, COSMM< bigint >,
	havas_best_remainder_ext< bigint, COSMM< bigint >, COSMM< bigint > > >;

template class havas_kernel< long, COSMM< long >, COSMM< bigint >,
	havas_euclidean_norm_ext< long, COSMM< long >, COSMM< bigint > > >;

template class havas_kernel< long, COSMM< long >, COSMM< bigint >,
	nf_conf3e< long, COSMM< long >, COSMM< bigint > > >;


template class kannan_kernel< bigint, COSMM< bigint >, COSMM< bigint >,
	Standard_normalization< bigint, COSMM< bigint >, COSMM< bigint > > >;

template class kannan_kernel< bigint, COSMM< bigint >, COSMM< bigint >,
	ChouCollins_normalization< bigint, COSMM< bigint >, COSMM< bigint > > >;

template class kannan_kernel< bigint, COSMM< bigint >, COSMM< bigint >,
	Jacobson_normalization< bigint, COSMM< bigint >, COSMM< bigint > > >;

template class normalization_kernel< bigint, COSMM< bigint >, COSMM< bigint > >;

template class Standard_normalization< bigint, COSMM< bigint >, COSMM< bigint > >;
template class ChouCollins_normalization< bigint, COSMM< bigint >, COSMM< bigint > >;
template class Jacobson_normalization< bigint, COSMM< bigint >, COSMM< bigint > >;

template class COSMM< long >;

//
// special for transformation matrix computation
//

template class havas_best_remainder_ext< long, COSMM< long >, RODMM< bigint > >;
template class havas_euclidean_norm_ext< long, COSMM< long >, RODMM< bigint > >;
template class nf_conf3e< long, COSMM< long >, RODMM< bigint > >;

template class havas_kernel< long, COSMM< long >, RODMM< bigint >,
	nf_conf3e< long, COSMM< long >, RODMM< bigint > > >;

template class havas_kernel< long, COSMM< long >, RODMM< bigint >,
	havas_best_remainder_ext< long, COSMM< long >, RODMM< bigint > > >;

template class havas_kernel< long, COSMM< long >, RODMM< bigint >,
	havas_euclidean_norm_ext< long, COSMM< long >, RODMM< bigint > > >;

template class normalization_kernel< long, COSMM< long >, RODMM< bigint > >;
template class normalization_kernel< long, COSMM< long >, COSMM< bigint > >;



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
