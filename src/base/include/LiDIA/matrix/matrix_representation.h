// -*- C++ -*-
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


#ifndef LIDIA_MATRIX_REPRESENTATION_H_GUARD_
#define LIDIA_MATRIX_REPRESENTATION_H_GUARD_



#ifndef LIDIA_BIGINT_H_GUARD_
# include	"LiDIA/bigint.h"
#endif
#ifndef LIDIA_BIGMOD_H_GUARD_
# include	"LiDIA/bigmod.h"
#endif
#ifndef LIDIA_ARRAY_FUNCTIONS_H_GUARD_
# include	"LiDIA/base/array_functions.h"
#endif
#ifndef LIDIA_MATRIX_FLAGS_H_GUARD_
# include	"LiDIA/matrix_flags.h"
#endif



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif



// type defines
#define matrix_representation MR

// base defines
#define DBMKex DBMK< T >
#define SBMKex SBMK< T >

// ring defines
#define DRMKex DRMK< T >
#define SRMKex SRMK< T >

// field defines
#define DFMKex DFMK< T >
#define SFMKex SFMK< T >


#define FRIEND_HNF(A) friend class A< T, RODMM< T >, RODMM< bigint > >; \
                      friend class A< T, CODMM< T >, CODMM< bigint > >; \
                      friend class A< T, ROSMM< T >, ROSMM< bigint > >; \
                      friend class A< T, COSMM< T >, COSMM< bigint > >; \
                      friend class havas_kernel< T, RODMM< T >, RODMM< bigint >, A< T, RODMM< T >, RODMM< bigint > > >; \
                      friend class havas_kernel< T, COSMM< T >, COSMM< bigint >, A< T, COSMM< T >, COSMM< bigint > > >

// base matrix
template< class T > class DBMK;
template< class T > class SBMK;

template< class T, class R1, class R2 > class BMA;

template< class T > class base_matrix;

// ring matrix
template< class T > class DRMK;
template< class T > class SRMK;

template< class T, class R1, class R2, class R3 > class RMA;

template< class T > class ring_matrix;

// field matrix
template< class T > class DFMK;
template< class T > class SFMK;

template< class T, class R1, class R2, class R3 > class FMA;

// user interface
template< class T > class matrix;
template< class T > class dense_matrix;


// matrix< bigint >
template< class T > class RODMM;
template< class T > class CODMM;
template< class T > class ROSMM;
template< class T > class COSMM;

template< class T, class R, class R1 > class hermite;
template< class T, class R, class R1 > class Bradley;
template< class T, class R, class R1 > class Iliopoulos;
template< class T, class R, class R1 > class opt;
template< class T, class R, class R1 > class Blankinship;
template< class T, class R, class R1 > class Best_remainder;
template< class T, class R, class R1 > class havas_best_remainder;
template< class T, class R, class R1 > class havas_euclidean_norm;
template< class T, class R, class R1 > class havas_minimal_norm;
template< class T, class R, class R1 > class havas_sorting_gcd;
template< class T, class R, class R1 > class havas_min_no_of_elements;

template< class T, class R, class R1 > class havas_min_abs_of_row_plus_minimal_euclidean_norm;
template< class T, class R, class R1 > class havas_min_abs_of_row_plus_minimal_norm;
template< class T, class R, class R1 > class havas_min_abs_of_row_plus_min_no_of_elements;

template< class T, class R, class R1 > class havas_minimal_norm_plus_min_abs_of_row;
template< class T, class R, class R1 > class havas_minimal_norm_plus_sorting_gcd;
template< class T, class R, class R1 > class havas_minimal_norm_plus_min_no_of_elements;

template< class T, class R, class R1 > class havas_sorting_gcd_plus_minimal_euclidean_norm;
template< class T, class R, class R1 > class havas_sorting_gcd_plus_minimal_norm;
template< class T, class R, class R1 > class havas_sorting_gcd_plus_min_no_of_elements;

template< class T, class R, class R1 > class havas_min_no_of_elements_plus_min_abs_of_row;
template< class T, class R, class R1 > class havas_min_no_of_elements_plus_minimal_euclidean_norm;
template< class T, class R, class R1 > class havas_min_no_of_elements_plus_minimal_norm;
template< class T, class R, class R1 > class havas_min_no_of_elements_plus_sorting_gcd;

template< class T, class R, class R1 > class Storjohann;
template< class T, class R, class R1 > class Heuristik;

template< class T, class R, class R1 > class havas_best_remainder_ext;
template< class T, class R, class R1 > class havas_euclidean_norm_ext;
template< class T, class R, class R1 > class nf_conf3e;

template< class T, class R, class R1 > class Standard_normalization;
template< class T, class R, class R1 > class ChouCollins_normalization;
template< class T, class R, class R1 > class Jacobson_normalization;

template< class T, class R, class R1, class C > class havas_kernel;
template< class T, class R, class R1, class C > class kannan_kernel;
template< class T, class R, class R1 > class normalization_kernel;

template< class T1, class T2, class T3 > class bigint_matrix_algorithms;
template< class T, class S, class M > class modular_arithmetic;
template< class T, class S, class M > class modular_bigint_matrix_algorithms;

// matrix< fp >
template< class T, class M > class dense_fp_matrix_kernel;
template< class T, class M > class sparse_fp_matrix_kernel;

template< class T, class M > class sparse_fp_matrix_algorithms;

template< class T , class M > class solver_kernel;


// data structure
template< class T >
struct matrix_representation
{
	lidia_size_t rows; // number of rows
	lidia_size_t columns; // number of columns

	lidia_size_t sparse_rows; // index of the first row in sparse rep.
	lidia_size_t sparse_columns; // index of the first column in sparse rep.

	T **value;

	lidia_size_t **index;

	lidia_size_t *value_counter;
	lidia_size_t *allocated;

	T Zero;

 	matrix_flags bitfield;
};



#undef DBMKex
#undef SBMKex

#undef DRMKex
#undef SRMKex

#undef DFMKex
#undef SFMKex

#undef matrix_representation



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif



#endif	// LIDIA_MATRIX_REPRESENTATION_H_GUARD_
