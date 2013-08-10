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
#include	"LiDIA/matrix/bigint_matrix_algorithms.h"
#include	"LiDIA/matrix/bigint_matrix_algorithms.cc"
#include	"LiDIA/matrix/modular_arithmetic.h"
#include	"LiDIA/matrix/modular_arithmetic.cc"
#include	"LiDIA/matrix/crt_and_prime_handling.h"
#include	"LiDIA/matrix/dense_fp_matrix_kernel.h"
#include	"LiDIA/matrix/dense_fp_matrix_kernel.cc"
#include	"LiDIA/matrix/sparse_fp_matrix_kernel.h"
#include	"LiDIA/matrix/sparse_fp_matrix_kernel.cc"
#include	"LiDIA/matrix/sparse_fp_matrix_algorithms.h"
#include	"LiDIA/matrix/sparse_fp_matrix_algorithms.cc"



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



//
// dense representation
//

template class modular_bigint_matrix_algorithms< DRMK< bigint >,
	dense_fp_matrix_kernel< long, MR< long > >,
	dense_fp_matrix_kernel< bigint, MR< bigint > > >;

template class modular_arithmetic< DRMK< bigint >,
	dense_fp_matrix_kernel< long, MR< long > >,
	dense_fp_matrix_kernel< bigint, MR< bigint > > >;

template class dense_fp_matrix_kernel< bigint, MR< bigint > >;
template class dense_fp_matrix_kernel< bigint, ring_matrix< bigint > >;

template class dense_fp_matrix_kernel< long, MR< long > >;
template class dense_fp_matrix_kernel< long, ring_matrix< long > >;

//
// sparse representation
//

template class modular_bigint_matrix_algorithms< SRMK< bigint >,
	sparse_fp_matrix_kernel< long, MR< long > >,
	sparse_fp_matrix_kernel< bigint, MR< bigint > > >;

template class modular_arithmetic< SRMK< bigint >,
	sparse_fp_matrix_kernel< long, MR< long > >,
	sparse_fp_matrix_kernel< bigint, MR< bigint > > >;

template class sparse_fp_matrix_kernel< bigint, MR< bigint > >;
template class sparse_fp_matrix_kernel< bigint, ring_matrix< bigint > >;

template class sparse_fp_matrix_kernel< long, MR< long > >;
template class sparse_fp_matrix_kernel< long, ring_matrix< long > >;

template class sparse_fp_matrix_algorithms< bigint, MR< bigint > >;
template class sparse_fp_matrix_algorithms< bigint, ring_matrix< bigint > >;

template class sparse_fp_matrix_algorithms< long, MR< long > >;
template class sparse_fp_matrix_algorithms< long, ring_matrix< long > >;



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
