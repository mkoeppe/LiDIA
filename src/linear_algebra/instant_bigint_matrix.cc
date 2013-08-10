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
#include	"LiDIA/matrix/crt_and_prime_handling.h"
#include	"LiDIA/matrix/dense_bigint_matrix_kernel.cc"
#include	"LiDIA/matrix/sparse_bigint_matrix_kernel.cc"
#include	"LiDIA/base_power_product.h"
#include	"LiDIA/base_power_product.cc"
#include	"LiDIA/base_ppair.h"
#include	"LiDIA/base_ppair.cc" // MM



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



template class dense_bigint_matrix_kernel< matrix< bigint > >;
template class sparse_bigint_matrix_kernel< matrix< bigint > >;

template class bigint_matrix_algorithms< SRMK< bigint >, SRMK< bigint >, SRMK< bigint > >;
template class bigint_matrix_algorithms< SRMK< bigint >, SRMK< bigint >, DRMK< bigint > >;
template class bigint_matrix_algorithms< SRMK< bigint >, DRMK< bigint >, SRMK< bigint > >;
template class bigint_matrix_algorithms< SRMK< bigint >, DRMK< bigint >, DRMK< bigint > >;

template class bigint_matrix_algorithms< DRMK< bigint >, SRMK< bigint >, SRMK< bigint > >;
template class bigint_matrix_algorithms< DRMK< bigint >, SRMK< bigint >, DRMK< bigint > >;
template class bigint_matrix_algorithms< DRMK< bigint >, DRMK< bigint >, SRMK< bigint > >;
template class bigint_matrix_algorithms< DRMK< bigint >, DRMK< bigint >, DRMK< bigint > >;

//
// base_power_product
//

template class base_power_product< ring_matrix< bigint >, int >;
template void swap (base_power_product< ring_matrix< bigint >, int > &,
                    base_power_product< ring_matrix< bigint >, int > &);
template std::istream & operator >> (std::istream &, base_power_product< ring_matrix< bigint >, int > &);
template std::ostream & operator << (std::ostream &, const base_power_product< ring_matrix< bigint >, int > &);

template class base_power_product< ring_matrix< long >, int >;
template void swap (base_power_product< ring_matrix< long >, int > &,
                    base_power_product< ring_matrix< long >, int > &);
template std::istream & operator >> (std::istream &, base_power_product< ring_matrix< long >, int > &);
template std::ostream & operator << (std::ostream &, const base_power_product< ring_matrix< long >, int > &);

//
// base_pair
//

template class base_ppair< ring_matrix< bigint >, int >;
template void swap (base_ppair< ring_matrix< bigint >, int > &,
                    base_ppair< ring_matrix< bigint >, int > &);
template std::istream & operator >> (std::istream &, base_ppair< ring_matrix< bigint >, int > &);
template std::ostream & operator << (std::ostream &, const base_ppair< ring_matrix< bigint >, int > &);

template class base_ppair< ring_matrix< long >, int >;
template void swap (base_ppair< ring_matrix< long >, int > &, base_ppair< ring_matrix< long >, int > &);
template std::istream & operator >> (std::istream &, base_ppair< ring_matrix< long >, int > &);
template std::ostream & operator << (std::ostream &, const base_ppair< ring_matrix< long >, int > &);



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
