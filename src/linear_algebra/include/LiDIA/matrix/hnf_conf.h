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


#ifndef LIDIA_HNF_CONF_H_GUARD_
#define LIDIA_HNF_CONF_H_GUARD_

#include	<fstream>
#ifndef LIDIA_HNF_KERNEL_H_GUARD_
# include	"LiDIA/matrix/hnf_kernel.h"
#endif



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif



////////////////////////////////////////////////////////////////////////////////////////////
#define schema_hnf(name) template< class T, class REP, class REP1 > \
class name                                                                                 \
{	                                                                                   \
                                                                                           \
private:                                                                                   \
                                                                                           \
  const REP init_modul; \
  const REP pivot_modul; \
  const REP elim_modul; \
                                                                                           \
  REP modul; \
                                                                                           \
public:                                                                                    \
                                                                                           \
  name (); \
  ~name (); \
                                                                                           \
  /* init */                                                                               \
  void init(MR< T > &, const T &); \
                                                                                           \
  /* mgcd */                                                                               \
  lidia_size_t mgcd(MR< T > &, lidia_size_t, lidia_size_t) const; \
  lidia_size_t mgcd(MR< T > &, matrix< bigint > &, lidia_size_t, lidia_size_t) const; \
                                                                                           \
  /* normalization */                                                                      \
  bool normalize_row(MR< T > &, lidia_size_t, lidia_size_t, lidia_size_t, \
                     lidia_size_t) const; \
                                                                                           \
  bool normalize_row(MR< T > &, matrix< bigint > &, lidia_size_t, lidia_size_t, \
                     lidia_size_t index, lidia_size_t len) const; \
};

////////////////////////////////////////////////////////////////////////////////////////////
#define schema_hnf_extended(name) template< class T, class REP, class REP1 > \
class name                                                                                 \
{	                                                                                   \
                                                                                           \
private:                                                                                   \
                                                                                           \
  const REP init_modul; \
  const REP pivot_modul; \
  const REP elim_modul; \
                                                                                           \
  REP modul; \
                                                                                           \
  T *max_array; \
  T BOUND; \
                                                                                           \
public:                                                                                    \
                                                                                           \
  name (); \
  ~name (); \
                                                                                           \
  /* init */                                                                               \
  void init(MR< T > &, const T &); \
                                                                                           \
  /* mgcd */                                                                               \
  lidia_size_t mgcd(MR< T > &, lidia_size_t, lidia_size_t) const; \
  lidia_size_t mgcd(MR< T > &, matrix< bigint > &, lidia_size_t, lidia_size_t) const; \
                                                                                           \
  /* normalization */                                                                      \
  bool normalize_row(MR< T > &, lidia_size_t, lidia_size_t, lidia_size_t, \
                     lidia_size_t) const; \
                                                                                           \
  bool normalize_row(MR< T > &, matrix< bigint > &, lidia_size_t, lidia_size_t, \
                     lidia_size_t index, lidia_size_t len) const; \
};


//
// Hermite normal form configurations
//

schema_hnf(hermite)
	schema_hnf(Bradley)
	schema_hnf(Iliopoulos)
	schema_hnf(opt)
	schema_hnf(Blankinship)
	schema_hnf(Best_remainder)
	schema_hnf(havas_best_remainder)
	schema_hnf(havas_euclidean_norm)
	schema_hnf(havas_minimal_norm)
	schema_hnf(havas_sorting_gcd)
	schema_hnf(havas_min_no_of_elements)

	schema_hnf(havas_min_abs_of_row_plus_minimal_euclidean_norm)
	schema_hnf(havas_min_abs_of_row_plus_minimal_norm)
	schema_hnf(havas_min_abs_of_row_plus_min_no_of_elements)

	schema_hnf(havas_minimal_norm_plus_min_abs_of_row)
	schema_hnf(havas_minimal_norm_plus_sorting_gcd)
	schema_hnf(havas_minimal_norm_plus_min_no_of_elements)

	schema_hnf(havas_sorting_gcd_plus_minimal_euclidean_norm)
	schema_hnf(havas_sorting_gcd_plus_minimal_norm)
	schema_hnf(havas_sorting_gcd_plus_min_no_of_elements)

	schema_hnf(havas_min_no_of_elements_plus_min_abs_of_row)
	schema_hnf(havas_min_no_of_elements_plus_minimal_euclidean_norm)
	schema_hnf(havas_min_no_of_elements_plus_minimal_norm)
	schema_hnf(havas_min_no_of_elements_plus_sorting_gcd)

	schema_hnf(Storjohann)
	schema_hnf(Heuristik)

	schema_hnf_extended(havas_best_remainder_ext)
	schema_hnf_extended(havas_euclidean_norm_ext)
	schema_hnf_extended(nf_conf3e)

/////////////////////////////////////////////////////////////////////////////////////////

#define schema_hnf_kannan(name) template< class T, class REP, class REP1 > \
class name                                                                              \
{	                                                                                \
private:                                                                                \
                                                                                        \
  const normalization_kernel< T, REP, REP1 > modul; \
                                                                                        \
public:                                                                                 \
                                                                                        \
  /* constructor */                                                                     \
                                                                                        \
  name (); \
                                                                                        \
  /* destructor */                                                                      \
                                                                                        \
  ~name (); \
                                                                                        \
  /* normalize */                                                                       \
                                                                                        \
  lidia_size_t *normalize(MR< T > &, lidia_size_t, lidia_size_t) const; \
  lidia_size_t *normalize(MR< T > &, matrix< bigint > &, lidia_size_t, \
                          lidia_size_t) const; \
};

//
// kannan_kernel: Hermite normal form configurations
//
	schema_hnf_kannan(Standard_normalization)
		schema_hnf_kannan(ChouCollins_normalization)
		schema_hnf_kannan(Jacobson_normalization)



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif



#endif	// LIDIA_HNF_CONF_H_GUARD_
