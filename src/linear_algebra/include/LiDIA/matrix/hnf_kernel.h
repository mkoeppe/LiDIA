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


#ifndef LIDIA_HNF_KERNEL_H_GUARD_
#define LIDIA_HNF_KERNEL_H_GUARD_


#ifndef LIDIA_MATH_VECTOR_H_GUARD_
# include	"LiDIA/math_vector.h"
#endif
#ifndef LIDIA_NORMALIZE_KERNEL_H_GUARD_
# include	"LiDIA/matrix/normalize_kernel.h"
#endif
#ifndef LIDIA_BIGFLOAT_H_GUARD_
# include	"LiDIA/bigfloat.h"
#endif



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif



template< class T, class REP, class REP1, class CONF >
class havas_kernel
{
private:

	CONF conf;
	const REP modul;
	const REP1 modul_TR;

public:

	//
	// constructor
	//

	havas_kernel() {}

	//
	// destructor
	//

	~havas_kernel() {}

	//
	// column step form
	//

	bool hnf_Z2_mod(MR< T > &, lidia_size_t &, lidia_size_t &, const T &mod);

	// ADDED BY MJJ
	int hnf_Z2(MR< T > &, trans_matrix &, matrix< bigint > &,
		   lidia_size_t &, lidia_size_t &, lidia_size_t, const T &B = 0);
	//

	//
	// hermite normal form computation and column step form
	//

	//
	// row elimination
	//

	bool mgcd(MR< T > &, lidia_size_t &, lidia_size_t &, const T &B = 0);
	bool mgcd(MR< T > &, matrix< bigint > &, lidia_size_t &, lidia_size_t &, const T &B = 0);

	bool hnf_Z1(MR< T > &, lidia_size_t &, lidia_size_t &, const T &B = 0);
	bool hnf_Z1(MR< T > &, matrix< bigint > &, lidia_size_t &, lidia_size_t &, const T &B = 0);

	bool hnf_Z2(MR< T > &, lidia_size_t &, lidia_size_t &, const T &B = 0);
	bool hnf_Z2(MR< T > &, matrix< bigint > &, lidia_size_t &, lidia_size_t &, const T &B = 0);

	bool stf(MR< T > &, lidia_size_t &, lidia_size_t &, const T &B = 0);
	bool stf(MR< T > &, matrix< bigint > &, lidia_size_t &, lidia_size_t &, const T &B = 0);
};



template< class T, class REP, class REP1, class CONF >
class kannan_kernel
{
private:

	const CONF conf_modul;
	const REP modul;
	const REP1 modul_TR;

public:

	//
	// constructor
	//

	kannan_kernel() {}

	//
	// destructor
	//

	~kannan_kernel() {}

	//
	// column step form
	//

	bool stf(MR< T > &, lidia_size_t &, lidia_size_t &, const T &B = 0);
	bool stf(MR< T > &, matrix< bigint > &, lidia_size_t &, lidia_size_t &,
		 const T &B = 0);

	// ADDED BY MJJ
	int stf(MR< T > &, trans_matrix &, matrix< bigint > &,
		lidia_size_t &, lidia_size_t &, lidia_size_t, const T &B = 0);
	//


	//
	// hermite normal form
	//

	bool hnf(MR< T > &, lidia_size_t &, lidia_size_t &, const T &B = 0);
	bool hnf(MR< T > &, matrix< bigint > &, lidia_size_t &, lidia_size_t &,
		 const T &B = 0);

	lidia_size_t *normalize(MR< T > &, lidia_size_t);
	lidia_size_t *normalize(MR< T > &, matrix< bigint > &, lidia_size_t);

	// ADDED BY MJJ
	lidia_size_t *normalize_cg(MR< T > &, lidia_size_t);
	lidia_size_t *normalize_cg(MR< T > &, matrix< bigint > &TR, lidia_size_t s);
	lidia_size_t *normalize_cg(MR< T > &, trans_matrix &TR,
				   matrix< bigint > &tran, lidia_size_t s);
	//

	lidia_size_t *normalize_cg_extended(MR< T > &, lidia_size_t);
};



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif



//#ifdef LIDIA_INCLUDE_CC
//# include	"LiDIA/hnf_kernel.cc"
//#endif



#endif	// LIDIA_HNF_KERNEL_H_GUARD_
