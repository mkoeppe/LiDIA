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


#ifndef LIDIA_NORMALIZE_KERNEL_H_GUARD_
#define LIDIA_NORMALIZE_KERNEL_H_GUARD_



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif



template< class T, class REP, class REP1 >
class normalization_kernel
{
private:

	const REP modul;
	const REP1 modul_TR;

public:

	//
	// constructor
	//

	normalization_kernel() {}

	//
	// destructor
	//

	~normalization_kernel() {}

	//
	// normalize
	//

	lidia_size_t *normalize_Std(MR< T > &,
				    lidia_size_t,
				    lidia_size_t) const;
	lidia_size_t *normalize_Std(MR< T > &,
				    matrix< bigint > &,
				    lidia_size_t,
				    lidia_size_t) const;

	lidia_size_t *normalize_ChouCollins(MR< T > &,
					    lidia_size_t,
					    lidia_size_t) const;
	lidia_size_t *normalize_ChouCollins(MR< T > &,
					    matrix< bigint > &TR,
					    lidia_size_t,
					    lidia_size_t) const;
	lidia_size_t *normalize_ChouCollins(MR< T > &,
					    trans_matrix &TR,
					    matrix< bigint > &tran,
					    lidia_size_t,
					    lidia_size_t) const;

	lidia_size_t *normalize_ChouCollins_extended(MR< T > &,
						     lidia_size_t,
						     lidia_size_t) const;

	lidia_size_t *normalizeMod_Std(MR< T > &,
				       lidia_size_t,
				       lidia_size_t) const;
	lidia_size_t *normalizeMod_Std(MR< T > &,
				       matrix< bigint > &,
				       lidia_size_t,
				       lidia_size_t) const;
	lidia_size_t *normalizeMod_Std(MR< T > &,
				       trans_matrix &,
				       matrix< bigint > &,
				       lidia_size_t s);

	lidia_size_t *normalizeMod_ChouCollins(MR< T > &,
					       lidia_size_t,
					       lidia_size_t) const;
	lidia_size_t *normalizeMod_ChouCollins(MR< T > &,
					       matrix< bigint > &,
					       lidia_size_t,
					       lidia_size_t) const;
	lidia_size_t *normalizeMod_ChouCollins(MR< T > &,
					       trans_matrix &,
					       matrix< bigint > &,
					       lidia_size_t s);

	lidia_size_t *normalizeHybrid_Std(MR< T > &,
					  lidia_size_t,
					  lidia_size_t) const;
	lidia_size_t *normalizeHybrid_Std(MR< T > &,
					  matrix< bigint > &,
					  lidia_size_t,
					  lidia_size_t) const;
	lidia_size_t *normalizeHybrid_Std(MR< T > &,
					  trans_matrix &,
					  matrix< bigint > &,
					  lidia_size_t s);

	lidia_size_t *normalizeHybrid_ChouCollins(MR< T > &,
						  lidia_size_t,
						  lidia_size_t) const;
	lidia_size_t *normalizeHybrid_ChouCollins(MR< T > &,
						  matrix< bigint > &,
						  lidia_size_t,
						  lidia_size_t) const;
	lidia_size_t *normalizeHybrid_ChouCollins(MR< T > &,
						  trans_matrix &,
						  matrix< bigint > &,
						  lidia_size_t s);
};



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif



//#ifdef LIDIA_INCLUDE_CC_GUARD_
//# include	"LiDIA/normalize_kernel.cc"
//#endif



#endif	// LIDIA_NORMALIZE_KERNEL_H_GUARD_
