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
//	Author	:
//	Changes	: See CVS log
//
//==============================================================================================


//
//  usage:
//
//	include all headers that are necessary to describe the type TYPE
//
//	define the type TYPE
//
//	include this file
//


#if defined (BASE_DENSE_POWSER) || defined (DENSE_POWSER)
# include	"LiDIA/finite_fields/base_dense_power_series.h"
# include	"LiDIA/finite_fields/base_dense_power_series.cc"
#endif

#if defined (DENSE_POWSER)
# include	"LiDIA/dense_power_series.h"
# include	"LiDIA/dense_power_series.cc"
#endif

#if defined (BASE_SPARSE_POWSER) || defined (SPARSE_POWSER)
# include	"LiDIA/finite_fields/base_sparse_power_series.h"
# include	"LiDIA/finite_fields/base_sparse_power_series.cc"
#endif

#if defined (COEFF_SPARSE_POWSER) || defined (SPARSE_POWSER)
# include	"LiDIA/finite_fields/coeff_sparse_power_series.h"
# include	"LiDIA/finite_fields/coeff_sparse_power_series.cc"
#endif

#if defined (SPARSE_POWSER)
# include	"LiDIA/sparse_power_series.h"
# include	"LiDIA/sparse_power_series.cc"
#endif



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



// ================== dense_power_series ==================

#if defined (BASE_DENSE_POWSER) || defined (DENSE_POWSER)
template class base_dense_power_series< TYPE >;
#endif

#if defined (DENSE_POWSER)
template class      dense_power_series< TYPE >;
#endif


// ================== sparse_power_series ==================

#if defined (BASE_SPARSE_POWSER) || defined (SPARSE_POWSER)
template class base_sparse_power_series< TYPE >;
#endif

#if defined (COEFF_SPARSE_POWSER) || defined (SPARSE_POWSER)
template class spc< TYPE >;
#endif

#if defined (SPARSE_POWSER)
template class      sparse_power_series< TYPE >;
#endif



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
