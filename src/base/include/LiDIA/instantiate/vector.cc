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
//	define either BASE_VECTOR, SORT_VECTOR, MATH_VECTOR, FILE_VECTOR or LIDIA_VECTOR
//
//	include this file
//


// =============== base_vector ==================

#ifdef BASE_VECTOR

# include	"LiDIA/base_vector.h"
# include	"LiDIA/base_vector.cc"



# ifdef LIDIA_NAMESPACE
namespace LiDIA {
# endif



template class vector_representation< TYPE >;
template class base_vector< TYPE >;



# ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# endif



#endif	// BASE_VECTOR



// =============== file_vector ==================

#ifdef FILE_VECTOR

# include	"LiDIA/file_vector.h"
# include	"LiDIA/file_vector.cc"



# ifdef LIDIA_NAMESPACE
namespace LiDIA {
# endif



template class file_vector< TYPE >;



# ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# endif



#endif	// FILE_VECTOR



// =============== math_vector ==================

#ifdef MATH_VECTOR

# include	"LiDIA/math_vector.h"
# include	"LiDIA/math_vector.cc"



# ifdef LIDIA_NAMESPACE
namespace LiDIA {
# endif



template class math_vector< TYPE >;



# ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# endif



#endif	// MATH_VECTOR



// =============== sort_vector ==================

#ifdef SORT_VECTOR

# include	"LiDIA/comparator.h"
# include	"LiDIA/sort_vector.h"
# include	"LiDIA/sort_vector.cc"



# ifdef LIDIA_NAMESPACE
namespace LiDIA {
# endif



template class comparator< TYPE >;
template class sort_vector< TYPE >;



# ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# endif



#endif	// SORT_VECTOR



// =============== lidia_vector ==================

#ifdef LIDIA_VECTOR

# include	"LiDIA/lidia_vector.h"



# ifdef LIDIA_NAMESPACE
namespace LiDIA {
# endif



template class lidia_vector< TYPE >;



# ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# endif



#endif	// LIDIA_VECTOR
