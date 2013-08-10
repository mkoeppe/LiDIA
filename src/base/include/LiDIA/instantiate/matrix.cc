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


//
//  usage:
//
//	include all headers that are necessary to describe the type TYPE
//
//	define the type TYPE
//
//	define either NORMAL, DENSE, SPARSE, or SUBCLASSES
//
//	define either BASE_MATRIX, RING_MATRIX, or FIELD_MATRIX
//
//	include this file
//


#ifdef BASE_MATRIX

# ifdef SUBCLASSES
#  include	"LiDIA/base_matrix.h"
#  include	"LiDIA/base_matrix.cc"
#  include	"LiDIA/matrix/base_matrix_algorithms.h"
#  include	"LiDIA/matrix/base_matrix_algorithms.cc"
#  include	"LiDIA/matrix/dense_base_matrix_kernel.h"
#  include	"LiDIA/matrix/dense_base_matrix_kernel.cc"
#  include	"LiDIA/matrix/sparse_base_matrix_kernel.h"
#  include	"LiDIA/matrix/sparse_base_matrix_kernel.cc"

#  define DBMKex DBMK< TYPE >
#  define SBMKex SBMK< TYPE >



#  ifdef LIDIA_NAMESPACE
namespace LiDIA {
#  endif



template class DBMKex;
template class SBMKex;

template class BMA< TYPE, DBMKex, DBMKex >;
template class BMA< TYPE, DBMKex, SBMKex >;
template class BMA< TYPE, SBMKex, DBMKex >;
template class BMA< TYPE, SBMKex, SBMKex >;



#  ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#  endif

#  undef DBMKex
#  undef SBMKex

# endif	// SUBCLASSES



# ifdef NORMAL
#  include	"LiDIA/base_matrix.h"
#  include	"LiDIA/base_matrix.cc"
#  include	"LiDIA/matrix/base_matrix_algorithms.h"
#  include	"LiDIA/matrix/base_matrix_algorithms.cc"
#  include	"LiDIA/matrix/dense_base_matrix_kernel.h"
#  include	"LiDIA/matrix/dense_base_matrix_kernel.cc"
#  include	"LiDIA/matrix/sparse_base_matrix_kernel.h"
#  include	"LiDIA/matrix/sparse_base_matrix_kernel.cc"

#  ifdef LIDIA_NAMESPACE
namespace LiDIA {
#  endif



template class base_matrix< TYPE >;

template std::ostream & operator << (std::ostream &, const base_matrix< TYPE > &);
template std::istream & operator >> (std::istream &, base_matrix< TYPE > &);
template void swap(base_matrix< TYPE > &, base_matrix< TYPE > &);
template void assign(base_matrix< TYPE > &, const base_matrix< TYPE > &);
template void diag(base_matrix< TYPE > &, const TYPE &, const TYPE &);
template void trans(base_matrix< TYPE > &, const base_matrix< TYPE > &);
template base_matrix< TYPE > trans(const base_matrix< TYPE > &);



#  ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#  endif

# endif	// NORMAL



# ifdef DENSE
#  include	"LiDIA/dense_base_matrix.h"
#  include	"LiDIA/dense_base_matrix.cc"
#  include	"LiDIA/matrix/base_matrix_algorithms.h"
#  include	"LiDIA/matrix/base_matrix_algorithms.cc"
#  include	"LiDIA/matrix/dense_base_matrix_kernel.h"
#  include	"LiDIA/matrix/dense_base_matrix_kernel.cc"

#  ifdef LIDIA_NAMESPACE
namespace LiDIA {
#  endif



template class dense_base_matrix< TYPE >;

template std::ostream & operator << (std::ostream &, const dense_base_matrix< TYPE > &);
template std::istream & operator >> (std::istream &, dense_base_matrix< TYPE > &);
template void swap(dense_base_matrix< TYPE > &, dense_base_matrix< TYPE > &);
template void assign(dense_base_matrix< TYPE > &, const dense_base_matrix< TYPE > &);
template void diag(dense_base_matrix< TYPE > &, const TYPE &, const TYPE &);
template void trans(dense_base_matrix< TYPE > &, const dense_base_matrix< TYPE > &);
template dense_base_matrix< TYPE > trans(const dense_base_matrix< TYPE > &);



#  ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#  endif

# endif	// DENSE



# ifdef SPARSE
#  include	"LiDIA/sparse_base_matrix.h"
#  include	"LiDIA/sparse_base_matrix.cc"
#  include	"LiDIA/matrix/base_matrix_algorithms.h"
#  include	"LiDIA/matrix/base_matrix_algorithms.cc"
#  include	"LiDIA/matrix/sparse_base_matrix_kernel.h"
#  include	"LiDIA/matrix/sparse_base_matrix_kernel.cc"

#  ifdef LIDIA_NAMESPACE
namespace LiDIA {
#  endif



template class sparse_base_matrix< TYPE >;

template std::ostream & operator << (std::ostream &, const sparse_base_matrix< TYPE > &);
template std::istream & operator >> (std::istream &, sparse_base_matrix< TYPE > &);
template void swap(sparse_base_matrix< TYPE > &, sparse_base_matrix< TYPE > &);
template void assign(sparse_base_matrix< TYPE > &, const sparse_base_matrix< TYPE > &);
template void diag(sparse_base_matrix< TYPE > &, const TYPE &, const TYPE &);
template void trans(sparse_base_matrix< TYPE > &, const sparse_base_matrix< TYPE > &);
template sparse_base_matrix< TYPE > trans(const sparse_base_matrix< TYPE > &);



#  ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#  endif

# endif	// SPARSE

#endif	// BASE_MATRIX



#ifdef RING_MATRIX

# ifdef SUBCLASSES
#  include	"LiDIA/ring_matrix.h"
#  include	"LiDIA/ring_matrix.cc"
#  include	"LiDIA/matrix/ring_matrix_algorithms.h"
#  include	"LiDIA/matrix/ring_matrix_algorithms.cc"
#  include	"LiDIA/matrix/dense_ring_matrix_kernel.h"
#  include	"LiDIA/matrix/dense_ring_matrix_kernel.cc"
#  include	"LiDIA/matrix/sparse_ring_matrix_kernel.h"
#  include	"LiDIA/matrix/sparse_ring_matrix_kernel.cc"

#  define DRMKex DRMK< TYPE >
#  define SRMKex SRMK< TYPE >

#  ifdef LIDIA_NAMESPACE
namespace LiDIA {
#  endif



template class DRMKex;
template class SRMKex;

template class RMA< TYPE, DRMKex, DRMKex, DRMKex >;
template class RMA< TYPE, DRMKex, DRMKex, SRMKex >;
template class RMA< TYPE, DRMKex, SRMKex, DRMKex >;
template class RMA< TYPE, SRMKex, DRMKex, DRMKex >;
template class RMA< TYPE, DRMKex, SRMKex, SRMKex >;
template class RMA< TYPE, SRMKex, DRMKex, SRMKex >;
template class RMA< TYPE, SRMKex, SRMKex, DRMKex >;
template class RMA< TYPE, SRMKex, SRMKex, SRMKex >;



#  ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#  endif

#  undef DRMKex
#  undef SRMKex

# endif	// SUBCLASSES



# ifdef NORMAL
#  include	"LiDIA/ring_matrix.h"
#  include	"LiDIA/ring_matrix.cc"
#  include	"LiDIA/matrix/ring_matrix_algorithms.h"
#  include	"LiDIA/matrix/ring_matrix_algorithms.cc"
#  include	"LiDIA/matrix/dense_ring_matrix_kernel.h"
#  include	"LiDIA/matrix/dense_ring_matrix_kernel.cc"
#  include	"LiDIA/matrix/sparse_ring_matrix_kernel.h"
#  include	"LiDIA/matrix/sparse_ring_matrix_kernel.cc"

#  ifdef LIDIA_NAMESPACE
namespace LiDIA {
#  endif



template class ring_matrix< TYPE >;

template void add(ring_matrix< TYPE > &, const ring_matrix< TYPE > &, const ring_matrix< TYPE > &);
template void add(ring_matrix< TYPE > &, const ring_matrix< TYPE > &, const TYPE &);
template void add(ring_matrix< TYPE > &, const TYPE &, const ring_matrix< TYPE > &);
template void subtract(ring_matrix< TYPE > &, const ring_matrix< TYPE > &, const ring_matrix< TYPE > &);
template void subtract(ring_matrix< TYPE > &, const ring_matrix< TYPE > &, const TYPE &);
template void subtract(ring_matrix< TYPE > &, const TYPE &, const ring_matrix< TYPE > &);
template void multiply(ring_matrix< TYPE > &, const ring_matrix< TYPE > &, const ring_matrix< TYPE > &);
template void multiply(ring_matrix< TYPE > &, const ring_matrix< TYPE > &, const TYPE &);
template void multiply(ring_matrix< TYPE > &, const TYPE &, const ring_matrix< TYPE > &);
template void compwise_multiply(ring_matrix< TYPE > &, const ring_matrix< TYPE > &, const ring_matrix< TYPE > &);
template void multiply(TYPE *&, const ring_matrix< TYPE > &, const TYPE *);
template void multiply(TYPE *&, const TYPE *, const ring_matrix< TYPE > &);
template void multiply(math_vector< TYPE > &, const ring_matrix< TYPE > &, const math_vector< TYPE > &);
template void multiply(math_vector< TYPE > &, const math_vector< TYPE > &, const ring_matrix< TYPE > &);
template void negate(ring_matrix< TYPE > &, const ring_matrix< TYPE > &);
template ring_matrix< TYPE > operator + (const ring_matrix< TYPE > &, const ring_matrix< TYPE > &);
template ring_matrix< TYPE > operator + (const ring_matrix< TYPE > &, const TYPE &);
template ring_matrix< TYPE > operator + (const TYPE &, const ring_matrix< TYPE > &);
template ring_matrix< TYPE > & operator += (ring_matrix< TYPE > &, const ring_matrix< TYPE > &);
template ring_matrix< TYPE > & operator += (ring_matrix< TYPE > &, const TYPE &);
template ring_matrix< TYPE > operator - (const ring_matrix< TYPE > &, const ring_matrix< TYPE > &);
template ring_matrix< TYPE > operator - (const ring_matrix< TYPE > &, const TYPE &);
template ring_matrix< TYPE > operator - (const TYPE &, const ring_matrix< TYPE > &);
template ring_matrix< TYPE > & operator -= (ring_matrix< TYPE > &, const ring_matrix< TYPE > &);
template ring_matrix< TYPE > & operator -= (ring_matrix< TYPE > &, const TYPE &);
template ring_matrix< TYPE > operator * (const ring_matrix< TYPE > &, const ring_matrix< TYPE > &);
template ring_matrix< TYPE > operator * (const ring_matrix< TYPE > &, const TYPE &);
template ring_matrix< TYPE > operator * (const TYPE &, const ring_matrix< TYPE > &);
template TYPE * operator * (const ring_matrix< TYPE > &, const TYPE *);
template TYPE * operator * (const TYPE *, const ring_matrix< TYPE > &);
template math_vector< TYPE > operator * (const ring_matrix< TYPE > &, const math_vector< TYPE > &);
template math_vector< TYPE > operator * (const math_vector< TYPE > &, const ring_matrix< TYPE > &);
template ring_matrix< TYPE > & operator *= (ring_matrix< TYPE > &, const ring_matrix< TYPE > &);
template ring_matrix< TYPE > & operator *= (ring_matrix< TYPE > &, const TYPE &);
template ring_matrix< TYPE > operator - (const ring_matrix< TYPE > &);
template bool operator == (const ring_matrix< TYPE > &, const ring_matrix< TYPE > &);
template bool equal(const ring_matrix< TYPE > &, const ring_matrix< TYPE > &);
template bool operator != (const ring_matrix< TYPE > &, const ring_matrix< TYPE > &);
template bool unequal(const ring_matrix< TYPE > &, const ring_matrix< TYPE > &);
template TYPE trace(const ring_matrix< TYPE > &);

template class MR< TYPE >;



#  ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#  endif

# endif	// NORMAL



# ifdef DENSE
#  include	"LiDIA/dense_ring_matrix.h"
#  include	"LiDIA/dense_ring_matrix.cc"
#  include	"LiDIA/matrix/dense_ring_matrix_kernel.h"
#  include	"LiDIA/matrix/dense_ring_matrix_kernel.cc"

#  ifdef LIDIA_NAMESPACE
namespace LiDIA {
#  endif



template class dense_ring_matrix< TYPE >;

template void add(dense_ring_matrix< TYPE > &, const dense_ring_matrix< TYPE > &, const dense_ring_matrix< TYPE > &);
template void add(dense_ring_matrix< TYPE > &, const dense_ring_matrix< TYPE > &, const TYPE &);
template void add(dense_ring_matrix< TYPE > &, const TYPE &, const dense_ring_matrix< TYPE > &);
template void subtract(dense_ring_matrix< TYPE > &, const dense_ring_matrix< TYPE > &, const dense_ring_matrix< TYPE > &);
template void subtract(dense_ring_matrix< TYPE > &, const dense_ring_matrix< TYPE > &, const TYPE &);
template void subtract(dense_ring_matrix< TYPE > &, const TYPE &, const dense_ring_matrix< TYPE > &);
template void multiply(dense_ring_matrix< TYPE > &, const dense_ring_matrix< TYPE > &, const dense_ring_matrix< TYPE > &);
template void multiply(dense_ring_matrix< TYPE > &, const dense_ring_matrix< TYPE > &, const TYPE &);
template void multiply(dense_ring_matrix< TYPE > &, const TYPE &, const dense_ring_matrix< TYPE > &);
template void compwise_multiply(dense_ring_matrix< TYPE > &, const dense_ring_matrix< TYPE > &, const dense_ring_matrix< TYPE > &);
template void multiply(TYPE *&, const dense_ring_matrix< TYPE > &, const TYPE *);
template void multiply(TYPE *&, const TYPE *, const dense_ring_matrix< TYPE > &);
template void multiply(math_vector< TYPE > &, const dense_ring_matrix< TYPE > &, const math_vector< TYPE > &);
template void multiply(math_vector< TYPE > &, const math_vector< TYPE > &, const dense_ring_matrix< TYPE > &);
template void negate(dense_ring_matrix< TYPE > &, const dense_ring_matrix< TYPE > &);
template dense_ring_matrix< TYPE > operator + (const dense_ring_matrix< TYPE > &, const dense_ring_matrix< TYPE > &);
template dense_ring_matrix< TYPE > operator + (const dense_ring_matrix< TYPE > &, const TYPE &);
template dense_ring_matrix< TYPE > operator + (const TYPE &, const dense_ring_matrix< TYPE > &);
template dense_ring_matrix< TYPE > & operator += (dense_ring_matrix< TYPE > &, const dense_ring_matrix< TYPE > &);
template dense_ring_matrix< TYPE > & operator += (dense_ring_matrix< TYPE > &, const TYPE &);
template dense_ring_matrix< TYPE > operator - (const dense_ring_matrix< TYPE > &, const dense_ring_matrix< TYPE > &);
template dense_ring_matrix< TYPE > operator - (const dense_ring_matrix< TYPE > &, const TYPE &);
template dense_ring_matrix< TYPE > operator - (const TYPE &, const dense_ring_matrix< TYPE > &);
template dense_ring_matrix< TYPE > & operator -= (dense_ring_matrix< TYPE > &, const dense_ring_matrix< TYPE > &);
template dense_ring_matrix< TYPE > & operator -= (dense_ring_matrix< TYPE > &, const TYPE &);
template dense_ring_matrix< TYPE > operator * (const dense_ring_matrix< TYPE > &, const dense_ring_matrix< TYPE > &);
template dense_ring_matrix< TYPE > operator * (const dense_ring_matrix< TYPE > &, const TYPE &);
template dense_ring_matrix< TYPE > operator * (const TYPE &, const dense_ring_matrix< TYPE > &);
template TYPE * operator * (const dense_ring_matrix< TYPE > &, const TYPE *);
template TYPE * operator * (const TYPE *, const dense_ring_matrix< TYPE > &);
template math_vector< TYPE > operator * (const dense_ring_matrix< TYPE > &, const math_vector< TYPE > &);
template math_vector< TYPE > operator * (const math_vector< TYPE > &, const dense_ring_matrix< TYPE > &);
template dense_ring_matrix< TYPE > & operator *= (dense_ring_matrix< TYPE > &, const dense_ring_matrix< TYPE > &);
template dense_ring_matrix< TYPE > & operator *= (dense_ring_matrix< TYPE > &, const TYPE &);
template dense_ring_matrix< TYPE > operator - (const dense_ring_matrix< TYPE > &);
template bool operator == (const dense_ring_matrix< TYPE > &, const dense_ring_matrix< TYPE > &);
template bool equal(const dense_ring_matrix< TYPE > &, const dense_ring_matrix< TYPE > &);
template bool operator != (const dense_ring_matrix< TYPE > &, const dense_ring_matrix< TYPE > &);
template bool unequal(const dense_ring_matrix< TYPE > &, const dense_ring_matrix< TYPE > &);
template TYPE trace(const dense_ring_matrix< TYPE > &);



#  ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#  endif

# endif	// DENSE



# ifdef SPARSE
#  include	"LiDIA/sparse_ring_matrix.h"
#  include	"LiDIA/sparse_ring_matrix.cc"
#  include	"LiDIA/matrix/ring_matrix_algorithms.h"
#  include	"LiDIA/matrix/ring_matrix_algorithms.cc"
#  include	"LiDIA/matrix/sparse_ring_matrix_kernel.h"
#  include	"LiDIA/matrix/sparse_ring_matrix_kernel.cc"

#  ifdef LIDIA_NAMESPACE
namespace LiDIA {
#  endif



template class sparse_ring_matrix< TYPE >;

template void add(sparse_ring_matrix< TYPE > &, const sparse_ring_matrix< TYPE > &, const sparse_ring_matrix< TYPE > &);
template void add(sparse_ring_matrix< TYPE > &, const sparse_ring_matrix< TYPE > &, const TYPE &);
template void add(sparse_ring_matrix< TYPE > &, const TYPE &, const sparse_ring_matrix< TYPE > &);
template void subtract(sparse_ring_matrix< TYPE > &, const sparse_ring_matrix< TYPE > &, const sparse_ring_matrix< TYPE > &);
template void subtract(sparse_ring_matrix< TYPE > &, const sparse_ring_matrix< TYPE > &, const TYPE &);
template void subtract(sparse_ring_matrix< TYPE > &, const TYPE &, const sparse_ring_matrix< TYPE > &);
template void multiply(sparse_ring_matrix< TYPE > &, const sparse_ring_matrix< TYPE > &, const sparse_ring_matrix< TYPE > &);
template void multiply(sparse_ring_matrix< TYPE > &, const sparse_ring_matrix< TYPE > &, const TYPE &);
template void multiply(sparse_ring_matrix< TYPE > &, const TYPE &, const sparse_ring_matrix< TYPE > &);
template void compwise_multiply(sparse_ring_matrix< TYPE > &, const sparse_ring_matrix< TYPE > &, const sparse_ring_matrix< TYPE > &);
template void multiply(TYPE *&, const sparse_ring_matrix< TYPE > &, const TYPE *);
template void multiply(TYPE *&, const TYPE *, const sparse_ring_matrix< TYPE > &);
template void multiply(math_vector< TYPE > &, const sparse_ring_matrix< TYPE > &, const math_vector< TYPE > &);
template void multiply(math_vector< TYPE > &, const math_vector< TYPE > &, const sparse_ring_matrix< TYPE > &);
template void negate(sparse_ring_matrix< TYPE > &, const sparse_ring_matrix< TYPE > &);
template sparse_ring_matrix< TYPE > operator + (const sparse_ring_matrix< TYPE > &, const sparse_ring_matrix< TYPE > &);
template sparse_ring_matrix< TYPE > operator + (const sparse_ring_matrix< TYPE > &, const TYPE &);
template sparse_ring_matrix< TYPE > operator + (const TYPE &, const sparse_ring_matrix< TYPE > &);
template sparse_ring_matrix< TYPE > & operator += (sparse_ring_matrix< TYPE > &, const sparse_ring_matrix< TYPE > &);
template sparse_ring_matrix< TYPE > & operator += (sparse_ring_matrix< TYPE > &, const TYPE &);
template sparse_ring_matrix< TYPE > operator - (const sparse_ring_matrix< TYPE > &, const sparse_ring_matrix< TYPE > &);
template sparse_ring_matrix< TYPE > operator - (const sparse_ring_matrix< TYPE > &, const TYPE &);
template sparse_ring_matrix< TYPE > operator - (const TYPE &, const sparse_ring_matrix< TYPE > &);
template sparse_ring_matrix< TYPE > & operator -= (sparse_ring_matrix< TYPE > &, const sparse_ring_matrix< TYPE > &);
template sparse_ring_matrix< TYPE > & operator -= (sparse_ring_matrix< TYPE > &, const TYPE &);
template sparse_ring_matrix< TYPE > operator * (const sparse_ring_matrix< TYPE > &, const sparse_ring_matrix< TYPE > &);
template sparse_ring_matrix< TYPE > operator * (const sparse_ring_matrix< TYPE > &, const TYPE &);
template sparse_ring_matrix< TYPE > operator * (const TYPE &, const sparse_ring_matrix< TYPE > &);
template TYPE * operator * (const sparse_ring_matrix< TYPE > &, const TYPE *);
template TYPE * operator * (const TYPE *, const sparse_ring_matrix< TYPE > &);
template math_vector< TYPE > operator * (const sparse_ring_matrix< TYPE > &, const math_vector< TYPE > &);
template math_vector< TYPE > operator * (const math_vector< TYPE > &, const sparse_ring_matrix< TYPE > &);
template sparse_ring_matrix< TYPE > & operator *= (sparse_ring_matrix< TYPE > &, const sparse_ring_matrix< TYPE > &);
template sparse_ring_matrix< TYPE > & operator *= (sparse_ring_matrix< TYPE > &, const TYPE &);
template sparse_ring_matrix< TYPE > operator - (const sparse_ring_matrix< TYPE > &);
template bool operator == (const sparse_ring_matrix< TYPE > &, const sparse_ring_matrix< TYPE > &);
template bool equal(const sparse_ring_matrix< TYPE > &, const sparse_ring_matrix< TYPE > &);
template bool operator != (const sparse_ring_matrix< TYPE > &, const sparse_ring_matrix< TYPE > &);
template bool unequal(const sparse_ring_matrix< TYPE > &, const sparse_ring_matrix< TYPE > &);
template TYPE trace(const sparse_ring_matrix< TYPE > &);



#  ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#  endif

# endif	// SPARSE

#endif	// RING_MATRIX



#ifdef FIELD_MATRIX

# ifdef SUBCLASSES
#  include	"LiDIA/field_matrix.h"
#  include	"LiDIA/field_matrix.cc"
#  include	"LiDIA/matrix/field_matrix_algorithms.h"
#  include	"LiDIA/matrix/field_matrix_algorithms.cc"
#  include	"LiDIA/matrix/dense_field_matrix_kernel.h"
#  include	"LiDIA/matrix/dense_field_matrix_kernel.cc"
#  include	"LiDIA/matrix/sparse_field_matrix_kernel.h"
#  include	"LiDIA/matrix/sparse_field_matrix_kernel.cc"

#  define DFMKex DFMK< TYPE >
#  define SFMKex SFMK< TYPE >

#  ifdef LIDIA_NAMESPACE
namespace LiDIA {
#  endif



template class DFMKex;
template class SFMKex;

template class FMA< TYPE, DFMKex, DFMKex, DFMKex >;
template class FMA< TYPE, DFMKex, DFMKex, SFMKex >;
template class FMA< TYPE, DFMKex, SFMKex, DFMKex >;
template class FMA< TYPE, SFMKex, DFMKex, DFMKex >;
template class FMA< TYPE, DFMKex, SFMKex, SFMKex >;
template class FMA< TYPE, SFMKex, DFMKex, SFMKex >;
template class FMA< TYPE, SFMKex, SFMKex, DFMKex >;
template class FMA< TYPE, SFMKex, SFMKex, SFMKex >;



#  ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#  endif

#  undef DFMKex
#  undef SFMKex

# endif	// SUBCLASSES



# ifdef NORMAL
#  include	"LiDIA/field_matrix.h"
#  include	"LiDIA/field_matrix.cc"
#  include	"LiDIA/matrix/field_matrix_algorithms.h"
#  include	"LiDIA/matrix/field_matrix_algorithms.cc"
#  include	"LiDIA/matrix/dense_field_matrix_kernel.h"
#  include	"LiDIA/matrix/dense_field_matrix_kernel.cc"
#  include	"LiDIA/matrix/sparse_field_matrix_kernel.h"
#  include	"LiDIA/matrix/sparse_field_matrix_kernel.cc"

#  ifdef LIDIA_NAMESPACE
namespace LiDIA {
#  endif



template class field_matrix< TYPE >;

template void divide(field_matrix< TYPE > &, const field_matrix< TYPE > &, const TYPE &);
template void compwise_divide(field_matrix< TYPE > &, const field_matrix< TYPE > &, const field_matrix< TYPE > &);
template field_matrix< TYPE > operator / (const field_matrix< TYPE > &, const TYPE &);
template field_matrix< TYPE > & operator /= (field_matrix< TYPE > &, const TYPE &);



#  ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#  endif

# endif	// NORMAL



# ifdef DENSE
#  include	"LiDIA/dense_field_matrix.h"
#  include	"LiDIA/dense_field_matrix.cc"
#  include	"LiDIA/matrix/dense_field_matrix_kernel.h"
#  include	"LiDIA/matrix/dense_field_matrix_kernel.cc"

#  ifdef LIDIA_NAMESPACE
namespace LiDIA {
#  endif



template class dense_field_matrix< TYPE >;

template void divide(dense_field_matrix< TYPE > &, const dense_field_matrix< TYPE > &, const TYPE &);
template void compwise_divide(dense_field_matrix< TYPE > &, const dense_field_matrix< TYPE > &, const dense_field_matrix< TYPE > &);
template dense_field_matrix< TYPE > operator / (const dense_field_matrix< TYPE > &, const TYPE &);
template dense_field_matrix< TYPE > & operator /= (dense_field_matrix< TYPE > &, const TYPE &);



#  ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#  endif

# endif	// DENSE



# ifdef SPARSE
#  include	"LiDIA/sparse_field_matrix.h"
#  include	"LiDIA/sparse_field_matrix.cc"
#  include	"LiDIA/matrix/field_matrix_algorithms.h"
#  include	"LiDIA/matrix/field_matrix_algorithms.cc"
#  include	"LiDIA/matrix/sparse_field_matrix_kernel.h"
#  include	"LiDIA/matrix/sparse_field_matrix_kernel.cc"

#  ifdef LIDIA_NAMESPACE
namespace LiDIA {
#  endif



template class sparse_field_matrix< TYPE >;

template void divide(sparse_field_matrix< TYPE > &, const sparse_field_matrix< TYPE > &, const TYPE &);
template void compwise_divide(sparse_field_matrix< TYPE > &, const sparse_field_matrix< TYPE > &, const sparse_field_matrix< TYPE > &);
template sparse_field_matrix< TYPE > operator / (const sparse_field_matrix< TYPE > &, const TYPE &);
template sparse_field_matrix< TYPE > & operator /= (sparse_field_matrix< TYPE > &, const TYPE &);



#  ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#  endif

# endif	// SPARSE

#endif	// FIELD_MATRIX
