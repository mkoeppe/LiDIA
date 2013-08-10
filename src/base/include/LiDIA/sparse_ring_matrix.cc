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


#ifndef LIDIA_SPARSE_RING_MATRIX_CC_GUARD_
#define LIDIA_SPARSE_RING_MATRIX_CC_GUARD_



#ifndef LIDIA_SPARSE_RING_MATRIX_H_GUARD_
# include	"LiDIA/sparse_ring_matrix.h"
#endif
#ifndef LIDIA_ERROR_H_GUARD_
# include	"LiDIA/error.h"
#endif
#include "LiDIA/precondition_error.h"



#ifdef LIDIA_NAMESPACE
# ifndef IN_NAMESPACE_LIDIA
namespace LiDIA {
# endif
#endif



//
// debug defines / error defines
//

extern const char *PRT;
extern const char *matrix_error_msg[];

#define DVALUE LDBL_MATRIX + 10          // Debug value
#define DMESSAGE "sparse_ring_matrix"    // Debug message
#define EMESSAGE matrix_error_msg        // Error message array

//
// debug level
//
//   0 : addtion
//   1 : subtraction
//   2 : multiplication
//   3 : division
//   4 : negation
//   5 : comparison
//   6 : trace
//

//////////////////////////////////
// BEGIN: arithmetic procedures //
//////////////////////////////////

//
// addition
//

template< class T >
void sparse_ring_matrix< T >::
add(const sparse_ring_matrix< T > &M, const sparse_ring_matrix< T > &N)
{
	//
	//       Task: RES.add(M, N)
	// =  > RES = M + N,
	// Conditions: M.rows == N.rows and
	//             M.columns == N.columns
	//    Version: 2.0
	//

	debug_handler_l(DMESSAGE, "add(const sparse_ring_matrix< T > &, "
			"const sparse_ring_matrix< T > &)", DVALUE);

	if (M.rows != N.rows || M.columns != N.columns)
		precondition_error_handler(M.rows, "M.rows", "M.rows == N.rows",
				    N.rows, "N.rows", "M.rows == N.rows",
				    M.columns, "M.columns", "M.columns == N.columns",
				    N.columns, "N.columns", "M.columns == N.columns",
				    "void sparse_ring_matrix< T >::"
				    "add(const sparse_ring_matrix< T > &M, "
				    "const sparse_ring_matrix< T > &N)",
				    DMESSAGE, EMESSAGE[4]);

	if (this->rows != N.rows)
		this->S_ring_modul.set_no_of_rows(*this, N.rows);
	if (this->columns != N.columns)
		this->S_ring_modul.set_no_of_columns(*this, N.columns);

	this->S_ring_modul.add(*this, M, N);
}



template< class T >
void sparse_ring_matrix< T >::
add(const sparse_ring_matrix< T > &M, const T &a)
{
	//
	//    Task: RES.add(M, a);
	// =  > RES.value[x][y] = M.value[x][y] + a,
	//             x = 0, ..., M.rows-1, y = 0, ..., M.columns-1
	// Version: 2.0
	//

	debug_handler_l(DMESSAGE, "add(const sparse_ring_matrix< T > &, const T &)",
			DVALUE);

	if (this->rows != M.rows)
		this->S_ring_modul.set_no_of_rows(*this, M.rows);
	if (this->columns != M.columns)
		this->S_ring_modul.set_no_of_columns(*this, M.columns);

	this->S_ring_modul.add(*this, M, a);
}



template< class T >
void sparse_ring_matrix< T >::
add(const T &a, const sparse_ring_matrix< T > &M)
{
	//
	//    Task: RES.add(a, M);
	// =  > RES.value[x][y] = a + M.value[x][y],
	//             x = 0, ..., M.rows-1, y = 0, ..., M.columns-1
	// Version: 2.0
	//

	debug_handler_l(DMESSAGE, "add(const T &, const sparse_ring_matrix< T > &)",
			DVALUE);

	if (this->rows != M.rows)
		this->S_ring_modul.set_no_of_rows(*this, M.rows);
	if (this->columns != M.columns)
		this->S_ring_modul.set_no_of_columns(*this, M.columns);

	this->S_ring_modul.add(*this, a, M);
}



//
// subtraction
//

template< class T >
void sparse_ring_matrix< T >::
subtract(const sparse_ring_matrix< T > &M, const sparse_ring_matrix< T > &N)
{
	//
	//       Task: RES.subtract(M, N)
	// =  > RES = M - N
	// Conditions: M.rows == N.rows and
	//             M.columns == N.columns
	//    Version: 2.0
	//

	debug_handler_l(DMESSAGE, "subtract(const sparse_ring_matrix< T > &, "
			"const sparse_ring_matrix< T > &)", DVALUE + 1);

	if (M.rows != N.rows || M.columns != N.columns)
		precondition_error_handler(M.rows, "M.rows", "M.rows == N.rows",
				    N.rows, "N.rows", "M.rows == N.rows",
				    M.columns, "M.columns", "M.columns == N.columns",
				    N.columns, "N.columns", "M.columns == N.columns",
				    "void sparse_ring_matrix< T >::"
				    "subtract(const sparse_ring_matrix< T > &M, "
				    "const sparse_ring_matrix< T > &N)",
				    DMESSAGE, EMESSAGE[4]);
	if (this->rows != N.rows)
		this->S_ring_modul.set_no_of_rows(*this, N.rows);
	if (this->columns != N.columns)
		this->S_ring_modul.set_no_of_columns(*this, N.columns);

	this->S_ring_modul.subtract(*this, M, N);
}



template< class T >
void sparse_ring_matrix< T >::
subtract(const sparse_ring_matrix< T > &M, const T &a)
{
	//
	//    Task: RES.subtract(M, a)
	// =  > RES.value[x][y] = M.value[x][y] - a,
	//             x = 0, ..., M.rows-1, y = 0, ..., M.columns-1
	// Version: 2.0
	//

	debug_handler_l(DMESSAGE,
			"subtract(const sparse_ring_matrix< T > &, const T &)",
			DVALUE + 1);

	if (this->rows != M.rows)
		this->S_ring_modul.set_no_of_rows(*this, M.rows);
	if (this->columns != M.columns)
		this->S_ring_modul.set_no_of_columns(*this, M.columns);

	this->S_ring_modul.subtract(*this, M, a);
}



template< class T >
void sparse_ring_matrix< T >::
subtract(const T &a, const sparse_ring_matrix< T > &M)
{
	//
	//    Task: RES.subtract(a, M)
	// =  > RES.value[x][y] = a - M.value[x][y],
	//             x = 0, ..., M.rows-1, y = 0, ..., M.columns-1
	// Version: 2.0
	//

	debug_handler_l(DMESSAGE, "subtract(const T &, const sparse_ring_matrix< T > &)",
			DVALUE + 1);

	if (this->rows != M.rows)
		this->S_ring_modul.set_no_of_rows(*this, M.rows);
	if (this->columns != M.columns)
		this->S_ring_modul.set_no_of_columns(*this, M.columns);

	this->S_ring_modul.subtract(*this, a, M);
}



//
// multiplication
//

template< class T >
void sparse_ring_matrix< T >::
multiply(const sparse_ring_matrix< T > &A, const sparse_ring_matrix< T > &B)
{
	//
	//       Task: RES.multiply(A, B)
	// =  > RES = A * B
	// Conditions: A.columns == B.rows
	//    Version: 2.0
	//

	debug_handler_l(DMESSAGE, "multiply(const sparse_ring_matrix< T > &, "
			"const sparse_ring_matrix< T > &)", DVALUE + 2);

	if (A.columns != B.rows)
		precondition_error_handler(A.columns, "A.columns", "A.columns == B.rows",
				    B.rows, "B.rows", "A.columns == B.rows",
				    "void sparse_ring_matrix< T >::"
				    "multiply(const sparse_ring_matrix< T > &A, "
				    "const sparse_ring_matrix< T > &B)",
				    DMESSAGE, EMESSAGE[4]);

	if (A.value != this->value && B.value != this->value) {
		if (this->rows != A.rows)
			this->S_ring_modul.set_no_of_rows(*this, A.rows);
		if (this->columns != B.columns)
			this->S_ring_modul.set_no_of_columns(*this, B.columns);
		this->SSS_ring_modul.multiply(*this, A, B);
	}
	else {
		sparse_ring_matrix< T > RES1(A.rows, B.columns);
		this->SSS_ring_modul.multiply(RES1, A, B);
		if (this->rows != A.rows)
			this->S_ring_modul.set_no_of_rows(*this, A.rows);
		if (this->columns != B.columns)
			this->S_ring_modul.set_no_of_columns(*this, B.columns);
		this->S_ring_modul.assign(*this, RES1);
	}
}



template< class T >
void sparse_ring_matrix< T >::
multiply(const sparse_ring_matrix< T > &A, const T &k)
{
	//
	//    Task: RES.multiply(A, k)
	// =  > RES.value[x][y] = A.value[x][y]*k,
	//             x = 0, ..., A.rows-1, y = 0, ..., A.columns-1
	// Version: 2.0
	//

	debug_handler_l(DMESSAGE, "multiply(const sparse_ring_matrix< T > &, const T &)",
			DVALUE + 2);

	if (this->rows != A.rows)
		this->S_ring_modul.set_no_of_rows(*this, A.rows);
	if (this->columns != A.columns)
		this->S_ring_modul.set_no_of_columns(*this, A.columns);

	this->SSS_ring_modul.multiply(*this, A, k);
}



template< class T >
void sparse_ring_matrix< T >::
multiply(const T &k, const sparse_ring_matrix< T > &A)
{
	//
	//    Task: RES.multiply(k, A)
	// =  > RES.value[x][y] = k*A.value[x][y],
	//             x = 0, ..., A.rows-1, y = 0, ..., A.columns-1
	// Version: 2.0
	//

	debug_handler_l(DMESSAGE, "multiply(const T &, const sparse_ring_matrix< T > &)",
			DVALUE + 2);

	if (this->rows != A.rows)
		this->S_ring_modul.set_no_of_rows(*this, A.rows);
	if (this->columns != A.columns)
		this->S_ring_modul.set_no_of_columns(*this, A.columns);

	this->SSS_ring_modul.multiply(*this, k, A);
}



template< class T >
void sparse_ring_matrix< T >::
compwise_multiply(const sparse_ring_matrix< T > &A, const sparse_ring_matrix< T > &B)
{
	//
	//    Task: RES.compwise_multiply(A, B)
	// =  > RES.value[x][y] = A.value[x][y]*B.value[x][y],
	//             x = 0, ..., A.rows-1, y = 0, ..., A.columns-1
	// Version: 2.0
	//

	debug_handler_l(DMESSAGE, "compwise_multiply(const sparse_ring_matrix< T > &, "
			"const sparse_ring_matrix< T > &)", DVALUE + 2);

	if (A.rows != B.rows || A.columns != B.columns)
		precondition_error_handler(A.rows, "A.rows", "A.rows == B.rows",
				    B.rows, "B.rows", "A.rows == B.rows",
				    A.columns, "A.columns", "A.columns == B.columns",
				    B.columns, "B.columns", "B.columns == B.columns",
				    "void sparse_ring_matrix< T >::"
				    "compwise_multiply(const sparse_ring_matrix< T > &A, "
				    "const sparse_ring_matrix< T > &B)",
				    DMESSAGE, EMESSAGE[4]);

	if (this->rows != A.rows)
		this->S_ring_modul.set_no_of_rows(*this, A.rows);
	if (this->columns != A.columns)
		this->S_ring_modul.set_no_of_columns(*this, A.columns);

	this->SSS_ring_modul.compwise_multiply(*this, A, B);
}



template< class T >
void sparse_ring_matrix< T >::
multiply_right(math_vector< T > &res, const math_vector< T > &v) const
{
	//
	//       Task: c[x] = A.value[x][0]*v[0]+...
	//                      +A.value[x][A.columns-1]*v[A.columns-1],
	//             x = 0, ..., A.rows-1
	// Conditions: v.size == A.columns
	//    Version: 2.0
	//

	debug_handler_l(DMESSAGE, "multiply_right(math_vector< T > &, "
			"const math_vector< T > &)", DVALUE + 2);

	if (v.size() != this->columns)
		precondition_error_handler(this->columns, "columns", "columns == v.size",
				    v.size(), "v.size", "columns == v.size",
				    "void sparse_ring_matrix< T >::"
				    "multiply_right(math_vector< T > &res, "
				    "const math_vector< T > &v) const",
				    DMESSAGE, EMESSAGE[1]);

	this->SSS_ring_modul.multiply_right(*this, res, v);
}



template< class T >
void sparse_ring_matrix< T >::
multiply_right(T *&c, const T *v) const
{
	//
	//       Task: c[x] = A.value[x][0]*v[0]+...
	//                           +A.value[x][A.columns-1]*v[A.columns-1],
	//             x = 0, ..., A.rows-1
	// Conditions: v != NULL
	//    Version: 2.0
	//

	debug_handler_l(DMESSAGE, "multiply_right(T *&, const T *)", DVALUE + 2);

	if (v == NULL)
		precondition_error_handler(PRT, "v", "v != NULL",
				    "void sparse_ring_matrix< T >::"
				    "multiply_right(T *&c, const T *v) const",
				    DMESSAGE, EMESSAGE[6]);

	this->SSS_ring_modul.multiply_right(*this, c, v);
}



template< class T >
void sparse_ring_matrix< T >::
multiply_left(math_vector< T > &res, const math_vector< T > &v) const
{
	//
	//       Task: c[x] = v[0]*A.value[0][x]+...
	//                         +v[A.rows-1]*A.value[A.rows-1][x],
	//             x = 0, ..., A.columns-1
	// Conditions: v.length == A.rows
	//    Version: 2.0
	//

	debug_handler_l(DMESSAGE, "multiply_left(math_vector< T > &, const math_vector< T > &)",
			DVALUE + 2);

	if (v.size() != this->rows)
		precondition_error_handler(this->rows, "rows", "rows == v.size",
				    v.size(), "v.size", "rows == v.size",
				    "void sparse_ring_matrix< T >::"
				    "multiply_left(math_vector< T > &res, "
				    "const math_vector< T > &v) const",
				    DMESSAGE, EMESSAGE[6]);

	this->SSS_ring_modul.multiply_left(*this, res, v);
}



template< class T >
void sparse_ring_matrix< T >::
multiply_left(T *&c, const T *v) const
{
	//
	//       Task: c[x] = v[0]*A.value[0][x]+...
	//                           +v[A.rows-1]*A.value[A.rows-1][x],
	//             x = 0, ..., A.columns-1
	// Conditions: v != NULL
	//    Version: 2.0
	//

	debug_handler_l(DMESSAGE, "multiply_left(T *&, const T *)", DVALUE + 2);

	if (v == NULL)
		precondition_error_handler(PRT, "v", "v != NULL",
				    "void sparse_ring_matrix< T >::"
				    "multiply_left(T *&c, const T *v) const",
				    DMESSAGE, EMESSAGE[6]);

	this->SSS_ring_modul.multiply_left(*this, c, v);
}



//
// negation
//

template< class T >
void sparse_ring_matrix< T >::
negate(const sparse_ring_matrix< T > &B)
{
	//
	//    Task: A.negate(B)
	// =  > A.value[x][y] = - B.value[x][y],
	//             x = 0, ..., rows-1, y = 0, ..., columns-1
	// Version: 2.0
	//

	debug_handler_l(DMESSAGE, "negate(const sparse_ring_matrix< T > &)",
			DVALUE + 4);

	if (this->rows != B.rows)
		this->S_ring_modul.set_no_of_rows(*this, B.rows);
	if (this->columns != B.columns)
		this->S_ring_modul.set_no_of_columns(*this, B.columns);

	this->SSS_ring_modul.negate(*this, B);
}



//////////////////////////////////
// END: arithmetical procedures //
//////////////////////////////////

//
// comparisons
//

template< class T >
bool sparse_ring_matrix< T >::
equal(const sparse_ring_matrix< T > &N) const
{
	//
	//    Task: A.equal(B) returns true, if A.value[x][y] == B.value[x][y],
	//          x = 0, ..., A.rows-1, y = 0, ..., A.columns-1
	//          and A.columns == B.columns and A.rows == B.rows
	// Version: 2.0
	//

	debug_handler_l(DMESSAGE, "equal(const sparse_ring_matrix< T > &)",
			DVALUE + 5);

	return this->SSS_ring_modul.equal(*this, N);
}



//
// trace
//

template< class T >
void sparse_ring_matrix< T >::
trace(T &tr) const
{
	//
	//       Task: A.trace(tr)
	// =  > tr = trace of matrix A
	// Conditions: A.columns == A.rows
	// Version: 2.0
	//

	debug_handler_l(DMESSAGE, "trace(T &)", DVALUE + 6);

	if (this->rows != this->columns)
		precondition_error_handler(this->rows, "rows", "rows == columns",
				    this->columns, "columns", "rows == columns",
				    "void sparse_ring_matrix< T >::"
				    "trace(T &tr) const",
				    DMESSAGE, EMESSAGE[7]);

	this->SSS_ring_modul.trace(*this, tr);
}



//
// size reduction
//

template< class T >
void sparse_ring_matrix< T >::
size_reduction()
{
	//
	//    Task: A.size_reduction()
	//          reduces the no of rows according to the structured gauss elimination
	// Version: 2.0
	//

	debug_handler_l(DMESSAGE, "size_reduction()", DVALUE + 6);

	this->SSS_ring_modul.size_reduction(*this);
}



#undef DVALUE
#undef DMESSAGE
#undef EMESSAGE



#ifdef LIDIA_NAMESPACE
# ifndef IN_NAMESPACE_LIDIA
}	// end of namespace LiDIA
# endif
#endif



#endif	// LIDIA_SPARSE_RING_MATRIX_CC_GUARD_
