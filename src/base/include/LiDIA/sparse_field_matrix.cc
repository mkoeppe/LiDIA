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


#ifndef LIDIA_SPARSE_FIELD_MATRIX_CC_GUARD_
#define LIDIA_SPARSE_FIELD_MATRIX_CC_GUARD_



#ifndef LIDIA_SPARSE_FIELD_MATRIX_H_GUARD_
# include	"LiDIA/sparse_field_matrix.h"
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
#define DMESSAGE "sparse_field_matrix"     // Debug message
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
// divide
//

template< class T >
void sparse_field_matrix< T >::
divide(const sparse_field_matrix< T > &A, const T &k)
{
	//
	//    Task: RES.divide(A, k)
	// =  > RES.value[x][y] = A.value[x][y] / k,
	//             x = 0, ..., A.rows-1, y = 0, ..., A.columns-1
	// Version: 2.0
	//

	debug_handler_l(DMESSAGE,
			"divide(const sparse_field_matrix< T > &, const T &)",
			DVALUE + 2);

	if (this->rows != A.rows)
		this->set_no_of_rows(A.rows);
	if (this->columns != A.columns)
		this->set_no_of_columns(A.columns);

	this->SSS_field_modul.divide(*this, A, k);
}



template< class T >
void sparse_field_matrix< T >::
compwise_divide(const sparse_field_matrix< T > &A, const sparse_field_matrix< T > &B)
{
	//
	//       Task: RES.compwise_divide(A, B)
	// =  > RES.value[x][y] = A.value[x][y] / B.value[x][y],
	//                x = 0, ..., A.rows-1, y = 0, ..., A.columns-1
	// Conditions: A.rows == A.rows and
	//             A.columns == B.columns
	//    Version: 2.0
	//

	debug_handler_l(DMESSAGE,
			"compwise_divide(const sparse_field_matrix< T > &, "
			"const sparse_field_matrix< T > &)",
			DVALUE + 3);

	if (A.rows != B.rows || A.columns != B.columns)
		precondition_error_handler(A.rows, "A.rows", "A.rows == B.rows",
				    B.rows, "B.rows", "A.rows == B.rows",
				    A.columns, "A.columns", "A.columns == B.columns",
				    B.columns, "B.columns", "B.columns == B.columns",
				    "void sparse_field_matrix< T >::"
				    "compwise_divide(const sparse_field_matrix< T > &A, "
				    "const sparse_field_matrix< T > &B)",
				    DMESSAGE, EMESSAGE[4]);

	if (this->rows != A.rows)
		this->S_field_modul.set_no_of_rows(*this, A.rows);
	if (this->columns != A.columns)
		this->S_field_modul.set_no_of_columns(*this, A.columns);

	this->SSS_field_modul.compwise_divide(*this, A, B);
}



//////////////////////////////////
// END: arithmetical procedures //
//////////////////////////////////



#undef DVALUE
#undef DMESSAGE
#undef EMESSAGE



#ifdef LIDIA_NAMESPACE
# ifndef IN_NAMESPACE_LIDIA
}	// end of namespace LiDIA
# endif
#endif



#endif	// LIDIA_SPARSE_FIELD_MATRIX_CC_GUARD_
