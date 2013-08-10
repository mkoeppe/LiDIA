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


#ifndef LIDIA_FIELD_MATRIX_CC_GUARD_
#define LIDIA_FIELD_MATRIX_CC_GUARD_



#ifndef LIDIA_FIELD_MATRIX_H_GUARD_
# include	"LiDIA/field_matrix.h"
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
#define DMESSAGE "field_matrix"     // Debug message
#define EMESSAGE matrix_error_msg        // Error message array

#define compute_switch(A, B, C, D, E) ((((A*10 + B)*10 + C)*10 + D)*10 + E)

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

//
// divide
//

template< class T >
void field_matrix< T >::
divide(const field_matrix< T > &M, const T &a)
{
	//
	//    Task: RES.divide(A, k)
	// =  > RES.value[x][y] = A.value[x][y] / k,
	//             x = 0, ..., A.rows-1, y = 0, ..., A.columns-1
	// Version: 2.0
	//

	debug_handler_l(DMESSAGE,
			"divide(const field_matrix< T > &, const T &)",
			DVALUE + 2);

	switch (this->bitfield.get_representation()) {
	case matrix_flags::dense_representation:
		if (this->rows != M.rows)
			D_field_modul.set_no_of_rows(*this, M.rows);
		if (this->columns != M.columns)
			D_field_modul.set_no_of_columns(*this, M.columns);

		switch (M.bitfield.get_representation()) {
		case matrix_flags::dense_representation:
			DDD_field_modul.divide(*this, M, a);
			break;
		case matrix_flags::sparse_representation:
			DSD_field_modul.divide(*this, M, a);
			break;
		}
		break;

	case matrix_flags::sparse_representation:
		if (this->rows != M.rows)
			S_field_modul.set_no_of_rows(*this, M.rows);
		if (this->columns != M.columns)
			S_field_modul.set_no_of_columns(*this, M.columns);

		switch (M.bitfield.get_representation()) {
		case matrix_flags::dense_representation:
			SDD_field_modul.divide(*this, M, a);
			break;
		case matrix_flags::sparse_representation:
			SSD_field_modul.divide(*this, M, a);
			break;
		}
		break;
	}
}



template< class T >
void field_matrix< T >::
compwise_divide(const field_matrix< T > &A, const field_matrix< T > &B)
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
			"compwise_divide(const field_matrix< T > &, "
			"const field_matrix< T > &)",
			DVALUE + 3);

	if (A.rows != B.rows || A.columns != B.columns)
		precondition_error_handler(A.rows, "A.rows", "A.rows == B.rows",
				    B.rows, "B.rows", "A.rows == B.rows",
				    A.columns, "A.columns", "A.columns == B.columns",
				    B.columns, "B.columns", "B.columns == B.columns",
				    "void field_matrix< T >::"
				    "compwise_divide(const field_matrix< T > &A, "
				    "const field_matrix< T > &B)",
				    DMESSAGE, EMESSAGE[4]);

	unsigned long SW = compute_switch(this->bitfield.get_representation(), A.bitfield.get_representation(),
					  B.bitfield.get_representation(), 0, 0);

	const unsigned long DR = matrix_flags::dense_representation;
	const unsigned long SR = matrix_flags::sparse_representation;

	switch (this->bitfield.get_representation()) {
	case matrix_flags::dense_representation:
		if (this->rows != A.rows)
			this->D_field_modul.set_no_of_rows(*this, A.rows);
		if (this->columns != B.columns)
			this->D_field_modul.set_no_of_columns(*this, B.columns);

		switch(SW) {
		case (compute_switch(DR, DR, DR, 0, 0)):
			this->DDD_field_modul.compwise_divide(*this, A, B);
			break;
		case (compute_switch(DR, DR, SR, 0, 0)):
			this->DDS_field_modul.compwise_divide(*this, A, B);
			break;
		case (compute_switch(DR, SR, DR, 0, 0)):
			this->DSD_field_modul.compwise_divide(*this, A, B);
			break;
		case (compute_switch(DR, SR, SR, 0, 0)):
			this->DSS_field_modul.compwise_divide(*this, A, B);
			break;
		}
		break;
	case matrix_flags::sparse_representation:
		if (this->rows != A.rows)
			this->S_field_modul.set_no_of_rows(*this, A.rows);
		if (this->columns != A.columns)
			this->S_field_modul.set_no_of_columns(*this, A.columns);

		switch(SW) {
		case (compute_switch(SR, DR, DR, 0, 0)):
			this->SDD_field_modul.compwise_divide(*this, A, B);
			break;
		case (compute_switch(SR, DR, SR, 0, 0)):
			this->SDS_field_modul.compwise_divide(*this, A, B);
			break;
		case (compute_switch(SR, SR, DR, 0, 0)):
			this->SSD_field_modul.compwise_divide(*this, A, B);
			break;
		case (compute_switch(SR, SR, SR, 0, 0)):
			this->SSS_field_modul.compwise_divide(*this, A, B);
			break;
		}
		break;
	}
}



#undef compute_switch

#undef DVALUE
#undef DMESSAGE
#undef EMESSAGE



#ifdef LIDIA_NAMESPACE
# ifndef IN_NAMESPACE_LIDIA
}	// end of namespace LiDIA
# endif
#endif



#endif	// LIDIA_FIELD_MATRIX_CC_GUARD_
