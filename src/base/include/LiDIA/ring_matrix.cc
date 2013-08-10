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


#ifndef LIDIA_RING_MATRIX_CC_GUARD_
#define LIDIA_RING_MATRIX_CC_GUARD_



#ifndef LIDIA_RING_MATRIX_H_GUARD_
# include	"LiDIA/ring_matrix.h"
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

#define DVALUE LDBL_MATRIX + 10    // Debug value
#define DMESSAGE "ring_matrix"     // Debug message
#define EMESSAGE matrix_error_msg  // Error message array

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

//////////////////////////////////
// BEGIN: arithmetic procedures //
//////////////////////////////////

//
// addition
//

template< class T >
void ring_matrix< T >::
add(const ring_matrix< T > &M, const ring_matrix< T > &N)
{
	//
	//       Task: RES.add(M, N)
	// =  > RES = M + N,
	// Conditions: M.rows == N.rows and
	//             M.columns == N.columns
	//    Version: 2.0
	//

	debug_handler_l(DMESSAGE, "add(const ring_matrix< T > &, "
			"const ring_matrix< T > &)", DVALUE);

	if (M.rows != N.rows || M.columns != N.columns)
		precondition_error_handler(M.rows, "M.rows", "M.rows == N.rows",
				    N.rows, "N.rows", "M.rows == N.rows",
				    M.columns, "M.columns", "M.columns == N.columns",
				    N.columns, "N.columns", "M.columns == N.columns",
				    "void ring_matrix< T >::"
				    "add(const ring_matrix< T > &M, "
				    "const ring_matrix< T > &N)",
				    DMESSAGE, EMESSAGE[4]);

	unsigned long SW = compute_switch(this->bitfield.get_representation(),
					  M.bitfield.get_representation(),
					  N.bitfield.get_representation(), 0, 0);

	const unsigned long DR = matrix_flags::dense_representation;
	const unsigned long SR = matrix_flags::sparse_representation;

	switch (this->bitfield.get_representation()) {
	case matrix_flags::dense_representation:
		if (this->rows != N.rows)
			this->D_ring_modul.set_no_of_rows(*this, N.rows);
		if (this->columns != N.columns)
			this->D_ring_modul.set_no_of_columns(*this, N.columns);

		switch(SW) {
		case (compute_switch(DR, DR, DR, 0, 0)):
			this->DDD_ring_modul.add(*this, M, N);
			break;
		case (compute_switch(DR, DR, SR, 0, 0)):
			this->DDS_ring_modul.add(*this, M, N);
			break;
		case (compute_switch(DR, SR, DR, 0, 0)):
			this->DSD_ring_modul.add(*this, M, N);
			break;
		case (compute_switch(DR, SR, SR, 0, 0)):
			this->DSS_ring_modul.add(*this, M, N);
			break;
		}
		break;
	case matrix_flags::sparse_representation:
		if (this->rows != N.rows)
			this->S_ring_modul.set_no_of_rows(*this, N.rows);
		if (this->columns != N.columns)
			this->S_ring_modul.set_no_of_columns(*this, N.columns);

		switch(SW) {
		case (compute_switch(SR, DR, DR, 0, 0)):
			this->SDD_ring_modul.add(*this, M, N);
			break;
		case (compute_switch(SR, DR, SR, 0, 0)):
			this->SDS_ring_modul.add(*this, M, N);
			break;
		case (compute_switch(SR, SR, DR, 0, 0)):
			this->SSD_ring_modul.add(*this, M, N);
			break;
		case (compute_switch(SR, SR, SR, 0, 0)):
			this->S_ring_modul.add(*this, M, N);
			break;
		}
		break;
	}
}



template< class T >
void ring_matrix< T >::
add(const ring_matrix< T > &M, const T &a)
{
	//
	//    Task: RES.add(M, a);
	// =  > RES.value[x][y] = M.value[x][y] + a,
	//             x = 0, ..., M.rows-1, y = 0, ..., M.columns-1
	// Version: 2.0
	//

	debug_handler_l(DMESSAGE, "add(const ring_matrix< T > &, const T &)",
			DVALUE);

	switch (this->bitfield.get_representation()) {
	case matrix_flags::dense_representation:
		if (this->rows != M.rows)
			this->D_ring_modul.set_no_of_rows(*this, M.rows);
		if (this->columns != M.columns)
			this->D_ring_modul.set_no_of_columns(*this, M.columns);

		switch (M.bitfield.get_representation()) {
		case matrix_flags::dense_representation:
			this->DDD_ring_modul.add(*this, M, a);
			break;
		case matrix_flags::sparse_representation:
			this->DSD_ring_modul.add(*this, M, a);
			break;
		}
		break;

	case matrix_flags::sparse_representation:
		if (this->rows != M.rows)
			this->S_ring_modul.set_no_of_rows(*this, M.rows);
		if (this->columns != M.columns)
			this->S_ring_modul.set_no_of_columns(*this, M.columns);

		switch (M.bitfield.get_representation()) {
		case matrix_flags::dense_representation:
			this->SDD_ring_modul.add(*this, M, a);
			break;
		case matrix_flags::sparse_representation:
			this->SSD_ring_modul.add(*this, M, a);
			break;
		}
		break;
	}
}



template< class T >
void ring_matrix< T >::
add(const T &a, const ring_matrix< T > &M)
{
	//
	//    Task: RES.add(a, M);
	// =  > RES.value[x][y] = a + M.value[x][y],
	//             x = 0, ..., M.rows-1, y = 0, ..., M.columns-1
	// Version: 2.0
	//

	debug_handler_l(DMESSAGE, "add(const T &, const ring_matrix< T > &)",
			DVALUE);

	switch (this->bitfield.get_representation()) {
	case matrix_flags::dense_representation:
		if (this->rows != M.rows)
			this->D_ring_modul.set_no_of_rows(*this, M.rows);
		if (this->columns != M.columns)
			this->D_ring_modul.set_no_of_columns(*this, M.columns);

		switch (M.bitfield.get_representation()) {
		case matrix_flags::dense_representation:
			this->DDD_ring_modul.add(*this, a, M);
			break;
		case matrix_flags::sparse_representation:
			this->DSD_ring_modul.add(*this, a, M);
			break;
		}
		break;

	case matrix_flags::sparse_representation:
		if (this->rows != M.rows)
			this->S_ring_modul.set_no_of_rows(*this, M.rows);
		if (this->columns != M.columns)
			this->S_ring_modul.set_no_of_columns(*this, M.columns);

		switch (M.bitfield.get_representation()) {
		case matrix_flags::dense_representation:
			this->SDD_ring_modul.add(*this, a, M);
			break;
		case matrix_flags::sparse_representation:
			this->SSD_ring_modul.add(*this, a, M);
			break;
		}
		break;
	}
}



//
// subtraction
//

template< class T >
void ring_matrix< T >::
subtract(const ring_matrix< T > &M, const ring_matrix< T > &N)
{
	//
	//       Task: RES.subtract(M, N)
	// =  > RES = M - N
	// Conditions: M.rows == N.rows and
	//             M.columns == N.columns
	//    Version: 2.0
	//

	debug_handler_l(DMESSAGE, "subtract(const ring_matrix< T > &, "
			"const ring_matrix< T > &)", DVALUE + 1);

	if (M.rows != N.rows || M.columns != N.columns)
		precondition_error_handler(M.rows, "M.rows", "M.rows == N.rows",
				    N.rows, "N.rows", "M.rows == N.rows",
				    M.columns, "M.columns", "M.columns == N.columns",
				    N.columns, "N.columns", "M.columns == N.columns",
				    "void ring_matrix< T >::"
				    "subtract(const ring_matrix< T > &M, "
				    "const ring_matrix< T > &N)",
				    DMESSAGE, EMESSAGE[4]);

	unsigned long SW = compute_switch(this->bitfield.get_representation(), M.bitfield.get_representation(),
					  N.bitfield.get_representation(), 0, 0);

	const unsigned long DR = matrix_flags::dense_representation;
	const unsigned long SR = matrix_flags::sparse_representation;

	switch (this->bitfield.get_representation()) {
	case matrix_flags::dense_representation:
		if (this->rows != N.rows)
			this->D_ring_modul.set_no_of_rows(*this, N.rows);
		if (this->columns != N.columns)
			this->D_ring_modul.set_no_of_columns(*this, N.columns);

		switch(SW) {
		case (compute_switch(DR, DR, DR, 0, 0)):
			this->DDD_ring_modul.subtract(*this, M, N);
			break;
		case (compute_switch(DR, DR, SR, 0, 0)):
			this->DDS_ring_modul.subtract(*this, M, N);
			break;
		case (compute_switch(DR, SR, DR, 0, 0)):
			this->DSD_ring_modul.subtract(*this, M, N);
			break;
		case (compute_switch(DR, SR, SR, 0, 0)):
			this->DSS_ring_modul.subtract(*this, M, N);
			break;
		}
		break;
	case matrix_flags::sparse_representation:
		if (this->rows != N.rows)
			this->S_ring_modul.set_no_of_rows(*this, N.rows);
		if (this->columns != N.columns)
			this->S_ring_modul.set_no_of_columns(*this, N.columns);

		switch(SW) {
		case (compute_switch(SR, DR, DR, 0, 0)):
			this->SDD_ring_modul.subtract(*this, M, N);
			break;
		case (compute_switch(SR, DR, SR, 0, 0)):
			this->SDS_ring_modul.subtract(*this, M, N);
			break;
		case (compute_switch(SR, SR, DR, 0, 0)):
			this->SSD_ring_modul.subtract(*this, M, N);
			break;
		case (compute_switch(SR, SR, SR, 0, 0)):
			this->S_ring_modul.subtract(*this, M, N);
			break;
		}
		break;
	}
}



template< class T >
void ring_matrix< T >::
subtract(const ring_matrix< T > &M, const T &a)
{
	//
	//    Task: RES.subtract(M, a)
	// =  > RES.value[x][y] = M.value[x][y] - a,
	//             x = 0, ..., M.rows-1, y = 0, ..., M.columns-1
	// Version: 2.0
	//

	debug_handler_l(DMESSAGE, "subtract(const ring_matrix< T > &, const T &)",
			DVALUE + 1);

	switch (this->bitfield.get_representation()) {
	case matrix_flags::dense_representation:
		if (this->rows != M.rows)
			this->D_ring_modul.set_no_of_rows(*this, M.rows);
		if (this->columns != M.columns)
			this->D_ring_modul.set_no_of_columns(*this, M.columns);

		switch (M.bitfield.get_representation()) {
		case matrix_flags::dense_representation:
			this->DDD_ring_modul.subtract(*this, M, a);
			break;
		case matrix_flags::sparse_representation:
			this->DSD_ring_modul.subtract(*this, M, a);
			break;
		}
		break;

	case matrix_flags::sparse_representation:
		if (this->rows != M.rows)
			this->S_ring_modul.set_no_of_rows(*this, M.rows);
		if (this->columns != M.columns)
			this->S_ring_modul.set_no_of_columns(*this, M.columns);

		switch (M.bitfield.get_representation()) {
		case matrix_flags::dense_representation:
			this->SDD_ring_modul.subtract(*this, M, a);
			break;
		case matrix_flags::sparse_representation:
			this->SSD_ring_modul.subtract(*this, M, a);
			break;
		}
		break;
	}
}



template< class T >
void ring_matrix< T >::
subtract(const T &a, const ring_matrix< T > &M)
{
	//
	//    Task: RES.subtract(a, M)
	// =  > RES.value[x][y] = a - M.value[x][y],
	//             x = 0, ..., M.rows-1, y = 0, ..., M.columns-1
	// Version: 2.0
	//

	debug_handler_l(DMESSAGE, "subtract(const T &, const ring_matrix< T > &)",
			DVALUE + 1);

	switch (this->bitfield.get_representation()) {
	case matrix_flags::dense_representation:
		if (this->rows != M.rows)
			this->D_ring_modul.set_no_of_rows(*this, M.rows);
		if (this->columns != M.columns)
			this->D_ring_modul.set_no_of_columns(*this, M.columns);

		switch (M.bitfield.get_representation()) {
		case matrix_flags::dense_representation:
			this->DDD_ring_modul.subtract(*this, a, M);
			break;
		case matrix_flags::sparse_representation:
			this->DSD_ring_modul.subtract(*this, a, M);
			break;
		}
		break;

	case matrix_flags::sparse_representation:
		if (this->rows != M.rows)
			this->S_ring_modul.set_no_of_rows(*this, M.rows);
		if (this->columns != M.columns)
			this->S_ring_modul.set_no_of_columns(*this, M.columns);

		switch (M.bitfield.get_representation()) {
		case matrix_flags::dense_representation:
			this->SDD_ring_modul.subtract(*this, a, M);
			break;
		case matrix_flags::sparse_representation:
			this->SSD_ring_modul.subtract(*this, a, M);
			break;
		}
		break;
	}
}



//
// multiplication
//

template< class T >
void ring_matrix< T >::
multiply(const ring_matrix< T > &M, const ring_matrix< T > &N)
{
	//
	//       Task: RES.multiply(A, B)
	// =  > RES = A * B
	// Conditions: A.columns == B.rows
	//    Version: 2.0
	//

	debug_handler_l(DMESSAGE, "multiply(const ring_matrix< T > &, "
			"const ring_matrix< T > &)", DVALUE + 2);

	if (M.columns != N.rows)
		precondition_error_handler(M.columns, "A.columns", "A.columns == B.rows",
				    N.rows, "B.rows", "A.columns == B.rows",
				    "void ring_matrix< T >::"
				    "multiply(const ring_matrix< T > &A, "
				    "const ring_matrix< T > &B)",
				    DMESSAGE, EMESSAGE[4]);

	unsigned long SW = compute_switch(this->bitfield.get_representation(), M.bitfield.get_representation(),
					  N.bitfield.get_representation(), 0, 0);

	const unsigned long DR = matrix_flags::dense_representation;
	const unsigned long SR = matrix_flags::sparse_representation;

	if (M.value != this->value && N.value != this->value) {
		switch (this->bitfield.get_representation()) {
		case matrix_flags::dense_representation:
			if (this->rows != M.rows)
				this->D_ring_modul.set_no_of_rows(*this, M.rows);
			if (this->columns != N.columns)
				this->D_ring_modul.set_no_of_columns(*this, N.columns);
			switch(SW) {
			case (compute_switch(DR, DR, DR, 0, 0)):
				this->DDD_ring_modul.multiply(*this, M, N);
				break;
			case (compute_switch(DR, DR, SR, 0, 0)):
				this->DDS_ring_modul.multiply(*this, M, N);
				break;
			case (compute_switch(DR, SR, DR, 0, 0)):
				this->DSD_ring_modul.multiply(*this, M, N);
				break;
			case (compute_switch(DR, SR, SR, 0, 0)):
				this->DSS_ring_modul.multiply(*this, M, N);
				break;
			}
			break;
		case matrix_flags::sparse_representation:
			if (this->rows != M.rows)
				this->S_ring_modul.set_no_of_rows(*this, M.rows);
			if (this->columns != N.columns)
				this->S_ring_modul.set_no_of_columns(*this, N.columns);
			switch(SW) {
			case (compute_switch(SR, DR, DR, 0, 0)):
				this->SDD_ring_modul.multiply(*this, M, N);
				break;
			case (compute_switch(SR, DR, SR, 0, 0)):
				this->SDS_ring_modul.multiply(*this, M, N);
				break;
			case (compute_switch(SR, SR, DR, 0, 0)):
				this->SSD_ring_modul.multiply(*this, M, N);
				break;
			case (compute_switch(SR, SR, SR, 0, 0)):
				this->SSS_ring_modul.multiply(*this, M, N);
				break;
			}
			break;
		}
	}
	else {
		switch (this->bitfield.get_representation()) {
		case matrix_flags::dense_representation:
		{
			ring_matrix< T > RES1(M.rows, N.columns, this->bitfield);
			switch(SW) {
			case (compute_switch(DR, DR, DR, 0, 0)):
				this->DDD_ring_modul.multiply(RES1, M, N);
				break;
			case (compute_switch(DR, DR, SR, 0, 0)):
				this->DDS_ring_modul.multiply(RES1, M, N);
				break;
			case (compute_switch(DR, SR, DR, 0, 0)):
				this->DSD_ring_modul.multiply(RES1, M, N);
				break;
			case (compute_switch(DR, SR, SR, 0, 0)):
				this->DSS_ring_modul.multiply(RES1, M, N);
				break;
			}
			if (this->rows != M.rows)
				this->D_ring_modul.set_no_of_rows(*this, M.rows);
			if (this->columns != N.columns)
				this->D_ring_modul.set_no_of_columns(*this, N.columns);
			this->D_ring_modul.assign(*this, RES1);
		}
		break;
		case matrix_flags::sparse_representation:
		{
			ring_matrix< T > RES1(M.rows, N.columns, this->bitfield);
			this->S_ring_modul.constructor(RES1, M.rows, N.columns);
			switch(SW) {
			case (compute_switch(SR, DR, DR, 0, 0)):
				this->SDD_ring_modul.multiply(RES1, M, N);
				break;
			case (compute_switch(SR, DR, SR, 0, 0)):
				this->SDS_ring_modul.multiply(RES1, M, N);
				break;
			case (compute_switch(SR, SR, DR, 0, 0)):
				this->SSD_ring_modul.multiply(RES1, M, N);
				break;
			case (compute_switch(SR, SR, SR, 0, 0)):
				this->SSS_ring_modul.multiply(RES1, M, N);
				break;
			}
			if (this->rows != M.rows)
				this->S_ring_modul.set_no_of_rows(*this, M.rows);
			if (this->columns != N.columns)
				this->S_ring_modul.set_no_of_columns(*this, N.columns);
			this->S_ring_modul.assign(*this, RES1);
		}
		break;
		}
	}
}



template< class T >
void ring_matrix< T >::
multiply(const ring_matrix< T > &M, const T &a)
{
	//
	//    Task: RES.multiply(A, k)
	// =  > RES.value[x][y] = A.value[x][y]*k,
	//             x = 0, ..., A.rows-1, y = 0, ..., A.columns-1
	// Version: 2.0
	//

	debug_handler_l(DMESSAGE,
			"multiply(const ring_matrix< T > &, const T &)",
			DVALUE + 2);

	switch (this->bitfield.get_representation()) {
	case matrix_flags::dense_representation:
		if (this->rows != M.rows)
			this->D_ring_modul.set_no_of_rows(*this, M.rows);
		if (this->columns != M.columns)
			this->D_ring_modul.set_no_of_columns(*this, M.columns);

		switch (M.bitfield.get_representation()) {
		case matrix_flags::dense_representation:
			this->DDD_ring_modul.multiply(*this, M, a);
			break;
		case matrix_flags::sparse_representation:
			this->DSD_ring_modul.multiply(*this, M, a);
			break;
		}
		break;

	case matrix_flags::sparse_representation:
		if (this->rows != M.rows)
			this->S_ring_modul.set_no_of_rows(*this, M.rows);
		if (this->columns != M.columns)
			this->S_ring_modul.set_no_of_columns(*this, M.columns);

		switch (M.bitfield.get_representation()) {
		case matrix_flags::dense_representation:
			this->SDD_ring_modul.multiply(*this, M, a);
			break;
		case matrix_flags::sparse_representation:
			this->SSD_ring_modul.multiply(*this, M, a);
			break;
		}
		break;
	}
}



template< class T >
void ring_matrix< T >::
multiply(const T &a, const ring_matrix< T > &M)
{
	//
	//    Task: RES.multiply(k, A)
	// =  > RES.value[x][y] = k*A.value[x][y],
	//             x = 0, ..., A.rows-1, y = 0, ..., A.columns-1
	// Version: 2.0
	//

	debug_handler_l(DMESSAGE, "multiply(const T &, const ring_matrix< T > &)",
			DVALUE + 2);

	switch (this->bitfield.get_representation()) {
	case matrix_flags::dense_representation:
		if (this->rows != M.rows)
			this->D_ring_modul.set_no_of_rows(*this, M.rows);
		if (this->columns != M.columns)
			this->D_ring_modul.set_no_of_columns(*this, M.columns);

		switch (M.bitfield.get_representation()) {
		case matrix_flags::dense_representation:
			this->DDD_ring_modul.multiply(*this, a, M);
			break;
		case matrix_flags::sparse_representation:
			this->DSD_ring_modul.multiply(*this, a, M);
			break;
		}
		break;

	case matrix_flags::sparse_representation:
		if (this->rows != M.rows)
			this->S_ring_modul.set_no_of_rows(*this, M.rows);
		if (this->columns != M.columns)
			this->S_ring_modul.set_no_of_columns(*this, M.columns);

		switch (M.bitfield.get_representation()) {
		case matrix_flags::dense_representation:
			this->SDD_ring_modul.multiply(*this, a, M);
			break;
		case matrix_flags::sparse_representation:
			this->SSD_ring_modul.multiply(*this, a, M);
			break;
		}
		break;
	}
}



template< class T >
void ring_matrix< T >::
compwise_multiply(const ring_matrix< T > &A, const ring_matrix< T > &B)
{
	//
	//    Task: RES.compwise_multiply(A, B)
	// =  > RES.value[x][y] = A.value[x][y]*B.value[x][y],
	//             x = 0, ..., A.rows-1, y = 0, ..., A.columns-1
	// Version: 2.0
	//

	debug_handler_l(DMESSAGE, "compwise_multiply(const ring_matrix< T > &, "
			"const ring_matrix< T > &)", DVALUE + 2);

	if (A.rows != B.rows || A.columns != B.columns)
		precondition_error_handler(A.rows, "A.rows", "A.rows == B.rows",
				    B.rows, "B.rows", "A.rows == B.rows",
				    A.columns, "A.columns", "A.columns == B.columns",
				    B.columns, "B.columns", "B.columns == B.columns",
				    "void ring_matrix< T >::"
				    "compwise_multiply(const ring_matrix< T > &A, "
				    "const ring_matrix< T > &B)",
				    DMESSAGE, EMESSAGE[4]);

	unsigned long SW = compute_switch(this->bitfield.get_representation(), A.bitfield.get_representation(),
					  B.bitfield.get_representation(), 0, 0);

	const unsigned long DR = matrix_flags::dense_representation;
	const unsigned long SR = matrix_flags::sparse_representation;

	switch (this->bitfield.get_representation()) {
	case matrix_flags::dense_representation:
		if (this->rows != A.rows)
			this->D_ring_modul.set_no_of_rows(*this, A.rows);
		if (this->columns != A.columns)
			this->D_ring_modul.set_no_of_columns(*this, B.columns);

		switch(SW) {
		case (compute_switch(DR, DR, DR, 0, 0)):
			this->DDD_ring_modul.compwise_multiply(*this, A, B);
			break;
		case (compute_switch(DR, DR, SR, 0, 0)):
			this->DDS_ring_modul.compwise_multiply(*this, A, B);
			break;
		case (compute_switch(DR, SR, DR, 0, 0)):
			this->DSD_ring_modul.compwise_multiply(*this, A, B);
			break;
		case (compute_switch(DR, SR, SR, 0, 0)):
			this->DSS_ring_modul.compwise_multiply(*this, A, B);
			break;
		}
		break;
	case matrix_flags::sparse_representation:
		if (this->rows != A.rows)
			this->S_ring_modul.set_no_of_rows(*this, A.rows);
		if (this->columns != A.columns)
			this->S_ring_modul.set_no_of_columns(*this, A.columns);

		switch(SW) {
		case (compute_switch(SR, DR, DR, 0, 0)):
			this->SDD_ring_modul.compwise_multiply(*this, A, B);
			break;
		case (compute_switch(SR, DR, SR, 0, 0)):
			this->SDS_ring_modul.compwise_multiply(*this, A, B);
			break;
		case (compute_switch(SR, SR, DR, 0, 0)):
			this->SSD_ring_modul.compwise_multiply(*this, A, B);
			break;
		case (compute_switch(SR, SR, SR, 0, 0)):
			this->SSS_ring_modul.compwise_multiply(*this, A, B);
			break;
		}
		break;
	}
}



template< class T >
void ring_matrix< T >::
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
		precondition_error_handler(this->columns, "this->columns", "columns == v.size",
				    v.size(), "v.size", "columns == v.size",
				    "void ring_matrix< T >::"
				    "multiply_right(math_vector< T > &res, "
				    "const math_vector< T > &v) const",
				    DMESSAGE, EMESSAGE[1]);

	switch (this->bitfield.get_representation()) {
	case matrix_flags::dense_representation:
	{
		lidia_size_t l = res.size();
		T *tmp1 = res.get_data_address(), *tmp2 = v.get_data_address();
		if (tmp1 == tmp2) {
			math_vector< T > res2(this->rows, this->rows);
			tmp1 = res2.get_data_address();
			this->D_ring_modul.multiply_right(*this, tmp1, tmp2);
			res = res2;
		}
		else {
			if (l != this->rows) {
				res.set_capacity(this->rows);
				res.set_size(this->rows);
			}
			tmp1 = res.get_data_address();
			this->D_ring_modul.multiply_right(*this, tmp1, tmp2);
		}
		break;
	}
	case matrix_flags::sparse_representation:
		this->SSS_ring_modul.multiply_right(*this, res, v);
		break;
	}

}



template< class T >
void ring_matrix< T >::
multiply_right(T *&res, const T *v) const
{
	//
	//       Task: c[x] = A.value[x][0]*v[0]+...
	//                           +A.value[x][A.columns-1]*v[A.columns-1],
	//             x = 0, ..., A.rows-1
	// Conditions: v != NULL
	//    Version: 2.0
	//

	debug_handler_l(DMESSAGE, "multiply_right(T *&, const T *)",
			DVALUE + 2);

	if (v == NULL)
		precondition_error_handler(PRT, "v", "v != NULL",
				    "void ring_matrix< T >::"
				    "multiply_right(T *&c, const T *v) const",
				    DMESSAGE, EMESSAGE[6]);

	switch (this->bitfield.get_representation()) {
	case matrix_flags::dense_representation:
		this->D_ring_modul.multiply_right(*this, res, v);
		break;
	case matrix_flags::sparse_representation:
		this->SSS_ring_modul.multiply_right(*this, res, v);
		break;
	}
}



template< class T >
void ring_matrix< T >::
multiply_left(math_vector< T > &res, const math_vector< T > &v) const
{
	//
	//       Task: c[x] = v[0]*A.value[0][x]+...
	//                         +v[A.rows-1]*A.value[A.rows-1][x],
	//             x = 0, ..., A.columns-1
	// Conditions: v.length == A.rows
	//    Version: 2.0
	//

	debug_handler_l(DMESSAGE,
			"multiply_left(math_vector< T > &, const math_vector< T > &)",
			DVALUE + 2);

	if (v.size() != this->rows)
		precondition_error_handler(this->rows, "rows", "rows == v.size",
				    v.size(), "v.size", "rows == v.size",
				    "void ring_matrix< T >::"
				    "multiply_left(math_vector< T > &res, "
				    "const math_vector< T > &v) const",
				    DMESSAGE, EMESSAGE[6]);

	switch (this->bitfield.get_representation()) {
	case matrix_flags::dense_representation:
	{
		lidia_size_t l = res.size();
		T *tmp1 = res.get_data_address(), *tmp2 = v.get_data_address();
		if (tmp1 == tmp2) {
			math_vector< T > res2(this->columns, this->columns);
			tmp1 = res2.get_data_address();

			this->D_ring_modul.multiply_left(*this, tmp1, tmp2);
			res = res2;
		}
		else {
			if (l != this->columns) {
				if (res.capacity() < this->columns)
					res.set_capacity(this->columns);
				if (l != this->columns)
					res.set_size(this->columns);
			}
			tmp1 = res.get_data_address();
			this->D_ring_modul.multiply_left(*this, tmp1, tmp2);
		}
		break;
	}
	case matrix_flags::sparse_representation:
		this->SSS_ring_modul.multiply_left(*this, res, v);
		break;
	}
}



template< class T >
void ring_matrix< T >::
multiply_left(T *&res, const T *v) const
{
	//
	//       Task: c[x] = v[0]*A.value[0][x]+...
	//                           +v[A.rows-1]*A.value[A.rows-1][x],
	//             x = 0, ..., A.columns-1
	// Conditions: v != NULL
	//    Version: 2.0
	//

	debug_handler_l(DMESSAGE, "multiply_left(T *&, const T *)",
			DVALUE + 2);

	if (v == NULL)
		precondition_error_handler(PRT, "v", "v != NULL",
				    "void ring_matrix< T >::"
				    "multiply_left(T *&c, const T *v) const",
				    DMESSAGE, EMESSAGE[6]);

	switch (this->bitfield.get_representation()) {
	case matrix_flags::dense_representation:
		this->D_ring_modul.multiply_left(*this, res, v);
		break;
	case matrix_flags::sparse_representation:
		this->SSS_ring_modul.multiply_left(*this, res, v);
		break;
	}
}



//
// negation
//

template< class T >
void ring_matrix< T >::
negate(const ring_matrix< T > &M)
{
	//
	//    Task: A.negate(B)
	// =  > A.value[x][y] = - B.value[x][y],
	//             x = 0, ..., rows-1, y = 0, ..., columns-1
	// Version: 2.0
	//

	debug_handler_l(DMESSAGE,
			"negate(const ring_matrix< T > &)",
			DVALUE + 4);

	switch (this->bitfield.get_representation()) {
	case matrix_flags::dense_representation:
		switch (M.bitfield.get_representation()) {
		case matrix_flags::dense_representation:
			this->DDD_ring_modul.negate(*this, M);
			break;
		case matrix_flags::sparse_representation:
			this->DSD_ring_modul.negate(*this, M);
			break;
		}
		break;

	case matrix_flags::sparse_representation:
		switch (M.bitfield.get_representation()) {
		case matrix_flags::dense_representation:
			this->SDD_ring_modul.negate(*this, M);
			break;
		case matrix_flags::sparse_representation:
			this->SSD_ring_modul.negate(*this, M);
			break;
		}
		break;
	}
}



//////////////////////////////////
// END: arithmetical procedures //
//////////////////////////////////

//
// comparisons
//

template< class T >
bool ring_matrix< T >::
equal(const ring_matrix< T > &M) const
{
	//
	//    Task: A.equal(B) returns true, if A.value[x][y] == B.value[x][y],
	//          x = 0, ..., A.rows-1, y = 0, ..., A.columns-1
	//          and A.columns == B.columns and A.rows == B.rows
	// Version: 2.0
	//

	debug_handler_l(DMESSAGE, "equal(const ring_matrix< T > &)",
			DVALUE + 5);

	switch (this->bitfield.get_representation()) {
	case matrix_flags::dense_representation:
		switch (M.bitfield.get_representation()) {

		case matrix_flags::dense_representation:
			return this->DDD_ring_modul.equal(*this, M);
		case matrix_flags::sparse_representation:
			return this->DSD_ring_modul.equal(*this, M);
		}
		break;

	case matrix_flags::sparse_representation:
		switch (M.bitfield.get_representation()) {
		case matrix_flags::dense_representation:
			return this->SDD_ring_modul.equal(*this, M);
		case matrix_flags::sparse_representation:
			return this->SSD_ring_modul.equal(*this, M);
		}
		break;
	}
	return false;
}



//
// trace
//

template< class T >
void ring_matrix< T >::
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
				    "void ring_matrix< T >::"
				    "trace(T &tr) const",
				    DMESSAGE, EMESSAGE[7]);

	switch (this->bitfield.get_representation()) {
	case matrix_flags::dense_representation:
		this->D_ring_modul.trace(*this, tr);
		break;
	case matrix_flags::sparse_representation:
		this->SSS_ring_modul.trace(*this, tr);
		break;
	}
}



//
// size reduction
//

template< class T >
void ring_matrix< T >::
size_reduction()
{
	//
	//    Task: A.size_reduction()
	//          reduces the no of rows according to the structured gauss elimination
	// Version: 2.0
	//

	debug_handler_l(DMESSAGE, "size_reduction()", DVALUE + 6);

	switch (this->bitfield.get_representation()) {
	case matrix_flags::dense_representation:
		this->DDD_ring_modul.size_reduction(*this);
		break;
	case matrix_flags::sparse_representation:
		this->SSS_ring_modul.size_reduction(*this);
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



#endif	// LIDIA_RING_MATRIX_CC_GUARD_
