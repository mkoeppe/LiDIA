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


#ifdef HAVE_CONFIG_H
# include	"config.h"
#endif
#include	"LiDIA/base_matrix.h"
#include	"LiDIA/modular_operations.inl"
#include "LiDIA/precondition_error.h"


#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



//
// debug defines / error defines
//

extern const char *matrix_error_msg[];
extern const char *PRT;

#define DVALUE   LDBL_MATRIX           // Debug value
#define DMESSAGE "base_matrix"         // Debug message
#define EMESSAGE matrix_error_msg      // Error message array

//
// debug level
//
//   0 : constructors / destructor
//   1 : input / output
//   2 : access functions
//   3 : split / compose
//   4 : swap functions
//   5 : structur functions
//   6 : assignments
//   7 : diagonal and transpose functions
//   8 : stream handling
//   9 : internal functions
//  10 : boolean functions
//

//
// constructors
//

base_matrix< bigmod >::
base_matrix()
{
	debug_handler_l(DMESSAGE, "base_matrix()", DVALUE);

	DD_base_modul_bigint.constructor(bigint_value);
	DD_base_modul_long.constructor(long_value);

	debug_handler_c(DMESSAGE, "base_matrix()", DVALUE,
			status_report());
}



base_matrix< bigmod >::
base_matrix(const bigint &mod)
{
	debug_handler_l(DMESSAGE, "base_matrix(const bigint &)", DVALUE);

	if (mod.bit_length() > bigint::bits_per_digit()) {
		BigintLong = true;
		bigint_modulus = mod;
		DD_base_modul_bigint.constructor(bigint_value);
	}
	else {
		BigintLong = false;
		mod.longify(long_modulus);
		DD_base_modul_long.constructor(long_value);
	}

	debug_handler_c(DMESSAGE, "base_matrix(const bigint &)", DVALUE,
			status_report());
}



base_matrix< bigmod >::
base_matrix(const bigint &mod, const matrix_flags &flag)
{
	debug_handler_l(DMESSAGE, "base_matrix(const bigint &, "
			"const matrix_flags &)", DVALUE);

	if (mod.bit_length() > bigint::bits_per_digit()) {
		BigintLong = true;
		bigint_modulus = mod;
		bigint_value.bitfield = flag;
		if (bigint_value.bitfield.get_representation() ==
		    matrix_flags::dense_representation)
			DD_base_modul_bigint.constructor(bigint_value);
		else
			SS_base_modul_bigint.constructor(bigint_value);
	}
	else {
		BigintLong = false;
		mod.longify(long_modulus);
		long_value.bitfield = flag;
		if (bigint_value.bitfield.get_representation() ==
		    matrix_flags::dense_representation)
			DD_base_modul_long.constructor(long_value);
		else
			SS_base_modul_long.constructor(long_value);
	}

	debug_handler_c(DMESSAGE, "base_matrix(const bigint &, "
			"const matrix_flags &)", DVALUE, status_report());
}



base_matrix< bigmod >::
base_matrix(lidia_size_t r, lidia_size_t c, const bigint &mod)
{
	debug_handler_l(DMESSAGE, "base_matrix(lidia_size_t, lidia_size_t, "
			"const bigint &)", DVALUE);

	if (r < 0 || c < 0)
		precondition_error_handler(r, "r", "r >= 0",
				    c, "c", "c >= 0",
				    "base_matrix< bigmod >::"
				    "base_matrix(lidia_size_t r, lidia_size_t c, "
				    "const bigint &)", DMESSAGE, EMESSAGE[0]);

	if (mod.bit_length() > bigint::bits_per_digit()) {
		BigintLong = true;
		bigint_modulus = mod;
		D_base_modul_bigint.constructor(bigint_value, r, c);
	}
	else {
		BigintLong = false;
		mod.longify(long_modulus);
		D_base_modul_long.constructor(long_value, r, c);
	}

	debug_handler_c(DMESSAGE, "base_matrix(lidia_size_t, lidia_size_t, "
			"const bigint &)", DVALUE, status_report());
}



base_matrix< bigmod >::
base_matrix(lidia_size_t r, lidia_size_t c, const bigint &mod,
	    const matrix_flags &flag)
{
	debug_handler_l(DMESSAGE, "base_matrix(lidia_size_t, lidia_size_t, "
			"const bigint &, const matrix_flags &)", DVALUE);

	if (r < 0 || c < 0)
		precondition_error_handler(r, "r", "r >= 0",
				    c, "c", "c >= 0",
				    "base_matrix< bigmod >::"
				    "base_matrix(lidia_size_t r, lidia_size_t c, "
				    "const bigint &, const matrix_flags &flag)",
				    DMESSAGE, EMESSAGE[0]);

	if (mod.bit_length() > bigint::bits_per_digit()) {
		BigintLong = true;
		bigint_modulus = mod;
		bigint_value.bitfield = flag;
		if (bigint_value.bitfield.get_representation() ==
		    matrix_flags::dense_representation)
			D_base_modul_bigint.constructor(bigint_value, r, c);
		else
			S_base_modul_bigint.constructor(bigint_value, r, c);
	}
	else {
		BigintLong = false;
		mod.longify(long_modulus);
		long_value.bitfield = flag;
		if (long_value.bitfield.get_representation() ==
		    matrix_flags::dense_representation)
			D_base_modul_long.constructor(long_value, r, c);
		else
			S_base_modul_long.constructor(long_value, r, c);
	}

	debug_handler_c(DMESSAGE, "base_matrix(lidia_size_t, lidia_size_t, "
			"const bigint &, const matrix_flags &)",
			DVALUE, status_report());
}



base_matrix< bigmod >::
base_matrix(const base_matrix< bigmod > &M)
{
	debug_handler_l(DMESSAGE, "base_matrix(const base_matrix< bigmod > &)",
			DVALUE);

	BigintLong = M.BigintLong;

	if (BigintLong) {
		bigint_modulus = M.bigint_modulus;
		bigint_value.bitfield = M.bigint_value.bitfield;

		switch (M.bigint_value.bitfield.get_representation()) {
		case matrix_flags::dense_representation:
			D_base_modul_bigint.constructor(bigint_value, M.bigint_value);
			break;
		case matrix_flags::sparse_representation:
			S_base_modul_bigint.constructor(bigint_value, M.bigint_value);
			break;
		}
	}
	else {
		long_modulus = M.long_modulus;
		long_value.bitfield = M.long_value.bitfield;

		switch (M.long_value.bitfield.get_representation()) {
		case matrix_flags::dense_representation:
			D_base_modul_long.constructor(long_value, M.long_value);
			break;
		case matrix_flags::sparse_representation:
			S_base_modul_long.constructor(long_value, M.long_value);
			break;
		}
	}
	debug_handler_c(DMESSAGE, "base_matrix(const base_matrix< bigmod > &)",
			DVALUE, status_report());
}



base_matrix< bigmod >::
base_matrix(const base_matrix< bigmod > &M, const matrix_flags &flag)
{
	debug_handler_l(DMESSAGE, "base_matrix(const base_matrix< bigmod > &, "
			"const matrix_flags &)", DVALUE);

	BigintLong = M.BigintLong;

	if (BigintLong) {
		bigint_modulus = M.bigint_modulus;
		bigint_value.bitfield = flag;

		switch (flag.get_representation()) {
		case matrix_flags::dense_representation:
			switch(M.bigint_value.bitfield.get_representation()) {
			case matrix_flags::dense_representation:
				DD_base_modul_bigint.constructor(bigint_value, M.bigint_value);
				break;
			case matrix_flags::sparse_representation:
				DS_base_modul_bigint.constructor(bigint_value, M.bigint_value);
				break;
			}
			break;
		case matrix_flags::sparse_representation:
			switch(M.bigint_value.bitfield.get_representation()) {
			case matrix_flags::dense_representation:
				SD_base_modul_bigint.constructor(bigint_value, M.bigint_value);
				break;
			case matrix_flags::sparse_representation:
				SS_base_modul_bigint.constructor(bigint_value, M.bigint_value);
				break;
			}
			break;
		}
	}
	else {
		long_modulus = M.long_modulus;
		long_value.bitfield = flag;

		switch (flag.get_representation()) {
		case matrix_flags::dense_representation:
			switch(M.long_value.bitfield.get_representation()) {
			case matrix_flags::dense_representation:
				DD_base_modul_long.constructor(long_value, M.long_value);
				break;
			case matrix_flags::sparse_representation:
				DS_base_modul_long.constructor(long_value, M.long_value);
				break;
			}
			break;
		case matrix_flags::sparse_representation:
			switch(M.long_value.bitfield.get_representation()) {
			case matrix_flags::dense_representation:
				SD_base_modul_long.constructor(long_value, M.long_value);
				break;
			case matrix_flags::sparse_representation:
				SS_base_modul_long.constructor(long_value, M.long_value);
				break;
			}
			break;
		}
	}

	debug_handler_c(DMESSAGE, "base_matrix(const base_matrix< bigmod > &, "
			"const matrix_flags &)", DVALUE, status_report());
}



//
// destructor
//

base_matrix< bigmod >::
~base_matrix()
{
	debug_handler_l(DMESSAGE, "~base_matrix()", DVALUE);

	if (BigintLong)
		switch (bigint_value.bitfield.get_representation()) {
		case matrix_flags::dense_representation:
			D_base_modul_bigint.destructor(bigint_value);
			break;
		case matrix_flags::sparse_representation:
			S_base_modul_bigint.destructor(bigint_value);
			break;
		}
	else
		switch (long_value.bitfield.get_representation()) {
		case matrix_flags::dense_representation:
			D_base_modul_long.destructor(long_value);
			break;
		case matrix_flags::sparse_representation:
			S_base_modul_long.destructor(long_value);
			break;
		}
}



//
// Input / Output
//

std::ostream & base_matrix< bigmod >::
write(std::ostream &out) const
{
	//
	//    Task: Output to stream
	// Version: 2.0
	//

	debug_handler_l(DMESSAGE, "print(std::ostream &)", DVALUE + 1);

	if (BigintLong)
		switch(bigint_value.bitfield.get_print_mode()) {
		case matrix_flags::beauty_mode:
			switch (bigint_value.bitfield.get_representation()) {
			case matrix_flags::dense_representation:
				D_base_modul_bigint.write_to_beauty(bigint_value, out);
				break;
			case matrix_flags::sparse_representation:
				S_base_modul_bigint.write_to_beauty(bigint_value, out);
				break;
			}
			break;
		case matrix_flags::lidia_mode:
			switch (bigint_value.bitfield.get_representation()) {
			case matrix_flags::dense_representation:
				D_base_modul_bigint.write_to_stream(bigint_value, out);
				break;
			case matrix_flags::sparse_representation:
				S_base_modul_bigint.write_to_stream(bigint_value, out);
				break;
			}
			break;
		case matrix_flags::gp_mode:
			switch (bigint_value.bitfield.get_representation()) {
			case matrix_flags::dense_representation:
				D_base_modul_bigint.write_to_gp(bigint_value, out);
				break;
			case matrix_flags::sparse_representation:
				S_base_modul_bigint.write_to_gp(bigint_value, out);
				break;
			}
			break;
		case matrix_flags::maple_mode:
			switch (bigint_value.bitfield.get_representation()) {
			case matrix_flags::dense_representation:
				D_base_modul_bigint.write_to_maple(bigint_value, out);
				break;
			case matrix_flags::sparse_representation:
				S_base_modul_bigint.write_to_maple(bigint_value, out);
				break;
			}
			break;
		case matrix_flags::mathematica_mode:
			switch (bigint_value.bitfield.get_representation()) {
			case matrix_flags::dense_representation:
				D_base_modul_bigint.write_to_mathematica(bigint_value, out);
				break;
			case matrix_flags::sparse_representation:
				S_base_modul_bigint.write_to_mathematica(bigint_value, out);
				break;
			}
			break;
		case matrix_flags::kash_mode:
			switch (bigint_value.bitfield.get_representation()) {
			case matrix_flags::dense_representation:
				D_base_modul_bigint.write_to_kash(bigint_value, out);
				break;
			case matrix_flags::sparse_representation:
				S_base_modul_bigint.write_to_kash(bigint_value, out);
				break;
			}
			break;
		case matrix_flags::latex_mode:
			switch (bigint_value.bitfield.get_representation()) {
			case matrix_flags::dense_representation:
				D_base_modul_bigint.write_to_latex(bigint_value, out);
				break;
			case matrix_flags::sparse_representation:
				S_base_modul_bigint.write_to_latex(bigint_value, out);
				break;
			}
			break;
		}
	else
		switch(long_value.bitfield.get_print_mode()) {
		case matrix_flags::beauty_mode:
			switch (long_value.bitfield.get_representation()) {
			case matrix_flags::dense_representation:
				D_base_modul_long.write_to_beauty(long_value, out);
				break;
			case matrix_flags::sparse_representation:
				S_base_modul_long.write_to_beauty(long_value, out);
				break;
			}
			break;
		case matrix_flags::lidia_mode:
			switch (long_value.bitfield.get_representation()) {
			case matrix_flags::dense_representation:
				D_base_modul_long.write_to_stream(long_value, out);
				break;
			case matrix_flags::sparse_representation:
				S_base_modul_long.write_to_stream(long_value, out);
				break;
			}
			break;
		case matrix_flags::gp_mode:
			switch (long_value.bitfield.get_representation()) {
			case matrix_flags::dense_representation:
				D_base_modul_long.write_to_gp(long_value, out);
				break;
			case matrix_flags::sparse_representation:
				S_base_modul_long.write_to_gp(long_value, out);
				break;
			}
			break;
		case matrix_flags::maple_mode:
			switch (long_value.bitfield.get_representation()) {
			case matrix_flags::dense_representation:
				D_base_modul_long.write_to_maple(long_value, out);
				break;
			case matrix_flags::sparse_representation:
				S_base_modul_long.write_to_maple(long_value, out);
				break;
			}
			break;
		case matrix_flags::mathematica_mode:
			switch (long_value.bitfield.get_representation()) {
			case matrix_flags::dense_representation:
				D_base_modul_long.write_to_mathematica(long_value, out);
				break;
			case matrix_flags::sparse_representation:
				S_base_modul_long.write_to_mathematica(long_value, out);
				break;
			}
			break;
		case matrix_flags::kash_mode:
			switch (long_value.bitfield.get_representation()) {
			case matrix_flags::dense_representation:
				D_base_modul_long.write_to_kash(long_value, out);
				break;
			case matrix_flags::sparse_representation:
				S_base_modul_long.write_to_kash(long_value, out);
				break;
			}
			break;
		case matrix_flags::latex_mode:
			switch (long_value.bitfield.get_representation()) {
			case matrix_flags::dense_representation:
				D_base_modul_long.write_to_latex(long_value, out);
				break;
			case matrix_flags::sparse_representation:
				S_base_modul_long.write_to_latex(long_value, out);
				break;
			}
			break;
		}
	return out;
}



std::istream & base_matrix< bigmod >::
read(std::istream &in)
{
	//
	//       Task: Input from stream
	// Conditions: first character in {[,{,1,2,3,4,5,6,7,8,9,M,a}
	//    Version: 2.0
	//

	debug_handler_l(DMESSAGE, "read(std::istream &)", DVALUE + 1);

	char c;
	in >> c;
	in.putback(c);

	switch(c) {
	case '1':
	case '2':
	case '3':
	case '4':
	case '5':
	case '6':
	case '7':
	case '8':
	case '9':
		switch (bigint_value.bitfield.get_representation()) {
		case matrix_flags::dense_representation:
			D_base_modul_bigint.read_from_stream(bigint_value, in);
			break;
		case matrix_flags::sparse_representation:
			S_base_modul_bigint.read_from_stream(bigint_value, in);
			break;
		}
		break;
	case '[':
		switch (bigint_value.bitfield.get_representation()) {
		case matrix_flags::dense_representation:
			D_base_modul_bigint.read_from_gp(bigint_value, in);
			break;
		case matrix_flags::sparse_representation:
			S_base_modul_bigint.read_from_gp(bigint_value, in);
			break;
		}
		break;
	case 'a':
		switch (bigint_value.bitfield.get_representation()) {
		case matrix_flags::dense_representation:
			D_base_modul_bigint.read_from_maple(bigint_value, in);
			break;
		case matrix_flags::sparse_representation:
			S_base_modul_bigint.read_from_maple(bigint_value, in);
			break;
		}
		break;
	case '{':
		switch (bigint_value.bitfield.get_representation()) {
		case matrix_flags::dense_representation:
			D_base_modul_bigint.read_from_mathematica(bigint_value, in);
			break;
		case matrix_flags::sparse_representation:
			S_base_modul_bigint.read_from_mathematica(bigint_value, in);
			break;
		}
		break;
	case 'M':
		switch (bigint_value.bitfield.get_representation()) {
		case matrix_flags::dense_representation:
			D_base_modul_bigint.read_from_kash(bigint_value, in);
			break;
		case matrix_flags::sparse_representation:
			S_base_modul_bigint.read_from_kash(bigint_value, in);
			break;
		}
		break;
	default:
		precondition_error_handler(c, "c", "Character c have to be in "
				    "{[, {, a, M, 1, 2, 3, 4, 5, 6, 7, 8, 9}",
				    "std::istream & base_matrix< bigmod >::"
				    "read(std::istream &in)",
				    DMESSAGE, EMESSAGE[2]);
	}
	return in;
}



/////////////////////////////
// BEGIN: access functions //
/////////////////////////////

//
// element access
//

void base_matrix< bigmod >::
sto(lidia_size_t x, lidia_size_t y, const bigmod &l)
{
	//
	//       Task: A.sto(x,y,e) stores e at position (x,y) in matrix A.
	// Conditions: 0 <= x < A.rows and
	//             0 <= y < A.columns
	//    Version: 2.0
	//

	debug_handler_l(DMESSAGE, "sto(lidia_size_t, lidia_size_t, const bigmod &)",
			DVALUE + 2);

	if (BigintLong) {
		if (x< 0 || x >= bigint_value.rows
		    || y< 0 || y >= bigint_value.columns
		    || bigint_modulus != l.modulus())
			precondition_error_handler(x, "x", "0 <= x < rows",
					    bigint_value.rows, "rows", "",
					    y, "y", "0 <= y < columns",
					    bigint_value.columns, "columns", "",
					    bigint_modulus, "modulus", "modulus == e.modulus",
					    l.modulus(), "e.modulus", "",
					    "void base_matrix< bigmod >::"
					    "sto(lidia_size_t x, lidia_size_t y, "
					    "const bigmod &e)", DMESSAGE, EMESSAGE[3]);

		bigint e_bigint = l.mantissa();
		LiDIA::best_remainder(e_bigint, e_bigint, bigint_modulus);

		switch(bigint_value.bitfield.get_representation()) {
		case matrix_flags::dense_representation:
			D_base_modul_bigint.sto(bigint_value, x, y, e_bigint);
			break;
		case matrix_flags::sparse_representation:
			S_base_modul_bigint.sto(bigint_value, x, y, e_bigint);
			break;
		}
	}
	else {
		if (x< 0 || x >= long_value.rows
		    || y< 0 || y >= long_value.columns
		    || bigint(long_modulus) != l.modulus())
			precondition_error_handler(x, "x", "0 <= x < rows",
					    long_value.rows, "rows", "",
					    y, "y", "0 <= y < columns",
					    long_value.columns, "columns", "",
					    long_modulus, "modulus", "modulus == e.modulus",
					    l.modulus(), "e.modulus", "",
					    "void base_matrix< bigmod >::"
					    "sto(lidia_size_t x, lidia_size_t y, "
					    "const bigmod &e)", DMESSAGE, EMESSAGE[3]);

		long e_long;
		(l.mantissa()).longify(e_long);
		LiDIA::best_remainder(e_long, e_long, long_modulus);

		switch(long_value.bitfield.get_representation()) {
		case matrix_flags::dense_representation:
			D_base_modul_long.sto(long_value, x, y, e_long);
			break;
		case matrix_flags::sparse_representation:
			S_base_modul_long.sto(long_value, x, y, e_long);
			break;
		}
	}
}



const bigmod base_matrix< bigmod >::
member(lidia_size_t x, lidia_size_t y) const
{
	//
	//       Task: A.member(x,y) = value at position (x,y) in matrix A
	// Conditions: 0 <= x < A.rows and
	//             0 <= y < A.columns
	//    Version: 2.0
	//

	debug_handler_l(DMESSAGE, "member(lidia_size_t, lidia_size_t)", DVALUE + 2);

	bigmod TMP;
	if (BigintLong) {
		if (x< 0 || y < 0 || x >= bigint_value.rows
		    || y >= bigint_value.columns)
			precondition_error_handler(x, "x", "0 <= x < rows",
					    bigint_value.rows, "rows", "",
					    y, "y", "0 <= y < columns",
					    bigint_value.columns, "columns", "",
					    "const bigmod & base_matrix< bigmod >::"
					    "member(lidia_size_t x, lidia_size_t y) const",
					    DMESSAGE, EMESSAGE[3]);

		TMP.set_modulus(bigint_modulus);
		TMP = (bigint_value.bitfield.get_representation() ==
		       matrix_flags::dense_representation)
			? D_base_modul_bigint.member(bigint_value, x, y)
			: S_base_modul_bigint.member(bigint_value, x, y);
	}
	else {
		if (x< 0 || y < 0 || x >= long_value.rows
		    || y >= long_value.columns)
			precondition_error_handler(x, "x", "0 <= x < rows",
					    long_value.rows, "rows", "",
					    y, "y", "0 <= y < columns",
					    long_value.columns, "columns", "",
					    "const bigmod & base_matrix< bigmod >::"
					    "member(lidia_size_t x, lidia_size_t y) const",
					    DMESSAGE, EMESSAGE[3]);

		TMP.set_modulus(long_modulus);
		TMP = (long_value.bitfield.get_representation() ==
		       matrix_flags::dense_representation)
			? D_base_modul_long.member(long_value, x, y)
			: S_base_modul_long.member(long_value, x, y);
	}
	return (const bigmod &)TMP;
}



//
// column access
//

void base_matrix< bigmod >::
sto_column(const bigmod *v, lidia_size_t l, lidia_size_t j, lidia_size_t from)
{
	//
	//       Task: A.sto_column(v,l,j,p) stores v[0],...,v[l-1] in column j
	//             of matrix A starting at position p.
	// Conditions: 0 <= j < A.columns and
	//             0 <= l <= A.rows and
	//             0 <= from < A.rows and
	//             from + l <= A.rows and
	//             v != NULL
	// Version: 2.0
	//

	debug_handler_l(DMESSAGE, "sto_column(const bigmod *, lidia_size_t, "
			"lidia_size_t, lidia_size_t)", DVALUE + 2);

	lidia_size_t i;
	bigmod TMP;
	if (BigintLong) {
		if (j >= bigint_value.columns || j < 0
		    || l > bigint_value.rows || l < 0
		    || from< 0 || from >= bigint_value.rows
		    || from + l > bigint_value.rows || v == NULL
		    || TMP.modulus() != bigint_modulus)
			precondition_error_handler(j, "j", "0 <= j < columns",
					    bigint_value.columns, "columns", "",
					    l, "l", "0 <= l <= rows and l+from <= rows",
					    from, "from",
					    "0 <= from < rows and l + from <= rows",
					    bigint_value.rows, "rows", "",
					    PRT, "v", " v != NULL",
					    bigint_modulus, "modulus",
					    "modulus == global modulus",
					    TMP.modulus(), "global modulus", "",
					    "void base_matrix< bigmod >::"
					    "sto_column(const bigmod *v, lidia_size_t l, "
					    "lidia_size_t j, lidia_size_t from)",
					    DMESSAGE, EMESSAGE[1]);

		bigint *tmp = new bigint[l];
		memory_handler(tmp, DM_BM, "sto_column :: "
			       "Error in memory allocation (tmp)");
		for (i = 0; i < l; i++)
			tmp[i] = v[i].mantissa();

		switch(bigint_value.bitfield.get_representation()) {
		case matrix_flags::dense_representation:
			D_base_modul_bigint.sto_column(bigint_value, tmp, l, j, from);
			break;
		case matrix_flags::sparse_representation:
			S_base_modul_bigint.sto_column(bigint_value, tmp, l, j, from);
			break;
		}
		delete[] tmp;
	}
	else {
		if (j >= long_value.columns || j < 0
		    || l > long_value.rows || l < 0
		    || from< 0 || from >= long_value.rows
		    || from + l > long_value.rows || v == NULL
		    || long_modulus != TMP.modulus())
			precondition_error_handler(j, "j", "0 <= j < columns",
					    long_value.columns, "columns", "",
					    l, "l", "0 <= l <= rows and l+from <= rows",
					    from, "from",
					    "0 <= from < rows and l + from <= rows",
					    long_value.rows, "rows", "",
					    PRT, "v", " v != NULL",
					    long_modulus, "modulus",
					    "modulus == global modulus",
					    TMP.modulus(), "global modulus", "",
					    "void base_matrix< bigmod >::"
					    "sto_column(const bigmod *v, lidia_size_t l, "
					    "lidia_size_t j, lidia_size_t from)",
					    DMESSAGE, EMESSAGE[1]);

		long *tmp = new long[l];
		memory_handler(tmp, DM_BM, "sto_column :: "
			       "Error in memory allocation (tmp)");
		for (i = 0; i < l; i++)
			(v[i].mantissa()).longify(tmp[i]);

		switch(long_value.bitfield.get_representation()) {
		case matrix_flags::dense_representation:
			D_base_modul_long.sto_column(long_value, tmp, l, j, from);
			break;
		case matrix_flags::sparse_representation:
			S_base_modul_long.sto_column(long_value, tmp, l, j, from);
			break;
		}
		delete[] tmp;
	}
}



void base_matrix< bigmod >::
sto_column_vector(const base_vector< bigmod > &v, lidia_size_t l,
		  lidia_size_t j, lidia_size_t from)
{
	//
	//       Task: A.sto_column_vector(v,l,j,p) stores v[0],...,v[l-1] in
	//             column j of matrix A starting at position p
	// Conditions: 0 <= j < A.columns and
	//             0 <= l <= min(v.size, A.rows) and
	//             0 <= from < A.rows and
	//             from + l <= A.rows
	//    Version: 2.0
	//

	debug_handler_l(DMESSAGE, "sto_column_vector(const base_vector< bigmod > &, "
			"lidia_size_t, lidia_size_t, lidia_size_t)", DVALUE + 2);

	lidia_size_t i;
	bigmod TMP;
	if (BigintLong) {
		lidia_size_t min = (bigint_value.rows < v.get_size()
				    ? bigint_value.rows : v.get_size());

		if (j >= bigint_value.columns || j< 0 || l < 0 || l > min
		    || from< 0 || from >= bigint_value.rows
		    || l + from > bigint_value.rows
		    || TMP.modulus() != bigint_modulus)
			precondition_error_handler(j , "j", "0 <= j < columns",
					    bigint_value.columns, "columns", "",
					    l, "l",
					    "0 < l <= min(v.size, rows) and from + l <= rows",
					    from, "from",
					    "0 <= from < rows and from + l <= rows",
					    bigint_value.rows, "rows", "",
					    v.get_size(), "v.size",
					    "0 < l <= min(v.size, rows)",
					    bigint_modulus, "modulus", "global modulus",
					    TMP.modulus(), "global modulus", "",
					    "void base_matrix< bigmod >::"
					    "sto_column_vector(const base_vector< bigmod > &v, "
					    "lidia_size_t l, lidia_size_t j, "
					    "lidia_size_t from)", DMESSAGE, EMESSAGE[3]);

		bigint *tmp = new bigint[l];
		memory_handler(tmp, DM_BM, "sto_column_vector :: "
			       "Error in memory allocation (tmp)");
		for (i = 0; i < l; i++)
			tmp[i] = v[i].mantissa();

		switch(bigint_value.bitfield.get_representation()) {
		case matrix_flags::dense_representation:
			D_base_modul_bigint.sto_column(bigint_value, tmp, l, j, from);
			break;
		case matrix_flags::sparse_representation:
			S_base_modul_bigint.sto_column(bigint_value, tmp, l, j, from);
			break;
		}
		delete[] tmp;
	}
	else {
		lidia_size_t min = (long_value.rows < v.get_size()
				    ? long_value.rows : v.get_size());

		if (j >= long_value.columns || j< 0 || l < 0 || l > min
		    || from< 0 || from >= long_value.rows
		    || l + from > long_value.rows
		    || bigint(long_modulus) != TMP.modulus())
			precondition_error_handler(j , "j", "0 <= j < columns",
					    long_value.columns, "columns", "",
					    l, "l",
					    "0 < l <= min(v.size, rows) and from + l <= rows",
					    from, "from",
					    "0 <= from < rows and from + l <= rows",
					    long_value.rows, "rows", "",
					    v.get_size(), "v.size",
					    "0 < l <= min(v.size, rows)",
					    bigint_modulus, "modulus", "global modulus",
					    TMP.modulus(), "global modulus", "",
					    "void base_matrix< bigmod >::"
					    "sto_column_vector(const base_vector< bigmod > &v, "
					    "lidia_size_t l, lidia_size_t j, "
					    "lidia_size_t from)", DMESSAGE, EMESSAGE[3]);

		long *tmp = new long[l];
		memory_handler(tmp, DM_BM, "sto_column_vector :: "
			       "Error in memory allocation (tmp)");
		for (i = 0; i < l; i++)
			(v[i].mantissa()).longify(tmp[i]);

		switch(long_value.bitfield.get_representation()) {
		case matrix_flags::dense_representation:
			D_base_modul_long.sto_column(long_value, tmp, l, j, from);
			break;
		case matrix_flags::sparse_representation:
			S_base_modul_long.sto_column(long_value, tmp, l, j, from);
			break;
		}
		delete[] tmp;
	}
}



void base_matrix< bigmod >::
get_column(bigmod *RES, lidia_size_t i) const
{
	//
	//       Task: A.get_column(RES,i);
	//             => RES[0] = A.value[0][i],...,
	//                RES[A.rows - 1] = A.value[A.rows - 1][i]
	// Conditions: 0 <= i < A.columns and
	//             RES != NULL
	//    Version: 2.0
	//

	debug_handler_l(DMESSAGE, "get_column(T *, lidia_size_t)", DVALUE + 2);

	if (BigintLong) {
		if (i >= bigint_value.columns || i < 0 || RES == NULL)
			precondition_error_handler(i, "i", "0 <= i < columns",
					    bigint_value.columns, "columns", "",
					    PRT, "RES", "RES != NULL",
					    "void base_matrix< bigmod >::"
					    "get_column(T *RES, lidia_size_t i) const",
					    DMESSAGE, EMESSAGE[1]);

		bigint *tmp = new bigint[bigint_value.rows];
		memory_handler(tmp, DM_BM, "get_column :: "
			       "Error in memory allocation (tmp)");

		switch(bigint_value.bitfield.get_representation()) {
		case matrix_flags::dense_representation:
			D_base_modul_bigint.get_column(bigint_value, tmp, i);
			break;
		case matrix_flags::sparse_representation:
			S_base_modul_bigint.get_column(bigint_value, tmp, i);
			break;
		}

		for (i = 0; i < bigint_value.rows; i++) {
			RES[i].set_modulus(bigint_modulus);
			RES[i] = tmp[i];
		}
		delete[] tmp;
	}
	else {
		if (i >= long_value.columns || i < 0 || RES == NULL)
			precondition_error_handler(i, "i", "0 <= i < columns",
					    long_value.columns, "columns", "",
					    PRT, "RES", "RES != NULL",
					    "void base_matrix< bigmod >::"
					    "get_column(T *RES, lidia_size_t i) const",
					    DMESSAGE, EMESSAGE[1]);

		long *tmp = new long[long_value.rows];
		memory_handler(tmp, DM_BM, "get_column :: "
			       "Error in memory allocation (tmp)");

		switch(long_value.bitfield.get_representation()) {
		case matrix_flags::dense_representation:
			D_base_modul_long.get_column(long_value, tmp, i);
			break;
		case matrix_flags::sparse_representation:
			S_base_modul_long.get_column(long_value, tmp, i);
			break;
		}
		for (i = 0; i < long_value.rows; i++) {
			RES[i].set_modulus(long_modulus);
			RES[i] = tmp[i];
		}
		delete[] tmp;
	}
}



void base_matrix< bigmod >::
get_column_vector(base_vector< bigmod > &RES, lidia_size_t i) const
{
	//
	//       Task: A.get_column_vector(RES,i);
	//             => RES[0] = A.value[0][i],...,
	//                RES[A.rows - 1] = A.value[A.rows - 1][i]
	// Conditions: 0 <= i < A.columns
	//    Version: 2.0
	//

	debug_handler_l(DMESSAGE, "get_column_vector(base_vector< bigmod > &, "
			"lidia_size_t)", DVALUE + 2);

	if (BigintLong) {
		if (i >= bigint_value.columns || i < 0)
			precondition_error_handler(i, "i", "0 <= i < columns",
					    bigint_value.columns, "columns", "",
					    "void base_matrix< bigmod >::"
					    "get_column_vector(base_vector< bigmod > &RES, "
					    "lidia_size_t i) const",
					    DMESSAGE, EMESSAGE[3]);

		if (RES.get_capacity() < bigint_value.rows)
			RES.set_capacity(bigint_value.rows);
		if (RES.size() != bigint_value.rows)
			RES.set_size(bigint_value.rows);

		bigint *tmp = new bigint[bigint_value.rows];
		memory_handler(tmp, DM_BM, "get_column_vector :: "
			       "Error in memory allocation (tmp)");

		switch(bigint_value.bitfield.get_representation()) {
		case matrix_flags::dense_representation:
			D_base_modul_bigint.get_column(bigint_value, tmp, i);
			break;
		case matrix_flags::sparse_representation:
			S_base_modul_bigint.get_column(bigint_value, tmp, i);
			break;
		}
		for (i = 0; i < bigint_value.rows; i++) {
			RES[i].set_modulus(bigint_modulus);
			RES[i] = tmp[i];
		}
		delete[] tmp;
	}
	else {
		if (i >= long_value.columns || i < 0)
			precondition_error_handler(i, "i", "0 <= i < columns",
					    long_value.columns, "columns", "",
					    "void base_matrix< bigmod >::"
					    "get_column_vector(base_vector< bigmod > &RES, "
					    "lidia_size_t i) const",
					    DMESSAGE, EMESSAGE[3]);

		if (RES.get_capacity() < long_value.rows)
			RES.set_capacity(long_value.rows);
		if (RES.size() != long_value.rows)
			RES.set_size(long_value.rows);

		long *tmp = new long[long_value.rows];
		memory_handler(tmp, DM_BM, "get_column_vector :: "
			       "Error in memory allocation (tmp)");

		switch(long_value.bitfield.get_representation()) {
		case matrix_flags::dense_representation:
			D_base_modul_long.get_column(long_value, tmp, i);
			break;
		case matrix_flags::sparse_representation:
			S_base_modul_long.get_column(long_value, tmp, i);
			break;
		}
		for (i = 0; i < long_value.rows; i++) {
			RES[i].set_modulus(long_modulus);
			RES[i] = tmp[i];
		}
		delete[] tmp;
	}
}



//
// row access
//

void base_matrix< bigmod >::
sto_row(const bigmod *v, lidia_size_t l, lidia_size_t i, lidia_size_t from)
{
	//
	//       Task: A.sto_row(v,l,j,p) stores v[0],...,v[l-1] in row i
	//             of matrix A starting at position p.
	// Conditions: 0 <= i < A.rows and
	//             0 <= l <= A.columns and
	//             0 <= from < A.columns and
	//             v != NULL and
	//             from + l <= A.rows
	//    Version: 2.0
	//

	debug_handler_l(DMESSAGE, "sto_row(const bigmod *, lidia_size_t, "
			"lidia_size_t, lidia_size_t)", DVALUE + 2);

	bigmod TMP;
	if (BigintLong) {
		if (i< 0 || i >= bigint_value.rows
		    || 0 > from || from >= bigint_value.columns
		    || l >= bigint_value.columns || l <= 0
		    || from + l > bigint_value.columns || v == NULL
		    || TMP.modulus() != bigint_modulus)
			precondition_error_handler(i, "i", "0 <= i < rows",
					    bigint_value.rows, "rows", "",
					    l, "l",
					    "0 <= l <= columns and l + from <= columns" ,
					    from, "from",
					    "0 <= from < columns and l + from <= columns",
					    bigint_value.columns, "columns", "",
					    PRT, "v", "v != NULL",
					    bigint_modulus, "modulus",
					    "modulus == global modulus",
					    TMP.modulus(), "global modulus", "",
					    "void base_matrix< bigmod >::"
					    "sto_row(const bigmod *v, lidia_size_t l, "
					    "lidia_size_t i, lidia_size_t)",
					    DMESSAGE, EMESSAGE[1]);

		bigint *tmp = new bigint[l];
		memory_handler(tmp, DM_BM, "sto_row :: "
			       "Error in memory allocation (tmp)");
		for (i = 0; i < l; i++)
			tmp[i] = v[i].mantissa();

		switch(bigint_value.bitfield.get_representation()) {
		case matrix_flags::dense_representation:
			D_base_modul_bigint.sto_row(bigint_value, tmp, l, i, from);
			break;
		case matrix_flags::sparse_representation:
			S_base_modul_bigint.sto_row(bigint_value, tmp, l, i, from);
			break;
		}
		delete[] tmp;
	}
	else {
		if (i< 0 || i >= long_value.rows
		    || 0 > from || from >= long_value.columns
		    || l >= long_value.columns || l <= 0
		    || from + l > long_value.columns || v == NULL
		    || bigint(long_modulus) != TMP.modulus())
			precondition_error_handler(i, "i", "0 <= i < rows",
					    long_value.rows, "rows", "",
					    l, "l",
					    "0 <= l <= columns and l + from <= columns" ,
					    from, "from",
					    "0 <= from < columns and l + from <= columns",
					    long_value.columns, "columns", "",
					    PRT, "v", "v != NULL",
					    long_modulus, "modulus",
					    "modulus == global modulus",
					    TMP.modulus(), "global modulus", "",
					    "void base_matrix< bigmod >::"
					    "sto_row(const bigmod *v, lidia_size_t l, "
					    "lidia_size_t i, lidia_size_t)",
					    DMESSAGE, EMESSAGE[1]);

		long *tmp = new long[l];
		memory_handler(tmp, DM_BM, "sto_row :: "
			       "Error in memory allocation (tmp)");
		for (i = 0; i < l; i++)
			(v[i].mantissa()).longify(tmp[i]);

		switch(long_value.bitfield.get_representation()) {
		case matrix_flags::dense_representation:
			D_base_modul_long.sto_row(long_value, tmp, l, i, from);
			break;
		case matrix_flags::sparse_representation:
			S_base_modul_long.sto_row(long_value, tmp, l, i, from);
			break;
		}
		delete[] tmp;
	}
}



void base_matrix< bigmod >::
sto_row_vector(const base_vector< bigmod > &v, lidia_size_t l,
	       lidia_size_t i, lidia_size_t from)
{
	//
	//       Task: A.sto_row(v,l,j,p) stores v[0],...,v[l-1] in row i
	//             of matrix A starting at position p.
	// Conditions: 0 <= i < A.rows and
	//             0 <= l <= min(v.size, A.columns) and
	//             0 <= from < A.rows and
	//             from + l <= A.columns
	//    Version: 2.0
	//

	debug_handler_l(DMESSAGE, "sto_row_vector(const base_vector< bigmod > &, "
			"lidia_size_t, lidia_size_t, lidia_size_t)", DVALUE + 2);

	bigmod TMP;
	if (BigintLong) {
		lidia_size_t min = (bigint_value.columns < v.get_size()
				    ? bigint_value.columns : v.get_size());

		if (i >= bigint_value.rows || i< 0 || l < 0 || l > min
		    || from< 0 || from >= bigint_value.columns
		    || l + from > bigint_value.columns
		    || bigint_modulus != TMP.modulus())
			precondition_error_handler(i, "i", "0 <= i < rows",
					    bigint_value.rows, "rows", "",
					    l, "l", "0 < l <= min(columns, v.size) "
					    "and from + l <= columns",
					    from, "from",
					    "0 <= from < columns and from + l <= columns",
					    bigint_value.columns, "columns", "",
					    v.size(), "v.size",
					    "0 < l <= min(columns, v.size)",
					    bigint_modulus, "modulus",
					    "modulus == global modulus",
					    TMP.modulus(), "global modulus", "",
					    "void base_matrix< bigmod >::"
					    "sto_row_vector(const base_vector< bigmod > &v, "
					    "lidia_size_t l, lidia_size_t i, "
					    "lidia_size_t from)", DMESSAGE, EMESSAGE[3]);

		bigint *tmp = new bigint[l];
		memory_handler(tmp, DM_BM, "sto_row_vector :: "
			       "Error in memory allocation (tmp)");
		for (i = 0; i < l; i++)
			tmp[i] = v[i].mantissa();

		switch(bigint_value.bitfield.get_representation()) {
		case matrix_flags::dense_representation:
			D_base_modul_bigint.sto_row(bigint_value, tmp, l, i, from);
			break;
		case matrix_flags::sparse_representation:
			S_base_modul_bigint.sto_row(bigint_value, tmp, l, i, from);
			break;
		}
		delete[] tmp;
	}
	else {
		lidia_size_t min = (long_value.columns < v.get_size()
				    ? long_value.columns : v.get_size());

		if (i >= long_value.rows || i< 0 || l < 0 || l > min
		    || from< 0 || from >= long_value.columns
		    || l + from > long_value.columns
		    || bigint(long_modulus) != TMP.modulus())
			precondition_error_handler(i, "i", "0 <= i < rows",
					    long_value.rows, "rows", "",
					    l, "l", "0 < l <= min(columns, v.size) "
					    "and from + l <= columns",
					    from, "from",
					    "0 <= from < columns and from + l <= columns",
					    long_value.columns, "columns", "",
					    v.size(), "v.size",
					    "0 < l <= min(columns, v.size)",
					    long_modulus, "modulus",
					    "modulus == global modulus",
					    TMP.modulus(), "global modulus", "",
					    "void base_matrix< bigmod >::"
					    "sto_row_vector(const base_vector< bigmod > &v, "
					    "lidia_size_t l, lidia_size_t i, "
					    "lidia_size_t from)", DMESSAGE, EMESSAGE[3]);

		long *tmp = new long[l];
		memory_handler(tmp, DM_BM, "sto_row_vector :: "
			       "Error in memory allocation (tmp)");
		for (i = 0; i < l; i++)
			(v[i].mantissa()).longify(tmp[i]);

		switch(long_value.bitfield.get_representation()) {
		case matrix_flags::dense_representation:
			D_base_modul_long.sto_row(long_value, tmp, l, i, from);
			break;
		case matrix_flags::sparse_representation:
			S_base_modul_long.sto_row(long_value, tmp, l, i, from);
			break;
		}
		delete[] tmp;
	}
}



void base_matrix< bigmod >::
get_row(bigmod *RES, lidia_size_t i) const
{
	//
	//       Task: A.get_row(RES,i);
	//             => RES[0] = A.value[i][0],...,
	//                RES[A.columns - 1] = A.value[i][A.columns - 1]
	// Conditions: 0 <= i < A.rows and
	//             RES != NULL
	//    Version: 2.0
	//

	debug_handler_l(DMESSAGE, "get_row(T *, lidia_size_t)", DVALUE + 2);

	if (BigintLong) {
		if (i >= bigint_value.rows || i < 0 || RES == NULL)
			precondition_error_handler(i, "i", "0 <= i < rows",
					    bigint_value.rows, "rows", "",
					    PRT, "RES", "RES != NULL",
					    "void base_matrix< bigmod >::"
					    "get_row(T *RES, lidia_size_t i) const",
					    DMESSAGE, EMESSAGE[3]);

		bigint *tmp = new bigint[bigint_value.columns];
		memory_handler(tmp, DM_BM, "get_row :: "
			       "Error in memory allocation (tmp)");

		switch(bigint_value.bitfield.get_representation()) {
		case matrix_flags::dense_representation:
			D_base_modul_bigint.get_row(bigint_value, tmp, i);
			break;
		case matrix_flags::sparse_representation:
			S_base_modul_bigint.get_row(bigint_value, tmp, i);
			break;
		}

		for (i = 0; i < bigint_value.columns; i++) {
			RES[i].set_modulus(bigint_modulus);
			RES[i] = tmp[i];
		}
		delete[] tmp;
	}
	else {
		if (i >= long_value.rows || i < 0 || RES == NULL)
			precondition_error_handler(i, "i", "0 <= i < rows",
					    long_value.rows, "rows", "",
					    PRT, "RES", "RES != NULL",
					    "void base_matrix< bigmod >::"
					    "get_row(T *RES, lidia_size_t i) const",
					    DMESSAGE, EMESSAGE[3]);

		long *tmp = new long[long_value.columns];
		memory_handler(tmp, DM_BM, "get_row :: "
			       "Error in memory allocation (tmp)");

		switch(long_value.bitfield.get_representation()) {
		case matrix_flags::dense_representation:
			D_base_modul_long.get_row(long_value, tmp, i);
			break;
		case matrix_flags::sparse_representation:
			S_base_modul_long.get_row(long_value, tmp, i);
			break;
		}

		for (i = 0; i < long_value.columns; i++) {
			RES[i].set_modulus(long_modulus);
			RES[i] = tmp[i];
		}
		delete[] tmp;
	}
}



void base_matrix< bigmod >::
get_row_vector(base_vector< bigmod > &RES, lidia_size_t i) const
{
	//
	//       Task: A.get_row_vector(RES, i);
	//             => RES[0] = A.value[i][0],...,
	//                RES[A.columns - 1] = A.value[i][A.columns - 1]
	// Conditions: 0 <= i < A.rows
	//    Version: 2.0
	//

	debug_handler_l(DMESSAGE, "get_row_vector(base_vector< bigmod > &, "
			"lidia_size_t)", DVALUE + 2);

	if (BigintLong) {
		if (i >= bigint_value.rows || i < 0)
			precondition_error_handler(i, "i", "0 <= i < rows",
					    bigint_value.rows, "rows", "",
					    "void base_matrix< bigmod >::"
					    "get_row_vector(base_vector< bigmod > &RES, "
					    "lidia_size_t i) const", DMESSAGE, EMESSAGE[3]);

		if (RES.get_capacity() < bigint_value.columns)
			RES.set_capacity(bigint_value.columns);
		if (RES.size() != bigint_value.columns)
			RES.set_size(bigint_value.columns);

		bigint *tmp = new bigint[bigint_value.columns];
		memory_handler(tmp, DM_BM, "get_row_vector :: "
			       "Error in memory allocation (tmp)");

		switch(bigint_value.bitfield.get_representation()) {
		case matrix_flags::dense_representation:
			D_base_modul_bigint.get_row(bigint_value, tmp, i);
			break;
		case matrix_flags::sparse_representation:
			S_base_modul_bigint.get_row(bigint_value, tmp, i);
			break;
		}
		for (i = 0; i < bigint_value.columns; i++) {
			RES[i].set_modulus(bigint_modulus);
			RES[i] = tmp[i];
		}
		delete[] tmp;
	}
	else {
		if (i >= long_value.rows || i < 0)
			precondition_error_handler(i, "i", "0 <= i < rows",
					    long_value.rows, "rows", "",
					    "void base_matrix< bigmod >::"
					    "get_row_vector(base_vector< bigmod > &RES, "
					    "lidia_size_t i) const", DMESSAGE, EMESSAGE[3]);

		if (RES.get_capacity() < long_value.columns)
			RES.set_capacity(long_value.columns);
		if (RES.size() != long_value.columns)
			RES.set_size(long_value.columns);

		long *tmp = new long[long_value.columns];
		memory_handler(tmp, DM_BM, "get_row_vector :: "
			       "Error in memory allocation (tmp)");

		switch(long_value.bitfield.get_representation()) {
		case matrix_flags::dense_representation:
			D_base_modul_long.get_row(long_value, tmp, i);
			break;
		case matrix_flags::sparse_representation:
			S_base_modul_long.get_row(long_value, tmp, i);
			break;
		}
		for (i = 0; i < long_value.columns; i++) {
			RES[i].set_modulus(long_modulus);
			RES[i] = tmp[i];
		}
		delete[] tmp;
	}
}



//
// array access
//

bigmod ** base_matrix< bigmod >::
get_data() const
{
	//
	//    Task: prt = A.get_value()
	//          => prt = pointer to a copy of array A.value
	// Version: 2.0
	//

	debug_handler_l(DMESSAGE, "get_value()", DVALUE + 2);

	if (BigintLong) {
		bigmod **copy = new bigmod *[bigint_value.rows];
		memory_handler(copy, DM_BM, "get_value :: "
			       "Error in memory allocation (copy)");

		lidia_size_t i, j, l;
		bigmod *tmp;
		bigint *tmp1;

		switch(bigint_value.bitfield.get_representation()) {
		case matrix_flags::dense_representation:
			for (i = 0; i < bigint_value.rows; i++) {
				tmp1 = bigint_value.value[i];
				tmp = new bigmod[bigint_value.columns];
				memory_handler(tmp, DM_BM, "get_value :: "
					       "Error in memory allocation (tmp)");
				for (j = 0; j < bigint_value.columns; j++) {
					tmp[j].set_modulus(bigint_modulus);
					tmp[j] = tmp1[j];
				}
				copy[i] = tmp;
			}
			break;
		case matrix_flags::sparse_representation:
			for (i = 0; i < bigint_value.rows; i++) {
				l = 0;
				tmp = new bigmod[bigint_value.columns];
				memory_handler(tmp, DM_BM, "get_value :: "
					       "Error in memory allocation (tmp)");
				for (j = 0; j < bigint_value.columns; j++) {
					tmp[j].set_modulus(bigint_modulus);
					if (l < bigint_value.value_counter[i]
					    && bigint_value.index[i][l] == j) {
						tmp[j] = bigint_value.value[i][l];
						l++;
					}
					else
						tmp[j] = bigint_value.Zero;
				}
				copy[i] = tmp;
			}
			break;
		}
		return copy;
	}
	else {
		bigmod **copy = new bigmod *[long_value.rows];
		memory_handler(copy, DM_BM, "get_value :: "
			       "Error in memory allocation (copy)");

		lidia_size_t i, j, l;
		bigmod *tmp;
		long *tmp1;

		switch(long_value.bitfield.get_representation()) {
		case matrix_flags::dense_representation:
			for (i = 0; i < long_value.rows; i++) {
				tmp1 = long_value.value[i];
				tmp = new bigmod[long_value.columns];
				memory_handler(tmp, DM_BM, "get_value :: "
					       "Error in memory allocation (tmp)");
				for (j = 0; j < long_value.columns; j++) {
					tmp[j].set_modulus(long_modulus);
					tmp[j] = tmp1[j];
				}
				copy[i] = tmp;
			}
			break;
		case matrix_flags::sparse_representation:
			for (i = 0; i < long_value.rows; i++) {
				l = 0;
				tmp = new bigmod[long_value.columns];
				memory_handler(tmp, DM_BM, "get_value :: "
					       "Error in memory allocation (tmp)");
				for (j = 0; j < long_value.columns; j++) {
					tmp[j].set_modulus(long_modulus);
					if (l < long_value.value_counter[i]
					    && long_value.index[i][l] == j) {
						tmp[j] = long_value.value[i][l];
						l++;
					}
					else
						tmp[j] = long_value.Zero;
				}
				copy[i] = tmp;
			}
			break;
		}
		return copy;
	}
}



void base_matrix< bigmod >::
set_data(const bigint **v, lidia_size_t r, lidia_size_t c)
{
	//
	//       Task: A.set_value(v,r,c)
	//             => A.value[i][j] = v[i][j]
	//             0 <= i < r and 0 <= j < c
	// Conditions: r >= 0 and
	//             c >= 0 and
	//             v != NULL
	//    Version: 2.0
	//

	debug_handler_l(DMESSAGE, "set_value(const bigint **, "
			"lidia_size_t, lidia_size_t)", DVALUE + 2);

	lidia_size_t len, len1;
	if (BigintLong) {
		if (r < 0 || c < 0 || v == NULL)
			precondition_error_handler(r, "r", "0 <= r",
					    c, "c", "0 <= c",
					    PRT, "v", "v != NULL",
					    "void base_matrix< bigmod >::"
					    "set_data(const bigmod **v, lidia_size_t r, "
					    "lidia_size_t c)", DMESSAGE, EMESSAGE[1]);

		switch(bigint_value.bitfield.get_representation()) {
		case matrix_flags::dense_representation:
			if (r != bigint_value.rows)
				D_base_modul_bigint.set_no_of_rows(bigint_value, r);
			if (c != bigint_value.columns)
				D_base_modul_bigint.set_no_of_columns(bigint_value, c);

			for (len1 = bigint_value.rows - 1; len1 >= 0; len1--)
				for (len = bigint_value.columns - 1; len >= 0; len--)
					LiDIA::best_remainder(bigint_value.value[len1][len], v[len1][len],
								 bigint_modulus);
			break;
		case matrix_flags::sparse_representation:
			if (r != bigint_value.rows)
				S_base_modul_bigint.set_no_of_rows(bigint_value, r);
			if (c != bigint_value.columns)
				S_base_modul_bigint.set_no_of_columns(bigint_value, c);

			bigint TMP;
			for (len1 = bigint_value.rows - 1; len1 >= 0; len1--)
				for (len = bigint_value.columns - 1; len >= 0; len--) {
					LiDIA::best_remainder(TMP, v[len1][len], bigint_modulus);
					S_base_modul_bigint.sto(bigint_value, len1, len, v[len1][len]);
				}
			break;
		}
	}
	else {
		if (r < 0 || c < 0 || v == NULL)
			precondition_error_handler(r, "r", "0 <= r",
					    c, "c", "0 <= c",
					    PRT, "v", "v != NULL",
					    "void base_matrix< bigmod >::"
					    "set_data(const bigmod **v, lidia_size_t r, "
					    "lidia_size_t c)", DMESSAGE, EMESSAGE[1]);

		switch(long_value.bitfield.get_representation()) {
		case matrix_flags::dense_representation:
			if (r != long_value.rows)
				D_base_modul_long.set_no_of_rows(long_value, r);
			if (c != long_value.columns)
				D_base_modul_long.set_no_of_columns(long_value, c);

			for (len1 = long_value.rows - 1; len1 >= 0; len1--)
				for (len = long_value.columns - 1; len >= 0; len--)
					LiDIA::best_remainder(long_value.value[len1][len], v[len1][len],
								 long_modulus);
			break;
		case matrix_flags::sparse_representation:
			if (r != long_value.rows)
				S_base_modul_long.set_no_of_rows(long_value, r);
			if (c != long_value.columns)
				S_base_modul_long.set_no_of_columns(long_value, c);

			long TMP;
			for (len1 = long_value.rows - 1; len1 >= 0; len1--)
				for (len = long_value.columns - 1; len >= 0; len--) {
					LiDIA::best_remainder(TMP, v[len1][len], long_modulus);
					S_base_modul_long.sto(long_value, len1, len, TMP);
				}
			break;
		}
	}
}



//
// set_modulus
//

void base_matrix< bigmod >::
set_modulus(const bigint &mod)
{
	if (mod.bit_length() > bigint::bits_per_digit()) {
		if (!BigintLong) {
			BigintLong = true;
			bigint_modulus = mod;
			bigint_value.bitfield = long_value.bitfield;
			if (bigint_value.bitfield.get_representation() ==
			    matrix_flags::dense_representation) {
				D_base_modul_bigint.constructor(bigint_value, long_value.rows, long_value.columns);
				D_base_modul_long.destructor(long_value);
			}
			else {
				S_base_modul_bigint.constructor(bigint_value, long_value.rows, long_value.columns);
				S_base_modul_long.destructor(long_value);
			}

		}
		else
			bigint_modulus = mod;
	}
	else {
		if (BigintLong) {
			BigintLong = false;
			mod.longify(long_modulus);
			long_value.bitfield = bigint_value.bitfield;
			if (long_value.bitfield.get_representation() ==
			    matrix_flags::dense_representation) {
				D_base_modul_long.constructor(long_value, bigint_value.rows, bigint_value.columns);
				D_base_modul_bigint.destructor(bigint_value);
			}
			else {
				S_base_modul_long.constructor(long_value, bigint_value.rows, bigint_value.columns);
				S_base_modul_bigint.destructor(bigint_value);
			}
		}
		else
			mod.longify(long_modulus);
	}
}



///////////////////////////
// END: access functions //
///////////////////////////

//
// insert_at
//

void base_matrix< bigmod >::
insert_at(lidia_size_t r, lidia_size_t c, const base_matrix< bigmod > &A,
	  lidia_size_t startr, lidia_size_t startc,
	  lidia_size_t l1, lidia_size_t l2)
{
	//
	//       Task: A.insert_at(r,c,B,startr,startc,l1,l2)
	//             stores the submatrix B[startr, startc]...
	//             B[startr+l1][startc+l2]
	//             in matrix A starting at position (r,c)
	// Conditions: 0 <= r < A.rows and
	//             0 <= c < A.columns and
	//             0 <= startr < B.rows and
	//             0 <= startc < B.columns and
	//             0 <= l1 <= B.rows and
	//             0 <= l2 <= B.columns and
	//             startr + l1 <= B.rows and
	//             startc + l2 <= B.columns and
	//             r + l1 <= A.rows and
	//             c + l2 <= A.columns
	//             v != NULL
	// Version: 2.0
	//

	bigint MOD1 = A.get_modulus();
	bigint MOD2 = get_modulus();

	if (BigintLong) {
		if (r< 0 || r >= bigint_value.rows || c < 0
		    || c >= bigint_value.columns
		    || startr< 0 || startr >= A.bigint_value.rows || startc < 0
		    || startc >= A.bigint_value.columns
		    || l1< 0 || l1 > A.bigint_value.rows || l2 < 0
		    || l2 > A.bigint_value.columns
		    || startr + l1 > A.bigint_value.rows
		    || startc + l2 > A.bigint_value.columns
		    || r + l1 > bigint_value.rows || c + l2 > bigint_value.columns
		    || BigintLong != A.BigintLong
		    || MOD1 != MOD2)
			precondition_error_handler(r, "r", "0 <= r < rows and r + l1 <= rows",
					    bigint_value.rows, "rows", "",
					    c, "c", "0 <= c < columns and c + l2 <= columns",
					    bigint_value.columns, "columns", "",
					    startr, "startr",
					    "0 <= startr < A.rows and startr + l1 <= A.rows",
					    A.bigint_value.rows, "A.rows", "",
					    startc, "startc",
					    "0 <= startc < A.columns "
					    "and startc + l2 <= A.columns",
					    A.bigint_value.columns, "A.columns", "",
					    l1, "l1",
					    "r + l1 <= rows and startr + l1 <= A.rows",
					    l2, "l2",
					    "c + l2 <= columns and startc + l2 <= A.columns",
					    MOD1, "A.modulus", "modulus == A.modulus",
					    MOD2, "modulus", "",
					    "void base_matrix< bigmod >::"
					    "insert_at(lidia_size_t r, lidia_size_t c, "
					    "const base_matrix< bigmod > &A, "
					    "lidia_size_t startr, lidia_size_t startc, "
					    "lidia_size_t l1, lidia_size_t l2)",
					    DMESSAGE, EMESSAGE[1]);

		switch(bigint_value.bitfield.get_representation()) {
		case matrix_flags::dense_representation:
			switch(A.bigint_value.bitfield.get_representation()) {
			case matrix_flags::dense_representation:
				DD_base_modul_bigint.insert_at(bigint_value, r, c,
							       A.bigint_value, startr, startc,
							       l1, l2);
				break;
			case matrix_flags::sparse_representation:
				DS_base_modul_bigint.insert_at(bigint_value, r, c,
							       A.bigint_value, startr, startc,
							       l1, l2);
				break;
			}
			break;
		case matrix_flags::sparse_representation:
			switch(A.bigint_value.bitfield.get_representation()) {
			case matrix_flags::dense_representation:
				SD_base_modul_bigint.insert_at(bigint_value, r, c,
							       A.bigint_value, startr, startc,
							       l1, l2);
				break;
			case matrix_flags::sparse_representation:
				SS_base_modul_bigint.insert_at(bigint_value, r, c,
							       A.bigint_value, startr, startc,
							       l1, l2);
				break;
			}
			break;
		}
	}
	else {
		if (r< 0 || r >= long_value.rows || c < 0
		    || c >= long_value.columns
		    || startr< 0 || startr >= A.long_value.rows || startc < 0
		    || startc >= A.long_value.columns
		    || l1< 0 || l1 > A.long_value.rows || l2 < 0
		    || l2 > A.long_value.columns
		    || startr + l1 > A.long_value.rows
		    || startc + l2 > A.long_value.columns
		    || r + l1 > long_value.rows || c + l2 > long_value.columns
		    || MOD1 != MOD2)
			precondition_error_handler(r, "r", "0 <= r < rows and r + l1 <= rows",
					    long_value.rows, "rows", "",
					    c, "c", "0 <= c < columns and c + l2 <= columns",
					    long_value.columns, "columns", "",
					    startr, "startr",
					    "0 <= startr < A.rows and startr + l1 <= A.rows",
					    A.long_value.rows, "A.rows", "",
					    startc, "startc",
					    "0 <= startc < A.columns "
					    "and startc + l2 <= A.columns",
					    A.long_value.columns, "A.columns", "",
					    l1, "l1",
					    "r + l1 <= rows and startr + l1 <= A.rows",
					    l2, "l2",
					    "c + l2 <= columns and startc + l2 <= A.columns",
					    MOD1, "A.modulus", "modulus == A.modulus",
					    MOD2, "modulus", "",
					    "void base_matrix< bigmod >::"
					    "insert_at(lidia_size_t r, lidia_size_t c, "
					    "const base_matrix< bigmod > &A, "
					    "lidia_size_t startr, lidia_size_t startc, "
					    "lidia_size_t l1, lidia_size_t l2)",
					    DMESSAGE, EMESSAGE[1]);

		switch(long_value.bitfield.get_representation()) {
		case matrix_flags::dense_representation:
			switch(A.long_value.bitfield.get_representation()) {
			case matrix_flags::dense_representation:
				DD_base_modul_long.insert_at(long_value, r, c,
							     A.long_value, startr, startc,
							     l1, l2);
				break;
			case matrix_flags::sparse_representation:
				DS_base_modul_long.insert_at(long_value, r, c,
							     A.long_value, startr, startc,
							     l1, l2);
				break;
			}
			break;
		case matrix_flags::sparse_representation:
			switch(A.long_value.bitfield.get_representation()) {
			case matrix_flags::dense_representation:
				SD_base_modul_long.insert_at(long_value, r, c,
							     A.long_value, startr, startc,
							     l1, l2);
				break;
			case matrix_flags::sparse_representation:
				SS_base_modul_long.insert_at(long_value, r, c,
							     A.long_value, startr, startc,
							     l1, l2);
				break;
			}
			break;
		}
	}
}



//
// insert_columns, insert_rows, remove_columns, remove_rows
//

void base_matrix< bigmod >::
insert_columns(lidia_size_t *ind, const bigint **news)
{
	//
	//    Task: A.insert_columns(ind, news);
	//          => insert the columns of array news in matrix A
	// Version: 2.0
	//

	debug_handler_l(DMESSAGE, "insert_columns(lidia_size_t *, const bigmod **)",
			DVALUE + 2);

	lidia_size_t i, j;
	if (BigintLong) {
		bigint **tmp = new bigint *[ind[0]];
		memory_handler(tmp, DM_BM, "insert_columns :: "
			       "Error in memory allocation (tmp)");
		for (i = 0; i < ind[0]; i++) {
			tmp[i] = new bigint[bigint_value.rows];
			memory_handler(tmp[i], DM_BM, "insert_columns :: "
				       "Error in memory allocation (tmp[i])");

			for (j = 0; j < bigint_value.rows; j++)
				LiDIA::best_remainder(tmp[i][j], news[i][j], bigint_modulus);
		}

		switch(bigint_value.bitfield.get_representation()) {
		case matrix_flags::dense_representation:
			D_base_modul_bigint.insert_columns(bigint_value, ind,
							   const_cast<const bigint **>(tmp));
			break;
		case matrix_flags::sparse_representation:
			S_base_modul_bigint.insert_columns(bigint_value, ind,
							   const_cast<const bigint **>(tmp));
			break;
		}
	}
	else {
		long **tmp = new long *[ind[0]];
		memory_handler(tmp, DM_BM, "insert_columns :: "
			       "Error in memory allocation (tmp)");
		for (i = 0; i < ind[0]; i++) {
			tmp[i] = new long[long_value.rows];
			memory_handler(tmp[i], DM_BM, "insert_columns :: "
				       "Error in memory allocation (tmp[i])");

			for (j = 0; j < long_value.rows; j++)
				LiDIA::best_remainder(tmp[i][j], news[i][j], long_modulus);
		}

		switch(long_value.bitfield.get_representation()) {
		case matrix_flags::dense_representation:
			D_base_modul_long.insert_columns(long_value, ind,
							 const_cast<const long **>(tmp));
			break;
		case matrix_flags::sparse_representation:
			S_base_modul_long.insert_columns(long_value, ind,
							 const_cast<const long **>(tmp));
			break;
		}
	}
}



void base_matrix< bigmod >::
remove_columns(lidia_size_t *rem)
{
	//
	//    Task: A.remove_column(rem);
	//          => remove columns
	// Version: 2.0
	//

	debug_handler_l(DMESSAGE, "remove_column(lidia_size_t *)", DVALUE + 2);

	if (BigintLong) {
		switch(bigint_value.bitfield.get_representation()) {
		case matrix_flags::dense_representation:
			D_base_modul_bigint.remove_columns(bigint_value, rem);
			break;
		case matrix_flags::sparse_representation:
			S_base_modul_bigint.remove_columns(bigint_value, rem);
			break;
		}
	}
	else {
		switch(long_value.bitfield.get_representation()) {
		case matrix_flags::dense_representation:
			D_base_modul_long.remove_columns(long_value, rem);
			break;
		case matrix_flags::sparse_representation:
			S_base_modul_long.remove_columns(long_value, rem);
			break;
		}
	}
}



void base_matrix< bigmod >::
insert_rows(lidia_size_t *ind, const bigint **news)
{
	//
	//    Task: A.insert_rows(ind, news);
	//          => insert the rows of array in matrix A
	// Version: 2.0
	//

	debug_handler_l(DMESSAGE, "insert_rows(lidia_size_t *, const bigmod **)",
			DVALUE + 2);

	lidia_size_t i, j;
	if (BigintLong) {
		bigint **tmp = new bigint *[ind[0]];
		memory_handler(tmp, DM_BM, "insert_rows :: "
			       "Error in memory allocation (tmp)");

		for (i = 0; i < ind[0]; i++) {
			tmp[i] = new bigint[bigint_value.columns];
			memory_handler(tmp[i], DM_BM, "insert_rows :: "
				       "Error in memory allocation (tmp[i])");

			for (j = 0; j < bigint_value.columns; j++)
				LiDIA::best_remainder(tmp[i][j], news[i][j], bigint_modulus);
		}

		switch(bigint_value.bitfield.get_representation()) {
		case matrix_flags::dense_representation:
			D_base_modul_bigint.insert_rows(bigint_value, ind,
							const_cast<const bigint **>(tmp));
			break;
		case matrix_flags::sparse_representation:
			S_base_modul_bigint.insert_rows(bigint_value, ind,
							const_cast<const bigint **>(tmp));
			break;
		}
	}
	else {
		long **tmp = new long *[ind[0]];
		memory_handler(tmp, DM_BM, "insert_rows :: "
			       "Error in memory allocation (tmp)");
		for (i = 0; i < ind[0]; i++) {
			tmp[i] = new long[long_value.columns];
			memory_handler(tmp[i], DM_BM, "insert_rows :: "
				       "Error in memory allocation (tmp[i])");

			for (j = 0; j < long_value.columns; j++)
				LiDIA::best_remainder(tmp[i][j], news[i][j], long_modulus);
		}

		switch(long_value.bitfield.get_representation()) {
		case matrix_flags::dense_representation:
			D_base_modul_long.insert_rows(long_value, ind,
						      const_cast<const long **>(tmp));
			break;
		case matrix_flags::sparse_representation:
			S_base_modul_long.insert_rows(long_value, ind,
						      const_cast<const long **>(tmp));
			break;
		}
	}
}



void base_matrix< bigmod >::
remove_rows(lidia_size_t *rem)
{
	//
	//       Task: A.remove_rows(rem);
	//             => remove rows
	//    Version: 2.0
	//

	debug_handler_l(DMESSAGE, "remove_rows(lidia_size_t *)", DVALUE + 2);

	if (BigintLong) {
		switch(bigint_value.bitfield.get_representation()) {
		case matrix_flags::dense_representation:
			D_base_modul_bigint.remove_rows(bigint_value, rem);
			break;
		case matrix_flags::sparse_representation:
			S_base_modul_bigint.remove_rows(bigint_value, rem);
			break;
		}
	}
	else {
		switch(long_value.bitfield.get_representation()) {
		case matrix_flags::dense_representation:
			D_base_modul_long.remove_rows(long_value, rem);
			break;
		case matrix_flags::sparse_representation:
			S_base_modul_long.remove_rows(long_value, rem);
			break;
		}
	}
}



//
// split functions
//

void base_matrix< bigmod >::
split_t(base_matrix< bigmod > &A, base_matrix< bigmod > &B,
	base_matrix< bigmod > &C, base_matrix< bigmod > &D) const
{
	//
	//       Task: E.split_t(A,B,C,D);
	//             => splits matrix E into the submatrices A,B,C,D
	//                    (A B)
	//             => E = (C D)
	// Conditions: A.columns, B.columns, C.columns, D.columns <= columns and
	//             A.rows, B.rows, C.rows, D.rows <= rows
	//    Version: 2.0
	//

	debug_handler_l(DMESSAGE, "split_t(base_matrix< bigmod > &, "
			"base_matrix< bigmod > &, "
			"base_matrix< bigmod > &, base_matrix< bigmod > &)",
			DVALUE + 3);

	if (BigintLong) {
		if (A.bigint_value.rows > bigint_value.rows
		    || A.bigint_value.columns > bigint_value.columns
		    || B.bigint_value.rows > bigint_value.rows
		    || B.bigint_value.columns > bigint_value.columns
		    || C.bigint_value.rows > bigint_value.rows
		    || C.bigint_value.columns > bigint_value.columns
		    || D.bigint_value.rows > bigint_value.rows
		    || D.bigint_value.columns > bigint_value.columns)
			precondition_error_handler(A.bigint_value.rows, "A.rows", "A.rows <= rows",
					    A.bigint_value.columns, "A.columns",
					    "A.columns <= columns",
					    B.bigint_value.rows, "B.rows", "B.rows <= rows",
					    B.bigint_value.columns, "B.columns",
					    "B.columns <= columns",
					    C.bigint_value.rows, "C.rows", "C.rows <= rows",
					    C.bigint_value.columns, "C.columns",
					    "C.columns <= columns",
					    D.bigint_value.rows, "D.rows", "D.rows <= rows",
					    D.bigint_value.columns, "D.columns",
					    "D.columns <= columns",
					    bigint_value.rows, "rows", "",
					    bigint_value.columns, "columns", "",
					    "void base_matrix< bigmod >::"
					    "split_t(base_matrix< bigmod > &A, "
					    "base_matrix< bigmod > &B, "
					    "base_matrix< bigmod > &C, "
					    "base_matrix< bigmod > &D) const",
					    DMESSAGE, EMESSAGE[4]);

		A.set_modulus(bigint_modulus);
		B.set_modulus(bigint_modulus);
		C.set_modulus(bigint_modulus);
		D.set_modulus(bigint_modulus);

		switch(bigint_value.bitfield.get_representation()) {
		case matrix_flags::dense_representation:
			if (A.bigint_value.bitfield.get_representation() ==
			    matrix_flags::dense_representation)
				DD_base_modul_bigint.insert_at(A.bigint_value, 0, 0, bigint_value,
							       0, 0, A.bigint_value.rows,
							       A.bigint_value.columns);
			else
				DS_base_modul_bigint.insert_at(A.bigint_value, 0, 0, bigint_value,
							       0, 0, A.bigint_value.rows,
							       A.bigint_value.columns);

			if (B.bigint_value.bitfield.get_representation() ==
			    matrix_flags::dense_representation)
				DD_base_modul_bigint.insert_at(B.bigint_value, 0, 0,
							       bigint_value, 0,
							       bigint_value.columns -
							       B.bigint_value.columns,
							       B.bigint_value.rows,
							       B.bigint_value.columns);
			else
				DS_base_modul_bigint.insert_at(B.bigint_value, 0, 0,
							       bigint_value, 0,
							       bigint_value.columns -
							       B.bigint_value.columns,
							       B.bigint_value.rows,
							       B.bigint_value.columns);

			if (C.bigint_value.bitfield.get_representation() ==
			    matrix_flags::dense_representation)
				DD_base_modul_bigint.insert_at(C.bigint_value, 0, 0, bigint_value,
							       bigint_value.rows -
							       C.bigint_value.rows,
							       0, C.bigint_value.rows,
							       C.bigint_value.columns);
			else
				DS_base_modul_bigint.insert_at(C.bigint_value, 0, 0, bigint_value,
							       bigint_value.rows -
							       C.bigint_value.rows,
							       0, C.bigint_value.rows,
							       C.bigint_value.columns);

			if (C.bigint_value.bitfield.get_representation() ==
			    matrix_flags::dense_representation)
				DD_base_modul_bigint.insert_at(D.bigint_value, 0, 0, bigint_value,
							       bigint_value.rows -
							       D.bigint_value.rows,
							       bigint_value.columns -
							       D.bigint_value.columns,
							       D.bigint_value.rows,
							       D.bigint_value.columns);
			else
				DS_base_modul_bigint.insert_at(D.bigint_value, 0, 0, bigint_value,
							       bigint_value.rows -
							       D.bigint_value.rows,
							       bigint_value.columns -
							       D.bigint_value.columns,
							       D.bigint_value.rows,
							       D.bigint_value.columns);
			break;
		case matrix_flags::sparse_representation:
			if (A.bigint_value.bitfield.get_representation() ==
			    matrix_flags::dense_representation)
				SD_base_modul_bigint.insert_at(A.bigint_value, 0, 0, bigint_value,
							       0, 0, A.bigint_value.rows,
							       A.bigint_value.columns);
			else
				SS_base_modul_bigint.insert_at(A.bigint_value, 0, 0, bigint_value,
							       0, 0, A.bigint_value.rows,
							       A.bigint_value.columns);

			if (B.bigint_value.bitfield.get_representation() ==
			    matrix_flags::dense_representation)
				SD_base_modul_bigint.insert_at(B.bigint_value, 0, 0, bigint_value,
							       0, bigint_value.columns -
							       B.bigint_value.columns,
							       B.bigint_value.rows,
							       B.bigint_value.columns);
			else
				SS_base_modul_bigint.insert_at(B.bigint_value, 0, 0, bigint_value,
							       0, bigint_value.columns -
							       B.bigint_value.columns,
							       B.bigint_value.rows,
							       B.bigint_value.columns);

			if (C.bigint_value.bitfield.get_representation() ==
			    matrix_flags::dense_representation)
				SD_base_modul_bigint.insert_at(C.bigint_value, 0, 0, bigint_value,
							       bigint_value.rows -
							       C.bigint_value.rows, 0,
							       C.bigint_value.rows,
							       C.bigint_value.columns);
			else
				SS_base_modul_bigint.insert_at(C.bigint_value, 0, 0, bigint_value,
							       bigint_value.rows -
							       C.bigint_value.rows, 0,
							       C.bigint_value.rows,
							       C.bigint_value.columns);

			if (C.bigint_value.bitfield.get_representation() ==
			    matrix_flags::dense_representation)
				SD_base_modul_bigint.insert_at(D.bigint_value, 0, 0, bigint_value,
							       bigint_value.rows -
							       D.bigint_value.rows,
							       bigint_value.columns -
							       D.bigint_value.columns,
							       D.bigint_value.rows,
							       D.bigint_value.columns);
			else
				SS_base_modul_bigint.insert_at(D.bigint_value, 0, 0, bigint_value,
							       bigint_value.rows -
							       D.bigint_value.rows,
							       bigint_value.columns -
							       D.bigint_value.columns,
							       D.bigint_value.rows,
							       D.bigint_value.columns);
			break;
		}
	}
	else {
		if (A.long_value.rows > long_value.rows
		    || A.long_value.columns > long_value.columns
		    || B.long_value.rows > long_value.rows
		    || B.long_value.columns > long_value.columns
		    || C.long_value.rows > long_value.rows
		    || C.long_value.columns > long_value.columns
		    || D.long_value.rows > long_value.rows
		    || D.long_value.columns > long_value.columns)
			precondition_error_handler(A.long_value.rows, "A.rows", "A.rows <= rows",
					    A.long_value.columns, "A.columns",
					    "A.columns <= columns",
					    B.long_value.rows, "B.rows", "B.rows <= rows",
					    B.long_value.columns, "B.columns",
					    "B.columns <= columns",
					    C.long_value.rows, "C.rows", "C.rows <= rows",
					    C.long_value.columns, "C.columns",
					    "C.columns <= columns",
					    D.long_value.rows, "D.rows", "D.rows <= rows",
					    D.long_value.columns, "D.columns",
					    "D.columns <= columns",
					    long_value.rows, "rows", "",
					    long_value.columns, "columns", "",
					    "void base_matrix< bigmod >::"
					    "split_t(base_matrix< bigmod > &A, "
					    "base_matrix< bigmod > &B, "
					    "base_matrix< bigmod > &C, "
					    "base_matrix< bigmod > &D) const",
					    DMESSAGE, EMESSAGE[4]);

		A.set_modulus(long_modulus);
		B.set_modulus(long_modulus);
		C.set_modulus(long_modulus);
		D.set_modulus(long_modulus);

		switch(long_value.bitfield.get_representation()) {
		case matrix_flags::dense_representation:
			if (A.long_value.bitfield.get_representation() ==
			    matrix_flags::dense_representation)
				DD_base_modul_long.insert_at(A.long_value, 0, 0, long_value,
							     0, 0, A.long_value.rows,
							     A.long_value.columns);
			else
				DS_base_modul_long.insert_at(A.long_value, 0, 0, long_value,
							     0, 0, A.long_value.rows,
							     A.long_value.columns);

			if (B.long_value.bitfield.get_representation() ==
			    matrix_flags::dense_representation)
				DD_base_modul_long.insert_at(B.long_value, 0, 0,
							     long_value, 0,
							     long_value.columns -
							     B.long_value.columns,
							     B.long_value.rows,
							     B.long_value.columns);
			else
				DS_base_modul_long.insert_at(B.long_value, 0, 0,
							     long_value, 0,
							     long_value.columns -
							     B.long_value.columns,
							     B.long_value.rows,
							     B.long_value.columns);

			if (C.long_value.bitfield.get_representation() ==
			    matrix_flags::dense_representation)
				DD_base_modul_long.insert_at(C.long_value, 0, 0, long_value,
							     long_value.rows -
							     C.long_value.rows,
							     0, C.long_value.rows,
							     C.long_value.columns);
			else
				DS_base_modul_long.insert_at(C.long_value, 0, 0, long_value,
							     long_value.rows -
							     C.long_value.rows,
							     0, C.long_value.rows,
							     C.long_value.columns);

			if (C.long_value.bitfield.get_representation() ==
			    matrix_flags::dense_representation)
				DD_base_modul_long.insert_at(D.long_value, 0, 0, long_value,
							     long_value.rows -
							     D.long_value.rows,
							     long_value.columns -
							     D.long_value.columns,
							     D.long_value.rows,
							     D.long_value.columns);
			else
				DS_base_modul_long.insert_at(D.long_value, 0, 0, long_value,
							     long_value.rows -
							     D.long_value.rows,
							     long_value.columns -
							     D.long_value.columns,
							     D.long_value.rows,
							     D.long_value.columns);
			break;
		case matrix_flags::sparse_representation:
			if (A.long_value.bitfield.get_representation() ==
			    matrix_flags::dense_representation)
				SD_base_modul_long.insert_at(A.long_value, 0, 0, long_value,
							     0, 0, A.long_value.rows,
							     A.long_value.columns);
			else
				SS_base_modul_long.insert_at(A.long_value, 0, 0, long_value,
							     0, 0, A.long_value.rows,
							     A.long_value.columns);

			if (B.long_value.bitfield.get_representation() ==
			    matrix_flags::dense_representation)
				SD_base_modul_long.insert_at(B.long_value, 0, 0, long_value,
							     0, long_value.columns -
							     B.long_value.columns,
							     B.long_value.rows,
							     B.long_value.columns);
			else
				SS_base_modul_long.insert_at(B.long_value, 0, 0, long_value,
							     0, long_value.columns -
							     B.long_value.columns,
							     B.long_value.rows,
							     B.long_value.columns);

			if (C.long_value.bitfield.get_representation() ==
			    matrix_flags::dense_representation)
				SD_base_modul_long.insert_at(C.long_value, 0, 0, long_value,
							     long_value.rows -
							     C.long_value.rows, 0,
							     C.long_value.rows,
							     C.long_value.columns);
			else
				SS_base_modul_long.insert_at(C.long_value, 0, 0, long_value,
							     long_value.rows -
							     C.long_value.rows, 0,
							     C.long_value.rows,
							     C.long_value.columns);

			if (C.long_value.bitfield.get_representation() ==
			    matrix_flags::dense_representation)
				SD_base_modul_long.insert_at(D.long_value, 0, 0, long_value,
							     long_value.rows -
							     D.long_value.rows,
							     long_value.columns -
							     D.long_value.columns,
							     D.long_value.rows,
							     D.long_value.columns);
			else
				SS_base_modul_long.insert_at(D.long_value, 0, 0, long_value,
							     long_value.rows -
							     D.long_value.rows,
							     long_value.columns -
							     D.long_value.columns,
							     D.long_value.rows,
							     D.long_value.columns);
			break;
		}
	}
}



void base_matrix< bigmod >::
split_h(base_matrix< bigmod > &A, base_matrix< bigmod > &B) const
{
	//
	//       Task: C.split_h(A,B);
	//             => splits matrix C into the submatrices A,B
	//             => C = (A B)
	// Conditions: A.columns, B.columns <= columns or
	//             A.rows, B.rows <= rows
	//    Version: 2.0
	//

	debug_handler_l(DMESSAGE, "split_h(base_matrix< bigmod > &, "
			"base_matrix< bigmod > &)", DVALUE + 3);

	if (BigintLong) {
		if (A.bigint_value.rows > bigint_value.rows
		    || A.bigint_value.columns > bigint_value.columns
		    || B.bigint_value.rows > bigint_value.rows
		    || B.bigint_value.columns > bigint_value.columns)
			precondition_error_handler(A.bigint_value.rows, "A.rows", "A.rows <= rows",
					    A.bigint_value.columns, "A.columns",
					    "A.columns <= coluns",
					    B.bigint_value.rows, "B.rows", "B.rows <= rows",
					    B.bigint_value.columns, "B.columns",
					    "B.columns <= coluns",
					    bigint_value.rows, "rows", "",
					    bigint_value.columns, "columns", "",
					    "void base_matrix< bigmod >::"
					    "split_h(base_matrix< bigmod > &A, "
					    "base_matrix< bigmod > &B) const",
					    DMESSAGE, EMESSAGE[4]);

		A.set_modulus(bigint_modulus);
		B.set_modulus(bigint_modulus);

		switch(bigint_value.bitfield.get_representation()) {
		case matrix_flags::dense_representation:
			if (A.bigint_value.bitfield.get_representation() ==
			    matrix_flags::dense_representation)
				DD_base_modul_bigint.insert_at(A.bigint_value, 0, 0, bigint_value,
							       0, 0, A.bigint_value.rows,
							       A.bigint_value.columns);
			else
				DS_base_modul_bigint.insert_at(A.bigint_value, 0, 0, bigint_value,
							       0, 0, A.bigint_value.rows,
							       A.bigint_value.columns);

			if (B.bigint_value.bitfield.get_representation() ==
			    matrix_flags::dense_representation)
				DD_base_modul_bigint.insert_at(B.bigint_value, 0, 0, bigint_value,
							       0, bigint_value.columns -
							       B.bigint_value.columns,
							       B.bigint_value.rows,
							       B.bigint_value.columns);
			else
				DS_base_modul_bigint.insert_at(B.bigint_value, 0, 0, bigint_value,
							       0, bigint_value.columns -
							       B.bigint_value.columns,
							       B.bigint_value.rows,
							       B.bigint_value.columns);
			break;
		case matrix_flags::sparse_representation:
			if (A.bigint_value.bitfield.get_representation() ==
			    matrix_flags::dense_representation)
				SD_base_modul_bigint.insert_at(A.bigint_value, 0, 0, bigint_value,
							       0, 0, A.bigint_value.rows,
							       A.bigint_value.columns);
			else
				SS_base_modul_bigint.insert_at(A.bigint_value, 0, 0, bigint_value,
							       0, 0, A.bigint_value.rows,
							       A.bigint_value.columns);

			if (B.bigint_value.bitfield.get_representation() ==
			    matrix_flags::dense_representation)
				SD_base_modul_bigint.insert_at(B.bigint_value, 0, 0, bigint_value,
							       0, bigint_value.columns -
							       B.bigint_value.columns,
							       B.bigint_value.rows,
							       B.bigint_value.columns);
			else
				SS_base_modul_bigint.insert_at(B.bigint_value, 0, 0, bigint_value,
							       0, bigint_value.columns -
							       B.bigint_value.columns,
							       B.bigint_value.rows,
							       B.bigint_value.columns);
			break;
		}
	}
	else {
		if (A.long_value.rows > long_value.rows
		    || A.long_value.columns > long_value.columns
		    || B.long_value.rows > long_value.rows
		    || B.long_value.columns > long_value.columns)
			precondition_error_handler(A.long_value.rows, "A.rows", "A.rows <= rows",
					    A.long_value.columns, "A.columns",
					    "A.columns <= coluns",
					    B.long_value.rows, "B.rows", "B.rows <= rows",
					    B.long_value.columns, "B.columns",
					    "B.columns <= coluns",
					    long_value.rows, "rows", "",
					    long_value.columns, "columns", "",
					    "void base_matrix< bigmod >::"
					    "split_h(base_matrix< bigmod > &A, "
					    "base_matrix< bigmod > &B) const",
					    DMESSAGE, EMESSAGE[4]);

		A.set_modulus(long_modulus);
		B.set_modulus(long_modulus);

		switch(long_value.bitfield.get_representation()) {
		case matrix_flags::dense_representation:
			if (A.long_value.bitfield.get_representation() ==
			    matrix_flags::dense_representation)
				DD_base_modul_long.insert_at(A.long_value, 0, 0, long_value,
							     0, 0, A.long_value.rows,
							     A.long_value.columns);
			else
				DS_base_modul_long.insert_at(A.long_value, 0, 0, long_value,
							     0, 0, A.long_value.rows,
							     A.long_value.columns);

			if (B.long_value.bitfield.get_representation() ==
			    matrix_flags::dense_representation)
				DD_base_modul_long.insert_at(B.long_value, 0, 0, long_value,
							     0, long_value.columns -
							     B.long_value.columns,
							     B.long_value.rows,
							     B.long_value.columns);
			else
				DS_base_modul_long.insert_at(B.long_value, 0, 0, long_value,
							     0, long_value.columns -
							     B.long_value.columns,
							     B.long_value.rows,
							     B.long_value.columns);
			break;
		case matrix_flags::sparse_representation:
			if (A.long_value.bitfield.get_representation() ==
			    matrix_flags::dense_representation)
				SD_base_modul_long.insert_at(A.long_value, 0, 0, long_value,
							     0, 0, A.long_value.rows,
							     A.long_value.columns);
			else
				SS_base_modul_long.insert_at(A.long_value, 0, 0, long_value,
							     0, 0, A.long_value.rows,
							     A.long_value.columns);

			if (B.long_value.bitfield.get_representation() ==
			    matrix_flags::dense_representation)
				SD_base_modul_long.insert_at(B.long_value, 0, 0, long_value,
							     0, long_value.columns -
							     B.long_value.columns,
							     B.long_value.rows,
							     B.long_value.columns);
			else
				SS_base_modul_long.insert_at(B.long_value, 0, 0, long_value,
							     0, long_value.columns -
							     B.long_value.columns,
							     B.long_value.rows,
							     B.long_value.columns);
			break;
		}
	}
}



void base_matrix< bigmod >::
split_h(bigmod *v, base_matrix< bigmod > &A) const
{
	//
	//       Task: C.split_h(v, A);
	//             => splits matrix C into the array v and the submatrix A
	//             => C = (v A)
	// Conditions: A.columns <= columns and
	//             A.rows <= rows and
	//             v != NULL
	//    Version: 2.0
	//

	debug_handler_l(DMESSAGE, "split_h(T *, base_matrix< bigmod > &)",
			DVALUE + 3);

	if (BigintLong) {
		if (A.bigint_value.rows > bigint_value.rows
		    || A.bigint_value.columns > bigint_value.columns || v == NULL)
			precondition_error_handler(A.bigint_value.rows, "A.rows", "A.rows <= rows",
					    A.bigint_value.columns, "A.columns",
					    "A.columns <= columns",
					    bigint_value.rows, "rows", "",
					    bigint_value.columns, "columns", "",
					    PRT, "v", "v != NULL",
					    "void base_matrix< bigmod >::"
					    "split_h(T *v, base_matrix< bigmod > &A) const",
					    DMESSAGE, EMESSAGE[4]);

		A.set_modulus(bigint_modulus);

		bigint *tmp = new bigint[bigint_value.rows];
		memory_handler(tmp, DMESSAGE, "split_h :: "
			       "Error in memory allocation (tmp)");

		switch(bigint_value.bitfield.get_representation()) {
		case matrix_flags::dense_representation:
			if (A.bigint_value.bitfield.get_representation() ==
			    matrix_flags::dense_representation)
				DD_base_modul_bigint.insert_at(A.bigint_value, 0, 0, bigint_value,
							       0, bigint_value.columns -
							       A.bigint_value.columns,
							       A.bigint_value.rows,
							       A.bigint_value.columns);
			else
				DS_base_modul_bigint.insert_at(A.bigint_value, 0, 0, bigint_value,
							       0, bigint_value.columns -
							       A.bigint_value.columns,
							       A.bigint_value.rows,
							       A.bigint_value.columns);
			D_base_modul_bigint.get_column(bigint_value, tmp, 0);
			break;
		case matrix_flags::sparse_representation:
			if (A.bigint_value.bitfield.get_representation() ==
			    matrix_flags::dense_representation)
				SD_base_modul_bigint.insert_at(A.bigint_value, 0, 0, bigint_value,
							       0, bigint_value.columns -
							       A.bigint_value.columns,
							       A.bigint_value.rows,
							       A.bigint_value.columns);
			else
				SS_base_modul_bigint.insert_at(A.bigint_value, 0, 0, bigint_value,
							       0, bigint_value.columns -
							       A.bigint_value.columns,
							       A.bigint_value.rows,
							       A.bigint_value.columns);
			S_base_modul_bigint.get_column(bigint_value, tmp, 0);
			break;
		}
		lidia_size_t i;
		for (i = 0; i < bigint_value.rows; i++) {
			v[i].set_modulus(bigint_modulus);
			v[i] = tmp[i];
		}
		delete[] tmp;
	}
	else {
		if (A.long_value.rows > long_value.rows
		    || A.long_value.columns > long_value.columns || v == NULL)
			precondition_error_handler(A.long_value.rows, "A.rows", "A.rows <= rows",
					    A.long_value.columns, "A.columns",
					    "A.columns <= columns",
					    long_value.rows, "rows", "",
					    long_value.columns, "columns", "",
					    PRT, "v", "v != NULL",
					    "void base_matrix< bigmod >::"
					    "split_h(T *v, base_matrix< bigmod > &A) const",
					    DMESSAGE, EMESSAGE[4]);

		A.set_modulus(long_modulus);

		long *tmp = new long[long_value.rows];
		memory_handler(tmp, DMESSAGE, "split_h :: "
			       "Error in memory allocation (tmp)");

		switch(long_value.bitfield.get_representation()) {
		case matrix_flags::dense_representation:
			if (A.long_value.bitfield.get_representation() ==
			    matrix_flags::dense_representation)
				DD_base_modul_long.insert_at(A.long_value, 0, 0, long_value,
							     0, long_value.columns -
							     A.long_value.columns,
							     A.long_value.rows,
							     A.long_value.columns);
			else
				DS_base_modul_long.insert_at(A.long_value, 0, 0, long_value,
							     0, long_value.columns -
							     A.long_value.columns,
							     A.long_value.rows,
							     A.long_value.columns);
			D_base_modul_long.get_column(long_value, tmp, 0);
			break;
		case matrix_flags::sparse_representation:
			if (A.long_value.bitfield.get_representation() ==
			    matrix_flags::dense_representation)
				SD_base_modul_long.insert_at(A.long_value, 0, 0, long_value,
							     0, long_value.columns -
							     A.long_value.columns,
							     A.long_value.rows,
							     A.long_value.columns);
			else
				SS_base_modul_long.insert_at(A.long_value, 0, 0, long_value,
							     0, long_value.columns -
							     A.long_value.columns,
							     A.long_value.rows,
							     A.long_value.columns);
			S_base_modul_long.get_column(long_value, tmp, 0);
			break;
		}
		lidia_size_t i;
		for (i = 0; i < long_value.rows; i++) {
			v[i].set_modulus(long_modulus);
			v[i] = tmp[i];
		}
		delete[] tmp;
	}
}



void base_matrix< bigmod >::
split_h(base_matrix< bigmod > &A, bigmod *v) const
{
	//
	//       Task: C.split_h(A, v);
	//             => splits matrix C into the submatrix A and the array v
	//             => C = (A v)
	// Conditions: A.columns <= columns and
	//             A.rows <= rows and
	//             v != NULL
	//    Version: 2.0
	//

	debug_handler_l(DMESSAGE, "split_h(base_matrix< bigmod > &, bigmod *)",
			DVALUE + 3);

	if (BigintLong) {
		if (A.bigint_value.rows > bigint_value.rows
		    || A.bigint_value.columns > bigint_value.columns || v == NULL)
			precondition_error_handler(A.bigint_value.rows, "A.rows", "A.rows <= rows",
					    A.bigint_value.columns, "A.columns",
					    "A.columns <= columns",
					    bigint_value.rows, "rows", "",
					    bigint_value.columns, "columns", "",
					    PRT, "v", "v != NULL",
					    "void base_matrix< bigmod >::"
					    "split_h(base_matrix< bigmod > &A, "
					    "bigmod *v) const", DMESSAGE, EMESSAGE[4]);

		A.set_modulus(bigint_modulus);

		bigint *tmp = new bigint[bigint_value.rows];
		memory_handler(tmp, DMESSAGE, "split_h :: "
			       "Error in memory allocation (tmp)");

		switch(bigint_value.bitfield.get_representation()) {
		case matrix_flags::dense_representation:
			if (A.bigint_value.bitfield.get_representation() ==
			    matrix_flags::dense_representation)
				DD_base_modul_bigint.insert_at(A.bigint_value, 0, 0, bigint_value,
							       0, 0, A.bigint_value.rows,
							       A.bigint_value.columns);
			else
				DS_base_modul_bigint.insert_at(A.bigint_value, 0, 0, bigint_value,
							       0, 0, A.bigint_value.rows,
							       A.bigint_value.columns);
			D_base_modul_bigint.get_column(bigint_value, tmp,
						       bigint_value.columns - 1);
			break;
		case matrix_flags::sparse_representation:
			if (A.bigint_value.bitfield.get_representation() ==
			    matrix_flags::dense_representation)
				SD_base_modul_bigint.insert_at(A.bigint_value, 0, 0, bigint_value,
							       0, 0, A.bigint_value.rows,
							       A.bigint_value.columns);
			else
				SS_base_modul_bigint.insert_at(A.bigint_value, 0, 0, bigint_value,
							       0, 0, A.bigint_value.rows,
							       A.bigint_value.columns);
			S_base_modul_bigint.get_column(bigint_value, tmp,
						       bigint_value.columns - 1);
			break;
		}
		lidia_size_t i;
		for (i = 0; i < bigint_value.rows; i++) {
			v[i].set_modulus(bigint_modulus);
			v[i] = tmp[i];
		}
		delete[] tmp;
	}
	else {
		if (A.long_value.rows > long_value.rows
		    || A.long_value.columns > long_value.columns || v == NULL)
			precondition_error_handler(A.long_value.rows, "A.rows", "A.rows <= rows",
					    A.long_value.columns, "A.columns",
					    "A.columns <= columns",
					    long_value.rows, "rows", "",
					    long_value.columns, "columns", "",
					    PRT, "v", "v != NULL",
					    "void base_matrix< bigmod >::"
					    "split_h(base_matrix< bigmod > &A, "
					    "bigmod *v) const", DMESSAGE, EMESSAGE[4]);

		A.set_modulus(long_modulus);

		long *tmp = new long[long_value.rows];
		memory_handler(tmp, DMESSAGE, "split_h :: "
			       "Error in memory allocation (tmp)");

		switch(long_value.bitfield.get_representation()) {
		case matrix_flags::dense_representation:
			if (A.long_value.bitfield.get_representation() ==
			    matrix_flags::dense_representation)
				DD_base_modul_long.insert_at(A.long_value, 0, 0, long_value,
							     0, 0, A.long_value.rows,
							     A.long_value.columns);
			else
				DS_base_modul_long.insert_at(A.long_value, 0, 0, long_value,
							     0, 0, A.long_value.rows,
							     A.long_value.columns);
			D_base_modul_long.get_column(long_value, tmp,
						     long_value.columns - 1);
			break;
		case matrix_flags::sparse_representation:
			if (A.long_value.bitfield.get_representation() ==
			    matrix_flags::dense_representation)
				SD_base_modul_long.insert_at(A.long_value, 0, 0, long_value,
							     0, 0, A.long_value.rows,
							     A.long_value.columns);
			else
				SS_base_modul_long.insert_at(A.long_value, 0, 0, long_value,
							     0, 0, A.long_value.rows,
							     A.long_value.columns);
			S_base_modul_long.get_column(long_value, tmp,
						     long_value.columns - 1);
			break;
		}
		lidia_size_t i;
		for (i = 0; i < long_value.rows; i++) {
			v[i].set_modulus(long_modulus);
			v[i] = tmp[i];
		}
		delete[] tmp;
	}
}



void base_matrix< bigmod >::
split_h(base_vector< bigmod > &v, base_matrix< bigmod > &A) const
{
	//
	//       Task: C.split_h(v, A);
	//             => splits matrix C into the vector v and the submatrix A
	//             => C = (v A)
	// Conditions: v.size <= rows and
	//             A.columns <= columns and
	//             A.rows <= rows
	//    Version: 2.0
	//

	debug_handler_l(DMESSAGE, "split_h(base_vector< bigmod > &, "
			"base_matrix< bigmod > &)", DVALUE + 3);

	lidia_size_t l = v.size();

	if (BigintLong) {
		if (A.bigint_value.rows > bigint_value.rows
		    || A.bigint_value.columns > bigint_value.columns
		    || l > bigint_value.rows)
			precondition_error_handler(A.bigint_value.rows, "A.rows", "A.rows <= rows",
					    A.bigint_value.columns, "A.columns",
					    "A.columns <= columns",
					    bigint_value.rows, "rows", "",
					    bigint_value.columns, "columns", "",
					    l, "v.size", "v.size <= rows",
					    "void base_matrix< bigmod >::"
					    "split_h(base_vector< bigmod > &v, "
					    "base_matrix< bigmod > &A) const", DMESSAGE,
					    EMESSAGE[4]);

		A.set_modulus(bigint_modulus);

		bigint *tmp = new bigint[bigint_value.rows];
		memory_handler(tmp, DMESSAGE, "split_h :: "
			       "Error in memory allocation (tmp)");

		switch(bigint_value.bitfield.get_representation()) {
		case matrix_flags::dense_representation:
			if (A.bigint_value.bitfield.get_representation() ==
			    matrix_flags::dense_representation)
				DD_base_modul_bigint.insert_at(A.bigint_value, 0, 0, bigint_value,
							       0, bigint_value.columns -
							       A.bigint_value.columns,
							       A.bigint_value.rows,
							       A.bigint_value.columns);
			else
				DS_base_modul_bigint.insert_at(A.bigint_value, 0, 0, bigint_value,
							       0, bigint_value.columns -
							       A.bigint_value.columns,
							       A.bigint_value.rows,
							       A.bigint_value.columns);
			D_base_modul_bigint.get_column(bigint_value, tmp, 0);
			break;
		case matrix_flags::sparse_representation:
			if (A.bigint_value.bitfield.get_representation() ==
			    matrix_flags::dense_representation)
				SD_base_modul_bigint.insert_at(A.bigint_value, 0, 0, bigint_value,
							       0, bigint_value.columns -
							       A.bigint_value.columns,
							       A.bigint_value.rows,
							       A.bigint_value.columns);
			else
				SS_base_modul_bigint.insert_at(A.bigint_value, 0, 0, bigint_value,
							       0, bigint_value.columns -
							       A.bigint_value.columns,
							       A.bigint_value.rows,
							       A.bigint_value.columns);
			S_base_modul_bigint.get_column(bigint_value, tmp, 0);
			break;
		}
		lidia_size_t i;
		for (i = 0; i < bigint_value.rows; i++) {
			v[i].set_modulus(bigint_modulus);
			v[i] = tmp[i];
		}
		delete[] tmp;
	}
	else {
		if (A.long_value.rows > long_value.rows
		    || A.long_value.columns > long_value.columns
		    || l > long_value.rows)
			precondition_error_handler(A.long_value.rows, "A.rows", "A.rows <= rows",
					    A.long_value.columns, "A.columns",
					    "A.columns <= columns",
					    long_value.rows, "rows", "",
					    long_value.columns, "columns", "",
					    l, "v.size", "v.size <= rows",
					    "void base_matrix< bigmod >::"
					    "split_h(base_vector< bigmod > &v, "
					    "base_matrix< bigmod > &A) const", DMESSAGE,
					    EMESSAGE[4]);

		A.set_modulus(long_modulus);

		long *tmp = new long[long_value.rows];
		memory_handler(tmp, DMESSAGE, "split_h :: "
			       "Error in memory allocation (tmp)");

		switch(long_value.bitfield.get_representation()) {
		case matrix_flags::dense_representation:
			if (A.long_value.bitfield.get_representation() ==
			    matrix_flags::dense_representation)
				DD_base_modul_long.insert_at(A.long_value, 0, 0, long_value,
							     0, long_value.columns -
							     A.long_value.columns,
							     A.long_value.rows,
							     A.long_value.columns);
			else
				DS_base_modul_long.insert_at(A.long_value, 0, 0, long_value,
							     0, long_value.columns -
							     A.long_value.columns,
							     A.long_value.rows,
							     A.long_value.columns);
			D_base_modul_long.get_column(long_value, tmp, 0);
			break;
		case matrix_flags::sparse_representation:
			if (A.long_value.bitfield.get_representation() ==
			    matrix_flags::dense_representation)
				SD_base_modul_long.insert_at(A.long_value, 0, 0, long_value,
							     0, long_value.columns -
							     A.long_value.columns,
							     A.long_value.rows,
							     A.long_value.columns);
			else
				SS_base_modul_long.insert_at(A.long_value, 0, 0, long_value,
							     0, long_value.columns -
							     A.long_value.columns,
							     A.long_value.rows,
							     A.long_value.columns);
			S_base_modul_long.get_column(long_value, tmp, 0);
			break;
		}
		lidia_size_t i;
		for (i = 0; i < long_value.rows; i++) {
			v[i].set_modulus(long_modulus);
			v[i] = tmp[i];
		}
		delete[] tmp;
	}
}



void base_matrix< bigmod >::
split_h(base_matrix< bigmod > &A, base_vector< bigmod > &v) const
{
	//
	//       Task: C.split_h(A, v);
	//             => splits matrix C into the submatrix A and the vector v
	//             => C = (A v)
	// Conditions: A.columns <= columns and
	//             A.rows <= rows and
	//             v.size <= rows
	//    Version: 2.0
	//

	debug_handler_l(DMESSAGE, "split_h(base_matrix< bigmod > &, "
			"base_vector< bigmod > &)", DVALUE + 3);

	lidia_size_t l = v.size();

	if (BigintLong) {
		if (A.bigint_value.rows > bigint_value.rows
		    || A.bigint_value.columns > bigint_value.columns
		    || l > bigint_value.rows)
			precondition_error_handler(A.bigint_value.rows, "A.rows", "A.rows <= rows",
					    A.bigint_value.columns, "A.columns",
					    "A.columns <= columns",
					    bigint_value.rows, "rows", "",
					    bigint_value.columns, "columns", "",
					    l, "v.size", "v.size <= rows",
					    "void base_matrix< bigmod >::"
					    "split_h(base_matrix< bigmod > &A, "
					    "base_vector< bigmod > &v) const",
					    DMESSAGE, EMESSAGE[4]);

		A.set_modulus(bigint_modulus);

		bigint *tmp = new bigint[bigint_value.rows];
		memory_handler(tmp, DMESSAGE, "split_h :: "
			       "Error in memory allocation (tmp)");

		switch(bigint_value.bitfield.get_representation()) {
		case matrix_flags::dense_representation:
			if (A.bigint_value.bitfield.get_representation() ==
			    matrix_flags::dense_representation)
				DD_base_modul_bigint.insert_at(A.bigint_value, 0, 0, bigint_value,
							       0, 0, A.bigint_value.rows,
							       A.bigint_value.columns);
			else
				DS_base_modul_bigint.insert_at(A.bigint_value, 0, 0, bigint_value,
							       0, 0, A.bigint_value.rows,
							       A.bigint_value.columns);
			D_base_modul_bigint.get_column(bigint_value, tmp,
						       bigint_value.columns - 1);
			break;
		case matrix_flags::sparse_representation:
			if (A.bigint_value.bitfield.get_representation() ==
			    matrix_flags::dense_representation)
				SD_base_modul_bigint.insert_at(A.bigint_value, 0, 0, bigint_value,
							       0, 0, A.bigint_value.rows,
							       A.bigint_value.columns);
			else
				SS_base_modul_bigint.insert_at(A.bigint_value, 0, 0, bigint_value,
							       0, 0, A.bigint_value.rows,
							       A.bigint_value.columns);
			S_base_modul_bigint.get_column(bigint_value, tmp,
						       bigint_value.columns - 1);
			break;
		}
		lidia_size_t i;
		for (i = 0; i < bigint_value.rows; i++) {
			v[i].set_modulus(bigint_modulus);
			v[i] = tmp[i];
		}
		delete[] tmp;
	}
	else {
		if (A.long_value.rows > long_value.rows
		    || A.long_value.columns > long_value.columns
		    || l > long_value.rows)
			precondition_error_handler(A.long_value.rows, "A.rows", "A.rows <= rows",
					    A.long_value.columns, "A.columns",
					    "A.columns <= columns",
					    long_value.rows, "rows", "",
					    long_value.columns, "columns", "",
					    l, "v.size", "v.size <= rows",
					    "void base_matrix< bigmod >::"
					    "split_h(base_matrix< bigmod > &A, "
					    "base_vector< bigmod > &v) const",
					    DMESSAGE, EMESSAGE[4]);

		A.set_modulus(long_modulus);

		long *tmp = new long[long_value.rows];
		memory_handler(tmp, DMESSAGE, "split_h :: "
			       "Error in memory allocation (tmp)");

		switch(long_value.bitfield.get_representation()) {
		case matrix_flags::dense_representation:
			if (A.long_value.bitfield.get_representation() ==
			    matrix_flags::dense_representation)
				DD_base_modul_long.insert_at(A.long_value, 0, 0, long_value,
							     0, 0, A.long_value.rows,
							     A.long_value.columns);
			else
				DS_base_modul_long.insert_at(A.long_value, 0, 0, long_value,
							     0, 0, A.long_value.rows,
							     A.long_value.columns);
			D_base_modul_long.get_column(long_value, tmp,
						     long_value.columns - 1);
			break;
		case matrix_flags::sparse_representation:
			if (A.long_value.bitfield.get_representation() ==
			    matrix_flags::dense_representation)
				SD_base_modul_long.insert_at(A.long_value, 0, 0, long_value,
							     0, 0, A.long_value.rows,
							     A.long_value.columns);
			else
				SS_base_modul_long.insert_at(A.long_value, 0, 0, long_value,
							     0, 0, A.long_value.rows,
							     A.long_value.columns);
			S_base_modul_long.get_column(long_value, tmp,
						     long_value.columns - 1);
			break;
		}
		lidia_size_t i;
		for (i = 0; i < long_value.rows; i++) {
			v[i].set_modulus(long_modulus);
			v[i] = tmp[i];
		}
		delete[] tmp;
	}
}



void base_matrix< bigmod >::
split_v(base_matrix< bigmod > &A, base_matrix< bigmod > &B) const
{
	//
	//       Task: C.split_v(A,B);
	//             => splits matrix C into the submatrices A,B
	//                    (A)
	//             => C = (B)
	// Conditions: A.columns <= columns and
	//             A.rows <= rows and
	//             B.columns <= columns and
	//             B.rows <= rows
	//    Version: 2.0
	//

	debug_handler_l(DMESSAGE, "split_v(base_matrix< bigmod > &, "
			"base_matrix< bigmod > &)", DVALUE + 3);

	if (BigintLong) {
		if (A.bigint_value.rows > bigint_value.rows
		    || A.bigint_value.columns > bigint_value.columns
		    || B.bigint_value.rows > bigint_value.rows
		    || B.bigint_value.columns > bigint_value.columns)
			precondition_error_handler(A.bigint_value.rows, "A.rows", "A.rows <= rows",
					    A.bigint_value.columns, "A.columns",
					    "A.columns <= columns",
					    B.bigint_value.rows, "B.rows", "B.rows <= rows",
					    B.bigint_value.columns, "B.columns",
					    "B.columns <= columns",
					    bigint_value.rows, "rows", "",
					    bigint_value.columns, "columns", "",
					    "void base_matrix< bigmod >::"
					    "split_v(base_matrix< bigmod > &A, "
					    "base_matrix< bigmod > &B) const",
					    DMESSAGE, EMESSAGE[4]);

		A.set_modulus(bigint_modulus);
		B.set_modulus(bigint_modulus);

		switch(bigint_value.bitfield.get_representation()) {
		case matrix_flags::dense_representation:
			if (A.bigint_value.bitfield.get_representation() ==
			    matrix_flags::dense_representation)
				DD_base_modul_bigint.insert_at(A.bigint_value, 0, 0, bigint_value,
							       0, 0, A.bigint_value.rows,
							       A.bigint_value.columns);
			else
				DS_base_modul_bigint.insert_at(A.bigint_value, 0, 0, bigint_value,
							       0, 0, A.bigint_value.rows,
							       A.bigint_value.columns);

			if (B.bigint_value.bitfield.get_representation() ==
			    matrix_flags::dense_representation)
				DD_base_modul_bigint.insert_at(B.bigint_value, 0, 0, bigint_value,
							       bigint_value.rows -
							       B.bigint_value.rows, 0,
							       B.bigint_value.rows,
							       B.bigint_value.columns);
			else
				DS_base_modul_bigint.insert_at(B.bigint_value, 0, 0, bigint_value,
							       bigint_value.rows -
							       B.bigint_value.rows, 0,
							       B.bigint_value.rows,
							       B.bigint_value.columns);
			break;
		case matrix_flags::sparse_representation:
			if (A.bigint_value.bitfield.get_representation() ==
			    matrix_flags::dense_representation)
				SD_base_modul_bigint.insert_at(A.bigint_value, 0, 0, bigint_value,
							       0, 0, A.bigint_value.rows,
							       A.bigint_value.columns);
			else
				SS_base_modul_bigint.insert_at(A.bigint_value, 0, 0, bigint_value,
							       0, 0, A.bigint_value.rows,
							       A.bigint_value.columns);

			if (B.bigint_value.bitfield.get_representation() ==
			    matrix_flags::dense_representation)
				SD_base_modul_bigint.insert_at(B.bigint_value, 0, 0, bigint_value,
							       0, bigint_value.columns -
							       B.bigint_value.columns,
							       B.bigint_value.rows,
							       B.bigint_value.columns);
			else
				SS_base_modul_bigint.insert_at(B.bigint_value, 0, 0, bigint_value,
							       0, bigint_value.columns -
							       B.bigint_value.columns,
							       B.bigint_value.rows,
							       B.bigint_value.columns);
			break;
		}
	}
	else {
		if (A.long_value.rows > long_value.rows
		    || A.long_value.columns > long_value.columns
		    || B.long_value.rows > long_value.rows
		    || B.long_value.columns > long_value.columns)
			precondition_error_handler(A.long_value.rows, "A.rows", "A.rows <= rows",
					    A.long_value.columns, "A.columns",
					    "A.columns <= columns",
					    B.long_value.rows, "B.rows", "B.rows <= rows",
					    B.long_value.columns, "B.columns",
					    "B.columns <= columns",
					    long_value.rows, "rows", "",
					    long_value.columns, "columns", "",
					    "void base_matrix< bigmod >::"
					    "split_v(base_matrix< bigmod > &A, "
					    "base_matrix< bigmod > &B) const",
					    DMESSAGE, EMESSAGE[4]);

		A.set_modulus(long_modulus);
		B.set_modulus(long_modulus);

		switch(long_value.bitfield.get_representation()) {
		case matrix_flags::dense_representation:
			if (A.long_value.bitfield.get_representation() ==
			    matrix_flags::dense_representation)
				DD_base_modul_long.insert_at(A.long_value, 0, 0, long_value,
							     0, 0, A.long_value.rows,
							     A.long_value.columns);
			else
				DS_base_modul_long.insert_at(A.long_value, 0, 0, long_value,
							     0, 0, A.long_value.rows,
							     A.long_value.columns);

			if (B.long_value.bitfield.get_representation() ==
			    matrix_flags::dense_representation)
				DD_base_modul_long.insert_at(B.long_value, 0, 0, long_value,
							     long_value.rows -
							     B.long_value.rows, 0,
							     B.long_value.rows,
							     B.long_value.columns);
			else
				DS_base_modul_long.insert_at(B.long_value, 0, 0, long_value,
							     long_value.rows -
							     B.long_value.rows, 0,
							     B.long_value.rows,
							     B.long_value.columns);
			break;
		case matrix_flags::sparse_representation:
			if (A.long_value.bitfield.get_representation() ==
			    matrix_flags::dense_representation)
				SD_base_modul_long.insert_at(A.long_value, 0, 0, long_value,
							     0, 0, A.long_value.rows,
							     A.long_value.columns);
			else
				SS_base_modul_long.insert_at(A.long_value, 0, 0, long_value,
							     0, 0, A.long_value.rows,
							     A.long_value.columns);

			if (B.long_value.bitfield.get_representation() ==
			    matrix_flags::dense_representation)
				SD_base_modul_long.insert_at(B.long_value, 0, 0, long_value,
							     0, long_value.columns -
							     B.long_value.columns,
							     B.long_value.rows,
							     B.long_value.columns);
			else
				SS_base_modul_long.insert_at(B.long_value, 0, 0, long_value,
							     0, long_value.columns -
							     B.long_value.columns,
							     B.long_value.rows,
							     B.long_value.columns);
			break;
		}
	}
}



void base_matrix< bigmod >::
split_v(bigmod *v, base_matrix< bigmod > &A) const
{
	//
	//       Task: C.split_v(v, A);
	//             => splits matrix C into the array v and the submatrix A
	//                    (v)
	//             => C = (A)
	// Conditions: A.columns <= columns and
	//             A.rows <= rows and
	//             v != NULL
	//    Version: 2.0
	//

	debug_handler_l(DMESSAGE, "split_v(bigmod *, base_matrix< bigmod > &)",
			DVALUE + 3);

	if (BigintLong) {
		if (A.bigint_value.rows > bigint_value.rows
		    || A.bigint_value.columns > bigint_value.columns || v == NULL)
			precondition_error_handler(A.bigint_value.rows, "A.rows", "A.rows <= rows",
					    A.bigint_value.columns, "A.columns",
					    "A.columns <= columns",
					    bigint_value.rows, "rows", "",
					    bigint_value.columns, "columns", "",
					    PRT, "v", "v != NULL",
					    "void base_matrix< bigmod >::"
					    "split_v(bigmod *v, "
					    "base_matrix< bigmod > &A) const",
					    DMESSAGE, EMESSAGE[4]);

		A.set_modulus(bigint_modulus);

		bigint *tmp = new bigint[bigint_value.columns];
		memory_handler(tmp, DMESSAGE, "split_v :: "
			       "Error in memory allocation (tmp)");

		switch(bigint_value.bitfield.get_representation()) {
		case matrix_flags::dense_representation:
			if (A.bigint_value.bitfield.get_representation() ==
			    matrix_flags::dense_representation)
				DD_base_modul_bigint.insert_at(A.bigint_value, 0, 0, bigint_value,
							       bigint_value.rows -
							       A.bigint_value.rows, 0,
							       A.bigint_value.rows,
							       A.bigint_value.columns);
			else
				DS_base_modul_bigint.insert_at(A.bigint_value, 0, 0, bigint_value,
							       bigint_value.rows -
							       A.bigint_value.rows, 0,
							       A.bigint_value.rows,
							       A.bigint_value.columns);
			D_base_modul_bigint.get_row(bigint_value, tmp, 0);
			break;
		case matrix_flags::sparse_representation:
			if (A.bigint_value.bitfield.get_representation() ==
			    matrix_flags::dense_representation)
				SD_base_modul_bigint.insert_at(A.bigint_value, 0, 0, bigint_value,
							       bigint_value.rows -
							       A.bigint_value.rows, 0,
							       A.bigint_value.rows,
							       A.bigint_value.columns);
			else
				SS_base_modul_bigint.insert_at(A.bigint_value, 0, 0, bigint_value,
							       bigint_value.rows -
							       A.bigint_value.rows, 0,
							       A.bigint_value.rows,
							       A.bigint_value.columns);
			S_base_modul_bigint.get_row(bigint_value, tmp, 0);
			break;
		}
		lidia_size_t i;
		for (i = 0; i < bigint_value.columns; i++) {
			v[i].set_modulus(bigint_modulus);
			v[i] = tmp[i];
		}
		delete[] tmp;
	}
	else {
		if (A.long_value.rows > long_value.rows
		    || A.long_value.columns > long_value.columns || v == NULL)
			precondition_error_handler(A.long_value.rows, "A.rows", "A.rows <= rows",
					    A.long_value.columns, "A.columns",
					    "A.columns <= columns",
					    long_value.rows, "rows", "",
					    long_value.columns, "columns", "",
					    PRT, "v", "v != NULL",
					    "void base_matrix< bigmod >::"
					    "split_v(bigmod *v, "
					    "base_matrix< bigmod > &A) const",
					    DMESSAGE, EMESSAGE[4]);

		A.set_modulus(long_modulus);

		long *tmp = new long[long_value.columns];
		memory_handler(tmp, DMESSAGE, "split_v :: "
			       "Error in memory allocation (tmp)");

		switch(long_value.bitfield.get_representation()) {
		case matrix_flags::dense_representation:
			if (A.long_value.bitfield.get_representation() ==
			    matrix_flags::dense_representation)
				DD_base_modul_long.insert_at(A.long_value, 0, 0, long_value,
							     long_value.rows -
							     A.long_value.rows, 0,
							     A.long_value.rows,
							     A.long_value.columns);
			else
				DS_base_modul_long.insert_at(A.long_value, 0, 0, long_value,
							     long_value.rows -
							     A.long_value.rows, 0,
							     A.long_value.rows,
							     A.long_value.columns);
			D_base_modul_long.get_row(long_value, tmp, 0);
			break;
		case matrix_flags::sparse_representation:
			if (A.long_value.bitfield.get_representation() ==
			    matrix_flags::dense_representation)
				SD_base_modul_long.insert_at(A.long_value, 0, 0, long_value,
							     long_value.rows -
							     A.long_value.rows, 0,
							     A.long_value.rows,
							     A.long_value.columns);
			else
				SS_base_modul_long.insert_at(A.long_value, 0, 0, long_value,
							     long_value.rows -
							     A.long_value.rows, 0,
							     A.long_value.rows,
							     A.long_value.columns);
			S_base_modul_long.get_row(long_value, tmp, 0);
			break;
		}
		lidia_size_t i;
		for (i = 0; i < long_value.columns; i++) {
			v[i].set_modulus(long_modulus);
			v[i] = tmp[i];
		}
		delete[] tmp;
	}
}



void base_matrix< bigmod >::
split_v(base_matrix< bigmod > &A, bigmod *v) const
{
	//
	//       Task: C.split_v(A,v);
	//             => splits matrix C into the submatrix A and the vector v
	//                    (A)
	//             => C = (v)
	// Conditions: A.columns <= columns and
	//             A.rows <= rows and
	//             v != NULL
	//    Version: 2.0
	//

	debug_handler_l(DMESSAGE, "split_v(base_matrix< bigmod > &, bigmod *)",
			DVALUE + 3);

	if (BigintLong) {
		if (A.bigint_value.rows > bigint_value.rows
		    || A.bigint_value.columns > bigint_value.columns || v == NULL)
			precondition_error_handler(A.bigint_value.rows, "A.rows", "A.rows <= rows",
					    A.bigint_value.columns, "A.columns",
					    "A.columsn <= columns",
					    bigint_value.rows, "rows", "",
					    bigint_value.columns, "columns", "",
					    PRT, "v", "v != NULL",
					    "void base_matrix< bigmod >::"
					    "split_v(base_matrix< bigmod > &A, "
					    "bigmod *v) const", DMESSAGE, EMESSAGE[4]);

		A.set_modulus(bigint_modulus);

		bigint *tmp = new bigint[bigint_value.columns];
		memory_handler(tmp, DMESSAGE, "split_v :: "
			       "Error in memory allocation (tmp)");

		switch(bigint_value.bitfield.get_representation()) {
		case matrix_flags::dense_representation:
			if (A.bigint_value.bitfield.get_representation() ==
			    matrix_flags::dense_representation)
				DD_base_modul_bigint.insert_at(A.bigint_value, 0, 0, bigint_value,
							       0, 0, A.bigint_value.rows,
							       A.bigint_value.columns);
			else
				DS_base_modul_bigint.insert_at(A.bigint_value, 0, 0, bigint_value,
							       0, 0, A.bigint_value.rows,
							       A.bigint_value.columns);
			D_base_modul_bigint.get_row(bigint_value, tmp,
						    bigint_value.rows - 1);
			break;
		case matrix_flags::sparse_representation:
			if (A.bigint_value.bitfield.get_representation() ==
			    matrix_flags::dense_representation)
				SD_base_modul_bigint.insert_at(A.bigint_value, 0, 0, bigint_value,
							       0, 0, A.bigint_value.rows,
							       A.bigint_value.columns);
			else
				SS_base_modul_bigint.insert_at(A.bigint_value, 0, 0, bigint_value,
							       0, 0, A.bigint_value.rows,
							       A.bigint_value.columns);
			S_base_modul_bigint.get_row(bigint_value, tmp,
						    bigint_value.rows - 1);
			break;
		}
		lidia_size_t i;
		for (i = 0; i < bigint_value.columns; i++) {
			v[i].set_modulus(bigint_modulus);
			v[i] = tmp[i];
		}
		delete[] tmp;
	}
	else {
		if (A.long_value.rows > long_value.rows
		    || A.long_value.columns > long_value.columns || v == NULL)
			precondition_error_handler(A.long_value.rows, "A.rows", "A.rows <= rows",
					    A.long_value.columns, "A.columns",
					    "A.columsn <= columns",
					    long_value.rows, "rows", "",
					    long_value.columns, "columns", "",
					    PRT, "v", "v != NULL",
					    "void base_matrix< bigmod >::"
					    "split_v(base_matrix< bigmod > &A, "
					    "bigmod *v) const", DMESSAGE, EMESSAGE[4]);

		A.set_modulus(long_modulus);

		long *tmp = new long[long_value.columns];
		memory_handler(tmp, DMESSAGE, "split_v :: "
			       "Error in memory allocation (tmp)");

		switch(long_value.bitfield.get_representation()) {
		case matrix_flags::dense_representation:
			if (A.long_value.bitfield.get_representation() ==
			    matrix_flags::dense_representation)
				DD_base_modul_long.insert_at(A.long_value, 0, 0, long_value,
							     0, 0, A.long_value.rows,
							     A.long_value.columns);
			else
				DS_base_modul_long.insert_at(A.long_value, 0, 0, long_value,
							     0, 0, A.long_value.rows,
							     A.long_value.columns);
			D_base_modul_long.get_row(long_value, tmp,
						  long_value.rows - 1);
			break;
		case matrix_flags::sparse_representation:
			if (A.long_value.bitfield.get_representation() ==
			    matrix_flags::dense_representation)
				SD_base_modul_long.insert_at(A.long_value, 0, 0, long_value,
							     0, 0, A.long_value.rows,
							     A.long_value.columns);
			else
				SS_base_modul_long.insert_at(A.long_value, 0, 0, long_value,
							     0, 0, A.long_value.rows,
							     A.long_value.columns);
			S_base_modul_long.get_row(long_value, tmp,
						  long_value.rows - 1);
			break;
		}
		lidia_size_t i;
		for (i = 0; i < long_value.columns; i++) {
			v[i].set_modulus(long_modulus);
			v[i] = tmp[i];
		}
		delete[] tmp;
	}
}



void base_matrix< bigmod >::
split_v(base_vector< bigmod > &v, base_matrix< bigmod > &A) const
{
	//
	//       Task: C.split_v(v, A);
	//             => splits matrix C into the vector v and the submatrix A
	//                    (v)
	//             => C = (A)
	// Conditions: A.columns <= columns and
	//             A.rows <= rows and
	//             v.size <= columns
	//    Version: 2.0
	//

	debug_handler_l(DMESSAGE, "split_v(base_vector< bigmod > &, "
			"base_matrix< bigmod > &)", DVALUE + 3);

	lidia_size_t l = v.size();

	if (BigintLong) {
		if (A.bigint_value.rows > bigint_value.rows
		    || A.bigint_value.columns > bigint_value.columns
		    || l > bigint_value.columns)
			precondition_error_handler(A.bigint_value.rows, "A.rows", "A.rows <= rows",
					    A.bigint_value.columns, "A.columns",
					    "A.columns <= columns",
					    bigint_value.rows, "rows", "",
					    bigint_value.columns, "columns", "",
					    l, "v.size", "v.size <= columns",
					    "void base_matrix< bigmod >::"
					    "split_v(base_vector< bigmod > &v, "
					    "base_matrix< bigmod > &A) const",
					    DMESSAGE, EMESSAGE[4]);

		A.set_modulus(bigint_modulus);

		bigint *tmp = new bigint[bigint_value.columns];
		memory_handler(tmp, DMESSAGE, "split_v :: "
			       "Error in memory allocation (tmp)");

		switch(bigint_value.bitfield.get_representation()) {
		case matrix_flags::dense_representation:
			if (A.bigint_value.bitfield.get_representation() ==
			    matrix_flags::dense_representation)
				DD_base_modul_bigint.insert_at(A.bigint_value, 0, 0, bigint_value,
							       bigint_value.rows -
							       A.bigint_value.rows, 0,
							       A.bigint_value.rows,
							       A.bigint_value.columns);
			else
				DS_base_modul_bigint.insert_at(A.bigint_value, 0, 0, bigint_value,
							       bigint_value.rows -
							       A.bigint_value.rows, 0,
							       A.bigint_value.rows,
							       A.bigint_value.columns);
			D_base_modul_bigint.get_row(bigint_value, tmp, 0);
			break;
		case matrix_flags::sparse_representation:
			if (A.bigint_value.bitfield.get_representation() ==
			    matrix_flags::dense_representation)
				SD_base_modul_bigint.insert_at(A.bigint_value, 0, 0, bigint_value,
							       bigint_value.rows -
							       A.bigint_value.rows, 0,
							       A.bigint_value.rows,
							       A.bigint_value.columns);
			else
				SS_base_modul_bigint.insert_at(A.bigint_value, 0, 0, bigint_value,
							       bigint_value.rows -
							       A.bigint_value.rows, 0,
							       A.bigint_value.rows,
							       A.bigint_value.columns);
			S_base_modul_bigint.get_row(bigint_value, tmp, 0);
			break;
		}
		lidia_size_t i;
		for (i = 0; i < bigint_value.columns; i++) {
			v[i].set_modulus(bigint_modulus);
			v[i] = tmp[i];
		}
		delete[] tmp;
	}
	else {
		if (A.long_value.rows > long_value.rows
		    || A.long_value.columns > long_value.columns
		    || l > long_value.columns)
			precondition_error_handler(A.long_value.rows, "A.rows", "A.rows <= rows",
					    A.long_value.columns, "A.columns",
					    "A.columns <= columns",
					    long_value.rows, "rows", "",
					    long_value.columns, "columns", "",
					    l, "v.size", "v.size <= columns",
					    "void base_matrix< bigmod >::"
					    "split_v(base_vector< bigmod > &v, "
					    "base_matrix< bigmod > &A) const",
					    DMESSAGE, EMESSAGE[4]);

		A.set_modulus(long_modulus);

		long *tmp = new long[long_value.columns];
		memory_handler(tmp, DMESSAGE, "split_v :: "
			       "Error in memory allocation (tmp)");

		switch(long_value.bitfield.get_representation()) {
		case matrix_flags::dense_representation:
			if (A.long_value.bitfield.get_representation() ==
			    matrix_flags::dense_representation)
				DD_base_modul_long.insert_at(A.long_value, 0, 0, long_value,
							     long_value.rows -
							     A.long_value.rows, 0,
							     A.long_value.rows,
							     A.long_value.columns);
			else
				DS_base_modul_long.insert_at(A.long_value, 0, 0, long_value,
							     long_value.rows -
							     A.long_value.rows, 0,
							     A.long_value.rows,
							     A.long_value.columns);
			D_base_modul_long.get_row(long_value, tmp, 0);
			break;
		case matrix_flags::sparse_representation:
			if (A.long_value.bitfield.get_representation() ==
			    matrix_flags::dense_representation)
				SD_base_modul_long.insert_at(A.long_value, 0, 0, long_value,
							     long_value.rows -
							     A.long_value.rows, 0,
							     A.long_value.rows,
							     A.long_value.columns);
			else
				SS_base_modul_long.insert_at(A.long_value, 0, 0, long_value,
							     long_value.rows -
							     A.long_value.rows, 0,
							     A.long_value.rows,
							     A.long_value.columns);
			S_base_modul_long.get_row(long_value, tmp, 0);
			break;
		}
		lidia_size_t i;
		for (i = 0; i < long_value.columns; i++) {
			v[i].set_modulus(long_modulus);
			v[i] = tmp[i];
		}
		delete[] tmp;
	}
}



void base_matrix< bigmod >::
split_v(base_matrix< bigmod > &A, base_vector< bigmod > &v) const
{
	//
	//       Task: C.split_v(A,B);
	//             => splits matrix C into the submatrix A and the vector v
	//                    (A)
	//             => C = (v)
	// Conditions: A.columns <= columns and
	//             A.rows <= rows and
	//             v.size <= columns
	//    Version: 2.0
	//

	debug_handler_l(DMESSAGE, "split_v(base_matrix< bigmod > &, "
			"base_vector< bigmod > &)", DVALUE + 3);

	lidia_size_t l = v.size();

	if (BigintLong) {
		if (A.bigint_value.rows > bigint_value.rows
		    || A.bigint_value.columns > bigint_value.columns
		    || l > bigint_value.columns)
			precondition_error_handler(A.bigint_value.rows, "A.rows", "A.rows <= rows",
					    A.bigint_value.columns, "A.columns",
					    "A.columns <= columns",
					    bigint_value.rows, "rows", "",
					    bigint_value.columns, "columns", "",
					    l, "v.size", "v.size <= columns",
					    "void base_matrix< bigmod >::"
					    "split_v(base_matrix< bigmod > &A, "
					    "base_vector< bigmod > &v) const",
					    DMESSAGE, EMESSAGE[4]);

		A.set_modulus(bigint_modulus);

		bigint *tmp = new bigint[bigint_value.columns];
		memory_handler(tmp, DMESSAGE, "split_v :: "
			       "Error in memory allocation (tmp)");

		switch(bigint_value.bitfield.get_representation()) {
		case matrix_flags::dense_representation:
			if (A.bigint_value.bitfield.get_representation() ==
			    matrix_flags::dense_representation)
				DD_base_modul_bigint.insert_at(A.bigint_value, 0, 0, bigint_value,
							       0, 0, A.bigint_value.rows,
							       A.bigint_value.columns);
			else
				DS_base_modul_bigint.insert_at(A.bigint_value, 0, 0, bigint_value,
							       0, 0, A.bigint_value.rows,
							       A.bigint_value.columns);
			D_base_modul_bigint.get_row(bigint_value, tmp,
						    bigint_value.rows - 1);
			break;
		case matrix_flags::sparse_representation:
			if (A.bigint_value.bitfield.get_representation() ==
			    matrix_flags::dense_representation)
				SD_base_modul_bigint.insert_at(A.bigint_value, 0, 0, bigint_value,
							       0, 0, A.bigint_value.rows,
							       A.bigint_value.columns);
			else
				SS_base_modul_bigint.insert_at(A.bigint_value, 0, 0, bigint_value,
							       0, 0, A.bigint_value.rows,
							       A.bigint_value.columns);
			S_base_modul_bigint.get_row(bigint_value, tmp,
						    bigint_value.rows - 1);
			break;
		}
		lidia_size_t i;
		for (i = 0; i < bigint_value.columns; i++) {
			v[i].set_modulus(bigint_modulus);
			v[i] = tmp[i];
		}
		delete[] tmp;
	}
	else {
		if (A.long_value.rows > long_value.rows
		    || A.long_value.columns > long_value.columns
		    || l > long_value.columns)
			precondition_error_handler(A.long_value.rows, "A.rows", "A.rows <= rows",
					    A.long_value.columns, "A.columns",
					    "A.columns <= columns",
					    long_value.rows, "rows", "",
					    long_value.columns, "columns", "",
					    l, "v.size", "v.size <= columns",
					    "void base_matrix< bigmod >::"
					    "split_v(base_matrix< bigmod > &A, "
					    "base_vector< bigmod > &v) const",
					    DMESSAGE, EMESSAGE[4]);

		A.set_modulus(long_modulus);

		long *tmp = new long[long_value.columns];
		memory_handler(tmp, DMESSAGE, "split_v :: "
			       "Error in memory allocation (tmp)");

		switch(long_value.bitfield.get_representation()) {
		case matrix_flags::dense_representation:
			if (A.long_value.bitfield.get_representation() ==
			    matrix_flags::dense_representation)
				DD_base_modul_long.insert_at(A.long_value, 0, 0, long_value,
							     0, 0, A.long_value.rows,
							     A.long_value.columns);
			else
				DS_base_modul_long.insert_at(A.long_value, 0, 0, long_value,
							     0, 0, A.long_value.rows,
							     A.long_value.columns);
			D_base_modul_long.get_row(long_value, tmp,
						  long_value.rows - 1);
			break;
		case matrix_flags::sparse_representation:
			if (A.long_value.bitfield.get_representation() ==
			    matrix_flags::dense_representation)
				SD_base_modul_long.insert_at(A.long_value, 0, 0, long_value,
							     0, 0, A.long_value.rows,
							     A.long_value.columns);
			else
				SS_base_modul_long.insert_at(A.long_value, 0, 0, long_value,
							     0, 0, A.long_value.rows,
							     A.long_value.columns);
			S_base_modul_long.get_row(long_value, tmp,
						  long_value.rows - 1);
			break;
		}
		lidia_size_t i;
		for (i = 0; i < long_value.columns; i++) {
			v[i].set_modulus(long_modulus);
			v[i] = tmp[i];
		}
		delete[] tmp;
	}
}



//
// compose functions
//

void base_matrix< bigmod >::
compose_t(const base_matrix< bigmod > &A, const base_matrix< bigmod > &B,
	  const base_matrix< bigmod > &C, const base_matrix< bigmod > &D)
{
	//
	//       Task: E.compose_t(A,B,C,D)
	//             => compose the submatrices A,B,C,D to E
	//                    (A B)
	//             => E = (C D)
	// Conditions: A.columns, B.compose, C.compose, D.compose <= columns and
	//             A.rows, B.rows, C.rows, D.rows <= rows
	//    Version: 2.0
	//

	debug_handler_l(DMESSAGE,
			"compose_t(const base_matrix< bigmod > &, "
			"const base_matrix< bigmod > &, "
			"const base_matrix< bigmod > &, "
			"const base_matrix< bigmod > &)", DVALUE + 3);

	bigint Amod = A.get_modulus();
	bigint Bmod = B.get_modulus();
	bigint Cmod = C.get_modulus();
	bigint Dmod = D.get_modulus();

	if (Amod != Bmod || Amod != Cmod || Amod != Dmod)
		precondition_error_handler(Amod, "A.modulus", "A.modulus == B.modulus",
				    Bmod, "B.modulus", "B.modulus == C.modulus",
				    Cmod, "C.modulus", "C.modulus == D.modulus",
				    Dmod, "D.modulus", "",
				    "compose_t(const base_matrix< bigmod > &A, "
				    "const base_matrix< bigmod > &B, "
				    "const base_matrix< bigmod > &C, "
				    "const base_matrix< bigmod > &D)",
				    DMESSAGE, EMESSAGE[4]);

	set_modulus(Amod);

	if (BigintLong) {
		if (A.bigint_value.rows > bigint_value.rows
		    || A.bigint_value.columns > bigint_value.columns
		    || B.bigint_value.rows > bigint_value.rows
		    || B.bigint_value.columns > bigint_value.columns
		    || C.bigint_value.rows > bigint_value.rows
		    || C.bigint_value.columns > bigint_value.columns
		    || D.bigint_value.rows > bigint_value.rows
		    || D.bigint_value.columns > bigint_value.columns)
			precondition_error_handler(A.bigint_value.rows, "A.rows", "A.rows <= rows",
					    A.bigint_value.columns, "A.columns",
					    "A.columns <= columns",
					    B.bigint_value.rows, "B.rows", "B.rows <= rows",
					    B.bigint_value.columns, "B.columns",
					    "B.columns <= columns",
					    C.bigint_value.rows, "C.rows", "C.rows <= rows",
					    C.bigint_value.columns, "C.columns",
					    "C.columns <= columns",
					    D.bigint_value.rows, "D.rows", "D.rows <= rows",
					    D.bigint_value.columns, "D.columns",
					    "D.columns <= columns",
					    bigint_value.rows, "rows", "",
					    bigint_value.columns, "columns", "",
					    "void base_matrix< bigmod >::"
					    "compose_t(const base_matrix< bigmod > &A, "
					    "const base_matrix< bigmod > &B, "
					    "const base_matrix< bigmod > &C, "
					    "const base_matrix< bigmod > &D)",
					    DMESSAGE, EMESSAGE[4]);

		switch(bigint_value.bitfield.get_representation()) {
		case matrix_flags::dense_representation:
			if (A.bigint_value.bitfield.get_representation() ==
			    matrix_flags::dense_representation)
				DD_base_modul_bigint.insert_at(bigint_value, 0, 0, A.bigint_value,
							       0, 0, A.bigint_value.rows,
							       A.bigint_value.columns);
			else
				DS_base_modul_bigint.insert_at(bigint_value, 0, 0, A.bigint_value,
							       0, 0, A.bigint_value.rows,
							       A.bigint_value.columns);

			if (B.bigint_value.bitfield.get_representation() ==
			    matrix_flags::dense_representation)
				DD_base_modul_bigint.insert_at(bigint_value, 0,
							       bigint_value.columns -
							       B.bigint_value.columns,
							       B.bigint_value, 0, 0,
							       B.bigint_value.rows,
							       B.bigint_value.columns);
			else
				DS_base_modul_bigint.insert_at(bigint_value, 0,
							       bigint_value.columns -
							       B.bigint_value.columns,
							       B.bigint_value, 0, 0,
							       B.bigint_value.rows,
							       B.bigint_value.columns);

			if (C.bigint_value.bitfield.get_representation() ==
			    matrix_flags::dense_representation)
				DD_base_modul_bigint.insert_at(bigint_value, bigint_value.rows -
							       C.bigint_value.rows, 0,
							       C.bigint_value, 0, 0,
							       C.bigint_value.rows,
							       C.bigint_value.columns);
			else
				DS_base_modul_bigint.insert_at(bigint_value, bigint_value.rows -
							       C.bigint_value.rows, 0,
							       C.bigint_value, 0, 0,
							       C.bigint_value.rows,
							       C.bigint_value.columns);

			if (D.bigint_value.bitfield.get_representation() ==
			    matrix_flags::dense_representation)
				DD_base_modul_bigint.insert_at(bigint_value, bigint_value.rows -
							       D.bigint_value.rows,
							       bigint_value.columns -
							       D.bigint_value.columns,
							       D.bigint_value, 0, 0,
							       D.bigint_value.rows,
							       D.bigint_value.columns);
			else
				DS_base_modul_bigint.insert_at(bigint_value, bigint_value.rows -
							       D.bigint_value.rows,
							       bigint_value.columns -
							       D.bigint_value.columns,
							       D.bigint_value, 0, 0,
							       D.bigint_value.rows,
							       D.bigint_value.columns);
			break;
		case matrix_flags::sparse_representation:
			if (A.bigint_value.bitfield.get_representation() ==
			    matrix_flags::dense_representation)
				SD_base_modul_bigint.insert_at(bigint_value, 0, 0, A.bigint_value,
							       0, 0, A.bigint_value.rows,
							       A.bigint_value.columns);
			else
				SS_base_modul_bigint.insert_at(bigint_value, 0, 0, A.bigint_value,
							       0, 0, A.bigint_value.rows,
							       A.bigint_value.columns);

			if (B.bigint_value.bitfield.get_representation() ==
			    matrix_flags::dense_representation)
				SD_base_modul_bigint.insert_at(bigint_value, 0,
							       bigint_value.columns -
							       B.bigint_value.columns,
							       B.bigint_value, 0, 0,
							       B.bigint_value.rows,
							       B.bigint_value.columns);
			else
				SS_base_modul_bigint.insert_at(bigint_value, 0,
							       bigint_value.columns -
							       B.bigint_value.columns,
							       B.bigint_value, 0, 0,
							       B.bigint_value.rows,
							       B.bigint_value.columns);

			if (C.bigint_value.bitfield.get_representation() ==
			    matrix_flags::dense_representation)
				SD_base_modul_bigint.insert_at(bigint_value, bigint_value.rows -
							       C.bigint_value.rows, 0,
							       C.bigint_value, 0, 0,
							       C.bigint_value.rows,
							       C.bigint_value.columns);
			else
				SS_base_modul_bigint.insert_at(bigint_value, bigint_value.rows -
							       C.bigint_value.rows, 0,
							       C.bigint_value, 0, 0,
							       C.bigint_value.rows,
							       C.bigint_value.columns);

			if (D.bigint_value.bitfield.get_representation() ==
			    matrix_flags::dense_representation)
				SD_base_modul_bigint.insert_at(bigint_value, bigint_value.rows -
							       D.bigint_value.rows,
							       bigint_value.columns -
							       D.bigint_value.columns,
							       D.bigint_value, 0, 0,
							       D.bigint_value.rows,
							       D.bigint_value.columns);
			else
				SS_base_modul_bigint.insert_at(bigint_value, bigint_value.rows -
							       D.bigint_value.rows,
							       bigint_value.columns -
							       D.bigint_value.columns,
							       D.bigint_value, 0, 0,
							       D.bigint_value.rows,
							       D.bigint_value.columns);
			break;
		}
	}
	else {
		if (A.long_value.rows > long_value.rows
		    || A.long_value.columns > long_value.columns
		    || B.long_value.rows > long_value.rows
		    || B.long_value.columns > long_value.columns
		    || C.long_value.rows > long_value.rows
		    || C.long_value.columns > long_value.columns
		    || D.long_value.rows > long_value.rows
		    || D.long_value.columns > long_value.columns)
			precondition_error_handler(A.long_value.rows, "A.rows", "A.rows <= rows",
					    A.long_value.columns, "A.columns",
					    "A.columns <= columns",
					    B.long_value.rows, "B.rows", "B.rows <= rows",
					    B.long_value.columns, "B.columns",
					    "B.columns <= columns",
					    C.long_value.rows, "C.rows", "C.rows <= rows",
					    C.long_value.columns, "C.columns",
					    "C.columns <= columns",
					    D.long_value.rows, "D.rows", "D.rows <= rows",
					    D.long_value.columns, "D.columns",
					    "D.columns <= columns",
					    long_value.rows, "rows", "",
					    long_value.columns, "columns", "",
					    "void base_matrix< bigmod >::"
					    "compose_t(const base_matrix< bigmod > &A, "
					    "const base_matrix< bigmod > &B, "
					    "const base_matrix< bigmod > &C, "
					    "const base_matrix< bigmod > &D)",
					    DMESSAGE, EMESSAGE[4]);

		switch(long_value.bitfield.get_representation()) {
		case matrix_flags::dense_representation:
			if (A.long_value.bitfield.get_representation() ==
			    matrix_flags::dense_representation)
				DD_base_modul_long.insert_at(long_value, 0, 0, A.long_value,
							     0, 0, A.long_value.rows,
							     A.long_value.columns);
			else
				DS_base_modul_long.insert_at(long_value, 0, 0, A.long_value,
							     0, 0, A.long_value.rows,
							     A.long_value.columns);

			if (B.long_value.bitfield.get_representation() ==
			    matrix_flags::dense_representation)
				DD_base_modul_long.insert_at(long_value, 0,
							     long_value.columns -
							     B.long_value.columns,
							     B.long_value, 0, 0,
							     B.long_value.rows,
							     B.long_value.columns);
			else
				DS_base_modul_long.insert_at(long_value, 0,
							     long_value.columns -
							     B.long_value.columns,
							     B.long_value, 0, 0,
							     B.long_value.rows,
							     B.long_value.columns);

			if (C.long_value.bitfield.get_representation() ==
			    matrix_flags::dense_representation)
				DD_base_modul_long.insert_at(long_value, long_value.rows -
							     C.long_value.rows, 0,
							     C.long_value, 0, 0,
							     C.long_value.rows,
							     C.long_value.columns);
			else
				DS_base_modul_long.insert_at(long_value, long_value.rows -
							     C.long_value.rows, 0,
							     C.long_value, 0, 0,
							     C.long_value.rows,
							     C.long_value.columns);

			if (D.long_value.bitfield.get_representation() ==
			    matrix_flags::dense_representation)
				DD_base_modul_long.insert_at(long_value, long_value.rows -
							     D.long_value.rows,
							     long_value.columns -
							     D.long_value.columns,
							     D.long_value, 0, 0,
							     D.long_value.rows,
							     D.long_value.columns);
			else
				DS_base_modul_long.insert_at(long_value, long_value.rows -
							     D.long_value.rows,
							     long_value.columns -
							     D.long_value.columns,
							     D.long_value, 0, 0,
							     D.long_value.rows,
							     D.long_value.columns);
			break;
		case matrix_flags::sparse_representation:
			if (A.long_value.bitfield.get_representation() ==
			    matrix_flags::dense_representation)
				SD_base_modul_long.insert_at(long_value, 0, 0, A.long_value,
							     0, 0, A.long_value.rows,
							     A.long_value.columns);
			else
				SS_base_modul_long.insert_at(long_value, 0, 0, A.long_value,
							     0, 0, A.long_value.rows,
							     A.long_value.columns);

			if (B.long_value.bitfield.get_representation() ==
			    matrix_flags::dense_representation)
				SD_base_modul_long.insert_at(long_value, 0,
							     long_value.columns -
							     B.long_value.columns,
							     B.long_value, 0, 0,
							     B.long_value.rows,
							     B.long_value.columns);
			else
				SS_base_modul_long.insert_at(long_value, 0,
							     long_value.columns -
							     B.long_value.columns,
							     B.long_value, 0, 0,
							     B.long_value.rows,
							     B.long_value.columns);

			if (C.long_value.bitfield.get_representation() ==
			    matrix_flags::dense_representation)
				SD_base_modul_long.insert_at(long_value, long_value.rows -
							     C.long_value.rows, 0,
							     C.long_value, 0, 0,
							     C.long_value.rows,
							     C.long_value.columns);
			else
				SS_base_modul_long.insert_at(long_value, long_value.rows -
							     C.long_value.rows, 0,
							     C.long_value, 0, 0,
							     C.long_value.rows,
							     C.long_value.columns);

			if (D.long_value.bitfield.get_representation() ==
			    matrix_flags::dense_representation)
				SD_base_modul_long.insert_at(long_value, long_value.rows -
							     D.long_value.rows,
							     long_value.columns -
							     D.long_value.columns,
							     D.long_value, 0, 0,
							     D.long_value.rows,
							     D.long_value.columns);
			else
				SS_base_modul_long.insert_at(long_value, long_value.rows -
							     D.long_value.rows,
							     long_value.columns -
							     D.long_value.columns,
							     D.long_value, 0, 0,
							     D.long_value.rows,
							     D.long_value.columns);
			break;
		}
	}
}



void base_matrix< bigmod >::
compose_h(const base_matrix< bigmod > &A, const base_matrix< bigmod > &B)
{
	//
	//       Task: C.compose_h(A,B);
	//             => compose the submatrices A,B to C
	//             => C = (A B)
	// Conditions: A.columns, B.columns <= columns and
	//             A.rows, B.rows <= rows
	//    Version: 2.0
	//

	debug_handler_l(DMESSAGE, "compose_h(const base_matrix< bigmod > &, "
			"const base_matrix< bigmod > &)", DVALUE + 3);

	bigint Amod = A.get_modulus();
	bigint Bmod = B.get_modulus();

	if (Amod != Bmod)
		precondition_error_handler(Amod, "A.modulus", "A.modulus == B.modulus",
				    Bmod, "B.modulus", "",
				    "compose_h(const base_matrix< bigmod > &A, "
				    "const base_matrix< bigmod > &B)",
				    DMESSAGE, EMESSAGE[4]);

	set_modulus(Amod);

	if (BigintLong) {
		if (A.bigint_value.rows > bigint_value.rows
		    || A.bigint_value.columns > bigint_value.columns
		    || B.bigint_value.rows > bigint_value.rows
		    || B.bigint_value.columns > bigint_value.columns)
			precondition_error_handler(A.bigint_value.rows, "A.rows", "A.rows <= rows",
					    A.bigint_value.columns, "A.columns",
					    "A.columns <= columns",
					    B.bigint_value.rows, "B.rows", "B.rows <= rows",
					    B.bigint_value.columns, "B.columns",
					    "B.columns <= columns",
					    bigint_value.rows, "rows", "",
					    bigint_value.columns, "columns", "",
					    "void base_matrix< bigmod >::"
					    "compose_h(const base_matrix< bigmod > &A, "
					    "const base_matrix< bigmod > &B)",
					    DMESSAGE, EMESSAGE[4]);

		switch(bigint_value.bitfield.get_representation()) {
		case matrix_flags::dense_representation:
			if (A.bigint_value.bitfield.get_representation() ==
			    matrix_flags::dense_representation)
				DD_base_modul_bigint.insert_at(bigint_value, 0, 0, A.bigint_value,
							       0, 0, A.bigint_value.rows,
							       A.bigint_value.columns);
			else
				DS_base_modul_bigint.insert_at(bigint_value, 0, 0, A.bigint_value,
							       0, 0, A.bigint_value.rows,
							       A.bigint_value.columns);

			if (B.bigint_value.bitfield.get_representation() ==
			    matrix_flags::dense_representation)
				DD_base_modul_bigint.insert_at(bigint_value, 0,
							       bigint_value.columns -
							       B.bigint_value.columns,
							       B.bigint_value, 0, 0,
							       B.bigint_value.rows,
							       B.bigint_value.columns);
			else
				DD_base_modul_bigint.insert_at(bigint_value, 0,
							       bigint_value.columns -
							       B.bigint_value.columns,
							       B.bigint_value, 0, 0,
							       B.bigint_value.rows,
							       B.bigint_value.columns);
			break;
		case matrix_flags::sparse_representation:
			if (A.bigint_value.bitfield.get_representation() ==
			    matrix_flags::dense_representation)
				SD_base_modul_bigint.insert_at(bigint_value, 0, 0,
							       A.bigint_value, 0, 0,
							       A.bigint_value.rows,
							       A.bigint_value.columns);
			else
				SS_base_modul_bigint.insert_at(bigint_value, 0, 0,
							       A.bigint_value, 0, 0,
							       A.bigint_value.rows,
							       A.bigint_value.columns);

			if (B.bigint_value.bitfield.get_representation() ==
			    matrix_flags::dense_representation)
				SD_base_modul_bigint.insert_at(bigint_value, 0,
							       bigint_value.columns -
							       B.bigint_value.columns,
							       B.bigint_value, 0, 0,
							       B.bigint_value.rows,
							       B.bigint_value.columns);
			else
				SD_base_modul_bigint.insert_at(bigint_value, 0,
							       bigint_value.columns -
							       B.bigint_value.columns,
							       B.bigint_value, 0, 0,
							       B.bigint_value.rows,
							       B.bigint_value.columns);
			break;
		}
	}
	else {
		if (A.long_value.rows > long_value.rows
		    || A.long_value.columns > long_value.columns
		    || B.long_value.rows > long_value.rows
		    || B.long_value.columns > long_value.columns)
			precondition_error_handler(A.long_value.rows, "A.rows", "A.rows <= rows",
					    A.long_value.columns, "A.columns",
					    "A.columns <= columns",
					    B.long_value.rows, "B.rows", "B.rows <= rows",
					    B.long_value.columns, "B.columns",
					    "B.columns <= columns",
					    long_value.rows, "rows", "",
					    long_value.columns, "columns", "",
					    "void base_matrix< bigmod >::"
					    "compose_h(const base_matrix< bigmod > &A, "
					    "const base_matrix< bigmod > &B)",
					    DMESSAGE, EMESSAGE[4]);

		switch(long_value.bitfield.get_representation()) {
		case matrix_flags::dense_representation:
			if (A.long_value.bitfield.get_representation() ==
			    matrix_flags::dense_representation)
				DD_base_modul_long.insert_at(long_value, 0, 0, A.long_value,
							     0, 0, A.long_value.rows,
							     A.long_value.columns);
			else
				DS_base_modul_long.insert_at(long_value, 0, 0, A.long_value,
							     0, 0, A.long_value.rows,
							     A.long_value.columns);

			if (B.long_value.bitfield.get_representation() ==
			    matrix_flags::dense_representation)
				DD_base_modul_long.insert_at(long_value, 0,
							     long_value.columns -
							     B.long_value.columns,
							     B.long_value, 0, 0,
							     B.long_value.rows,
							     B.long_value.columns);
			else
				DD_base_modul_long.insert_at(long_value, 0,
							     long_value.columns -
							     B.long_value.columns,
							     B.long_value, 0, 0,
							     B.long_value.rows,
							     B.long_value.columns);
			break;
		case matrix_flags::sparse_representation:
			if (A.long_value.bitfield.get_representation() ==
			    matrix_flags::dense_representation)
				SD_base_modul_long.insert_at(long_value, 0, 0,
							     A.long_value, 0, 0,
							     A.long_value.rows,
							     A.long_value.columns);
			else
				SS_base_modul_long.insert_at(long_value, 0, 0,
							     A.long_value, 0, 0,
							     A.long_value.rows,
							     A.long_value.columns);

			if (B.long_value.bitfield.get_representation() ==
			    matrix_flags::dense_representation)
				SD_base_modul_long.insert_at(long_value, 0,
							     long_value.columns -
							     B.long_value.columns,
							     B.long_value, 0, 0,
							     B.long_value.rows,
							     B.long_value.columns);
			else
				SD_base_modul_long.insert_at(long_value, 0,
							     long_value.columns -
							     B.long_value.columns,
							     B.long_value, 0, 0,
							     B.long_value.rows,
							     B.long_value.columns);
			break;
		}
	}
}



void base_matrix< bigmod >::
compose_h(const bigmod *v, const base_matrix< bigmod > &A)
{
	//
	//       Task: C.compose_h(v, A);
	//             => compose the submatrix A and the vector v to matrix C
	//             => C = (v A)
	// Conditions: A.columns <= columns and
	//             A.rows <= rows and
	//             v != NULL
	//    Version: 2.0
	//

	debug_handler_l(DMESSAGE, "compose_h(const bigmod *, "
			"const base_matrix< bigmod > &)", DVALUE + 3);

	if (BigintLong) {
		if (A.bigint_value.rows > bigint_value.rows
		    || A.bigint_value.columns > bigint_value.columns || v == NULL)
			precondition_error_handler(A.bigint_value.rows, "A.rows", "A.rows <= rows",
					    A.bigint_value.columns, "A.columns",
					    "A.columns <= columns",
					    bigint_value.rows, "rows", "",
					    bigint_value.columns, "columns", "",
					    PRT, "v", "v != NULL",
					    "void base_matrix< bigmod >::"
					    "compose_h(const bigmod *v, "
					    "const base_matrix< bigmod > &A)",
					    DMESSAGE, EMESSAGE[4]);

		bigint *tmp = new bigint[bigint_value.rows];
		memory_handler(tmp, DMESSAGE, "compose_h :: "
			       "Error in memory allocation (tmp)");
		lidia_size_t i;
		for (i = 0; i < bigint_value.rows; i++)
			tmp[i] = v[i].mantissa();

		switch(bigint_value.bitfield.get_representation()) {
		case matrix_flags::dense_representation:
			D_base_modul_bigint.insert_column_at(bigint_value, 0, tmp,
							     bigint_value.rows);
			if (A.bigint_value.bitfield.get_representation() ==
			    matrix_flags::dense_representation)
				DD_base_modul_bigint.insert_at(bigint_value, 0,
							       bigint_value.columns -
							       A.bigint_value.columns,
							       A.bigint_value, 0, 0,
							       A.bigint_value.rows,
							       A.bigint_value.columns);
			else
				DS_base_modul_bigint.insert_at(bigint_value, 0,
							       bigint_value.columns -
							       A.bigint_value.columns,
							       A.bigint_value, 0, 0,
							       A.bigint_value.rows,
							       A.bigint_value.columns);
			break;
		case matrix_flags::sparse_representation:
			S_base_modul_bigint.insert_column_at(bigint_value, 0, tmp,
							     bigint_value.rows);
			if (A.bigint_value.bitfield.get_representation() ==
			    matrix_flags::dense_representation)
				SD_base_modul_bigint.insert_at(bigint_value, 0,
							       bigint_value.columns -
							       A.bigint_value.columns,
							       A.bigint_value, 0, 0,
							       A.bigint_value.rows,
							       A.bigint_value.columns);
			else
				SS_base_modul_bigint.insert_at(bigint_value, 0,
							       bigint_value.columns -
							       A.bigint_value.columns,
							       A.bigint_value, 0, 0,
							       A.bigint_value.rows,
							       A.bigint_value.columns);
			break;
		}
		delete[] tmp;
	}
	else {
		if (A.long_value.rows > long_value.rows
		    || A.long_value.columns > long_value.columns || v == NULL)
			precondition_error_handler(A.long_value.rows, "A.rows", "A.rows <= rows",
					    A.long_value.columns, "A.columns",
					    "A.columns <= columns",
					    long_value.rows, "rows", "",
					    long_value.columns, "columns", "",
					    PRT, "v", "v != NULL",
					    "void base_matrix< bigmod >::"
					    "compose_h(const bigmod *v, "
					    "const base_matrix< bigmod > &A)",
					    DMESSAGE, EMESSAGE[4]);

		long *tmp = new long[long_value.rows];
		memory_handler(tmp, DMESSAGE, "compose_h :: "
			       "Error in memory allocation (tmp)");
		lidia_size_t i;
		for (i = 0; i < long_value.rows; i++)
			(v[i].mantissa()).longify(tmp[i]);

		switch(long_value.bitfield.get_representation()) {
		case matrix_flags::dense_representation:
			D_base_modul_long.insert_column_at(long_value, 0, tmp,
							   long_value.rows);




			if (A.long_value.bitfield.get_representation() ==
			    matrix_flags::dense_representation)
				DD_base_modul_long.insert_at(long_value, 0,
							     long_value.columns -
							     A.long_value.columns,
							     A.long_value, 0, 0,
							     A.long_value.rows,
							     A.long_value.columns);
			else
				DS_base_modul_long.insert_at(long_value, 0,
							     long_value.columns -
							     A.long_value.columns,
							     A.long_value, 0, 0,
							     A.long_value.rows,
							     A.long_value.columns);
			break;
		case matrix_flags::sparse_representation:
			S_base_modul_long.insert_column_at(long_value, 0, tmp,
							   long_value.rows);
			if (A.long_value.bitfield.get_representation() ==
			    matrix_flags::dense_representation)
				SD_base_modul_long.insert_at(long_value, 0,
							     long_value.columns -
							     A.long_value.columns,
							     A.long_value, 0, 0,
							     A.long_value.rows,
							     A.long_value.columns);
			else
				SS_base_modul_long.insert_at(long_value, 0,
							     long_value.columns -
							     A.long_value.columns,
							     A.long_value, 0, 0,
							     A.long_value.rows,
							     A.long_value.columns);
			break;
		}
		delete[] tmp;
	}
}



void base_matrix< bigmod >::
compose_h(const base_matrix< bigmod > &A, const bigmod *v)
{
	//
	//       Task: C.compose_h(A, v);
	//             => compose the submatrix A and vector v to matrix C
	//             => C = (A v)
	// Conditions: A.columns <= columns and
	//             A.rows <= rows and
	//             v != NULL
	//    Version: 2.0
	//

	debug_handler_l(DMESSAGE, "compose_h(const base_matrix< bigmod > &, "
			"const bigmod *)", DVALUE + 3);

	if (BigintLong) {
		if (A.bigint_value.rows > bigint_value.rows
		    || A.bigint_value.columns > bigint_value.columns || v == NULL)
			precondition_error_handler(A.bigint_value.rows, "A.rows", "A.rows <= rows",
					    A.bigint_value.columns, "A.columns",
					    "A.columns <= columns",
					    bigint_value.rows, "rows", "",
					    bigint_value.columns, "columns", "",
					    PRT, "v", "v != NULL",
					    "void base_matrix< bigmod >::"
					    "compose_h(const base_matrix< bigmod > &A, "
					    "const bigmod *v)", DMESSAGE, EMESSAGE[4]);

		bigint *tmp = new bigint[bigint_value.rows];
		memory_handler(tmp, DMESSAGE, "compose_h :: "
			       "Error in memory allocation (tmp)");
		lidia_size_t i;
		for (i = 0; i < bigint_value.rows; i++)
			tmp[i] = v[i].mantissa();

		switch(bigint_value.bitfield.get_representation()) {
		case matrix_flags::dense_representation:
			if (A.bigint_value.bitfield.get_representation() ==
			    matrix_flags::dense_representation)
				DD_base_modul_bigint.insert_at(bigint_value, 0, 0,
							       A.bigint_value, 0, 0,
							       A.bigint_value.rows,
							       A.bigint_value.columns);
			else
				DS_base_modul_bigint.insert_at(bigint_value, 0, 0,
							       A.bigint_value, 0, 0,
							       A.bigint_value.rows,
							       A.bigint_value.columns);
			D_base_modul_bigint.insert_column_at(bigint_value,
							     bigint_value.columns - 1, tmp,
							     bigint_value.rows);
			break;
		case matrix_flags::sparse_representation:
			if (A.bigint_value.bitfield.get_representation() ==
			    matrix_flags::dense_representation)
				SD_base_modul_bigint.insert_at(bigint_value, 0, 0,
							       A.bigint_value, 0, 0,
							       A.bigint_value.rows,
							       A.bigint_value.columns);
			else
				SS_base_modul_bigint.insert_at(bigint_value, 0, 0,
							       A.bigint_value, 0, 0,
							       A.bigint_value.rows,
							       A.bigint_value.columns);
			D_base_modul_bigint.insert_column_at(bigint_value,
							     bigint_value.columns - 1, tmp,
							     bigint_value.rows);
			break;
		}
		delete[] tmp;
	}
	else {
		if (A.long_value.rows > long_value.rows
		    || A.long_value.columns > long_value.columns || v == NULL)
			precondition_error_handler(A.long_value.rows, "A.rows", "A.rows <= rows",
					    A.long_value.columns, "A.columns",
					    "A.columns <= columns",
					    long_value.rows, "rows", "",
					    long_value.columns, "columns", "",
					    PRT, "v", "v != NULL",
					    "void base_matrix< bigmod >::"
					    "compose_h(const base_matrix< bigmod > &A, "
					    "const bigmod *v)", DMESSAGE, EMESSAGE[4]);

		long *tmp = new long[long_value.rows];
		memory_handler(tmp, DMESSAGE, "compose_h :: "
			       "Error in memory allocation (tmp)");
		lidia_size_t i;
		for (i = 0; i < long_value.rows; i++)
			(v[i].mantissa()).longify(tmp[i]);

		switch(long_value.bitfield.get_representation()) {
		case matrix_flags::dense_representation:
			if (A.long_value.bitfield.get_representation() ==
			    matrix_flags::dense_representation)
				DD_base_modul_long.insert_at(long_value, 0, 0,
							     A.long_value, 0, 0,
							     A.long_value.rows,
							     A.long_value.columns);
			else
				DS_base_modul_long.insert_at(long_value, 0, 0,
							     A.long_value, 0, 0,
							     A.long_value.rows,
							     A.long_value.columns);
			D_base_modul_long.insert_column_at(long_value,
							   long_value.columns - 1, tmp,
							   long_value.rows);
			break;
		case matrix_flags::sparse_representation:
			if (A.long_value.bitfield.get_representation() ==
			    matrix_flags::dense_representation)
				SD_base_modul_long.insert_at(long_value, 0, 0,
							     A.long_value, 0, 0,
							     A.long_value.rows,
							     A.long_value.columns);
			else
				SS_base_modul_long.insert_at(long_value, 0, 0,
							     A.long_value, 0, 0,
							     A.long_value.rows,
							     A.long_value.columns);
			D_base_modul_long.insert_column_at(long_value,
							   long_value.columns - 1, tmp,
							   long_value.rows);
			break;
		}
		delete[] tmp;
	}
}



void base_matrix< bigmod >::
compose_h(const base_vector< bigmod > &v, const base_matrix< bigmod > &A)
{
	//
	//       Task: C.compose_h(v, A);
	//             => compose the submatrix A and the vector v to matrix C
	//             => C = (v A)
	// Conditions: A.columns <= columns and
	//             A.rows <= rows and
	//             v.size > rows
	//    Version: 2.0
	//

	debug_handler_l(DMESSAGE, "compose_h(const base_vector< bigmod > &, "
			"const base_matrix< bigmod > &)", DVALUE + 3);

	lidia_size_t l = v.size();

	if (BigintLong) {
		if (A.bigint_value.rows > bigint_value.rows
		    || A.bigint_value.columns > bigint_value.columns
		    || l > bigint_value.rows)
			precondition_error_handler(A.bigint_value.rows, "A.rows", "A.rows <= rows",
					    A.bigint_value.columns, "A.columns",
					    "A.columns <= columns",
					    bigint_value.rows, "rows", "",
					    bigint_value.columns, "columns", "",
					    l, "v.size", "v.size <= rows",
					    "void base_matrix< bigmod >::"
					    "compose_h(const base_vector< bigmod > &v, "
					    "const base_matrix< bigmod > &A)",
					    DMESSAGE, EMESSAGE[4]);

		bigint *tmp = new bigint[bigint_value.rows];
		memory_handler(tmp, DMESSAGE, "compose_h :: "
			       "Error in memory allocation (tmp)");
		lidia_size_t i;
		for (i = 0; i < bigint_value.rows; i++)
			tmp[i] = v[i].mantissa();

		switch(bigint_value.bitfield.get_representation()) {
		case matrix_flags::dense_representation:
			D_base_modul_bigint.insert_column_at(bigint_value, 0, tmp, l);
			if (A.bigint_value.bitfield.get_representation() ==
			    matrix_flags::dense_representation)
				DD_base_modul_bigint.insert_at(bigint_value, 0,
							       bigint_value.columns -
							       A.bigint_value.columns,
							       A.bigint_value, 0, 0,
							       A.bigint_value.rows,
							       A.bigint_value.columns);
			else
				DS_base_modul_bigint.insert_at(bigint_value, 0,
							       bigint_value.columns -
							       A.bigint_value.columns,
							       A.bigint_value, 0, 0,
							       A.bigint_value.rows,
							       A.bigint_value.columns);
			break;
		case matrix_flags::sparse_representation:
			S_base_modul_bigint.insert_column_at(bigint_value, 0, tmp, l);
			if (A.bigint_value.bitfield.get_representation() ==
			    matrix_flags::dense_representation)
				SD_base_modul_bigint.insert_at(bigint_value, 0,
							       bigint_value.columns -
							       A.bigint_value.columns,
							       A.bigint_value, 0, 0,
							       A.bigint_value.rows,
							       A.bigint_value.columns);
			else
				SS_base_modul_bigint.insert_at(bigint_value, 0,
							       bigint_value.columns -
							       A.bigint_value.columns,
							       A.bigint_value, 0, 0,
							       A.bigint_value.rows,
							       A.bigint_value.columns);
			break;
		}
		delete[] tmp;
	}
	else {
		if (A.long_value.rows > long_value.rows
		    || A.long_value.columns > long_value.columns
		    || l > long_value.rows)
			precondition_error_handler(A.long_value.rows, "A.rows", "A.rows <= rows",
					    A.long_value.columns, "A.columns",
					    "A.columns <= columns",
					    long_value.rows, "rows", "",
					    long_value.columns, "columns", "",
					    l, "v.size", "v.size <= rows",
					    "void base_matrix< bigmod >::"
					    "compose_h(const base_vector< bigmod > &v, "
					    "const base_matrix< bigmod > &A)",
					    DMESSAGE, EMESSAGE[4]);

		long *tmp = new long[long_value.rows];
		memory_handler(tmp, DMESSAGE, "compose_h :: "
			       "Error in memory allocation (tmp)");
		lidia_size_t i;
		for (i = 0; i < long_value.rows; i++)
			(v[i].mantissa()).longify(tmp[i]);

		switch(long_value.bitfield.get_representation()) {
		case matrix_flags::dense_representation:
			D_base_modul_long.insert_column_at(long_value, 0, tmp, l);
			if (A.long_value.bitfield.get_representation() ==
			    matrix_flags::dense_representation)
				DD_base_modul_long.insert_at(long_value, 0,
							     long_value.columns -
							     A.long_value.columns,
							     A.long_value, 0, 0,
							     A.long_value.rows,
							     A.long_value.columns);
			else
				DS_base_modul_long.insert_at(long_value, 0,
							     long_value.columns -
							     A.long_value.columns,
							     A.long_value, 0, 0,
							     A.long_value.rows,
							     A.long_value.columns);
			break;
		case matrix_flags::sparse_representation:
			S_base_modul_long.insert_column_at(long_value, 0, tmp, l);
			if (A.long_value.bitfield.get_representation() ==
			    matrix_flags::dense_representation)
				SD_base_modul_long.insert_at(long_value, 0,
							     long_value.columns -
							     A.long_value.columns,
							     A.long_value, 0, 0,
							     A.long_value.rows,
							     A.long_value.columns);
			else
				SS_base_modul_long.insert_at(long_value, 0,
							     long_value.columns -
							     A.long_value.columns,
							     A.long_value, 0, 0,
							     A.long_value.rows,
							     A.long_value.columns);
			break;
		}
		delete[] tmp;
	}
}



void base_matrix< bigmod >::
compose_h(const base_matrix< bigmod > &A, const base_vector< bigmod > &v)
{
	//
	//       Task: C.compose_h(A, v);
	//             => compose the submatrix A and vector v to matrix C
	//             => C = (A v)
	// Conditions: A.columns <= columns and
	//             A.rows <= rows and
	//             v.size <= rows
	//    Version: 2.0
	//

	debug_handler_l(DMESSAGE, "compose_h(const base_matrix< bigmod > &, "
			"base_vector< bigmod > &)", DVALUE + 3);

	lidia_size_t l = v.size();

	if (BigintLong) {
		if (A.bigint_value.rows > bigint_value.rows
		    || A.bigint_value.columns > bigint_value.columns
		    || l > bigint_value.rows)
			precondition_error_handler(A.bigint_value.rows, "A.rows", "A.rows <= rows",
					    A.bigint_value.columns, "A.columns",
					    "A.columns <= columns",
					    bigint_value.rows, "rows", "",
					    bigint_value.columns, "columns", "",
					    l, "v.size", "v.size <= rows",
					    "void base_matrix< bigmod >::"
					    "compose_h(const base_matrix< bigmod > &A, "
					    "const base_vector< bigmod > &v)",
					    DMESSAGE, EMESSAGE[4]);

		bigint *tmp = new bigint[bigint_value.rows];
		memory_handler(tmp, DMESSAGE, "compose_h :: "
			       "Error in memory allocation (tmp)");
		lidia_size_t i;
		for (i = 0; i < bigint_value.rows; i++)
			tmp[i] = v[i].mantissa();

		switch(bigint_value.bitfield.get_representation()) {
		case matrix_flags::dense_representation:
			if (A.bigint_value.bitfield.get_representation() ==
			    matrix_flags::dense_representation)
				DD_base_modul_bigint.insert_at(bigint_value, 0, 0, A.bigint_value,
							       0, 0, A.bigint_value.rows,
							       A.bigint_value.columns);
			else
				DS_base_modul_bigint.insert_at(bigint_value, 0, 0, A.bigint_value,
							       0, 0, A.bigint_value.rows,
							       A.bigint_value.columns);
			D_base_modul_bigint.insert_column_at(bigint_value,
							     bigint_value.columns - 1, tmp,
							     l);
			break;
		case matrix_flags::sparse_representation:
			if (A.bigint_value.bitfield.get_representation() ==
			    matrix_flags::dense_representation)
				SD_base_modul_bigint.insert_at(bigint_value, 0, 0, A.bigint_value,
							       0, 0, A.bigint_value.rows,
							       A.bigint_value.columns);
			else
				SS_base_modul_bigint.insert_at(bigint_value, 0, 0, A.bigint_value,
							       0, 0, A.bigint_value.rows,
							       A.bigint_value.columns);
			D_base_modul_bigint.insert_column_at(bigint_value,
							     bigint_value.columns - 1, tmp,
							     l);
			break;
		}
		delete[] tmp;
	}
	else {
		if (A.long_value.rows > long_value.rows
		    || A.long_value.columns > long_value.columns
		    || l > long_value.rows)
			precondition_error_handler(A.long_value.rows, "A.rows", "A.rows <= rows",
					    A.long_value.columns, "A.columns",
					    "A.columns <= columns",
					    long_value.rows, "rows", "",
					    long_value.columns, "columns", "",
					    l, "v.size", "v.size <= rows",
					    "void base_matrix< bigmod >::"
					    "compose_h(const base_matrix< bigmod > &A, "
					    "const base_vector< bigmod > &v)",
					    DMESSAGE, EMESSAGE[4]);

		long *tmp = new long[long_value.rows];
		memory_handler(tmp, DMESSAGE, "compose_h :: "
			       "Error in memory allocation (tmp)");
		lidia_size_t i;
		for (i = 0; i < long_value.rows; i++)
			(v[i].mantissa()).longify(tmp[i]);

		switch(long_value.bitfield.get_representation()) {
		case matrix_flags::dense_representation:
			if (A.long_value.bitfield.get_representation() ==
			    matrix_flags::dense_representation)
				DD_base_modul_long.insert_at(long_value, 0, 0, A.long_value,
							     0, 0, A.long_value.rows,
							     A.long_value.columns);
			else
				DS_base_modul_long.insert_at(long_value, 0, 0, A.long_value,
							     0, 0, A.long_value.rows,
							     A.long_value.columns);
			D_base_modul_long.insert_column_at(long_value,
							   long_value.columns - 1, tmp,
							   l);
			break;
		case matrix_flags::sparse_representation:
			if (A.long_value.bitfield.get_representation() ==
			    matrix_flags::dense_representation)
				SD_base_modul_long.insert_at(long_value, 0, 0, A.long_value,
							     0, 0, A.long_value.rows,
							     A.long_value.columns);
			else
				SS_base_modul_long.insert_at(long_value, 0, 0, A.long_value,
							     0, 0, A.long_value.rows,
							     A.long_value.columns);
			D_base_modul_long.insert_column_at(long_value,
							   long_value.columns - 1, tmp,
							   l);
			break;
		}
		delete[] tmp;
	}
}



void base_matrix< bigmod >::
compose_v(const base_matrix< bigmod > &A, const base_matrix< bigmod > &B)
{
	//
	//       Task: C.compose_v(A,B);
	//             => compose the submatrices A,B to C
	//                    (A)
	//             => C = (B)
	// Conditions: A.columns, B.columns <= columns and
	//             A.rows, B.rows <= rows
	//    Version: 2.0
	//

	debug_handler_l(DMESSAGE, "compose_v(const base_matrix< bigmod > &, "
			"const base_matrix< bigmod > &)", DVALUE + 3);

	bigint Amod = A.get_modulus();
	bigint Bmod = B.get_modulus();

	if (Amod != Bmod)
		precondition_error_handler(Amod, "A.modulus", "A.modulus == B.modulus",
				    Bmod, "B.modulus", "",
				    "compose_v(const base_matrix< bigmod > &A, "
				    "const base_matrix< bigmod > &B)",
				    DMESSAGE, EMESSAGE[4]);

	set_modulus(Amod);

	if (BigintLong) {
		if (A.bigint_value.rows > bigint_value.rows
		    || A.bigint_value.columns > bigint_value.columns
		    || B.bigint_value.rows > bigint_value.rows
		    || B.bigint_value.columns > bigint_value.columns)
			precondition_error_handler(A.bigint_value.rows, "A.rows", "A.rows <= rows",
					    A.bigint_value.columns, "A.columns",
					    "A.columns <= columns",
					    B.bigint_value.rows, "B.rows", "B.rows <= rows",
					    B.bigint_value.columns, "B.columns",
					    "B.columns <= columns",
					    bigint_value.rows, "rows", "",
					    bigint_value.columns, "columns", "",
					    "void base_matrix< bigmod >::"
					    "compose_v(const base_matrix< bigmod > &A, "
					    "const base_matrix< bigmod > &B)",
					    DMESSAGE, EMESSAGE[4]);

		switch(bigint_value.bitfield.get_representation()) {
		case matrix_flags::dense_representation:
			if (A.bigint_value.bitfield.get_representation() ==
			    matrix_flags::dense_representation)
				DD_base_modul_bigint.insert_at(bigint_value, 0, 0, A.bigint_value,
							       0, 0, A.bigint_value.rows,
							       A.bigint_value.columns);
			else
				DS_base_modul_bigint.insert_at(bigint_value, 0, 0, A.bigint_value,
							       0, 0, A.bigint_value.rows,
							       A.bigint_value.columns);

			if (B.bigint_value.bitfield.get_representation() ==
			    matrix_flags::dense_representation)
				DD_base_modul_bigint.insert_at(bigint_value, bigint_value.rows -
							       B.bigint_value.rows, 0,
							       B.bigint_value, 0, 0,
							       B.bigint_value.rows,
							       B.bigint_value.columns);
			else
				DS_base_modul_bigint.insert_at(bigint_value, bigint_value.rows -
							       B.bigint_value.rows, 0,
							       B.bigint_value, 0, 0,
							       B.bigint_value.rows,
							       B.bigint_value.columns);
			break;
		case matrix_flags::sparse_representation:
			if (A.bigint_value.bitfield.get_representation() ==
			    matrix_flags::dense_representation)
				SD_base_modul_bigint.insert_at(bigint_value, 0, 0,
							       A.bigint_value, 0, 0,
							       A.bigint_value.rows,
							       A.bigint_value.columns);
			else
				SS_base_modul_bigint.insert_at(bigint_value, 0, 0, A.bigint_value,
							       0, 0, A.bigint_value.rows,
							       A.bigint_value.columns);

			if (B.bigint_value.bitfield.get_representation() ==
			    matrix_flags::dense_representation)
				SD_base_modul_bigint.insert_at(bigint_value, bigint_value.rows -
							       B.bigint_value.rows, 0,
							       B.bigint_value, 0, 0,
							       B.bigint_value.rows,
							       B.bigint_value.columns);
			else
				SS_base_modul_bigint.insert_at(bigint_value, bigint_value.rows -
							       B.bigint_value.rows, 0,
							       B.bigint_value, 0, 0,
							       B.bigint_value.rows,
							       B.bigint_value.columns);
			break;
		}
	}
	else {
		if (A.long_value.rows > long_value.rows
		    || A.long_value.columns > long_value.columns
		    || B.long_value.rows > long_value.rows
		    || B.long_value.columns > long_value.columns)
			precondition_error_handler(A.long_value.rows, "A.rows", "A.rows <= rows",
					    A.long_value.columns, "A.columns",
					    "A.columns <= columns",
					    B.long_value.rows, "B.rows", "B.rows <= rows",
					    B.long_value.columns, "B.columns",
					    "B.columns <= columns",
					    long_value.rows, "rows", "",
					    long_value.columns, "columns", "",
					    "void base_matrix< bigmod >::"
					    "compose_v(const base_matrix< bigmod > &A, "
					    "const base_matrix< bigmod > &B)",
					    DMESSAGE, EMESSAGE[4]);

		switch(long_value.bitfield.get_representation()) {
		case matrix_flags::dense_representation:
			if (A.long_value.bitfield.get_representation() ==
			    matrix_flags::dense_representation)
				DD_base_modul_long.insert_at(long_value, 0, 0, A.long_value,
							     0, 0, A.long_value.rows,
							     A.long_value.columns);
			else
				DS_base_modul_long.insert_at(long_value, 0, 0, A.long_value,
							     0, 0, A.long_value.rows,
							     A.long_value.columns);

			if (B.long_value.bitfield.get_representation() ==
			    matrix_flags::dense_representation)
				DD_base_modul_long.insert_at(long_value, long_value.rows -
							     B.long_value.rows, 0,
							     B.long_value, 0, 0,
							     B.long_value.rows,
							     B.long_value.columns);
			else
				DS_base_modul_long.insert_at(long_value, long_value.rows -
							     B.long_value.rows, 0,
							     B.long_value, 0, 0,
							     B.long_value.rows,
							     B.long_value.columns);
			break;
		case matrix_flags::sparse_representation:
			if (A.long_value.bitfield.get_representation() ==
			    matrix_flags::dense_representation)
				SD_base_modul_long.insert_at(long_value, 0, 0,
							     A.long_value, 0, 0,
							     A.long_value.rows,
							     A.long_value.columns);
			else
				SS_base_modul_long.insert_at(long_value, 0, 0, A.long_value,
							     0, 0, A.long_value.rows,
							     A.long_value.columns);

			if (B.long_value.bitfield.get_representation() ==
			    matrix_flags::dense_representation)
				SD_base_modul_long.insert_at(long_value, long_value.rows -
							     B.long_value.rows, 0,
							     B.long_value, 0, 0,
							     B.long_value.rows,
							     B.long_value.columns);
			else
				SS_base_modul_long.insert_at(long_value, long_value.rows -
							     B.long_value.rows, 0,
							     B.long_value, 0, 0,
							     B.long_value.rows,
							     B.long_value.columns);
			break;
		}
	}
}



void base_matrix< bigmod >::
compose_v(const bigmod *v, const base_matrix< bigmod > &A)
{
	//
	//       Task: C.compose_v(v, A);
	//             => compose the vector v and submatrix A to matrix C
	//                    (v)
	//             => C = (A)
	// Conditions: A.columns <= columns and
	//             A.rows <= rows and
	//             v != NULL
	//    Version: 2.0
	//

	debug_handler_l(DMESSAGE, "compose_v(const bigmod *, "
			"const base_matrix< bigmod > &)", DVALUE + 3);

	if (BigintLong) {
		if (A.bigint_value.rows > bigint_value.rows
		    || A.bigint_value.columns > bigint_value.columns || v == NULL)
			precondition_error_handler(A.bigint_value.rows, "A.rows", "A.rows <= rows",
					    A.bigint_value.columns, "A.columns",
					    "A.columns <= columns",
					    bigint_value.rows, "rows", "",
					    bigint_value.columns, "columns", "",
					    PRT, "v", "v != NULL",
					    "void base_matrix< bigmod >::"
					    "compose_v(const bigmod *v, "
					    "const base_matrix< bigmod > &A)",
					    DMESSAGE, EMESSAGE[4]);

		bigint *tmp = new bigint[bigint_value.columns];
		memory_handler(tmp, DMESSAGE, "compose_v :: "
			       "Error in memory allocation (tmp)");
		lidia_size_t i;
		for (i = 0; i < bigint_value.columns; i++)
			tmp[i] = v[i].mantissa();

		switch(bigint_value.bitfield.get_representation()) {
		case matrix_flags::dense_representation:
			D_base_modul_bigint.insert_row_at(bigint_value, 0, tmp,
							  bigint_value.columns);
			if (A.bigint_value.bitfield.get_representation() ==
			    matrix_flags::dense_representation)
				DD_base_modul_bigint.insert_at(bigint_value, bigint_value.rows -
							       A.bigint_value.rows, 0,
							       A.bigint_value, 0, 0,
							       A.bigint_value.rows,
							       A.bigint_value.columns);
			else
				DS_base_modul_bigint.insert_at(bigint_value, bigint_value.rows -
							       A.bigint_value.rows, 0,
							       A.bigint_value, 0, 0,
							       A.bigint_value.rows,
							       A.bigint_value.columns);
			break;
		case matrix_flags::sparse_representation:
			S_base_modul_bigint.insert_row_at(bigint_value, 0, tmp,
							  bigint_value.columns);
			if (A.bigint_value.bitfield.get_representation() ==
			    matrix_flags::dense_representation)
				SD_base_modul_bigint.insert_at(bigint_value, bigint_value.rows -
							       A.bigint_value.rows, 0,
							       A.bigint_value, 0, 0,
							       A.bigint_value.rows,
							       A.bigint_value.columns);
			else
				SS_base_modul_bigint.insert_at(bigint_value, bigint_value.rows -
							       A.bigint_value.rows, 0,
							       A.bigint_value, 0, 0,
							       A.bigint_value.rows,
							       A.bigint_value.columns);
			break;
		}
		delete[] tmp;
	}
	else {
		if (A.long_value.rows > long_value.rows
		    || A.long_value.columns > long_value.columns || v == NULL)
			precondition_error_handler(A.long_value.rows, "A.rows", "A.rows <= rows",
					    A.long_value.columns, "A.columns",
					    "A.columns <= columns",
					    long_value.rows, "rows", "",
					    long_value.columns, "columns", "",
					    PRT, "v", "v != NULL",
					    "void base_matrix< bigmod >::"
					    "compose_v(const bigmod *v, "
					    "const base_matrix< bigmod > &A)",
					    DMESSAGE, EMESSAGE[4]);

		long *tmp = new long[long_value.columns];
		memory_handler(tmp, DMESSAGE, "compose_v :: "
			       "Error in memory allocation (tmp)");
		lidia_size_t i;
		for (i = 0; i < long_value.columns; i++)
			(v[i].mantissa()).longify(tmp[i]);

		switch(long_value.bitfield.get_representation()) {
		case matrix_flags::dense_representation:
			D_base_modul_long.insert_row_at(long_value, 0, tmp,
							long_value.columns);
			if (A.long_value.bitfield.get_representation() ==
			    matrix_flags::dense_representation)
				DD_base_modul_long.insert_at(long_value, long_value.rows -
							     A.long_value.rows, 0,
							     A.long_value, 0, 0,
							     A.long_value.rows,
							     A.long_value.columns);
			else
				DS_base_modul_long.insert_at(long_value, long_value.rows -
							     A.long_value.rows, 0,
							     A.long_value, 0, 0,
							     A.long_value.rows,
							     A.long_value.columns);
			break;
		case matrix_flags::sparse_representation:
			S_base_modul_long.insert_row_at(long_value, 0, tmp,
							long_value.columns);
			if (A.long_value.bitfield.get_representation() ==
			    matrix_flags::dense_representation)
				SD_base_modul_long.insert_at(long_value, long_value.rows -
							     A.long_value.rows, 0,
							     A.long_value, 0, 0,
							     A.long_value.rows,
							     A.long_value.columns);
			else
				SS_base_modul_long.insert_at(long_value, long_value.rows -
							     A.long_value.rows, 0,
							     A.long_value, 0, 0,
							     A.long_value.rows,
							     A.long_value.columns);
			break;
		}
		delete[] tmp;
	}
}



void base_matrix< bigmod >::
compose_v(const base_matrix< bigmod > &A, const bigmod *v)
{
	//
	//       Task: C.compose_v(A, v);
	//             => compose the submatrix A and the vector v to matrix C
	//                    (A)
	//             => C = (v)
	// Conditions: A.columns <= columns and
	//             A.rows <= rows and
	//             v != NULL
	//    Version: 2.0
	//

	debug_handler_l(DMESSAGE, "compose_v(const base_matrix< bigmod > &, "
			"const bigmod *)", DVALUE + 3);

	if (BigintLong) {
		if (A.bigint_value.rows > bigint_value.rows
		    || A.bigint_value.columns > bigint_value.columns || v == NULL)
			precondition_error_handler(A.bigint_value.rows, "A.rows", "A.rows <= rows",
					    A.bigint_value.columns, "A.columns",
					    "A.columns <= columns",
					    bigint_value.rows, "rows", "",
					    bigint_value.columns, "columns", "",
					    PRT, "v", "v != NULL",
					    "void base_matrix< bigmod >::"
					    "compose_v(const base_matrix< bigmod > &A, "
					    "const bigmod *v)", DMESSAGE, EMESSAGE[4]);

		bigint *tmp = new bigint[bigint_value.columns];
		memory_handler(tmp, DMESSAGE, "compose_v :: "
			       "Error in memory allocation (tmp)");
		lidia_size_t i;
		for (i = 0; i < bigint_value.columns; i++)
			tmp[i] = v[i].mantissa();

		switch(bigint_value.bitfield.get_representation()) {
		case matrix_flags::dense_representation:
			if (A.bigint_value.bitfield.get_representation() ==
			    matrix_flags::dense_representation)
				DD_base_modul_bigint.insert_at(bigint_value, 0, 0, A.bigint_value,
							       0, 0, A.bigint_value.rows,
							       A.bigint_value.columns);
			else
				DS_base_modul_bigint.insert_at(bigint_value, 0, 0, A.bigint_value,
							       0, 0, A.bigint_value.rows,
							       A.bigint_value.columns);
			D_base_modul_bigint.insert_row_at(bigint_value,
							  bigint_value.rows - 1, tmp,
							  bigint_value.columns);
			break;
		case matrix_flags::sparse_representation:
			if (A.bigint_value.bitfield.get_representation() ==
			    matrix_flags::dense_representation)
				SD_base_modul_bigint.insert_at(bigint_value, 0, 0, A.bigint_value,
							       0, 0, A.bigint_value.rows,
							       A.bigint_value.columns);
			else
				SS_base_modul_bigint.insert_at(bigint_value, 0, 0, A.bigint_value,
							       0, 0, A.bigint_value.rows,
							       A.bigint_value.columns);
			S_base_modul_bigint.insert_row_at(bigint_value,
							  bigint_value.rows - 1, tmp,
							  bigint_value.columns);
			break;
		}
		delete[] tmp;
	}
	else {
		if (A.long_value.rows > long_value.rows
		    || A.long_value.columns > long_value.columns || v == NULL)
			precondition_error_handler(A.long_value.rows, "A.rows", "A.rows <= rows",
					    A.long_value.columns, "A.columns",
					    "A.columns <= columns",
					    long_value.rows, "rows", "",
					    long_value.columns, "columns", "",
					    PRT, "v", "v != NULL",
					    "void base_matrix< bigmod >::"
					    "compose_v(const base_matrix< bigmod > &A, "
					    "const bigmod *v)", DMESSAGE, EMESSAGE[4]);

		long *tmp = new long[long_value.columns];
		memory_handler(tmp, DMESSAGE, "compose_v :: "
			       "Error in memory allocation (tmp)");
		lidia_size_t i;
		for (i = 0; i < long_value.columns; i++)
			(v[i].mantissa()).longify(tmp[i]);

		switch(long_value.bitfield.get_representation()) {
		case matrix_flags::dense_representation:
			if (A.long_value.bitfield.get_representation() ==
			    matrix_flags::dense_representation)
				DD_base_modul_long.insert_at(long_value, 0, 0, A.long_value,
							     0, 0, A.long_value.rows,
							     A.long_value.columns);
			else
				DS_base_modul_long.insert_at(long_value, 0, 0, A.long_value,
							     0, 0, A.long_value.rows,
							     A.long_value.columns);
			D_base_modul_long.insert_row_at(long_value,
							long_value.rows - 1, tmp,
							long_value.columns);
			break;
		case matrix_flags::sparse_representation:
			if (A.long_value.bitfield.get_representation() ==
			    matrix_flags::dense_representation)
				SD_base_modul_long.insert_at(long_value, 0, 0, A.long_value,
							     0, 0, A.long_value.rows,
							     A.long_value.columns);
			else
				SS_base_modul_long.insert_at(long_value, 0, 0, A.long_value,
							     0, 0, A.long_value.rows,
							     A.long_value.columns);
			S_base_modul_long.insert_row_at(long_value,
							long_value.rows - 1, tmp,
							long_value.columns);
			break;
		}
		delete[] tmp;
	}
}



void base_matrix< bigmod >::
compose_v(const base_vector< bigmod > &v, const base_matrix< bigmod > &A)
{
	//
	//       Task: C.compose_v(v, A);
	//             => compose the vector v and submatrix A to matrix C
	//                    (v)
	//             => C = (A)
	// Conditions: A.columns <= columns and
	//             A.rows <= rows and
	//             v.size <= columns
	//    Version: 2.0
	//

	debug_handler_l(DMESSAGE, "compose_v(const base_vector< bigmod > &, "
			"const base_matrix< bigmod > &)", DVALUE + 3);

	lidia_size_t l = v.size();

	if (BigintLong) {
		if (A.bigint_value.rows > bigint_value.rows
		    || A.bigint_value.columns > bigint_value.columns
		    || l > bigint_value.columns)
			precondition_error_handler(A.bigint_value.rows, "A.rows", "A.rows <= rows",
					    A.bigint_value.columns, "A.columns",
					    "A.columns <= columns",
					    bigint_value.rows, "rows", "",
					    bigint_value.columns, "columns", "",
					    l, "v.size", "v.size <= coluns",
					    "void base_matrix< bigmod >::"
					    "compose_v(const base_vector< bigmod > &v, "
					    "const base_matrix< bigmod > &A)",
					    DMESSAGE, EMESSAGE[4]);

		bigint *tmp = new bigint[bigint_value.columns];
		memory_handler(tmp, DMESSAGE, "compose_v :: "
			       "Error in memory allocation (tmp)");
		lidia_size_t i;
		for (i = 0; i < bigint_value.columns; i++)
			tmp[i] = v[i].mantissa();

		switch(bigint_value.bitfield.get_representation()) {
		case matrix_flags::dense_representation:
			D_base_modul_bigint.insert_row_at(bigint_value, 0, tmp, l);
			if (A.bigint_value.bitfield.get_representation() ==
			    matrix_flags::dense_representation)
				DD_base_modul_bigint.insert_at(bigint_value, bigint_value.rows -
							       A.bigint_value.rows, 0,
							       A.bigint_value, 0, 0,
							       A.bigint_value.rows,
							       A.bigint_value.columns);
			else
				DS_base_modul_bigint.insert_at(bigint_value, bigint_value.rows -
							       A.bigint_value.rows, 0,
							       A.bigint_value, 0, 0,
							       A.bigint_value.rows,
							       A.bigint_value.columns);
			break;
		case matrix_flags::sparse_representation:
			S_base_modul_bigint.insert_row_at(bigint_value, 0, tmp, l);
			if (A.bigint_value.bitfield.get_representation() ==
			    matrix_flags::dense_representation)
				SD_base_modul_bigint.insert_at(bigint_value, bigint_value.rows -
							       A.bigint_value.rows, 0,
							       A.bigint_value, 0, 0,
							       A.bigint_value.rows,
							       A.bigint_value.columns);
			else
				SS_base_modul_bigint.insert_at(bigint_value, bigint_value.rows -
							       A.bigint_value.rows, 0,
							       A.bigint_value, 0, 0,
							       A.bigint_value.rows,
							       A.bigint_value.columns);
			break;
		}
		delete[] tmp;
	}
	else {
		if (A.long_value.rows > long_value.rows
		    || A.long_value.columns > long_value.columns
		    || l > long_value.columns)
			precondition_error_handler(A.long_value.rows, "A.rows", "A.rows <= rows",
					    A.long_value.columns, "A.columns",
					    "A.columns <= columns",
					    long_value.rows, "rows", "",
					    long_value.columns, "columns", "",
					    l, "v.size", "v.size <= coluns",
					    "void base_matrix< bigmod >::"
					    "compose_v(const base_vector< bigmod > &v, "
					    "const base_matrix< bigmod > &A)",
					    DMESSAGE, EMESSAGE[4]);

		long *tmp = new long[long_value.columns];
		memory_handler(tmp, DMESSAGE, "compose_v :: "
			       "Error in memory allocation (tmp)");
		lidia_size_t i;
		for (i = 0; i < long_value.columns; i++)
			(v[i].mantissa()).longify(tmp[i]);

		switch(long_value.bitfield.get_representation()) {
		case matrix_flags::dense_representation:
			D_base_modul_long.insert_row_at(long_value, 0, tmp, l);
			if (A.long_value.bitfield.get_representation() ==
			    matrix_flags::dense_representation)
				DD_base_modul_long.insert_at(long_value, long_value.rows -
							     A.long_value.rows, 0,
							     A.long_value, 0, 0,
							     A.long_value.rows,
							     A.long_value.columns);
			else
				DS_base_modul_long.insert_at(long_value, long_value.rows -
							     A.long_value.rows, 0,
							     A.long_value, 0, 0,
							     A.long_value.rows,
							     A.long_value.columns);
			break;
		case matrix_flags::sparse_representation:
			S_base_modul_long.insert_row_at(long_value, 0, tmp, l);
			if (A.long_value.bitfield.get_representation() ==
			    matrix_flags::dense_representation)
				SD_base_modul_long.insert_at(long_value, long_value.rows -
							     A.long_value.rows, 0,
							     A.long_value, 0, 0,
							     A.long_value.rows,
							     A.long_value.columns);
			else
				SS_base_modul_long.insert_at(long_value, long_value.rows -
							     A.long_value.rows, 0,
							     A.long_value, 0, 0,
							     A.long_value.rows,
							     A.long_value.columns);
			break;
		}
		delete[] tmp;
	}
}



void base_matrix< bigmod >::
compose_v(const base_matrix< bigmod > &A, const base_vector< bigmod > &v)
{
	//
	//       Task: C.compose_v(A, v);
	//             => compose the submatrix A and the vector v to matrix C
	//                    (A)
	//             => C = (v)
	// Conditions: A.columns <= columns and
	//             A.rows <= rows and
	//             v.size <= columns
	//    Version: 2.0
	//

	debug_handler_l(DMESSAGE, "compose_v(const base_matrix< bigmod > &, "
			"const base_vector< bigmod > &)", DVALUE + 3);

	lidia_size_t l = v.size();

	if (BigintLong) {
		if (A.bigint_value.rows > bigint_value.rows
		    || A.bigint_value.columns > bigint_value.columns
		    || l > bigint_value.columns)
			precondition_error_handler(A.bigint_value.rows, "A.rows", "A.rows <= rows",
					    A.bigint_value.columns, "A.columns",
					    "A.columns <= columns",
					    bigint_value.rows, "rows", "",
					    bigint_value.columns, "columns", "",
					    l, "v.size", "v.size <= columns",
					    "void base_matrix< bigmod >::"
					    "compose_v(const base_matrix< bigmod > &A, "
					    "const base_vector< bigmod > &v)",
					    DMESSAGE, EMESSAGE[4]);

		bigint *tmp = new bigint[bigint_value.columns];
		memory_handler(tmp, DMESSAGE, "compose_v :: "
			       "Error in memory allocation (tmp)");
		lidia_size_t i;
		for (i = 0; i < bigint_value.columns; i++)
			tmp[i] = v[i].mantissa();

		switch(bigint_value.bitfield.get_representation()) {
		case matrix_flags::dense_representation:
			if (A.bigint_value.bitfield.get_representation() ==
			    matrix_flags::dense_representation)
				DD_base_modul_bigint.insert_at(bigint_value, 0, 0, A.bigint_value,
							       0, 0, A.bigint_value.rows,
							       A.bigint_value.columns);
			else
				DS_base_modul_bigint.insert_at(bigint_value, 0, 0, A.bigint_value,
							       0, 0, A.bigint_value.rows,
							       A.bigint_value.columns);
			D_base_modul_bigint.insert_row_at(bigint_value,
							  bigint_value.rows - 1, tmp, l);
			break;
		case matrix_flags::sparse_representation:
			if (A.bigint_value.bitfield.get_representation() ==
			    matrix_flags::dense_representation)
				SD_base_modul_bigint.insert_at(bigint_value, 0, 0, A.bigint_value,
							       0, 0, A.bigint_value.rows,
							       A.bigint_value.columns);
			else
				SS_base_modul_bigint.insert_at(bigint_value, 0, 0, A.bigint_value,
							       0, 0, A.bigint_value.rows,
							       A.bigint_value.columns);
			S_base_modul_bigint.insert_row_at(bigint_value,
							  bigint_value.rows - 1, tmp, l);
			break;
		}
		delete[] tmp;
	}
	else {
		if (A.long_value.rows > long_value.rows
		    || A.long_value.columns > long_value.columns
		    || l > long_value.columns)
			precondition_error_handler(A.long_value.rows, "A.rows", "A.rows <= rows",
					    A.long_value.columns, "A.columns",
					    "A.columns <= columns",
					    long_value.rows, "rows", "",
					    long_value.columns, "columns", "",
					    l, "v.size", "v.size <= columns",
					    "void base_matrix< bigmod >::"
					    "compose_v(const base_matrix< bigmod > &A, "
					    "const base_vector< bigmod > &v)",
					    DMESSAGE, EMESSAGE[4]);

		long *tmp = new long[long_value.columns];
		memory_handler(tmp, DMESSAGE, "compose_v :: "
			       "Error in memory allocation (tmp)");
		lidia_size_t i;
		for (i = 0; i < long_value.columns; i++)
			(v[i].mantissa()).longify(tmp[i]);

		switch(long_value.bitfield.get_representation()) {
		case matrix_flags::dense_representation:
			if (A.long_value.bitfield.get_representation() ==
			    matrix_flags::dense_representation)
				DD_base_modul_long.insert_at(long_value, 0, 0, A.long_value,
							     0, 0, A.long_value.rows,
							     A.long_value.columns);
			else
				DS_base_modul_long.insert_at(long_value, 0, 0, A.long_value,
							     0, 0, A.long_value.rows,
							     A.long_value.columns);
			D_base_modul_long.insert_row_at(long_value,
							long_value.rows - 1, tmp, l);
			break;
		case matrix_flags::sparse_representation:
			if (A.long_value.bitfield.get_representation() ==
			    matrix_flags::dense_representation)
				SD_base_modul_long.insert_at(long_value, 0, 0, A.long_value,
							     0, 0, A.long_value.rows,
							     A.long_value.columns);
			else
				SS_base_modul_long.insert_at(long_value, 0, 0, A.long_value,
							     0, 0, A.long_value.rows,
							     A.long_value.columns);
			S_base_modul_long.insert_row_at(long_value,
							long_value.rows - 1, tmp, l);
			break;
		}
		delete[] tmp;
	}
}



//
// exchange functions / swap functions
//

void base_matrix< bigmod >::
swap(base_matrix< bigmod > &B)
{
	//
	//    Task: A.swap(B) exchanges matrix A with matrix B.
	// Version: 2.0
	//

	debug_handler_l(DMESSAGE, "swap(base_matrix< bigmod > &)", DVALUE + 4);

	LiDIA::swap(bigint_modulus, B.bigint_modulus);
	LiDIA::swap(long_modulus, B.long_modulus);

	bool SW = BigintLong;
	BigintLong = B.BigintLong;
	B.BigintLong = SW;

	//////////////////////////
	// bigint part
	/////////////////////////

	// swap no_of_rows
	LiDIA::swap(bigint_value.rows, B.bigint_value.rows);

	// swap no_of_columns
	LiDIA::swap(bigint_value.columns, B.bigint_value.columns);

	// swap no_of_sparse_rows
	LiDIA::swap(bigint_value.sparse_rows, B.bigint_value.sparse_rows);

	// swap no_of_sparse_columns
	LiDIA::swap(bigint_value.sparse_columns, B.bigint_value.sparse_columns);

	// swap value_counter
	lidia_size_t *tmp2 = bigint_value.value_counter;
	bigint_value.value_counter = B.bigint_value.value_counter;
	B.bigint_value.value_counter = tmp2;

	// swap allocated
	tmp2 = bigint_value.allocated;
	bigint_value.allocated = B.bigint_value.allocated;
	B.bigint_value.allocated = tmp2;

	// swap values
	bigint **tmp = bigint_value.value;
	bigint_value.value = B.bigint_value.value;
	B.bigint_value.value = tmp;

	// swap index
	lidia_size_t **tmp1 = bigint_value.index;
	bigint_value.index = B.bigint_value.index;
	B.bigint_value.index = tmp1;

	//////////////////////////
	// long part
	/////////////////////////

	// swap no_of_rows
	LiDIA::swap(long_value.rows, B.long_value.rows);

	// swap no_of_columns
	LiDIA::swap(long_value.columns, B.long_value.columns);

	// swap no_of_sparse_rows
	LiDIA::swap(long_value.sparse_rows, B.long_value.sparse_rows);

	// swap no_of_sparse_columns
	LiDIA::swap(long_value.sparse_columns, B.long_value.sparse_columns);

	// swap value_counter
	tmp2 = long_value.value_counter;
	long_value.value_counter = B.long_value.value_counter;
	B.long_value.value_counter = tmp2;

	// swap allocated
	tmp2 = long_value.allocated;
	long_value.allocated = B.long_value.allocated;
	B.long_value.allocated = tmp2;

	// swap values
	long **tmp3 = long_value.value;
	long_value.value = B.long_value.value;
	B.long_value.value = tmp3;

	// swap index
	tmp1 = long_value.index;
	long_value.index = B.long_value.index;
	B.long_value.index = tmp1;

	// swap bitfields
	LiDIA::swap(bigint_value.bitfield, B.bigint_value.bitfield);
	LiDIA::swap(long_value.bitfield, B.long_value.bitfield);
}



void base_matrix< bigmod >::
swap_columns(lidia_size_t i, lidia_size_t j)
{
	//
	//       Task: A.swap_columns(i,j) exchanges column i with column j
	//             in matrix A.
	// Conditions: 0 <= i < A.columns and
	//             0 <= j < A.columns
	//    Version: 2.0
	//

	debug_handler_l(DMESSAGE, "swap_columns(lidia_size_t , lidia_size_t)",
			DVALUE + 4);

	if (BigintLong) {
		if (i< 0 || i >= bigint_value.columns
		    || j< 0 || j >= bigint_value.columns)
			precondition_error_handler(i, "i", "0 <= i < columns",
					    j, "j", "0 <= j < columns",
					    bigint_value.columns, "columns", "",
					    "void base_matrix< bigmod >::"
					    "swap_columns(lidia_size_t i, lidia_size_t j)",
					    DMESSAGE, EMESSAGE[1]);

		switch(bigint_value.bitfield.get_representation()) {
		case matrix_flags::dense_representation:
			D_base_modul_bigint.swap_columns(bigint_value, i, j);
			break;
		case matrix_flags::sparse_representation:
			SD_base_modul_bigint.swap_columns(bigint_value, i, j);
			break;
		}
	}
	else {
		if (i< 0 || i >= long_value.columns
		    || j< 0 || j >= long_value.columns)
			precondition_error_handler(i, "i", "0 <= i < columns",
					    j, "j", "0 <= j < columns",
					    long_value.columns, "columns", "",
					    "void base_matrix< bigmod >::"
					    "swap_columns(lidia_size_t i, lidia_size_t j)",
					    DMESSAGE, EMESSAGE[1]);

		switch(long_value.bitfield.get_representation()) {
		case matrix_flags::dense_representation:
			D_base_modul_long.swap_columns(long_value, i, j);
			break;
		case matrix_flags::sparse_representation:
			SD_base_modul_long.swap_columns(long_value, i, j);
			break;
		}
	}
}



void base_matrix< bigmod >::
swap_rows(lidia_size_t i, lidia_size_t j)
{
	//
	//       Task: A.swap_rows(i,j) exchanges row i with row j
	//             in matrix A.
	// Conditions: 0 <= i < A.rows and
	//             0 <= j < A.rows
	//    Version: 2.0
	//

	debug_handler_l(DMESSAGE, "swap_rows(lidia_size_t , lidia_size_t)",
			DVALUE + 4);

	if (BigintLong) {
		if (i< 0 || i >= bigint_value.rows
		    || j< 0 || j >= bigint_value.rows)
			precondition_error_handler(i, "i", "0 <= i < rows",
					    j, "j", "0 <= j < rows",
					    bigint_value.rows, "rows", "",
					    "void base_matrix< bigmod >::"
					    "swap_rows(lidia_size_t i, lidia_size_t j)",
					    DMESSAGE, EMESSAGE[1]);

		switch(bigint_value.bitfield.get_representation()) {
		case matrix_flags::dense_representation:
			D_base_modul_bigint.swap_rows(bigint_value, i, j);
			break;
		case matrix_flags::sparse_representation:
			SD_base_modul_bigint.swap_rows(bigint_value, i, j);
			break;
		}
	}
	else {
		if (i< 0 || i >= long_value.rows
		    || j< 0 || j >= long_value.rows)
			precondition_error_handler(i, "i", "0 <= i < rows",
					    j, "j", "0 <= j < rows",
					    long_value.rows, "rows", "",
					    "void base_matrix< bigmod >::"
					    "swap_rows(lidia_size_t i, lidia_size_t j)",
					    DMESSAGE, EMESSAGE[1]);

		switch(long_value.bitfield.get_representation()) {
		case matrix_flags::dense_representation:
			D_base_modul_long.swap_rows(long_value, i, j);
			break;
		case matrix_flags::sparse_representation:
			SD_base_modul_long.swap_rows(long_value, i, j);
			break;
		}
	}
}



//
// structur functions
//

void base_matrix< bigmod >::
set_no_of_rows(lidia_size_t r)
{
	//
	//       Task: A.set_no_of_rows(r) sets number of rows of matrix A to r.
	// Conditions: r >= 0
	//    Version: 2.0
	//

	debug_handler_l(DMESSAGE, "set_no_of_rows(lidia_size_t)", DVALUE + 5);

	if (r < 0)
		precondition_error_handler(r, "r", "r >= 0",
				    "void base_matrix< bigmod >::"
				    "set_no_of_rows(lidia_size_t r)",
				    DMESSAGE, EMESSAGE[1]);

	if (BigintLong) {
		switch(bigint_value.bitfield.get_representation()) {
		case matrix_flags::dense_representation:
			D_base_modul_bigint.set_no_of_rows(bigint_value, r);
			break;
		case matrix_flags::sparse_representation:
			S_base_modul_bigint.set_no_of_rows(bigint_value, r);
			break;
		}
	}
	else {
		switch(long_value.bitfield.get_representation()) {
		case matrix_flags::dense_representation:
			D_base_modul_long.set_no_of_rows(long_value, r);
			break;
		case matrix_flags::sparse_representation:
			S_base_modul_long.set_no_of_rows(long_value, r);
			break;
		}
	}
}



void base_matrix< bigmod >::
set_no_of_columns(lidia_size_t c)
{
	//
	//       Task: A.set_no_of_columns(c) sets number of
	//             columns of matrix A to c.
	// Conditions: c >= 0
	//    Version: 2.0
	//

	debug_handler_l(DMESSAGE, "set_no_of_columns(lidia_size_t)", DVALUE + 5);

	if (c < 0)
		precondition_error_handler(c, "c", "c >= 0",
				    "void base_matrix< bigmod >::"
				    "set_no_of_columns(lidia_size_t c)",
				    DMESSAGE, EMESSAGE[1]);

	if (BigintLong) {
		switch(bigint_value.bitfield.get_representation()) {
		case matrix_flags::dense_representation:
			D_base_modul_bigint.set_no_of_columns(bigint_value, c);
			break;
		case matrix_flags::sparse_representation:
			S_base_modul_bigint.set_no_of_columns(bigint_value, c);
			break;
		}
	}
	else {
		switch(long_value.bitfield.get_representation()) {
		case matrix_flags::dense_representation:
			D_base_modul_long.set_no_of_columns(long_value, c);
			break;
		case matrix_flags::sparse_representation:
			S_base_modul_long.set_no_of_columns(long_value, c);
			break;
		}
	}
	debug_handler_c(DMESSAGE, "set_no_of_columns(lidia_size_t)", DVALUE + 5,
			status_report());
}



void base_matrix< bigmod >::
resize(lidia_size_t r, lidia_size_t c)
{
	//
	//       Task: A.resize(r,c)
	//             => sets number of columns of matrix A to c and
	//                number of rows of matrix A to r.
	// Conditions: r >= 0 and
	//             c >= 0
	//    Version: 2.0
	//

	debug_handler_l(DMESSAGE, "resize(lidia_size_t, lidia_size_t)", DVALUE + 5);

	if (c < 0 || r < 0)
		precondition_error_handler(r, "r", "r >= 0",
				    c, "c", "c >= 0",
				    "void base_matrix< bigmod >::"
				    "resize(lidia_size_t r, lidia_size_t c)",
				    DMESSAGE, EMESSAGE[1]);

	if (BigintLong) {
		switch(bigint_value.bitfield.get_representation()) {
		case matrix_flags::dense_representation:
			DD_base_modul_bigint.resize(bigint_value, r, c);
			break;
		case matrix_flags::sparse_representation:
			SD_base_modul_bigint.resize(bigint_value, r, c);
			break;
		}
	}
	else {
		switch(long_value.bitfield.get_representation()) {
		case matrix_flags::dense_representation:
			DD_base_modul_long.resize(long_value, r, c);
			break;
		case matrix_flags::sparse_representation:
			SD_base_modul_long.resize(long_value, r, c);
			break;
		}
	}

	debug_handler_c(DMESSAGE, "resize(lidia_size_t, lidia_size_t)", DVALUE + 5,
			status_report());
}



void base_matrix< bigmod >::
kill()
{
	//
	//    Task: A.kill()
	//          => sets number of columns of matrix A to 0 and
	//             number of rows of matrix A to 0.
	// Version: 2.0
	//

	debug_handler_l(DMESSAGE, "kill()", DVALUE + 5);

	if (BigintLong) {
		switch(bigint_value.bitfield.get_representation()) {
		case matrix_flags::dense_representation:
			D_base_modul_bigint.kill(bigint_value);
			break;
		case matrix_flags::sparse_representation:
			SD_base_modul_bigint.kill(bigint_value);
			break;
		}
	}
	else {
		switch(long_value.bitfield.get_representation()) {
		case matrix_flags::dense_representation:
			D_base_modul_long.kill(long_value);
			break;
		case matrix_flags::sparse_representation:
			SD_base_modul_long.kill(long_value);
			break;
		}
	}

	debug_handler_c(DMESSAGE, "kill()", DVALUE + 5, status_report());
}



//
// assignments
//

void base_matrix< bigmod >::
assign(const base_matrix< bigmod > &M)
{
	//
	//    Task: A.assign(B);
	//          => A.value[x][y] = B.value[x][y],
	//             x=0,...,A.rows-1, y=0,...,A.columns-1
	//          => A.rows = B.rows and A.columns = B.columns
	// Version: 2.0
	//

	debug_handler_l(DMESSAGE, "assign(const base_matrix< bigmod > &)",
			DVALUE + 6);

	BigintLong = M.BigintLong;

	if (BigintLong) {
		set_representation(M.bigint_value.bitfield.get_representation());

		switch(bigint_value.bitfield.get_representation()) {
		case matrix_flags::dense_representation:
			if (bigint_value.rows != M.bigint_value.rows)
				D_base_modul_bigint.set_no_of_rows(bigint_value,
								   M.bigint_value.rows);
			if (bigint_value.columns != M.bigint_value.columns)
				D_base_modul_bigint.set_no_of_columns(bigint_value,
								      M.bigint_value.columns);

			DD_base_modul_bigint.assign(bigint_value, M.bigint_value);
			break;
		case matrix_flags::sparse_representation:
			if (bigint_value.rows != M.bigint_value.rows)
				S_base_modul_bigint.set_no_of_rows(bigint_value,
								   M.bigint_value.rows);
			if (bigint_value.columns != M.bigint_value.columns)
				S_base_modul_bigint.set_no_of_columns(bigint_value,
								      M.bigint_value.columns);

			SS_base_modul_bigint.assign(bigint_value, M.bigint_value);
			break;
		}
	}
	else {
		set_representation(M.long_value.bitfield.get_representation());

		switch(long_value.bitfield.get_representation()) {
		case matrix_flags::dense_representation:
			if (long_value.rows != M.long_value.rows)
				D_base_modul_long.set_no_of_rows(long_value,
								 M.long_value.rows);
			if (long_value.columns != M.long_value.columns)
				D_base_modul_long.set_no_of_columns(long_value,
								    M.long_value.columns);

			DD_base_modul_long.assign(long_value, M.long_value);
			break;
		case matrix_flags::sparse_representation:
			if (long_value.rows != M.long_value.rows)
				S_base_modul_long.set_no_of_rows(long_value,
								 M.long_value.rows);
			if (long_value.columns != M.long_value.columns)
				S_base_modul_long.set_no_of_columns(long_value,
								    M.long_value.columns);

			SS_base_modul_long.assign(long_value, M.long_value);
			break;
		}
	}
}



//
// diagonal function
//

void base_matrix< bigmod >::
diag(const bigmod &a, const bigmod &b)
{
	//
	//    Task: A.diag(a,b);
	//          => A.value[i][i] = a, i=0,...,min(columns,rows)
	//          => A.value[i][j] = b, i=0,...,rows and j=0,...,columns
	//             and i != j
	// Version: 2.0
	//

	debug_handler_l(DMESSAGE, "diag(const bigmod &, const bigmod &)",
			DVALUE + 7);

	if (BigintLong) {
		bigint a1, b1;
		LiDIA::best_remainder(a1, a.mantissa(), bigint_modulus);
		LiDIA::best_remainder(b1, b.mantissa(), bigint_modulus);

		switch(bigint_value.bitfield.get_representation()) {
		case matrix_flags::dense_representation:
			D_base_modul_bigint.diag(bigint_value, a1, b1);
			break;
		case matrix_flags::sparse_representation:
			SS_base_modul_bigint.diag(bigint_value, a1, b1);
			break;
		}
	}
	else {
		long a1, b1;
		LiDIA::best_remainder(a1, a.mantissa(), long_modulus);
		LiDIA::best_remainder(b1, b.mantissa(), long_modulus);

		switch(long_value.bitfield.get_representation()) {
		case matrix_flags::dense_representation:
			D_base_modul_long.diag(long_value, a1, b1);
			break;
		case matrix_flags::sparse_representation:
			SS_base_modul_long.diag(long_value, a1, b1);
			break;
		}
	}
}



//
// transpose function
//

base_matrix< bigmod > base_matrix< bigmod >::
trans() const
{
	//
	//    Task: AT = A.trans();
	//          => AT.value[x][y] = A.value[y][x],
	//             x=0,...,A.columns-1, y=0,...,A.rows-1
	// Version: 2.0
	//

	debug_handler_l(DMESSAGE, "trans()", DVALUE + 7);

	if (BigintLong) {
		base_matrix< bigmod > TRANS(bigint_value.columns, bigint_value.rows,
					    bigint_modulus, bigint_value.bitfield);

		switch(bigint_value.bitfield.get_representation()) {
		case matrix_flags::dense_representation:
			DD_base_modul_bigint.trans(TRANS.bigint_value, bigint_value);
			break;
		case matrix_flags::sparse_representation:
			SS_base_modul_bigint.trans(TRANS.bigint_value, bigint_value);
			break;
		}
		return TRANS;
	}
	else {
		base_matrix< bigmod > TRANS(long_value.columns, long_value.rows,
					    bigint(long_modulus), long_value.bitfield);

		switch(long_value.bitfield.get_representation()) {
		case matrix_flags::dense_representation:
			DD_base_modul_long.trans(TRANS.long_value, long_value);
			break;
		case matrix_flags::sparse_representation:
			SS_base_modul_long.trans(TRANS.long_value, long_value);
			break;
		}
		return TRANS;
	}
}



void base_matrix< bigmod >::
trans(const base_matrix< bigmod > &A)
{
	//
	//    Task: AT.trans(A);
	//          => AT.value[x][y] = A.value[y][x],
	//             x=0,...,A.columns-1, y=0,...,A.rows-1
	//          => AT.rows = A.columns and AT.columns = A.rows
	// Version: 2.0
	//

	debug_handler_l(DMESSAGE,
			"trans(const base_matrix< bigmod > &)", DVALUE + 7);

	if (BigintLong) {
		set_representation(A.bigint_value.bitfield.get_representation());

		switch(bigint_value.bitfield.get_representation()) {
		case matrix_flags::dense_representation:
			DD_base_modul_bigint.trans(bigint_value, A.bigint_value);
			break;
		case matrix_flags::sparse_representation:
			SS_base_modul_bigint.trans(bigint_value, A.bigint_value);
			break;
		}
	}
	else {
		set_representation(A.long_value.bitfield.get_representation());

		switch(long_value.bitfield.get_representation()) {
		case matrix_flags::dense_representation:
			DD_base_modul_long.trans(long_value, A.long_value);
			break;
		case matrix_flags::sparse_representation:
			SS_base_modul_long.trans(long_value, A.long_value);
			break;
		}
	}
}



//
// stream handling - LiDIA
//

void base_matrix< bigmod >::
write_to_beauty(std::ostream &out) const
{
	//
	//    Task: A.write_to_stream(out);
	//          => writes matrix A to stream out in beauty - format.
	// Version: 2.0
	//

	debug_handler_l(DMESSAGE, "write_to_beauty(std::ostream &)", DVALUE + 8);

	if (BigintLong) {
		switch(bigint_value.bitfield.get_representation()) {
		case matrix_flags::dense_representation:
			D_base_modul_bigint.write_to_beauty(bigint_value, out);
			break;
		case matrix_flags::sparse_representation:
			S_base_modul_bigint.write_to_beauty(bigint_value, out);
			break;
		}
	}
	else {
		switch(long_value.bitfield.get_representation()) {
		case matrix_flags::dense_representation:
			D_base_modul_long.write_to_beauty(long_value, out);
			break;
		case matrix_flags::sparse_representation:
			S_base_modul_long.write_to_beauty(long_value, out);
			break;
		}
	}
}



void base_matrix< bigmod >::
write_to_stream(std::ostream &out) const
{
	//
	//    Task: A.write_to_stream(out);
	//          => writes matrix A to stream out in LiDIA - format.
	// Version: 2.0
	//

	debug_handler_l(DMESSAGE, "write_to_stream(std::ostream &)", DVALUE + 8);

	if (BigintLong) {
		switch(bigint_value.bitfield.get_representation()) {
		case matrix_flags::dense_representation:
			D_base_modul_bigint.write_to_stream(bigint_value, out);
			break;
		case matrix_flags::sparse_representation:
			S_base_modul_bigint.write_to_stream(bigint_value, out);
			break;
		}
	}
	else {
		switch(long_value.bitfield.get_representation()) {
		case matrix_flags::dense_representation:
			D_base_modul_long.write_to_stream(long_value, out);
			break;
		case matrix_flags::sparse_representation:
			S_base_modul_long.write_to_stream(long_value, out);
			break;
		}
	}
}



void base_matrix< bigmod >::
read_from_stream(std::istream &in)
{
	//
	//    Task: A.read_from_stream(in);
	//          => reads matrix A in LiDIA - format from stream in .
	// Version: 2.0
	//

	debug_handler_l(DMESSAGE, "read_from_stream(std::istream &)", DVALUE + 8);

	if (BigintLong) {
		switch(bigint_value.bitfield.get_representation()) {
		case matrix_flags::dense_representation:
			D_base_modul_bigint.read_from_stream(bigint_value, in);
			break;
		case matrix_flags::sparse_representation:
			S_base_modul_bigint.read_from_stream(bigint_value, in);
			break;
		}
	}
	else {
		switch(long_value.bitfield.get_representation()) {
		case matrix_flags::dense_representation:
			D_base_modul_long.read_from_stream(long_value, in);
			break;
		case matrix_flags::sparse_representation:
			S_base_modul_long.read_from_stream(long_value, in);
			break;
		}
	}
}



//
// stream handling - MATHEMATICA
//

void base_matrix< bigmod >::
write_to_mathematica(std::ostream &out) const
{
	//
	//    Task: A.write_to_mathematica(out);
	//          => writes matrix A to stream out in mathematica format
	// Version: 2.0
	//

	debug_handler_l(DMESSAGE, "write_to_mathematica(std::ostream &)", DVALUE + 8);

	if (BigintLong) {
		switch(bigint_value.bitfield.get_representation()) {
		case matrix_flags::dense_representation:
			D_base_modul_bigint.write_to_mathematica(bigint_value, out);
			break;
		case matrix_flags::sparse_representation:
			S_base_modul_bigint.write_to_mathematica(bigint_value, out);
			break;
		}
	}
	else {
		switch(long_value.bitfield.get_representation()) {
		case matrix_flags::dense_representation:
			D_base_modul_long.write_to_mathematica(long_value, out);
			break;
		case matrix_flags::sparse_representation:
			S_base_modul_long.write_to_mathematica(long_value, out);
			break;
		}
	}
}



void base_matrix< bigmod >::
read_from_mathematica(std::istream &in)
{
	//
	//    Task: A.read_from_mathematica(dz);
	//          => reads matrix A from stream dz in mathematica format
	// Version: 2.0
	//

	debug_handler_l(DMESSAGE, "read_from_mathematica(std::istream &dz)", DVALUE + 8);

	if (BigintLong) {
		switch(bigint_value.bitfield.get_representation()) {
		case matrix_flags::dense_representation:
			D_base_modul_bigint.read_from_mathematica(bigint_value, in);
			break;
		case matrix_flags::sparse_representation:
			S_base_modul_bigint.read_from_mathematica(bigint_value, in);
			break;
		}
	}
	else {
		switch(long_value.bitfield.get_representation()) {
		case matrix_flags::dense_representation:
			D_base_modul_long.read_from_mathematica(long_value, in);
			break;
		case matrix_flags::sparse_representation:
			S_base_modul_long.read_from_mathematica(long_value, in);
			break;
		}
	}
}



//
// stream handling - MAPLE
//

void base_matrix< bigmod >::
write_to_maple(std::ostream &out) const
{
	//
	//    Task: A.write_to_maple(out);
	//          => writes matrix A to stream out in maple format
	// Version: 2.0
	//

	debug_handler_l(DMESSAGE, "write_to_maple(std::ostream &)", DVALUE + 8);

	if (BigintLong) {
		switch(bigint_value.bitfield.get_representation()) {
		case matrix_flags::dense_representation:
			D_base_modul_bigint.write_to_maple(bigint_value, out);
			break;
		case matrix_flags::sparse_representation:
			S_base_modul_bigint.write_to_maple(bigint_value, out);
			break;
		}
	}
	else {
		switch(long_value.bitfield.get_representation()) {
		case matrix_flags::dense_representation:
			D_base_modul_long.write_to_maple(long_value, out);
			break;
		case matrix_flags::sparse_representation:
			S_base_modul_long.write_to_maple(long_value, out);
			break;
		}
	}
}



void base_matrix< bigmod >::
read_from_maple(std::istream &in)
{
	//
	//    Task: A.read_from_maple(dz);
	//          => reads matrix A from stream dz in maple format
	// Version: 2.0
	//

	debug_handler_l(DMESSAGE, "read_from_maple(std::istream &dz)", DVALUE + 8);

	if (BigintLong) {
		switch(bigint_value.bitfield.get_representation()) {
		case matrix_flags::dense_representation:
			D_base_modul_bigint.read_from_maple(bigint_value, in);
			break;
		case matrix_flags::sparse_representation:
			S_base_modul_bigint.read_from_maple(bigint_value, in);
			break;
		}
	}
	else {
		switch(long_value.bitfield.get_representation()) {
		case matrix_flags::dense_representation:
			D_base_modul_long.read_from_maple(long_value, in);
			break;
		case matrix_flags::sparse_representation:
			S_base_modul_long.read_from_maple(long_value, in);
			break;
		}
	}
}



//
// stream handling - PARI
//

void base_matrix< bigmod >::
write_to_gp(std::ostream &out) const
{
	//
	//    Task: A.write_to_gp(out);
	//          => writes matrix A to stream out in PARI format
	// Version: 2.0
	//

	debug_handler_l(DMESSAGE, "write_to_gp(std::ostream &)", DVALUE + 8);

	if (BigintLong) {
		switch(bigint_value.bitfield.get_representation()) {
		case matrix_flags::dense_representation:
			D_base_modul_bigint.write_to_gp(bigint_value, out);
			break;
		case matrix_flags::sparse_representation:
			S_base_modul_bigint.write_to_gp(bigint_value, out);
			break;
		}
	}
	else {
		switch(long_value.bitfield.get_representation()) {
		case matrix_flags::dense_representation:
			D_base_modul_long.write_to_gp(long_value, out);
			break;
		case matrix_flags::sparse_representation:
			S_base_modul_long.write_to_gp(long_value, out);
			break;
		}
	}
}



void base_matrix< bigmod >::
read_from_gp(std::istream &in)
{
	//
	//    Task: A.read_from_gp(dz);
	//          => reads matrix A from stream dz in PARI format
	// Version: 2.0
	//

	debug_handler_l(DMESSAGE, "read_from_gp(std::istream &)", DVALUE + 8);

	if (BigintLong) {
		switch(bigint_value.bitfield.get_representation()) {
		case matrix_flags::dense_representation:
			D_base_modul_bigint.read_from_gp(bigint_value, in);
			break;
		case matrix_flags::sparse_representation:
			S_base_modul_bigint.read_from_gp(bigint_value, in);
			break;
		}
	}
	else {
		switch(long_value.bitfield.get_representation()) {
		case matrix_flags::dense_representation:
			D_base_modul_long.read_from_gp(long_value, in);
			break;
		case matrix_flags::sparse_representation:
			S_base_modul_long.read_from_gp(long_value, in);
			break;
		}
	}
}



//
// stream handling - KANT
//

void base_matrix< bigmod >::
write_to_kash(std::ostream &out) const
{
	//
	//    Task: A.write_to_kash(out);
	//          => writes matrix A to stream out in kash format
	// Version: 2.0
	//

	debug_handler_l(DMESSAGE, "write_to_kash(std::ostream &)", DVALUE + 8);

	if (BigintLong) {
		switch(bigint_value.bitfield.get_representation()) {
		case matrix_flags::dense_representation:
			D_base_modul_bigint.write_to_kash(bigint_value, out);
			break;
		case matrix_flags::sparse_representation:
			S_base_modul_bigint.write_to_kash(bigint_value, out);
			break;
		}
	}
	else {
		switch(long_value.bitfield.get_representation()) {
		case matrix_flags::dense_representation:
			D_base_modul_long.write_to_kash(long_value, out);
			break;
		case matrix_flags::sparse_representation:
			S_base_modul_long.write_to_kash(long_value, out);
			break;
		}
	}
}



void base_matrix< bigmod >::
read_from_kash(std::istream &in)
{
	//
	//    Task: A.read_from_kash(dz);
	//          => reads matrix A from stream dz in kash format
	// Version: 2.0
	//

	debug_handler_l(DMESSAGE, "read_from_kash(std::istream &dz)", DVALUE + 8);

	if (BigintLong) {
		switch(bigint_value.bitfield.get_representation()) {
		case matrix_flags::dense_representation:
			D_base_modul_bigint.read_from_kash(bigint_value, in);
			break;
		case matrix_flags::sparse_representation:
			S_base_modul_bigint.read_from_kash(bigint_value, in);
			break;
		}
	}
	else {
		switch(long_value.bitfield.get_representation()) {
		case matrix_flags::dense_representation:
			D_base_modul_long.read_from_kash(long_value, in);
			break;
		case matrix_flags::sparse_representation:
			S_base_modul_long.read_from_kash(long_value, in);
			break;
		}
	}
}



//
// stream handling - LaTeX
//

void base_matrix< bigmod >::
write_to_latex(std::ostream &out) const
{
	//
	//    Task: A.write_to_latex(out);
	//          => writes matrix A to stream out in LaTeX format
	// Version: 2.0
	//

	debug_handler_l(DMESSAGE, "write_to_latex(std::ostream &)", DVALUE + 8);

	if (BigintLong) {
		switch(bigint_value.bitfield.get_representation()) {
		case matrix_flags::dense_representation:
			D_base_modul_bigint.write_to_latex(bigint_value, out);
			break;
		case matrix_flags::sparse_representation:
			S_base_modul_bigint.write_to_latex(bigint_value, out);
			break;
		}
	}
	else {
		switch(long_value.bitfield.get_representation()) {
		case matrix_flags::dense_representation:
			D_base_modul_long.write_to_latex(long_value, out);
			break;
		case matrix_flags::sparse_representation:
			S_base_modul_long.write_to_latex(long_value, out);
			break;
		}
	}
}



//
// boolean functions
//

bool base_matrix< bigmod >::
is_column_zero(lidia_size_t c) const
{
	//
	//       Task: A.is_column_zero(c) == true <=> A.value[x][c] == 0,
	//             x=0,...,rows-1
	// Conditions: 0 <= c < A.columns
	//    Version: 2.0
	//

	debug_handler_l(DMESSAGE, "is_column_zero(lidia_size_t)", DVALUE + 10);

	if (BigintLong) {
		if (c< 0 || c >= bigint_value.columns)
			precondition_error_handler(c, "c", "0 <= c < columns",
					    bigint_value.columns , "columns", "",
					    "bool base_matrix< bigmod >::"
					    "is_column_zero(lidia_size_t c) const",
					    DMESSAGE, EMESSAGE[3]);

		return ((bigint_value.bitfield.get_representation() ==
			 matrix_flags::dense_representation)
			? D_base_modul_bigint.is_column_zero(bigint_value, c)
			: S_base_modul_bigint.is_column_zero(bigint_value, c));
	}
	else {
		if (c< 0 || c >= long_value.columns)
			precondition_error_handler(c, "c", "0 <= c < columns",
					    long_value.columns , "columns", "",
					    "bool base_matrix< bigmod >::"
					    "is_column_zero(lidia_size_t c) const",
					    DMESSAGE, EMESSAGE[3]);

		return ((long_value.bitfield.get_representation() ==
			 matrix_flags::dense_representation)
			? D_base_modul_long.is_column_zero(long_value, c)
			: S_base_modul_long.is_column_zero(long_value, c));
	}
}



bool base_matrix< bigmod >::
is_row_zero(lidia_size_t r) const
{
	//
	//       Task: A.is_row_zero(r) == true <=> A.value[r][x] == 0,
	//             x=0,...,columns-1
	// Conditions: 0 <= r < A.rows
	//    Version: 2.0
	//

	debug_handler_l(DMESSAGE, "is_row_zero(lidia_size_t)", DVALUE + 10);

	if (BigintLong) {
		if (r< 0 || r >= bigint_value.rows)
			precondition_error_handler(r, "r", "0 <= r < rows",
					    bigint_value.rows, "rows", "",
					    "bool base_matrix< bigmod >::"
					    "is_row_zero(lidia_size_t r) const",
					    DMESSAGE, EMESSAGE[3]);

		return (bigint_value.bitfield.get_representation() ==
			matrix_flags::dense_representation)
			? D_base_modul_bigint.is_row_zero(bigint_value, r)
			: S_base_modul_bigint.is_row_zero(bigint_value, r);
	}
	else {
		if (r< 0 || r >= long_value.rows)
			precondition_error_handler(r, "r", "0 <= r < rows",
					    long_value.rows, "rows", "",
					    "bool base_matrix< bigmod >::"
					    "is_row_zero(lidia_size_t r) const",
					    DMESSAGE, EMESSAGE[3]);

		return (long_value.bitfield.get_representation() ==
			matrix_flags::dense_representation)
			? D_base_modul_long.is_row_zero(long_value, r)
			: S_base_modul_long.is_row_zero(long_value, r);
	}
}



bool base_matrix< bigmod >::
is_matrix_zero() const
{
	//
	//    Task: A.is_matrix_zero() == true <=> A.value[x][y] == 0,
	//          x=0,...,rows-1 and y=0,...,columns-1
	// Version: 2.0
	//

	debug_handler_l(DMESSAGE, "is_matrix_zero()", DVALUE + 10);

	if (BigintLong) {
		return (bigint_value.bitfield.get_representation() ==
			matrix_flags::dense_representation)
			? D_base_modul_bigint.is_matrix_zero(bigint_value)
			: S_base_modul_bigint.is_matrix_zero(bigint_value);
	}
	else {
		return (long_value.bitfield.get_representation() ==
			matrix_flags::dense_representation)
			? D_base_modul_long.is_matrix_zero(long_value)
			: S_base_modul_long.is_matrix_zero(long_value);
	}
}



//
// internal functions
//

void base_matrix< bigmod >::
change_representation(unsigned long art)
{
	debug_handler_l(DMESSAGE, "change_representation(unsigned long)",
			DVALUE + 9);

	lidia_size_t i, j, k;
	if (BigintLong) {
		switch (bigint_value.bitfield.get_representation()) {
		case matrix_flags::dense_representation:
		{
			switch (art) {
			case matrix_flags::sparse_representation:
			{
				bigint_value.index = new lidia_size_t *[bigint_value.rows];
				memory_handler(bigint_value.index, DMESSAGE,
					       "change_representation(unsigned long) :: "
					       "Error in memory allocation (index)");

				bigint_value.value_counter
					= new lidia_size_t[bigint_value.rows];
				memory_handler(value_counter, DMESSAGE,
					       "change_representation(unsigned long) :: "
					       "Error in memory allocation (value_counter)");

				bigint_value.allocated = new lidia_size_t[bigint_value.rows];
				memory_handler(bigint_value.allocated, DMESSAGE,
					       "change_representation(unsigned long) :: "
					       "Error in memory allocation (allocated)");

				for (i = 0; i < bigint_value.rows; i++) {
					// set allocated
					bigint_value.allocated[i] = bigint_value.columns;

					// set value_counter
					register lidia_size_t size = bigint_value.columns;
					for (j = 0; j < bigint_value.columns; j++)
						if (bigint_value.value[i][j] == bigint_value.Zero)
							size--;

					bigint_value.value_counter[i] = size;

					// set index
					bigint_value.index[i]
						= new lidia_size_t[bigint_value.columns];
					memory_handler(bigint_value.index[i], DMESSAGE,
						       "change_representation(unsigned long) :: "
						       "Error in memory allocation (index[i])");

					register lidia_size_t p = 0;
					for (k = 0; k < bigint_value.columns; k++)
						if (bigint_value.value[i][k] != bigint_value.Zero) {
							LiDIA::swap(bigint_value.value[i][p],
								       bigint_value.value[i][k]);
							bigint_value.index[i][p] = k;
							p++;
						}
				}
				break;
			}
			default:
				break;
			}
			break;
		}
		case matrix_flags::sparse_representation:
		{
			switch (art) {
			case matrix_flags::dense_representation:
			{
				bigint *tmp = NULL;
				for (i = 0; i < bigint_value.rows; i++) {
					tmp = new bigint[bigint_value.columns];
					memory_handler(tmp, DMESSAGE,
						       "change_representation(unsigned long) :: "
						       "Error in memory allocation (tmp)");

					k = 0;
					for (j = 0; j < bigint_value.columns; j++) {
						if (k < bigint_value.value_counter[i]
						    && bigint_value.index[i][k] == j) {
							tmp[j] = bigint_value.value[i][k];
							k++;
						}
						else {
							tmp[j] = bigint_value.Zero;
						}
					}
					delete[] bigint_value.value[i];
					bigint_value.value[i] = tmp;
					delete[] bigint_value.index[i];
				}
				delete[] bigint_value.value_counter;
				bigint_value.value_counter = NULL;

				delete[] bigint_value.allocated;
				bigint_value.allocated = NULL;

				delete[] bigint_value.index;
				bigint_value.index = NULL;

				break;
			}
			default:
				break;
			}
		}
		break;
		}
	}
	else {
		switch (long_value.bitfield.get_representation()) {
		case matrix_flags::dense_representation:
		{
			switch (art) {
			case matrix_flags::sparse_representation:
			{
				long_value.index = new lidia_size_t *[long_value.rows];
				memory_handler(long_value.index, DMESSAGE,
					       "change_representation(unsigned long) :: "
					       "Error in memory allocation (index)");

				long_value.value_counter
					= new lidia_size_t[long_value.rows];
				memory_handler(value_counter, DMESSAGE,
					       "change_representation(unsigned long) :: "
					       "Error in memory allocation (value_counter)");

				long_value.allocated = new lidia_size_t[long_value.rows];
				memory_handler(long_value.allocated, DMESSAGE,
					       "change_representation(unsigned long) :: "
					       "Error in memory allocation (allocated)");

				for (i = 0; i < long_value.rows; i++) {
					// set allocated
					long_value.allocated[i] = long_value.columns;

					// set value_counter
					register lidia_size_t size = long_value.columns;
					for (j = 0; j < long_value.columns; j++)
						if (long_value.value[i][j] == long_value.Zero)
							size--;

					long_value.value_counter[i] = size;

					// set index
					long_value.index[i]
						= new lidia_size_t[long_value.columns];
					memory_handler(long_value.index[i], DMESSAGE,
						       "change_representation(unsigned long) :: "
						       "Error in memory allocation (index[i])");

					register lidia_size_t p = 0;
					for (k = 0; k < long_value.columns; k++)
						if (long_value.value[i][k] != long_value.Zero) {
							LiDIA::swap(long_value.value[i][p],
								       long_value.value[i][k]);
							long_value.index[i][p] = k;
							p++;
						}
				}
				break;
			}
			default:
				break;
			}
			break;
		}
		case matrix_flags::sparse_representation:
		{
			switch (art) {
			case matrix_flags::dense_representation:
			{
				long *tmp = NULL;
				for (i = 0; i < long_value.rows; i++) {
					tmp = new long[long_value.columns];
					memory_handler(tmp, DMESSAGE,
						       "change_representation(unsigned long) :: "
						       "Error in memory allocation (tmp)");

					k = 0;
					for (j = 0; j < long_value.columns; j++) {
						if (k < long_value.value_counter[i]
						    && long_value.index[i][k] == j) {
							tmp[j] = long_value.value[i][k];
							k++;
						}
						else {
							tmp[j] = long_value.Zero;
						}
					}
					delete[] long_value.value[i];
					long_value.value[i] = tmp;
					delete[] long_value.index[i];
				}
				delete[] long_value.value_counter;
				long_value.value_counter = NULL;

				delete[] long_value.allocated;
				long_value.allocated = NULL;

				delete[] long_value.index;
				long_value.index = NULL;

				break;
			}
			default:
				break;
			}
		}
		break;
		}
	}
}



void base_matrix< bigmod >::
status_report()
{
	debug_handler_l(DMESSAGE, "status_report()", DVALUE + 9);

	lidia_size_t i;
	if (BigintLong) {
		std::cout << "configuration (bigint): "
			  << "\n\t rows = " << bigint_value.rows
			  << "\n\t columns = " << bigint_value.columns
			  << "\n\t sparse_rows = " << bigint_value.sparse_rows
			  << "\n\t sparse_columns = " << bigint_value.sparse_columns << std::endl;

		if (bigint_value.bitfield.get_representation() ==
		    matrix_flags::sparse_representation)
			for (i = 0; i < bigint_value.rows; i++) {
				std::cout << "row " << i << ": allocated = "
					  << bigint_value.allocated[i]
					  << " value_counter = " << bigint_value.value_counter[i]
					  << std::endl;
			}
	}
	else {
		std::cout << "configuration (long): "
			  << "\n\t rows = " << bigint_value.rows
			  << "\n\t columns = " << bigint_value.columns
			  << "\n\t sparse_rows = " << bigint_value.sparse_rows
			  << "\n\t sparse_columns = " << bigint_value.sparse_columns << std::endl;

		if (bigint_value.bitfield.get_representation() ==
		    matrix_flags::sparse_representation)
			for (i = 0; i < bigint_value.rows; i++) {
				std::cout << "row " << i << ": allocated = "
					  << bigint_value.allocated[i]
					  << " value_counter = " << bigint_value.value_counter[i]
					  << std::endl;
			}
	}
}



template std::ostream & operator << (std::ostream &, const base_matrix< bigmod > &);
template std::istream & operator >> (std::istream &, base_matrix< bigmod > &);
template void swap(base_matrix< bigmod > &, base_matrix< bigmod > &);
template void assign(base_matrix< bigmod > &, const base_matrix< bigmod > &);
template void diag(base_matrix< bigmod > &, const bigmod &, const bigmod &);
template void trans(base_matrix< bigmod > &, const base_matrix< bigmod > &);
template base_matrix< bigmod > trans(const base_matrix< bigmod > &);



#undef DVALUE
#undef DMESSAGE
#undef EMESSAGE



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
