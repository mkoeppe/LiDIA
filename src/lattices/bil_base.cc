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
//	Author	: Werner Backes (WB), Thorsten Lauer (TL)
//	Changes	: See CVS log
//
//==============================================================================================


#ifdef HAVE_CONFIG_H
# include	"config.h"
#endif
#include	"LiDIA/bigint_lattice.h"



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



//
// Set special lattice characteristics
//
void bigint_lattice::set_red_orientation_columns()
{
	debug_handler("bigint_lattice", "set_red_orentation_columns()");
	bitfield.lattice_mode &= !REDUCE_ROWS;
}



void bigint_lattice::set_red_orientation_rows()
{
	debug_handler("bigint_lattice", "set_red_orentation_rows()");
	bitfield.lattice_mode |= REDUCE_ROWS;
}



void bigint_lattice::set_gram_flag()
{
	debug_handler("bigint_lattice", "set_gram_flag()");
	bitfield.info_mode |= GRAM_MATRIX;
}



void bigint_lattice::set_basis_flag()
{
	debug_handler("bigint_lattice", "set_basis_flag()");
	if (bitfield.lattice_mode & REDUCE_ROWS)
		bitfield.structure_mode |= ROWS_LININD;
	else
		bitfield.structure_mode |= COLUMNS_LININD;
}



void bigint_lattice::delete_gram_flag()
{
	debug_handler("bigint_lattice", "delete_gram_flag()");
	bitfield.info_mode &= !GRAM_MATRIX;
}



void bigint_lattice::delete_basis_flag()
{
	debug_handler("bigint_lattice", "delete_basis_flag()");
	if (bitfield.lattice_mode & REDUCE_ROWS)
		bitfield.structure_mode &= !ROWS_LININD;
	else
		bitfield.structure_mode &= !COLUMNS_LININD;
}



//
// get the special lattice characteristics
//
bool bigint_lattice::get_red_orientation()
{
	debug_handler("bigint_lattice", "get_red_orentation()");
	return ((bitfield.lattice_mode & REDUCE_ROWS) ? false : true);
	//  return (!((bool )(bitfield.lattice_mode & REDUCE_ROWS)));
}



bool bigint_lattice::get_basis_flag()
{
	debug_handler("bigint_lattice", "get_basis_flag()");
	if (bitfield.lattice_mode & REDUCE_ROWS)
		return ((bitfield.structure_mode & ROWS_LININD) ? true : false);
	else
		return ((bitfield.structure_mode & COLUMNS_LININD) ? true : false);
}



bool bigint_lattice::get_gram_flag()
{
	debug_handler("bigint_lattice", "get_gram_flag()");
	return ((bitfield.structure_mode & GRAM_MATRIX) ? true : false);
}



bool bigint_lattice::chk_basis()
{
	debug_handler("bigint_lattice", "chk_basis()");
	if (bitfield.lattice_mode & REDUCE_ROWS)
		return ((bitfield.structure_mode & ROWS_LININD) ? true : false);
	else
		return ((bitfield.structure_mode & COLUMNS_LININD) ? true : false);
}



bool bigint_lattice::chk_gram()
{
	debug_handler("bigint_lattice", "chk_gram()");
	return ((bitfield.info_mode & GRAM_MATRIX) ? true : false);
	//  return ((bool )(bitfield.info_mode & GRAM_MATRIX));
}



bool bigint_lattice::chk_trans()
{
	debug_handler("bigint_lattice", "chk_trans()");
	if (chk_gram())
		return (false);
	else {
		if (((bitfield.lattice_mode & REDUCE_ROWS) && (bitfield.structure_mode & COLUMN_ORIENTED)) ||
		    (!(bitfield.lattice_mode & REDUCE_ROWS) && !(bitfield.structure_mode & COLUMN_ORIENTED)))
			return (true);
		else
			return (false);
	}
	//    return (((bool )(bitfield.lattice_mode & REDUCE_ROWS)) ==
	//	      ((bool )(bitfield.structure_mode & COLUMN_ORIENTED)));

}



bool bigint_lattice::chk_reduce_columns()
{
	debug_handler("bigint_lattice", "chk_reduce_columns()");
	return ((bitfield.lattice_mode & REDUCE_ROWS) ? false : true);
	//  return(!((bool )(bitfield.lattice_mode & REDUCE_ROWS)));
}



//
// Dimension checking
//
bool bigint_lattice::chk_mlll_dim()
{
	debug_handler("bigint_lattice", "chk_mlll_dim()");
	return (rows+static_cast<lidia_size_t>(chk_reduce_columns()) ==
		columns+static_cast<lidia_size_t>(!chk_reduce_columns()));
}



bool bigint_lattice::chk_lll_dim()
{
	debug_handler("bigint_lattice", "chk_lll_dim()");
	bool chk = false;
	if (chk_gram())
		chk = (rows == columns);
	else {
		if (chk_reduce_columns())
			chk = (rows >= columns);
		else
			chk = (columns >= rows);
		if (!chk_basis())
			chk = true;
	}
	return (chk);
}



//
// Other checkings
//
void bigint_lattice::chk_corr_param(double& y)
{
	debug_handler("bigint_lattice", "chk_corr_param(y)");
	if ((y > 1.0) || (y <= 0.5))
		y = 0.99;
}



void bigint_lattice::chk_corr_param(sdigit& nom, sdigit& denom)
{
	debug_handler("bigint_lattice", "chk_corr_param(nom, denom)");
	double y;
	y = static_cast<double>(nom)/static_cast<double>(denom);
	if ((y > 1.0) || (y <= 0.5)) {
		nom = 99;
		denom = 100;
	}
}



//
// Schnorr - Euchner
//
void bigint_lattice::lll_schnorr_euchner_orig(double y, lattice_info& li,
					      sdigit x_factor)
{
	debug_handler("bigint_lattice", "lll_schnorr_euchner_orig(y, li, factor)");
	ALG_DEF_BASIS(bigint, vector_op, Normal) alg_basis;
	ALG_DEF_GENSYS(bigint, vector_op, Normal) alg_gensys;

	if (!(chk_lll_dim()))
		lidia_error_handler("bigint_lattice", "lll_schnorr_euchner_orig(y, li, factor) ::"
				    "wrong dimension (no basis)");
	dense_alg< bigint > da;
	da.b.y = y;
	chk_corr_param(da.b.y);
	da.b.alg_trans = false;
	x_factor = ((x_factor < 1) ? 1 : x_factor);
	da.d.bit_prec = x_factor*DOUBLE_MANTISSA_BITS;
	Tr_dense_create(da);
	if (chk_basis()) {
		if (chk_gram())
			ALG_CALL(alg_basis, lll_gram, da, li, x_factor)
				else
					ALG_CALL(alg_basis, lll, da, li, x_factor)
						}
	else {
		if (chk_gram())
			ALG_CALL(alg_gensys, lll_gram, da, li, x_factor)
				else
					ALG_CALL(alg_gensys, lll, da, li, x_factor)
						}
	Tr_dense_extract(da);
}



void bigint_lattice::lll_schnorr_euchner_orig(ring_matrix< bigint > & T, double y,
			                      lattice_info& li, sdigit x_factor)
{
	debug_handler("bigint_lattice", "lll_schnorr_euchner_orig(T, y, li, factor)");
	ALG_DEF_BASIS(bigint, vector_op, Normal) alg_basis;
	ALG_DEF_GENSYS(bigint, vector_op, Normal) alg_gensys;

	if (!(chk_lll_dim()))
		lidia_error_handler("bigint_lattice", "lll_schnorr_euchner_orig(T, y, li, "
				    " factor) :: wrong dimension (no basis)");
	dense_alg< bigint > da;
	da.b.y = y;
	chk_corr_param(da.b.y);
	da.b.alg_trans = true;
	da.s.TMatrix = &T;
	x_factor = ((x_factor < 1) ? 1 : x_factor);
	da.d.bit_prec = x_factor*DOUBLE_MANTISSA_BITS;
	Tr_dense_create(da);
	if (chk_basis()) {
		if (chk_gram())
			ALG_CALL(alg_basis, lll_gram, da, li, x_factor)
				else
					ALG_CALL(alg_basis, lll, da, li, x_factor)
						}
	else {
		if (chk_gram())
			ALG_CALL(alg_gensys, lll_gram, da, li, x_factor)
				else
					ALG_CALL(alg_gensys, lll, da, li, x_factor)
						}
	Tr_dense_extract(da);
}



void bigint_lattice::lll_schnorr_euchner_fact(double y, lattice_info& li,
					      sdigit x_factor)
{
	debug_handler("bigint_lattice", "lll_schnorr_euchner_fact(y, li, factor)");
	ALG_DEF_BASIS(bigint, vector_op, VariationI) alg_basis;
	ALG_DEF_GENSYS(bigint, vector_op, VariationI) alg_gensys;

	if (!(chk_lll_dim()))
		lidia_error_handler("bigint_lattice", "lll_schnorr_euchner_fact(y, li, factor) ::"
				    "wrong dimension (no basis)");
	dense_alg< bigint > da;
	da.b.y = y;
	chk_corr_param(da.b.y);
	da.b.alg_trans = false;
	x_factor = ((x_factor < 1) ? 1 : x_factor);
	da.d.bit_prec = x_factor*DOUBLE_MANTISSA_BITS;
	Tr_dense_create(da);
	da.d.bit_factor = TrD_search_factor(da);
	if (chk_basis()) {
		if (chk_gram())
			ALG_CALL(alg_basis, lll_gram, da, li, x_factor)
				else
					ALG_CALL(alg_basis, lll, da, li, x_factor)
						}
	else {
		if (chk_gram())
			ALG_CALL(alg_gensys, lll_gram, da, li, x_factor)
				else
					ALG_CALL(alg_gensys, lll, da, li, x_factor)
						}
	Tr_dense_extract(da);
}



void bigint_lattice::lll_schnorr_euchner_fact(ring_matrix< bigint > & T, double y,
			                      lattice_info& li, sdigit x_factor)
{
	debug_handler("bigint_lattice", "lll_schnorr_euchner_fact(T, y, li, factor)");
	ALG_DEF_BASIS(bigint, vector_op, VariationI) alg_basis;
	ALG_DEF_GENSYS(bigint, vector_op, VariationI) alg_gensys;

	if (!(chk_lll_dim()))
		lidia_error_handler("bigint_lattice", "lll_schnorr_euchner_fact(T, y, li, "
				    " factor) :: wrong dimension (no basis)");
	dense_alg< bigint > da;
	da.b.y = y;
	chk_corr_param(da.b.y);
	da.b.alg_trans = true;
	da.s.TMatrix = &T;
	x_factor = ((x_factor < 1) ? 1 : x_factor);
	da.d.bit_prec = x_factor*DOUBLE_MANTISSA_BITS;
	Tr_dense_create(da);
	da.d.bit_factor = TrD_search_factor(da);
	if (chk_basis()) {
		if (chk_gram())
			ALG_CALL(alg_basis, lll_gram, da, li, x_factor)
				else
					ALG_CALL(alg_basis, lll, da, li, x_factor)
						}
	else {
		if (chk_gram())
			ALG_CALL(alg_gensys, lll_gram, da, li, x_factor)
				else
					ALG_CALL(alg_gensys, lll, da, li, x_factor)
						}
	Tr_dense_extract(da);
}



//
// Schnorr - Euchner for user defined Scalar Product
//
void bigint_lattice::lll_schnorr_euchner_orig(double y, lattice_info& li,
			                      user_SP SP, sdigit x_factor)
{
	debug_handler("bigint_lattice", "lll_schnorr_euchner_orig(y, li, SP, factor)");
	ALG_DEF_BASIS(bigint, vector_op_SP, Normal) alg_basis;
	ALG_DEF_GENSYS(bigint, vector_op_SP, Normal) alg_gensys;

	ALG_POINTER(alg_basis, SP, bin);
	ALG_POINTER(alg_gensys, SP, bin);

	if (!(chk_lll_dim()))
		lidia_error_handler("bigint_lattice", "lll_schnorr_euchner_orig(y, li, SP, "
				    " factor) :: wrong dimension (no basis)");
	dense_alg< bigint > da;
	da.b.y = y;
	chk_corr_param(da.b.y);
	da.b.alg_trans = false;
	x_factor = ((x_factor < 1) ? 1 : x_factor);
	da.d.bit_prec = x_factor*DOUBLE_MANTISSA_BITS;
	Tr_dense_create(da);
	if (chk_basis()) {
		if (chk_gram())
			ALG_CALL(alg_basis, lll_gram, da, li, x_factor)
				else
					ALG_CALL(alg_basis, lll, da, li, x_factor)
						}
	else {
		if (chk_gram())
			ALG_CALL(alg_gensys, lll_gram, da, li, x_factor)
				else
					ALG_CALL(alg_gensys, lll, da, li, x_factor)
						}
	Tr_dense_extract(da);
}



void bigint_lattice::lll_schnorr_euchner_orig(ring_matrix< bigint > & T, double y,
			                      lattice_info& li, user_SP SP,
					      sdigit x_factor)
{
	debug_handler("bigint_lattice", "lll_schnorr_euchner_orig(T, y, li, SP, factor)");
	ALG_DEF_BASIS(bigint, vector_op_SP, Normal) alg_basis;
	ALG_DEF_GENSYS(bigint, vector_op_SP, Normal) alg_gensys;

	ALG_POINTER(alg_basis, SP, bin);
	ALG_POINTER(alg_gensys, SP, bin);

	if (!(chk_lll_dim()))
		lidia_error_handler("bigint_lattice", "lll_schnorr_euchner_orig(T, y, li, SP, "
				    " factor) :: wrong dimension (no basis)");
	dense_alg< bigint > da;
	da.b.y = y;
	chk_corr_param(da.b.y);
	da.b.alg_trans = true;
	da.s.TMatrix = &T;
	x_factor = ((x_factor < 1) ? 1 : x_factor);
	da.d.bit_prec = x_factor*DOUBLE_MANTISSA_BITS;
	Tr_dense_create(da);
	if (chk_basis()) {
		if (chk_gram())
			ALG_CALL(alg_basis, lll_gram, da, li, x_factor)
				else
					ALG_CALL(alg_basis, lll, da, li, x_factor)
						}
	else {
		if (chk_gram())
			ALG_CALL(alg_gensys, lll_gram, da, li, x_factor)
				else
					ALG_CALL(alg_gensys, lll, da, li, x_factor)
						}
	Tr_dense_extract(da);
}



void bigint_lattice::lll_schnorr_euchner_fact(double y, lattice_info& li,
			                      user_SP SP, sdigit x_factor)
{
	debug_handler("bigint_lattice", "lll_schnorr_euchner_fact(y, li, SP, factor)");
	ALG_DEF_BASIS(bigint, vector_op_SP, VariationI) alg_basis;
	ALG_DEF_GENSYS(bigint, vector_op_SP, VariationI) alg_gensys;

	ALG_POINTER(alg_basis, SP, bin);
	ALG_POINTER(alg_gensys, SP, bin);

	if (!(chk_lll_dim()))
		lidia_error_handler("bigint_lattice", "lll_schnorr_euchner_fact(y, li, SP, "
				    "factor) :: wrong dimension (no basis)");
	dense_alg< bigint > da;
	da.b.y = y;
	chk_corr_param(da.b.y);
	da.b.alg_trans = false;
	x_factor = ((x_factor < 1) ? 1 : x_factor);
	da.d.bit_prec = x_factor*DOUBLE_MANTISSA_BITS;
	Tr_dense_create(da);
	da.d.bit_factor = TrD_search_factor(da);
	if (chk_basis()) {
		if (chk_gram())
			ALG_CALL(alg_basis, lll_gram, da, li, x_factor)
				else
					ALG_CALL(alg_basis, lll, da, li, x_factor)
						}
	else {
		if (chk_gram())
			ALG_CALL(alg_gensys, lll_gram, da, li, x_factor)
				else
					ALG_CALL(alg_gensys, lll, da, li, x_factor)
						}
	Tr_dense_extract(da);
}



void bigint_lattice::lll_schnorr_euchner_fact(ring_matrix< bigint > & T, double y,
			                      lattice_info& li, user_SP SP,
					      sdigit x_factor)
{
	debug_handler("bigint_lattice", "lll_schnorr_euchner_fact(T, y, li, SP, factor)");
	ALG_DEF_BASIS(bigint, vector_op_SP, VariationI) alg_basis;
	ALG_DEF_GENSYS(bigint, vector_op_SP, VariationI) alg_gensys;

	ALG_POINTER(alg_basis, SP, bin);
	ALG_POINTER(alg_gensys, SP, bin);

	if (!(chk_lll_dim()))
		lidia_error_handler("bigint_lattice", "lll_schnorr_euchner_fact(T, y, li, SP, "
				    "factor) :: wrong dimension (no basis)");
	dense_alg< bigint > da;
	da.b.y = y;
	chk_corr_param(da.b.y);
	da.b.alg_trans = true;
	da.s.TMatrix = &T;
	x_factor = ((x_factor < 1) ? 1 : x_factor);
	da.d.bit_prec = x_factor*DOUBLE_MANTISSA_BITS;
	Tr_dense_create(da);
	da.d.bit_factor = TrD_search_factor(da);
	if (chk_basis()) {
		if (chk_gram())
			ALG_CALL(alg_basis, lll_gram, da, li, x_factor)
				else
					ALG_CALL(alg_basis, lll, da, li, x_factor)
						}
	else {
		if (chk_gram())
			ALG_CALL(alg_gensys, lll_gram, da, li, x_factor)
				else
					ALG_CALL(alg_gensys, lll, da, li, x_factor)
						}
	Tr_dense_extract(da);
}



//
// Benne de Weger
//
void bigint_lattice::lll_benne_de_weger(sdigit y_nom, sdigit y_denom,
					lattice_info& li)
{
	debug_handler("bigint_lattice", "lll_benne_de_weger(y_nom, y_denom, li)");

	if (!chk_basis() || !chk_lll_dim() || chk_gram())
		lidia_error_handler("bigint_lattice", "lll_benne_de_weger(y_nom, y_denom, "
				    "li) :: not implemented for gensys or gram");

	dense_alg< bigint > da;
	da.b.y_nom = y_nom;
	da.b.y_denom = y_denom;
	chk_corr_param(da.b.y_nom, da.b.y_denom);
	da.b.alg_trans = false;
	Tr_dense_create(da);
	TrD_lll(da, li);
	Tr_dense_extract(da);
}



void bigint_lattice::lll_benne_de_weger(ring_matrix< bigint > & T,
					sdigit y_nom, sdigit y_denom,
					lattice_info& li)
{
	debug_handler("bigint_lattice", "lll_benne_de_weger(T, y_nom, y_denom, li)");
	if (!chk_basis() || !chk_lll_dim() || chk_gram())
		lidia_error_handler("bigint_lattice", "lll_benne_de_weger(T, y_nom, y_denom, "
				    " li) :: not implemented for gensys or gram");

	dense_alg< bigint > da;
	da.b.y_nom = y_nom;
	da.b.y_denom = y_denom;
	chk_corr_param(da.b.y_nom, da.b.y_denom);
	da.b.alg_trans = false;
	da.s.TMatrix = &T;
	Tr_dense_create(da);
	TrD_lll_trans(da, li);
	Tr_dense_extract(da);
	if (da.b.transpose)
		LiDIA::multiply(*this, *this, T);
	else
		LiDIA::multiply(*this, T, *this);
}



//
// Buchmann - Kessler
//
void bigint_lattice::buchmann_kessler(ring_matrix< bigint > & T, double y,
				      lattice_info& li)
{
	debug_handler("bigint_lattice", "buchmann_kessler(T, y, li)");
	if (chk_gram())
		lidia_error_handler("bigint_lattice", "buchmann_kessler(T, y, "
				    " li) :: not avaidable for gram");

	dense_alg< bigint > da;
	da.b.y = y;
	chk_corr_param(da.b.y);
	da.b.alg_trans = true;
	da.s.TMatrix = &T;
	Tr_dense_create(da);
	TrD_buchmann_kessler(da, li);
	Tr_dense_extract(da);
	if (da.b.transpose)
		LiDIA::multiply(*this, *this, T);
	else
		LiDIA::multiply(*this, T, *this);
}



//
// Modified lll
//
void bigint_lattice::mlll(double y, bigint*& v, lattice_info& li)
{
	debug_handler("bigint_lattice", "mlll(y, v, li)");

	if (!chk_mlll_dim())
		lidia_error_handler("bigint_lattice", "mlll(y, v, li) ::"
				    " wrong dimension");

	dense_alg< bigint > da;
	da.b.y = y;
	chk_corr_param(da.b.y);
	da.b.alg_trans = false;
	Tr_dense_create(da);
	v = TrD_mlll_bfl(da, li);
	Tr_dense_extract(da);
}



void bigint_lattice::mlll(double y, base_vector< bigint > & bv, lattice_info& li)
{
	debug_handler("bigint_lattice", "mlll(y, bv, li)");
	if (!chk_mlll_dim())
		lidia_error_handler("bigint_lattice", "mlll(y, bv, li) ::"
				    " wrong dimension");

	dense_alg< bigint > da;
	bigint *v;
	da.b.y = y;
	chk_corr_param(da.b.y);
	da.b.alg_trans = false;
	Tr_dense_create(da);
	v = TrD_mlll_bfl(da, li);
	base_vector< bigint > temp(v, da.b.rows);
	bv = temp;
	Tr_dense_extract(da);
}



void bigint_lattice::close_vector(const base_vector< bigint > & v,
				  base_vector< bigint > & cv,
				  sdigit x_factor)
{
	debug_handler("bigint_lattice", "close_vector(v, cv)");
	ALG_DEF_BASIS(bigint, vector_op, Normal) alg_basis;
	bigint_lattice cA(*this);
	dense_alg< bigint > da;
	lattice_info li;
	p_vector< bigint > vector;
	bigfloat sqrtBmax;
	bigint C, B, Bmax;

	if ((!chk_basis()) || (chk_gram()))
		lidia_error_handler("bigint_lattice", "close_vector(c, cv) :: "
				    "not implemented for gram or gensys !");
	cA.set_no_of_rows(get_no_of_rows()+1);
	cA.set_no_of_columns(get_no_of_columns()+1);
	da.b.alg_trans = false;
	x_factor = ((x_factor < 1) ? 1 : x_factor);
	da.d.bit_prec = x_factor*DOUBLE_MANTISSA_BITS;
	da.b.y = 0.99;
	cA.Tr_dense_create(da);
	vector.vectsize = da.b.columns;
	//
	// Berechne C
	//
	Bmax.assign_zero();
	for (lidia_size_t i = 0; i < da.b.rows; i++) {
		vector.scalprod(B, da.s.value[i], da.s.value[i]);
		if (B > Bmax)
			Bmax.assign(B);
	}
	sqrtBmax.assign(Bmax);
	sqrt(sqrtBmax, sqrtBmax);
	sqrtBmax.divide_by_2();
	LiDIA::multiply(sqrtBmax, sqrtBmax, sqrt(bigfloat(3.0)));
	ceil(C, sqrtBmax);
	//
	// Insert vector
	//
	if (v.size() != da.b.columns-1)
		lidia_error_handler("bigint_lattice", "close_vector(c, cv) :: "
				    "illegal size of vector !");
	for (lidia_size_t i = 0; i < da.b.columns-1; i++)
		da.s.value[da.b.rows-1][i].assign(v[i]);
	da.s.value[da.b.rows-1][da.b.columns-1].assign(C);
	ALG_CALL(alg_basis, lll, da, li, x_factor)
		for (lidia_size_t i = 0; i < da.b.columns-1; i++)
			LiDIA::subtract(da.s.value[da.b.rows-1][i], v[i], da.s.value[da.b.rows-1][i]);
	cv.set_size(da.b.columns-1);
	cv.set_data(da.s.value[da.b.rows-1], da.b.columns-1);
}



//
// Interface to Algorithms
//
// friend functions
//
bigint_lattice lll_schnorr_euchner_orig(const bigint_lattice& A, double y,
		                        lattice_info& li, sdigit factor)
{
	bigint_lattice tA(A);
	tA.lll_schnorr_euchner_orig(y, li, factor);
	return(tA);
}



bigint_lattice lll_schnorr_euchner_orig(const bigint_lattice& A,
		                        ring_matrix< bigint > & T, double y,
		                        lattice_info& li, sdigit factor)
{
	bigint_lattice tA(A);
	tA.lll_schnorr_euchner_orig(T, y, li, factor);
	return(tA);
}



bigint_lattice lll_schnorr_euchner_orig(const bigint_lattice& A, double y,
		                        sdigit factor)
{
	bigint_lattice tA(A);
	tA.lll_schnorr_euchner_orig(y, factor);
	return(tA);
}



bigint_lattice lll_schnorr_euchner_orig(const bigint_lattice& A,
	                                ring_matrix< bigint > & T, double y,
		                        sdigit factor)
{
	bigint_lattice tA(A);
	tA.lll_schnorr_euchner_orig(T, y, factor);
	return(tA);
}



//
// user defined Scalar Product
//
bigint_lattice lll_schnorr_euchner_orig(const bigint_lattice& A, double y,
		                        lattice_info& li, user_SP SP,
					sdigit factor)
{
	bigint_lattice tA(A);
	tA.lll_schnorr_euchner_orig(y, li, SP, factor);
	return(tA);
}



bigint_lattice lll_schnorr_euchner_orig(const bigint_lattice& A,
					ring_matrix< bigint > & T, double y,
					lattice_info& li, user_SP SP, sdigit factor)
{
	bigint_lattice tA(A);
	tA.lll_schnorr_euchner_orig(T, y, li, SP, factor);
	return(tA);
}



bigint_lattice lll_schnorr_euchner_orig(const bigint_lattice& A, double y,
		                        user_SP SP, sdigit factor)
{
	bigint_lattice tA(A);
	tA.lll_schnorr_euchner_orig(y, SP, factor);
	return(tA);
}



bigint_lattice lll_schnorr_euchner_orig(const bigint_lattice& A,
	                                ring_matrix< bigint > & T, double y,
		                        user_SP SP, sdigit factor)
{
	bigint_lattice tA(A);
	tA.lll_schnorr_euchner_orig(T, y, SP, factor);
	return(tA);
}



//
// Variation of Schnorr - Euchner
//
// friend functions
//
bigint_lattice lll_schnorr_euchner_fact(const bigint_lattice& A, double y,
		                        lattice_info& li, sdigit factor)
{
	bigint_lattice tA(A);
	tA.lll_schnorr_euchner_fact(y, li, factor);
	return(tA);
}



bigint_lattice lll_schnorr_euchner_fact(const bigint_lattice& A,
		                        ring_matrix< bigint > & T, double y,
		                        lattice_info& li, sdigit factor)
{
	bigint_lattice tA(A);
	tA.lll_schnorr_euchner_fact(T, y, li, factor);
	return(tA);
}



bigint_lattice lll_schnorr_euchner_fact(const bigint_lattice& A, double y,
		                        sdigit factor)
{
	bigint_lattice tA(A);
	tA.lll_schnorr_euchner_fact(y, factor);
	return(tA);
}



bigint_lattice lll_schnorr_euchner_fact(const bigint_lattice& A,
		                        ring_matrix< bigint > & T, double y,
		                        sdigit factor)
{
	bigint_lattice tA(A);
	tA.lll_schnorr_euchner_fact(T, y, factor);
	return(tA);
}



//
// user defined Scalar Product
//
bigint_lattice lll_schnorr_euchner_fact(const bigint_lattice& A, double y,
		                        lattice_info& li, user_SP SP,
					sdigit factor)
{
	bigint_lattice tA(A);
	tA.lll_schnorr_euchner_fact(y, li, SP, factor);
	return(tA);
}



bigint_lattice lll_schnorr_euchner_fact(const bigint_lattice& A,
		                        ring_matrix< bigint > & T, double y,
		                        lattice_info& li, user_SP SP, sdigit factor)
{
	bigint_lattice tA(A);
	tA.lll_schnorr_euchner_fact(T, y, li, SP, factor);
	return(tA);
}



bigint_lattice lll_schnorr_euchner_fact(const bigint_lattice& A, double y,
		                        user_SP SP, sdigit factor)
{
	bigint_lattice tA(A);
	tA.lll_schnorr_euchner_fact(y, SP, factor);
	return(tA);
}



bigint_lattice lll_schnorr_euchner_fact(const bigint_lattice& A,
	                                ring_matrix< bigint > & T, double y,
		                        user_SP SP, sdigit factor)
{
	bigint_lattice tA(A);
	tA.lll_schnorr_euchner_fact(T, y, SP, factor);
	return(tA);
}



//
// Benne de Weger
//
// friend functions
//
bigint_lattice lll_benne_de_weger(const bigint_lattice& A,
				  sdigit y_nom, sdigit y_denom,
			  	  lattice_info& li)
{
	bigint_lattice tA(A);
	tA.lll_benne_de_weger(y_nom, y_denom, li);
	return(tA);
}



bigint_lattice lll_benne_de_weger(const bigint_lattice& A,
				  ring_matrix< bigint > & T,
				  sdigit y_nom, sdigit y_denom,
				  lattice_info& li)
{
	bigint_lattice tA(A);
	tA.lll_benne_de_weger(T, y_nom, y_denom, li);
	return(tA);
}



bigint_lattice lll_benne_de_weger(const bigint_lattice& A,
				  sdigit y_nom, sdigit y_denom)
{
	bigint_lattice tA(A);
	tA.lll_benne_de_weger(y_nom, y_denom);
	return(tA);
}



bigint_lattice lll_benne_de_weger(const bigint_lattice& A,
				  ring_matrix< bigint > & T,
				  sdigit y_nom, sdigit y_denom)
{
	bigint_lattice tA(A);
	tA.lll_benne_de_weger(T, y_nom, y_denom);
	return(tA);
}



//
// Buchmann - Kessler
//
bigint_lattice buchmann_kessler(const bigint_lattice& A,
			        ring_matrix< bigint > & T,
			        double y, lattice_info& li)
{
	bigint_lattice tA(A);
	tA.buchmann_kessler(T, y, li);
	return(tA);
}



bigint_lattice buchmann_kessler(const bigint_lattice& A,
			        ring_matrix< bigint > & T, double y)
{
	bigint_lattice tA(A);
	tA.buchmann_kessler(T, y);
	return(tA);
}



//
// Modified lll
//
// friend functions
//
bigint_lattice mlll(const bigint_lattice& A, double y,
		    bigint*& v, lattice_info& li)
{
	bigint_lattice tA(A);
	tA.mlll(y, v, li);
	return(tA);
}



bigint_lattice mlll(const bigint_lattice& A, double y,
		    base_vector< bigint > & v, lattice_info& li)
{
	bigint_lattice tA(A);
	tA.mlll(y, v, li);
	return(tA);
}



bigint_lattice mlll(const bigint_lattice& A, double y, bigint*& v)
{
	bigint_lattice tA(A);
	tA.mlll(y, v);
	return(tA);
}



bigint_lattice mlll(const bigint_lattice& A, double y, base_vector< bigint > & v)
{
	bigint_lattice tA(A);
	tA.mlll(y, v);
	return(tA);
}



//
// Other Things
//
// flag - checkings
//
bool bigint_lattice::check_basis()
{
	debug_handler("bigint_lattice", "check_basis()");
	dense_alg< bigint > da;
	bool flag;

	da.b.alg_trans = false;
	Tr_dense_create(da);
	flag = TrD_check_basis(da);
	Tr_dense_extract(da);
	if (flag)
		set_basis_flag();
	return(flag);
}



bool bigint_lattice::check_gram()
{
	debug_handler("bigint_lattice", "check_gram()");
	dense_alg< bigint > da;
	bool flag;

	da.b.alg_trans = false;
	Tr_dense_create(da);
	flag = TrD_check_gram(da);
	Tr_dense_extract(da);
	if (flag)
		set_gram_flag();
	return(flag);
}



//
// lll - checkings
//
bool bigint_lattice::lll_check(sdigit nom, sdigit denom)
{
	debug_handler("bigint_lattice", "lll_check(nom, denom)");
	dense_alg< bigint > da;
	bool red_flag;

	if (chk_gram())
		lidia_error_handler("bigint_lattice", "lll_check(nom, denom) "
				    ":: not avaidable for gram");

	da.b.y = static_cast<double>(nom)/static_cast<double>(denom);
	da.b.alg_trans = false;
	if ((da.b.y > 1.0) || (da.b.y < 0.5)) {
		lidia_warning_handler("bigint_lattice", "lll_check(nom, denom) :: "
				      "no allowed y for schnorr - euchner - lll");
		return(false);
	}
	Tr_dense_create(da);
	red_flag = TrD_lll_check(da, nom, denom);
	Tr_dense_extract(da);
	return(red_flag);
}



void bigint_lattice::lll_check_search(sdigit& nom, sdigit& denom)
{
	debug_handler("bigint_lattice", "lll_check_search(nom, denom)");
	dense_alg< bigint > da;

	if (chk_gram())
		lidia_error_handler("bigint_lattice", "lll_check_search(nom, denom) "
				    ":: not avaidable for gram");

	da.b.alg_trans = false;
	Tr_dense_create(da);
	TrD_lll_check_search(da, nom, denom);
	Tr_dense_extract(da);
}



//
// Creating needed structure for dense lattice algorithms
//
// Reducing 2 matrix assignments in performing transposition
// 2*rows*columns real bigint (> 4 Bytes) assignments to
// 2*rows*columns*3 pointer assignments (4 Bytes)
//
void bigint_lattice::Tr_dense_create(dense_alg< bigint > & da)
{
	debug_handler("bigint_lattice", "Tr_dense_create(da)");

	// Hier soll spaeter die Abfrage nach dem entspr. Bit hin !!
	if (!(da.b.transpose = chk_trans())) {
		//
		// Copy pointer only, no trans needed
		//
		da.b.rows = rows;
		da.b.columns = columns;
		if (da.b.alg_trans)
			da.b.real_columns = da.b.rows+da.b.columns;
		else
			da.b.real_columns = da.b.columns;
		da.b.rank = da.b.rows;
		resize(da.b.rows, da.b.real_columns);
		da.s.value = value;
		//
		// If trans - version, concate identic matrix
		//
		for (lidia_size_t i = da.b.columns; i < da.b.real_columns; i++)
			da.s.value[i-da.b.columns][i].assign_one();
	}
	else {
		//
		// Transpose and swap matrix, structure of bigint_lattice will be
		// destroyed
		//
		da.b.rows = columns;
		da.b.columns = rows;
		if (da.b.alg_trans)
			da.b.real_columns = da.b.rows+da.b.columns;
		else
			da.b.real_columns = da.b.columns;
		da.b.rank = da.b.rows;

		//
		// Allocate space for transposed structure
		//
		da.s.value = new bigint*[da.b.rows];
		memory_handler(da.s.value, "bigint_lattice", "Tr_dense_create(da) :: "
			       "not enough memory !");
		da.s.delvalue = new bigint[da.b.rows*da.b.real_columns];
		memory_handler(da.s.delvalue, "bigint_lattice", "Tr_dense_create(da) ::"
			       "not enough memory !");

		for (lidia_size_t i = 0; i < da.b.rows; i++)
			da.s.value[i] = &da.s.delvalue[i*da.b.real_columns];
		for (lidia_size_t i = 0; i < da.b.rows; i++)
			for (lidia_size_t j = 0; j < da.b.columns; j++)
				LiDIA::swap(value[j][i], da.s.value[i][j]);
		//
		// If trans - version, concate identic matrix
		//
		for (lidia_size_t i = da.b.columns; i < da.b.real_columns; i++)
			da.s.value[i-da.b.columns][i].assign_one();
	}
}



//
// Extracting lattice after performing a dense lattice algorithm
// (See Tr_dense_create() above)
//
void bigint_lattice::Tr_dense_extract(const dense_alg< bigint > & da)
{
	debug_handler("bigint_lattice", "Tr_dense_create(da)");
	bigint **Taddr;

	// Hier soll spaeter die Abfrage nach dem entspr. Bit hin !!
	if (!da.b.transpose) {
		value = da.s.value;
		//
		// Trans matrix computed ->
		// extract and store into da.TMatrix
		//
		if (da.b.alg_trans) {
			(da.s.TMatrix)->resize(da.b.rows, da.b.rows);
			Taddr = (da.s.TMatrix)->get_data_address();
			for (lidia_size_t i = 0; i < da.b.rows; i++)
				for (lidia_size_t j = 0; j < da.b.rows; j++)
					LiDIA::swap(Taddr[i][j], da.s.value[i][j+da.b.columns]);
		}
		resize(da.b.rows, da.b.columns);
	}
	else {
		resize(da.b.columns, da.b.rows);
		for (lidia_size_t i = 0; i < da.b.rows; i++)
			for (lidia_size_t j = 0; j < da.b.columns; j++)
				LiDIA::swap(da.s.value[i][j], value[j][i]);
		//
		// Trans matrix computed ->
		// extract and store transposed into da.TMatrix
		//
		if (da.b.alg_trans) {
			(da.s.TMatrix)->resize(da.b.rows, da.b.rows);
			Taddr = (da.s.TMatrix)->get_data_address();
			for (lidia_size_t i = 0; i < da.b.rows; i++)
				for (lidia_size_t j = 0; j < da.b.rows; j++)
					LiDIA::swap(Taddr[j][i], da.s.value[i][j+da.b.columns]);
		}
		//
		// Free storage allocated by Tr_dense_create !!
		// Every Tr_dense_create is followed by a Tr_dense_extract
		//
		delete[] da.s.delvalue;
		delete[] da.s.value;
	}
}



void bigint_lattice::randomize_vectors()
{
	debug_handler("bigint_lattice", "randomize_vectors()");
	//
	// First step is to transpose the lattice
	// Second is to generate a permutation of the rows (using random - functions)
	// Third to transpose again
	//
	dense_alg< bigint > da;
	da.b.alg_trans = false;
	Tr_dense_create(da);
	TrD_randomize_vectors(da);
	Tr_dense_extract(da);
}



void bigint_lattice::sort_big_vectors()
{
	debug_handler("bigint_lattice", "sort_big_vectors()");
	//
	// First step is to transpose the lattice
	// Second is to sort the rows by length (biggest first)
	// Third to transpose again
	//
	dense_alg< bigint > da;
	da.b.alg_trans = false;
	Tr_dense_create(da);
	TrD_sort_vectors(da, 1);
	Tr_dense_extract(da);
}



void bigint_lattice::sort_small_vectors()
{
	debug_handler("bigint_lattice", "sort_small_vectors()");
	//
	// First step is to transpose the lattice
	// Second is to sort the rows by length (smallest first)
	// Third to transpose again
	//
	dense_alg< bigint > da;
	da.b.alg_trans = false;
	Tr_dense_create(da);
	TrD_sort_vectors(da, -1);
	Tr_dense_extract(da);
}



void bigint_lattice::sort_vectors(bin_cmp_func cpf)
{
	debug_handler("bigint_lattice", "sort_vectors(cpf)");
	//
	// First step is to transpose the lattice
	// Second is to sort the rows by the given compare functions cpf
	// (See header file for more information)
	// Third to transpose again
	//
	dense_alg< bigint > da;
	da.b.alg_trans = false;
	Tr_dense_create(da);
	TrD_sort_vectors(da, cpf);
	Tr_dense_extract(da);
}



//
// Routine needed by another variation of the schnorr - euchner - lll
// try to reduced lattices with bigger entries with less precision
// for the approximation (try to use doubles)
//
sdigit bigint_lattice::TrD_search_factor(const dense_alg< bigint > & da)
{
	debug_handler("bigint_lattice", "TrD_search_factor(da)");
	sdigit min_bits = da.s.value[0][0].bit_length();
	sdigit max_bits = 0;
	sdigit bi_bit_len;

	for (lidia_size_t i = 0; i < da.b.rows; i++)
		for (lidia_size_t j = 0; j < da.b.columns; j++) {
			if ((bi_bit_len = da.s.value[i][j].bit_length()) < min_bits)
				min_bits = bi_bit_len;
			if (bi_bit_len > max_bits)
				max_bits = bi_bit_len;
		}
	return ((max_bits+min_bits) << 1);
}



//
// Algorithm`s
//
//
// Real implementation of the algorithms
// They are working on the transposed lattice
//
void bigint_lattice::TrD_randomize_vectors(dense_alg< bigint > & da)
{
	debug_handler("bigint_lattice", "TrD_randomize_vectors(da)");
	random_generator rg;
	char *bitvect;
	sdigit *perm;
	sdigit ran;
	bigint **temp;


	//
	// Allocate memory
	//
	bitvect = new char[da.b.rows];
	memory_handler(bitvect, "bigint_lattice", "TrD_randomize_vectors(da) :: "
		       "not enough memory !");
	perm = new sdigit[da.b.rows];
	memory_handler(perm, "bigint_lattice", "TrD_randomize_vectors(da) :: "
		       "not enough memory !");
	temp = new bigint*[da.b.rows];
	memory_handler(temp, "bigint_lattice", "TrD_randomize_vectors(da) :: "
		       "not enough memory !");

	//
	// Clear bitvector
	//
	for (lidia_size_t i = 0; i < da.b.rows; bitvect[i++] = 0);
	//
	// Generate a permutation
	// Try rows times to find valid value
	//
	for (lidia_size_t i = 0; i < da.b.rows; i++) {
		for (lidia_size_t j = 0; j < da.b.rows; j++) {
			rg >> ran;
			ran = ran%da.b.rows;
			if (bitvect[ran] == 0) {
				bitvect[ran] = 1;
				perm[i] = ran;
				break;
			}
			else
				if (j == rows-1) {
					for (lidia_size_t l = 0; l < da.b.rows; l++)
						if (bitvect[l] == 0) {
							perm[i] = l;
							bitvect[l] = 1;
							break;
						}
					break;
				}
		}
	}

	//
	// Perform permutation on lattice
	//
	for (lidia_size_t i = 0; i < da.b.rows; i++)
		temp[perm[i]] = da.s.value[i];
	for (lidia_size_t i = 0; i < da.b.rows; i++)
		da.s.value[i] = temp[i];

	//
	// Free allocated storage
	//
	delete[] perm;
	delete[] bitvect;
	delete[] temp;
}



void bigint_lattice::TrD_sort_vectors(dense_alg< bigint > & da, sdigit vgl)
{
	debug_handler("bigint_lattice", "TrD_sort_vectors(da, vgl)");
	p_vector< bigint > vector;
	bigint *quads;

	//
	// Allocate memory for scalar product of the vectors
	//
	quads = new bigint[da.b.rows];
	memory_handler(quads, "bigint_lattice", "TrD_sort_vectors(da, vgl) :: "
		       "not enough memory !");

	//
	// Compute scalar products
	//
	vector.vectsize = da.b.columns;
	for (lidia_size_t i = 0; i < da.b.rows; i++)
		vector.scalprod(quads[i], da.s.value[i], da.s.value[i]);
	//
	// sort by scalar products ("length^2")
	// vgl = -1   ->   smallest vector first
	// vgl =  1   ->   biggest vector first
	//
	for (lidia_size_t i = 0; i < da.b.rows-1; i++)
		for (lidia_size_t j = i+1; j < da.b.rows; j++) {
			if (quads[i].compare(quads[j]) == vgl) {
				vector.swap(da.s.value[i], da.s.value[j]);
				LiDIA::swap(quads[i], quads[j]);
			}
		}
	//
	// Free allocated storage
	//
	delete[] quads;
}



void bigint_lattice::TrD_sort_vectors(dense_alg< bigint > & da, bin_cmp_func cpf)
{
	debug_handler("bigint_lattice", "TrD_sort_vectors(da, cpf)");
	p_vector< bigint > vector;
	//
	// Sort vectors by specified compare function cpf
	// Perform bubble sort (for typical sizes of the
	// lattice quicksort not nessesary)
	//
	vector.vectsize = da.b.columns;
	for (lidia_size_t i = 0; i < da.b.rows-1; i++)
		for (lidia_size_t j = i+1; j < da.b.rows; j++)
			if (cpf(da.s.value[i], da.s.value[j], da.b.columns) > 0)
				vector.swap(da.s.value[i], da.s.value[j]);
}



bool bigint_lattice::TrD_gram_schmidt_orth_bdw(const dense_alg< bigint > & da,
					       bigint** my, bigint** gso,
					       bigint* v)
{
	debug_handler("bigint_lattice", "TrD_gram_schmidt_orth_bdw(da, my, gso, v)");

	bigint* tempPbin0;
	bigint* tempPbin1;
	p_vector< bigint > vector;

	if (da.b.rows > da.b.columns)
		return(false);

	tempPbin0 = new bigint[da.b.rows+da.b.columns+1];
	memory_handler(tempvect0, "bigint_lattice", "TrD_gram_schmidt_orth_bdw(da, my, "
		       " gso, v) :: not enough memory !");
	tempPbin1 = &tempPbin0[da.b.rows+1];

	vector.vectsize = da.b.rows+1;
	vector.assign_zero(tempPbin0);
	tempPbin0[0].assign_one();
	vector.vectsize = da.b.columns;

	for (lidia_size_t i = 0; i < da.b.rows; i++) {
		vector.assign(gso[i], da.s.value[i]);
		for (lidia_size_t j = 0; j < i; j++) {
			vector.scalprod(my[j][i], da.s.value[i], gso[j]);
			for (lidia_size_t cx = 0; cx < da.b.columns; cx++) {
				LiDIA::multiply(gso[i][cx], tempPbin0[j+1], gso[i][cx]);
				LiDIA::multiply(tempPbin1[cx], my[j][i], gso[j][cx]);
				LiDIA::subtract(gso[i][cx], gso[i][cx], tempPbin1[cx]);
			}
			if (tempPbin0[j].is_zero()) {
				delete[] tempPbin0;
				return(false);
			}
			for (lidia_size_t m = 0; m < da.b.columns; m++)
				LiDIA::divide(gso[i][m], gso[i][m], tempPbin0[j]);
		}
		vector.scalprod(tempPbin0[i+1], gso[i], gso[i]);
		if (tempPbin0[i].is_zero()) {
			delete[] tempPbin0;
			return(false);
		}
		else
			LiDIA::divide(tempPbin0[i+1], tempPbin0[i+1], tempPbin0[i]);
	}
	vector.vectsize = da.b.rows+1;
	vector.assign(v, tempPbin0);
	delete[] tempPbin0;
	return(true);
}



bool bigint_lattice::TrD_lll_check(const dense_alg< bigint > & da,
				   sdigit y_nom, sdigit y_denom)
{
	debug_handler("bigint_lattice", "Tr_lll_check(da, y_nom, y_denom)");
	bigint tempbin0;
	bigint tempbin1;
	bigint tempbin2;
	bigint* tempPbin0;
	bigint** my;
	bigint* mydel;
	bigint** gso;
	p_vector< bigint > vector;

	my = new bigint*[da.b.rows*2];
	memory_handler(my, "bigint_lattice", "Tr_lll_check(da, y_nom, y_denom) :: "
		       "not enough memory !");
	mydel = new bigint[da.b.rows*(da.b.rows+da.b.columns)];
	memory_handler(mydel, "bigint_lattice", "Tr_lll_check(da, y_nom, y_denom) :: "
		       "not enough memory !");
	gso = &my[da.b.rows];
	for (lidia_size_t i = 0; i < da.b.rows; i++) {
		my[i] = &mydel[i*da.b.rows];
		gso[i] = &mydel[da.b.rows*da.b.rows+i*da.b.columns];
	}

	tempPbin0 = new bigint[da.b.rows+1];
	memory_handler(tempvect0, "bigint_lattice", "Tr_lll_check(da, y_nom, "
		       "y_denom) :: not enough memory !");

	if (!TrD_gram_schmidt_orth_bdw(da, my, gso, tempPbin0)) {
		delete[] mydel;
		delete[] my;
		delete[] tempPbin0;
		lidia_error_handler("bigint_lattice", "Tr_lll_check(da, y_nom, "
				    "y_denom) :: lattice is no basis !");
		return (false);
	}

	for (lidia_size_t i = 0; i < da.b.rows; i++)
		for (lidia_size_t j = 0; j < i; j++) {
			tempbin0.assign(abs(my[j][i]));
			tempbin0.multiply_by_2();
			if (tempbin0 > tempPbin0[j+1]) {
				delete[] mydel;
				delete[] my;
				delete[] tempPbin0;
				return (false);
			}
		}
	vector.vectsize = da.b.columns;
	for (lidia_size_t i = 1; i < da.b.rows; i++) {
		LiDIA::square(tempbin0, tempPbin0[i]);
		LiDIA::multiply(tempbin1, tempbin0, y_nom);

		LiDIA::square(tempbin0, my[i-1][i]);
		LiDIA::multiply(tempbin2, tempbin0, y_denom);
		LiDIA::subtract(tempbin0, tempbin1, tempbin2);
		LiDIA::multiply(tempbin1, tempPbin0[i-1], tempPbin0[i+1]);
		LiDIA::multiply(tempbin2, tempbin1, y_denom);

		if (tempbin2 < tempbin0) {
			delete[] mydel;
			delete[] my;
			delete[] tempPbin0;
			return(false);
		}

	}
	delete[] mydel;
	delete[] my;
	delete[] tempPbin0;
	return (true);
}



void bigint_lattice::TrD_lll_check_search(const dense_alg< bigint > & da,
					  sdigit& a, sdigit& b)
{
	debug_handler("bigint_lattice", "Tr_lll_check_search(da, a, b)");
	bool success;
	sdigit param_nom = 3;
	sdigit param_den = 4;
	sdigit upper_nom = 4;
	sdigit upper_den = 4;
	sdigit downer_nom = 2;
	sdigit downer_den = 4;
	bigint tempbin0;
	bigint tempbin1;
	bigint tempbin2;
	bigint tempbin3;
	bigint* tempPbin0;
	bigint** my;
	bigint* mydel;
	bigint** gso;
	p_vector< bigint > vector;


	my = new bigint*[da.b.rows*2];
	memory_handler(my, "bigint_lattice", "Tr_lll_check_search(da, a, b) :: "
		       "not enough memory !");
	mydel = new bigint[da.b.rows*(da.b.rows+da.b.columns)];
	memory_handler(mydel, "bigint_lattice", "Tr_lll_check_search(da, a, b) :: "
		       "not enough memory !");
	gso = &my[da.b.rows];
	for (lidia_size_t i = 0; i < da.b.rows; i++) {
		my[i] = &mydel[i*da.b.rows];
		gso[i] = &mydel[da.b.rows*da.b.rows+i*da.b.columns];
	}

	tempPbin0 = new bigint[da.b.rows+1];
	memory_handler(tempPbin0, "bigint_lattice", "Tr_lll_check_search(da, "
		       "a, b) :: not enough memory !");

	if (!TrD_gram_schmidt_orth_bdw(da, my, gso, tempPbin0)) {
		delete[] mydel;
		delete[] my;
		delete[] tempPbin0;
		lidia_error_handler("bigint_lattice", "Tr_lll_check_search(da, a, b) :: "
				    "lattice is no basis !");
		return;
	}

	a = 0;
	b = 1;
	for (lidia_size_t i = 0; i < da.b.rows; i++)
		for (lidia_size_t j = 0; j < i; j++) {
			tempbin1.assign(abs(my[j][i]));
			tempbin1.multiply_by_2();
			if (tempbin1 > tempPbin0[j+1]) {
				delete[] mydel;
				delete[] my;
				delete[] tempPbin0;
				return;
			}
		}
	vector.vectsize = da.b.columns;
	for (lidia_size_t j = 0; j < 20; j++)   // 2^-20
	{
		success = true;
		for (lidia_size_t i = 1; i < da.b.rows; i++) {
			LiDIA::square(tempbin0, tempPbin0[i]);
			LiDIA::multiply(tempbin1, tempbin0, param_nom);

			LiDIA::square(tempbin0, my[i-1][i]);
			LiDIA::multiply(tempbin2, tempbin0, param_den);
			LiDIA::subtract(tempbin0, tempbin1, tempbin2);
			LiDIA::multiply(tempbin1, tempPbin0[i-1], tempPbin0[i+1]);
			LiDIA::multiply(tempbin2, tempbin1, param_den);

			if (tempbin2 < tempbin0) {
				success = false;
				break;
			}

		}
		if (success) {
			a = param_nom;
			b = param_den;

			downer_nom = param_nom;
			downer_den = param_den << 1;
			upper_den <<= 1;
			param_nom = (downer_nom+upper_nom);
			param_den <<= 1;
			upper_nom <<= 1;
			downer_nom <<= 1;
		}
		else {
			upper_nom = param_nom;
			upper_den = param_den << 1;
			downer_den <<= 1;
			param_nom = (upper_nom+downer_nom);
			param_den <<= 1;
			upper_nom <<= 1;
			downer_nom <<= 1;
		}
	}
	delete[] mydel;
	delete[] my;
	delete[] tempPbin0;
}



bool bigint_lattice::TrD_check_basis(const dense_alg< bigint > & da)
{
	debug_handler("bigint_lattice", "TrD_check_basis(da)");
	bigint tempbin0;
	matrix< bigint > Tbin0(da.b.rows, da.b.columns, const_cast<const bigint **>(da.s.value));
	matrix< bigint > Tbin1;

	if (da.b.rows > da.b.columns)
		return (false);
	if (da.b.rows == da.b.columns) {
		tempbin0.assign(Tbin0.det());
		if (tempbin0.is_zero())
			return(false);
		else
			return(true);
	}
	// Tbin1.assign(Tbin0*Tbin0.trans());
	Tbin1.assign(Tbin0.trans());
	LiDIA::multiply(Tbin1, Tbin0, Tbin1);
	// A * A^T


	tempbin0.assign(Tbin1.det());
	if (tempbin0.is_zero())
		return(false);
	else
		return(true);

	return(false);
}



bool bigint_lattice::TrD_check_gram(const dense_alg< bigint > & da)
{
	debug_handler("bigint_lattice", "TrD_check_gram(da)");
	bigint tempbin0;
	matrix< bigint > Tbin0(da.b.rows, da.b.columns, const_cast<const bigint **>(da.s.value));


	//
	// quadratic
	//
	if (da.b.rows != da.b.columns)
		return(false);

	//
	// symmetric
	//
	for (lidia_size_t i = 0; i < da.b.rows; i++)
		for (lidia_size_t j = i; j < da.b.columns; j++)
			if (da.s.value[i][j] != da.s.value[j][i])
				return(false);
	//
	// positiv definit, funktioniert nicht bei gensys
	//
	for (lidia_size_t i = da.b.rows; i > 0; i--) {
		Tbin0.set_no_of_rows(i);
		Tbin0.set_no_of_columns(i);
		tempbin0.assign(Tbin0.det());
		//      if (tempbin0.is_le_zero())
		if (tempbin0.is_lt_zero())
			return(false);
	}
	return(true);
}



sdigit bigint_lattice::compute_read_precision(const dense_alg< bigint > & da)
{
	debug_handler("bigint_lattice", "compute_read_precision(da)");
	//
	// Give the digit size of longest value
	//
	sdigit new_prec;
	sdigit read_prec = 0;
	for (lidia_size_t x = 0; x < da.b.rows; ++x)
		for (lidia_size_t y = 0; y < da.b.columns; ++y)
			if (read_prec < (new_prec = static_cast<sdigit>(da.s.value[x][y].bit_length()/
									std::log(10.0)+1)))
				read_prec = new_prec;
	return (read_prec);
}



sdigit bigint_lattice::compute_precision(const dense_alg< bigint > & da)
{
	debug_handler("bigint_lattice", "compute_precision(da)");
	bigfloat alpha, zweipotq;
	sdigit read_prec = compute_read_precision(da);
	sdigit n2 = da.b.rows;
	sdigit prec;
	alpha_compute(da, alpha);
	zwei_pot_q_compute(da, zweipotq, n2, alpha);
	read_prec = 2*read_prec+da.b.rows-1;
	prec = static_cast<sdigit>(2*(((zweipotq.bit_length()+1)*
				       std::log(2.0)/std::log(10.0)+1)+read_prec)+(da.b.columns+da.b.rows)-1);
	return (prec);
}



void bigint_lattice::gamma_compute(bigfloat& g, sdigit l)
{
	debug_handler("bigint_lattice", "gamma_compute(g, l)");
	bigfloat ha[] = {1, 4.0/3.0, 2, 4, 8, 64.0/3.0, 64, 256};
	//
	// Computation of the l-th Hermite - Constant gamma,
	//
	bigfloat lg;
	if ((l > 0) && (l < 9))
		g.assign(ha[l-1]);
	else {
		lg.assign(0.75);
		LiDIA::log(lg, lg);
		g.assign(static_cast<sdigit>(l * (l-1)));
		g.divide_by_2();
		LiDIA::multiply(lg, lg, g);
		LiDIA::exp(g, lg);
	}
}



void bigint_lattice::alpha_compute(const dense_alg< bigint > & da, bigfloat& alpha)
{
	//
	// alpha = max{vect[1].l2_norm(),...,vect[columns].l2_norm()}
	// norm = vect[i].l2_norm(), 0 <= i < columns
	//
	debug_handler("bigint_lattice", "alpha_compute(da, alpha)");
	bigint tempbin0;
	bigint bi_alpha;
	p_vector< bigint > vector;
	vector.vectsize = da.b.columns;
	vector.scalprod(bi_alpha, da.s.value[0], da.s.value[0]);
	for (lidia_size_t i = 1; i < da.b.rows; ++i) {
		vector.scalprod(tempbin0, da.s.value[i], da.s.value[i]);
		if (bi_alpha.compare(tempbin0) < 0)
			bi_alpha.assign(tempbin0);
	}
	//
	// calculating sqrt of l2_norm
	//
	alpha.assign(bi_alpha);
	LiDIA::sqrt(alpha, alpha);
}



void bigint_lattice::zwei_pot_q_compute(const dense_alg< bigint > & da,
					bigfloat& zweipotq,
					sdigit& n2, bigfloat& alpha)
{
	debug_handler("bigint_lattice", "zwei_pot_q_compute(da, zwpq, n2, alpha)");
	sdigit beta = da.b.rows;
	bigint tempbin0;
	bigfloat tempbfl0;
	bigfloat tempbfl1;
	bigfloat tempbfl2;

	// Computation of beta = min (A.columns, A.rows)
	if (da.b.columns < da.b.rows)
		beta = da.b.columns;

	tempbfl0.assign(n2);
	LiDIA::sqrt(tempbfl0, tempbfl0);
	tempbfl0.divide_by_2();
	LiDIA::multiply(tempbfl0, rows, tempbfl0);
	tempbfl1.assign(rows);
	LiDIA::sqrt(tempbfl1, tempbfl1);
	LiDIA::add(tempbfl0, tempbfl0, tempbfl1);

	LiDIA::log(tempbfl1, alpha);
	LiDIA::multiply(tempbfl1, tempbfl1, beta);
	LiDIA::exp(tempbfl1, tempbfl1);
	gamma_compute(tempbfl2 , beta);

	LiDIA::sqrt(tempbfl2, tempbfl2);
	LiDIA::multiply(tempbfl2, tempbfl2, tempbfl1);
	LiDIA::multiply(tempbfl2, tempbfl2, tempbfl0);
	tempbin0.assign(n2);
	LiDIA::multiply(tempbfl0, tempbfl0, rows);
	LiDIA::sqrt(tempbfl0, tempbin0);
	LiDIA::add(tempbfl0, tempbfl0, 2);
	LiDIA::multiply(zweipotq, tempbfl0, tempbfl2);

	tempbfl0.assign(beta + rows);
	tempbfl0.divide_by_2();
	LiDIA::subtract(tempbfl0, tempbfl0, 1);
	tempbfl0.multiply_by_2();
	LiDIA::exp(tempbfl0, tempbfl0);
	LiDIA::multiply(tempbfl0, zweipotq, tempbfl0);

	tempbfl2.assign(beta + 1);
	tempbfl2.divide_by_2();
	tempbfl1.assign(columns);
	LiDIA::log(tempbfl1, tempbfl1);
	LiDIA::multiply(tempbfl2, tempbfl2, tempbfl1);
	LiDIA::exp(tempbfl2, tempbfl2);
	LiDIA::divide(tempbfl0, tempbfl0, tempbfl2);
	LiDIA::ceil(zweipotq, tempbfl0);
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
