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
#include	"LiDIA/bigfloat_lattice.h"



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



bigint_lattice bigfloat_lattice::CallBiL;

//
// Set special lattice characteristics
//
void bigfloat_lattice::set_red_orientation_columns()
{
	debug_handler("bigfloat_lattice", "set_red_orentation_columns()");
	bitfield.lattice_mode &= !REDUCE_ROWS;
}



void bigfloat_lattice::set_red_orientation_rows()
{
	debug_handler("bigfloat_lattice", "set_red_orentation_rows()");
	bitfield.lattice_mode |= REDUCE_ROWS;
}



void bigfloat_lattice::set_gram_flag()
{
	debug_handler("bigfloat_lattice", "set_gram_flag()");
	bitfield.info_mode |= GRAM_MATRIX;
}



void bigfloat_lattice::set_basis_flag()
{
	debug_handler("bigfloat_lattice", "set_basis_flag()");
	if (bitfield.lattice_mode & REDUCE_ROWS)
		bitfield.structure_mode |= ROWS_LININD;
	else
		bitfield.structure_mode |= COLUMNS_LININD;
}



void bigfloat_lattice::delete_gram_flag()
{
	debug_handler("bigfloat_lattice", "delete_gram_flag()");
	bitfield.info_mode &= !GRAM_MATRIX;
}



void bigfloat_lattice::delete_basis_flag()
{
	debug_handler("bigfloat_lattice", "delete_basis_flag()");
	if (bitfield.lattice_mode & REDUCE_ROWS)
		bitfield.structure_mode &= !ROWS_LININD;
	else
		bitfield.structure_mode &= !COLUMNS_LININD;
}



//
// get the special lattice characteristics
//
bool bigfloat_lattice::get_red_orientation()
{
	debug_handler("bigfloat_lattice", "get_red_orentation()");
	return ((bitfield.lattice_mode&REDUCE_ROWS)?false:true);
	//  return (!((bool )(bitfield.lattice_mode&REDUCE_ROWS)));
}



bool bigfloat_lattice::get_basis_flag()
{
	debug_handler("bigfloat_lattice", "get_basis_flag()");
	if (bitfield.lattice_mode & REDUCE_ROWS)
		return ((bitfield.structure_mode&ROWS_LININD)?true:false);
	else
		return ((bitfield.structure_mode&COLUMNS_LININD)?true:false);
}



bool bigfloat_lattice::get_gram_flag()
{
	debug_handler("bigfloat_lattice", "get_gram_flag()");
	return ((bitfield.structure_mode&GRAM_MATRIX)?true:false);
}



bool bigfloat_lattice::chk_basis()
{
	debug_handler("bigfloat_lattice", "chk_basis()");
	if (bitfield.lattice_mode & REDUCE_ROWS)
		return ((bitfield.structure_mode&ROWS_LININD)?true:false);
	else
		return ((bitfield.structure_mode&COLUMNS_LININD)?true:false);
}



bool bigfloat_lattice::chk_gram()
{
	debug_handler("bigfloat_lattice", "chk_gram()");
	return ((bitfield.info_mode&GRAM_MATRIX)?true:false);
	//  return ((bool )(bitfield.info_mode&GRAM_MATRIX));
}



bool bigfloat_lattice::chk_trans()
{
	debug_handler("bigfloat_lattice", "chk_trans()");
	if (chk_gram())
		return (false);
	else {
		if (((bitfield.lattice_mode&REDUCE_ROWS) && (bitfield.structure_mode&COLUMN_ORIENTED)) ||
		    (!(bitfield.lattice_mode&REDUCE_ROWS) && !(bitfield.structure_mode&COLUMN_ORIENTED)))
			return (true);
		else
			return (false);
	}
	//    return (((bool )(bitfield.lattice_mode&REDUCE_ROWS)) ==
	//	      ((bool )(bitfield.structure_mode&COLUMN_ORIENTED)));

}



bool bigfloat_lattice::chk_reduce_columns()
{
	debug_handler("bigfloat_lattice", "chk_reduce_columns()");
	return ((bitfield.lattice_mode&REDUCE_ROWS)?false:true);
	//  return(!((bool )(bitfield.lattice_mode&REDUCE_ROWS)));
}



//
// Dimension checking
//
bool bigfloat_lattice::chk_mlll_dim()
{
	debug_handler("bigfloat_lattice", "chk_mlll_dim()");
	return (rows+static_cast<lidia_size_t>(chk_reduce_columns()) ==
		columns+static_cast<lidia_size_t>(!chk_reduce_columns()));
}



bool bigfloat_lattice::chk_lll_dim()
{
	debug_handler("bigfloat_lattice", "chk_lll_dim()");
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
void bigfloat_lattice::chk_corr_param(double& y)
{
	debug_handler("bigfloat_lattice", "chk_corr_param(y)");
	if ((y > 1.0) || (y <= 0.5))
		y = 0.99;
}



void bigfloat_lattice::chk_corr_param(sdigit& nom, sdigit& denom)
{
	debug_handler("bigfloat_lattice", "chk_corr_param(nom, denom)");
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
void bigfloat_lattice::lll_schnorr_euchner_orig(double y, lattice_info& li,
						sdigit x_factor)
{
	debug_handler("bigfloat_lattice", "lll_schnorr_euchner_orig(y, li, factor)");
	ALG_DEF_BASIS(bigfloat, vector_op, Normal) alg_basis;
	ALG_DEF_GENSYS(bigfloat, vector_op, Normal) alg_gensys;

	if (!(chk_lll_dim()))
		lidia_error_handler("bigfloat_lattice", "lll_schnorr_euchner_orig(y, li, "
				    "factor) :: wrong dimension (no basis)");
	dense_alg< bigfloat > da;
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



void bigfloat_lattice::lll_schnorr_euchner_orig(ring_matrix< bigint > & T, double y,
			                        lattice_info& li, sdigit x_factor)
{
	debug_handler("bigfloat_lattice", "lll_schnorr_euchner_orig(T, y, li, factor)");
	field_matrix< bigfloat > tT;
	lll_schnorr_euchner_orig(tT, y, li, x_factor);
	conv_bfl_bin(tT, T);
}



void bigfloat_lattice::lll_schnorr_euchner_orig(field_matrix< bigfloat > & T, double y,
			                        lattice_info& li, sdigit x_factor)
{
	debug_handler("bigfloat_lattice", "lll_schnorr_euchner_orig(T, y, li, factor)");
	ALG_DEF_BASIS(bigfloat, vector_op, Normal) alg_basis;
	ALG_DEF_GENSYS(bigfloat, vector_op, Normal) alg_gensys;

	if (!(chk_lll_dim()))
		lidia_error_handler("bigfloat_lattice", "lll_schnorr_euchner_orig(T, y, li, "
				    "factor) :: wrong dimension (no basis)");
	dense_alg< bigfloat > da;
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



void bigfloat_lattice::lll_schnorr_euchner_fact(double y, lattice_info& li,
						sdigit x_factor)
{
	debug_handler("bigfloat_lattice", "lll_schnorr_euchner_fact(y, li, factor)");
	ALG_DEF_BASIS(bigint, vector_op, VariationI) alg_basis;
	ALG_DEF_GENSYS(bigint, vector_op, VariationI) alg_gensys;

	if (!(chk_lll_dim()))
		lidia_error_handler("bigfloat_lattice", "lll_schnorr_euchner_fact(y, li, "
				    "factor) :: wrong dimension (no basis)");
	dense_alg< bigint > da;
	da.b.y = y;
	chk_corr_param(da.b.y);
	da.b.alg_trans = false;
	x_factor = ((x_factor < 1) ? 1 : x_factor);
	da.d.bit_prec = x_factor*DOUBLE_MANTISSA_BITS;
	Tr_dense_create_f(da);
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
	Tr_dense_extract_f(da);
}



void bigfloat_lattice::lll_schnorr_euchner_fact(ring_matrix< bigint > & T, double y,
			                        lattice_info& li, sdigit x_factor)
{
	debug_handler("bigfloat_lattice", "lll_schnorr_euchner_fact(T, y, li, factor)");
	ALG_DEF_BASIS(bigint, vector_op, VariationI) alg_basis;
	ALG_DEF_GENSYS(bigint, vector_op, VariationI) alg_gensys;

	if (!(chk_lll_dim()))
		lidia_error_handler("bigfloat_lattice", "lll_schnorr_euchner_fact(T, y, li, "
				    "factor) :: wrong dimension (no basis)");
	dense_alg< bigint > da;
	da.b.y = y;
	chk_corr_param(da.b.y);
	da.b.alg_trans = true;
	da.s.TMatrix = &T;
	x_factor = ((x_factor < 1) ? 1 : x_factor);
	da.d.bit_prec = x_factor*DOUBLE_MANTISSA_BITS;
	Tr_dense_create_f(da);
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
	Tr_dense_extract_f(da);
}



void bigfloat_lattice::lll_schnorr_euchner_fact(field_matrix< bigfloat > & T, double y,
			                        lattice_info& li, sdigit x_factor)
{
	debug_handler("bigfloat_lattice", "lll_schnorr_euchner_fact(T, y, li, factor)");
	ring_matrix< bigint > tT;
	lll_schnorr_euchner_fact(tT, y, li, x_factor);
	conv_bin_bfl(tT, T);
}



//
// Schnorr - Euchner for user defined Scalar Product
//
void bigfloat_lattice::lll_schnorr_euchner_orig(double y, lattice_info& li,
			                        user_SP SP, sdigit x_factor)
{
	debug_handler("bigfloat_lattice", "lll_schnorr_euchner_orig(y, li, SP, factor)");
	ALG_DEF_BASIS(bigfloat, vector_op_SP, Normal) alg_basis;
	ALG_DEF_GENSYS(bigfloat, vector_op_SP, Normal) alg_gensys;

	ALG_POINTER(alg_basis, SP, bfl);
	ALG_POINTER(alg_gensys, SP, bfl);

	if (!(chk_lll_dim()))
		lidia_error_handler("bigfloat_lattice", "lll_schnorr_euchner_orig(y, li, SP, "
				    "factor) :: wrong dimension (no basis)");
	dense_alg< bigfloat > da;
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



void bigfloat_lattice::lll_schnorr_euchner_orig(ring_matrix< bigint > & T, double y,
			                        lattice_info& li, user_SP SP,
						sdigit x_factor)
{
	debug_handler("bigfloat_lattice", "lll_schnorr_euchner_orig(T, y, li, SP, factor)");
	field_matrix< bigfloat > tT;
	lll_schnorr_euchner_orig(tT, y, li, SP, x_factor);
	conv_bfl_bin(tT, T);
}



void bigfloat_lattice::lll_schnorr_euchner_orig(field_matrix< bigfloat > & T, double y,
			                        lattice_info& li, user_SP SP,
						sdigit x_factor)
{
	debug_handler("bigfloat_lattice", "lll_schnorr_euchner_orig(T, y, li, SP, factor)");
	ALG_DEF_BASIS(bigfloat, vector_op_SP, Normal) alg_basis;
	ALG_DEF_GENSYS(bigfloat, vector_op_SP, Normal) alg_gensys;

	ALG_POINTER(alg_basis, SP, bfl);
	ALG_POINTER(alg_gensys, SP, bfl);

	if (!(chk_lll_dim()))
		lidia_error_handler("bigfloat_lattice", "lll_schnorr_euchner_orig(T, y, li, SP, "
				    "factor) :: wrong dimension (no basis)");
	dense_alg< bigfloat > da;
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



void bigfloat_lattice::lll_schnorr_euchner_fact(double y, lattice_info& li,
			                        user_SP SP, sdigit x_factor)
{
	debug_handler("bigfloat_lattice", "lll_schnorr_euchner_fact(y, li, SP, factor)");
	ALG_DEF_BASIS(bigint, vector_op_SP, VariationI) alg_basis;
	ALG_DEF_GENSYS(bigint, vector_op_SP, VariationI) alg_gensys;

	ALG_POINTER(alg_basis, SP, bin);
	ALG_POINTER(alg_gensys, SP, bin);

	if (!(chk_lll_dim()))
		lidia_error_handler("bigfloat_lattice", "lll_schnorr_euchner_fact(y, li, SP, "
				    "factor) :: wrong dimension (no basis)");
	dense_alg< bigint > da;
	da.b.y = y;
	chk_corr_param(da.b.y);
	da.b.alg_trans = false;
	x_factor = ((x_factor < 1) ? 1 : x_factor);
	da.d.bit_prec = x_factor*DOUBLE_MANTISSA_BITS;
	Tr_dense_create_f(da);
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
	Tr_dense_extract_f(da);
}



void bigfloat_lattice::lll_schnorr_euchner_fact(ring_matrix< bigint > & T, double y,
			                        lattice_info& li, user_SP SP,
						sdigit x_factor)
{
	debug_handler("bigfloat_lattice", "lll_schnorr_euchner_fact(T, y, li, SP, factor)");
	ALG_DEF_BASIS(bigint, vector_op_SP, VariationI) alg_basis;
	ALG_DEF_GENSYS(bigint, vector_op_SP, VariationI) alg_gensys;

	ALG_POINTER(alg_basis, SP, bin);
	ALG_POINTER(alg_gensys, SP, bin);

	if (!(chk_lll_dim()))
		lidia_error_handler("bigfloat_lattice", "lll_schnorr_euchner_fact(T, y, li, SP, "
				    "factor) :: wrong dimension (no basis)");
	dense_alg< bigint > da;
	da.b.y = y;
	chk_corr_param(da.b.y);
	da.b.alg_trans = true;
	da.s.TMatrix = &T;
	x_factor = ((x_factor < 1) ? 1 : x_factor);
	da.d.bit_prec = x_factor*DOUBLE_MANTISSA_BITS;
	Tr_dense_create_f(da);
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
	Tr_dense_extract_f(da);
}



void bigfloat_lattice::lll_schnorr_euchner_fact(field_matrix< bigfloat > & T, double y,
			                        lattice_info& li, user_SP SP,
						sdigit x_factor)
{
	debug_handler("bigfloat_lattice", "lll_schnorr_euchner_fact(T, y, li, SP, factor)");
	ring_matrix< bigint > tT;
	lll_schnorr_euchner_fact(tT, y, li, SP, x_factor);
	conv_bin_bfl(tT, T);
}



//
// Buchmann - Kessler
//
void bigfloat_lattice::buchmann_kessler(ring_matrix< bigint > & T, double y,
				        lattice_info& li)
{
	debug_handler("bigint_lattice", "buchmann_kessler(T, y, li)");
	field_matrix< bigfloat > tT;
	buchmann_kessler(tT, y, li);
	conv_bfl_bin(tT, T);
}



void bigfloat_lattice::buchmann_kessler(field_matrix< bigfloat > & T, double y,
				        lattice_info& li)
{
	debug_handler("bigint_lattice", "buchmann_kessler(T, y, li)");
	if (chk_gram())
		lidia_error_handler("bigint_lattice", "buchmann_kessler(T, y, "
				    " li) :: not avaidable for gram");

	dense_alg< bigfloat > da;
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
void bigfloat_lattice::mlll(double y, bigint*& v, lattice_info& li)
{
	debug_handler("bigfloat_lattice", "mlll(y, v, li)");

	if (!chk_mlll_dim())
		lidia_error_handler("bigfloat_lattice", "mlll(y, v, li) ::"
				    " wrong dimension");

	dense_alg< bigfloat > da;
	da.b.y = y;
	chk_corr_param(da.b.y);
	da.b.alg_trans = false;
	Tr_dense_create(da);
	v = TrD_mlll_bfl(da, li);
	Tr_dense_extract(da);
}



void bigfloat_lattice::mlll(double y, bigfloat*& v, lattice_info& li)
{
	debug_handler("bigfloat_lattice", "mlll(y, v, li)");
	base_vector< bigint > bv;
	mlll(y, bv, li);
	v = new bigfloat[bv.size()];
	for (lidia_size_t i = 0; i < bv.size(); i++)
		v[i].assign(bv[i]);
}



void bigfloat_lattice::mlll(double y, base_vector< bigint > & bv, lattice_info& li)
{
	debug_handler("bigfloat_lattice", "mlll(y, bv, li)");
	if (!chk_mlll_dim())
		lidia_error_handler("bigfloat_lattice", "mlll(y, bv, li) ::"
				    " wrong dimension");

	dense_alg< bigfloat > da;
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



void bigfloat_lattice::mlll(double y, base_vector< bigfloat > & bv,
			    lattice_info& li)
{
	debug_handler("bigfloat_lattice", "mlll(y, bv, li)");
	base_vector< bigint > v;

	mlll(y, v, li);
	for (lidia_size_t i = 0; i < v.size(); i++)
		bv[i].assign(v[i]);
}



void bigfloat_lattice::close_vector(const base_vector< bigfloat > & v,
				    base_vector< bigfloat > & cv,
				    sdigit x_factor)
{
	debug_handler("bigfloat_lattice", "close_vector(v, cv)");
	ALG_DEF_BASIS(bigfloat, vector_op, Normal) alg_basis;
	bigfloat_lattice cA(*this);
	dense_alg< bigfloat > da;
	lattice_info li;
	p_vector< bigfloat > vector;
	bigfloat sqrtBmax;
	bigfloat C, B, Bmax;

	if ((!chk_basis()) || (chk_gram()))
		lidia_error_handler("bigfloat_lattice", "close_vector(c, cv) :: "
				    "not implemented for gram or gensys !");
	cA.set_no_of_rows(get_no_of_rows()+1);
	cA.set_no_of_columns(get_no_of_columns()+1);
	da.b.alg_trans = false;
	da.b.y = 0.99;
	x_factor = ((x_factor < 1) ? 1 : x_factor);
	da.d.bit_prec = x_factor*DOUBLE_MANTISSA_BITS;
	cA.Tr_dense_create(da);
	vector.vectsize = da.b.columns;
	//
	// Berechne C vereinfacht
	// aus sqrt(3)/2 wird 1, wegen >=
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
		lidia_error_handler("bigfloat_lattice", "close_vector(c, cv) :: "
				    "illegal size of vector !");
	for (lidia_size_t i = 0; i < da.b.columns-1; i++)
		da.s.value[da.b.rows-1][i].assign(v[i]);
	da.s.value[da.b.rows-1][da.b.columns-1].assign(C);
	ALG_CALL(alg_basis, lll, da, li, x_factor)
		for (lidia_size_t i = 0; i < da.b.columns-1; i++)
			LiDIA::subtract(da.s.value[da.b.rows-1][i], v[i], da.s.value[da.b.rows-1][i]);
	cv.set_size(da.b.columns-1);
	cv.set_data(da.s.value[da.b.rows-1], da.b.columns-1);
	cA.Tr_dense_extract(da);
}



//
// Interface to Algorithms
//
// friend functions
//
bigfloat_lattice lll_schnorr_euchner_orig(const bigfloat_lattice& A, double y,
		                          lattice_info& li, sdigit factor)
{
	bigfloat_lattice tA(A);
	tA.lll_schnorr_euchner_orig(y, li, factor);
	return(tA);
}



bigfloat_lattice lll_schnorr_euchner_orig(const bigfloat_lattice& A,
		                          ring_matrix< bigint > & T, double y,
		                          lattice_info& li, sdigit factor)
{
	bigfloat_lattice tA(A);
	tA.lll_schnorr_euchner_orig(T, y, li, factor);
	return(tA);
}



bigfloat_lattice lll_schnorr_euchner_orig(const bigfloat_lattice& A,
		                          field_matrix< bigfloat > & T, double y,
		                          lattice_info& li, sdigit factor)
{
	bigfloat_lattice tA(A);
	tA.lll_schnorr_euchner_orig(T, y, li, factor);
	return(tA);
}



bigfloat_lattice lll_schnorr_euchner_orig(const bigfloat_lattice& A, double y,
		                          sdigit factor)
{
	bigfloat_lattice tA(A);
	tA.lll_schnorr_euchner_orig(y, factor);
	return(tA);
}



bigfloat_lattice lll_schnorr_euchner_orig(const bigfloat_lattice& A,
	                                  ring_matrix< bigint > & T, double y,
		                          sdigit factor)
{
	bigfloat_lattice tA(A);
	tA.lll_schnorr_euchner_orig(T, y, factor);
	return(tA);
}



bigfloat_lattice lll_schnorr_euchner_orig(const bigfloat_lattice& A,
					  field_matrix< bigfloat > & T, double y,
		                          sdigit factor)
{
	bigfloat_lattice tA(A);
	tA.lll_schnorr_euchner_orig(T, y, factor);
	return(tA);
}



//
// user defined Scalar Product
//
bigfloat_lattice lll_schnorr_euchner_orig(const bigfloat_lattice& A, double y,
		                          lattice_info& li, user_SP SP,
					  sdigit factor)
{
	bigfloat_lattice tA(A);
	tA.lll_schnorr_euchner_orig(y, li, SP, factor);
	return(tA);
}



bigfloat_lattice lll_schnorr_euchner_orig(const bigfloat_lattice& A,
		                          ring_matrix< bigint > & T, double y,
		                          lattice_info& li, user_SP SP,
					  sdigit factor)
{
	bigfloat_lattice tA(A);
	tA.lll_schnorr_euchner_orig(T, y, li, SP, factor);
	return(tA);
}



bigfloat_lattice lll_schnorr_euchner_orig(const bigfloat_lattice& A,
		                          field_matrix< bigfloat > & T, double y,
		                          lattice_info& li, user_SP SP,
					  sdigit factor)
{
	bigfloat_lattice tA(A);
	tA.lll_schnorr_euchner_orig(T, y, li, SP, factor);
	return(tA);
}



bigfloat_lattice lll_schnorr_euchner_orig(const bigfloat_lattice& A, double y,
		                          user_SP SP, sdigit factor)
{
	bigfloat_lattice tA(A);
	tA.lll_schnorr_euchner_orig(y, SP, factor);
	return(tA);
}



bigfloat_lattice lll_schnorr_euchner_orig(const bigfloat_lattice& A,
	                                  ring_matrix< bigint > & T, double y,
		                          user_SP SP, sdigit factor)
{
	bigfloat_lattice tA(A);
	tA.lll_schnorr_euchner_orig(T, y, SP, factor);
	return(tA);
}



bigfloat_lattice lll_schnorr_euchner_orig(const bigfloat_lattice& A,
	                                  field_matrix< bigfloat > & T, double y,
		                          user_SP SP, sdigit factor)
{
	bigfloat_lattice tA(A);
	tA.lll_schnorr_euchner_orig(T, y, SP, factor);
	return(tA);
}



//
// Variation of Schnorr - Euchner
//
// friend functions
//
bigfloat_lattice lll_schnorr_euchner_fact(const bigfloat_lattice& A, double y,
		                          lattice_info& li, sdigit factor)
{
	bigfloat_lattice tA(A);
	tA.lll_schnorr_euchner_fact(y, li, factor);
	return(tA);
}



bigfloat_lattice lll_schnorr_euchner_fact(const bigfloat_lattice& A,
		                          ring_matrix< bigint > & T, double y,
		                          lattice_info& li, sdigit factor)
{
	bigfloat_lattice tA(A);
	tA.lll_schnorr_euchner_fact(T, y, li, factor);
	return(tA);
}



bigfloat_lattice lll_schnorr_euchner_fact(const bigfloat_lattice& A,
		                          field_matrix< bigfloat > & T, double y,
		                          lattice_info& li, sdigit factor)
{
	bigfloat_lattice tA(A);
	tA.lll_schnorr_euchner_fact(T, y, li, factor);
	return(tA);
}



bigfloat_lattice lll_schnorr_euchner_fact(const bigfloat_lattice& A, double y,
		                          sdigit factor)
{
	bigfloat_lattice tA(A);
	tA.lll_schnorr_euchner_fact(y, factor);
	return(tA);
}



bigfloat_lattice lll_schnorr_euchner_fact(const bigfloat_lattice& A,
		                          ring_matrix< bigint > & T, double y,
		                          sdigit factor)
{
	bigfloat_lattice tA(A);
	tA.lll_schnorr_euchner_fact(T, y, factor);
	return(tA);
}



bigfloat_lattice lll_schnorr_euchner_fact(const bigfloat_lattice& A,
		                          field_matrix< bigfloat > & T, double y,
		                          sdigit factor)
{
	bigfloat_lattice tA(A);
	tA.lll_schnorr_euchner_fact(T, y, factor);
	return(tA);
}



//
// user defined Scalar Product
//
bigfloat_lattice lll_schnorr_euchner_fact(const bigfloat_lattice& A, double y,
		                          lattice_info& li, user_SP SP,
					  sdigit factor)
{
	bigfloat_lattice tA(A);
	tA.lll_schnorr_euchner_fact(y, li, SP, factor);
	return(tA);
}



bigfloat_lattice lll_schnorr_euchner_fact(const bigfloat_lattice& A,
		                          ring_matrix< bigint > & T, double y,
		                          lattice_info& li, user_SP SP,
					  sdigit factor)
{
	bigfloat_lattice tA(A);
	tA.lll_schnorr_euchner_fact(T, y, li, SP, factor);
	return(tA);
}



bigfloat_lattice lll_schnorr_euchner_fact(const bigfloat_lattice& A,
		                          field_matrix< bigfloat > & T, double y,
		                          lattice_info& li, user_SP SP,
					  sdigit factor)
{
	bigfloat_lattice tA(A);
	tA.lll_schnorr_euchner_fact(T, y, li, SP, factor);
	return(tA);
}



bigfloat_lattice lll_schnorr_euchner_fact(const bigfloat_lattice& A, double y,
		                          user_SP SP, sdigit factor)
{
	bigfloat_lattice tA(A);
	tA.lll_schnorr_euchner_fact(y, SP, factor);
	return(tA);
}



bigfloat_lattice lll_schnorr_euchner_fact(const bigfloat_lattice& A,
	                                  ring_matrix< bigint > & T, double y,
		                          user_SP SP, sdigit factor)
{
	bigfloat_lattice tA(A);
	tA.lll_schnorr_euchner_fact(T, y, SP, factor);
	return(tA);
}



bigfloat_lattice lll_schnorr_euchner_fact(const bigfloat_lattice& A,
	                                  field_matrix< bigfloat > & T, double y,
		                          user_SP SP, sdigit factor)
{
	bigfloat_lattice tA(A);
	tA.lll_schnorr_euchner_fact(T, y, SP, factor);
	return(tA);
}



//
// Buchmann - Kessler
//
bigfloat_lattice buchmann_kessler(const bigfloat_lattice& A,
			          ring_matrix< bigint > & T,
			          double y, lattice_info& li)
{
	bigfloat_lattice tA(A);
	tA.buchmann_kessler(T, y, li);
	return(tA);
}



bigfloat_lattice buchmann_kessler(const bigfloat_lattice& A,
			          field_matrix< bigfloat > & T,
			          double y, lattice_info& li)
{
	bigfloat_lattice tA(A);
	tA.buchmann_kessler(T, y, li);
	return(tA);
}



bigfloat_lattice buchmann_kessler(const bigfloat_lattice& A,
			          ring_matrix< bigint > & T, double y)
{
	bigfloat_lattice tA(A);
	tA.buchmann_kessler(T, y);
	return(tA);
}



bigfloat_lattice buchmann_kessler(const bigfloat_lattice& A,
			          field_matrix< bigfloat > & T, double y)
{
	bigfloat_lattice tA(A);
	tA.buchmann_kessler(T, y);
	return(tA);
}



//
// Modified lll
//
// friend functions
//
bigfloat_lattice mlll(const bigfloat_lattice& A, double y,
		      bigint*& v, lattice_info& li)
{
	bigfloat_lattice tA(A);
	tA.mlll(y, v, li);
	return(tA);
}



bigfloat_lattice mlll(const bigfloat_lattice& A, double y,
		      bigfloat*& v, lattice_info& li)
{
	bigfloat_lattice tA(A);
	tA.mlll(y, v, li);
	return(tA);
}



bigfloat_lattice mlll(const bigfloat_lattice& A, double y,
		      base_vector< bigint > & v, lattice_info& li)
{
	bigfloat_lattice tA(A);
	tA.mlll(y, v, li);
	return(tA);
}



bigfloat_lattice mlll(const bigfloat_lattice& A, double y,
		      base_vector< bigfloat > & v, lattice_info& li)
{
	bigfloat_lattice tA(A);
	tA.mlll(y, v, li);
	return(tA);
}



bigfloat_lattice mlll(const bigfloat_lattice& A, double y, bigint*& v)
{
	bigfloat_lattice tA(A);
	tA.mlll(y, v);
	return(tA);
}



bigfloat_lattice mlll(const bigfloat_lattice& A, double y, bigfloat*& v)
{
	bigfloat_lattice tA(A);
	tA.mlll(y, v);
	return(tA);
}



bigfloat_lattice mlll(const bigfloat_lattice& A, double y,
		      base_vector< bigint > & v)
{
	bigfloat_lattice tA(A);
	tA.mlll(y, v);
	return(tA);
}



bigfloat_lattice mlll(const bigfloat_lattice& A, double y,
		      base_vector< bigfloat > & v)
{
	bigfloat_lattice tA(A);
	tA.mlll(y, v);
	return(tA);
}



//
// Other Things
//
// flag - checkings
//
bool bigfloat_lattice::check_basis()
{
	debug_handler("bigfloat_lattice", "check_basis()");
	bigfloat_lattice tBfL(*this);
	dense_alg< bigint > da;
	bool flag;

	da.b.alg_trans = false;
	tBfL.Tr_dense_create_f(da);
	flag = CallBiL.TrD_check_basis(da);
	tBfL.Tr_dense_extract_f(da);
	if (flag)
		set_basis_flag();
	return(flag);
}



bool bigfloat_lattice::check_gram()
{
	debug_handler("bigfloat_lattice", "check_gram()");
	bigfloat_lattice tBfL(*this);
	dense_alg< bigint > da;
	bool flag;

	da.b.alg_trans = false;
	tBfL.Tr_dense_create_f(da);
	flag = CallBiL.TrD_check_gram(da);
	tBfL.Tr_dense_extract_f(da);
	if (flag)
		set_gram_flag();
	return(flag);
}



//
// lll - checkings
//
bool bigfloat_lattice::lll_check(sdigit nom, sdigit denom)
{
	debug_handler("bigfloat_lattice", "lll_check(nom, denom)");
	bigfloat_lattice tBfl(*this);
	dense_alg< bigint > da;
	bool red_flag;

	if (chk_gram())
		lidia_error_handler("bigint_lattice", "lll_check(nom, denom) "
				    ":: not avaidable for gram");

	da.b.y = static_cast<double>(nom)/static_cast<double>(denom);
	da.b.alg_trans = false;
	if ((da.b.y > 1.0) || (da.b.y < 0.5)) {
		lidia_warning_handler("bigfloat_lattice", "lll_check(nom, denom) :: "
				      "no allowed y for schnorr - euchner - lll");
		return(false);
	}
	tBfl.Tr_dense_create_f(da);
	red_flag = CallBiL.TrD_lll_check(da, nom, denom);
	tBfl.Tr_dense_extract_f(da);
	return(red_flag);
}



void bigfloat_lattice::lll_check_search(sdigit& nom, sdigit& denom)
{
	debug_handler("bigfloat_lattice", "lll_check_search(nom, denom)");
	bigfloat_lattice tBfl(*this);
	dense_alg< bigint > da;

	if (chk_gram())
		lidia_error_handler("bigfloat_lattice", "lll_check_search(nom, denom) "
				    ":: not avaidable for gram");

	da.b.alg_trans = false;
	tBfl.Tr_dense_create_f(da);
	CallBiL.TrD_lll_check_search(da, nom, denom);
	tBfl.Tr_dense_extract_f(da);
}



//
// Creating needed structure for dense lattice algorithms
//
// Reducing 2 matrix assignments in performing transposition
// 2*rows*columns real bigfloat (> 4 Bytes) assignments to
// 2*rows*columns*3 pointer assignments (4 Bytes)
//
void bigfloat_lattice::Tr_dense_create(dense_alg< bigfloat > & da)
{
	debug_handler("bigfloat_lattice", "Tr_dense_create(da)");

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
		// Transpose and swap matrix, structure of bigfloat_lattice will be
		// destroyed
		//
		da.b.transpose = true;
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
		da.s.value = new bigfloat*[da.b.rows];
		memory_handler(da.s.value, "bigfloat_lattice", "Tr_dense_create(da) :: "
			       "not enough memory !");
		da.s.delvalue = new bigfloat[da.b.rows*da.b.real_columns];
		memory_handler(da.s.delvalue, "bigfloat_lattice", "Tr_dense_create(da) ::"
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
void bigfloat_lattice::Tr_dense_extract(const dense_alg< bigfloat > & da)
{
	debug_handler("bigfloat_lattice", "Tr_dense_create(da)");
	bigfloat **Taddr;

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



void bigfloat_lattice::Tr_dense_create_f(dense_alg< bigint > & da)
{
	debug_handler("bigfloat_lattice", "Tr_dense_create_f(da)");
	sdigit min_exp;

	//
	// Search for the exponent
	//
	if (value[0][0].mantissa().is_zero())
		min_exp = 0;
	else
		min_exp = value[0][0].exponent();
	for (lidia_size_t i = 0; i < rows; i++)
		for (lidia_size_t j = 0; j < columns; j++)
			if (!(value[i][j].mantissa().is_zero()))
				min_exp = (value[i][j].exponent() < min_exp)?
					value[i][j].exponent():min_exp;

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
		//
		// Allocate memory for bigint lattice
		//
		da.s.value = new bigint*[da.b.rows];
		memory_handler(da.s.value, "bigfloat_lattice", "Tr_dense_create_f(da) :: "
			       "not enough memory !");
		da.s.delvalue = new bigint[da.b.rows*da.b.real_columns];
		memory_handler(da.s.delvalue, "bigfloat_lattice", "Tr_dense_create_f(da) ::"
			       "not enough memory !");

		for (lidia_size_t i = 0; i < da.b.rows; i++)
			da.s.value[i] = &da.s.delvalue[i*da.b.real_columns];

		//
		// Shift the matrix and store factor
		//
		if (min_exp < 0) {
			da.d.bit_factor = -min_exp;
			for (lidia_size_t i = 0; i < da.b.rows; i++)
				for (lidia_size_t j = 0; j < da.b.columns; j++) {
					da.s.value[i][j].assign(value[i][j].mantissa());
					if (!(da.s.value[i][j].is_zero()))
						shift_left(da.s.value[i][j], da.s.value[i][j],
							   (value[i][j].exponent()-min_exp));
				}
		}
		else {
			da.d.bit_factor = 0;
			for (lidia_size_t i = 0; i < rows; i++)
				for (lidia_size_t j = 0; j < columns; j++) {
					da.s.value[i][j].assign(value[i][j].mantissa());
					if (!(da.s.value[i][j].is_zero()))
						shift_left(da.s.value[i][j], da.s.value[i][j],
							   value[i][j].exponent());
				}
		}

		//
		// If trans - version, concate identic matrix
		//
		for (lidia_size_t i = da.b.columns; i < da.b.real_columns; i++)
			da.s.value[i-da.b.columns][i].assign_one();
	}
	else {
		//
		// Transpose and swap matrix, structure of bigfloat_lattice will be
		// destroyed
		//
		da.b.transpose = true;
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
		memory_handler(da.s.value, "bigfloat_lattice", "Tr_dense_create_f(da) :: "
			       "not enough memory !");
		da.s.delvalue = new bigint[da.b.rows*da.b.real_columns];
		memory_handler(da.s.delvalue, "bigfloat_lattice", "Tr_dense_create_f(da) ::"
			       "not enough memory !");

		for (lidia_size_t i = 0; i < da.b.rows; i++)
			da.s.value[i] = &da.s.delvalue[i*da.b.real_columns];

		//
		// Shift the matrix and store factor
		//
		if (min_exp < 0) {
			da.d.bit_factor = -min_exp;
			for (lidia_size_t i = 0; i < da.b.rows; i++)
				for (lidia_size_t j = 0; j < da.b.columns; j++) {
					da.s.value[i][j].assign(value[j][i].mantissa());
					if (!(da.s.value[i][j].is_zero()))
						shift_left(da.s.value[i][j], da.s.value[i][j],
							   (value[j][i].exponent()-min_exp));
				}
		}
		else {
			da.d.bit_factor = 0;
			for (lidia_size_t i = 0; i < rows; i++)
				for (lidia_size_t j = 0; j < columns; j++) {
					da.s.value[i][j].assign(value[j][i].mantissa());
					if (!(da.s.value[i][j].is_zero()))
						shift_left(da.s.value[i][j], da.s.value[i][j],
							   value[j][i].exponent());
				}
		}


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
void bigfloat_lattice::Tr_dense_extract_f(const dense_alg< bigint > & da)
{
	debug_handler("bigfloat_lattice", "Tr_dense_create_f(da)");
	bigint **Taddr;

	// Hier soll spaeter die Abfrage nach dem entspr. Bit hin !!
	if (!da.b.transpose) {
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
		for (lidia_size_t i = 0; i < da.b.rows; i++)
			for (lidia_size_t j = 0; j < da.b.columns; j++) {
				value[i][j].assign(da.s.value[i][j]);
				shift_right(value[i][j], value[i][j], da.d.bit_factor);
			}
	}
	else {
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
		resize(da.b.columns, da.b.rows);
		for (lidia_size_t i = 0; i < da.b.rows; i++)
			for (lidia_size_t j = 0; j < da.b.columns; j++) {
				value[j][i].assign(da.s.value[i][j]);
				shift_right(value[j][i], value[j][i], da.d.bit_factor);
			}
	}
	//
	// Free storage allocated by Tr_dense_create_f !!
	// Every Tr_dense_create_f is followed by a Tr_dense_extract_f
	//
	delete[] da.s.delvalue;
	delete[] da.s.value;
}



//
// Conversions
//
void bigfloat_lattice::conv_bfl_bin(field_matrix< bigfloat > & Mbfl,
			            ring_matrix< bigint > & Mbin)
{
	debug_handler("bigfloat_lattice", "conv_bfl_bin(Mbfl, Mbin)");
	bigint **MbinAddr;
	bigfloat **MbflAddr;
	lidia_size_t Mbflrows;
	lidia_size_t Mbflcolumns;

	Mbin.resize(Mbflrows = Mbfl.get_no_of_rows(),
		    Mbflcolumns = Mbfl.get_no_of_columns());
	MbinAddr = Mbin.get_data_address();
	MbflAddr = Mbfl.get_data_address();
	for (lidia_size_t i = 0; i < Mbflrows; i++)
		for (lidia_size_t j = 0; j < Mbflcolumns; j++)
			MbflAddr[i][j].bigintify(MbinAddr[i][j]);
}



void bigfloat_lattice::conv_bin_bfl(ring_matrix< bigint > & Mbin,
			            field_matrix< bigfloat > & Mbfl)
{
	debug_handler("bigfloat_lattice", "conv_bin_bfl(Mbin, Mbfl)");
	bigint **MbinAddr;
	bigfloat **MbflAddr;
	lidia_size_t Mbinrows;
	lidia_size_t Mbincolumns;

	Mbfl.resize(Mbinrows = Mbin.get_no_of_rows(),
		    Mbincolumns = Mbin.get_no_of_columns());
	MbinAddr = Mbin.get_data_address();
	MbflAddr = Mbfl.get_data_address();
	for (lidia_size_t i = 0; i < Mbinrows; i++)
		for (lidia_size_t j = 0; j < Mbincolumns; j++)
			MbflAddr[i][j].assign(MbinAddr[i][j]);
}



//
// Vector positions
//
void bigfloat_lattice::randomize_vectors()
{
	debug_handler("bigfloat_lattice", "randomize_vectors()");
	//
	// First step is to transpose the lattice
	// Second is to generate a permutation of the rows (using random - functions)
	// Third to transpose again
	//
	dense_alg< bigfloat > da;
	da.b.alg_trans = false;
	Tr_dense_create(da);
	TrD_randomize_vectors(da);
	Tr_dense_extract(da);
}



void bigfloat_lattice::sort_big_vectors()
{
	debug_handler("bigfloat_lattice", "sort_big_vectors()");
	//
	// First step is to transpose the lattice
	// Second is to sort the rows by length (biggest first)
	// Third to transpose again
	//
	dense_alg< bigfloat > da;
	da.b.alg_trans = false;
	Tr_dense_create(da);
	TrD_sort_vectors(da, 1);
	Tr_dense_extract(da);
}



void bigfloat_lattice::sort_small_vectors()
{
	debug_handler("bigfloat_lattice", "sort_small_vectors()");
	//
	// First step is to transpose the lattice
	// Second is to sort the rows by length (smallest first)
	// Third to transpose again
	//
	dense_alg< bigfloat > da;
	da.b.alg_trans = false;
	Tr_dense_create(da);
	TrD_sort_vectors(da, -1);
	Tr_dense_extract(da);
}



void bigfloat_lattice::sort_vectors(bfl_cmp_func cpf)
{
	debug_handler("bigfloat_lattice", "sort_vectors(cpf)");
	//
	// First step is to transpose the lattice
	// Second is to sort the rows by the given compare functions cpf
	// (See header file for more information)
	// Third to transpose again
	//
	dense_alg< bigfloat > da;
	da.b.alg_trans = false;
	Tr_dense_create(da);
	TrD_sort_vectors(da, cpf);
	Tr_dense_extract(da);
}



//
// Algorithm`s
//
//
// Real implementation of the algorithms
// They are working on the transposed lattice
//
void bigfloat_lattice::TrD_randomize_vectors(dense_alg< bigfloat > & da)
{
	debug_handler("bigfloat_lattice", "TrD_randomize_vectors(da)");
	random_generator rg;
	char *bitvect;
	sdigit *perm;
	sdigit ran;
	bigfloat **temp;


	//
	// Allocate memory
	//
	bitvect = new char[rows];
	memory_handler(bitvect, "bigfloat_lattice", "TrD_randomize_vectors(da) :: "
		       "not enough memory !");
	perm = new sdigit[rows];
	memory_handler(perm, "bigfloat_lattice", "TrD_randomize_vectors(da) :: "
		       "not enough memory !");
	temp = new bigfloat*[rows];
	memory_handler(temp, "bigfloat_lattice", "TrD_randomize_vectors(da) :: "
		       "not enough memory !");

	//
	// Clear bitvector
	//
	for (lidia_size_t i = 0; i < rows; bitvect[i++] = 0);
	//
	// Generate a permutation
	// Try rows times to find valid value
	//
	for (lidia_size_t i = 0; i < rows; i++) {
		for (lidia_size_t j = 0; j < rows; j++) {
			rg >> ran;
			ran = static_cast<sdigit>(ran%rows);
			if (bitvect[ran] == 0) {
				bitvect[ran] = 1;
				perm[i] = ran;
				break;
			}
			else
				if (j == rows-1) {
					for (lidia_size_t l = 0; l < rows; l++)
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
	for (lidia_size_t i = 0; i < rows; i++)
		temp[perm[i]] = da.s.value[i];
	for (lidia_size_t i = 0; i < rows; i++)
		da.s.value[i] = temp[i];

	//
	// Free allocated storage
	//
	delete[] perm;
	delete[] bitvect;
	delete[] temp;
}



void bigfloat_lattice::TrD_sort_vectors(dense_alg< bigfloat > & da, sdigit vgl)
{
	debug_handler("bigfloat_lattice", "TrD_sort_vectors(da, vgl)");
	p_vector< bigfloat > vector;
	bigfloat *quads;

	//
	// Allocate memory for scalar product of the vectors
	//
	quads = new bigfloat[rows];
	memory_handler(quads, "bigfloat_lattice", "TrD_sort_vectors(da, vgl) :: "
		       "not enough memory !");

	//
	// Compute scalar products
	//
	vector.vectsize = da.b.columns;
	for (lidia_size_t i = 0; i < rows; i++)
		vector.scalprod(quads[i], da.s.value[i], da.s.value[i]);
	//
	// sort by scalar products ("length^2")
	// vgl = -1   ->   smallest vector first
	// vgl =  1   ->   biggest vector first
	//
	for (lidia_size_t i = 0; i < rows-1; i++)
		for (lidia_size_t j = i+1; j < rows; j++) {
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



void bigfloat_lattice::TrD_sort_vectors(dense_alg< bigfloat > & da,
					bfl_cmp_func cpf)
{
	debug_handler("bigfloat_lattice", "TrD_sort_vectors(da, cpf)");
	p_vector< bigfloat > vector;
	//
	// Sort vectors by specified compare function cpf
	// Perform bubble sort (for typical sizes of the
	// lattice quicksort not nessesary)
	//
	vector.vectsize = da.b.columns;
	for (lidia_size_t i = 0; i < rows-1; i++)
		for (lidia_size_t j = i+1; j < rows; j++)
			if (cpf(da.s.value[i], da.s.value[j], da.b.columns) > 0)
				vector.swap(da.s.value[i], da.s.value[j]);
}



sdigit bigfloat_lattice::compute_read_precision(const dense_alg< bigfloat > & da)
{
	debug_handler("bigfloat_lattice", "compute_read_precision(da)");
	//
	// Give the digit size of longest value
	//
	sdigit new_prec;
	sdigit read_prec = 0;
	for (lidia_size_t x = 0; x < rows; ++x)
		for (lidia_size_t y = 0; y < columns; ++y)
			if (read_prec < (new_prec = static_cast<sdigit>(da.s.value[x][y].bit_length()/
									std::log(10.0)+1)))
				read_prec = new_prec;
	return (read_prec);
}



sdigit bigfloat_lattice::compute_precision(const dense_alg< bigfloat > & da)
{
	debug_handler("bigfloat_lattice", "compute_precision(da)");
	bigfloat alpha, zweipotq;
	sdigit read_prec = compute_read_precision(da);
	sdigit n2 = rows;
	sdigit prec;
	alpha_compute(da, alpha);
	zwei_pot_q_compute(da, zweipotq, n2, alpha);
	read_prec = 2*read_prec+rows-1;
	prec = static_cast<sdigit>(2*(((zweipotq.bit_length()+1)*
				       std::log(2.0)/std::log(10.0)+1)+read_prec)+(columns+rows)-1);
	return (prec);
}



void bigfloat_lattice::gamma_compute(bigfloat& g, sdigit l)
{
	debug_handler("bigfloat_lattice", "gamma_compute(g, l)");
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



void bigfloat_lattice::alpha_compute(const dense_alg< bigfloat > & da,
				     bigfloat& alpha)
{
	//
	// alpha = max{vect[1].l2_norm(),...,vect[columns].l2_norm()}
	// norm = vect[i].l2_norm(), 0 <= i < columns
	//
	debug_handler("bigfloat_lattice", "alpha_compute(da, alpha)");
	bigfloat tempbfl0;
	p_vector< bigfloat > vector;
	vector.vectsize = da.b.columns;
	vector.scalprod(alpha, da.s.value[0], da.s.value[0]);
	for (lidia_size_t i = 1; i < rows; ++i) {
		vector.scalprod(tempbfl0, da.s.value[i], da.s.value[i]);
		if (alpha.compare(tempbfl0) < 0)
			alpha.assign(tempbfl0);
	}
	//
	// calculating sqrt of l2_norm
	//
	LiDIA::sqrt(alpha, alpha);
}



void bigfloat_lattice::zwei_pot_q_compute(const dense_alg< bigfloat > & da,
        				  bigfloat& zweipotq,
					  sdigit& n2, bigfloat& alpha)
{
	debug_handler("bigfloat_lattice", "zwei_pot_q_compute(da, zwpq, n2, alpha)");
	sdigit beta = da.b.rows;
	bigfloat tempbin0;
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
