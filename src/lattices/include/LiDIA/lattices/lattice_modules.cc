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
//	Changes	: see CVS log
//
//==============================================================================================


#ifndef LIDIA_LATTICE_MODULES_CC_GUARD_
#define LIDIA_LATTICE_MODULES_CC_GUARD_



/*
 * g++-2.95.2 occasionally does not instantiate inline specializations
 * of template class member functions when compiling with
 * -fno-implicit-templates (though that should not affect inline
 * definitions).  To work around this, we control the instantiations
 * ourselves, using the pragmas `interface' and `implementation',
 * which have the additional advantage of avoiding the duplication of
 * instantiations.
 */
#if __GNUC__ && __GNUC__ < 3
#pragma interface "lattice_modules.cc"
#endif



#ifdef LIDIA_NAMESPACE
# ifndef IN_NAMESPACE_LIDIA
namespace LiDIA {
# endif
#endif



//
// Code for
// 	     * base_modules < {bigfloat, bigint}, {double, xdouble, bigfloat} >
//
//
// bigint / double
//
template <>
inline void
prec_modules< bigint, double >::prec_startup (dense_alg< bigint > &)
{
}



template <>
inline void
prec_modules< bigint, double >::prec_update (dense_alg< bigint > & da)
{
	da.d.cut_bit_prec = (1+(da.d.bit_prec/MANTISSA_CUT))*MANTISSA_CUT;
	da.d.approx_prec = static_cast<sdigit>(static_cast<double>(da.d.bit_prec)*std::log(2.0)/std::log(10.0));
}



template <>
inline void
prec_modules< bigint, double >::prec_exact (const dense_alg< bigint > &)
{
}



template <>
inline void
prec_modules< bigint, double >::prec_approx (const dense_alg< bigint > &)
{
}



template <>
inline void
prec_modules< bigint, double >::prec_correct (const dense_alg< bigint > &)
{
}



template <>
inline void
prec_modules< bigint, double >::prec_restore (const dense_alg< bigint > &)
{
}


//
// bigint / xdouble
//
template <>
inline void
prec_modules< bigint, xdouble >::prec_startup (dense_alg< bigint > &)
{
}



template <>
inline void
prec_modules< bigint, xdouble >::prec_update (dense_alg< bigint > & da)
{
	da.d.cut_bit_prec = (1+(da.d.bit_prec/MANTISSA_CUT))*MANTISSA_CUT;
	da.d.approx_prec = static_cast<sdigit>(static_cast<double>(da.d.bit_prec)*std::log(2.0)/std::log(10.0));
}



template <>
inline void
prec_modules< bigint, xdouble >::prec_exact (const dense_alg< bigint > &)
{
}



template <>
inline void
prec_modules< bigint, xdouble >::prec_approx (const dense_alg< bigint > &)
{
}



template <>
inline void
prec_modules< bigint, xdouble >::prec_correct (const dense_alg< bigint > &)
{
}



template <>
inline void
prec_modules< bigint, xdouble >::prec_restore (const dense_alg< bigint > &)
{
}



//
// bigint / bigfloat
//
template <>
inline void
prec_modules< bigint, bigfloat >::prec_startup (dense_alg< bigint > & da)
{
	da.d.old_prec = bigfloat::get_precision();
	da.d.approx_prec = static_cast<sdigit>(static_cast<double>(da.d.bit_prec)*std::log(2.0)/std::log(10.0));
}



template <>
inline void
prec_modules< bigint, bigfloat >::prec_update (dense_alg< bigint > & da)
{
	da.d.cut_bit_prec = (1+(da.d.bit_prec/MANTISSA_CUT))*MANTISSA_CUT;
	da.d.approx_prec = static_cast<sdigit>(static_cast<double>(da.d.bit_prec)*std::log(2.0)/std::log(10.0));
}



template <>
inline void
prec_modules< bigint, bigfloat >::prec_exact (const dense_alg< bigint > &)
{
}



template <>
inline void
prec_modules< bigint, bigfloat >::prec_approx (const dense_alg< bigint > & da)
{
	bigfloat::set_precision(da.d.approx_prec);
}



template <>
inline void
prec_modules< bigint, bigfloat >::prec_correct (const dense_alg< bigint > &)
{
}



template <>
inline void
prec_modules< bigint, bigfloat >::prec_restore (const dense_alg< bigint > & da)
{
	bigfloat::set_precision(da.d.old_prec);
}



//
// bigfloat / double
//
template <>
inline void
prec_modules< bigfloat, double >::prec_startup (dense_alg< bigfloat > & da)
{
	//
	// Compute precision of read lattice
	//
	da.d.old_prec = bigfloat::get_precision();
	da.d.read_prec = 0;
	for (lidia_size_t i = 0; i < da.b.rows; ++i)
		for (lidia_size_t j = 0; j < da.b.columns; ++j)
			if (da.d.read_prec < (tempsdt = static_cast<sdigit>(da.s.value[i][j].bit_length()/
									    std::log(10.0)+1)))
				da.d.read_prec = tempsdt;
	da.d.exact_prec = da.d.read_prec*((da.b.columns > da.b.rows)?
					  da.b.columns:da.b.rows);
	da.d.exact_prec *= Magic_Precision_Factor;
}



template <>
inline void
prec_modules< bigfloat, double >::prec_update (dense_alg< bigfloat > & da)
{
	da.d.cut_bit_prec = (1+(da.d.bit_prec/MANTISSA_CUT))*MANTISSA_CUT;
	da.d.approx_prec = static_cast<sdigit>(static_cast<double>(da.d.bit_prec)*std::log(2.0)/std::log(10.0));
}



template <>
inline void
prec_modules< bigfloat, double >::prec_approx (const dense_alg< bigfloat > &)
{
}



template <>
inline void
prec_modules< bigfloat, double >::prec_exact (const dense_alg< bigfloat > & da)
{
	bigfloat::set_precision(da.d.exact_prec);
}



template <>
inline void
prec_modules< bigfloat, double >::prec_correct (const dense_alg< bigfloat > & da)
{
	bigfloat::set_precision(da.b.columns+(static_cast<sdigit>(DOUBLE_MANTISSA_BITS*std::log(2.0)/
							      std::log(10.0)+1)*2));
}



template <>
inline void
prec_modules< bigfloat, double >::prec_restore (const dense_alg< bigfloat > & da)
{
	bigfloat::set_precision(da.d.old_prec);
}



//
// bigfloat / xdouble
//
template <>
inline void
prec_modules< bigfloat, xdouble >::prec_startup (dense_alg< bigfloat > & da)
{
	da.d.old_prec = bigfloat::get_precision();
	//
	// Compute precision of read lattice
	//
	da.d.read_prec = 0;
	for (lidia_size_t i = 0; i < da.b.rows; ++i)
		for (lidia_size_t j = 0; j < da.b.columns; ++j)
			if (da.d.read_prec < (tempsdt = static_cast<sdigit>(da.s.value[i][j].bit_length()/
									    std::log(10.0)+1)))
				da.d.read_prec = tempsdt;
	da.d.exact_prec = da.d.read_prec*((da.b.columns > da.b.rows)?
					  da.b.columns:da.b.rows);
	da.d.exact_prec *= Magic_Precision_Factor;
}



template <>
inline void
prec_modules< bigfloat, xdouble >::prec_update (dense_alg< bigfloat > & da)
{
	da.d.cut_bit_prec = (1+(da.d.bit_prec/MANTISSA_CUT))*MANTISSA_CUT;
	da.d.approx_prec = static_cast<sdigit>(static_cast<double>(da.d.bit_prec)*std::log(2.0)/std::log(10.0));
}



template <>
inline void
prec_modules< bigfloat, xdouble >::prec_approx (const dense_alg< bigfloat > &)
{
}



template <>
inline void
prec_modules< bigfloat, xdouble >::prec_exact (const dense_alg< bigfloat > & da)
{
	bigfloat::set_precision(da.d.exact_prec);
}



template <>
inline void
prec_modules< bigfloat, xdouble >::prec_correct (const dense_alg< bigfloat > & da)
{
	bigfloat::set_precision(da.b.columns+(static_cast<sdigit>(DOUBLE_MANTISSA_BITS*std::log(2.0)/
							      std::log(10.0)+1)*4));
}



template <>
inline void
prec_modules< bigfloat, xdouble >::prec_restore (const dense_alg< bigfloat > & da)
{
	bigfloat::set_precision(da.d.old_prec);
}



//
// bigfloat / bigfloat
//
template <>
inline void
prec_modules< bigfloat, bigfloat >::prec_startup (dense_alg< bigfloat > & da)
{
	da.d.old_prec = bigfloat::get_precision();
	da.d.approx_prec = static_cast<sdigit>(static_cast<double>(da.d.bit_prec)*std::log(2.0)/std::log(10.0));
//
// Compute precision of read lattice
//
	da.d.read_prec = 0;
	for (lidia_size_t i = 0; i < da.b.rows; ++i)
		for (lidia_size_t j = 0; j < da.b.columns; ++j)
			if (da.d.read_prec < (tempsdt = static_cast<sdigit>(da.s.value[i][j].bit_length()/
									    std::log(10.0)+1)))
				da.d.read_prec = tempsdt;
	da.d.exact_prec = da.d.read_prec*((da.b.columns > da.b.rows)?
					  da.b.columns:da.b.rows);
	da.d.exact_prec *= Magic_Precision_Factor;
}



template <>
inline void
prec_modules< bigfloat, bigfloat >::prec_update (dense_alg< bigfloat > & da)
{
	da.d.cut_bit_prec = (1+(da.d.bit_prec/MANTISSA_CUT))*MANTISSA_CUT;
	da.d.approx_prec = static_cast<sdigit>(static_cast<double>(da.d.bit_prec)*std::log(2.0)/std::log(10.0));
}



template <>
inline void
prec_modules< bigfloat, bigfloat >::prec_approx (const dense_alg< bigfloat > & da)
{
	bigfloat::set_precision(da.d.approx_prec);
}



template <>
inline void
prec_modules< bigfloat, bigfloat >::prec_exact (const dense_alg< bigfloat > & da)
{
	bigfloat::set_precision(da.d.exact_prec);
}



template <>
inline void
prec_modules< bigfloat, bigfloat >::prec_correct (const dense_alg< bigfloat > & da)
{
	bigfloat::set_precision(da.b.columns+(static_cast<sdigit>(da.d.bit_prec*std::log(2.0)/
							      std::log(10.0)+1)*2));
}



template <>
inline void
prec_modules< bigfloat, bigfloat >::prec_restore (const dense_alg< bigfloat > & da)
{
	bigfloat::set_precision(da.d.old_prec);
}



//
// End of
// 	     * prec_modules < {bigfloat, bigint}, {double, xdouble, bigfloat} >
//

//
// Code for
// 	     * base_modules < bigint, {double, xdouble, bigfloat}, Normal >
// 	     * basis_modules < bigint, {double, xdouble, bigfloat}, Normal >
// 	     * gensys_modules < bigint, {double, xdouble, bigfloat}, Normal >
//
template <>
inline void
base_modules< bigint, double, Normal >::E_convert_value_A (dense_alg< bigint > &,
							   double& d,
							   const bigint& bin)
{
	d = dbl(bin);
}



template <>
inline bool
base_modules< bigint, double, Normal >::A_convert_value_E_bound (dense_alg< bigint > &,
								 bigint& bin,
								 const double& dbl,
								 const double& bound)
{
	if (fabs(dbl) > bound) {
		tempbfl.assign(dbl);
		tempbfl.bigintify(bin);
		return (true);
	}
	else {
		bin.assign(static_cast<sdigit>(dbl));
		return (false);
	}
}



template <>
inline void
base_modules< bigint, xdouble, Normal >::E_convert_value_A (dense_alg< bigint > &,
							    xdouble& xd,
							    const bigint& bin)
{
	xd = xdbl(bin);
}



template <>
inline bool
base_modules< bigint, xdouble, Normal >::A_convert_value_E_bound (dense_alg< bigint > &,
								  bigint& bin,
								  const xdouble& xdbl,
								  const xdouble& bound)
{
	tempbfl.assign(xdbl);
	tempbfl.bigintify(bin);
	if (fabs(xdbl) > bound)
		return (true);
	return (false);
}



template <>
inline void
base_modules< bigint, bigfloat, Normal >::E_convert_value_A (dense_alg< bigint > & da,
							     bigfloat& bfl,
							     const bigint& bin)
{
	bi_bit_len = bin.bit_length();
	if (bi_bit_len > da.d.cut_bit_prec) {
		bit_diff = bi_bit_len-da.d.cut_bit_prec;
		shift_right(tempbin, bin, bit_diff);
		bfl.assign(tempbin);
		shift_left(bfl, bfl, bit_diff);
	}
	else
		bfl.assign(bin);
}



template <>
inline bool
base_modules< bigint, bigfloat, Normal >::A_convert_value_E_bound (dense_alg< bigint > &,
								   bigint& bin,
								   const bigfloat& bfl,
								   const bigfloat& bound)
{
	bfl.bigintify(bin);
	if (abs(bfl).compare(bound) > 0)
		return (true);
	return (false);
}



template <>
inline void
basis_modules< bigint, double, Normal >::E_convert_lattice_A (dense_alg< bigint > & da,
							      double** dblvalue)
{
	for (lidia_size_t i = 0; i < da.b.rows; i++)
		for (lidia_size_t j = 0; j < da.b.columns; j++)
			dblvalue[i][j] = dbl(da.s.value[i][j]);
}



template <>
inline bool
basis_modules< bigint, double, Normal >::E_convert_vector_A (dense_alg< bigint > & da,
							     double **dblvalue,
							     const lidia_size_t k)
{
	for (lidia_size_t i = 0; i < da.b.columns; i++)
		dblvalue[k][i] = dbl(da.s.value[k][i]);
	return(false);
}



template <>
inline void
basis_modules< bigint, xdouble, Normal >::E_convert_lattice_A (dense_alg< bigint > & da,
							       xdouble** xdblvalue)
{
	for (lidia_size_t i = 0; i < da.b.rows; i++)
		for (lidia_size_t j = 0; j < da.b.columns; j++)
			xdblvalue[i][j] = xdbl(da.s.value[i][j]);
}



template <>
inline bool
basis_modules< bigint, xdouble, Normal >::E_convert_vector_A (dense_alg< bigint > & da,
							      xdouble **xdblvalue,
							      const lidia_size_t k)
{
	for (lidia_size_t i = 0; i < da.b.columns; i++)
		xdblvalue[k][i] = xdbl(da.s.value[k][i]);
	return(false);
}



template <>
inline void
basis_modules< bigint, bigfloat, Normal >::E_convert_lattice_A (dense_alg< bigint > & da,
								bigfloat** bflvalue)
{
	for (lidia_size_t i = 0; i < da.b.rows; i++)
		for (lidia_size_t j = 0; j < da.b.columns; j++) {
			bi_bit_len = da.s.value[i][j].bit_length();
			if (bi_bit_len > da.d.cut_bit_prec) {
				bit_diff = bi_bit_len-da.d.cut_bit_prec;
				shift_right(tempbin, da.s.value[i][j], bit_diff);
				bflvalue[i][j].assign(tempbin);
				shift_left(bflvalue[i][j], bflvalue[i][j], bit_diff);
			}
			else
				bflvalue[i][j].assign(da.s.value[i][j]);
		}
}



template <>
inline bool
basis_modules< bigint, bigfloat, Normal >::E_convert_vector_A (dense_alg< bigint > & da,
							       bigfloat** bflvalue,
							       const lidia_size_t k)
{
	for (lidia_size_t i = 0; i < da.b.columns; i++) {
		bi_bit_len = da.s.value[k][i].bit_length();
		if (bi_bit_len > da.d.cut_bit_prec) {
			bit_diff = bi_bit_len-da.d.cut_bit_prec;
			shift_right(tempbin, da.s.value[k][i], bit_diff);
			bflvalue[k][i].assign(tempbin);
			shift_left(bflvalue[k][i], bflvalue[k][i], bit_diff);
		}
		else
			bflvalue[k][i].assign(da.s.value[k][i]);
	}
	return (false);
}



template <>
inline void
gensys_modules< bigint, double, Normal >::E_convert_lattice_A (dense_alg< bigint > & da,
							       double** dblvalue)
{
	for (lidia_size_t i = 0; i < da.b.rank;) {
		is_zero = true;
		for (lidia_size_t j = 0; j < da.b.columns; j++)
			if ((!is_zero) || (!da.s.value[i][j].is_zero())) {
				dblvalue[i][j] = dbl(da.s.value[i][j]);
				is_zero = false;
			}
			else
				dblvalue[i][j] = 0.0;
		if (is_zero) {
			da.b.rank--;
			for (lidia_size_t j = i; j < da.b.rank; j++) {
				tempP = static_cast<void *>(dblvalue[j]);
				dblvalue[j] = dblvalue[j+1];
				dblvalue[j+1] = static_cast<double *>(tempP);
				tempP = static_cast<void *>(da.s.value[j]);
				da.s.value[j] = da.s.value[j+1];
				da.s.value[j+1] = static_cast<bigint *>(tempP);
			}
		}
		else
			i++;
	}
}



template <>
inline bool
gensys_modules< bigint, double, Normal >::E_convert_vector_A (dense_alg< bigint > & da,
							      double** dblvalue,
							      const lidia_size_t k)
{
	is_zero = true;
	for (lidia_size_t i = 0; i < da.b.columns; i++)
		if ((!is_zero) || (!da.s.value[k][i].is_zero())) {
			dblvalue[k][i] = dbl(da.s.value[k][i]);
			is_zero = false;
		}
		else
			dblvalue[k][i] = 0.0;
	if (is_zero) {
		da.b.rank--;
		for (lidia_size_t i = k; i < da.b.rank; i++) {
			tempP = static_cast<void *>(dblvalue[i]);
			dblvalue[i] = dblvalue[i+1];
			dblvalue[i+1] = static_cast<double *>(tempP);
			tempP = static_cast<void *>(da.s.value[i]);
			da.s.value[i] = da.s.value[i+1];
			da.s.value[i+1] = static_cast<bigint *>(tempP);
		}
		return(true);
	}
	return(false);
}



template <>
inline void
gensys_modules< bigint, xdouble, Normal >::E_convert_lattice_A (dense_alg< bigint > & da,
								xdouble** xdblvalue)
{
	for (lidia_size_t i = 0; i < da.b.rank;) {
		is_zero = true;
		for (lidia_size_t j = 0; j < da.b.columns; j++)
			if ((!is_zero) || (!da.s.value[i][j].is_zero())) {
				xdblvalue[i][j] = xdbl(da.s.value[i][j]);
				is_zero = false;
			}
			else
				xdblvalue[i][j] = 0.0;
		if (is_zero) {
			da.b.rank--;
			for (lidia_size_t j = i; j < da.b.rank; j++) {
				tempP = static_cast<void *>(xdblvalue[j]);
				xdblvalue[j] = xdblvalue[j+1];
				xdblvalue[j+1] = static_cast<xdouble *>(tempP);
				tempP = static_cast<void *>(da.s.value[j]);
				da.s.value[j] = da.s.value[j+1];
				da.s.value[j+1] = static_cast<bigint *>(tempP);
			}
		}
		else
			i++;
	}
}



template <>
inline bool
gensys_modules< bigint, xdouble, Normal >::E_convert_vector_A (dense_alg< bigint > & da,
							       xdouble** xdblvalue,
							       const lidia_size_t k)
{
	is_zero = true;
	for (lidia_size_t i = 0; i < da.b.columns; i++)
		if ((!is_zero) || (!da.s.value[k][i].is_zero())) {
			xdblvalue[k][i] = xdbl(da.s.value[k][i]);
			is_zero = false;
		}
		else
			xdblvalue[k][i] = 0.0;
	if (is_zero) {
		da.b.rank--;
		for (lidia_size_t i = k; i < da.b.rank; i++) {
			tempP = static_cast<void *>(xdblvalue[i]);
			xdblvalue[i] = xdblvalue[i+1];
			xdblvalue[i+1] = static_cast<xdouble *>(tempP);
			tempP = static_cast<void *>(da.s.value[i]);
			da.s.value[i] = da.s.value[i+1];
			da.s.value[i+1] = static_cast<bigint *>(tempP);
		}
		return(true);
	}
	return(false);
}



template <>
inline void
gensys_modules< bigint, bigfloat, Normal >::E_convert_lattice_A (dense_alg< bigint > & da,
								 bigfloat** bflvalue)
{
	for (lidia_size_t i = 0; i < da.b.rank;) {
		is_zero = true;
		for (lidia_size_t j = 0; j < da.b.columns; j++)
			if ((!is_zero) || (!da.s.value[i][j].is_zero())) {
				bi_bit_len = da.s.value[i][j].bit_length();
				if (bi_bit_len > da.d.cut_bit_prec) {
					bit_diff = bi_bit_len-da.d.cut_bit_prec;
					shift_right(tempbin, da.s.value[i][j], bit_diff);
					bflvalue[i][j].assign(tempbin);
					shift_left(bflvalue[i][j], bflvalue[i][j], bit_diff);
				}
				else
					bflvalue[i][j].assign(da.s.value[i][j]);
				is_zero = false;
			}
			else
				bflvalue[i][j].assign_zero();
		if (is_zero) {
			da.b.rank--;
			for (lidia_size_t j = i; j < da.b.rank; j++) {
				tempP = static_cast<void *>(bflvalue[j]);
				bflvalue[j] = bflvalue[j+1];
				bflvalue[j+1] = static_cast<bigfloat *>(tempP);
				tempP = static_cast<void *>(da.s.value[j]);
				da.s.value[j] = da.s.value[j+1];
				da.s.value[j+1] = static_cast<bigint *>(tempP);
			}
		}
		else
			i++;
	}
}



template <>
inline bool
gensys_modules< bigint, bigfloat, Normal >::E_convert_vector_A (dense_alg< bigint > & da,
								bigfloat** bflvalue,
								const lidia_size_t k)
{
	is_zero = true;
	for (lidia_size_t i = 0; i < da.b.columns; i++)
		if ((!is_zero) || (!da.s.value[k][i].is_zero())) {
			bi_bit_len = da.s.value[k][i].bit_length();
			if (bi_bit_len > da.d.cut_bit_prec) {
				bit_diff = bi_bit_len-da.d.cut_bit_prec;
				shift_right(tempbin, da.s.value[k][i], bit_diff);
				bflvalue[k][i].assign(tempbin);
				shift_left(bflvalue[k][i], bflvalue[k][i], bit_diff);
			}
			else
				bflvalue[k][i].assign(da.s.value[k][i]);
			is_zero = false;
		}
		else
			bflvalue[k][i].assign_zero();
	if (is_zero) {
		da.b.rank--;
		for (lidia_size_t i = k; i < da.b.rank; i++) {
			tempP = static_cast<void *>(bflvalue[i]);
			bflvalue[i] = bflvalue[i+1];
			bflvalue[i+1] = static_cast<bigfloat *>(tempP);
			tempP = static_cast<void *>(da.s.value[i]);
			da.s.value[i] = da.s.value[i+1];
			da.s.value[i+1] = static_cast<bigint *>(tempP);
		}
		return(true);
	}
	return(false);
}



//
// End of
// 	     * base_modules < bigint, {double, xdouble, bigfloat}, Normal >
// 	     * basis_modules < bigint, {double, xdouble, bigfloat}, Normal >
// 	     * gensys_modules < bigint, {double, xdouble, bigfloat}, Normal >
//

//
// Code for
// 	     * base_modules < bigint, {double, xdouble, bigfloat}, VariationI >
// 	     * basis_modules < bigint, {double, xdouble, bigfloat}, VariationI >
// 	     * gensys_modules < bigint, {double, xdouble, bigfloat}, VariationI >
//
template <>
inline void
base_modules< bigint, double, VariationI >::E_convert_value_A (dense_alg< bigint > & da,
							       double& d,
							       const bigint& bin)
{
	tempbfl.assign(bin);
	shift_right(tempbfl, tempbfl, da.d.bit_factor);
	tempbfl.doublify(d);
}



template <>
inline bool
base_modules< bigint, double, VariationI >::A_convert_value_E_bound (dense_alg< bigint > &,
								     bigint& bin,
								     const double& dbl,
								     const double& bound)
{
	if (fabs(dbl) > bound) {
		tempbfl.assign(dbl);
		tempbfl.bigintify(bin);
		return (true);
	}
	else {
		bin.assign(static_cast<sdigit>(dbl));
		return (false);
	}
}



template <>
inline void
base_modules< bigint, xdouble, VariationI >::E_convert_value_A (dense_alg< bigint > & da,
								xdouble& xd,
								const bigint& bin)
{
	tempbfl.assign(bin);
	shift_right(tempbfl, tempbfl, da.d.bit_factor);
	tempbfl.xdoublify(xd);
}



template <>
inline bool
base_modules< bigint, xdouble, VariationI >::A_convert_value_E_bound (dense_alg< bigint > &,
								      bigint& bin,
								      const xdouble& xdbl,
								      const xdouble& bound)
{
	tempbfl.assign(xdbl);
	tempbfl.bigintify(bin);
	if (fabs(xdbl) > bound)
		return (true);
	return (false);
}



template <>
inline void
base_modules< bigint, bigfloat, VariationI >::E_convert_value_A (dense_alg< bigint > & da,
								 bigfloat& bfl,
								 const bigint& bin)
{
	bi_bit_len = bin.bit_length();
	if (bi_bit_len > da.d.cut_bit_prec) {
		bit_diff = bi_bit_len-da.d.cut_bit_prec;
		shift_right(tempbin, bin, bit_diff);
		bfl.assign(tempbin);
		shift_left(bfl, bfl, bit_diff);
	}
	else
		bfl.assign(bin);
	shift_right(bfl, bfl, da.d.bit_factor);
}



template <>
inline bool
base_modules< bigint, bigfloat, VariationI >::A_convert_value_E_bound (dense_alg< bigint > &,
								       bigint& bin,
								       const bigfloat& bfl,
								       const bigfloat& bound)
{
	bfl.bigintify(bin);
	if (abs(bfl).compare(bound) > 0)
		return (true);
	return (false);
}



template <>
inline void
basis_modules< bigint, double, VariationI >::E_convert_lattice_A (dense_alg< bigint > & da,
								  double** dblvalue)
{
	for (lidia_size_t i = 0; i < da.b.rows; i++)
		for (lidia_size_t j = 0; j < da.b.columns; j++) {
			tempbfl.assign(da.s.value[i][j]);
			shift_right(tempbfl, tempbfl, da.d.bit_factor);
			tempbfl.doublify(dblvalue[i][j]);
		}
}



template <>
inline bool
basis_modules< bigint, double, VariationI >::E_convert_vector_A (dense_alg< bigint > & da,
								 double **dblvalue,
								 const lidia_size_t k)
{
	for (lidia_size_t i = 0; i < da.b.columns; i++) {
		tempbfl.assign(da.s.value[k][i]);
		shift_right(tempbfl, tempbfl, da.d.bit_factor);
		tempbfl.doublify(dblvalue[k][i]);
	}
	return(false);
}



template <>
inline void
basis_modules< bigint, xdouble, VariationI >::E_convert_lattice_A (dense_alg< bigint > & da,
								   xdouble** xdblvalue)
{
	for (lidia_size_t i = 0; i < da.b.rows; i++)
		for (lidia_size_t j = 0; j < da.b.columns; j++) {
			tempbfl.assign(da.s.value[i][j]);
			shift_right(tempbfl, tempbfl, da.d.bit_factor);
			tempbfl.xdoublify(xdblvalue[i][j]);
		}
}



template <>
inline bool
basis_modules< bigint, xdouble, VariationI >::E_convert_vector_A (dense_alg< bigint > & da,
								  xdouble **xdblvalue,
								  const lidia_size_t k)
{
	for (lidia_size_t i = 0; i < da.b.columns; i++) {
		tempbfl.assign(da.s.value[k][i]);
		shift_right(tempbfl, tempbfl, da.d.bit_factor);
		tempbfl.xdoublify(xdblvalue[k][i]);
	}
	return(false);
}



template <>
inline void
basis_modules< bigint, bigfloat, VariationI >::E_convert_lattice_A (dense_alg< bigint > & da,
								    bigfloat** bflvalue)
{
	for (lidia_size_t i = 0; i < da.b.rows; i++)
		for (lidia_size_t j = 0; j < da.b.columns; j++) {
			bi_bit_len = da.s.value[i][j].bit_length();
			if (bi_bit_len > da.d.cut_bit_prec) {
				bit_diff = bi_bit_len-da.d.cut_bit_prec;
				shift_right(tempbin, da.s.value[i][j], bit_diff);
				bflvalue[i][j].assign(tempbin);
				shift_left(bflvalue[i][j], bflvalue[i][j], bit_diff);
			}
			else
				bflvalue[i][j].assign(da.s.value[i][j]);
			shift_right(bflvalue[i][j], bflvalue[i][j], da.d.bit_factor);
		}
}



template <>
inline bool
basis_modules< bigint, bigfloat, VariationI >::E_convert_vector_A (dense_alg< bigint > & da,
								   bigfloat** bflvalue,
								   const lidia_size_t k)
{
	for (lidia_size_t i = 0; i < da.b.columns; i++) {
		bi_bit_len = da.s.value[k][i].bit_length();
		if (bi_bit_len > da.d.cut_bit_prec) {
			bit_diff = bi_bit_len-da.d.cut_bit_prec;
			shift_right(tempbin, da.s.value[k][i], bit_diff);
			bflvalue[k][i].assign(tempbin);
			shift_left(bflvalue[k][i], bflvalue[k][i], bit_diff);
		}
		else
			bflvalue[k][i].assign(da.s.value[k][i]);
		shift_right(bflvalue[k][i], bflvalue[k][i], da.d.bit_factor);
	}
	return (false);
}



template <>
inline void
gensys_modules< bigint, double, VariationI >::E_convert_lattice_A (dense_alg< bigint > & da,
								   double** dblvalue)
{
	for (lidia_size_t i = 0; i < da.b.rank;) {
		is_zero = true;
		for (lidia_size_t j = 0; j < da.b.columns; j++)
			if ((!is_zero) || (!da.s.value[i][j].is_zero())) {
				tempbfl.assign(da.s.value[i][j]);
				shift_right(tempbfl, tempbfl, da.d.bit_factor);
				tempbfl.doublify(dblvalue[i][j]);
				is_zero = false;
			}
			else
				dblvalue[i][j] = 0.0;
		if (is_zero) {
			da.b.rank--;
			for (lidia_size_t j = i; j < da.b.rank; j++) {
				tempP = static_cast<void *>(dblvalue[j]);
				dblvalue[j] = dblvalue[j+1];
				dblvalue[j+1] = static_cast<double *>(tempP);
				tempP = static_cast<void *>(da.s.value[j]);
				da.s.value[j] = da.s.value[j+1];
				da.s.value[j+1] = static_cast<bigint *>(tempP);
			}
		}
		else
			i++;
	}
}



template <>
inline bool
gensys_modules< bigint, double, VariationI >::E_convert_vector_A (dense_alg< bigint > & da,
								  double** dblvalue,
								  const lidia_size_t k)
{
	is_zero = true;
	for (lidia_size_t i = 0; i < da.b.columns; i++)
		if ((!is_zero) || (!da.s.value[k][i].is_zero())) {
			tempbfl.assign(da.s.value[k][i]);
			shift_right(tempbfl, tempbfl, da.d.bit_factor);
			tempbfl.doublify(dblvalue[k][i]);
			is_zero = false;
		}
		else
			dblvalue[k][i] = 0.0;
	if (is_zero) {
		da.b.rank--;
		for (lidia_size_t i = k; i < da.b.rank; i++) {
			tempP = static_cast<void *>(dblvalue[i]);
			dblvalue[i] = dblvalue[i+1];
			dblvalue[i+1] = static_cast<double *>(tempP);
			tempP = static_cast<void *>(da.s.value[i]);
			da.s.value[i] = da.s.value[i+1];
			da.s.value[i+1] = static_cast<bigint *>(tempP);
		}
		return(true);
	}
	return(false);
}



template <>
inline void
gensys_modules< bigint, xdouble, VariationI >::E_convert_lattice_A (dense_alg< bigint > & da,
								    xdouble** xdblvalue)
{
	for (lidia_size_t i = 0; i < da.b.rank;) {
		is_zero = true;
		for (lidia_size_t j = 0; j < da.b.columns; j++)
			if ((!is_zero) || (!da.s.value[i][j].is_zero())) {
				tempbfl.assign(da.s.value[i][j]);
				shift_right(tempbfl, tempbfl, da.d.bit_factor);
				tempbfl.xdoublify(xdblvalue[i][j]);
				is_zero = false;
			}
			else
				xdblvalue[i][j] = 0.0;
		if (is_zero) {
			da.b.rank--;
			for (lidia_size_t j = i; j < da.b.rank; j++) {
				tempP = static_cast<void *>(xdblvalue[j]);
				xdblvalue[j] = xdblvalue[j+1];
				xdblvalue[j+1] = static_cast<xdouble *>(tempP);
				tempP = static_cast<void *>(da.s.value[j]);
				da.s.value[j] = da.s.value[j+1];
				da.s.value[j+1] = static_cast<bigint *>(tempP);
			}
		}
		else
			i++;
	}
}



template <>
inline bool
gensys_modules< bigint, xdouble, VariationI >::E_convert_vector_A (dense_alg< bigint > & da,
								   xdouble** xdblvalue,
								   const lidia_size_t k)
{
	is_zero = true;
	for (lidia_size_t i = 0; i < da.b.columns; i++)
		if ((!is_zero) || (!da.s.value[k][i].is_zero())) {
			tempbfl.assign(da.s.value[k][i]);
			shift_right(tempbfl, tempbfl, da.d.bit_factor);
			tempbfl.xdoublify(xdblvalue[k][i]);
			is_zero = false;
		}
		else
			xdblvalue[k][i] = 0.0;
	if (is_zero) {
		da.b.rank--;
		for (lidia_size_t i = k; i < da.b.rank; i++) {
			tempP = static_cast<void *>(xdblvalue[i]);
			xdblvalue[i] = xdblvalue[i+1];
			xdblvalue[i+1] = static_cast<xdouble *>(tempP);
			tempP = static_cast<void *>(da.s.value[i]);
			da.s.value[i] = da.s.value[i+1];
			da.s.value[i+1] = static_cast<bigint *>(tempP);
		}
		return(true);
	}
	return(false);
}



template <>
inline void
gensys_modules< bigint, bigfloat, VariationI >::E_convert_lattice_A (dense_alg< bigint > & da,
								     bigfloat** bflvalue)
{
	for (lidia_size_t i = 0; i < da.b.rank;) {
		is_zero = true;
		for (lidia_size_t j = 0; j < da.b.columns; j++)
			if ((!is_zero) || (!da.s.value[i][j].is_zero())) {
				bi_bit_len = da.s.value[i][j].bit_length();
				if (bi_bit_len > da.d.cut_bit_prec) {
					bit_diff = bi_bit_len-da.d.cut_bit_prec;
					shift_right(tempbin, da.s.value[i][j], bit_diff);
					bflvalue[i][j].assign(tempbin);
					shift_left(bflvalue[i][j], bflvalue[i][j], bit_diff);
				}
				else
					bflvalue[i][j].assign(da.s.value[i][j]);
				shift_right(bflvalue[i][j], bflvalue[i][j], da.d.bit_factor);
				is_zero = false;
			}
			else
				bflvalue[i][j].assign_zero();
		if (is_zero) {
			da.b.rank--;
			for (lidia_size_t j = i; j < da.b.rank; j++) {
				tempP = static_cast<void *>(bflvalue[j]);
				bflvalue[j] = bflvalue[j+1];
				bflvalue[j+1] = static_cast<bigfloat *>(tempP);
				tempP = static_cast<void *>(da.s.value[j]);
				da.s.value[j] = da.s.value[j+1];
				da.s.value[j+1] = static_cast<bigint *>(tempP);
			}
		}
		else
			i++;
	}
}



template <>
inline bool
gensys_modules< bigint, bigfloat, VariationI >::E_convert_vector_A (dense_alg< bigint > & da,
								    bigfloat** bflvalue,
								    const lidia_size_t k)
{
	is_zero = true;
	for (lidia_size_t i = 0; i < da.b.columns; i++)
		if ((!is_zero) || (!da.s.value[k][i].is_zero())) {
			bi_bit_len = da.s.value[k][i].bit_length();
			if (bi_bit_len > da.d.cut_bit_prec) {
				bit_diff = bi_bit_len-da.d.cut_bit_prec;
				shift_right(tempbin, da.s.value[k][i], bit_diff);
				bflvalue[k][i].assign(tempbin);
				shift_left(bflvalue[k][i], bflvalue[k][i], bit_diff);
			}
			else
				bflvalue[k][i].assign(da.s.value[k][i]);
			shift_right(bflvalue[k][i], bflvalue[k][i], da.d.bit_factor);
			is_zero = false;
		}
		else
			bflvalue[k][i].assign_zero();
	if (is_zero) {
		da.b.rank--;
		for (lidia_size_t i = k; i < da.b.rank; i++) {
			tempP = static_cast<void *>(bflvalue[i]);
			bflvalue[i] = bflvalue[i+1];
			bflvalue[i+1] = static_cast<bigfloat *>(tempP);
			tempP = static_cast<void *>(da.s.value[i]);
			da.s.value[i] = da.s.value[i+1];
			da.s.value[i+1] = static_cast<bigint *>(tempP);
		}
		return(true);
	}
	return(false);
}



//
// End of
// 	     * base_modules < bigint, {double, xdouble, bigfloat}, VariationI >
// 	     * basis_modules < bigint, {double, xdouble, bigfloat}, VariationI >
// 	     * gensys_modules < bigint, {double, xdouble, bigfloat}, VariationI >
//

//
// Code for
// 	     * base_modules < bigfloat, {double, xdouble, bigfloat}, Normal >
// 	     * basis_modules < bigfloat, {double, xdouble, bigfloat}, Normal >
// 	     * gensys_modules < bigfloat, {double, xdouble, bigfloat}, Normal >
//
template <>
inline void
base_modules< bigfloat, double, Normal >::E_convert_value_A (dense_alg< bigfloat > &,
							     double& d,
							     const bigfloat& bfl)
{
	bfl.doublify(d);
}



template <>
inline bool
base_modules< bigfloat, double, Normal >::A_convert_value_E_bound (dense_alg< bigfloat > &,
								   bigfloat& bfl,
								   const double& dbl,
								   const double& bound)
{
	if (fabs(dbl) > bound) {
		bfl.assign(dbl);
		return (true);
	}
	else {
		bfl.assign(static_cast<sdigit>(dbl));
		return (false);
	}
}



template <>
inline void
base_modules< bigfloat, xdouble, Normal >::E_convert_value_A (dense_alg< bigfloat > &,
							      xdouble& xd,
							      const bigfloat& bfl)
{
	bfl.xdoublify(xd);
}



template <>
inline bool
base_modules< bigfloat, xdouble, Normal >::A_convert_value_E_bound (dense_alg< bigfloat > &,
								    bigfloat& bfl,
								    const xdouble& xdbl,
								    const xdouble& bound)
{
	bfl.assign(xdbl);
	if (fabs(xdbl) > bound)
		return (true);
	return (false);
}



template <>
inline void
base_modules< bigfloat, bigfloat, Normal >::E_convert_value_A (dense_alg< bigfloat > & da,
							       bigfloat& bfla,
							       const bigfloat& bfl)
{
	bigfloat::set_precision(da.d.approx_prec);
	mant.assign(bfl.mantissa());
	expo = bfl.exponent();
	bi_bit_len = mant.bit_length();
	if (bi_bit_len > da.d.cut_bit_prec) {
		bit_diff = bi_bit_len-da.d.cut_bit_prec;
		shift_right(mant, mant, bit_diff);
		bfla.assign(mant);
		shift_left(bfla, bfla, bit_diff);
	}
	else
		bfla.assign(mant);
	shift_left(bfla, bfla, expo);
	bigfloat::set_precision(da.d.exact_prec);
}



template <>
inline bool
base_modules< bigfloat, bigfloat, Normal >::A_convert_value_E_bound (dense_alg< bigfloat > &,
								     bigfloat& bfla,
								     const bigfloat& bfl,
								     const bigfloat& bound)
{
	bfla.assign(bfl);
	if (abs(bfl).compare(bound) > 0)
		return (true);
	return (false);
}



template <>
inline void
basis_modules< bigfloat, double, Normal >::E_convert_lattice_A (dense_alg< bigfloat > & da,
								double** dblvalue)
{
	for (lidia_size_t i = 0; i < da.b.rows; i++)
		for (lidia_size_t j = 0; j < da.b.columns; j++)
			da.s.value[i][j].doublify(dblvalue[i][j]);
}



template <>
inline bool
basis_modules< bigfloat, double, Normal >::E_convert_vector_A (dense_alg< bigfloat > & da,
							       double **dblvalue,
							       const lidia_size_t k)
{
	for (lidia_size_t i = 0; i < da.b.columns; i++)
		da.s.value[k][i].doublify(dblvalue[k][i]);
	return(false);
}



template <>
inline void
basis_modules< bigfloat, xdouble, Normal >::E_convert_lattice_A (dense_alg< bigfloat > & da,
								 xdouble** xdblvalue)
{
	for (lidia_size_t i = 0; i < da.b.rows; i++)
		for (lidia_size_t j = 0; j < da.b.columns; j++)
			da.s.value[i][j].xdoublify(xdblvalue[i][j]);
}



template <>
inline bool
basis_modules< bigfloat, xdouble, Normal >::E_convert_vector_A (dense_alg< bigfloat > & da,
								xdouble **xdblvalue,
								const lidia_size_t k)
{
	for (lidia_size_t i = 0; i < da.b.columns; i++)
		da.s.value[k][i].xdoublify(xdblvalue[k][i]);
	return(false);
}



template <>
inline void
basis_modules< bigfloat, bigfloat, Normal >::E_convert_lattice_A (dense_alg< bigfloat > & da,
								  bigfloat** bflvalue)
{
	bigfloat::set_precision(da.d.approx_prec);
	for (lidia_size_t i = 0; i < da.b.rows; i++)
		for (lidia_size_t j = 0; j < da.b.columns; j++) {
			mant.assign(da.s.value[i][j].mantissa());
			expo = da.s.value[i][j].exponent();
			bi_bit_len = mant.bit_length();
			if (bi_bit_len > da.d.cut_bit_prec) {
				bit_diff = bi_bit_len-da.d.cut_bit_prec;
				shift_right(mant, mant, bit_diff);
				bflvalue[i][j].assign(mant);
				shift_left(bflvalue[i][j], bflvalue[i][j], bit_diff);
			}
			else
				bflvalue[i][j].assign(mant);
			shift_left(bflvalue[i][j], bflvalue[i][j], expo);
		}
	bigfloat::set_precision(da.d.exact_prec);
}



template <>
inline bool
basis_modules< bigfloat, bigfloat, Normal >::E_convert_vector_A (dense_alg< bigfloat > & da,
								 bigfloat** bflvalue,
								 const lidia_size_t k)
{
	bigfloat::set_precision(da.d.approx_prec);
	for (lidia_size_t i = 0; i < da.b.columns; i++) {
		mant.assign(da.s.value[k][i].mantissa());
		expo = da.s.value[k][i].exponent();
		bi_bit_len = mant.bit_length();
		if (bi_bit_len > da.d.cut_bit_prec) {
			bit_diff = bi_bit_len-da.d.cut_bit_prec;
			shift_right(mant, mant, bit_diff);
			bflvalue[k][i].assign(mant);
			shift_left(bflvalue[k][i], bflvalue[k][i], bit_diff);
		}
		else
			bflvalue[k][i].assign(mant);
		shift_left(bflvalue[k][i], bflvalue[k][i], expo);
	}
	bigfloat::set_precision(da.d.exact_prec);
	return (false);
}



template <>
inline void
gensys_modules< bigfloat, double, Normal >::E_convert_lattice_A (dense_alg< bigfloat > & da,
								 double** dblvalue)
{
	for (lidia_size_t i = 0; i < da.b.rank;) {
		is_zero = true;
		for (lidia_size_t j = 0; j < da.b.columns; j++)
			if ((!is_zero) || (!da.s.value[i][j].is_zero())) {
				da.s.value[i][j].doublify(dblvalue[i][j]);
				is_zero = false;
			}
			else
				dblvalue[i][j] = 0.0;
		if (is_zero) {
			da.b.rank--;
			for (lidia_size_t j = i; j < da.b.rank; j++) {
				tempP = static_cast<void *>(dblvalue[j]);
				dblvalue[j] = dblvalue[j+1];
				dblvalue[j+1] = static_cast<double *>(tempP);
				tempP = static_cast<void *>(da.s.value[j]);
				da.s.value[j] = da.s.value[j+1];
				da.s.value[j+1] = static_cast<bigfloat *>(tempP);
			}
		}
		else
			i++;
	}
}



template <>
inline bool
gensys_modules< bigfloat, double, Normal >::E_convert_vector_A (dense_alg< bigfloat > & da,
								double** dblvalue,
								const lidia_size_t k)
{
	is_zero = true;
	for (lidia_size_t i = 0; i < da.b.columns; i++)
		if ((!is_zero) || (!da.s.value[k][i].is_zero())) {
			da.s.value[k][i].doublify(dblvalue[k][i]);
			is_zero = false;
		}
		else
			dblvalue[k][i] = 0.0;
	if (is_zero) {
		da.b.rank--;
		for (lidia_size_t i = k; i < da.b.rank; i++) {
			tempP = static_cast<void *>(dblvalue[i]);
			dblvalue[i] = dblvalue[i+1];
			dblvalue[i+1] = static_cast<double *>(tempP);
			tempP = static_cast<void *>(da.s.value[i]);
			da.s.value[i] = da.s.value[i+1];
			da.s.value[i+1] = static_cast<bigfloat *>(tempP);
		}
		return(true);
	}
	return(false);
}



template <>
inline void
gensys_modules< bigfloat, xdouble, Normal >::E_convert_lattice_A (dense_alg< bigfloat > & da,
								  xdouble** xdblvalue)
{
	for (lidia_size_t i = 0; i < da.b.rank;) {
		is_zero = true;
		for (lidia_size_t j = 0; j < da.b.columns; j++)
			if ((!is_zero) || (!da.s.value[i][j].is_zero())) {
				da.s.value[i][j].xdoublify(xdblvalue[i][j]);
				is_zero = false;
			}
			else
				xdblvalue[i][j] = 0.0;
		if (is_zero) {
			da.b.rank--;
			for (lidia_size_t j = i; j < da.b.rank; j++) {
				tempP = static_cast<void *>(xdblvalue[j]);
				xdblvalue[j] = xdblvalue[j+1];
				xdblvalue[j+1] = static_cast<xdouble *>(tempP);
				tempP = static_cast<void *>(da.s.value[j]);
				da.s.value[j] = da.s.value[j+1];
				da.s.value[j+1] = static_cast<bigfloat *>(tempP);
			}
		}
		else
			i++;
	}
}



template <>
inline bool
gensys_modules< bigfloat, xdouble, Normal >::E_convert_vector_A (dense_alg< bigfloat > & da,
								 xdouble** xdblvalue,
								 const lidia_size_t k)
{
	is_zero = true;
	for (lidia_size_t i = 0; i < da.b.columns; i++)
		if ((!is_zero) || (!da.s.value[k][i].is_zero())) {
			da.s.value[k][i].xdoublify(xdblvalue[k][i]);
			is_zero = false;
		}
		else
			xdblvalue[k][i] = 0.0;
	if (is_zero) {
		da.b.rank--;
		for (lidia_size_t i = k; i < da.b.rank; i++) {
			tempP = static_cast<void *>(xdblvalue[i]);
			xdblvalue[i] = xdblvalue[i+1];
			xdblvalue[i+1] = static_cast<xdouble *>(tempP);
			tempP = static_cast<void *>(da.s.value[i]);
			da.s.value[i] = da.s.value[i+1];
			da.s.value[i+1] = static_cast<bigfloat *>(tempP);
		}
		return(true);
	}
	return(false);
}



template <>
inline void
gensys_modules< bigfloat, bigfloat, Normal >::E_convert_lattice_A (dense_alg< bigfloat > & da,
								   bigfloat** bflvalue)
{
	bigfloat::set_precision(da.d.approx_prec);
	for (lidia_size_t i = 0; i < da.b.rank;) {
		is_zero = true;
		for (lidia_size_t j = 0; j < da.b.columns; j++)
			if ((!is_zero) || (!da.s.value[i][j].is_zero())) {
				mant.assign(da.s.value[i][j].mantissa());
				expo = da.s.value[i][j].exponent();
				bi_bit_len = mant.bit_length();
				if (bi_bit_len > da.d.cut_bit_prec) {
					bit_diff = bi_bit_len-da.d.cut_bit_prec;
					shift_right(mant, mant, bit_diff);
					bflvalue[i][j].assign(mant);
					shift_left(bflvalue[i][j], bflvalue[i][j], bit_diff);
				}
				else
					bflvalue[i][j].assign(mant);
				shift_left(bflvalue[i][j], bflvalue[i][j], expo);
				is_zero = false;
			}
			else
				bflvalue[i][j].assign_zero();
		if (is_zero) {
			da.b.rank--;
			for (lidia_size_t j = i; j < da.b.rank; j++) {
				tempP = static_cast<void *>(bflvalue[j]);
				bflvalue[j] = bflvalue[j+1];
				bflvalue[j+1] = static_cast<bigfloat *>(tempP);
				tempP = static_cast<void *>(da.s.value[j]);
				da.s.value[j] = da.s.value[j+1];
				da.s.value[j+1] = static_cast<bigfloat *>(tempP);
			}
		}
		else
			i++;
	}
	bigfloat::set_precision(da.d.exact_prec);
}



template <>
inline bool
gensys_modules< bigfloat, bigfloat, Normal >::E_convert_vector_A (dense_alg< bigfloat > & da,
								  bigfloat** bflvalue,
								  const lidia_size_t k)
{
	bigfloat::set_precision(da.d.approx_prec);
	is_zero = true;
	for (lidia_size_t i = 0; i < da.b.columns; i++)
		if ((!is_zero) || (!da.s.value[k][i].is_zero())) {
			mant.assign(da.s.value[k][i].mantissa());
			expo = da.s.value[k][i].exponent();
			bi_bit_len = mant.bit_length();
			if (bi_bit_len > da.d.cut_bit_prec) {
				bit_diff = bi_bit_len-da.d.cut_bit_prec;
				shift_right(mant, mant, bit_diff);
				bflvalue[k][i].assign(mant);
				shift_left(bflvalue[k][i], bflvalue[k][i], bit_diff);
			}
			else
				bflvalue[k][i].assign(mant);
			shift_left(bflvalue[k][i], bflvalue[k][i], expo);
			is_zero = false;
		}
		else
			bflvalue[k][i].assign_zero();
	if (is_zero) {
		da.b.rank--;
		for (lidia_size_t i = k; i < da.b.rank; i++) {
			tempP = static_cast<void *>(bflvalue[i]);
			bflvalue[i] = bflvalue[i+1];
			bflvalue[i+1] = static_cast<bigfloat *>(tempP);
			tempP = static_cast<void *>(da.s.value[i]);
			da.s.value[i] = da.s.value[i+1];
			da.s.value[i+1] = static_cast<bigfloat *>(tempP);
		}
		bigfloat::set_precision(da.d.exact_prec);
		return(true);
	}
	bigfloat::set_precision(da.d.exact_prec);
	return(false);
}



//
// End of
// 	     * base_modules < bigfloat, {double, xdouble, bigfloat}, Normal >
// 	     * basis_modules < bigfloat, {double, xdouble, bigfloat}, Normal >
// 	     * gensys_modules < bigfloat, {double, xdouble, bigfloat}, Normal >
//



#ifdef LIDIA_NAMESPACE
# ifndef IN_NAMESPACE_LIDIA
}	// end of namespace LiDIA
# endif
#endif



#endif	// LIDIA_LATTICE_MODULES_CC_GUARD_
