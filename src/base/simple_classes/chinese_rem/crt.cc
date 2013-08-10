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
//	$Id: crt.cc,v 2.6 2004/06/15 10:19:26 lidiaadm Exp $
//
//	Author	: Frank J. Lehmann (FL), Thomas Pfahler (TPF)
//	Changes	: See CVS log
//
//==============================================================================================


#ifdef HAVE_CONFIG_H
# include	"config.h"
#endif
#include	"LiDIA/bigint.h"
#include	"LiDIA/bigmod.h"
#include	"LiDIA/base_vector.h"
#include	"LiDIA/crt_table.h"
#include	"LiDIA/crt.h"



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



//***************************************************************
//
//					class crt
//
//***************************************************************



const char crt::INIT = 0;
const char crt::SINGLE = 1;
const char crt::POINTER = 2;
const char crt::VECTOR = 4;
const char crt::MATRIX = 8;

const char crt::INTEGER = 1;
const char crt::MODULAR = 2;



crt::crt ()
{
	debug_handler("crt", "crt()");

	m = NULL;

	H = d = 0;
	Hv = NULL;
	dv = NULL;

	ix = NULL;

	size = used = 0;

	mode = INIT;
}



crt::crt (crt_table &M)
{
	debug_handler("crt", "crt(crt_table &)");

	if (M.num_primes == 0)
		lidia_error_handler("crt", "crt(crt_table)::argument crt_table must be initialized before");

	m = NULL;

	H = d = 0;
	Hv = NULL;
	dv = NULL;

	ix = NULL;

	init (M);
}



void
crt::clear ()
{
	debug_handler("crt", "clear ()");

	if (m != NULL) {
		m->reference_counter --;
		m = NULL;
	}

	H = d = 0;

	clear_pointer();
	clear_vector();
	clear_matrix();

	if (ix != NULL) {
		delete [] ix;
		ix = NULL;
	}

	used = 0;

	mode = INIT;
}



void
crt::reset ()
{
	debug_handler("crt", "reset ()");

	lidia_size_t i;

	H = d = 0;

	clear_pointer();
	clear_vector();
	clear_matrix();


	if (m != NULL) // i.e. ix != NULL and ix has appropriate size
	{
		for (i = 0; i < m->num_primes; i++) ix[i] = 0;
	}

	used = 0;

	mode = INIT;
}



void
crt::init (crt_table & M)
{
	debug_handler("crt", "init (crt_table &)");

	lidia_size_t i;

	if (M.num_primes == 0)
		lidia_error_handler("crt", "init(crt_table &)::argument crt_table must be initialized before");


	clear ();

	m = &M;
	m->reference_counter ++;

	ix = new unsigned char [ m->num_primes ];
	for (i = 0; i < m->num_primes; i++) ix[i] = 0;

}



void
crt::init_single ()
{
	debug_handler("crt", "init_single ()");

	H = d = 0;
}



void
crt::init_pointer (lidia_size_t sz)
{
	debug_handler("crt", "init_pointer (lidia_size_t)");

	lidia_size_t i;

	if (sz < 0) {
		lidia_error_handler("crt", "init_pointer (lidia_size_t)::negative length");
	}
	else {
		if (Hv != NULL)
			delete[] Hv;
		if (dv != NULL)
			delete[] dv;

		Hv = new bigint[sz];
		dv = new double[sz];

		if (Hv == NULL)
			lidia_error_handler("crt", "init_pointer(lidia_size_t)::out of memory (Hv)");
		if (dv == NULL)
			lidia_error_handler("crt", "init_pointer(lidia_size_t)::out of memory (dv)");

		size = sz;

		for (i = 0; i < sz; i++) {
			Hv[i] = 0;
			dv[i] = 0;
		}
	}
}



void
crt::clear_pointer ()
{
	debug_handler("crt", "clear_pointer()");

	if (Hv != NULL)
		delete[] Hv;
	if (dv != NULL)
		delete[] dv;

	Hv = NULL;
	dv = NULL;
	size = 0;
}



void
crt::init_vector (lidia_size_t sz)
{
	debug_handler("crt", "init_vector (lidia_size_t)");

	lidia_size_t i;

	if (sz < 0) {
		lidia_error_handler("crt", "init_vector (lidia_size_t)::negative length");
	}
	else {
		Hvec.set_capacity (sz);
		dvec.set_capacity (sz);

		for (i = 0; i < sz; i++) {
			Hvec[i] = 0;
			dvec[i] = 0;
		}
	}
}



void
crt::clear_vector ()
{
	debug_handler("crt", "clear_vector ()");

	Hvec.kill ();
	dvec.kill ();
}



void
crt::init_matrix (lidia_size_t row, lidia_size_t col)
{

	debug_handler("crt", "init_matrix (lidia_size_t, lidia_size_t)");

	lidia_size_t r, c;
	double null = 0;
	bigint zero(0);

	if (row < 0 || col < 0) {
		lidia_error_handler("crt", "init_vector (lidia_size_t, lidia_size_t)::negative matrix dimension");
	}
	else {
		Hmat.resize (row, col);
		dmat.resize (row, col);

		for (r = 0; r < row; r ++)
			for (c = 0; c < col; c ++) {
				Hmat.sto (r, c , zero);
				dmat.sto (r, c , null);
			}
	}

}



void
crt::clear_matrix ()
{

	debug_handler("crt", "clear_matrix ()");

	Hmat.resize (1, 1);
	dmat.resize (1, 1);

}



crt::~crt()
{
	debug_handler("crt", "~crt ()");

	// only do part of what clear would do!
	if (m != NULL)
		m->reference_counter --;

	clear_pointer();

	if (ix != NULL)
		delete [] ix;
}



sdigit
crt::get_prime (lidia_size_t index) const
{
	debug_handler("crt", "get_prime(lidia_size_t)");

	if (index< 0 || index >= m->num_primes)
		lidia_error_handler("crt", "get_prime(lidia_size_t)::index out of range");

	return (m->primes[index]);
}



bool
crt::test_mode (char md)
{
	debug_handler("crt", "test_mode(char)");

	bool ok;

	if (mode == md) {
		ok = true;
	}
	else if (mode == INIT) {
		mode = md;
		ok = true;
	}
	else {
		ok = false;
	}

	return ok;
}



//***************************************************************
//
//                 reduction for INTEGER mode
//
//***************************************************************

void
crt::reduce (sdigit &small_value, const bigint &big, lidia_size_t index) const
{
	debug_handler("crt", "reduce(sdigit&, bigint&, lidia_size_t ");

	if (index< 0 || index > m->num_primes)
		lidia_error_handler("crt", "reduce(sdigit&, bigint&, lidia_size_t)::index out of range");
#if 0
	if (m->mode != INTEGER)
		lidia_error_handler("crt", "reduce(sdigit&, bigint&, lidia_size_t)::table not in INTEGER mode");
#endif
	m->SINGLE_REMAINDER(small_value, big, index);
}



void
crt::reduce (sdigit* small_vec, const bigint* big_vec, lidia_size_t bsize, lidia_size_t index) const
{
	debug_handler("crt", "reduce(sdigit*, bigint*, lidia_size_t, lidia_size_t) ");

	if (index< 0 || index > m->num_primes)
		lidia_error_handler("crt", "reduce(sdigit*, bigint*, lidia_size_t, lidia_size_t)::index out of range");
#if 0
	if (m->mode != INTEGER)
		lidia_error_handler("crt", "reduce(sdigit*, bigint*, lidia_size_t, lidia_size_t)::table not in INTEGER mode");
#endif
	m->MULTI_REMAINDER(small_vec, big_vec, bsize, index);
}



void
crt::reduce (base_vector< sdigit > & small_vec,
	     const base_vector< bigint > & big_vec, lidia_size_t index) const
{
	debug_handler("crt", "reduce(base_vector< sdigit > &, const base_vector< bigint > &, lidia_size_t) ");

	if (index< 0 || index > m->num_primes)
		lidia_error_handler("crt", "reduce(base_vector< sdigit > &, const base_vector< bigint > &, lidia_size_t)::index out of range");

#if 0
	if (m->mode != INTEGER)
		lidia_error_handler("crt", "reduce(base_vector< sdigit > &, const base_vector< bigint > &, lidia_size_t)::table not in INTEGER mode");
#endif

	lidia_size_t     bsize;
	sdigit *small_ptr;
	bigint *big_ptr;

	bsize = big_vec.size();
	small_vec.set_capacity (bsize);
	small_vec.set_size     (bsize);

	small_ptr = small_vec.get_data_address();
	big_ptr = big_vec.get_data_address();

	m->MULTI_REMAINDER(small_ptr, big_ptr, bsize , index);
}



void
crt::reduce (base_matrix< sdigit > & small_mat,
	     const base_matrix< bigint > & big_mat, lidia_size_t index) const
{
	debug_handler("crt", "reduce(base_matrix< sdigit > &, const base_matrix< bigint > &, lidia_size_t) ");

	if (index< 0 || index > m->num_primes)
		lidia_error_handler("crt", "reduce(base_matrix< sdigit > &, const base_matrix< bigint > &, lidia_size_t)::index out of range");

#if 0
	if (m->mode != INTEGER)
		lidia_error_handler("crt", "reduce(base_matrix< sdigit > &, const base_matrix< bigint > &, lidia_size_t)::table not in INTEGER mode");
#endif

	lidia_size_t    brow, bcol, i;

	base_vector< bigint > bvec;
	bigint *bptr;

	base_vector< sdigit > svec;
	sdigit *sptr;

	brow = big_mat.get_no_of_rows   ();
	bcol = big_mat.get_no_of_columns ();

	small_mat.resize (brow, bcol);

	svec.set_capacity (bcol);
	svec.set_size     (bcol);
	sptr = svec.get_data_address ();

	bvec.set_capacity (bcol);
	bvec.set_size     (bcol);
	bptr = bvec.get_data_address ();

	if (sptr == NULL)
		lidia_error_handler ("crt", "reduce(base_matrix< sdigit > &, const base_matrix< bigint > &, lidia_size_t)::out of memory");

	for (i = 0; i < brow; i++) {
		big_mat.get_row_vector (bvec , i);
		m->MULTI_REMAINDER(sptr, bptr, bcol , index);

		small_mat.sto_row_vector (svec , bcol , i);
	}

	// free ( sptr );
}



//***************************************************************
//
//                 reduction for MODULAR mode
//
//***************************************************************

void
crt::reduce (sdigit &small_value, const bigmod &big, lidia_size_t index) const
{
	debug_handler("crt", "reduce(sdigit&, bigmod&, lidia_size_t ");

	if (index< 0 || index > m->num_primes)
		lidia_error_handler("crt", "reduce(sdigit&, bigmod&, lidia_size_t)::index out of range");

#if 0
	if (m->mode != MODULAR)
		lidia_error_handler("crt", "reduce(sdigit&, bigmod&, lidia_size_t)::table not in MODULAR mode");
#endif

	m->SINGLE_REMAINDER(small_value, big, index);
}



void
crt::reduce (sdigit* small_vec, const bigmod* big_vec, lidia_size_t bsize, lidia_size_t index) const
{
	debug_handler("crt", "reduce(sdigit*, bigmod*, lidia_size_t, lidia_size_t) ");

	if (index< 0 || index > m->num_primes)
		lidia_error_handler("crt", "reduce(sdigit*, bigmod*, lidia_size_t, lidia_size_t)::index out of range");

#if 0
	if (m->mode != MODULAR)
		lidia_error_handler("crt", "reduce(sdigit*, bigmod*, lidia_size_t, lidia_size_t)::table not in MODULAR mode");
#endif

	m->MULTI_REMAINDER(small_vec, big_vec, bsize, index);
}



void
crt::reduce (base_vector< sdigit > & small_vec,
	     const base_vector< bigmod > & big_vec, lidia_size_t index) const
{
	debug_handler("crt", "reduce(base_vector< sdigit > &, const base_vector< bigmod > &, lidia_size_t) ");

	if (index< 0 || index > m->num_primes)
		lidia_error_handler("crt", "reduce(base_vector< sdigit > &, const base_vector< bigmod > &, lidia_size_t)::index out of range");

#if 0
	if (m->mode != MODULAR)
		lidia_error_handler("crt", "reduce(base_vector< sdigit > &, const base_vector< bigmod > &, lidia_size_t)::table not in MODULAR mode");
#endif

	lidia_size_t     bsize;
	sdigit *small_ptr;
	bigmod *big_ptr;

	bsize = big_vec.size();
	small_vec.set_capacity (bsize);
	small_vec.set_size     (bsize);

	small_ptr = small_vec.get_data_address();
	big_ptr = big_vec.get_data_address();

	m->MULTI_REMAINDER(small_ptr, big_ptr, bsize , index);
}



void
crt::reduce (base_matrix< sdigit > & small_mat,
	     const base_matrix< bigmod > & big_mat, lidia_size_t index) const
{
	debug_handler("crt", "reduce(base_matrix< sdigit > &, const base_matrix< bigmod > &, lidia_size_t) ");

	if (index< 0 || index > m->num_primes)
		lidia_error_handler("crt", "reduce(base_matrix< sdigit > &, const base_matrix< bigmod > &, lidia_size_t)::index out of range");

#if 0
	if (m->mode != MODULAR)
		lidia_error_handler("crt", "reduce(base_matrix< sdigit > &, const base_matrix< bigmod > &, lidia_size_t)::table not in MODULAR mode");
#endif

	lidia_size_t    brow, bcol, i;

	base_vector< bigmod > bvec;
	bigmod *bptr;

	base_vector< sdigit > svec;
	sdigit *sptr;

	brow = big_mat.get_no_of_rows   ();
	bcol = big_mat.get_no_of_columns ();

	small_mat.resize (brow, bcol);

	svec.set_capacity (bcol);
	svec.set_size     (bcol);
	sptr = svec.get_data_address();

	bvec.set_capacity (bcol);
	bvec.set_size     (bcol);
	bptr = bvec.get_data_address ();

	if (sptr == NULL)
		lidia_error_handler ("crt", "reduce(base_matrix< sdigit > &, const base_matrix< bigmod > &, lidia_size_t)::out of memory");

	for (i = 0; i < brow; i++) {
		big_mat.get_row_vector (bvec , i);
		m->MULTI_REMAINDER(sptr, bptr, bcol , index);

		small_mat.sto_row_vector (svec , bcol , i);
	}

	// free ( sptr );
}



//***************************************************************
//
//         combination for both, INTEGER and MODULAR mode
//
//***************************************************************

void
crt::combine (const sdigit & small_value, lidia_size_t index)
{
	debug_handler("crt", "combine(const sdigit&, lidia_size_t)");

	if (m == NULL)
		lidia_error_handler("crt", "combine(const sdigit&, lidia_size_t)::no reference to crt_table yet");

	if (index< 0 || index > m->num_primes)
		lidia_error_handler("crt", "combine(const sdigit&, lidia_size_t)::index out of range");

	if (! test_mode (SINGLE))
		lidia_error_handler("crt", "combine(const sdigit&, lidia_size_t)::object may not be used with single instances");


	bigint tmp;

	if (! used)
		init_single ();

	if (ix[index] == 0) {
		multiply(tmp, m->coeff_mod_p[index], small_value);
		add     (H, H, tmp);

		d += (static_cast<double>(small_value)) * m->x[index];

		ix[index]++;
		used++;
	}
	else {
		lidia_error_handler("crt", "combine(const sdigit&, lidia_size_t)::prime[lidia_size_t] used twice");
	}
}



void
crt::combine (const sdigit* small_vec, lidia_size_t l, lidia_size_t index)
{
	debug_handler("crt", "combine(const sdigit*, lidia_size_t, lidia_size_t)");

	if (m == NULL)
		lidia_error_handler("crt", "combine(const sdigit*, lidia_size_t, lidia_size_t)::no reference to crt_table yet");

	if (index< 0 || index > m->num_primes)
		lidia_error_handler("crt", "combine(const sdigit*, lidia_size_t, lidia_size_t)::index out of range");

	if (used && l != size)
		lidia_error_handler("crt", "combine(const sdigit*, lidia_size_t, lidia_size_t)::incompatible length of pointer object");

	if (! test_mode (POINTER))
		lidia_error_handler("crt", "combine(const sdigit*, lidia_size_t, lidia_size_t)::object may not be used with pointers");

	lidia_size_t i;
	bigint tmp;

	if (! used)
		init_pointer (l);

	if (ix[index] == 0) {
		for (i = 0; i < l; i++) {
			multiply(tmp, m->coeff_mod_p[index], small_vec[i]);
			add     (Hv[i], Hv[i], tmp);

			dv[i] += static_cast<double>(small_vec[i]) * m->x[index];
		}

		ix[index]++;
		used++;
	}
	else {
		lidia_error_handler("crt", "combine(const sdigit*, lidia_size_t, lidia_size_t)::prime[lidia_size_t] used twice");
	}
}



void
crt::combine (const udigit* small_vec, lidia_size_t l, lidia_size_t index)
{
	// This is the (udigit*) version of combine.
	//
	// IMPORTANT NOTE:
	// ---------------
	// We need small_vec[i] = (long(small_vec[i])) for all i, i.e. the high
	// order bit of small_vec[i] must not be set (otherwise we get wrong results
	// for "multiply(...)" below).
	// This _should_ not be a problem, though, because prime[index] is of type
	// long, and we assume small_vec[i] reduced modulo prime[index].

	debug_handler("crt", "combine(const udigit*, lidia_size_t, lidia_size_t)");

	if (m == NULL)
		lidia_error_handler("crt", "combine(const udigit*, lidia_size_t, lidia_size_t)::no reference to crt_table yet");

	if (index< 0 || index > m->num_primes)
		lidia_error_handler("crt", "combine(const udigit*, lidia_size_t, lidia_size_t)::index out of range");

	if (used && l != size)
		lidia_error_handler("crt", "combine(const udigit*, lidia_size_t, lidia_size_t)::incompatible length of pointer object");

	if (! test_mode (POINTER))
		lidia_error_handler("crt", "combine(const udigit*, lidia_size_t, lidia_size_t)::object may not be used with pointers");

	lidia_size_t i;
	bigint tmp;

	if (! used)
		init_pointer (l);

	if (ix[index] == 0) {
		for (i = 0; i < l; i++) {
			multiply(tmp, m->coeff_mod_p[index], small_vec[i]);
			add     (Hv[i], Hv[i], tmp);

			dv[i] += static_cast<double>(small_vec[i]) * m->x[index];
		}

		ix[index]++;
		used++;
	}
	else {
		lidia_error_handler("crt", "combine(const udigit*, lidia_size_t, lidia_size_t)::prime[lidia_size_t] used twice");
	}
}



void
crt::combine (const base_vector< sdigit > & small_vec, lidia_size_t index)
{
	debug_handler("crt", "combine(const base_vector< sdigit > &, lidia_size_t)");

	if (m == NULL)
		lidia_error_handler("crt", "combine(const base_vector< sdigit > &, lidia_size_t)::no reference to crt_table yet");

	if (index< 0 || index > m->num_primes)
		lidia_error_handler("crt", "combine(const base_vector< sdigit > &, lidia_size_t)::index out of range");

	if (used && small_vec.size() != Hvec.size())
		lidia_error_handler("crt", "combine(const base_vector< sdigit > &, lidia_size_t)::incompatible length of vector");

	if (! test_mode (VECTOR))
		lidia_error_handler("crt", "combine(const base_vector< sdigit > &, lidia_size_t)::object may not be used with vectors");

	lidia_size_t i;
	lidia_size_t l;
	bigint tmp;


	l = small_vec.size();

	if (! used)
		init_vector (l);

#if 0
	base_vector< bigint > & Hptr = Hvec;
	base_vector< double > & dptr = dvec;
	const base_vector< sdigit > & sptr = small_vec;
#else
	bigint * Hptr = Hvec.get_data_address();
	double * dptr = dvec.get_data_address();
	sdigit * sptr = small_vec.get_data_address();
#endif


	if (ix[index] == 0) {
		for (i = 0; i < l; i++) {
			multiply(tmp, m->coeff_mod_p[index], sptr[i]);
			add     (Hptr[i], Hptr[i], tmp);

			dptr[i] += static_cast<double>(sptr[i]) * m->x[index];
		}

		ix[index]++;
		used++;
	}
	else {
		lidia_error_handler("crt", "combine(const base_vector< sdigit > &, lidia_size_t, lidia_size_t)::prime[lidia_size_t] used twice");
	}
}



void
crt::combine (const base_matrix< sdigit > & small_mat, lidia_size_t index)
{
	debug_handler("crt", "combine(const base_matrix< sdigit > &, lidia_size_t)");

	lidia_size_t srow, scol;
	lidia_size_t Hrow, Hcol;

	srow = small_mat.get_no_of_rows   ();
	scol = small_mat.get_no_of_columns ();
	Hrow = Hmat.get_no_of_rows   ();
	Hcol = Hmat.get_no_of_columns ();

	if (m == NULL)
		lidia_error_handler("crt", "combine(const base_matrix< sdigit > &, lidia_size_t)::no reference to crt_table yet");

	if (index< 0 || index > m->num_primes)
		lidia_error_handler("crt", "combine(const base_matrix< sdigit > &, lidia_size_t)::index out of range");

	if (used && (srow != Hrow || scol != Hcol))
		lidia_error_handler("crt", "combine(const base_matrix< sdigit > &, lidia_size_t)::incompatible size of matrix");

	if (! test_mode (MATRIX))
		lidia_error_handler("crt", "combine(const base_matrix< sdigit > &, lidia_size_t)::object may not be used with matrices");

	lidia_size_t i, j;
	sdigit stmp;
	double dtmp;
	bigint btmp, tmp;

	if (! used) {
		init_matrix (srow, scol);
	}

	if (ix[index] == 0) {
		for (i = 0; i < srow; i++) {
			for (j = 0; j < scol; j++) {
				stmp = small_mat(i, j);
				btmp = Hmat(i, j);
				dtmp = dmat(i, j);

				multiply(tmp, m->coeff_mod_p[index], stmp);
				add     (btmp, btmp, tmp);

				dtmp += static_cast<double>(stmp) * m->x[index];

				Hmat.sto (i, j, btmp);
				dmat.sto (i, j, dtmp);
			}
		}

		ix[index]++;
		used++;
	}
	else {
		lidia_error_handler("crt", "combine(const base_matrix< sdigit > &, lidia_size_t, lidia_size_t)::prime[lidia_size_t] used twice");
	}
}



//***************************************************************
//
//		get_result for INTEGER mode
//
//***************************************************************

void
crt::get_result(bigint & big)
{
	debug_handler("crt", "get_result(bigint&)");

	if (m == NULL)
		lidia_error_handler("crt", "get_result(bigint&)::no reference to crt_table yet");

	if (m->mode != INTEGER)
		lidia_error_handler("crt", "get_result(bigint&)::table not in INTEGER mode");

	if (! test_mode (SINGLE))
		lidia_error_handler("crt", "get_result(bigint&)::object may not be used with pointers");

	lidia_size_t    i;
	bigint prod, prod_over_2;

	if (used == m->num_primes) {
		big.assign(d + 0.5); //assignment bigint < - double
		multiply (big, big, m->minus_M_mod_p);
		add      (big, big, H);

		m->NORMALIZE (big, m->M_over_2, m->M);
	}
	else {
		//reduce modulo product of used primes

		warning_handler("crt", "get_result(bigint&)::not all primes have been used");

		prod.assign_one();
		for (i = 0; i < m->num_primes; i++) {
			if (ix[i] == 1)
				multiply (prod, prod, m->primes[i]);
		}

		add (prod_over_2 , prod , 1UL);
		shift_right (prod_over_2, prod_over_2 , 1);

		remainder (big, H, prod);

		m->NORMALIZE (big, prod_over_2, prod);
	}

	reset();
}



void
crt::get_result(bigint *& big_vec, lidia_size_t & l)
{
	debug_handler("crt", "get_result(bigint*&, lidia_size_t &)");

	if (m == NULL)
		lidia_error_handler("crt", "get_result(bigint*& , lidia_size_t&)::no reference to crt_table yet");

	if (m->mode != INTEGER)
		lidia_error_handler("crt", "get_result(bigint*& , lidia_size_t&)::table not in INTEGER mode");

	if (! test_mode (POINTER))
		lidia_error_handler("crt", "get_result(bigint*& , lidia_size_t&)::object may not be used with pointers");

	lidia_size_t    i;
	bigint prod, prod_over_2;

	if (l != size) {
		// reallocate memory for big_vec

		if (l > 0) delete [] big_vec;

		l = size;
		big_vec = new bigint [size];

		if (big_vec == NULL)
			lidia_error_handler("crt", "get_result(bigint*& , lidia_size_t&)::out of memory (big_vec)");
	}

	if (used == m->num_primes) {
		for (i = 0; i < l; i++) {
			big_vec[i] = (dv[i] + 0.5); //assignment bigint < - double
			multiply (big_vec[i], big_vec[i], m->minus_M_mod_p);
			add      (big_vec[i], big_vec[i], Hv[i]);

			m->NORMALIZE (big_vec[i], m->M_over_2, m->M);
		}
	}
	else {
		//reduce modulo product of used primes

		warning_handler("crt", "get_result(bigint*, lidia_size_t)::not all primes have been used");

		prod.assign_one();
		for (i = 0; i < m->num_primes; i++) {
			if (ix[i] == 1)
				multiply (prod, prod, m->primes[i]);
		}

		add (prod_over_2 , prod , 1UL);
		shift_right (prod_over_2, prod_over_2 , 1);

		for (i = 0; i < l; i++) {
			remainder (big_vec[i], Hv[i], prod);
			m->NORMALIZE (big_vec[i], prod_over_2, prod);
		}
	}

	reset();
}



void
crt::get_result(base_vector< bigint > & big_vec)
{
	debug_handler("crt", "get_result(base_vector< bigint > &)");

	if (m == NULL)
		lidia_error_handler("crt", "get_result(base_vector< bigint > &)::no reference to crt_table yet");

	if (m->mode != INTEGER)
		lidia_error_handler("crt", "get_result(base_vector< bigint > &)::table not in INTEGER mode");

	if (! test_mode (VECTOR))
		lidia_error_handler("crt", "get_result(base_vector< bigint > &)::object may not be used with pointers");

	lidia_size_t    i;
	lidia_size_t    l;
	bigint prod, prod_over_2;

	l = Hvec.size();

	big_vec.set_capacity (l);
	big_vec.set_size     (l);

#if 0
	base_vector< bigint > & bptr = big_vec;
	base_vector< bigint > & Hptr = Hvec;
	base_vector< double > & dptr = dvec;
#else
	bigint * bptr = big_vec.get_data_address();
	bigint * Hptr = Hvec.get_data_address();
	double * dptr = dvec.get_data_address();
#endif


	l = Hvec.size();

	if (used == m->num_primes) {
		for (i = 0; i < l; i++) {
			bptr[i] = (dptr[i] + 0.5); //assignment bigint < - double
			multiply (bptr[i], bptr[i], m->minus_M_mod_p);
			add      (bptr[i], bptr[i], Hptr[i]);

			m->NORMALIZE (bptr[i], m->M_over_2, m->M);
		}
	}
	else {
		//reduce modulo product of used primes

		warning_handler("crt", "get_result(base_vector< bigint > &)::not all primes have been used");

		prod.assign_one();
		for (i = 0; i < m->num_primes; i++) {
			if (ix[i] == 1)
				multiply (prod, prod, m->primes[i]);
		}

		add (prod_over_2 , prod , 1UL);
		shift_right (prod_over_2, prod_over_2 , 1);

		for (i = 0; i < l; i++) {
			remainder (bptr[i], Hptr[i], prod);
			m->NORMALIZE (bptr[i], prod_over_2, prod);
		}
	}

	reset();
}



void
crt::get_result(base_matrix< bigint > & big_mat)
{
	debug_handler("crt", "get_result(base_matrix< bigint > &)");

	if (m == NULL)
		lidia_error_handler("crt", "get_result(base_matrix< bigint > &)::no reference to crt_table yet");

	if (m->mode != INTEGER)
		lidia_error_handler("crt", "get_result(base_matrix< bigint > &)::table not in INTEGER mode");

	if (! test_mode (MATRIX))
		lidia_error_handler("crt", "get_result(base_matrix< bigint > &)::object must be used with martrix");

	lidia_size_t   i, j;

	double dtmp;
	bigint Htmp, btmp;
	bigint tmp;
	bigint prod, prod_over_2;

	lidia_size_t Hrow, Hcol;

	Hrow = Hmat.get_no_of_rows   ();
	Hcol = Hmat.get_no_of_columns ();

	big_mat.resize (Hrow, Hcol);

	if (used == m->num_primes) {
		for (i = 0; i < Hrow; i++) {
			for (j = 0; j < Hcol; j++) {
				Htmp = Hmat(i, j);
				dtmp = dmat(i, j);

				tmp = (dtmp + 0.5); //assignment bigint < - double
				multiply (tmp, tmp, m->minus_M_mod_p);
				add      (tmp, tmp, Htmp);

				m->NORMALIZE (tmp, m->M_over_2, m->M);

				big_mat.sto (i, j , tmp);
			}
		}
	}
	else {
		//reduce modulo product of used primes

		warning_handler("crt", "get_result(base_matrix< bigint > &)::not all primes have been used");

		prod.assign_one();
		for (i = 0; i < m->num_primes; i++) {
			if (ix[i] == 1)
				multiply (prod, prod, m->primes[i]);
		}

		add (prod_over_2 , prod , 1UL);
		shift_right (prod_over_2, prod_over_2 , 1);

		for (i = 0; i < Hrow; i++) {
			for (j = 0; j < Hcol; j++) {
				tmp = Hmat(i, j);

				remainder (tmp, tmp, prod);
				m->NORMALIZE (tmp, prod_over_2, prod);

				Hmat.sto (i, j , tmp);
			}
		}
	}

	reset();
}



//***************************************************************
//
//		get_result for MODULAR mode
//
//***************************************************************

void
crt::get_result(bigmod & big)
{
	debug_handler("crt", "get_result(bigmod&)");

	if (m == NULL)
		lidia_error_handler("crt", "get_result(bigmod&)::no reference to crt_table yet");

	if (m->mode != MODULAR)
		lidia_error_handler("crt", "get_result(bigmod&)::table not in MODULAR mode");

	if (! test_mode (SINGLE))
		lidia_error_handler("crt", "get_result(bigmod&)::object may not be used with pointers");

	bigint t;

	if (used == m->num_primes) {
		t.assign(d + 0.5); //assignment bigint < - double
		multiply (t, t, m->minus_M_mod_p);
		add      (t, t, H);

		big.assign (t);
	}
	else {
		lidia_error_handler("crt", "get_result(bigmod&)::not all primes have been used");
	}

	reset();
}



void
crt::get_result (bigmod *& big_vec, lidia_size_t & l)
{
	debug_handler("crt", "get_result(bigmod*, lidia_size_t)");

	if (m == NULL)
		lidia_error_handler("crt", "get_result(bigmod*, lidia_size_t)::no reference to crt_table yet");

	if (m->mode != MODULAR)
		lidia_error_handler("crt", "get_result(bigmod*, lidia_size_t)::table not in MODULAR mode");

	if (! test_mode (POINTER))
		lidia_error_handler("crt", "get_result(bigmod*, lidia_size_t)::object may not be used with pointers");

	lidia_size_t    i;
	bigint t;


	if (l != size) {
		// reallocate memory for big_vec

		if (l > 0) delete [] big_vec;

		l = size;
		big_vec = new bigmod [size];

		if (big_vec == NULL)
			lidia_error_handler("crt", "get_result(bigmod*& , lidia_size_t&)::out of memory (big_vec)");
	}


	if (used == m->num_primes) {
		for (i = 0; i < l; i++) {
			t.assign(dv[i] + 0.5); //assignment bigint < - double
			multiply (t, t, m->minus_M_mod_p);
			add      (t, t, Hv[i]);

			big_vec[i].assign(t); //assignment bigmod < - bigint
		}
	}
	else {
		lidia_error_handler("crt", "get_result(bigmod*, lidia_size_t)::not all primes have been used");
	}

	reset();
}



void
crt::get_result(base_vector< bigmod > & big_vec)
{
	debug_handler("crt", "get_result(base_vector< bigmod > &)");

	if (m == NULL)
		lidia_error_handler("crt", "get_result(base_vector< bigmod > &)::no reference to crt_table yet");

	if (m->mode != MODULAR)
		lidia_error_handler("crt", "get_result(base_vector< bigmod > &)::table not in MODULAR mode");

	if (! test_mode (VECTOR))
		lidia_error_handler("crt", "get_result(base_vector< bigmod > &)::object may not be used with pointers");

	lidia_size_t    i;
	lidia_size_t    l;
	bigint  t;

	l = Hvec.size();

	big_vec.set_capacity (l);
	big_vec.set_size     (l);

#if 0
	base_vector< bigmod > & bptr = big_vec;
	base_vector< bigint > & Hptr = Hvec;
	base_vector< double > & dptr = dvec;
#else
	bigmod * bptr = big_vec.get_data_address();
	bigint * Hptr = Hvec.get_data_address();
	double * dptr = dvec.get_data_address();
#endif


	l = Hvec.size();

	if (used == m->num_primes) {
		for (i = 0; i < l; i++) {
			t.assign(dptr[i] + 0.5); //  assignment bigint < - double
			multiply (t, t, m->minus_M_mod_p);
			add      (t, t, Hptr[i]);

			bptr[i].assign (t);
		}
	}
	else {
		lidia_error_handler("crt", "get_result(base_vector< bigint > &, lidia_size_t)::not all primes have been used");
	}

	reset();
}



void
crt::get_result(base_matrix< bigmod > & big_mat)
{
	debug_handler("crt", "get_result(base_matrix< bigmod > &, lidia_size_t)");

	if (m == NULL)
		lidia_error_handler("crt", "get_result(base_matrix< bigmod > &)::no reference to crt_table yet");

	if (m->mode != MODULAR)
		lidia_error_handler("crt", "get_result(base_matrix< bigmod > &)::table not in MODULAR mode");

	if (! test_mode (MATRIX))
		lidia_error_handler("crt", "get_result(base_matrix< bigmod > &)::object must be used with martrix");

	lidia_size_t   i, j;

	double dtmp;
	bigmod btmp;
	bigint Htmp, tmp;

	lidia_size_t Hrow, Hcol;

	Hrow = Hmat.get_no_of_rows    ();
	Hcol = Hmat.get_no_of_columns ();

	big_mat.resize (Hrow, Hcol);

	if (used == m->num_primes) {
		for (i = 0; i < Hrow; i++) {
			for (j = 0; j < Hcol; j++) {
				Htmp = Hmat(i, j);
				dtmp = dmat(i, j);

				tmp = (dtmp + 0.5); //assignment bigint < - double
				multiply (tmp, tmp, m->minus_M_mod_p);
				add      (tmp, tmp, Htmp);

				big_mat.sto (i, j , tmp);
			}
		}
	}
	else {
		lidia_error_handler("crt", "get_result(base_matrix< bigmod > &)::not all primes have been used");
	}

	reset();
}



lidia_size_t
crt::number_of_primes() const
{
	if (m == NULL)
		lidia_error_handler("crt", "number_of_prime()::no reference to crt_table yet");

	return m->num_primes;
}



lidia_size_t
crt::how_much_to_use (const bigint &B) const
{
	if (m == NULL)
		lidia_error_handler("crt", "how_much_to_use(const bigint&)::no reference to crt_table yet");

	return m->how_much_to_use(B);
}



void
crt::info () const
{
	lidia_size_t i;

	std::cout << "==================crt==================\n" << std::endl;

	if (m == NULL) {
		std::cout << " No connection to crt_table object !!! \n" << std::endl;
	}
	else {
		std::cout << "     Mode                  : ";
		switch (mode) {
		case INIT    :  std::cout << "   INIT   "; break;
		case SINGLE  :  std::cout << "   SINGLE "; break;
		case POINTER :  std::cout << "   POINTER"; break;
		case VECTOR  :  std::cout << "   VECTOR "; break;
		case MATRIX  :  std::cout << "   MATRIX "; break;
		default      :  std::cout << "   unknown"; break;
		}
		std::cout << "\n";
		std::cout << "     Number of used primes :    " << used << std::endl;
		std::cout << "\n";
		std::cout << "     Current storage       : ";
		switch (mode) {
		case INIT    :

			std::cout << "   empty  ";

			break;

		case SINGLE  :

			std::cout << "   " << H << " " << d;

			break;

		case POINTER :

			std::cout << "   { ";
			for (i = 0; i < size; i++)
				std::cout << Hv[i] << " ";
			std::cout << " } " << std::endl;
			std::cout << "                               { ";
			for (i = 0; i < size; i++)
				std::cout << dv[i] << " ";
			std::cout << " } " << std::endl;

			break;

		case VECTOR  :

			std::cout << "   " << Hvec << std::endl;
			std::cout << "                                " << dvec << std::endl;

			break;


		case MATRIX :

			std::cout << "   " << Hmat << std::endl;
			std::cout << "                                " << dmat << std::endl;

			break;

		}
	}

	std::cout << "\n==================***==================\n" << std::endl;

}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
