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


#ifndef LIDIA_HNF_KERNEL_CC_GUARD_
#define LIDIA_HNF_KERNEL_CC_GUARD_


#ifndef LIDIA_INFO_H_GUARD_
# include	"LiDIA/info.h"
#endif
#ifndef LIDIA_HNF_CONF_H_GUARD_
# include	"LiDIA/matrix/hnf_conf.h"
#endif
#ifndef LIDIA_HNF_KERNEL_H_GUARD_
# include	"LiDIA/matrix/hnf_kernel.h"
#endif

#if defined LIDIA_COMPUTE_TIME_INFO || defined LIDIA_COMPUTE_RUNTIME_INFO
# include	<fstream>
#endif
#ifdef LIDIA_COMPUTE_TIME_INFO
# ifndef LIDIA_TIMER_H_GUARD_
#  include	"LiDIA/timer.h"
# endif
# ifndef LIDIA_BIGFLOAT_H_GUARD_
#  include	"LiDIA/bigfloat.h"
# endif
#endif



#ifdef LIDIA_NAMESPACE
# ifndef IN_NAMESPACE_LIDIA
namespace LiDIA {
# endif
#endif



#define DMESSAGE "hnf_kernel"


////////////////////////////////////////////////////////////////////////////////
// HNF computation
////////////////////////////////////////////////////////////////////////////////

//
// mgcd
//

template< class T, class REP, class REP1, class CONF >
inline bool havas_kernel< T, REP, REP1, CONF >::
mgcd(MR< T > &A, lidia_size_t &startr, lidia_size_t &startc, const T &BOUND)
{
	// Init
	conf.init(A, BOUND);
	--startc;
	--startr;

	conf.mgcd(A, startr, startc);
	return true;
}



template< class T, class REP, class REP1, class CONF >
inline bool havas_kernel< T, REP, REP1, CONF >::
mgcd(MR< T > &A, matrix< bigint > &TR, lidia_size_t &startr,
     lidia_size_t &startc, const T &BOUND)
{
	// Init
	conf.init(A, BOUND);
	--startr;
	--startc;

	conf.mgcd(A, TR, startr, startc);
	return true;
}



//
// row elimination
//

template< class T, class REP, class REP1, class CONF >
bool havas_kernel< T, REP, REP1, CONF >::
hnf_Z1(MR< T > &A, lidia_size_t &startr, lidia_size_t &startc, const T &BOUND)
{

	// Init
	conf.init(A, BOUND);

#ifdef LIDIA_COMPUTE_RUNTIME_INFO
	std::ofstream log1_ofs("maximaler_Eintrag.dat");
	std::ofstream log2_ofs("Eintragsdichte.dat");
	bigint Bound = 10;
	power(Bound, Bound, 100);
#endif

#ifdef LIDIA_COMPUTE_TIME_INFO
	std::ofstream log3_ofs("Time_pro_Iteration.dat");
	timer t;
#endif

	// row elimination
	for (--startc, --startr; startr >= 0 && startc >= 0; startc--, startr--) {
		lidia_xinfo_handler("stf", startr, A.rows);

#ifdef LIDIA_COMPUTE_RUNTIME_INFO
		T MAX = 0, Durch = 0;
		lidia_size_t no_of_elements;
		modul.kennwerte(A, MAX, no_of_elements, Durch);

		no_of_elements = 0;
		for (lidia_size_t i = 0; i <= startc; i++)
			no_of_elements += A.value_counter[i];

		log1_ofs << A.rows - startr << " " << MAX << std::endl;
		if (startr >= 0 && startc >= 0)
			log2_ofs << A.rows - startr << " " << bigfloat(100 * no_of_elements)/bigfloat((startr+1) * (startc+1)) << std::endl;

#if 0
		if (MAX > Bound) {
			for (lidia_size_t i = A.rows - startr + 1; i < A.rows; i++) {
				log1_ofs << i << " " << 0 << std::endl;
				log2_ofs << i << " " << 0 << std::endl;
			}
			log1_ofs.close();
			log2_ofs.close();
			return true;
		}
#endif
#endif

#ifdef LIDIA_COMPUTE_TIME_INFO
		t.start_timer();
#endif

		lidia_size_t st = conf.mgcd(A, startr, startc);

		if (st == -1) {
			startr++;
			startc++;
			return false;
		}
		else if (st == 0) // All elements are zero
			startc++;
		else {
			if (modul.member(A, startr, startc) < A.Zero)
				modul.negate_column(A, startc, startr);

			if (!conf.normalize_row(A, startr, A.columns - 1, startc, A.rows - 1)) {
				startr++;
				startc++;
				return false;
			}
		}
#ifdef LIDIA_COMPUTE_TIME_INFO
		t.stop_timer();
		log3_ofs << A.rows - 1 - startr << " " << " " << t.user_time() << std::endl;
#endif
	}
	startr++;
	startc++;

#ifdef LIDIA_COMPUTE_TIME_INFO
	log3_ofs.close();
#endif

#ifdef LIDIA_COMPUTE_RUNTIME_INFO
	log1_ofs.close();
	log2_ofs.close();
#endif

	return true;
}



template< class T, class REP, class REP1, class CONF >
bool havas_kernel< T, REP, REP1, CONF >::
hnf_Z1(MR< T > &A, matrix< bigint > &TR, lidia_size_t &startr,
       lidia_size_t &startc, const T &BOUND)
{
	// Init
	conf.init(A, BOUND);


	// row elimination
	for (--startc, --startr; startr >= 0 && startc >= 0; startc--, startr--) {
		lidia_xinfo_handler("stf", startr, A.rows);

		lidia_size_t st = conf.mgcd(A, TR, startr, startc);

		if (st == -1) {
			startr++;
			startc++;
			return false;
		}
		else if (st == 0) // All elements are zero
			startc++;
		else {
			if (modul.member(A, startr, startc) < A.Zero) {
				modul.negate_column(A, startc, startr);

				for (lidia_size_t j = 0; j < A.columns; j++)
					TR.sto(j, startc, -TR.member(j, startc));
			}

			if (!conf.normalize_row(A, TR, startr, A.columns - 1, startc, A.rows - 1)) {
				startr++;
				startc++;
				return false;
			}
		}
	}
	startr++;
	startc++;
	return true;
}



template< class T, class REP, class REP1, class CONF >
bool havas_kernel< T, REP, REP1, CONF >::
hnf_Z2(MR< T > &A, lidia_size_t &startr, lidia_size_t &startc, const T &BOUND)
{

#ifdef LIDIA_COMPUTE_TIME_INFO
	std::ofstream log3_ofs("Time_pro_Iteration.dat");
	timer t;
#endif

	// Init
	conf.init(A, BOUND);
	for (--startc, --startr; startr >= 0 && startc >= 0; startc--, startr--) {
		lidia_xinfo_handler("stf", startr, A.rows);

#ifdef LIDIA_COMPUTE_TIME_INFO
		t.start_timer();
#endif

		// row elimination
		lidia_size_t st = conf.mgcd(A, startr, startc);

		// stop computation
		if (st == -1) {
			startr++;
			startc++;
			return false;
		}
		else if (st == 1) {
			if (modul.member(A, startr, startc) < A.Zero)
				modul.negate_column(A, startc, startr);
		}

#ifdef LIDIA_COMPUTE_TIME_INFO
		t.stop_timer();
		log3_ofs << A.rows - 1 - startr << " " << " " << t.user_time() << std::endl;
#endif

	}

#ifdef LIDIA_COMPUTE_TIME_INFO
	log3_ofs.close();
#endif

	startr++;
	startc++;
	return true;
}



template< class T, class REP, class REP1, class CONF >
bool havas_kernel< T, REP, REP1, CONF >::
hnf_Z2(MR< T > &A, matrix< bigint > &TR, lidia_size_t &startr, lidia_size_t &startc, const T &BOUND)
{
	// Init
	conf.init(A, BOUND);

	// row elimination
	for (--startc, --startr; startr >= 0 && startc >= 0; startc--, startr--) {
		lidia_xinfo_handler("stf", startr, A.rows);

		lidia_size_t st = conf.mgcd(A, TR, startr, startc);

		if (st == -1) {
			startr++;
			startc++;
			return false;
		}
		else if (st == 1) {
			if (modul.member(A, startr, startc) < A.Zero) {
				modul.negate_column(A, startc, startr);

				for (lidia_size_t j = 0; j < A.columns; j++)
					TR.sto(j, startc, -TR.member(j, startc));
			}
		}
	}
	startr++;
	startc++;
	return true;
}



template< class T, class REP, class REP1, class CONF >
inline bool havas_kernel< T, REP, REP1, CONF >::
stf(MR< T > &A, lidia_size_t &startr, lidia_size_t &startc, const T &BOUND)
{
	// Init
	conf.init(A, BOUND);

	// row elimination
	for (--startc, --startr; startr >= 0 && startc >= 0; startc--, startr--) {
		lidia_xinfo_handler("stf", startr, A.rows);

		lidia_size_t st = conf.mgcd(A, startr, startc);
		if (st == -1) {
			startr++;
			startc++;
			return false;
		}
		else if (st == 0) // All elements are zero
			startc++;
		else {
			if (modul.member(A, startr, startc) < A.Zero)
				modul.negate_column(A, startc, startr);
		}
	}
	startr++;
	startc++;
	return true;
}



template< class T, class REP, class REP1, class CONF >
bool havas_kernel< T, REP, REP1, CONF >::
stf(MR< T > &A, matrix< bigint > &TR, lidia_size_t &startr,
    lidia_size_t &startc, const T &BOUND)
{
	// Init
	conf.init(A, BOUND);

	// row elimination
	for (--startc, --startr; startr >= 0 && startc >= 0; startc--, startr--) {
		lidia_xinfo_handler("stf", startr, A.rows);
		lidia_size_t st = conf.mgcd(A, TR, startr, startc);

		if (st == -1) {
			startr++;
			startc++;
			return false;
		}
		else if (st == 0) // All elements are zero
			startc++;
		else {
			if (modul.member(A, startr, startc) < A.Zero) {
				modul.negate_column(A, startc, startr);
				modul_TR.negate_column(TR, startc, TR.get_no_of_rows() - 1);
			}
		}
	}
	startr++;
	startc++;
	return true;
}



#if 0
template< class T, class REP, class REP1, class CONF >
bool havas_kernel< T, REP, REP1, CONF >::
hnf_Z2_mod(MR< T > &A, lidia_size_t &startr, lidia_size_t &startc, const T &DET)
{

#ifdef LIDIA_COMPUTE_TIME_INFO
	std::ofstream log3_ofs("Time_pro_Iteration.dat");
	timer t;
#endif

	// Init
	//conf.init(A, BOUND);
	for (--startc, --startr; startr >= 0 && startc >= 0; startc--, startr--) {
		lidia_xinfo_handler("stf", startr, A.rows);

#ifdef LIDIA_COMPUTE_TIME_INFO
		t.start_timer();
#endif

		// row elimination
		lidia_size_t st = conf.mgcd(A, startr, startc);
		modul.remainder(A, DET);

		// stop computation
		if (st == -1) {
			startr++;
			startc++;
			return false;
		}
		else if (st == 1) {
			if (modul.member(A, startr, startc) < A.Zero)
				modul.negate_column_mod(A, startc, startr, DET);
		}

#ifdef LIDIA_COMPUTE_TIME_INFO
		t.stop_timer();
		log3_ofs << A.rows - 1 - startr << " " << " " << t.user_time() << std::endl;
#endif

	}

#ifdef LIDIA_COMPUTE_TIME_INFO
	log3_ofs.close();
#endif

	startr++;
	startc++;
	return true;
}
#endif



// ADDED BY MJJ

template< class T, class REP, class REP1, class CONF >
int havas_kernel< T, REP, REP1, CONF >::
hnf_Z2(MR< T > &A, trans_matrix & TR, matrix< bigint > & tran,
       lidia_size_t &startr, lidia_size_t &startc, lidia_size_t rowsorig,
       const T &BOUND)
{
	lidia_size_t n = 2;
	lidia_size_t dense_stop;

	if ((startr <= 500) &&
	    (tran.get_representation() == matrix_flags::sparse_representation))
		return -1;

	if (rowsorig >= 4000)
		dense_stop = 1400;
	else if (rowsorig >= 3000)
		dense_stop = 800;
	else if (rowsorig >= 2000)
		dense_stop = 600;
	else if (rowsorig >= 800)
		dense_stop = 500;
	else
		dense_stop = -1;


	tran.reset();
	tran.resize(startc, startc);
	tran.diag(1, 0);

	conf.init(A, BOUND);

	// Elimination
	for (--startc, --startr; startr >= 0 && startc >= 0; startc--, startr--) {
		lidia_xinfo_handler("stf", startr, A.rows);

		lidia_size_t st = conf.mgcd(A, tran, startr, startc);
		if (st == -1) {
			startr++;
			startc++;
			return 0;
		}
		else if (st == 0) {
			if (modul.member(A, startr, startc) < A.Zero) {
				modul.negate_column(A, startc, startr);
				modul_TR.negate_column(tran, startc, tran.get_no_of_rows() - 1);
			}
		}

		if (tran.get_representation() == matrix_flags::sparse_representation) {
			if (startr <= dense_stop) {
				startr++;
				startc++;
				return -1;
			}
		}
		else if ((A.bitfield.get_representation() == matrix_flags::sparse_representation) &&
			 ((startr <= (dense_stop-200)))) {
			startr++;
			startc++;
			return 0;
		}

		if (startr == (rowsorig/n)) {
			TR.store_matrix(tran);

			tran.resize(startc+1, startc+1);
			tran.diag(1, 0);

			++n;
		}
	}

	startr++;
	startc++;
	return 1;
}



////////////////////////////////////////////////////////////////////////////////
// kannan_kernel
////////////////////////////////////////////////////////////////////////////////

//
// column step form
//

template< class T, class REP, class REP1, class CONF >
bool kannan_kernel< T, REP, REP1, CONF >::
hnf(MR< T > &A, lidia_size_t &startr, lidia_size_t &startc, const T &BOUND)
{
	normalization_kernel< T, REP, REP > normalize_modul;

#ifdef LIDIA_COMPUTE_RUNTIME_INFO
	std::ofstream log1_ofs("maximaler_Eintrag.dat");
	std::ofstream log2_ofs("Eintragsdichte.dat");
	bigint Bound = 10;
	power(Bound, Bound, 100);
#endif

#ifdef LIDIA_COMPUTE_TIME_INFO
	std::ofstream log3_ofs("Time_pro_Iteration.dat");
	timer t;
#endif

	lidia_size_t pos1 = 0, pos2 = 0;

	T TMP;
	T RES0, RES1, RES2, TMP1, x, y;
	lidia_size_t i, j, l, k;

	for (i = startc - 2, k = startr - 2; i >= 0; i--, k--) {

#ifdef LIDIA_COMPUTE_RUNTIME_INFO
		T MAX = 0, Durch = 0;
		lidia_size_t no_of_elements;
		modul.kennwerte(A, MAX, no_of_elements, Durch);

		no_of_elements = 0;
		for (lidia_size_t f = 0; f <= startc; f++)
			no_of_elements += A.value_counter[f];

		log1_ofs << A.rows - k << " " << MAX << std::endl;
		if (i >= 0 && k >= 0)
			log2_ofs << A.rows - k << " " << bigfloat(100 * no_of_elements)/bigfloat((k+1) * (i+1)) << std::endl;
#endif

#ifdef LIDIA_COMPUTE_TIME_INFO
		t.start_timer();
#endif

		lidia_xinfo_handler("stf", k, A.rows);
		for (j = startr - 1, l = startc - 1; j >= 0 && l > i; j--, l--) {
			if (modul.member(A, j, i) != A.Zero) {
				if (modul.member(A, j, l) == A.Zero)
					modul.swap_columns(A, i, l);
				else {
					RES2 = xgcd(RES0, RES1, modul.member(A, j, i), modul.member(A, j, l));
					x = modul.member(A, j, l) / RES2;
					y = modul.member(A, j , i) / RES2;

					for (lidia_size_t len = 0; len <= j; len++) {
						TMP = modul.member(A, len, i);
						TMP1 = modul.member(A, len, l);

						modul.sto(A, len, i, TMP*x - TMP1*y);
						modul.sto(A, len, l, TMP*RES0 + TMP1*RES1);
					}

				}
			}
		}
		if (k >= 0 && i >= 0) {
			pos1 = k;
			pos2 = i;
			if (modul.member(A, pos1, pos2) < A.Zero) {
				for (lidia_size_t len = 0; len <= pos1; len++)
					modul.sto(A, len, pos2, -modul.member(A, len, pos2));
			}
		}
		conf_modul.normalize(A, pos1, pos2);

#ifdef LIDIA_COMPUTE_TIME_INFO
		t.stop_timer();
		log3_ofs << A.rows - 1 - k << " " << " " << t.user_time() << std::endl;
#endif

	}

#ifdef LIDIA_COMPUTE_TIME_INFO
	log3_ofs.close();
#endif

#ifdef LIDIA_COMPUTE_RUNTIME_INFO
	log1_ofs.close();
	log2_ofs.close();
#endif

	return true;
}



#undef DMESSAGE



#ifdef LIDIA_NAMESPACE
# ifndef IN_NAMESPACE_LIDIA
}	// end of namespace LiDIA
# endif
#endif



#endif	// LIDIA_HNF_KERNEL_CC_GUARD_
