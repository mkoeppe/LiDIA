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
// using `only` bigints
//
// Buchmann Kessler
//
// * Tr_buchmann_kessler(T)
//


//
// the modified lll-algorithm
//
bigint*
bigint_lattice::TrD_mlll_bfl (dense_alg< bigint > &da, lattice_info &li)
{
	debug_handler("bigint_lattice", "TrD_mlll_bfl(da, li, v)");
	lidia_size_t j, k = da.b.columns, l, m;
	sdigit prec;
	sdigit old_prec;

	bigint tempbin0;
	bigint *tempvectbin0;
	bigint *tempvectbin1;
	bigfloat *tempvectbfl0;
	bigfloat *tempvectbfl1;
	bigfloat *B;
	bigfloat halb(0.5);
	bigfloat y_bfl(da.b.y);
	bigfloat tempbfl0, tempbfl1, tempbfl2, Mu, Bz;
	bigfloat *mydel;
	bigfloat **my, **bs;
	bigint *Trdel;
	bigint **Tr;
	bigint *v;
	p_vector< bigfloat > vector_bfl;
	p_vector< bigint > vector_bin;

//
// Allocating needed storage
//
// Lattices
//
	Tr = new bigint*[da.b.rows];
	memory_handler(T, "bigint_lattice", "Tr_mlll_bfl(v) :: "
		       "not enough memory !");
	Trdel = new bigint[da.b.rows * (da.b.rows + 1) + da.b.columns];
	memory_handler(Trdel, "bigint_lattice", "Tr_mlll_bfl(v) :: "
                       "not enough memory !");

	my = new bigfloat*[2 * da.b.rows];
	memory_handler(my, "bigint_lattice", "Tr_mlll_bfl(v) :: "
		       "not enough memory !");
	mydel = new bigfloat[da.b.rows * ((da.b.rows + da.b.columns) + 3)];
	memory_handler(mydel, "bigint_lattice", "Tr_mlll_bfl(v) :: "
                       "not enough memory !");
	bs = &my[da.b.rows];

//
// Restore and set precision
//
	old_prec = bigfloat::get_precision();
	prec = compute_read_precision(da);
	halb.set_precision(da.b.rows * prec);
	li.lll.reduction_steps = 0;
	li.lll.correction_steps = 0;
	li.lll.swaps = 0;
	for (lidia_size_t i = 0; i < da.b.rows; i++) {
		Tr[i] = &Trdel[da.b.rows * i];
		Tr[i][i].assign_one();
		my[i] = &mydel[da.b.rows * i];
		bs[i] = &mydel[da.b.rows * da.b.rows + i * da.b.columns];
	}


//
// Vectors
// Setting pointers
//
	tempvectbfl0 = &mydel[(da.b.rows + da.b.columns) * da.b.rows];
	tempvectbfl1 = &tempvectbfl0[da.b.rows];
	B = &tempvectbfl0[2 * da.b.rows];
	tempvectbin0 = &Trdel[da.b.rows * da.b.rows];
	tempvectbin1 = &tempvectbin0[da.b.rows];

	vector_bin.vectsize = da.b.columns;
	vector_bfl.vectsize = da.b.columns;

	for (lidia_size_t i = 0; i < k + 1; i++) {
		for (lidia_size_t h = 0; h < i; h++) {
			my[h][i].assign_zero();
			for (lidia_size_t cx = 0; cx < da.b.columns; cx++) {
				tempbfl2.assign(da.s.value[i][cx]);
				LiDIA::multiply(tempbfl2, tempbfl2, bs[h][cx]);
				LiDIA::add(my[h][i], my[h][i], tempbfl2);
			}
//
//          for mixed datatypes
//          vector_bfl.scalprod(my[h][i], bs[h], da.s.value[i]);
//
			LiDIA::divide(my[h][i], my[h][i], B[h]);
		}
		vector_bfl.assign_zero(tempvectbfl1);
		for (lidia_size_t h = 0; h < i; h++) {
			vector_bfl.scalmul(tempvectbfl0, my[h][i], bs[h]);
			vector_bfl.add(tempvectbfl1, tempvectbfl1, tempvectbfl0);
		}
		for (lidia_size_t cx = 0; cx < da.b.columns; cx++) {
			tempbfl2.assign(da.s.value[i][cx]);
			LiDIA::subtract(bs[i][cx], tempbfl2, tempvectbfl1[cx]);
		}
//
//      Mixed datatypes
//      vector_bfl.subtract(bs[i], da.s.value[i], tempvectbfl1);
//
		vector_bfl.scalprod(B[i], bs[i], bs[i]);
	}

	m = 1;
	l = 0;

	while (true) {
//
// Reduction Step
//
		if (abs(my[l][m]).compare(halb) > 0) {
			li.lll.reduction_steps++;
			round(tempbfl0, my[l][m]);
			tempbfl0.bigintify(tempbin0);
			vector_bin.scalmul(tempvectbin0, tempbin0, da.s.value[l]);
			vector_bin.subtract(da.s.value[m], da.s.value[m], tempvectbin0);
			vector_bin.vectsize = da.b.rows;
			vector_bin.scalmul(tempvectbin0, tempbin0, Tr[l]);
			vector_bin.subtract(Tr[m], Tr[m], tempvectbin0);
			vector_bin.vectsize = da.b.columns;

			tempbfl1.assign(tempbin0);
			LiDIA::subtract(my[l][m], my[l][m], tempbfl1);
			for (lidia_size_t h = 0; h < l; h++) {
				LiDIA::multiply(tempbfl0, tempbfl1, my[h][l]);
				LiDIA::subtract(my[h][m], my[h][m], tempbfl0);
			}
		}


		j = 0;
		while ((j < k) && (da.s.value[m][j].is_zero()))
			j++;

//
// Exiting - condition
//
		if (j == k) {
			for (lidia_size_t i = m; i < k; i++) {
				vector_bin.assign(da.s.value[i], da.s.value[i + 1]);
			}
			vector_bin.assign_zero(da.s.value[da.b.rows - 1]);
			vector_bin.vectsize = da.b.rows;
//
// Alloc storage for vector of relations
//
			v = new bigint[da.b.rows];
			vector_bin.assign(v, Tr[m]);
			vector_bin.vectsize = da.b.columns;
//
// Free allocated storage
// and restore precision
//
			da.b.rank = da.b.rows - 1;
			bigfloat::set_precision(old_prec);
			delete[] Trdel;
			delete[] Tr;
			delete[] mydel;
			delete[] my;
			return (v);
		}


		if (l >= m - 1) {
			LiDIA::square(tempbfl1, my[m - 1][m]);
			LiDIA::subtract(tempbfl0, y_bfl, tempbfl1);
			LiDIA::multiply(tempbfl0, tempbfl0, B[m - 1]);

			if (B[m].compare(tempbfl0) < 0) {
				Mu.assign(my[m - 1][m]);
				LiDIA::square(tempbfl0, Mu);
				LiDIA::multiply(tempbfl0, tempbfl0, B[m - 1]);
				LiDIA::add(Bz, B[m], tempbfl0);
				if (!Bz.is_approx_zero()) {
					LiDIA::divide(tempbfl0, B[m - 1], Bz);
					LiDIA::multiply(my[m - 1][m], Mu, tempbfl0);
					LiDIA::multiply(B[m], B[m], tempbfl0);
					for (lidia_size_t i = m + 1; i < k + 1; i++) {
						LiDIA::multiply(tempbfl0, Mu, my[m][i]);
						LiDIA::subtract(tempbfl1, my[m - 1][i], tempbfl0);
						LiDIA::multiply(tempbfl0, tempbfl1, my[m - 1][m]);
						LiDIA::add(my[m - 1][i], my[m][i], tempbfl0);
						my[m][i].assign(tempbfl1);
					}
				}
				B[m - 1].assign(Bz);
				li.lll.swaps++;
				vector_bin.swap(da.s.value[m - 1], da.s.value[m]);
				vector_bin.swap(Tr[m - 1], Tr[m]);
				for (j = 0; j <= m - 2; j++) {
					tempbfl1.assign(my[j][m - 1]);
					my[j][m - 1].assign(my[j][m]);
					my[j][m].assign(tempbfl1);
				}
				if (m > 1) {
					m--;
				}
				l = m - 1;
				continue;
			}
		}
		l--;
		if (l < 0) {
			m++;
			l = m - 1;
		}
	}
	return (NULL); // Impossible
}



//
// Buchmann - Kessler version for generating systems
//
void
bigint_lattice::TrD_buchmann_kessler (dense_alg< bigint > & da, lattice_info& li)
{
	debug_handler("bigint_lattice", "TrD_buchmann_kessler(da, li) [1]");
	bigint *help;
	bigint *rel;
	bigint *tempPbin;
	bigint **oldvalue = da.s.value;
	bigint *olddelvalue = da.s.delvalue;
	bigint **tempPPbin;
	bigint *deltempPPbin;
	bigfloat vor, rechts;
	bigfloat zweipotq, alpha;
	bigfloat norm1, norm2;
	bigint bi_norm1, bi_norm2;
	bigint zwpq;
	sdigit n2 = da.b.rows;
	sdigit prec;
	p_vector< bigint > vector;
	lll_kernel_fu< bigint, bigfloat, vector_op < bigint, bigfloat >,
                basis_modules< bigint, bigfloat, Normal > > alg;

//
// Compute bigint approximation of lattice
// to use the schnorr - euchner version of lll
//
	prec = compute_precision(da);
	bigfloat::set_precision(prec);
	alpha_compute(da, alpha);
	zwei_pot_q_compute(da, zweipotq, n2, alpha);
	zweipotq.bigintify(zwpq);

//
// Allocate memory
//
	da.s.value = new bigint*[da.b.rows];
	memory_handler(da.s.value, "bigint_lattice", "Tr_buchmann_kessler(da, li) :: "
		       "not enough memory !");
	da.s.delvalue = new bigint[(da.b.rows + da.b.columns) * da.b.rows];
	memory_handler(da.s.delvalue, "bigint_lattice", "Tr_buchmann_kessler(da, li) :: "
		       "not enough memory !");
	tempPPbin = new bigint*[da.b.rows];
	memory_handler(tempPPbin, "bigint_lattice", "Tr_buchmann_kessler(da, li) :: "
		       "not enough memory !");
	deltempPPbin = new bigint[da.b.rows * da.b.columns + 2*da.b.columns + da.b.rows];
	memory_handler(deltempPPbin, "bigint_lattice", "Tr_buchmann_kessler(da, li) :: "
		       "not enough memory !");
//
// Set pointer
//
	for (lidia_size_t i = 0; i < da.b.rows; i++) {
		da.s.value[i] = &da.s.delvalue[i * (da.b.columns + da.b.rows)];
		tempPPbin[i] = &deltempPPbin[i * da.b.columns];
	}
	help = &deltempPPbin[da.b.rows * da.b.columns];
	rel = &help[da.b.columns];
	tempPbin = &rel[da.b.rows];

//
// Prepare matrix for buchmann  - kessler
//
	for (lidia_size_t i = 0; i < da.b.rows; ++i) {
		da.s.value[i][i + da.b.columns].assign_one(); // 1, wenn i = j, 0 sonst
		for (lidia_size_t j = 0; j < da.b.columns; ++j) {
			LiDIA::multiply(da.s.value[i][j], zwpq, oldvalue[i][j]);
			tempPPbin[i][j].assign(da.s.value[i][j]); // merke Kopie
		}
	}

//
// Compute needed Precision for approximation bigfloats
//
	da.b.columns += da.b.rows;
	prec = compute_read_precision(da);
	da.d.bit_prec = static_cast<sdigit>(static_cast<double>(prec) * (std::log(10.0)/std::log(2.0)));
	da.d.bit_prec = ((da.d.bit_prec/200) + 1) * DOUBLE_MANTISSA_BITS;

//
// Perform lll
//
	debug_handler("bigint_lattice", "Tr_buchmann_kessler(da, li) [2]");
	alg.lll(da, li);
	debug_handler("bigint_lattice", "Tr_buchmann_kessler(da, li) [3]");
	da.b.columns -= da.b.rows;

//
// Check rank of the lattice
//
	lidia_size_t l = 0;
	do {
		vector.vectsize = da.b.columns;
		vector.assign_zero(help); // Initializes help with the zero - vector
		for (lidia_size_t j = 0; j < da.b.rows; ++j) {
			rel[j].assign(da.s.value[l][j + da.b.columns]);
			vector.scalmul(tempPbin, rel[j], tempPPbin[j]);
			vector.add(help, help, tempPbin);
		}
		vector.scalprod(bi_norm2, help, help);
		norm2.assign(bi_norm2);
		vector.vectsize = da.b.rows;
//
// vector.l1_norm_bin(bi_norm1,rel);
//
		bi_norm1.assign_zero();
		for (lidia_size_t j = 0; j < da.b.rows; j++) {
			LiDIA::add(bi_norm1, bi_norm1, abs(rel[j]));
		}

		norm1.assign(bi_norm1);
		sqrt(norm1, norm1);
		++l;
		LiDIA::divide(vor, n2, da.b.rows);
		vor.multiply_by_2();
		sqrt(vor, vor);
		LiDIA::multiply(vor, vor, norm1);
		LiDIA::divide(rechts, bigfloat(std::sqrt(static_cast<double>(n2))), bigfloat(2.0));
		LiDIA::multiply(rechts, rechts, norm1);
	} while ((zweipotq.compare(vor) > 0) &&
		 (norm2.compare(rechts) <= 0) && (l < da.b.rows));
	debug_handler("bigint_lattice", "Tr_buchmann_kessler(da, li) [4]");
	if (l >= da.b.rows) {
		warning_handler("bigint_lattice", "Tr_buchmann_kessler(da, li) :: "
				"lattice of dimension 1");
	}

	da.b.rank = da.b.rows - l + 1; // da.b.rows is the dimension of the lattice
	li.lll.rank = da.b.rank;
//
// Store transformation lattice in math_matrix< bigint >
//
//  for (lidia_size_t i = 0; i < da.b.rows; i++)
//    for (lidia_size_t j = da.b.columns; j < da.b.columns + da.b.rows - da.b.rank; j++)
//      da.s.value[i][j].assign_zero();
	for (lidia_size_t i = 0; i < da.b.rows - da.b.rank; i++) {
		for (lidia_size_t j = 0; j < da.b.rows - i - 1; j++) {
			vector.swap(da.s.value[j], da.s.value[j + 1]);
		}
	}

//
// Restore original matrix !
//
	for (lidia_size_t i = 0; i < da.b.rows; i++) {
		for (lidia_size_t j = 0; j < da.b.columns; j++) {
			LiDIA::swap(da.s.value[i][j], oldvalue[i][j]);
		}
	}

//
// Free allocated storage
//
	delete[] oldvalue;
	delete[] olddelvalue;
	delete[] tempPPbin;
	delete[] deltempPPbin;
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
