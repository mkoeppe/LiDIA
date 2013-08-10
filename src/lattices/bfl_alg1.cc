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



//
// the modified lll-algorithm
//
bigint *bigfloat_lattice::TrD_mlll_bfl(dense_alg< bigfloat > & da,
				       lattice_info& li)
{
	debug_handler("bigfloat_lattice", "Tr_mlll_bfl(da, li)");
	lidia_size_t j, k = columns, l, m;
	sdigit old_prec;
	sdigit prec;
	bigint tempbin0;
	bigint *tempvectbin0;
	bigint *Trdel;
	bigint **Tr;
	bigint *v;
	bigfloat tempbfl0;
	bigfloat tempbfl1;
	bigfloat halb(0.5);
	bigfloat y_bfl(da.b.y);
	bigfloat Mu, Bz;
	bigfloat *mydel;
	bigfloat **my, **bs;
	bigfloat *tempvectbfl0;
	bigfloat *tempvectbfl1;
	bigfloat *tempvectbfl2;
	bigfloat *B;
	p_vector< bigint > vector_bin;
	p_vector< bigfloat > vector_bfl;

	//
	// Allocating needed storage
	//
	// Lattices
	//

	Tr = new bigint*[da.b.rows];
	memory_handler(T, "bigfloat_lattice", "Tr_mlll_bfl(da, li) :: "
		       "not enough memory !");
	my = new bigfloat*[2*da.b.rows];
	memory_handler(my, "bigfloat_lattice", "Tr_mlll_bfl(da, li) :: "
		       "not enough memory !");

	bs = &my[da.b.rows];
	mydel = new bigfloat[da.b.rows*(da.b.rows+da.b.columns)+
			    3*da.b.rows+da.b.columns];
	memory_handler(mydel, "bigfloat_lattice", "Tr_mlll_bfl(da, li) :: "
		       "not enough memory !");
	Trdel = new bigint[da.b.rows*da.b.rows+da.b.rows];
	memory_handler(Trdel, "bigfloat_lattice", "Tr_mlll_bfl(da, li) :: "
		       "not enough memory !");

	//
	// Restore and set precision
	//
	//
	old_prec = bigfloat::get_precision();
	prec = compute_read_precision(da);
	bigfloat::set_precision(da.b.rows*prec);
	li.lll.reduction_steps = 0;
	li.lll.correction_steps = 0;
	li.lll.swaps = 0;
	for (lidia_size_t i = 0; i < da.b.rows; i++) {
		Tr[i] = &Trdel[i*da.b.rows];
		Tr[i][i].assign_one();
		my[i] = &mydel[i*da.b.rows];
		bs[i] = &mydel[da.b.rows*da.b.rows+i*da.b.columns];
	}


	//
	// Vectors
	// Setting pointers
	//
	tempvectbfl1 = &mydel[da.b.rows*(da.b.rows+da.b.columns)];
	tempvectbfl2 = &tempvectbfl1[da.b.rows];
	B = &tempvectbfl1[2*da.b.rows];

	tempvectbin0 = &Trdel[da.b.rows*da.b.rows];
	tempvectbfl0 = &tempvectbfl1[3*da.b.rows];

	vector_bin.vectsize = da.b.columns;
	vector_bfl.vectsize = da.b.columns;


	for (lidia_size_t i = 0; i < k+1; i++) {
		for (lidia_size_t h = 0; h < i; h++) {
			my[h][i].assign_zero();
			vector_bfl.scalprod(my[h][i], bs[h], value[i]);
			LiDIA::divide(my[h][i], my[h][i], B[h]);
		}
		vector_bfl.assign_zero(tempvectbfl2);
		for (lidia_size_t h = 0; h < i; h++) {
			vector_bfl.scalmul(tempvectbfl1, my[h][i], bs[h]);
			vector_bfl.add(tempvectbfl2, tempvectbfl2, tempvectbfl1);
		}
		vector_bfl.subtract(bs[i], value[i], tempvectbfl2);
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
			LiDIA::round(tempbfl0, my[l][m]);
			tempbfl0.bigintify(tempbin0);
			vector_bfl.scalmul(tempvectbfl0, tempbfl0, value[l]);
			vector_bfl.subtract(value[m], value[m], tempvectbfl0);
			vector_bin.vectsize = da.b.rows;
			vector_bin.scalmul(tempvectbin0, tempbin0, Tr[l]);
			vector_bin.subtract(Tr[m], Tr[m], tempvectbin0);
			vector_bin.vectsize = da.b.columns;

			LiDIA::subtract(my[l][m], my[l][m], tempbfl0);
			for (lidia_size_t h = 0; h < l; h++) {
				LiDIA::multiply(tempbfl1, tempbfl0, my[h][l]);
				LiDIA::subtract(my[h][m], my[h][m], tempbfl1);
			}
		}


		j = 0;
		while ((j < k) && (value[m][j].is_approx_zero())) {
			j++;
		}

		//
		// Exiting - condition
		//
		if (j == k) {
			for (lidia_size_t i = m; i < k; i++)
				vector_bfl.assign(value[i], value[i+1]);
			vector_bfl.assign_zero(value[rows-1]);
			vector_bin.vectsize = da.b.rows;
			//
			// Alloc storage for vector of relations
			//
			v = new bigint[da.b.rows];
			vector_bin.assign(v, Tr[m]);
			vector_bin.vectsize = da.b.columns;
			bigfloat::set_precision(old_prec);
			delete[] mydel;
			delete[] my;
			delete[] Trdel;
			delete[] Tr;
			return (v);
		}


		if (l >= m-1) {
			LiDIA::square(tempbfl1, my[m-1][m]);
			LiDIA::subtract(tempbfl0, y_bfl, tempbfl1);
			LiDIA::multiply(tempbfl0, tempbfl0, B[m-1]);

			if (B[m].compare(tempbfl0) < 0) {
				Mu.assign(my[m-1][m]);
				LiDIA::square(tempbfl0, Mu);
				LiDIA::multiply(tempbfl0, tempbfl0, B[m-1]);
				LiDIA::add(Bz, B[m], tempbfl0);
				if (!Bz.is_approx_zero()) {
					LiDIA::divide(tempbfl0, B[m-1], Bz);
					LiDIA::multiply(my[m-1][m], Mu, tempbfl0);
					LiDIA::multiply(B[m], B[m], tempbfl0);
					for (lidia_size_t i = m+1; i < k+1; i++) {
						LiDIA::multiply(tempbfl0, Mu, my[m][i]);
						LiDIA::subtract(tempbfl1, my[m-1][i], tempbfl0);
						LiDIA::multiply(tempbfl0, tempbfl1, my[m-1][m]);
						LiDIA::add(my[m-1][i], my[m][i], tempbfl0);
						my[m][i].assign(tempbfl1);
					}
				}
				B[m-1].assign(Bz);
				li.lll.swaps++;
				vector_bfl.swap(value[m-1], value[m]);
				vector_bin.swap(Tr[m-1], Tr[m]);
				for (j = 0; j <= m-2; j++) {
					tempbfl1.assign(my[j][m-1]);
					my[j][m-1].assign(my[j][m]);
					my[j][m].assign(tempbfl1);
				}
				if (m > 1)
					m--;
				l = m-1;
				continue;
			}
		}
		l--;
		if (l < 0) {
			m++;
			l = m-1;
		}
	}
	return(NULL);
}



//
// Buchmann - Kessler version for generating systems
//
void bigfloat_lattice::TrD_buchmann_kessler(dense_alg< bigfloat > & da,
					    lattice_info& li)
{
	debug_handler("bigfloat_lattice", "TrD_buchmann_kessler(da, li) [1]");
	bigfloat *help;
	bigfloat *rel;
	bigfloat *tempPbfl;
	bigfloat **oldvalue = da.s.value;
	bigfloat *olddelvalue = da.s.delvalue;
	bigfloat **tempPPbfl;
	bigfloat *deltempPPbfl;
	bigfloat vor, rechts;
	bigfloat zweipotq, alpha;
	bigfloat norm1, norm2;
	bigint bi_norm1, bi_norm2;
	sdigit n2 = da.b.rows;
	sdigit prec;
	p_vector< bigfloat > vector;
	lll_kernel_fu< bigfloat, bigfloat, vector_op < bigfloat, bigfloat >,
                basis_modules< bigfloat, bigfloat, Normal > > alg;

	//
	// Compute bigint approximation of lattice
	// to use the schnorr - euchner version of lll
	//
	prec = compute_precision(da);
	bigfloat::set_precision(prec);
	alpha_compute(da, alpha);
	zwei_pot_q_compute(da, zweipotq, n2, alpha);

	//
	// Allocate memory
	//
	da.s.value = new bigfloat*[da.b.rows];
	memory_handler(da.s.value, "bigfloat_lattice", "Tr_buchmann_kessler(da, li) :: "
		       "not enough memory !");
	da.s.delvalue = new bigfloat[(da.b.rows+da.b.columns)*da.b.rows];
	memory_handler(da.s.delvalue, "bigfloat_lattice", "Tr_buchmann_kessler(da, li) :: "
		       "not enough memory !");
	tempPPbfl = new bigfloat*[da.b.rows];
	memory_handler(tempPPbfl, "bigfloat_lattice", "Tr_buchmann_kessler(da, li) :: "
		       "not enough memory !");
	deltempPPbfl = new bigfloat[da.b.rows*da.b.columns+2*da.b.columns+da.b.rows];
	memory_handler(deltempPPbfl, "bigfloat_lattice", "Tr_buchmann_kessler(da, li) :: "
		       "not enough memory !");
	//
	// Set pointer
	//
	for (lidia_size_t i = 0; i < da.b.rows; i++) {
		da.s.value[i] = &da.s.delvalue[i*(da.b.columns+da.b.rows)];
		tempPPbfl[i] = &deltempPPbfl[i*da.b.columns];
	}
	help = &deltempPPbfl[da.b.rows*da.b.columns];
	rel = &help[da.b.columns];
	tempPbfl = &rel[da.b.rows];

	//
	// Prepare matrix for buchmann  - kessler
	//
	for (lidia_size_t i = 0; i < da.b.rows; ++i) {
		da.s.value[i][i+da.b.columns].assign_one(); // 1, wenn i = j, 0 sonst
		for (lidia_size_t j = 0; j < da.b.columns; ++j) {
			LiDIA::multiply(da.s.value[i][j], zweipotq, oldvalue[i][j]);
			tempPPbfl[i][j].assign(da.s.value[i][j]); // merke Kopie
		}
	}

	//
	// Compute needed Precision for approximation bigfloats
	//
	prec = compute_read_precision(da);
	da.d.bit_prec = static_cast<sdigit>(static_cast<double>(prec)*(std::log(10.0)/std::log(2.0)));
	da.d.bit_prec = ((da.d.bit_prec/200)+1)*DOUBLE_MANTISSA_BITS;

	//
	// Perform lll
	//
	da.b.columns += da.b.rows;
	debug_handler("bigfloat_lattice", "Tr_buchmann_kessler(da, li) [2]");
	alg.lll(da, li);
	debug_handler("bigfloat_lattice", "Tr_buchmann_kessler(da, li) [3]");
	da.b.columns -= da.b.rows;

	//
	// Check rank of the lattice
	//
	lidia_size_t l = 0;
	do {
		vector.vectsize = da.b.columns;
		vector.assign_zero(help); // Initializes help with the zero - vector
		for (lidia_size_t j = 0; j < da.b.rows; ++j) {
			rel[j].assign(da.s.value[l][j+da.b.columns]);
			vector.scalmul(tempPbfl, rel[j], tempPPbfl[j]);
			vector.add(help, help, tempPbfl);
		}
		vector.scalprod(norm2, help, help);
		vector.vectsize = da.b.rows;
		//
		// vector.l1_norm(norm1,rel);
		//
		norm1.assign_zero();
		for (lidia_size_t j = 0; j < da.b.rows; j++)
			LiDIA::add(norm1, norm1, abs(rel[j]));

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
	debug_handler("bigfloat_lattice", "Tr_buchmann_kessler(da, li) [4]");
	if (l >= da.b.rows)
		warning_handler("bigfloat_lattice", "Tr_buchmann_kessler(da, li) :: "
				"lattice of dimension 1");

	da.b.rank = da.b.rows-l+1; // da.b.rows is the dimension of the lattice
	li.lll.rank = da.b.rank;
	//
	// Store transformation lattice in math_matrix< bigint >
	//
	//  for (lidia_size_t i=0;i<da.b.rows;i++)
	//    for (lidia_size_t j=da.b.columns;j<da.b.columns+da.b.rows-da.b.rank;j++)
	//      da.s.value[i][j].assign_zero();
	for (lidia_size_t i = 0; i < da.b.rows-da.b.rank; i++)
		for (lidia_size_t j = 0; j < da.b.rows-i-1; j++)
			vector.swap(da.s.value[j], da.s.value[j+1]);


	//
	// Restore original matrix !
	//
	for (lidia_size_t i = 0; i < da.b.rows; i++)
		for (lidia_size_t j = 0; j < da.b.columns; j++)
			LiDIA::swap(da.s.value[i][j], oldvalue[i][j]);

	//
	// Free allocated storage
	//
	delete[] oldvalue;
	delete[] olddelvalue;
	delete[] tempPPbfl;
	delete[] deltempPPbfl;
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
