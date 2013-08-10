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
//	Author	: Michael Jacobson (MJ)
//	Changes	: See CVS log
//
//==============================================================================================


#ifdef HAVE_CONFIG_H
# include	"config.h"
#endif
#include	"LiDIA/quadratic_order.h"
#include	"LiDIA/qi_class.h"
#include	"LiDIA/qi_class_real.h"



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



///////////////////////
// PRIVATE FUNCTIONS //
///////////////////////


//
// quadratic_order::rBJT
//
// Task:
//      computes R using the BJT algorithm
//

void
quadratic_order::rBJT()
{
	debug_handler("quadratic_order", "rBJT");

	bigfloat y, usqr, u, temp;
	long upper;
	bigint oa, temp2;
	qi_class_real A, C, D, *F, UU;

	UU.assign_one();

	if (prin_list.no_of_elements() > 0) {
		// if a list of principal ideals already exists, use it
		A.assign(prin_list.last_entry());
		ceil(temp, A.get_distance());
		temp.bigintify(temp2);
		if (temp2.is_odd())
			inc(temp2);
		u.assign(temp2);
	}
	else {
		// get hash table size
		temp = ceil(power(bigfloat(Delta), bigfloat(0.25)));
		temp.multiply_by_2();
		multiply(temp, temp, bigfloat(1.33));
		temp.longify(upper);
		prin_list.initialize(static_cast<lidia_size_t>(upper));
		prin_list.set_key_function(&qi_class_real_key);

		// compute initial step-width = ceil((log(Delta)+3)/4)
		add(temp, log(bigfloat(Delta)), 3);
		shift_right(temp, temp, 2);
		ceil(temp, temp);
		temp.bigintify(temp2);
		if (temp2.is_odd())
			inc(temp2);
		u.assign(temp2);

		A.assign_one();
	}

	D.assign_one();
	y.assign(u);

	// oa = (D - b^2) / 4a
	oa.assign(abs(A.get_c()));

	while (R.is_zero()) {
		// compute more baby steps
		R = A.make_list(oa, u, prin_list);
		if (R.compare(min_log) < 0)
			R.assign_zero();

		// compute giant steps to u^2
		if (R.is_zero()) {
			inverse(C, A);
			C.adjust_abs(u);
			square(usqr, u);
		}
		while ((y.compare(usqr) < 0) && (R.is_zero())) {
			multiply_real(D, D, C);
			D.adjust_abs(y);

			F = prin_list.search(D);
			if (F) {
				// found D in list:  R = y-r
				subtract(R, D.get_distance(), F->get_distance());
				R.absolute_value();
				if (R.compare(min_log) < 0)
					R.assign_zero();
			}

			add(y, y, u);
		}

		u.multiply_by_2();
	}
}



//
// quadratic_order::rshanks
//
// Task:
//      computes the regulator using Shanks' algorithm
//

void
quadratic_order::rshanks()
{
	debug_handler("quadratic_order", "rshanks");

	bigfloat y, usqr, u, hR, temp, E, Aval, Fval;
	long upper, OQ;
	bigint oa, temp2;
	qi_class_real A, C, CINV, BB, D, *F, UU;

	UU.assign_one();

	// compute approximation of h
	OQ = get_optimal_Q_cnum();
	E = estimate_L1(OQ);
	E *= sqrt(bigfloat(Delta)) / (bigfloat(2.0));

	Aval = estimate_L1_error(OQ);
	Fval = exp(Aval) - bigfloat(1.0);
	temp = bigfloat(1.0) - exp(-Aval);
	if (temp > Fval)
		Fval = temp;
	multiply(Fval, Fval, E);
	sqrt(Fval, Fval);
	ceil(temp2, Fval);
	if (temp2.is_odd())
		inc(temp2);
	u.assign(temp2);

	power(temp, bigfloat(Delta), bigfloat(0.2));
	ceil(temp, temp);

	if (prin_list.no_of_elements() > 0) {
		// if a list of principal ideals already exists, use it
		A.assign(prin_list.last_entry());
		if (A.get_distance() > u) {
			ceil(temp, A.get_distance());
			temp.bigintify(temp2);
			if (temp2.is_odd())
				inc(temp2);
			u.assign(temp2);
		}
	}
	else {
		// get hash table size
		temp.multiply_by_2();
		multiply(temp, temp, bigfloat(1.33));
		temp.longify(upper);
		prin_list.initialize(static_cast<lidia_size_t>(upper));
		prin_list.set_key_function(&qi_class_real_key);

		// compute initial step-width = ceil((log(Delta)+3)/4)
		add(temp, log(bigfloat(Delta)), 3);
		shift_right(temp, temp, 2);
		ceil(temp, temp);
		temp.bigintify(temp2);
		if (temp2.is_odd())
			inc(temp2);
		if (bigfloat(temp2) > u)
			u.assign(temp2);

		A.assign_one();
		prin_list.hash(A);
	}

	y.assign_zero();
	D.assign(nearest(UU, E));
	BB.assign(D);

	// oa = (D - b^2) / 4a
	oa.assign(abs(A.get_c()));

	while (R.is_zero()) {
		// compute more baby steps
		R = A.make_list(oa, u, prin_list);
		if (R.compare(min_log) < 0)
			R.assign_zero();

		if (D.is_one()) {
			hR.assign(D.get_distance());
			if (hR.compare(min_log) >= 0)
				R = find_hstar(hR, E, u, Fval);
			else
				hR.assign_zero();
		}

		// compute giant steps to u^2
		if (R.is_zero()) {
			C.assign(A);
			C.adjust_pos(u);

			inverse(CINV, C);
			CINV.adjust_abs(u);
			square(usqr, u);
		}

		while ((y.compare(usqr) < 0) && (R.is_zero())) {
			subtract(temp, E, y);
			D.adjust_pos(temp);
			F = prin_list.search(D);
			if (F) {
				// found BINV in list
				subtract(hR, D.get_distance(), F->get_distance());
				if (hR.compare(min_log) >= 0)
					R = find_hstar(hR, E, u, Fval);
				else
					hR.assign_zero();
			}

			if (R.is_zero()) {
				add(temp, E, y);
				BB.adjust_pos(temp);

				F = prin_list.search(BB);
				if (F) {
					// found B in list
					subtract(hR, BB.get_distance(), F->get_distance());
					if (hR.compare(min_log) >= 0)
						R = find_hstar(hR, E, u, Fval);
					else
						hR.assign_zero();
				}
			}

			if (R.is_zero()) {
				multiply_real(D, D, CINV);
				multiply_real(BB, BB, C);
				add(y, y, u);
			}
		}

		u.multiply_by_2();
	}


//  	subtract(temp, E, hR);
//  	temp.absolute_value();
//  	divide(temp, Fval*Fval, temp);
//  	std::cout << Delta << "   " << temp << "   " << E << "   " << hR << "   " << Fval*Fval << "   " << abs(E-hR) << std::endl;

}



//
// quadratic_order::find_hstar
//
// Task:
//      computes the regulator given a multiple of it.
//

bigfloat
quadratic_order::find_hstar(const bigfloat & hR, const bigfloat & E,
                            const bigfloat & u, const bigfloat & Lval)
{
	debug_handler("quadratic_order", "find_hstar");

	bigfloat Eval, Rval, y, Bval, temp, dist, Lv2;
	bigint QQ, P;
	qi_class_real A, B, C, CINV, *F, UU, TT;
	base_vector< qi_class_real > Ilist;
	lidia_size_t i;
	long upper;

	Ilist.set_mode(EXPAND);
	UU.assign_one();

	Rval.assign_zero();
	if (prin_list.no_of_elements() > 0)
		C = prin_list.last_entry();
	else
		C.assign_one();
	divide(Bval, E, sqrt(Lval));

	// check whether R < E / sqrt(dist(C))
	y.assign_zero();
	B.assign(C);
	while ((y.compare(Bval) <= 0) && (Rval.is_zero())) {
		add(y, y, u);
		multiply_real(B, B, C);
		B.adjust_pos(y);

		F = prin_list.search(B);
		if (F) {
			subtract(Rval, B.get_distance(), F->get_distance());
			Rval.absolute_value();
			if (Rval.compare(min_log) < 0)
				Rval.assign_zero();
		}
	}

	// if R > E / sqrt(dist(C)), compute h*
	if (Rval.is_zero()) {
		Rval.assign(hR);

		// Eval = (E + Lval^2) / Bval
		square(temp, Lval);
		add(temp, temp, E);
		divide(Eval, temp, Bval);
		if (Eval.compare(hR) > 0)
			Eval.assign(hR);

		// compute table of powers of C
		shift_right(dist, hR, 1);
		temp.assign(C.get_distance());
		Lv2.assign(temp);
		i = 0;
		A.assign(C);
		Ilist[0] = A;
		while (A.get_distance() < dist) {
			++i;
			temp.multiply_by_2();
			square(A, A);
			A.adjust_pos(temp);
			Ilist[i] = A;
		}

		Eval.longify(upper);
		if (upper > 2)
			if (!(upper & 1))
				++upper;
		if ((static_cast<long>(PL.get_upper_bound()) < upper) || (PL.get_lower_bound() != 2))
			PL.resize(2, upper);
		P = PL.get_last_prime(upper);

		A.assign_one();

		while (P.is_gt_zero()) {
			divide(dist, Rval, bigfloat(P));
			subtract(temp, dist, A.get_distance());
			divide(temp, temp, Lv2);
			floor(temp, temp);
			temp.bigintify(QQ);
			inc(QQ);
			i = 0;
			while (QQ.is_gt_zero()) {
				if (QQ.is_odd())
					multiply_real(A, A, Ilist[i]);
				QQ.divide_by_2();
				++i;
			}

			while (A.get_distance() <= dist)
				A.rho();
			add(temp, dist, C.get_distance());
			while (A.get_distance() >= temp)
				A.inverse_rho();

			F = prin_list.search(A);
			if (F) {
				subtract(temp, A.get_distance(), F->get_distance());
				subtract(temp, temp, dist);
				temp.absolute_value();
				TT = nearest(UU, temp);
				if (TT.get_distance() > temp)
					TT.inverse_rho();
				if (TT.is_one()) {
					divide(Rval, Rval, bigfloat(P));
					A.assign_one();
				}
				else {
					P = PL.get_prev_prime();
				}
			}
			else {
				P = PL.get_prev_prime();
			}
		}
	}

	return Rval;
}



//
// quadratic_order::cns_imag
//
// Task:
//      computes the class number using the analytic class number formula and
//      a divisor of h.
//

void
quadratic_order::cns_imag()
{
	debug_handler("quadratic_order", "cns_imag");

	bigfloat nFI, temp, temp1, temp2, F, A, k, B1, B2, val, bfone;
	bigint h1;

	// compute approximation of h
	Q = get_optimal_Q_cnum();
	nFI.assign(estimate_L1(Q));
	nFI *= sqrt(bigfloat(-Delta)) / Pi();
	if (Delta == -3)
		nFI *= bigfloat(3.0);
	if (Delta == -4)
		nFI *= bigfloat(2.0);
	nFI.bigintify(h);

	A = estimate_L1_error(Q);
	h1.assign_one();

	// compute divisor of h, if necessary
	if (h > 3) {
		temp1 = exp(A) - bigfloat(1.0);
		temp2 = bigfloat(1.0) - exp(-A);
		if (temp1 > temp2)
			F = temp1;
		else
			F = temp2;
		cgs_imag(nFI, F, h1);
		nFI /= h1;
		nFI.bigintify(h);
	}

	// compute sufficiently accurate estimate of L(1,X) to prove correctness
	bfone = bigfloat(-1.0);
	k = nFI - bigfloat(h);
	B1 = nFI*exp(A);
	B2 = abs(A) + B1 - nFI;
	val = bigfloat(h) - B2 - floor(B2 + bigfloat(h));
	while (val <= bfone) {
		Q += 10000;
		nFI.assign(estimate_L1(Q));
		nFI *= sqrt(bigfloat(-Delta)) / (Pi()*bigfloat(h1));
		nFI.bigintify(h);
		A = estimate_L1_error(Q);
		k = nFI - bigfloat(h);
		B1 = nFI*exp(A);
		B2 = abs(A) + B1 - nFI;
		val = bigfloat(h) - B2 - floor(B2 + bigfloat(h));
	}

	add(B2, B2, h);
	temp.assign(floor(B2));
	temp.bigintify(h);
	h *= h1;
}



//
// quadratic_order::cns_real
//
// Task:
//      computes the class number using the analytic class number formula and
//      a divisor of h.
//

void
quadratic_order::cns_real()
{
	debug_handler("quadratic_order", "cns_real");

	bigfloat nFI, temp, temp1, temp2, F, A, val;
	bigint h1;

	// compute regulator, if needed
	if (R.is_zero())
		regulator();

	// exit if class number is computed via regulator routine
	if (!h.is_zero())
		return;

	// compute estimate of h
	Q = get_optimal_Q_cnum();
	nFI.assign(estimate_L1(Q));
	nFI *= sqrt(bigfloat(Delta)) / (R*bigfloat(2.0));
	nFI.bigintify(h);

	A = estimate_L1_error(Q);
	val.assign(h);
	inc(val);
	temp = bigfloat(h) + abs(nFI - bigfloat(h));
	val = log(val/temp);
	h1.assign_one();

	// compute divisor of h, if necessary
	if ((h > 3) && (A >= val)) {
		temp1 = exp(A) - bigfloat(1.0);
		temp2 = bigfloat(1.0) - exp(-A);
		if (temp1 > temp2)
			F = temp1;
		else
			F = temp2;
		cgs_real(nFI, h1);
		nFI /= h1;
		nFI.bigintify(h);

		val.assign(h);
		inc(val);
		temp = bigfloat(h) + abs(nFI - bigfloat(h));
		val = log(val/temp);
	}

	// compute sufficiently accurate estimate of L(1,X) to prove correctness
	while (A >= val) {
		Q += 5000;
		nFI.assign(estimate_L1(Q));
		nFI *= sqrt(bigfloat(Delta)) / (R*bigfloat(h1 << 1));
		nFI.bigintify(h);

		A = estimate_L1_error(Q);
		val.assign(h);
		inc(val);
		temp = bigfloat(h) + abs(nFI - bigfloat(h));
		val = log(val/temp);
	}

	h *= h1;
}



//
// quadratic_order::set_transformations
//
// Task:
//      computes the transformation matrix needed to compute a set of
//      generators.
//

void
quadratic_order::set_transformations()
{
	debug_handler("quadratic_order", "set_transformations");

	matrix< bigint > Imat(1, 1), temp(1, 1);
	lidia_size_t dif, numFB, numCL;

	numFB = contributors.size();
	numCL = CL.size();
	dif = numFB - numCL;

	Imat.resize(numFB, numFB);
	Imat.diag(1, 0);
	UINV.reginvimage(U, Imat);
	UINV.set_no_of_rows(numFB);

	// compute inverse
	if (dif) {
		// remove unneccessary rows of U and cols of UINV
		temp.resize(numCL, numFB);
		Imat.resize(dif, numFB);
		U.split_v(Imat, temp);
		U.assign(temp);

		temp.resize(numFB, numCL);
		Imat.resize(numFB, dif);
		UINV.split_h(Imat, temp);
		UINV.assign(temp);
	}
}



//
// quadratic_order::cgBJT_imag
//
// Task:
//      computes the class group using the BJT algorithm.
//

void
quadratic_order::cgBJT_imag()
{
	debug_handler("quadratic_order", "cgBJT_imag");

	lidia_size_t pidx;
	bigint y, usqr, det, Bjj;
	long s, u, r, q, upper, Bj;
	lidia_size_t i, j, k, crank, numRpr, numQ, curr_index, idx;
	qi_class G, A, B, C, D, E, HI, Gq, GBj;
	matrix< bigint > Bmat, junk(1, 1);
	base_vector< long > Rvec, Qvec;
	ideal_node *Inode;
	indexed_hash_table< ideal_node > QT, RT;
	bigfloat temp, hstar, nFI;

	Rvec.set_mode(EXPAND);
	Qvec.set_mode(EXPAND);
	Rvec.reset();
	Qvec.reset();

	// initialize sets R and Q
	temp = ceil(power(bigfloat(-Delta), bigfloat(0.25)));
	temp.longify(upper);

	QT.initialize(static_cast<lidia_size_t>(upper) << 1);
	QT.set_key_function(&ideal_node_key);
	B.assign_one();
	QT.hash(ideal_node(B, 0));

	RT.initialize(static_cast<lidia_size_t>(upper) << 1);
	RT.set_key_function(&ideal_node_key);
	RT.hash(ideal_node(B, 0));

	crank = 0;
	det.assign_one();
	h.assign_one();
	pidx = 0;

	// compute estimate of h
	Q = get_optimal_Q();
	nFI.assign(estimate_L1(Q));
	nFI *= sqrt(bigfloat(-Delta)) / Pi();
	if (Delta == -3)
		nFI *= bigfloat(3.0);
	if (Delta == -4)
		nFI *= bigfloat(2.0);
	hstar = nFI / sqrt(bigfloat(2.0));

	// while h < hstar, we only have a subgroup
	while (bigfloat(h) < hstar) {
		// get prime ideal
		while (!generate_prime_ideal(G, bigint(SPL[pidx])))
			++pidx;
		fact_base[crank] = G;
		contributors[crank] = crank;
		s = 1;

		// compute initial step-width
		temp.assign(floor(sqrt(nFI)));
		temp.longify(u);
		if (u == 0) u = 2;
		if ((u & 1) == 1)
			++u;
		y.assign(u);

		A.assign_one();
		power_imag(B, G, u);
		C.assign(B);
		HI.assign(inverse(G));

		numRpr = RT.no_of_elements();
		curr_index = numRpr;
		numQ = QT.no_of_elements();

		// if G is one, we have no new information
		if (G.is_one())
			Bjj.assign_one();
		else
			Bjj.assign_zero();

		// check whether current ideal is in previously generated subgroup
		if ((Bjj.is_zero()) && (crank > 0)) {
			for (i = 0; i < numQ; ++i) {
				D.assign(QT[i].get_A());
				multiply_imag(E, D, G);
				Inode = RT.search(ideal_node(E, 0));
				if ((E.is_one()) || (Inode)) {
					Bjj.assign_one();
					break;
				}
			}
		}

		while (Bjj.is_zero()) {
			// compute more baby steps
			for (r = s; r <= u; ++r) {
				multiply_imag(A, A, HI);

				if ((s == 1) && (r > 1)) {
					for (k = 0; k < numRpr; ++k) {
						D.assign(RT[k].get_A());
						multiply_imag(E, D, A);

						// check if relation already found
						Inode = QT.search(ideal_node(E, 0));
						if (Inode) {
							q = Inode->get_index();
							Bjj.assign(r);
							decode_vector(Bmat, Bjj, k, q, Rvec, Qvec, numRpr, numQ);
							break;
						}
						else {
							RT.hash(ideal_node(E, curr_index));
							++curr_index;
						}
					}
					if (!Bjj.is_zero())
						break;
				}
				else {
					for (k = 0; k < numRpr; ++k) {
						D.assign(RT[k].get_A());
						multiply_imag(E, D, A);

						RT.hash(ideal_node(E, curr_index));
						++curr_index;
					}
				}
			}

			// compute giant steps to u^2
			square(usqr, u);
			while ((y.compare(usqr) < 0) && (Bjj.is_zero())) {
				// search
				for (i = 0; i < numQ; ++i) {
					E.assign(QT[i].get_A());
					multiply_imag(D, E, B);

					Inode = RT.search(ideal_node(D, 0));
					if (Inode) {
						r = Inode->get_index();
						q = i;
						add(Bjj, y, (r/numRpr));
						if (Bjj > 1) {
							r %= numRpr;
							decode_vector(Bmat, Bjj, r, q, Rvec, Qvec, numRpr, numQ);
						}
						break;
					}
				}

				if (Bjj.is_zero()) {
					// not found, take another giant step
					add(y, y, u);
					multiply_imag(B, B, C);
				}
			}

			// double u
			if (Bjj.is_zero()) {
				s = u+1;
				u <<= 1;
				square_imag(C, C);
			}
		}

		if (!Bjj.is_one()) {
			h *= Bjj;
			divide(nFI, nFI, bigfloat(Bjj));
			++crank;
		}

		if (bigfloat(h) < hstar) {
			temp = ceil(sqrt(bigfloat(Bjj)));
			temp.longify(Bj);

			// compute new R' (remove entries with too large exponents)
			det *= Bj;
			det.longify(upper);
			idx = static_cast<lidia_size_t>(upper);
			numRpr = RT.no_of_elements();
			for (i = numRpr-1; i >= idx; --i)
				RT.remove_from(i);

			if (!Bjj.is_one()) {
				Rvec[crank-1] = Bj;
				Qvec[crank-1] = 1;

				// compute new Q
				numQ = QT.no_of_elements();
				curr_index = numQ;
				power_imag(GBj, G, Bj);
				Gq.assign(GBj);
				for (i = 1; i < Bj; ++i) {
					if (Bjj.compare(i*Bj) <= 0)
						break;
					for (k = 0; k < numQ; ++k) {
						E.assign(QT[k].get_A());
						multiply_imag(D, E, Gq);
						QT.hash(ideal_node(D, curr_index));
						++curr_index;
					}
					multiply_imag(Gq, Gq, GBj);
					++Qvec[crank-1];
				}
			}

			++pidx;
		}
	}

	// compute structure
//  RHNF.assign(Bmat);
	if (h == 1)
		CL[0] = 1;
	else {
		// delete all rows and columns with diagonal 1
		Bmat.snf_havas(U, junk);
		i = 0;
		for (j = 0; j < crank; ++j) {
			Bjj.assign(Bmat.member(j, j));
			if (!Bjj.is_one())
				CL[i++] = Bjj;
		}
	}
}



//
// quadratic_order::cgBJT_real
//
// Task:
//      computes the class group using the BJT algorithm.
//

void
quadratic_order::cgBJT_real()
{
	debug_handler("quadratic_order", "cgBJT_real");

	lidia_size_t pidx;
	bigint y, usqr, det, Bjj, ht;
	long s, u, r, q, upper, Bj;
	lidia_size_t i, j, k, crank, numRpr, numQ, curr_index, idx, numRreps;
	qi_class G, A, Aprime, B, C, D, E, HI, Gq, GBj;
	qi_class_real Areal, F, Gstep, *FS;
	matrix< bigint > Bmat, junk(1, 1);
	base_vector< long > Rvec, Qvec;
	ideal_node *Inode;
	indexed_hash_table< ideal_node > QT, RT;
	base_vector< qi_class > Rreps;
	bigfloat temp, hstar, nFI, sqReg, val, Aval, GStepWidth;

	Q = get_optimal_Q_cnum();
	nFI.assign(estimate_L1(Q));

	// compute regulator, if needed
	if (R.is_zero())
		regulator();

	// exit if class number is computed via regulator routine
	if (!h.is_zero())
		return;

	crank = 0;
	det.assign_one();
	h.assign_one();
	pidx = 0;

	// compute estimate of h
	nFI *= sqrt(bigfloat(Delta)) / (R*bigfloat(2.0));
	nFI.bigintify(ht);
	hstar = nFI / sqrt(bigfloat(2.0));

	// test if h is a small prime
	Aval = estimate_L1_error(Q);
	val.assign(ht);
	inc(val);
	temp = bigfloat(ht) + abs(nFI - bigfloat(ht));
	val = log(val/temp);
	if ((Aval < val) && (ht < bigint(8)) &&
	    (ht != bigint(4)) && (ht != bigint(6))) {
		h.assign(ht);
		CL[0] = h;
		// find generator
		while (!generate_prime_ideal(G, SPL[pidx]))
			++pidx;
		while ((h > 1) && (G.is_principal())) {
			++pidx;
			while (!generate_prime_ideal(G, SPL[pidx]))
				++pidx;
		}
		Bmat.resize(1, 1);
		Bmat.sto(0, 0, h);
		fact_base[0] = G;
		U.sto(0, 0, bigint(1));
		UINV.sto(0, 0, bigint(1));
		if (!h.is_one())
			gens[0] = G;
	}

	if (bigfloat(h) < hstar) {
		Rvec.set_mode(EXPAND);
		Qvec.set_mode(EXPAND);
		Rreps.set_mode(EXPAND);
		Rvec.reset();
		Qvec.reset();
		Rreps.reset();

		// initialize sets R and Q
		temp = ceil(power(bigfloat(Delta), bigfloat(0.25)));
		temp.longify(upper);

		QT.initialize(static_cast<lidia_size_t>(upper) << 1);
		QT.set_key_function(&ideal_node_key);
		B.assign_one();
		QT.hash(ideal_node(B, 0));

		RT.initialize(static_cast<lidia_size_t>(upper) << 1);
		RT.set_key_function(&ideal_node_key);
		RT.hash(ideal_node(B, 0));

		Rreps[0] = B;

		// compute giant step for equivalence testing
		sqrt(sqReg, R);
		F.assign_one();
		Gstep.assign(nearest(F, sqReg));
		if (Gstep.get_distance() > sqReg)
			Gstep.inverse_rho();
		F.assign(prin_list.last_entry());
		if (Gstep.is_one()) {
			while (!F.is_one()) {
				F.rho();
				prin_list.hash(F);
			}
		}
		else  {
			while ((F.get_distance() < Gstep.get_distance()) && (!F.is_one())) {
				F.rho();
				prin_list.hash(F);
			}
			if (!F.is_one()) {
				F.rho();
				prin_list.hash(F);
			}
			if (F.is_one())
				Gstep.assign(F);
		}
		floor(GStepWidth, Gstep.get_distance());
	}

	// while h < hstar, we only have a subgroup
	while (bigfloat(h) < hstar) {
		// get prime ideal
		while (!generate_prime_ideal(G, bigint(SPL[pidx])))
			++pidx;
		fact_base[crank] = G;
		contributors[crank] = crank;

		s = u = 1;
		y.assign_one();

		A.assign_one();
		B.assign(G);
		C.assign(B);
		HI.assign(inverse(G));

		numQ = QT.no_of_elements();
		curr_index = numRpr = numRreps = Rreps.size();

		Bjj.assign_zero();

		// check whether current ideal is in previously generated subgroup
		for (i = 0; i < numQ; ++i) {
			D.assign(QT[i].get_A());
			multiply_real(E, D, G);

			if (Gstep.is_one()) {
				FS = prin_list.search(qi_class_real(E));
				if (FS)
					Bjj.assign_one();
				else {
					Inode = RT.search(ideal_node(E, 0));
					if (Inode)
						Bjj.assign_one();
				}
			}
			else {
				F.assign(E, 0.0);
				temp.assign_zero();
				while ((F.get_distance() <= R) && (Bjj.is_zero())) {
					FS = prin_list.search(F);
					if (FS)
						Bjj.assign_one();
					else {
						Inode = RT.search(ideal_node(qi_class(F), 0));
						if (Inode)
							Bjj.assign_one();
						else {
							add(temp, temp, GStepWidth);
							multiply_real(F, F, Gstep);
							F.adjust_pos(temp);
						}
					}
				}
			}

			if (Bjj.is_one())
				break;
		}

		while (Bjj.is_zero()) {
			// compute more baby steps
			for (r = s; r <= u; ++r) {
				multiply_real(A, A, HI);

				for (k = 0; k < numRreps; ++k) {
					D.assign(Rreps[k]);
					multiply_real(E, D, A);
					Rreps[curr_index] = E;
					RT.hash(ideal_node(E, curr_index));

					if (Gstep.is_one()) {
						// store each cycle of ideals
						E.make_cycle(curr_index, RT);
					}
					else {
						// store ideals with distance < sqReg
						Areal.assign(E);
						Areal.make_list(sqReg, curr_index, RT);
					}

					++curr_index;
				}
			}

			// compute giant steps to u^2
			square(usqr, u);
			while ((y.compare(usqr) <= 0) && (Bjj.is_zero())) {
				// search
				for (i = 0; i < numQ; ++i) {
					D.assign(QT[i].get_A());
					multiply_real(E, D, B);

					if (Gstep.is_one()) {
						FS = prin_list.search(qi_class_real(E));
						if (FS) {
							r = 0;
							q = i;
							add(Bjj, y, (r/numRpr));
							r %= numRpr;
							decode_vector(Bmat, Bjj, r, q, Rvec, Qvec, numRpr, numQ);
						}
						else {
							Inode = RT.search(ideal_node(E, 0));
							if (Inode) {
								r = Inode->get_index();
								q = i;
								add(Bjj, y, (r/numRpr));
								r %= numRpr;
								decode_vector(Bmat, Bjj, r, q, Rvec, Qvec, numRpr, numQ);
							}
						}
					}
					else {
						F.assign(E, 0.0);
						temp.assign_zero();
						while ((F.get_distance() <= R) && (Bjj.is_zero())) {
							FS = prin_list.search(F);
							if (FS) {
								r = 0;
								q = i;
								add(Bjj, y, (r/numRpr));
								r %= numRpr;
								decode_vector(Bmat, Bjj, r, q, Rvec, Qvec, numRpr, numQ);
							}
							else {
								Inode = RT.search(ideal_node(qi_class(F), 0));
								if (Inode) {
									r = Inode->get_index();
									q = i;
									add(Bjj, y, (r/numRpr));
									r %= numRpr;
									decode_vector(Bmat, Bjj, r, q, Rvec, Qvec, numRpr, numQ);
								}
								else {
									add(temp, temp, GStepWidth);
									multiply_real(F, F, Gstep);
									F.adjust_pos(temp);
								}
							}
						}
					}

					if (Bjj.is_gt_zero())
						break;
				}

				if (Bjj.is_zero()) {
					// not found, take another giant step
					add(y, y, u);
					multiply_real(B, B, C);
				}
			}

			// double u
			if (Bjj.is_zero()) {
				s = u+1;
				u <<= 1;
				square_real(C, C);
			}
		}

		if (!Bjj.is_one()) {
			h *= Bjj;
			++crank;
		}

		if (bigfloat(h) < hstar) {
			temp = ceil(sqrt(bigfloat(Bjj)));
			temp.longify(Bj);

			// compute new R' (remove entries with too large exponents)
			det *= Bj;
			det.longify(upper);
			idx = static_cast<lidia_size_t>(upper);
			numRpr = RT.no_of_elements();
			i = numRpr-1;
			while (RT[i].get_index() >= idx) {
				RT.remove_from(i);
				--i;
			}
			Rreps.set_size(static_cast<lidia_size_t>(idx));

			if (!Bjj.is_one()) {
				Rvec[crank-1] = Bj;
				Qvec[crank-1] = 1;

				// compute new Q
				numQ = QT.no_of_elements();
				curr_index = numQ;
				power_real(GBj, G, Bj);
				Gq.assign(GBj);
				for (i = 1; i < Bj; ++i) {
					if (Bjj.compare(i*Bj) <= 0)
						break;
					for (k = 0; k < numQ; ++k) {
						E.assign(QT[k].get_A());
						multiply_real(D, E, Gq);
						QT.hash(ideal_node(D, curr_index));
						++curr_index;
					}
					multiply_real(Gq, Gq, GBj);
					++Qvec[crank-1];
				}
			}

			++pidx;
		}
	}

	// compute structure
//  RHNF.assign(Bmat);
	if (h == 1)
		CL[0] = 1;
	else {
		Bmat.snf_havas(U, junk);
		i = 0;
		for (j = 0; j < crank; ++j) {
			Bjj.assign(Bmat.member(j, j));
			if (!Bjj.is_one())
				CL[i++] = Bjj;
		}
	}

}



//
// quadratic_order::cgs_imag
//
// Task:
//      computes the class group using Shanks' algorithm.
//

void
quadratic_order::cgs_imag(const bigfloat & tFI, const bigfloat & tF,
                          bigint & h1)
{
	debug_handler("quadratic_order", "cgs_imag");

	lidia_size_t pidx;
	bigint y, usqr, det, Bjj, ht;
	long s, u, r, q, upper, Bj;
	lidia_size_t i, j, k, crank, numRpr, numR, numQ, curr_index, idx;
	matrix< bigint > Bmat, junk(1, 1);
	base_vector< long > Rvec, Qvec;
	qi_class G, A, B, Bcomp, C, D, E, HI, Gq, GBj, CINV;
	ideal_node *Inode;
	indexed_hash_table< ideal_node > QT, RT;
	bigfloat temp, temp1, temp2, hstar, nFI, F;
	rational_factorization htfact;

	Rvec.set_mode(EXPAND);
	Qvec.set_mode(EXPAND);
	Rvec.reset();
	Qvec.reset();

	// initialize sets R and Q
	temp = ceil(power(bigfloat(-Delta), bigfloat(0.25)));
	temp.longify(upper);

	QT.initialize(static_cast<lidia_size_t>(upper) << 1);
	QT.set_key_function(&ideal_node_key);
	B.assign_one();
	QT.hash(ideal_node(B, 0));

	RT.initialize(static_cast<lidia_size_t>(upper) << 1);
	RT.set_key_function(&ideal_node_key);
	RT.hash(ideal_node(B, 0));

	crank = 0;
	det.assign_one();
	h1.assign_one();
	pidx = 0;

	if (tFI.is_zero()) {
		// compute entire class group
		Q = get_optimal_Q_cnum();
		nFI.assign(estimate_L1(Q));
		nFI *= sqrt(bigfloat(-Delta)) / Pi();
		if (Delta == -3)
			nFI *= bigfloat(3.0);
		if (Delta == -4)
			nFI *= bigfloat(2.0);
		hstar = nFI / sqrt(bigfloat(2.0));

		temp = estimate_L1_error(Q);
		temp1 = exp(temp) - bigfloat(1.0);
		temp2 = bigfloat(1.0) - exp(-temp);
		if (temp1 > temp2)
			F = temp1;
		else
			F = temp2;
	}
	else {
		// compute divisor
		nFI.assign(tFI);
		hstar = power(bigfloat(-Delta), bigfloat(0.4));
		temp = nFI / 3;
		if (temp < hstar)
			hstar = temp;
		F.assign(tF);
	}

	while (bigfloat(h1) < hstar) {
		// get prime ideal
		while (!generate_prime_ideal(G, bigint(SPL[pidx])))
			++pidx;
		fact_base[crank] = G;
		contributors[crank] = crank;

		numRpr = RT.no_of_elements();
		curr_index = numRpr;
		numQ = QT.no_of_elements();

		s = 1;

		if (G.is_one())
			Bjj.assign_one();
		else
			Bjj.assign_zero();

		if (Bjj.is_zero()) {
			HI.assign(inverse(G));
			A.assign_one();

			nFI.bigintify(ht);
			temp = sqrt(nFI*F + abs(nFI - bigfloat(ht)));
			temp.longify(u);

			if (u == 0)
				u = 2;
			if ((u & 1) == 1)
				++u;
			y.assign_zero();

			power_imag(B, G, ht);

			if (crank == 0) {
				// check whether G^ht = (1)
				if (B.is_one()) {
					htfact.assign(ht);
					htfact.factor();
					Bjj = G.order_mult(ht, htfact);
					if (tFI.is_zero()) {
						Bmat.resize(1, 1);
						Bmat.sto(0, 0, Bjj);
					}
				}
			}
			else {
				// check whether G^ht is in currently generated subgroup
				for (i = 0; i < numQ; ++i) {
					D.assign(QT[i].get_A());
					multiply_imag(E, D, B);

					if (E.is_one())
						r = 0;
					else {
						Inode = RT.search(ideal_node(E, 0));
						if (Inode)
							r = Inode->get_index();
						else
							r = -1;
					}
					if (r >= 0) {
						// determine smallest divisor of ht
						q = i;
						Bjj = exact_power(ht, G, RT, QT, numRpr, r, q);
						if ((!Bjj.is_one()) && (tFI.is_zero()))
							decode_vector(Bmat, Bjj, r, q, Rvec, Qvec, numRpr, numQ);
						break;
					}
				}
			}

			if (Bjj.is_zero()) {
				power_imag(C, G, u);
				CINV = inverse(C);
				Bcomp.assign(B);
			}
			else
				u = 0;
		}

		while (Bjj.is_zero()) {
			// compute more baby steps
			for (r = s; r <= u; ++r) {
				multiply_imag(A, A, HI);

				if ((s == 1) && (r > 1)) {
					for (k = 0; k < numRpr; ++k) {
						D.assign(RT[k].get_A());
						multiply_imag(E, D, A);

						// check if relation already found
						Inode = QT.search(ideal_node(E, 0));
						if (Inode) {
							q = Inode->get_index();
							Bjj.assign(r);
							if (tFI.is_zero())
								decode_vector(Bmat, Bjj, k, q, Rvec, Qvec, numRpr, numQ);
							break;
						}
						else {
							RT.hash(ideal_node(E, curr_index));
							++curr_index;
						}
					}
					if (!Bjj.is_zero())
						break;
				}
				else {
					for (k = 0; k < numRpr; ++k) {
						D.assign(RT[k].get_A());
						multiply_imag(E, D, A);

						RT.hash(ideal_node(E, curr_index));
						++curr_index;
					}
				}
			}

			// compute giant steps to u^2
			square(usqr, u);
			while ((y.compare(usqr) < 0) && (Bjj.is_zero())) {
				// search
				for (i = 0; i < numQ; ++i) {
					E.assign(QT[i].get_A());
					multiply_imag(D, E, Bcomp);

					Inode = RT.search(ideal_node(D, 0));
					if (Inode) {
						r = Inode->get_index();
						subtract(Bjj, ht, y);
						add(Bjj, Bjj, (r/numRpr));
						if (!Bjj.is_zero()) {
							// determine smallest divisor of Bjj
							ht.assign(Bjj);
							if (crank == 0) {
								htfact.assign(Bjj);
								htfact.factor();
								Bjj = G.order_mult(Bjj, htfact);
								if (tFI.is_zero()) {
									Bmat.resize(1, 1);
									Bmat.sto(0, 0, Bjj);
								}
							}
							else {
								q = i;
								r %= numRpr;
								Bjj = exact_power(Bjj, G, RT, QT, numRpr, r, q);
								if ((!Bjj.is_one()) && (tFI.is_zero()))
									decode_vector(Bmat, Bjj, r, q, Rvec, Qvec, numRpr, numQ);
							}
							break;
						}
					}

					multiply_imag(D, E, B);
					Inode = RT.search(ideal_node(D, 0));
					if (Inode) {
						r = Inode->get_index();
						add(Bjj, ht, y);
						add(Bjj, Bjj, (r/numRpr));
						// determine smallest divisor of Bjj
						ht.assign(Bjj);
						if (crank == 0) {
							htfact.assign(Bjj);
							htfact.factor();
							Bjj = G.order_mult(Bjj, htfact);
							if (tFI.is_zero()) {
								Bmat.resize(1, 1);
								Bmat.sto(0, 0, Bjj);
							}
						}
						else {
							q = i;
							r %= numRpr;
							Bjj = exact_power(Bjj, G, RT, QT, numRpr, r, q);
							if ((!Bjj.is_one()) && (tFI.is_zero()))
								decode_vector(Bmat, Bjj, r, q, Rvec, Qvec, numRpr, numQ);
						}
						break;
					}
				}

				if (Bjj.is_zero()) {
					// not found, take another giant step
					add(y, y, u);
					multiply_imag(B, B, C);
					multiply_imag(Bcomp, Bcomp, CINV);
				}
			}

			// double u
			if (Bjj.is_zero()) {
				s = u+1;
				u <<= 1;
				square_imag(C, C);
				CINV = inverse(C);
			}
		}

		if (!Bjj.is_one()) {
			h1 *= Bjj;
			divide(nFI, nFI, bigfloat(Bjj));
			++crank;
		}

		if (bigfloat(h1) < hstar) {
			temp = ceil(sqrt(bigfloat(Bjj)));
			temp.longify(Bj);

			// compute new R' (remove entries with too large exponents)
			det *= Bj;
			det.longify(upper);
			idx = static_cast<lidia_size_t>(upper);
			numR = RT.no_of_elements();
			if (idx < numR) {
				// removing entries
				for (i = numR-1; i >= idx; --i)
					RT.remove_from(i);
			}
			else {
				// adding entries
				idx = static_cast<lidia_size_t>(u)+1;
				for (r = idx; r < Bj; ++r) {
					multiply_imag(A, A, HI);
					for (k = 0; k < numRpr; ++k) {
						D.assign(RT[k].get_A());
						multiply_imag(E, D, A);
						RT.hash(ideal_node(E, curr_index));
						++curr_index;
					}
				}
			}

			if (!Bjj.is_one()) {
				Rvec[crank-1] = Bj;
				Qvec[crank-1] = 1;

				// compute new Q
				numQ = QT.no_of_elements();
				curr_index = numQ;
				power_imag(GBj, G, Bj);
				Gq.assign(GBj);
				for (i = 1; i < Bj; ++i) {
					if (Bjj.compare(i*Bj) <= 0)
						break;
					for (k = 0; k < numQ; ++k) {
						E.assign(QT[k].get_A());
						multiply_imag(D, E, Gq);
						QT.hash(ideal_node(D, curr_index));
						++curr_index;
					}
					multiply_imag(Gq, Gq, GBj);
					++Qvec[crank-1];
				}
			}

			++pidx;
		}
	}

//  RHNF.assign(Bmat);
	if (tFI.is_zero()) {
		// compute structure
		h.assign(h1);
		if (h == 1)
			CL[0] = 1;
		else {
			Bmat.snf_havas(U, junk);
			i = 0;
			for (j = 0; j < crank; ++j) {
				Bjj.assign(Bmat.member(j, j));
				if (!Bjj.is_one())
					CL[i++] = Bjj;
			}
		}
	}
}



//
// quadratic_order::cgs_real
//
// Task:
//      computes the class group using Shanks' algorithm.
//

void
quadratic_order::cgs_real(const bigfloat & tFI, bigint & h1)
{
	debug_handler("quadratic_order", "cgs_real");

	lidia_size_t pidx;
	bigint y, usqr, det, Bjj, ht;
	long s, u, r, q, upper, Bj;
	lidia_size_t i, j, k, crank, numRpr, numQ, curr_index, idx, numRreps;
	matrix< bigint > Bmat, junk(1, 1);
	base_vector< long > Rvec, Qvec;
	qi_class G, A, B, Bcomp, C, D, E, HI, Gq, GBj, CINV, Aprime;
	qi_class_real Areal, FIDL, Gstep, *FS;
	ideal_node *Inode;
	indexed_hash_table< ideal_node > QT, RT;
	rational_factorization htfact;
	base_vector< qi_class > Rreps;
	bigfloat temp, temp1, temp2, hstar, nFI, sqReg, Aval, val, GStepWidth;

	crank = 0;
	det.assign_one();
	h1.assign_one();
	pidx = 0;

	if (tFI.is_zero()) {
		// compute entire class group
		Q = get_optimal_Q_cnum();
		nFI.assign(estimate_L1(Q));

		// compute regulator, if needed
		if (R.is_zero())
			regulator();

		// exit if class number is computed via regulator routine
		if (!h.is_zero())
			return;

		nFI *= sqrt(bigfloat(Delta)) / (R*bigfloat(2.0));
		hstar = nFI / sqrt(bigfloat(2.0));
		nFI.bigintify(ht);
	}
	else {
		// compute divisor
		nFI.assign(tFI);
		nFI.bigintify(ht);
		hstar = power(bigfloat(Delta), bigfloat(0.4));
		temp = nFI / 3;
		if (temp < hstar)
			hstar = temp;
	}

	if (tFI.is_zero()) {
		// check whether h is a small prime
		Aval = estimate_L1_error(Q);
		val.assign(ht);
		inc(val);
		temp = ht + abs(nFI - bigfloat(ht));
		val = log(val/temp);
		if ((Aval < val) && (ht < bigint(8)) &&
		    (ht != bigint(4)) && (ht != bigint(6))) {
			h.assign(ht);
			CL[0] = h;
			// find generator
			while (!generate_prime_ideal(G, SPL[pidx]))
				++pidx;
			while ((h > 1) && (G.is_principal())) {
				++pidx;
				while (!generate_prime_ideal(G, SPL[pidx]))
					++pidx;
			}
			Bmat.resize(1, 1);
			Bmat.sto(0, 0, h);
			fact_base[0] = G;
			U.sto(0, 0, bigint(1));
			UINV.sto(0, 0, bigint(1));
			if (!h.is_one())
				gens[0] = G;
		}
	}

	if (bigfloat(h1) < hstar) {
		Rvec.set_mode(EXPAND);
		Qvec.set_mode(EXPAND);
		Rreps.set_mode(EXPAND);
		Rvec.reset();
		Qvec.reset();
		Rreps.reset();

		// initialize sets R and Q
		temp = ceil(power(bigfloat(Delta), bigfloat(0.2)));
		temp.longify(upper);

		QT.initialize(static_cast<lidia_size_t>(upper) << 1);
		QT.set_key_function(&ideal_node_key);
		B.assign_one();
		QT.hash(ideal_node(B, 0));

		RT.initialize(static_cast<lidia_size_t>(upper) << 1);
		RT.set_key_function(&ideal_node_key);
		RT.hash(ideal_node(B, 0));

		Rreps[0] = B;

		// compute giant step for equivalence testing
		sqrt(sqReg, R);
		FIDL.assign_one();
		Gstep.assign(nearest(FIDL, sqReg));
		if (Gstep.get_distance() > sqReg)
			Gstep.inverse_rho();
		FIDL.assign(prin_list.last_entry());
		if (Gstep.is_one()) {
			while (!FIDL.is_one()) {
				FIDL.rho();
				prin_list.hash(FIDL);
			}
		}
		else {
			while ((FIDL.get_distance() < Gstep.get_distance()) &&
			       (!FIDL.is_one())) {
				FIDL.rho();
				prin_list.hash(FIDL);
			}
			if (!FIDL.is_one()) {
				FIDL.rho();
				prin_list.hash(FIDL);
			}
			if (FIDL.is_one())
				Gstep.assign(FIDL);
		}
		floor(GStepWidth, Gstep.get_distance());
	}

	while (bigfloat(h1) < hstar) {
		// get prime ideal
		while (!generate_prime_ideal(G, bigint(SPL[pidx])))
			++pidx;
		fact_base[crank] = G;
		contributors[crank] = crank;

		s = u = 1;
		y.assign_zero();

		A.assign_one();
		HI.assign(inverse(G));

		numQ = QT.no_of_elements();
		curr_index = numRpr = numRreps = Rreps.size();

		Bjj.assign_zero();

		HI.assign(inverse(G));
		A.assign_one();

		power_real(B, G, ht);

		if (crank == 0) {
			// check whether G^ht = (1)
			if (Gstep.is_one()) {
				FS = prin_list.search(qi_class_real(B));
				if (FS) {
					htfact.assign(ht);
					htfact.factor();
					Bjj = G.order_mult(ht, htfact);
					if (tFI.is_zero()) {
						Bmat.resize(1, 1);
						Bmat.sto(0, 0, Bjj);
					}
				}
			}
			else {
				FIDL.assign(B, 0.0);
				temp.assign_zero();
				while ((FIDL.get_distance() <= R) && (Bjj.is_zero())) {
					FS = prin_list.search(FIDL);
					if (FS) {
						htfact.assign(ht);
						htfact.factor();
						Bjj = G.order_mult(ht, htfact);
						if (tFI.is_zero()) {
							Bmat.resize(1, 1);
							Bmat.sto(0, 0, Bjj);
						}
					}
					else {
						add(temp, temp, GStepWidth);
						multiply_real(FIDL, FIDL, Gstep);
						FIDL.adjust_pos(temp);
					}
				}
			}
		}
		else {
			// check whether G^ht is in currently generated subgroup
			r = -1;
			for (i = 0; i < numQ; ++i) {
				D.assign(QT[i].get_A());
				multiply_real(E, D, B);

				if (Gstep.is_one()) {
					FS = prin_list.search(qi_class_real(E));
					if (FS)
						r = 0;
					else {
						Inode = RT.search(ideal_node(E, 0));
						if (Inode)
							r = Inode->get_index();
						else
							r = -1;
					}
				}
				else {
					FIDL.assign(E, 0.0);
					temp.assign_zero();
					while ((FIDL.get_distance() <= R) && (r < 0)) {
						FS = prin_list.search(FIDL);
						if (FS)
							r = 0;
						else {
							Inode = RT.search(ideal_node(qi_class(FIDL), 0));
							if (Inode)
								r = Inode->get_index();
							else {
								add(temp, temp, GStepWidth);
								multiply_real(FIDL, FIDL, Gstep);
								FIDL.adjust_pos(temp);
							}
						}
					}
				}

				if (r >= 0) {
					// determine smallest divisor of ht
					q = i;
					Bjj = exact_power_real(ht, G, RT, QT, numRpr, r, q, Gstep);
					if ((!Bjj.is_one()) && (tFI.is_zero()))
						decode_vector(Bmat, Bjj, r, q, Rvec, Qvec, numRpr, numQ);
					break;
				}
			}
		}

		if (Bjj.is_zero()) {
			C.assign(G);
			CINV = inverse(C);
			Bcomp.assign(B);
		}
		else
			u = 0;

		while (Bjj.is_zero()) {
			// compute more baby steps
			for (r = s; r <= u; ++r) {
				multiply_real(A, A, HI);
				for (k = 0; k < numRreps; ++k) {
					D.assign(Rreps[k]);
					multiply_real(E, D, A);
					Rreps[curr_index] = E;
					RT.hash(ideal_node(E, curr_index));

					if (Gstep.is_one()) {
						// store each cycle of ideals
						E.make_cycle(curr_index, RT);
					}
					else {
						// store ideals with distance < sqReg
						Areal.assign(E);
						Areal.make_list(sqReg, curr_index, RT);
					}

					++curr_index;
				}
			}

			// compute giant steps to u^2
			square(usqr, u);
			while ((y.compare(usqr) < 0) && (Bjj.is_zero())) {
				// search
				for (i = 0; i < numQ; ++i) {
					E.assign(QT[i].get_A());
					multiply_real(D, E, Bcomp);
					r = -1;

					if (Gstep.is_one()) {
						FS = prin_list.search(qi_class_real(D));
						if (FS)
							r = 0;
						else {
							Inode = RT.search(ideal_node(D, 0));
							if (Inode)
								r = Inode->get_index();
						}
					}
					else {
						FIDL.assign(D, 0.0);
						temp.assign_zero();
						while ((FIDL.get_distance() <= R) && (r < 0)) {
							FS = prin_list.search(FIDL);
							if (FS)
								r = 0;
							else {
								Inode = RT.search(ideal_node(qi_class(FIDL), 0));
								if (Inode)
									r = Inode->get_index();
								else {
									add(temp, temp, GStepWidth);
									multiply_real(FIDL, FIDL, Gstep);
									FIDL.adjust_pos(temp);
								}
							}
						}
					}

					if (r >= 0) {
						subtract(Bjj, ht, y);
						add(Bjj, Bjj, (r/numRpr));
						if (!Bjj.is_zero()) {
							// determine smallest divisor of Bjj
							ht.assign(Bjj);
							if (crank == 0) {
								htfact.assign(Bjj);
								htfact.factor();
								Bjj = G.order_mult(Bjj, htfact);
								if (tFI.is_zero()) {
									Bmat.resize(1, 1);
									Bmat.sto(0, 0, Bjj);
								}
							}
							else {
								q = i;
								r %= numRpr;
								Bjj = exact_power_real(Bjj, G, RT, QT, numRpr, r, q, Gstep);
								if ((!Bjj.is_one()) && (tFI.is_zero()))
									decode_vector(Bmat, Bjj, r, q, Rvec, Qvec, numRpr, numQ);
							}
							break;
						}
					}

					multiply_real(D, E, B);
					r = -1;

					if (Gstep.is_one()) {
						FS = prin_list.search(qi_class_real(D));
						if (FS)
							r = 0;
						else {
							Inode = RT.search(ideal_node(D, 0));
							if (Inode)
								r = Inode->get_index();
						}
					}
					else {
						FIDL.assign(D, 0.0);
						temp.assign_zero();
						while ((FIDL.get_distance() <= R) && (r < 0)) {
							FS = prin_list.search(FIDL);
							if (FS)
								r = 0;
							else {
								Inode = RT.search(ideal_node(qi_class(FIDL), 0));
								if (Inode)
									r = Inode->get_index();
								else {
									add(temp, temp, GStepWidth);
									multiply_real(FIDL, FIDL, Gstep);
									FIDL.adjust_pos(temp);
								}
							}
						}
					}

					if (r >= 0) {
						add(Bjj, ht, y);
						add(Bjj, Bjj, (r/numRpr));

						// determine smallest divisor of Bjj
						ht.assign(Bjj);
						if (crank == 0) {
							htfact.assign(Bjj);
							htfact.factor();
							Bjj = G.order_mult(Bjj, htfact);
							if (tFI.is_zero()) {
								Bmat.resize(1, 1);
								Bmat.sto(0, 0, Bjj);
							}
						}
						else {
							q = i;
							r %= numRpr;
							Bjj = exact_power_real(Bjj, G, RT, QT, numRpr, r, q, Gstep);
							if ((!Bjj.is_one()) && (tFI.is_zero()))
								decode_vector(Bmat, Bjj, r, q, Rvec, Qvec, numRpr, numQ);
						}
						break;
					}
				}

				if (Bjj.is_zero()) {
					// not found, take another giant step
					add(y, y, u);
					multiply_real(B, B, C);
					multiply_real(Bcomp, Bcomp, CINV);
				}
			}

			// double u
			if (Bjj.is_zero()) {
				s = u+1;
				u <<= 1;
				square_real(C, C);
				CINV = inverse(C);
			}
		}

		if (!Bjj.is_one()) {
			h1 *= Bjj;
			divide(ht, ht, Bjj);
			++crank;
		}

		if (bigfloat(h1) < hstar) {
			temp = ceil(sqrt(bigfloat(Bjj)));
			temp.longify(Bj);

			// compute new R' (remove entries with too large exponents)
			det *= Bj;
			det.longify(upper);
			idx = static_cast<lidia_size_t>(upper);
			if (idx < curr_index) {
				numRpr = RT.no_of_elements();
				// removing entries
				i = numRpr-1;
				while (RT[i].get_index() >= idx) {
					RT.remove_from(i);
					--i;
				}
				Rreps.set_size(static_cast<lidia_size_t>(idx));
			}
			else {
				// adding entries
				idx = static_cast<lidia_size_t>(u)+1;
				for (r = idx; r < Bj; ++r) {
					multiply_real(A, A, HI);
					for (k = 0; k < numRreps; ++k) {
						D.assign(Rreps[k]);
						multiply_real(E, D, A);
						Rreps[curr_index] = E;
						RT.hash(ideal_node(E, curr_index));

						if (Gstep.is_one()) {
							// store each cycle of ideals
							E.make_cycle(curr_index, RT);
						}
						else {
							// store ideals with distance < sqReg
							Areal.assign(E);
							Areal.make_list(sqReg, curr_index, RT);
						}

						++curr_index;
					}
				}
			}


			if (!Bjj.is_one()) {
				Rvec[crank-1] = Bj;
				Qvec[crank-1] = 1;

				// compute new Q
				numQ = QT.no_of_elements();
				curr_index = numQ;
				power_real(GBj, G, Bj);
				Gq.assign(GBj);
				for (i = 1; i < Bj; ++i) {
					if (Bjj.compare(i*Bj) <= 0)
						break;
					for (k = 0; k < numQ; ++k) {
						E.assign(QT[k].get_A());
						multiply_real(D, E, Gq);
						QT.hash(ideal_node(D, curr_index));
						++curr_index;
					}
					multiply_real(Gq, Gq, GBj);
					++Qvec[crank-1];
				}
			}

			++pidx;
		}
	}

//  RHNF.assign(Bmat);
	if (tFI.is_zero()) {
		// compute structure
		h.assign(h1);
		if (h == 1)
			CL[0] = 1;
		else {
			Bmat.snf_havas(U, junk);
			i = 0;
			for (j = 0; j < crank; ++j) {
				Bjj.assign(Bmat.member(j, j));
				if (!Bjj.is_one())
					CL[i++] = Bjj;
			}
		}
	}
}



//
// FIX
//

//
// quadratic_order::cgh_imag
//
// Task:
//      computes the class group given the class number.
//

void
quadratic_order::cgh_imag()
{
	debug_handler("quadratic_order", "cgh_imag");
}



//
// FIX
//

//
// quadratic_order::cgh_real
//
// Task:
//      computes the class group given the class number.
//

void
quadratic_order::cgh_real()
{
	debug_handler("quadratic_order", "cgh_real");
}



//
// quadratic_order::factor_imag
//
// Task:
//      factors the discriminant using the 2-Sylow subgroup
//

void
quadratic_order::factor_imag()
{
	debug_handler("quadratic_order", "factor_imag");

	lidia_size_t i, j, numgens, numfacts;
	long pwr;
	bigint curr_fact, val;
	qi_class A;
	base_vector< qi_class > S2gens;
	rational_factorization new_comp;
	sort_vector< bigint > facts;

	// if class group has not been computed, factor using integer fact. alg.
	if (!is_CL_computed()) {
		disc_fact.assign(Delta);
		disc_fact.factor();
	}


	if (disc_fact != Delta) {
		facts.set_mode(EXPAND);
		facts.reset();
		S2gens.set_mode(EXPAND);
		S2gens.reset();
		disc_fact.assign(-1);

		// compute generators of 2-Sylow subgroup
		generators();
		numgens = gens.size();
		j = 0;
		for (i = 0; i < numgens; ++i) {
			if (CL[i].is_even())
				power(S2gens[j++], gens[i], CL[i]/2);
		}
		numgens = j;

		if (numgens == 0)
			// Delta is prime
			disc_fact.assign(Delta);
		else {
			// compute divisors of the discriminant
			numfacts = 0;
			if (Delta.is_even())
				facts[numfacts++] = bigint(2);
			for (i = 0; i < numgens; ++i) {
				A = S2gens[i];
				if (A.get_b().is_zero()) {
					// Delta = -4ac
					facts[numfacts++] = abs(A.get_a());
					facts[numfacts++] = abs(A.get_c());
				}
				else if (A.get_a() == A.get_b()) {
					// Delta = b(b-4c)
					facts[numfacts++] = abs(A.get_b());
					facts[numfacts++] = abs(A.get_b() - (A.get_c() << 2));
				}
				else {
					// Delta = (b-2a)(b+2a)
					facts[numfacts++] = abs(A.get_b() - (A.get_a() << 1));
					facts[numfacts++] = abs(A.get_b() + (A.get_a() << 1));
				}
			}

			// refine divisors
			facts.sort();
			for (i = 0; i < numfacts-1; ++i) {
				curr_fact = facts[i];
				if (!curr_fact.is_one()) {
					for (j = i+1; j < numfacts; ++j)
						while ((facts[j] % curr_fact) == 0)
							facts[j] /= curr_fact;
				}
			}

			// construct factorization
			val.assign(-Delta);
			for (i = 0; i < numfacts; ++i) {
				curr_fact = facts[i];
				if (curr_fact > 1) {
					new_comp.assign(curr_fact);

					pwr = 0;
					while ((val % curr_fact) == 0) {
						++pwr;
						val /= curr_fact;
					}

					for (j = 0; j < pwr; ++j)
						disc_fact = disc_fact*new_comp;
				}
			}
			if (!val.is_one()) {
				new_comp.assign(val);
				disc_fact = disc_fact*new_comp;
			}
		}
	}


	if (!disc_fact.is_prime_factorization())
		disc_fact.factor();


	if (disc_fact != Delta) {
		std::cout << "Delta = " << Delta << std::endl;
		std::cout << disc_fact << std::endl;
		lidia_error_handler("quadratic_order", "factor_imag() - factorization not "
				    "equal to the discriminant.  Please report!");
		return;
	}
}



//
//
// quadratic_order::factor_real
//
// Task:
//      factors the discriminant using the 2-Sylow subgroup
//

void
quadratic_order::factor_real()
{
	debug_handler("quadratic_order", "factor_real");

	lidia_size_t i, j, numgens, numfacts;
	long pwr;
	bigint curr_fact, val;
	qi_class_real A;
	base_vector< qi_class_real > S2gens;
	rational_factorization new_comp;
	sort_vector< bigint > facts;

	// if class group has not been computed, factor using integer fact. alg.
	if (!is_CL_computed()) {
		disc_fact.assign(Delta);
		disc_fact.factor();
	}


	if (disc_fact != Delta) {
		facts.set_mode(EXPAND);
		facts.reset();
		S2gens.set_mode(EXPAND);
		S2gens.reset();
		disc_fact.assign(1);

		numfacts = 0;

		// test whethere there is an ambiguous principal ideal (at R/2)
		A.assign_one();
		find_ambiguous(A, facts, numfacts);

		// compute generators of 2-Sylow subgroup
		if (h.is_even()) {
			generators();
			numgens = gens.size();
			j = 0;
			for (i = 0; i < numgens; ++i) {
				if (CL[i].is_even())
					power(S2gens[j++], qi_class_real(gens[i]), CL[i]/2);
			}
			numgens = j;

			if (Delta.is_even())
				facts[numfacts++] = bigint(2);

			// Don't do this until principality testing is faster!!!
			// find ambiguous ideals in each ambiguous class

//  			for (i = 0; i < numgens; ++i) {
//  				A = S2gens[i];
//  				find_ambiguous(A, facts, numfacts);
//  			}
		}

		// refine divisors
		facts.sort();
		for (i = 0; i < numfacts-1; ++i) {
			curr_fact = facts[i];
			if (!curr_fact.is_one()) {
				for (j = i+1; j < numfacts; ++j)
					while ((facts[j] % curr_fact) == 0)
						facts[j] /= curr_fact;
			}
		}

		// construct factorization
		val.assign(Delta);
		for (i = 0; i < numfacts; ++i) {
			curr_fact = facts[i];
			if (curr_fact > 1) {
				new_comp.assign(curr_fact);

				pwr = 0;
				while ((val % curr_fact) == 0) {
					++pwr;
					val /= curr_fact;
				}

				for (j = 0; j < pwr; ++j)
					disc_fact = disc_fact*new_comp;
			}
		}
		if (!val.is_one()) {
			new_comp.assign(val);
			disc_fact = disc_fact*new_comp;
		}

		if (info > 4)
			std::cout << "BEFORE:  " << disc_fact << std::endl;

		if (!disc_fact.is_prime_factorization())
			disc_fact.factor();

		if (info > 4)
			std::cout << "AFTER:   " << disc_fact << std::endl;

		if (disc_fact != Delta) {
			std::cout << "\nDelta = " << Delta << std::endl;
			std::cout << disc_fact << std::endl;
			lidia_error_handler("quadratic_order", "factor_real() - factorization not "
					    "equal to the discriminant.  Please report!");
			return;
		}
	}
}


//
// quadratic_order::find_ambiguous
//
// Task:
//      factors the discriminant using the 2-Sylow subgroup
//

void
quadratic_order::find_ambiguous(qi_class_real & AMB,
                                sort_vector< bigint > & facts,
                                int & numfacts)
{
	debug_handler("quadratic_order", "find_ambiguous");

	qi_class_real A, B, AINV;
	bigfloat R2, dist;
	bool found;

	R2.assign(regulator());
	R2.divide_by_2();

	if (info > 4) {
		std::cout << "\nIn find_ambiguous:  AMB = " << AMB << std::endl;
		std::cout << "R/2 = " << R2 << std::endl;
	}

	if (AMB.is_one()) {
		// just check at distance R/2
		A.assign_one();
		A.assign(nearest(A, R2));
		if (info > 4)
			std::cout << "at R/2:  " << A << std::endl;

		// test A,rho(A), and inverse_rho(A)
		found = false;
		if (((Delta % A.get_a()) == 0) && (!A.is_one())) {
			facts[numfacts++] = A.get_a();
			facts[numfacts++] = Delta / A.get_a();
			if (info > 4)
				std::cout << "factor:  " << A.get_a() << std::endl;
			found = true;
		}

		if (!found) {
			apply_rho(B, A);
			if (info > 4)
				std::cout << "rho:  " << B << std::endl;
			if (((Delta % B.get_a()) == 0) && (!B.is_one())) {
				facts[numfacts++] = B.get_a();
				facts[numfacts++] = Delta / B.get_a();
				if (info > 4)
					std::cout << "factor:  " << A.get_a() << std::endl;
				found = true;
			}
		}

		if (!found) {
			apply_inverse_rho(B, A);
			if (info > 4)
				std::cout << "rho^-1:  " << B << std::endl;
			if (((Delta % B.get_a()) == 0) && (!B.is_one())) {
				facts[numfacts++] = B.get_a();
				facts[numfacts++] = Delta / B.get_a();
				if (info > 4)
					std::cout << "factor:  " << A.get_a() << std::endl;
				found = true;
			}
		}
	}
	else {
		// not principal class - find ambiguous ideal halfway between A and A^-1
		inverse(AINV, AMB);
		AINV.is_equivalent(AMB, dist);
		if (info > 4)
			std::cout << "d(A, A^-1) = " << dist << std::endl;

		dist.divide_by_2();
		A.assign(nearest(AMB, dist));
		if (info > 4)
			std::cout << "at dist/2 = " << A << std::endl;

		// test A,rho(A), and inverse_rho(A)
		found = false;
		if (((Delta % A.get_a()) == 0) && (!A.is_one())) {
			facts[numfacts++] = A.get_a();
			facts[numfacts++] = Delta / A.get_a();
			if (info > 4)
				std::cout << "factor:  " << A.get_a() << std::endl;
			found = true;
		}

		if (!found) {
			apply_rho(B, A);
			if (info > 4)
				std::cout << "rho:  " << B << std::endl;
			if (((Delta % B.get_a()) == 0) && (!B.is_one())) {
				facts[numfacts++] = B.get_a();
				facts[numfacts++] = Delta / B.get_a();
				if (info > 4)
					std::cout << "factor:  " << A.get_a() << std::endl;
				found = true;
			}
		}

		if (!found) {
			apply_inverse_rho(B, A);
			if (info > 4)
				std::cout << "rho^-1:  " << B << std::endl;
			if (((Delta % B.get_a()) == 0) && (!B.is_one())) {
				facts[numfacts++] = B.get_a();
				facts[numfacts++] = Delta / B.get_a();
				if (info > 4)
					std::cout << "factor:  " << A.get_a() << std::endl;
				found = true;
			}
		}


		// test ideal at distance dist/2 + R/2 from AMB
		add(dist, dist, R2);
		A.assign(nearest(AMB, dist));
		if (info > 4)
			std::cout << "at dist/2 + R/2 = " << A << std::endl;

		// test A,rho(A), and inverse_rho(A)
		found = false;
		if (((Delta % A.get_a()) == 0) && (!A.is_one())) {
			facts[numfacts++] = A.get_a();
			facts[numfacts++] = Delta / A.get_a();
			if (info > 4)
				std::cout << "factor:  " << A.get_a() << std::endl;
			found = true;
		}

		if (!found) {
			apply_rho(B, A);
			if (info > 4)
				std::cout << "rho:  " << B << std::endl;
			if (((Delta % B.get_a()) == 0) && (!B.is_one())) {
				facts[numfacts++] = B.get_a();
				facts[numfacts++] = Delta / B.get_a();
				if (info > 4)
					std::cout << "factor:  " << A.get_a() << std::endl;
				found = true;
			}
		}

		if (!found) {
			apply_inverse_rho(B, A);
			if (info > 4)
				std::cout << "rho^-1:  " << B << std::endl;
			if (((Delta % B.get_a()) == 0) && (!B.is_one())) {
				facts[numfacts++] = B.get_a();
				facts[numfacts++] = Delta / B.get_a();
				if (info > 4)
					std::cout << "factor:  " << A.get_a() << std::endl;
				found = true;
			}
		}
	}
}



//
// quadratic_order::exact_power
//
// Task:
//      computes the exact power of G that is in a subgroup, given a multiple.
//

bigint
quadratic_order::exact_power(const bigint & omult, const qi_class & G,
                             indexed_hash_table< ideal_node > & RT,
                             indexed_hash_table< ideal_node > & QT,
                             lidia_size_t numRpr, long & r, long & q)
{
	debug_handler("quadratic_order", "exact_power");

	lidia_size_t i, j, num, numQ;
	long oldr, oldq;
	bigint ord, p, pwr, ex, oldmult;
	qi_class A, D, E;
	rational_factorization pfact;
	bool found_rel;
	ideal_node *Inode;

	oldr = r;
	oldq = q;
	oldmult = omult;

	numQ = QT.no_of_elements();
	pfact.assign(omult);
	pfact.factor();

	ord.assign_one();
	num = pfact.no_of_comp();
	for (i = 0; i < num; ++i) {
		A.assign(G);
		p.assign(pfact.base(i));
		ex = pfact.exponent(i);

		pwr.assign(omult);
		for (j = 0; j < ex; ++j)
			divide(pwr, pwr, p);

		if (pwr > 1)
			power_imag(A, A, pwr);

		found_rel = false;
		while (!found_rel) {
			for (j = 0; j < numQ; ++j) {
				r = -1;
				D.assign(QT[j].get_A());
				multiply_imag(E, D, A);

				if (E.is_one())
					r = 0;
				else {
					Inode = RT.search(ideal_node(E, 0));
					if (Inode)
						r = Inode->get_index();
				}
				if ((r >= 0) && (r < numRpr)) {
					found_rel = true;
					break;
				}
			}
			if (!found_rel) {
				multiply(ord, ord, p);
				power_imag(A, A, p);
			}
		}
	}

	// ord is the smallest power of G contained in the current subgroup
	if (ord == oldmult) {
		r = oldr;
		q = oldq;
	}
	else {
		power(A, G, ord);
		r = -1;
		for (j = 0; j < numQ; ++j) {
			D.assign(QT[j].get_A());
			multiply_imag(E, D, A);
			if (E.is_one())
				r = 0;
			else {
				Inode = RT.search(ideal_node(E, 0));
				if (Inode)
					r = Inode->get_index();
			}
			if ((r >= 0) && (r < numRpr)) {
				q = j;
				break;
			}
			else
				r = -1;
		}
	}

	return ord;
}



//
// quadratic_order::exact_power_real
//
// Task:
//      computes the exact power of G that is in a subgroup, given a multiple.
//

bigint
quadratic_order::exact_power_real(const bigint & omult, const qi_class & G,
                                  indexed_hash_table< ideal_node > & RT,
                                  indexed_hash_table< ideal_node > & QT,
                                  lidia_size_t numRpr, long & r, long & q,
                                  qi_class_real & Gstep)
{
	debug_handler("quadratic_order", "exact_power_real");

	lidia_size_t i, j, num, numQ;
	long oldr, oldq;
	bigint ord, p, pwr, ex, oldmult;
	qi_class A, D, E;
	qi_class_real F, *FS;
	rational_factorization pfact;
	bool found_rel;
	ideal_node *Inode;
	bigfloat temp, GStepWidth;

	oldr = r;
	oldq = q;
	oldmult = omult;
	floor(GStepWidth, Gstep.get_distance());

	numQ = QT.no_of_elements();
	pfact.assign(omult);
	pfact.factor();

	ord.assign_one();
	num = pfact.no_of_comp();
	for (i = 0; i < num; ++i) {
		A.assign(G);
		p.assign(pfact.base(i));
		ex = pfact.exponent(i);

		pwr.assign(omult);
		for (j = 0; j < ex; ++j)
			divide(pwr, pwr, p);

		if (pwr > 1)
			power_real(A, A, pwr);

		found_rel = false;
		while (!found_rel) {
			for (j = 0; j < numQ; ++j) {
				r = -1;
				D.assign(QT[j].get_A());
				multiply_real(E, D, A);

				if (Gstep.is_one()) {
					FS = prin_list.search(qi_class_real(E));
					if (FS)
						r = 0;
					else {
						Inode = RT.search(ideal_node(E, 0));
						if (Inode)
							r = Inode->get_index();
					}
				}
				else {
					F.assign(E, 0.0);
					temp.assign_zero();
					while ((F.get_distance() <= R) && (r < 0)) {
						FS = prin_list.search(F);
						if (FS)
							r = 0;
						else {
							Inode = RT.search(ideal_node(qi_class(F), 0));
							if (Inode)
								r = Inode->get_index();
							else {
								add(temp, temp, GStepWidth);
								multiply_real(F, F, Gstep);
								F.adjust_pos(temp);
							}
						}
					}
				}

				if ((r >= 0) && (r < numRpr)) {
					found_rel = true;
					break;
				}
			}

			if (!found_rel) {
				multiply(ord, ord, p);
				power_real(A, A, p);
			}
		}
	}

	// ord is the smallest power of G contained in the current subgroup
	if (ord == oldmult) {
		r = oldr;
		q = oldq;
	}
	else {
		power(A, G, ord);
		r = -1;
		for (j = 0; j < numQ; ++j) {
			D.assign(QT[j].get_A());
			multiply_real(E, D, A);

			if (Gstep.is_one()) {
				FS = prin_list.search(qi_class_real(E));
				if (FS)
					r = 0;
				else {
					Inode = RT.search(ideal_node(E, 0));
					if (Inode)
						r = Inode->get_index();
				}
			}
			else {
				F.assign(E, 0.0);
				temp.assign_zero();
				while ((F.get_distance() <= R) && (r < 0)) {
					FS = prin_list.search(F);
					if (FS)
						r = 0;
					else {
						Inode = RT.search(ideal_node(qi_class(F), 0));
						if (Inode)
							r = Inode->get_index();
						else {
							add(temp, temp, GStepWidth);
							multiply_real(F, F, Gstep);
							F.adjust_pos(temp);
						}
					}
				}
			}

			if ((r >= 0) && (r < numRpr)) {
				q = j;
				break;
			}
			else
				r = -1;
		}
	}

	return ord;
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
