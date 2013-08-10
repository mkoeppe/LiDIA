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
#include	"LiDIA/bigrational.h"
#include	"LiDIA/bigfloat_matrix.h"
#include	"LiDIA/bigfloat_lattice.h"
#include	"LiDIA/bigint_lattice.h"



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



void quad_erg(math_matrix< bigrational > &, const math_matrix< bigint > &);
void diag_min(bigint &, const math_matrix< bigint > &);
long special_diag_min(base_matrix< bigrational > &);
void sort_columns(base_matrix< bigint > &, base_matrix< bigfloat > &);
void auszaehlen_jurk(bigrational &, math_matrix< bigint > &, const math_matrix< bigint > &, const math_matrix< bigrational > &, int alle);
void cho_from_quad(base_matrix< bigfloat > &, const base_matrix< bigrational > &);
int is_null(const math_vector< bigint > &);



#define VGET_P(a, b, t) (a = (b).get_data_address())
#define VGET_L(a, b, t) (a = (b).size())
#define VGET_C(a, b, t) (a = (b).capacity())

#define MGET_R(a, b, t) (a = (b).get_no_of_rows())
#define MGET_C(a, b, t) (a = (b).get_no_of_columns())
#define MGET_P(a, b, t) (a = (b).get_data_address())



// #endif

// Algorithmus 1 berechnet die Matrix Q fuer die quadratisch Form
// aus einer positiv definiten Matrix A (obere Dreiecksmatrix)


void
quad_erg(math_matrix< bigrational > &Q, const math_matrix< bigint > &A)
{
	debug_handler("bigint_lattice", "in function quad_erg(Q, A)");

	lidia_size_t n;
	MGET_R(n, Q, bigrational);

	if (n != Q.get_no_of_columns() || n != A.get_no_of_rows() || n != A.get_no_of_columns())
		lidia_error_handler("bigint_lattice", "function quad_erg(Q, A) - "
				    "Matrices have incorrect number of rows or columns");

	lidia_size_t i, j, l;
	bigrational diag, **q;
	MGET_P(q, Q, bigrational);

	for (i = 0; i < n; i++) {
		for (j = i; j < n; j++) {
			q[i][j].assign(A.member(i, j));
		}
	}

	for (i = 0; i < n; i++) {
		diag.assign(q[i][i]); //diag = Q[i][i]

		for (j = i + 1; j < n; j++) {
			q[j][i].assign(q[i][j]); //Q[j][i] = Q[i][j]
			divide(q[i][j], q[i][j], diag); //Q[i][j] = Q[i][j] / Q[i][i]
		}

		for (j = i + 1; j < n; j++) {
			for (l = j; l < n; l++) {
				multiply(diag, q[j][i], q[i][l]);
				subtract(q[j][l], q[j][l], diag); //Q[j][l] = Q[j][l] - Q[j][i] * Q[i][l]
			}
		}
	}
}



void
diag_min(bigint &erg, const math_matrix< bigint > &A)
{
	debug_handler("bigint_lattice", "in function diag_min(erg, A)");

	long j, rows, columns, tmp = 0;
	bigint **a;
	MGET_R(rows, A, bigint);
	MGET_C(columns, A, bigint);
	MGET_P(a, A, bigint);

	if (rows < columns) {
		j = rows;
	}
	else {
		j = columns;
	}

	for (register long i = 1; i < j; i++) {
		if (a[tmp][tmp].compare(a[i][i]) > 0) {
			tmp = i;
		}
	}
	erg.assign(a[tmp][tmp]);
}



// kleinstes l bestimmen, so dass A[0][0] <= A[j][j] fuer alle j
// mit l <= j <= k, wobei l >= 1
long
special_diag_min(base_matrix< bigrational > &A)
{
	debug_handler("bigint_lattice", "in function special_diag_min(A)");

	long j, rows;
	bigrational **a;
	MGET_R(rows, A, bigrational);
	MGET_P(a, A, bigrational);

	j = rows;

	while (--j > 0) {
		if (a[0][0].compare(a[j][j]) >= 0) {
			return (j+1);
		}
	}
	return 1;
}



void
cho_from_quad(base_matrix< bigfloat > &R, const base_matrix< bigrational > &Q)
{
	long n;
	MGET_R(n, Q, bigrational);

	if (n != Q.get_no_of_columns() || n != R.get_no_of_rows() || n != R.get_no_of_columns()) {
		lidia_error_handler("bigint_lattice", "function cho_from_quad(R, Q) - "
				    "Matrices have incorrect number of rows or columns");
	}
	long i, j;
	bigfloat **r;
	bigrational **q;

	MGET_P(q, Q, bigrational);
	MGET_P(r, R, bigfloat);

	for (i = 0; i < n; i++) {
		for (j = 0; j < i; j++) {
			r[i][j].assign_zero();
		}

		sqrt(r[i][i], q[i][i]);

		for (j = i+1; j < n; j++) {
			multiply(r[i][j], r[i][i], q[i][j]);
		}
	}
}



//B sei Nullmatrix
void
sort_columns(base_matrix< bigint > &B, base_matrix< bigfloat > &R)
{
	debug_handler("bigint_lattice", "in function sort_columns(R)");

	const long columns = R.get_no_of_columns();

	B.resize(columns, columns); // Anpassen von Zeilen und Spalten, falls noetig

	struct sort {
		long pos;
		bigfloat lenght;
	};

	sort *A = new sort[columns]; // Anzahl der Spalten
	bigfloat len;
	math_vector< bigfloat > tmp1;
	tmp1.set_capacity(columns);

	long i = 0, j, akt;

	while (i < columns) {
		j = i;
		akt = i-1;
		R.get_column_vector(tmp1, i);
		len = tmp1.sum_of_squares(); // PT

		while (akt >= 0) {
			if ((A[akt].lenght).compare(len) < 0) {
				(A[j].lenght).assign(A[akt].lenght);
				A[j].pos = A[akt].pos;
				akt--;
				j--;
			}
			else {
				break;
			}
		}
		A[j].pos = i;
		(A[j].lenght).assign(len);
		i++;
	}
	for (i = 0; i < columns; i++) {
		B.sto(A[i].pos, i, 1);
	}
        delete[] A;
}



void
auszaehlen_jurk(bigrational &c, math_matrix< bigint > &RES,
		const math_matrix< bigint > &C, const math_matrix< bigrational > &Q, int alle)
{
	debug_handler("bigint_lattice", "auszaehlen_jurk(c, RES, C, Q, alle)");
	lidia_size_t k, cr;
	MGET_R(k, Q, bigrational);
	MGET_R(cr, C, bigint);
	long i = k-1, erg, j, anzahl = 0;
	int  flag = 1;
	bigint ganz, ganz2, *O, *Z, *Y;
	bigfloat z, tmpfloat;
	bigrational helpf, tmp, length, **q, *S, *T;
	math_vector< bigrational > Stmp(k, k), Ttmp(k, k);
	math_vector< bigint > Otmp(k, k), Ztmp(k, k), Ytmp(cr, cr);

	VGET_P(S, Stmp, bigrational);
	VGET_P(T, Ttmp, bigrational);
	VGET_P(O, Otmp, bigint);
	VGET_P(Z, Ztmp, bigint);
	VGET_P(Y, Ytmp, bigint);

	RES.set_no_of_columns(cr);

	MGET_P(q, Q, bigrational);

	length.assign(c);

	do {             // 1.do-while-Schleife
		do {             // 2.do-while-Schleife
			if (flag) {       // flag gleich 1: (2) fortfahren
				subtract(helpf, length, T[i]);
				if (helpf.is_lt_zero()) {
					i++;
				}
				else {
					divide(tmp, helpf, q[i][i]);
					z.assign(tmp);
					sqrt(tmpfloat, z);
					ceil(z, tmpfloat);
					z.bigintify(ganz);
					ceil(ganz2, bigfloat(S[i])); // MM
					add(ganz2, ganz, ganz2);
					inc(ganz2);
					negate(Ztmp[i], ganz2); // untere Schranke der i-ten Komponente ist erster Wert

					subtract(Otmp[i], ganz, floor(S[i])); // obere Schranke der i-ten Komponente
				} // else
			}// if (flag)

			do {       // 3.do-while-Schleife
				inc(Ztmp[i]);
				if (O[i].compare(Z[i]) >= 0) {
					break;
				}
				else {
					i++;
				}
			} while (true); // Ende der 3.do-while-Schleife

			if (i == 0) {
				break; // (7)-->(10)
			}

			flag = 1;
			erg = i--; // naechstkleinere Komponente

			tmp.assign_zero();

			for (j = erg; j < k; j++) {      // S berechnen
				multiply(helpf, q[i][j], Z[j]);
				add(tmp, tmp, helpf);
			}
			Stmp[i].assign(tmp);

			add(helpf, Z[erg], S[erg]); // T[i] berechnen
			square(helpf, helpf);
			multiply(helpf, helpf, q[erg][erg]);
			add(Ttmp[i], helpf, T[erg]);
		} while (true); // Ende der 2.do-while-Schleife

		if (is_null(Ztmp)) {    // (10) Ausgabe + Ende
			c.assign(length);
			return;
		}
		else {         // kein Nullvektor der Koeffizienten
			add(helpf, Z[0], S[0]); // T berechnen
			square(helpf, helpf);

			multiply(helpf, helpf, q[0][0]);
			add(helpf, helpf, T[0]);

			erg = helpf.compare(length);
			if (erg < 1) {
				multiply(Ytmp, C, Ztmp); // Loesung transformieren
				if (erg < 0 && !alle) {    // neuer Vektor kleiner und nicht alle kleineren Vektoren gesucht
					length.assign(helpf); // neue kuerzeste Laenge
					anzahl = 1;
					RES.set_no_of_rows(anzahl);
					RES.sto_row_vector(Ytmp, Ytmp.size(), 0);
				}
				else {   // erg == 0: weitere Loesung bei gleicher Laenge
					RES.set_no_of_rows(anzahl+1);
					RES.sto_row_vector(Ytmp, Ytmp.size(), anzahl);
					anzahl++;
				}
			}//if (erg < 1)
		} //else
		flag = 0; // flag gleich 0: (5) fortfahren
	} while (true); // Ende der 1.do-while-Schleife
}



int
is_null (const math_vector< bigint > &V)
{
	lidia_size_t l;
	bigint *v;
	VGET_L(l, V, bigint);
	VGET_P(v, V, bigint);

	while (--l >= 0) {
		if (!(v[l].is_zero())) {
			return 0;
		}
	}

	return 1;
}



lidia_size_t
bigint_lattice::shortest_vector (math_matrix< bigint > & RES, lidia_size_t prec)
{
	debug_handler("bigint_lattice", "shortest_vector(RES, T, prec)");
	if (check_basis()) {
		long dim, old_prec, i, j;
		int l;
		bool mtrans, gf;

		bigint len_int;
		bigrational len_rat;

		math_matrix< bigint > A(columns, columns), T(columns, columns);
		math_matrix< bigint > B(rows, columns), M(*this);
		math_matrix< bigrational > Q(columns, columns);

		mtrans = chk_trans();
		if (mtrans) {
			M.trans();
		}

		gf = check_gram();
		if (gf) {
			for (i = 0; i < columns; i++) {
				T.sto(i, i, 1);
			}

			A.assign(M);
		}
		else {
			// LLL-Reduktion
			bigint_lattice L(M);
			L.lll(T, 0.99);

			LiDIA::multiply(B, M, T);
			// Erstellen der quadratischen Form
			LiDIA::multiply(A, (ring_matrix< bigint > ) LiDIA::trans(B), B);
		}

		diag_min(len_int, A);

		// Cholesky-Zerlegung
		quad_erg(Q, A);

		dim = special_diag_min(Q);

		if (dim == 1) {              // fertig
			T.set_no_of_columns(1);
			LiDIA::multiply(RES, M, T);
			if (mtrans) {
				RES.trans();
			}
			Q[0][0].intify(l);
			return l;
		}

		// weitere Deklarationen
		bigint_matrix U(dim, dim), S(dim, dim);
		math_matrix< bigfloat > R(dim, dim), INV(dim, dim);

		old_prec = bigfloat::get_precision();
		bigfloat::set_precision(prec);

		Q.set_no_of_rows(dim);
		Q.set_no_of_columns(dim);

		cho_from_quad(R, Q);
		invert(INV, R);

		LiDIA::trans(INV, INV);

		bigfloat_lattice F(INV);
		F.set_basis_flag();
		F.lll(U, 0.99); // INV ist lll-reduziert

		for (i = 0; i < dim; i++) {
			for (j = 0; j < dim; j++) {
				R.sto(i, j, U[i][j]);
			}
		}

		LiDIA::multiply(INV, INV, R);

		sort_columns(S, INV);

		LiDIA::multiply(U, U, S);

		S.adj(U); // S: Inverse von U

		LiDIA::trans(U, S);

		T.set_no_of_columns(dim);

		LiDIA::multiply(T, T, U); //T: Produkt aller Transformationsmatrizen bzgl B
		LiDIA::multiply(B, M, T);

		// Erstellen der quadratischen Form
//		LiDIA::multiply(A,::trans(B),B);
		LiDIA::multiply(A, S, B);

		// Berechnung der quadratischen Ergaenzung
		quad_erg(Q, A);

		dim = special_diag_min(Q);

		if (dim == 1) {             // fertig
			T.set_no_of_columns(1);
			LiDIA::multiply(RES, M, T);
			if (mtrans) {
				RES.trans();
			}
			Q[0][0].intify(l);
			return l;
		}

		Q.set_no_of_rows(dim);
		Q.set_no_of_columns(dim);
		T.set_no_of_columns(dim);

		// Auszaehlen
		len_rat.assign(len_int);
		auszaehlen_jurk(len_rat, RES, T, Q, 0);

		if (gf) {
			RES.resize(1, 1);
			RES.sto(0, 0, 0);
		}
		else {
			LiDIA::multiply(RES, M, (ring_matrix< bigint > ) LiDIA::trans(RES));
			if (mtrans) {
				RES.trans();
			}
		}

		len_rat.intify(l);

		bigfloat::set_precision(old_prec);

		return l;
	}
	return 1; // PT
}



lidia_size_t
shortest_vector(const bigint_lattice & L, math_matrix< bigint > & RES, lidia_size_t prec)
{
	bigint_lattice X(L);
	return X.shortest_vector(RES, prec);
}



void
bigint_lattice::short_vector_coeff(math_matrix< bigint > & RES,
				   const bigrational & limit,
				   lidia_size_t prec)
{
	debug_handler("bigint_lattice", "short_vector_coeff(RES, limit, prec)");
	if (check_basis()) {
		long dim = columns, old_prec, i, j;
		bool mtrans;

		bigint len_int;
		bigrational len_rat;

		math_matrix< bigint > A(columns, columns), T(columns, columns);
		math_matrix< bigint > B(rows, columns), M(*this);
		math_matrix< bigrational > Q(columns, columns);
		bigint_matrix U(columns, columns), S(columns, columns);
		math_matrix< bigfloat > R(columns, columns), INV(columns, columns);

		mtrans = chk_trans();
		if (mtrans)
			M.trans();

		// LLL-Reduktion
		bigint_lattice L(M);
		L.lll(T, 0.99);
		LiDIA::multiply(B, M, T);

		// Erstellen der quadratischen Form
		LiDIA::multiply(A, (ring_matrix< bigint > ) LiDIA::trans(B), B);

		// Cholesky-Zerlegung
		quad_erg(Q, A);

		old_prec = bigfloat::get_precision();
		bigfloat::set_precision(prec);

		cho_from_quad(R, Q);
		invert(INV, R);
		LiDIA::trans(INV, INV);

		bigfloat_lattice F(INV);
		F.set_basis_flag();
		F.lll(U, 0.99); // INV ist lll-reduziert

		for (i = 0; i < dim; i++) {
			for (j = 0; j < dim; j++) {
				R.sto(i, j, U[i][j]);
			}
		}
		LiDIA::multiply(INV, INV, R);
		sort_columns(S, INV);
		LiDIA::multiply(U, U, S);
		S.adj(U); // S: Inverse von U
		LiDIA::trans(U, S);
		LiDIA::multiply(T, T, U); //T: Produkt aller Transformationsmatrizen bzgl B
		LiDIA::multiply(B, M, T);

		// Erstellen der quadratischen Form
		LiDIA::multiply(A, (ring_matrix< bigint > ) LiDIA::trans(B), B);

		// Berechnung der quadratischen Ergaenzung
		quad_erg(Q, A);

		// Auszaehlen
		len_rat.assign(limit);
		auszaehlen_jurk(len_rat, RES, T, Q, 1);
		if (mtrans) {
			RES.trans();
		}
		bigfloat::set_precision(old_prec);
	}
}



void
short_vector_coeff (const bigint_lattice & L, math_matrix< bigint > & RES,
		    const bigrational & limit, lidia_size_t prec)
{
	bigint_lattice X(L);
	X.short_vector_coeff(RES, limit, prec);
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
