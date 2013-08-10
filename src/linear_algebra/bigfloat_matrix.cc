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
//	Author	: Thorsten Lauer (TL)
//	Changes	: See CVS log
//
//==============================================================================================


#ifdef HAVE_CONFIG_H
# include	"config.h"
#endif
#include	"LiDIA/bigfloat_matrix.h"
#include	"LiDIA/bigfloat_int.h"



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



void
bigfloat_matrix::
gauss_jordan(base_vector< bigfloat > & x, const base_vector< bigfloat > & b) const
{
	debug_handler("bigfloat_matrix", "gauss_jordan_pc(x, b)");
	if (rows != b.size())
		lidia_error_handler("bigfloat_lattice", "gauss_jordan_pc(x, b) :: illegal vector size");
	if (rows > columns)
		lidia_error_handler("bigfloat_lattice", "gauss_jordan_pc(x, b) :: "
				    "rows have to be less or equal columns");

	long i, j, k, r = 0, c = 0, *index = new long[columns]; // MM
	bigfloat_int factor, tmp, *help0, *help1, *helpr, *l, **p;
	bigfloat max, temp;

	x.set_capacity(rows);

	help1 = new bigfloat_int[columns];
	memory_handler(help1, "bigfloat_lattice", "gauss_jordan_pc(x, b) :: not enough memory !");

	l = new bigfloat_int[rows];
	memory_handler(l, "bigfloat_lattice", "gauss_jordan_pc(x, b) :: not enough memory !");

	for (i = 0; i < rows; i++)
		l[i].assign(b[i]);
	for (i = 0; i < columns; i++)
		index[i] = i;

	p = new bigfloat_int*[rows];
	memory_handler(p, "bigfloat_matrix", "gauss_jordan(b) :: not enough memory");
	for (i = 0; i < rows; i++) {
		p[i] = new bigfloat_int[columns];
		memory_handler(p[i], "bigfloat_matrix", "gauss_jordan(b) :: not enough memory");
	}

	for (i = 0; i < rows; i++)
		for (j = 0; j < columns; j++)
			p[i][j].assign(value[i][j]);

	// Begin

	for (i = 0; i < rows; i++) {
		// search for the biggest absolute value in column

		max.assign_zero();
		k = i;
		for (j = i; j < rows; j++)
			for (k = i; k < rows; k++) {
				p[j][k].approx(temp);
				temp = abs(temp);
				if (max < temp) {
					r = j;
					c = k;
					max.assign(temp);
				}
			}
		if (r != i) {
			helpr = p[i];
			p[i] = p[r];
			p[r] = helpr;
			LiDIA::swap(l[i], l[r]);
		}
		if (c != i) {
			for (j = 0; j < rows; j++) {
				factor.assign(p[j][c]);
				p[j][c].assign(p[j][i]);
				p[j][i].assign(factor);
			}
			k = index[i];
			index[i] = index[c];
			index[c] = k;
		}

		// is it zero ?  then finish

		p[i][i].approx(temp);
		temp = abs(temp);
		if (temp.is_approx_zero())
			break;

		// else continue

		helpr = p[i];
		for (j = i+1; j < rows; j++) {
			LiDIA::divide(factor, p[j][i], p[i][i]);
			help0 = p[j];
			for (k = i; k < columns; k++) {
				LiDIA::multiply(help1[k], helpr[k], factor);
				LiDIA::subtract(help0[k], help0[k], help1[k]);
			}
			LiDIA::multiply(tmp, l[i], factor);
			LiDIA::subtract(l[j], l[j], tmp);
		}
		for (j = 0; j < i; j++) {
			LiDIA::divide(factor, p[j][i], p[i][i]);
			help0 = p[j];
			for (k = i; k < columns; k++) {
				LiDIA::multiply(help1[k], helpr[k], factor);
				LiDIA::subtract(help0[k], help0[k], help1[k]);
			}
			LiDIA::multiply(tmp, l[i], factor);
			LiDIA::subtract(l[j], l[j], tmp);
		}
	}

	// check if the equation-system is solvable

	for (j = i; j < rows; j++) {
		l[j].approx(temp);
		if (!temp.is_approx_zero()) {
			for (i = 0; i < rows; i++)
				x[i].assign_zero();
			delete[] index; // MM
			return; // no solution, return zero-vector
		}
	}

	// else compute x in Ax=b

	for (j = 0; j < i; j++) {
		LiDIA::divide(factor, l[j], p[j][j]);
		factor.approx(x[index[j]]);
	}

	for (j = i; j < columns; j++)
		x[index[j]].assign_zero();

	for (i = 0; i < rows; i++)
		delete[] p[i];
	delete[] p;

	delete[] help1;
	delete[] l;
	delete[] index; // MM
}



base_vector< bigfloat >
bigfloat_matrix::gauss_jordan(const base_vector< bigfloat > & b) const
{
	debug_handler("bigfloat_matrix", "gauss_jordan_pc(b)");
	if (rows != b.size())
		lidia_error_handler("bigfloat_lattice", "gauss_jordan_pc(b) :: illegal vector size");
	if (rows > columns)
		lidia_error_handler("bigfloat_lattice", "gauss_jordan_pc(b) :: "
				    "rows have to be less or equal columns");

	base_vector< bigfloat > x(columns, vector_flags(vector_flags::fixed));
	long i, j, k, r = 0, c = 0, *index = new long[columns];
	bigfloat_int factor, tmp, *help0, *help1, *helpr, *l, **p;
	bigfloat max, temp;

	help1 = new bigfloat_int[columns];
	memory_handler(help1, "bigfloat_lattice", "gauss_jordan_pc(b) :: not enough memory !");

	l = new bigfloat_int[rows];
	memory_handler(l, "bigfloat_lattice", "gauss_jordan_pc(b) :: not enough memory !");

	for (i = 0; i < rows; i++)
		l[i].assign(b[i]);
	for (i = 0; i < columns; i++)
		index[i] = i;

	p = new bigfloat_int*[rows];
	memory_handler(p, "bigfloat_matrix", "gauss_jordan(b) :: not enough memory");
	for (i = 0; i < rows; i++) {
		p[i] = new bigfloat_int[columns];
		memory_handler(p[i], "bigfloat_matrix", "gauss_jordan(b) :: not enough memory");
	}

	for (i = 0; i < rows; i++)
		for (j = 0; j < columns; j++)
			p[i][j].assign(value[i][j]);

	// Begin

	for (i = 0; i < rows; i++) {
		// search for the biggest absolute value in column

		max.assign_zero();
		k = i;
		for (j = i; j < rows; j++)
			for (k = i; k < rows; k++) {
				p[j][k].approx(temp);
				temp = abs(temp);
				if (max < temp) {
					r = j;
					c = k;
					max.assign(temp);
				}
			}
		if (r != i) {
			helpr = p[i];
			p[i] = p[r];
			p[r] = helpr;
			LiDIA::swap(l[i], l[r]);
		}
		if (c != i) {
			for (j = 0; j < rows; j++) {
				factor.assign(p[j][c]);
				p[j][c].assign(p[j][i]);
				p[j][i].assign(factor);
			}
			k = index[i];
			index[i] = index[c];
			index[c] = k;
		}

		// is it zero ?  then finish

		p[i][i].approx(temp);
		temp = abs(temp);
		if (temp.is_approx_zero())
			break;

		// else continue

		helpr = p[i];
		for (j = i+1; j < rows; j++) {
			LiDIA::divide(factor, p[j][i], p[i][i]);
			help0 = p[j];
			for (k = i; k < columns; k++) {
				LiDIA::multiply(help1[k], helpr[k], factor);
				LiDIA::subtract(help0[k], help0[k], help1[k]);
			}
			LiDIA::multiply(tmp, l[i], factor);
			LiDIA::subtract(l[j], l[j], tmp);
		}
		for (j = 0; j < i; j++) {
			LiDIA::divide(factor, p[j][i], p[i][i]);
			help0 = p[j];
			for (k = i; k < columns; k++) {
				LiDIA::multiply(help1[k], helpr[k], factor);
				LiDIA::subtract(help0[k], help0[k], help1[k]);
			}
			LiDIA::multiply(tmp, l[i], factor);
			LiDIA::subtract(l[j], l[j], tmp);
		}
	}

	// check if the equation-system is solvable

	for (j = i; j < rows; j++) {
		l[j].approx(temp);
		if (!temp.is_approx_zero()) {
			delete [] index; // MM
			return x; // no solution, return zero-vector
		}
	}

	// else compute x in Ax=b

	for (j = 0; j < i; j++) {
		LiDIA::divide(factor, l[j], p[j][j]);
		factor.approx(x[index[j]]);
	}

	for (j = i; j < columns; j++)
		x[index[j]].assign_zero();

	for (i = 0; i < rows; i++)
		delete[] p[i];
	delete[] p;

	delete[] help1;
	delete[] l;
	delete[] index; // MM
	// return it

	return x;
}



void
bigfloat_matrix::cholesky(bigfloat_matrix & A) const
{
	debug_handler("bigfloat_matrix", "cholesky()");
	if (rows != columns)
		lidia_error_handler("bigfloat_matrix", "cholesky() :: matrix is not symmetrical");

	bigfloat_int l, **p, *helpr;
	bigfloat max, temp;
	long n, k, i, r, c, j, *index = new long[columns];

	p = new bigfloat_int*[rows];
	memory_handler(p, "bigfloat_matrix", "cholesky() :: not enough memory");
	for (n = 0; n < rows; n++) {
		p[n] = new bigfloat_int[columns];
		memory_handler(p[n], "bigfloat_matrix", "cholesky() :: not enough memory");
		index[n] = n;
	}

	for (n = 0; n < rows; n++) {
		max.assign_zero();
		r = c = n;
		for (j = n; j < rows; j++)
			for (k = n; k < columns; k++) {
				p[j][k].approx(temp);
				temp = abs(temp);
				if (max < temp) {
					r = j;
					c = k;
					max.assign(temp);
				}
			}
		if (c != n) {
			helpr = p[n];
			p[n] = p[c];
			p[c] = helpr;
		}
		if (r != n) {
			for (i = 0; i < rows; i++) {
				l = p[i][r];
				p[i][r] = p[i][n];
				p[i][n] = l;
			}
			k = index[n];
			index[n] = index[r];
			index[r] = k;
		}


		(p[n][n]).assign(value[n][n]);
		for (k = 0; k < n; k++) {
			LiDIA::square(l, p[n][k]);
			LiDIA::subtract(p[n][n], p[n][n], l);
		}
		LiDIA::sqrt(p[n][n], p[n][n]);
		for (i = n+1; i < rows; i++) {
			(p[i][n]).assign(value[i][n]);
			for (k = 0; k < n; k++) {
				LiDIA::multiply(l, p[i][k], p[n][k]);
				LiDIA::subtract(p[i][n], p[i][n], l);
			}
			LiDIA::divide(p[i][n], p[i][n], p[n][n]);
		}
	}

	A.resize(rows, columns);
	for (n = 0; n < rows; n++)
		for (k = 0; k < rows; k++)
			p[n][index[k]].approx(A.value[n][k]);

	for (i = 0; i < rows; i++)
		delete[] p[i];
	delete[] p;
	delete[] index; // MM
}



void
bigfloat_matrix::mod_cholesky(bigfloat_matrix & Q) const
{
	debug_handler("bigfloat_matrix", "mod_cholesky(Q)");
	if (rows != columns)
		lidia_error_handler("bigfloat_lattice", "mod_cholesky(Q) :: rows != columns");

	long i, j, s, l;
	bigfloat_int temp, **p;

	Q.resize(rows, columns);

	p = new bigfloat_int*[rows];
	for (i = 0; i < rows; i++)
		p[i] = new bigfloat_int[columns];

	for (i = 0; i < rows; i++)
		for (j = i; j < columns; j++)
			p[i][j].assign(value[i][j]);

	for (i = 0; i < rows; i++) {
		for (j = i+1; j < rows; j++) {
			p[j][i].assign(p[i][j]);
			LiDIA::divide(p[i][j], p[i][j], p[i][i]);
		}

		for (s = i+1; s < rows; s++) {
			for (l = s; l < columns; l++) {
				LiDIA::multiply(temp, p[s][i], p[i][l]);
				LiDIA::subtract(p[s][l], p[s][l], temp);
			}
		}
	}

	for (i = 0; i < rows; i++)
		for (j = 0; j < columns; j++)
			(p[i][j]).approx(Q.value[i][j]);

	for (i = 0; i < rows; i++)
		delete[] p[i];
	delete[] p;
}



void
bigfloat_matrix::invert(bigfloat_matrix & A) const
{
	debug_handler("bigfloat_matrix", "invert(A)");
	if (rows != columns)
		lidia_error_handler("bigfloat_lattice", "invert(A) :: rows != columns");

	long i, j, k, c = 0, r = 0, *index = new long[columns];
	bigfloat temp, max;
	bigfloat_int tempi, **p, *helpr, *help0, factor;
	bigfloat_matrix I(rows, columns), M(rows, columns+columns);

	p = new bigfloat_int*[rows];
	for (i = 0; i < rows; i++)
		p[i] = new bigfloat_int[columns+columns];

	for (i = 0; i < rows; i++) {
		for (j = 0; j < columns; j++)
			p[i][j].assign(value[i][j]);
		(p[i][i+columns]).assign_one();
		index[i] = i;
	}

	// Begin

	for (i = 0; i < rows; i++) {
		// search for the biggest absolute value

		max.assign_zero();
		k = i;
		for (j = i; j < rows; j++)
			for (k = i; k < rows; k++) {
				p[j][k].approx(temp);
				temp = abs(temp);
				if (max < temp) {
					r = j;
					c = k;
					max.assign(temp);
				}
			}
		if (r != i) {
			helpr = p[i];
			p[i] = p[r];
			p[r] = helpr;
		}
		if (c != i) {
			for (k = 0; k < rows; k++) {
				LiDIA::swap(p[k][i], p[k][c]);
			}
			k = index[i];
			index[i] = index[c];
			index[c] = k;
		}

		// is it zero ?  then finish

		p[i][i].approx(temp);
		if (temp.is_approx_zero())
			break;

		// else continue

		p[i][i].Approx();
		helpr = p[i];
		for (j = i+1; j < rows; j++) {
			LiDIA::divide(factor, p[j][i], p[i][i]);
			help0 = p[j];
			for (k = i; k < columns+columns; k++) {
				LiDIA::multiply(tempi, helpr[k], factor);
				LiDIA::subtract(help0[k], help0[k], tempi);
			}
		}
		for (j = 0; j < i; j++) {
			LiDIA::divide(factor, p[j][i], p[i][i]);
			help0 = p[j];
			for (k = i; k < columns+columns; k++) {
				LiDIA::multiply(tempi, helpr[k], factor);
				LiDIA::subtract(help0[k], help0[k], tempi);
			}
		}
	}

	for (j = 0; j < rows; j++)
		for (i = 0; i < columns; i++)
			LiDIA::divide(p[j][i+columns], p[j][i+columns], p[j][j]);

	for (j = 0; j < rows; j++) {
		while (index[j] != j) {
			helpr = p[j];
			p[j] = p[index[j]];
			p[index[j]] = helpr;
			i = index[index[j]];
			index[index[j]] = index[j];
			index[j] = i;
		}
	}

	A.resize(rows, columns);
	for (i = 0; i < rows; i++)
		for (j = 0; j < columns; j++)
			p[i][columns+j].approx(A.value[i][j]);

	for (i = 0; i < rows; i++)
		delete[] p[i];
	delete[] p;
	delete[] index; // MM
}


void
bigfloat_matrix::LR(bigfloat_matrix & L, bigfloat_matrix & R) const
{
	debug_handler("bigfloat_matrix", "LR(L, R)");
	if (rows != columns)
		lidia_error_handler("bigfloat_matrix", "LR(L, R) :: rows != columns");
	lidia_size_t i, j, k;
	bigfloat_int t, **l, **r, **a;
	bigfloat temp;

	L.resize(rows, columns);
	R.resize(rows, columns);

	l = new bigfloat_int*[rows];
	for (i = 0; i < rows; i++)
		l[i] = new bigfloat_int[columns];

	r = new bigfloat_int*[rows];
	for (i = 0; i < rows; i++)
		r[i] = new bigfloat_int[columns];

	a = new bigfloat_int*[rows];
	for (i = 0; i < rows; i++)
		a[i] = new bigfloat_int[columns];

	for (i = 0; i < rows; i++)
		for (j = 0; j < columns; j++) {
			l[i][j].assign_zero();
			r[i][j].assign_zero();
			a[i][j].assign(value[i][j]);
		}

	for (j = 0; j < columns; j++) {
		for (k = 0; k < j; k++)
			for (i = k+1; i <= j; i++) {
				LiDIA::multiply(t, a[i][k], a[k][j]);
				LiDIA::subtract(a[i][j], a[i][j], t);
			}

		for (k = 0; k < j; k++)
			for (i = j+1; i < rows; i++) {
				LiDIA::multiply(t, a[i][k], a[k][j]);
				LiDIA::subtract(a[i][j], a[i][j], t);
			}

		for (i = j+1; i < rows; i++)
			LiDIA::divide(a[i][j], a[i][j], a[j][j]);
	}

	for (i = 0; i < rows; i++) {
		for (j = 0; j < i; j++)
			(a[i][j]).approx(L.value[i][j]);
		for (; j < columns; j++)
			(a[i][j]).approx(R.value[i][j]);
		(L.value[i][i]).assign_one();
	}

	for (i = 0; i < rows; i++) {
		delete[] l[i];
		delete[] r[i];
		delete[] a[i];
	}
	delete[] l;
	delete[] r;
	delete[] a;
}


void
bigfloat_matrix::QR(bigfloat_matrix & Q, bigfloat_matrix & R) const
{
	debug_handler("bigfloat_matrix", "QR(Q, R)");
	if (rows != columns)
		lidia_error_handler("bigfloat_matrix", "QR(Q, R) :: rows != columns");
	lidia_size_t i, j, k;
	bigfloat_int a, a2, c, **q, **r;
	bigfloat temp;
	base_vector< bigfloat_int > u(rows, vector_flags(vector_flags::fixed));
	base_vector< bigfloat_int > v(rows, vector_flags(vector_flags::fixed));
	base_vector< bigfloat_int > w(columns, vector_flags(vector_flags::fixed));
	base_vector< bigfloat_int > x(columns, vector_flags(vector_flags::fixed));

	Q.resize(rows, columns);
	R.resize(rows, columns);

	q = new bigfloat_int*[rows];
	for (i = 0; i < rows; i++)
		q[i] = new bigfloat_int[columns];

	for (i = 0; i < rows; i++) {
		for (j = 0; j < columns; j++)
			q[i][j].assign_zero();
		q[i][i].assign_one();
	}

	r = new bigfloat_int*[rows];
	for (i = 0; i < rows; i++)
		r[i] = new bigfloat_int[columns];

	for (i = 0; i < rows; i++)
		for (j = 0; j < columns; j++)
			r[i][j].assign(value[i][j]);

	for (j = 0; j < columns-1; j++) {
		a2.assign_zero();

		for (k = j; k < rows; k++) {
			square(a, r[k][j]);
			LiDIA::add(a2, a2, a);
		}

		a2.approx(temp);
		if (temp.is_approx_zero())
			lidia_error_handler("bigfloat_matrix", "QR(Q, R) :: Matrix was singular");
		a = sqrt(a2);

		LiDIA::multiply(c, a, abs(r[j][j]));
		LiDIA::add(c, c, a2);
		c.invert();

		if (r[j][j].is_lt_zero())
			a.negate();

		for (k = j; k < rows; k++)
			v[k].assign(r[k][j]);

		LiDIA::add(v[j], v[j], a);

		for (k = j; k < rows; k++)
			LiDIA::multiply(u[k], v[k], c);

		for (k = j; k < columns; k++) {
			w[k].assign_zero();
			x[k].assign_zero();
			for (i = j; i < rows; i++) {
				LiDIA::multiply(c, u[i], r[i][k]);
				LiDIA::add(w[k], w[k], c);
				LiDIA::multiply(c, u[i], q[i][k]);
				LiDIA::add(x[k], x[k], c);
			}
		}

		for (k = j; k < columns; k++)
			for (i = j; i < rows; i++) {
				LiDIA::multiply(c, v[i], w[k]);
				LiDIA::subtract(r[i][k], r[i][k], c);
				LiDIA::multiply(c, v[i], x[k]);
				LiDIA::subtract(q[i][k], q[i][k], c);
			}

		for (k = 0; k < j; k++) {
			x[k].assign_zero();
			for (i = j; i < rows; i++) {
				LiDIA::multiply(c, u[i], q[i][k]);
				LiDIA::add(x[k], x[k], c);
			}
		}

		for (k = 0; k < j; k++)
			for (i = j; i < rows; i++) {
				LiDIA::multiply(c, v[i], x[k]);
				LiDIA::subtract(q[i][k], q[i][k], c);
			}
	}

	for (i = 1; i < rows; i++)
		for (j = 0; j < i; j++)
			r[i][j].assign_zero();

	for (i = 0; i < rows; i++)
		for (j = 0; j < columns; j++)
			(r[i][j]).approx(R.value[i][j]);

	for (i = 0; i < rows; i++)
		for (j = 0; j < columns; j++)
			(q[i][j]).approx(Q.value[j][i]);

	for (i = 0; i < rows; i++) {
		delete[] q[i];
		delete[] r[i];
	}
	delete[] q;
	delete[] r;
}

int
bigfloat_matrix::equal(const bigfloat_matrix & A) const
{
	debug_handler("bigfloat_matrix", "equal(A)");
	if (rows != A.rows || columns != A.columns) {
		return 0;
	}
	lidia_size_t i = 0, r, c;
	bigfloat erg;
	for (r = 0; r < rows; r++)
		for (c = 0; c < columns; c++) {
			LiDIA::subtract(erg, value[r][c], A.value[r][c]);
			erg = abs(erg);
			if (erg.is_approx_zero())
				i++;
		}
	return (i == rows*columns);
}


void
bigfloat_matrix::randomize()
{
	debug_handler("bigfloat_matrix", "randomize()");
	for (int l = 0; l < rows; l++)
		for (int k = 0; k < columns; k++) {
			value[k][l].assign(bigfloat(LiDIA::randomize(bigint(10000000))));
		}
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
