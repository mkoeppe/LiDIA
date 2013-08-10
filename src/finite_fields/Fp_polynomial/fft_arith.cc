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
//	Author	: Victor Shoup (VS), Thomas Pfahler (TPf)
//	Changes	: See CVS log
//
//==============================================================================================


#ifdef HAVE_CONFIG_H
# include	"config.h"
#endif
#include	"LiDIA/Fp_polynomial.h"
#include	"LiDIA/finite_fields/Fp_polynomial_fft.h"
#include	"LiDIA/crt.h"



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



//***************************************************************
//	This File contains the implementation of the functions
//	-	fft_mul, fft_sqr
//	-	fft_rem, fft_div, fft_div_rem (and copy_reverse)
//	-	newton_inv (and class poly_mod_rep)
//	- 	build_from_roots (and two auxiliary functions)
//***************************************************************



//***************************************************************
//			fft_mul, fft_sqr
//***************************************************************

void fft_mul(Fp_polynomial& x, const Fp_polynomial& a, const Fp_polynomial& b)
{
	debug_handler("fftmul.c", "fft_mul(Fp_polynomial&, Fp_polynomial&, Fp_polynomial&)");

	lidia_size_t deg_a = a.degree();
	lidia_size_t deg_b = b.degree();

	if (deg_a < 0 || deg_b < 0) {
		//zero polynomial has degree -1
		x.set_modulus(a); //assigns zero
		return;
	}

	lidia_size_t d = deg_a + deg_b;
	lidia_size_t k = next_power_of_two(d+1);

	fft_data F(a.modulus(), k);

	lidia_size_t num_primes = F.number_of_primes();
	if (num_primes == 1) {
		fft_rep Ra(F), Rb(F);
		Ra.to_fft_rep(a);
		Rb.to_fft_rep(b);
		multiply(Ra, Ra, Rb);
		Ra.from_fft_rep(x, 0, d);
	}
	else {
		modular_fft_rep R1(F), R2(F); // more space efficient version
		lidia_size_t index;
		for (index = 0; index < num_primes; index++) {
			R1.to_modular_fft_rep(a, index);
			R2.to_modular_fft_rep(b, index);
			multiply(R1, R1, R2, index);
			R1.from_modular_fft_rep(0, d, index);
		}
		R1.get_result(x, 0, d);
	}
}



void fft_sqr(Fp_polynomial& x, const Fp_polynomial& a)
{
	debug_handler("fftmul.c", "fft_sqr(Fp_polynomial&, Fp_polynomial&)");

	lidia_size_t deg_a = a.degree();
	if (deg_a < 0) {
		x.set_modulus(a); //assigns zero
		return;
	}

	lidia_size_t d = 2*deg_a;
	lidia_size_t k = next_power_of_two(d+1);

	fft_data F(a.modulus(), k);

	lidia_size_t num_primes = F.number_of_primes();
	if (num_primes == 1) {
		fft_rep Ra(F);
		Ra.to_fft_rep(a);
		multiply(Ra, Ra, Ra);
		Ra.from_fft_rep(x, 0, d);
	}
	else {
		modular_fft_rep R1(F); // more space efficient version
		lidia_size_t index;
		for (index = 0; index < num_primes; index++) {
			R1.to_modular_fft_rep(a, index);
			multiply(R1, R1, R1, index);
			R1.from_modular_fft_rep(0, d, index);
		}
		R1.get_result(x, 0, d);
	}
}



//***************************************************************
//			fft_rem, -div, -div_rem
//***************************************************************

void fft_rem(Fp_polynomial& r, const Fp_polynomial& a, const Fp_polynomial& b)
{
	debug_handler("fftrem.c", "fft_rem (Fp_polynomial&, Fp_polynomial&, Fp_polynomial&)");

	lidia_size_t deg_b = b.degree(), deg_a = a.degree();
	if (deg_a < deg_b) {
		r.assign(a);
		return;
	}

	lidia_size_t m = deg_a - deg_b + 1;
	Fp_polynomial P1, P2, P3;
	copy_reverse(P3, b, 0, deg_b);
	invert(P2, P3, m);
	copy_reverse(P1, P2, 0, m-1);

	lidia_size_t k = next_power_of_two(2*m-1);
	lidia_size_t l = next_power_of_two(deg_b);
	lidia_size_t index;
	lidia_size_t mx = comparator< lidia_size_t >::max(k, l);

	fft_data F(a.modulus(), mx);

	lidia_size_t num_primes = F.number_of_primes();
	if (num_primes == 1) {
		fft_rep Ra(F), Rb(F);
		Ra.set_size(k);
		Rb.set_size(k);
		Ra.to_fft_rep(P1);
		Rb.to_fft_rep(a, deg_b, deg_a);
		multiply(Ra, Ra, Rb);
		Ra.from_fft_rep(P3, deg_a-deg_b, 2*(deg_a-deg_b));

		Ra.set_size(l);
		Rb.set_size(l);
		Ra.to_fft_rep(b, 0, b.degree());
		Rb.to_fft_rep(P3);
		multiply(Ra, Ra, Rb);
		Ra.from_fft_rep(P3, 0, deg_b-1);
	}
	else {
		modular_fft_rep R1(F); // more space efficient version
		modular_fft_rep R2(F);

		R1.set_size(k);
		R2.set_size(k);
		for (index = 0; index < num_primes; index++) {
			R1.to_modular_fft_rep(P1, index);
			R2.to_modular_fft_rep(a, deg_b, deg_a, index);
			multiply(R1, R1, R2, index);
			R1.from_modular_fft_rep(deg_a-deg_b, 2*(deg_a-deg_b), index);
		}
		R1.get_result(P3, deg_a-deg_b, 2*(deg_a-deg_b));

		R1.set_size(l);
		R2.set_size(l);
		for (index = 0; index < num_primes; index++) {
			R1.to_modular_fft_rep(b, 0, b.degree(), index);
			R2.to_modular_fft_rep(P3, index);
			multiply(R1, R1, R2, index);
			R1.from_modular_fft_rep(0, deg_b-1, index);
		}
		R1.get_result(P3, 0, deg_b-1);
	}

	lidia_size_t L = 1 << l;
	cyclic_reduce(P2, a, L);
	trunc(r, P2, deg_b);
	subtract(r, r, P3);
}



void fft_div(Fp_polynomial& q, const Fp_polynomial& a, const Fp_polynomial& b)
{
	debug_handler("fftrem.c", "fft_div (Fp_polynomial&, Fp_polynomial&, Fp_polynomial&)");

	lidia_size_t deg_b = b.degree(), deg_a = a.degree();
	if (deg_a < deg_b) {
		q.set_modulus(b); //assigns zero
		return;
	}

	lidia_size_t m = deg_a - deg_b + 1;
	Fp_polynomial P1, P2, P3;
	copy_reverse(P3, b, 0, deg_b);
	invert(P2, P3, m);
	copy_reverse(P1, P2, 0, m-1);

	lidia_size_t k = next_power_of_two(2*m-1);

	fft_data F(a.modulus(), k);

	lidia_size_t num_primes = F.number_of_primes();
	if (num_primes == 1) {
		fft_rep Ra(F), Rb(F);
		Ra.set_size(k);
		Rb.set_size(k);
		Ra.to_fft_rep(P1);
		Rb.to_fft_rep(a, deg_b, deg_a);
		multiply(Ra, Ra, Rb);
		Ra.from_fft_rep(q, deg_a-deg_b, 2*(deg_a-deg_b));
	}
	else {
		modular_fft_rep R1(F); // more space efficient version
		modular_fft_rep R2(F);

		lidia_size_t index;
		for (index = 0; index < num_primes; index++) {
			R1.to_modular_fft_rep(P1, index);
			R2.to_modular_fft_rep(a, deg_b, deg_a, index);

			multiply(R1, R1, R2, index);
			R1.from_modular_fft_rep(deg_a-deg_b, 2*(deg_a-deg_b), index);
		}
		R1.get_result(q, deg_a-deg_b, 2*(deg_a-deg_b));
	}
}



void fft_div_rem(Fp_polynomial& q, Fp_polynomial& r, const Fp_polynomial& a, const Fp_polynomial& b)
{
	debug_handler("fftrem.c", "fft_div_rem (Fp_polynomial&, Fp_polynomial&, Fp_polynomial&, Fp_polynomial&)");

	lidia_size_t deg_b = b.degree(), deg_a = a.degree();
	if (deg_a < deg_b) {
		q.set_modulus(b); //assigns zero
		r.assign(a);
		return;
	}

	Fp_polynomial P1, P2, P3;
	copy_reverse(P3, b, 0, deg_b);
	invert(P2, P3, deg_a-deg_b+1);
	copy_reverse(P1, P2, 0, deg_a-deg_b);

	lidia_size_t k = next_power_of_two(2*(deg_a-deg_b)+1);
	lidia_size_t l = next_power_of_two(deg_b);
	lidia_size_t index;
	lidia_size_t mx = comparator< lidia_size_t >::max(k, l);

	fft_data F(a.modulus(), mx);

	lidia_size_t num_primes = F.number_of_primes();
	if (num_primes == 1) {
		fft_rep Ra(F), Rb(F);
		Ra.set_size(k);
		Rb.set_size(k);
		Ra.to_fft_rep(P1);
		Rb.to_fft_rep(a, deg_b, deg_a);
		multiply(Ra, Ra, Rb);
		Ra.from_fft_rep(P3, deg_a-deg_b, 2*(deg_a-deg_b));

		Ra.set_size(l);
		Rb.set_size(l);
		Ra.to_fft_rep(b, 0, b.degree());
		Rb.to_fft_rep(P3);
		multiply(Ra, Ra, Rb);
		Ra.from_fft_rep(P1, 0, deg_b-1);
	}
	else {
		modular_fft_rep R1(F); // more space efficient version
		modular_fft_rep R2(F);

		R1.set_size(k);
		R2.set_size(k);
		for (index = 0; index < num_primes; index++) {
			R1.to_modular_fft_rep(P1, index);
			R2.to_modular_fft_rep(a, deg_b, deg_a, index);
			multiply(R1, R1, R2, index);
			R1.from_modular_fft_rep(deg_a-deg_b, 2*(deg_a-deg_b), index);
		}
		R1.get_result(P3, deg_a-deg_b, 2*(deg_a-deg_b));

		R1.set_size(l);
		R2.set_size(l);
		for (index = 0; index < num_primes; index++) {
			R1.to_modular_fft_rep(b, index);
			R2.to_modular_fft_rep(P3, index);
			multiply(R1, R1, R2, index);
			R1.from_modular_fft_rep(0, deg_b-1, index);
		}
		R1.get_result(P1, 0, deg_b-1);
	}

	lidia_size_t L = 1 << l;
	cyclic_reduce(P2, a, L);
	trunc(r, P2, deg_b);
	subtract(r, r, P1);

	q.assign(P3);
}



// x[0..hi-lo+1] = reverse(a[lo..hi]), with zero fill
// input may not alias output
void copy_reverse(Fp_polynomial& x, const Fp_polynomial& a, lidia_size_t lo, lidia_size_t hi)
{
	debug_handler("fftrem.c", "copy_reverse(Fp_polynomial&, Fp_polynomial&, lidia_size_t, lidia_size_t) ");
	lidia_size_t i, j, n, m;

	n = hi-lo+1;
	m = a.degree()+1; // = a.c_length

	x.set_modulus(a);
	x.set_degree(n-1);

	const bigint* ap = a.coeff;
	bigint* xp = x.coeff;

	for (i = 0; i < n; i++) {
		j = hi-i;
		if (j< 0 || j >= m)
			xp[i].assign_zero();
		else
			xp[i] = ap[j];
	}

	x.remove_leading_zeros();
}



//***************************************************************
//		    newton_inv, class poly_mod_rep
//***************************************************************


class poly_mod_rep
{

	// This data structure holds unconvoluted modular representations
	// of polynomials
	// used only in function newton_inv

	lidia_size_t num_primes;
	lidia_size_t size;
	fft_prime_t **tbl;
	fft_data fd;

	poly_mod_rep();
	poly_mod_rep(const poly_mod_rep&); // disable

public:
	poly_mod_rep(const Fp_polynomial&, lidia_size_t, lidia_size_t,
		     const fft_data&);

	~poly_mod_rep();

	friend class modular_fft_rep;
};


poly_mod_rep::poly_mod_rep(const Fp_polynomial& a, lidia_size_t lo,
			   lidia_size_t hi, const fft_data &F) :
	fd(F)
{
	debug_handler("poly_mod_rep", "poly_mod_rep (Fp_polynomial&, lidia_size_t, lidia_size_t, fft_data&)");

	if (lo < 0 || hi < 0 || F.node == 0) {
		lidia_error_handler("poly_mod_rep", "poly_mod_rep(...)::bad args");
		return;
	}

	num_primes = F.number_of_primes();

	hi = comparator< lidia_size_t >::min(hi, a.degree());
	size = comparator< lidia_size_t >::max(hi-lo+1, 0);
	lidia_size_t i;

	tbl = new fft_prime_t*[num_primes];
	memory_handler(tbl, "poly_mod_rep", "poly_mod_rep(...)::"
		       "Error in memory allocation (tbl)");
	for (i = 0; i < num_primes; i++) {
		tbl[i] = new fft_prime_t[size];
		memory_handler(tbl[i], "poly_mod_rep", "poly_mod_rep(...)::"
			       "Error in memory allocation (tbl[i])");
	}

#if 0
	crt help(*F.crt_node->ct);
	for (i = 0; i < num_primes; i++)
		help.reduce(tbl[i], &a.coeff[lo], size, i);
#else
	for (i = 0; i < num_primes; i++) {
		fft_prime_t q = fd.crt_node->ct->get_prime(i);
		fft_prime_t *yp = tbl[i];
		const bigint *ap = &a.coeff[lo];

		for (lidia_size_t j = 0; j < size; j++)
			remainder(yp[j], ap[j], q);
	}
#endif
}



poly_mod_rep::~poly_mod_rep()
{
	debug_handler("poly_mod_rep", "~poly_mod_rep ()");

	if (tbl)	//should never be zero
		for (lidia_size_t i = 0; i < num_primes; i++)
			delete[] tbl[i];
	delete[] tbl;
}



void modular_fft_rep::to_modular_fft_rep(const poly_mod_rep &a,
					 lidia_size_t lo, lidia_size_t hi, lidia_size_t index)
	// converts coefficients lo..hi to a 2^k-point fft_rep.
	// must have hi-lo+1 < 2^k
{
	debug_handler("modular_fft_rep", "to_modular_fft_rep(poly_mod_rep&, lidia_size_t, lidia_size_t, lidia_size_t)");

	if (fd != a.fd || k < 0 || lo < 0) {
		lidia_error_handler("modular_fft_rep", "to_modular_fft_rep"
				    "(poly_mod_rep&, lidia_size_t, lidia_size_t, lidia_size_t)::"
				    "bad args");
		return;
	}

	if (hi > a.size - 1)
		hi = a.size - 1;

	lidia_size_t K = 1 << k;
	lidia_size_t j, m = comparator< lidia_size_t >::max(hi-lo + 1, 0);

	if (m > K) {
		lidia_error_handler("modular_fft_rep", "to_modular_fft_rep"
				    "(poly_mod_rep&, lidia_size_t, lidia_size_t, lidia_size_t)::"
				    "hi-lo+1 is too large");
		return;
	}

	fft_prime_t *ap = (m == 0 ? 0 : a.tbl[index]);

	for (j = 0; j < m; j++)
		stat_vec[j] = ap[lo+j];
	for (; j < K; j++)
		stat_vec[j] = 0;

	bool ok;
	ok = fd.node->prime[index].evaluate(vec, 0, stat_vec, K-1, k);
	if (!ok)
		lidia_error_handler("modular_fft_rep", "to_modular_fft_rep"
				    "(poly_mod_rep&, ...)::re-init of FFT prime");
}



void newton_inv(Fp_polynomial& x, const Fp_polynomial& a, lidia_size_t m)
{
	debug_handler("fftrem.c", "newton_inv(Fp_polynomial&, Fp_polynomial&, lidia_size_t)");

	x.set_modulus(a);
	x.set_degree(m-1);

	lidia_size_t index, t, k;
	const bigint & p = a.modulus();
	lidia_size_t crov = Fp_polynomial::crossovers.log2_newton_crossover(p);

	plain_inv(x, a, (1 << crov));
	t = next_power_of_two(m);

	fft_data F(p, t);
	fft_rep R1(F);
	modular_fft_rep R2(F);
	Fp_polynomial P1;
	P1.set_max_degree(m/2 - 1);

	lidia_size_t a_len = comparator< lidia_size_t >::min(m, a.c_length);

	poly_mod_rep a_rep(a, 0, a_len-1, F);

	t = crov;
	k = 1 << t;

	lidia_size_t num_primes = F.number_of_primes();
	while (k < m) {
		lidia_size_t l = comparator< lidia_size_t >::min(2*k, m);

		R1.set_size(t+1);
		R2.set_size(t+1);

		R1.to_fft_rep(x);
		for (index = 0; index < num_primes; index++) {
			R2.to_modular_fft_rep(a_rep, 0, l-1, index);
			multiply(R2, R1, R2, index);
			R2.from_modular_fft_rep(k, l-1, index);
		}
		R2.get_result(P1, k, l-1);

		R2.set_size(t+1);
		for (index = 0; index < num_primes; index++) {
			R2.to_modular_fft_rep(P1, index);
			multiply(R2, R1, R2, index);
			R2.from_modular_fft_rep(0, l-k-1, index);
		}
		R2.get_result(P1, 0, l-k-1);

		x.set_degree(l-1);

		lidia_size_t i, y_len = P1.c_length;
		for (i = k; i < l; i++) {
			if (i-k >= y_len)
				(x.coeff[i]).assign_zero();
			else
				NegateMod(x.coeff[i], P1.coeff[i-k], p);
		}
		x.remove_leading_zeros();

		t++;
		k = l;
	}
}



//***************************************************************
//				build_from_roots
//***************************************************************

static
void iter_build(bigint* a, lidia_size_t n, const bigint& p)
{
	debug_handler("Fp_polynomial", "iter_build(bigint*, lidia_size_t, const bigint&)");
	lidia_size_t i, k;
	bigint b, t;

	if (n <= 0) return;

	NegateMod(a[0], a[0], p);

	for (k = 1; k < n; k++) {
		NegateMod(b, a[k], p);
		AddMod(a[k], b, a[k-1], p);
		for (i = k-1; i > 0; i--) {
			MulMod(t, a[i], b, p);
			AddMod(a[i], t, a[i-1], p);
		}
		MulMod(a[0], a[0], b, p);
	}
}



static
void mul_build(bigint* x, const bigint* a, const bigint* b, lidia_size_t n, const bigint& p)
{
	debug_handler("Fp_polynomial", "mul_build(bigint*, const bigint*, const bigint*, lidia_size_t, const bigint&)");
	bigint t, accum; //static
	lidia_size_t i, j, jmin, jmax;

	lidia_size_t d = 2*n-1;

	for (i = 0; i <= d; i++) {
		jmin = comparator< lidia_size_t >::max(0, i-(n-1));
		jmax = comparator< lidia_size_t >::min(n-1, i);
		accum.assign_zero();
		for (j = jmin; j <= jmax; j++) {
			multiply(t, a[j], b[i-j]);
			add(accum, accum, t);
		}
		if (i >= n) {
			add(accum, accum, a[i-n]);
			add(accum, accum, b[i-n]);
		}

		Remainder(x[i], accum, p);
	}
}



// computes the polynomial (X-a[0]) ... (X-a[n-1]), where n = a.length()
void Fp_polynomial::build_from_roots(const base_vector< bigint > & a)
{
	debug_handler("Fp_polynomial", "build_from_roots(base_vector< bigint > &)");

	const bigint& p = modulus();
	if (p.is_zero()) {
		lidia_error_handler("Fp_polynomial", "build_from_roots"
				    "(base_vector< bigint > &)::modulus was not set");
		return;
	}

	lidia_size_t n = a.size();
	if (n == 0) {
		this->assign_one();
		return;
	}

	lidia_size_t crov = Fp_polynomial::crossovers.fftmul_crossover(p);
	lidia_size_t k0 = next_power_of_two(crov);
	lidia_size_t crossover = 1 << k0;

	if (n <= crossover) {
		set_max_degree(n);
		assign(a, p);
		iter_build(coeff, n, p);
		//		set_degree(n);
		set_coefficient(n);
		return;
	}

	lidia_size_t k = next_power_of_two(n);
	lidia_size_t m = 1 << k;
	lidia_size_t i, j, index, l, width;

	Fp_polynomial b;

	b.assign(a, p);
	b.set_degree(m);

	for (i = n; i < m; i++)
		b[i].assign_zero();

	b[m].assign_one();

	bigint t1, one(1);

	bigint* g = new bigint[crossover];
	memory_handler(g, "Fp_polynomial", "build_from_roots(base_vector< bigint > &"
		       ")::Error in memory allocation");
	bigint* h = new bigint[crossover];
	memory_handler(h, "Fp_polynomial", "build_from_roots(base_vector< bigint > &"
		       ")::Error in memory allocation");
	bigint *tmp;


	for (i = 0; i < m; i += crossover) {
		for (j = 0; j < crossover; j++)
			LiDIA::NegateMod(g[j], b[i+j], p);

		if (k0 > 0) {
			for (j = 0; j < crossover; j += 2) {
				MulMod(t1, g[j], g[j+1], p);
				AddMod(g[j+1], g[j], g[j+1], p);
				g[j].assign(t1);
			}
		}

		for (l = 1; l < k0; l++) {
			width = 1 << l;
			for (j = 0; j < crossover; j += 2*width)
				mul_build(&h[j], &g[j], &g[j+width], width, p);
			tmp = g; g = h; h = tmp;
		}

		for (j = 0; j < crossover; j++)
			b[i+j].assign(g[j]);
	}

	fft_data F(p, k);
	modular_fft_rep R1(F), R2(F);

	for (l = k0; l < k; l++) {
		width = 1 << l;
		for (i = 0; i < m; i += 2*width) {
			R1.set_size(l+1);
			R2.set_size(l+1);
			for (index = 0; index < F.number_of_primes(); index++) {
				swap(one, b.coeff[i+width]);
				//t1 = b[i+width]; (b[i+width]).assign_one();

				R1.to_modular_fft_rep(b, i, i+width, index);

				swap(one, b.coeff[i+width]);
				swap(b.coeff[i+2*width], one);
				//b[i+width] = t1; t1 = b[i+2*width]; b[i+2*width].assign_one();

				R2.to_modular_fft_rep(b, i+width, i+2*width, index);

				swap(b.coeff[i+2*width], one);
				//b[i+2*width] = t1;

				multiply(R1, R1, R2, index);
				R1.from_modular_fft_rep(0, 2*width-1, index);
			}
			R1.get_result_ptr(&b.coeff[i], 0, 2*width-1);
			SubMod(b.coeff[i], b.coeff[i], one, p);
			//subtract(b[i], b[i], one);
		}
	}

	set_degree(n);
	lidia_size_t delta = m-n;
	for (i = 0; i <= n; i++)
		coeff[i].assign(b[i+delta]);

	delete [] h; // VM
	delete [] g; // VM

	// no need to normalize
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
