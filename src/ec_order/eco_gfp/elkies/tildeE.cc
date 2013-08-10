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
//	Author	:Frank Lehmann(FL), Volker Mueller (VM), Markus Maurer (MM)
//	Changes	: See CVS log
//
//==============================================================================================



#ifdef HAVE_CONFIG_H
# include	"config.h"
#endif
#include	"LiDIA/eco_prime.h"



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



int eco_prime::tildeEA (ff_element & tildeEa, ff_element & tildeEb,
                        ff_element & P1, const ff_element & A_tau,
                        const ff_element & j_ltau, const ff_element & Ea,
                        const ff_element & Eb)
{
	// powers of A(tau), j(tau), j(ltau)

	// base_vector< ff_element > Avec(l+2, 0), jvec(v+1, 0), jlvec(v+1, 0);
	// base_vector< ff_element > &A  = Avec;
	// base_vector< ff_element > &j  = jvec;
	// base_vector< ff_element > &jl = jlvec;

	ff_element * j, *jl, *A;

	j = new ff_element[v+1];
	jl = new ff_element[v+1];
	A = new ff_element[l+2];

	ff_element DF1, DF2;
	ff_element DJ1, DJ2;
	ff_element DFF1, DFF2;
	ff_element DJJ1, DJJ2;
	ff_element DFJ1, DFJ2;

	ff_element E4_tau, E6_tau;
	ff_element E4_ltau, E6_ltau;
	ff_element A_tau_1;
	ff_element j_tau, j_tau_1;
	ff_element j_ltau_1;
	ff_element E2_star;
	ff_element delta_tau, delta_ltau;
	ff_element W1, W2;

	ff_element tmp1, tmp2, tmp3;
	ff_element l_r, tmp;

	udigit r, k;
	int  rc;

	// determine E4(tau) und E6(tau)

	divide (E4_tau, Ea, ff_element(-3));
	divide (E6_tau, Eb, ff_element(-2));


	// determine j-invariant

	compute_jinv (j_tau, Ea, Eb);


	// j'(tau) = - j * E6 / E4

	divide   (j_tau_1, E6_tau, E4_tau);
	multiply (j_tau_1, j_tau_1, j_tau);
	negate   (j_tau_1, j_tau_1);


	// computation of DF1, DJ1, DFFi, DJJi, DFJi (i=1, 2)

	// compute powers of A(tau), j(tau), j(ltau)

	A[0]  = 1; A[1]  = A_tau;
	j[0]  = 1; j[1]  = j_tau;
	jl[0] = 1; jl[1] = j_ltau;

	for (r = 2; r <= l+1; r++) {
		multiply (A[r], A[r-1], A[1]);
	}
	for (k = 2; k <= v; k++) {
		multiply (j[k], j[k-1], j[1]);
	}
	for (k = 2; k <= v; k++) {
		multiply (jl[k], jl[k-1], jl[1]);
	}


	// initialize file for reading the cofficients

	rc = meq_prime::reset ();
	if (rc) {
		lidia_error_handler("eco_prime", "tildeEA()::Unable to reset file");
		delete [] j; delete[] jl; delete[] A;
		return 0;
	}

	DF1  = DJ1  =        0;
	DFF1 = DJJ1 = DFJ1 = 0;
	DFF2 = DJJ2 = DFJ2 = 0;


	for (r = 0; r <= l+1; r++) {
		rc = read_row (p);

		if (rc) {
			lidia_error_handler("eco_prime", "tildeEA()::Unable to reset file");
			delete [] j; delete[] jl; delete[] A;
			return 0;
		}

		if (r > 0) {
			tmp1 = 0;

			for (k = 0; k <= v; k++) {
				multiply (tmp2, row[k], j[k]);
				add      (tmp1, tmp1, tmp2);
			}
			multiply (tmp2, A[r-1], r);
			multiply (tmp1, tmp1, tmp2);
			add      (DF1, DF1, tmp1);
		}

		if (r > 0) {
			tmp1 = 0;

			for (k = 0; k <= v; k++) {
				multiply (tmp2, row[k], jl[k]);
				add      (tmp1, tmp1, tmp2);
			}
			multiply (tmp2, A[r-1], r);
			multiply (tmp1, tmp1, tmp2);
			add      (DF2, DF2, tmp1);
		}

		tmp1 = 0;

		for (k = 1; k <= v; k++) {
			multiply (tmp2, row[k], ff_element(k));
			multiply (tmp2, tmp2, j[k-1]);
			add      (tmp1, tmp1, tmp2);
		}
		multiply (tmp1, tmp1, A[r]);
		add      (DJ1, DJ1, tmp1);


		tmp1 = 0;

		for (k = 1; k <= v; k++) {
			multiply (tmp2, row[k], ff_element(k));
			multiply (tmp2, tmp2, jl[k-1]);
			add      (tmp1, tmp1, tmp2);
		}
		multiply (tmp1, tmp1, A[r]);
		add      (DJ2, DJ2, tmp1);


		if (r > 1) {
			tmp1 = 0;

			for (k = 0; k <= v; k++) {
				multiply (tmp2, row[k], j[k]);
				add      (tmp1, tmp1, tmp2);
			}
			multiply (tmp2, A[r-2], r*(r-1));
			multiply (tmp1, tmp1, tmp2);
			add      (DFF1, DFF1, tmp1);
		}

		if (r > 1) {
			tmp1 = 0;

			for (k = 0; k <= v; k++) {
				multiply (tmp2, row[k], jl[k]);
				add      (tmp1, tmp1, tmp2);
			}
			multiply (tmp2, A[r-2], r*(r-1));
			multiply (tmp1, tmp1, tmp2);
			add      (DFF2, DFF2, tmp1);
		}

		tmp1 = 0;

		for (k = 2; k <= v; k++) {
			multiply (tmp2, row[k], ff_element(k * (k-1)));
			multiply (tmp2, tmp2, j[k-2]);
			add      (tmp1, tmp1, tmp2);
		}
		multiply (tmp1, tmp1, A[r]);
		add      (DJJ1, DJJ1, tmp1);


		tmp1 = 0;

		for (k = 2; k <= v; k++) {
			multiply (tmp2, row[k], ff_element(k * (k-1)));
			multiply (tmp2, tmp2, jl[k-2]);
			add      (tmp1, tmp1, tmp2);
		}
		multiply (tmp1, tmp1, A[r]);
		add      (DJJ2, DJJ2, tmp1);


		if (r > 0) {
			tmp1 = 0;

			for (k = 1; k <= v; k++) {
				multiply (tmp2, row[k], ff_element(k));
				multiply (tmp2, tmp2, j[k-1]);
				add      (tmp1, tmp1, tmp2);
			}
			multiply (tmp2, A[r-1], r);
			multiply (tmp1, tmp1, tmp2);
			add      (DFJ1, DFJ1, tmp1);
		}

		if (r > 0) {
			tmp1 = 0;

			for (k = 1; k <= v; k++) {
				multiply (tmp2, row[k], ff_element(k));
				multiply (tmp2, tmp2, jl[k-1]);
				add      (tmp1, tmp1, tmp2);
			}
			multiply (tmp2, A[r-1], r);
			multiply (tmp1, tmp1, tmp2);
			add      (DFJ2, DFJ2, tmp1);
		}

	} // end for r

	// test whether a denominator is zero

	if (DF1.is_zero() || DF2.is_zero() || DJ1.is_zero() || DJ2.is_zero()) {
		lidia_error_handler("eco_prime", "tildeEA::division by zero");
		delete[] j; delete[] jl; delete[] A;
		return 1;
	}


	// A'(tau) = - j'(tau) * DJ1 / DF1

	multiply (A_tau_1, j_tau_1, DJ1);
	divide   (A_tau_1, A_tau_1, DF1);
	negate   (A_tau_1, A_tau_1);


	// j'(ltau) = - A'(tau) * DF2 / (l * DJ2)

	multiply (j_ltau_1, A_tau_1, DF2);
	multiply (tmp1, DJ2, l);
	divide   (j_ltau_1, j_ltau_1, tmp1);
	negate   (j_ltau_1, j_ltau_1);


	// W1 = -1/DF1 * [ A' DFF1 + 2 j' DFJ1 - j' DJJ1 DF1 / DJ1 ]

	divide   (tmp1, DF1, DJ1);
	multiply (tmp2, j_tau_1, DJJ1);
	multiply (W1, tmp1, tmp2);

	multiply (tmp1, j_tau_1, DFJ1);
	add      (tmp1, tmp1, tmp1);
	subtract (W1, tmp1, W1);

	multiply (tmp1, A_tau_1, DFF1);
	add      (W1, W1, tmp1);

	divide (W1, W1, DF1);
	negate (W1, W1);


	// W2 = -1/DF2 * [ A' DFF2 + 2 l j'(l) DFJ2 - l j'(l) DJJ2 DF2 / DJ2 ]

	divide   (tmp1, DF2, DJ2);
	multiply (tmp2, j_ltau_1, l);
	multiply (tmp2, tmp2, DJJ2);
	multiply (W2, tmp1, tmp2);

	multiply (tmp1, j_ltau_1, 2*l);
	multiply (tmp2, tmp1, DFJ2);
	subtract (W2, tmp2, W2);

	multiply (tmp1, A_tau_1, DFF2);
	add      (W2, W2, tmp1);

	divide  (W2, W2, DF2);
	negate  (W2, W2);


	// E4(ltau) = j_ltau_1 ^2 / (j_ltau * (j_ltau - 1728))

	if (j_ltau.is_zero()) {
		lidia_error_handler("eco_prime", "tildeEA::division by zero");
		delete[] j; delete[]jl; delete[] A;
		return 1;
	}

	multiply   (tmp, j_ltau_1, j_ltau_1);

	tmp2 = ff_element(1728);
	subtract (tmp3, j_ltau, tmp2);
	multiply (tmp3, tmp3, j_ltau);

	divide   (E4_ltau, tmp, tmp3);


	// E6(ltau) = - E4(ltau) * j_ltau_1 / j_ltau

	multiply (tmp, E4_ltau, j_ltau_1);
	divide   (tmp2, tmp, j_ltau);
	negate   (E6_ltau, tmp2);


	// delta(tau) = E4^3/j_tau

	multiply (tmp, E4_tau, E4_tau);
	multiply (tmp, tmp, E4_tau);
	divide   (delta_tau, tmp, j_tau);


	// delta(l*tau) = E4(l*tau)^3/j(l*tau)

	multiply (tmp, E4_ltau, E4_ltau);
	multiply (tmp, tmp, E4_ltau);
	divide   (delta_ltau, tmp, j_ltau);


	// test whether a denominator is equal to zero

	if (j_ltau_1.is_zero() || j_tau_1.is_zero() || delta_tau.is_zero() || delta_ltau.is_zero()) {
		lidia_error_handler("eco_prime", "tildeEA::division by zero");
		delete[] j; delete [] jl; delete[] A;
		return 1;
	}

	// E2_star = 6*(W2-W1)
        //         + 3*(l*j(lt)*E4(lt) / j'(lt) - j(t)*E4(t) / j'(t))
        //         + 4*(l*E4(lt)*E6(lt)^2 / (delta(lt) * j'(lt)) - E4(t)*E6(t)^2 / (delta(t) * j'(t)))

	multiply (tmp, j_ltau, l);
	divide   (tmp2, E4_ltau, j_ltau_1);
	multiply (tmp, tmp, tmp2);

	multiply (tmp2, j_tau, E4_tau);
	divide   (tmp3, tmp2, j_tau_1);

	subtract (tmp2, tmp, tmp3);
	multiply (tmp2, tmp2, 3);

	subtract (tmp, W2, W1);
	multiply (tmp, tmp, 6);

	add       (tmp, tmp, tmp2);

	multiply (tmp2, E4_ltau, 4 * l);
	multiply (tmp2, tmp2, E6_ltau);
	multiply (tmp2, tmp2, E6_ltau);

	multiply (tmp3, delta_ltau, j_ltau_1);

	divide   (tmp2, tmp2, tmp3);

	add      (tmp, tmp, tmp2);

	multiply (tmp2, E4_tau, 4);
	multiply (tmp2, tmp2, E6_tau);
	multiply (tmp2, tmp2, E6_tau);

	multiply (tmp3, delta_tau, j_tau_1);

	divide   (tmp2, tmp2, tmp3);

	subtract (E2_star, tmp, tmp2);


	// P1 = - l/2 * E2_star

	tmp  =  l;
	tmp2 =  2;

	divide   (tmp3, tmp, tmp2);
	multiply (tmp2, tmp3, E2_star);
	negate   (P1, tmp2);


	// E~ = (A~, B~),
	// A~ == -3 * l^4 * E4(l*tau),
	// B~ == -2 * l^6 * E6(l*tau)

	tmp = l;

	power  (l_r, tmp, 4);
	multiply (tildeEa, ff_element(-3), E4_ltau);
	multiply (tildeEa, tildeEa, l_r);

	power    (l_r, tmp, 6);
	multiply (tildeEb, ff_element(-2), E6_ltau);
	multiply (tildeEb, tildeEb, l_r);

	delete[] j; delete[] jl; delete[] A;

	return rc;
}



int eco_prime::tildeEf (ff_element & tildeEa, ff_element & tildeEb,
                        ff_element & P1, const ff_element & g_tau,
	                const ff_element & Ea, const ff_element & Eb)
{
  int rc;
  udigit r, k;
  ff_element * g, *j;

  g = new ff_element[l+2];
  j = new ff_element[v+1];
  
  ff_element E4_tau, E6_tau;
  ff_element E4_ltau, E6_ltau;
  ff_element f_tau;
  ff_element g_tau_1, f_tau_1;
  ff_element g_tau_1num, g_tau_1den;
  ff_element j_ltau_1num, j_ltau_1den;
  ff_element j_tau, j_ltau;
  ff_element j_tau_1, j_ltau_1;
  ff_element E2_star;
  ff_element E0_tau, E0_tau_1;
  ff_element delta_tau, delta_ltau;
  
  ff_element DF, DJ, DFF, DJJ;
  ff_element DFJ;
  ff_element tmp, tmp2, tmp3;
  
  ff_element sum, sum_k, sum_k2;
  
  // compute values for  E4(tau) and E6(tau)
  
  divide(E4_tau, Ea, ff_element(-3));
  divide(E6_tau, Eb, ff_element(-2));
  
  if (E4_tau.is_zero()) 
    lidia_error_handler ("eco_prime", "tildeEf()::E4_tau == 0,"
			 " division by zero");

  // compute value for j-invariant
  
  compute_jinv (j_tau, Ea, Eb);
  
  // first derivation of j at tau
  // j'(tau) = - j(tau) * E6(tau) / E4(tau)
  
  divide   (j_tau_1, E6_tau, E4_tau);
  multiply (j_tau_1, j_tau_1, j_tau);
  negate   (j_tau_1, j_tau_1);
  
  
  // initialize vectors g and j
  // g[r] = g(tau)^r  j[k] = j(tau)^k
  
  g[0] = ff_element(1); g[1] = g_tau;
  j[0] = ff_element(1); j[1] = j_tau;
  
  for (r = 2; r <= l+1; r++) 
    multiply (g[r], g[r-1], g_tau);

  for (k = 2; k <= v; k++) 
    multiply (j[k], j[k-1], j_tau);

  // reset file to read coeff. of modular equation
  
  rc = meq_prime::reset ();
  
  if (rc) 
    lidia_error_handler ("eco_prime", "tildeEf()::Unable to reset file");

  // computation of several derivations
  // for notation see [Ma94]
  
  
  DFF = DJJ = DFJ = ff_element(0);
  DF  = DJ  = ff_element(0);
  g_tau_1num = g_tau_1den = ff_element(0);

  for (r = 0; r <= l+1; r++) {
    rc = read_row (p);
    
    if (rc) 
      lidia_error_handler ("eco_prime", "tildeEA()::Unable to read row");

	// compute value for g'(tau)
	
	for (k = 1; k <= v; k++) {
	  multiply (tmp2, row[k], ff_element(k));
	  multiply (tmp2, tmp2, g[r]);
	  multiply (tmp2, tmp2, j[k-1]);
	  
	  add      (g_tau_1num, g_tau_1num, tmp2);
	}
    
    if (r > 0)
      for (k = 0; k <= v; k++) {
	multiply (tmp2, row[k], ff_element(r));
	multiply (tmp2, tmp2, g[r-1]);
	multiply (tmp2, tmp2, j[k]);
	
	add      (g_tau_1den, g_tau_1den, tmp2);
      }
    
    
    // sum      = \sum_{k=0}^{v}     a_{r, k} g(tau)^{r} j(tau)^{k}
    // sum_k    = \sum_{k=1}^{v} k   a_{r, k} g(tau)^{r} j(tau)^{k}
    // sum_k2   = \sum_{k=1}^{v} k^2 a_{r, k} g(tau)^{r} j(tau)^{k}
    
    sum    = ff_element(0);
    sum_k  = ff_element(0);
    sum_k2 = ff_element(0);
    
    for (k = 0; k <= v; k++) {
      multiply (tmp, row[k], j[k]); 
      
      add      (sum, sum, tmp);
      
      multiply (tmp2, k, tmp);
      add      (sum_k, sum_k, tmp2);
      
      multiply (tmp2, k, tmp2);
      add      (sum_k2, sum_k2, tmp2);
    }
    
    multiply (sum, sum, g[r]);
    multiply (sum_k, sum_k, g[r]);
    multiply (sum_k2, sum_k2, g[r]);
    
    
    // compute value for DFF
    
    if (r > 0) { 
      multiply (tmp, r*r, sum);
      add      (DFF, DFF, tmp);
    }
    
    
    // compute value for DJJ
    
    add  (DJJ, DJJ, sum_k2);
    
    
    // compute value for DF
    
    if (r > 0) { 
      multiply (tmp, r, sum);
      add      (DF, DF, tmp);
    }
    
    // compute value for DJ
    
    add  (DJ, DJ, sum_k);

    // compute value for DFJ
    
    multiply (tmp, r, sum_k);
    add      (DFJ, DFJ, tmp);
  }
  
  
  // compute value for f'(tau)
  // f'(tau) = - j'(tau) * \frac{g_tau_1num}{g_tau_1den}
  
  negate   (g_tau_1, j_tau_1);
  multiply (g_tau_1, g_tau_1, g_tau_1num);
  divide   (g_tau_1, g_tau_1, g_tau_1den);
  
  
  // compute value for E2_star = (E6*DJ)/(E4*DF)
  //                          = -(s/12) E2*(tau) in [Ma94]
  
  multiply (E2_star, E6_tau, DJ);
  multiply (tmp, E4_tau, DF);
  
  if (tmp.is_zero())
    lidia_error_handler ("eco_prime", "tildeEf()::DF or E4(tau) is zero");

  divide (E2_star, E2_star, tmp);
  
  
  // compute value for E0_tau = DF/DJ
  //                         = -(s/12) E0(tau) in [Ma94]
  
  if (DJ.is_zero())
    lidia_error_handler ("eco_prime", "tildeEf()::DF is zero");
  
  divide (E0_tau, DF, DJ);
  
  
  // compute value for E0_tau_1
  //  = [ DFF*E2_star - 2*DFJ*E6/E4 + DJJ*E0*E6/E4 ] / DJ
  
  divide   (tmp, E6_tau, E4_tau);
  multiply (tmp2, E0_tau, tmp);
  multiply (E0_tau_1, DJJ, tmp2);
  
  multiply (tmp2, 2, DFJ);
  multiply (tmp2, tmp2, tmp);
  subtract (E0_tau_1, E0_tau_1, tmp2);
  
  multiply (tmp2, DFF, E2_star);
  add      (E0_tau_1, E0_tau_1, tmp2);
  
  if (DJ.is_zero())
    lidia_error_handler ("eco_prime", "tildeEf()::DJ is zero");
  
  divide (E0_tau_1, E0_tau_1, DJ);
  
  
  // compute value of E4_ltau
  // = [ E4 + 12/s*E2_star* {12*E0'/E0 + 6*E4^2/E6 - 4*E6/E4 + 12/s*E2_star} ] / l^2
  
  multiply (tmp2, ff_element(4), tmp); // tmp = E6/E4
  negate   (E4_ltau, tmp2); // then: E4_ltau  =  -4*E6/E4
  
  if (tmp.is_zero())
    lidia_error_handler ("eco_prime", "tildeEf()::E6 is zero");
  
  invert   (tmp2, tmp);
  multiply (tmp2, tmp2, E4_tau);
  multiply (tmp2, tmp2, ff_element(6));
  
  add      (E4_ltau, E4_ltau, tmp2); // E4_ltau +=  6*E4^2/E6
  
  
  if (E0_tau.is_zero())
    lidia_error_handler ("eco_prime", "tildeEf()::E0_tau is zero");
  
  divide   (tmp2, E0_tau_1, E0_tau);
  multiply (tmp, ff_element(12), tmp2);
  
  add      (E4_ltau, E4_ltau, tmp); // E4_ltau +=  12*E0'/E0
  
  multiply (tmp, ff_element(12/s), E2_star);
  
  add      (E4_ltau, E4_ltau, tmp); // E4_ltau +=  12/s*E2_star
  multiply (E4_ltau, E4_ltau, tmp); // E4_ltau *=  12/s*E2_star
  add      (E4_ltau, E4_ltau, E4_tau); // E4_ltau +=  E4_tau
  
  tmp = ff_element(l*l);
  
  if (tmp.is_zero())
    lidia_error_handler ("eco_prime", "tildeEf()::l^2 is zero");

  divide (E4_ltau, E4_ltau, tmp); // E4_ltau /= l^2
  
  
  // compute value for delta_tau = E4^3/j_tau
  
  multiply (tmp, E4_tau, E4_tau);
  multiply (tmp, tmp, E4_tau);
  
  if (j_tau.is_zero())
    lidia_error_handler ("eco_prime", "tildeEf()::j_tau is zero");
  
  divide (delta_tau, tmp, j_tau);
  
  
  // compute value for delta_ltau = delta_tau/ f_tau^{12/s}
  // with f_tau = l^s / g_tau
  //            => delta_ltau = delta_tau * g_tau^{12/s} / l^12
  
  power     (tmp, g_tau, static_cast<long>(12/s));
  power     (tmp2, ff_element(l), static_cast<long>(12));

  multiply   (delta_ltau, delta_tau, tmp);
  
  if (tmp2.is_zero())
    lidia_error_handler ("eco_prime", "tildeEf()::l^12 is zero");
  
  divide   (delta_ltau, delta_ltau, tmp2);
  
  
  // compute value for j_ltau = E4_ltau^3 / delta_ltau
  
  multiply (tmp, E4_ltau, E4_ltau);
  multiply (tmp, tmp, E4_ltau);
  
  if (delta_ltau.is_zero())
    lidia_error_handler ("eco_prime", "tildeEf()::delta_ltau is zero");
  
  divide (j_ltau, tmp, delta_ltau);
  
  
  // perform the transformation tau -> -1/ltau
  // thus, we see: g(tau) -> l^s/g(tau) = f(tau)
  //              j(tau) -> j(l tau)
  // then, f'(tau) = - f(tau) g'(tau) / g(tau)
  
  power  (tmp, ff_element(l), static_cast<long>(s));
  divide (f_tau, tmp, g_tau);
  
  multiply (tmp, f_tau, g_tau_1);
  divide   (tmp, tmp, g_tau);
  negate   (f_tau_1, tmp);
  
  
  // initialize vectors g and j
  // g[r] = f(tau)^r  j[k] = j(l tau)^k
  
  
  g[0] = ff_element(1); g[1] = f_tau;
  j[0] = ff_element(1); j[1] = j_ltau;

  for (r = 2; r <= l+1; r++) 
    multiply (g[r], g[r-1], f_tau);

  for (k = 2; k <= v; k++) 
    multiply (j[k], j[k-1], j_ltau);


  // reset file to read coeff. of modular equation
  
  rc = meq_prime::reset ();
  
  
  // compute value for j_ltau_1
  
  j_ltau_1num = j_ltau_1den = ff_element(0);

  for (r = 0; r <= l+1; r++) {
    rc = read_row (p);
    
    if (rc)
      lidia_error_handler ("eco_prime", "tildeEA()::Unable to read row");
    
    // compute value for j'(l tau)

    if (r > 0)
      for (k = 0; k <= v; k++) {
	multiply (tmp2, ff_element(row[k]), ff_element(r));
	multiply (tmp2, tmp2, g[r-1]);
	multiply (tmp2, tmp2, j[k]);

	add      (j_ltau_1num, j_ltau_1num, tmp2);
      }
    
    for (k = 1; k <= v; k++) {
      multiply (tmp2, ff_element(row[k]), ff_element(k));
      multiply (tmp2, tmp2, g[r]);
      multiply (tmp2, tmp2, j[k-1]);

      add      (j_ltau_1den, j_ltau_1den, tmp2);
    }	 
  }
  
  multiply (j_ltau_1num, j_ltau_1num, f_tau_1);
  multiply (j_ltau_1den, j_ltau_1den, ff_element(- static_cast<long>(l)));

  if (j_ltau_1den.is_zero())
    lidia_error_handler ("eco_prime", "tildeEf::denominator of j'(l tau) is zero");
  
  divide (j_ltau_1, j_ltau_1num, j_ltau_1den);
  
  
  // compute value for E6(l tau)
  // E6(l tau) = - E4(ltau) * j'(l tau) / j(l tau)
  
  multiply   (tmp, E4_ltau, j_ltau_1);
  
  if (j_ltau.is_zero())
    lidia_error_handler ("eco_prime", "tildeEf::j(l tau) is zero");
  
  divide (tmp2, tmp, j_ltau);
  negate (E6_ltau, tmp2);
  
  
  // compute parameters P1, tildeEa, and tildeEb
  // P1 = -l/2 E2*(tau) = l/2 12/s E2_star
  // tildeEa == -3 * l^4 * E4(l tau)
  // tildeEb == -2 * l^6 * E6(l tau)
  
  tmp = tmp2 = ff_element(l);
  
  square(tmp, tmp); square(tmp, tmp);
  power (tmp2, tmp2, static_cast<long>(6));
  
  multiply (P1, ff_element(l * 6 / s), E2_star);
  multiply (tildeEa, ff_element(-3) * tmp, E4_ltau);
  multiply (tildeEb, ff_element(-2) * tmp2, E6_ltau);

  delete[] g; delete[] j;
  
  return 0;
}


#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
