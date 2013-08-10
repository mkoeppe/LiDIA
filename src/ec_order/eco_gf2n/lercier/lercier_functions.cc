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
//      $Id$
//
//	Author	:
//	Changes	: See CVS log
//
//==============================================================================================



#ifdef HAVE_CONFIG_H
# include	"config.h"
#endif
#include	"LiDIA/eco_gf2n.h"
#include	"LiDIA/sort_vector.h"
#include	"LiDIA/timer.h"
#ifdef DEBUG
#include	<cassert>
#endif



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif


const int length_for_guessed = 300;
const int max_number_of_guessed = 4;

const char eco_gf2n::PLAIN_LERCIER = 0;
const char eco_gf2n::WITH_VELU = 1;

//****************************************************************
//  This function computes the vector of p_k as in Lercier's paper.
//  The notation is similar to that paper, see algorithm there.
//
// compute p_k[i] (0 <= i <=d) using Lercier´s algorithm          
// 
// strategy = eco_gf2n::PLAIN_LERCIER (without optimizing a la Velu)       
//          = eco_gf2n::WITH_VELU     (use practical improvement based on Velu)
//****************************************************************

void eco_gf2n::compute_sqrtP(gf2n *& p_k, const power_tab<gf2n> & alpha_pow,
			     const gf2n & sqrt4_alpha, 
			     const gf2n & sqrt4_beta, 
			     const gf2n & beta,
			     lidia_size_t d, char strategy)
{
  mv_poly * p_mv;       // coefficients p_i of P(X)=\sum_{i=0}^{d} p_i^2 X^i
  ff_element temp, temp2;
  ff_element alpha; 
  lidia_size_t Anz_bin = 0;   
  lidia_size_t i = 1, j, k;
  
  p_mv = new mv_poly[d+1];        
  alpha_pow.get_power(alpha, 1);

  p_k[d].assign_one();              // the highest coefficient
  p_mv[d] = mv_term(p_k[d], 0);  

  add(p_k[d-1], alpha, beta);      // the second highest coefficient
  p_mv[d-1] = mv_term(p_k[d-1], 0);

  if (d >= 2)                      // coefficient p[d-2]
    {
      square(temp, p_k[d-1]);
      square(temp, temp);
      add(p_k[d-2], temp, alpha*p_k[d-1]);
      if (d & 1)
	{
	  alpha_pow.get_power(temp, 2);
	  add(p_k[d-2], p_k[d-2], temp);
	}
      p_mv[d-2] = mv_term(p_k[d-2], 0);
    }

  if (d > 2)    // coefficient p_k[0]
    {
      alpha_pow.get_power(temp, 2*d - 1);
      add (p_k[0], alpha, p_k[d-1]);
      multiply(temp, p_k[0], temp);
      sqrt(p_k[0], sqrt(temp));
      p_mv[0] = mv_term(p_k[0], 0);
    }

  if (d <= 3)
    {
      delete[] p_mv;
      return;
    }

  //---------------------------------------------------------
  //               * *  PHASE I  * * 
  //
  
     
  ff_element b_K;             
  mv_poly    c_K;
  mv_poly    C;
  mv_poly    Tr;
  mv_poly    sum1, sum2, temp_mp1, temp_mp2;
  ff_element temp_pow;
  lidia_size_t K = 1;
  ff_element sqrt_alpha, sqrt_beta;
  

  square(sqrt_alpha, sqrt4_alpha);
  square(sqrt_beta, sqrt4_beta);


  while (1) 
    {
      if (strategy == eco_gf2n::WITH_VELU && K >= 3 && d >= 7 && (K & 1))
	{
	  // Velu Improvement is possible!
	  // compute p_mv[K] without introducing new binary variable 

	  velu_improvement(p_mv, K, d, alpha_pow, sqrt_alpha, 
			   sqrt_beta);
	}
      else
	{
	  // compute p_mv[K] with Lercier's base algorithm:
	  // compute b_K ..
	  
	  alpha_pow.get_power(temp, 2*K);
	  multiply(temp, temp, sqrt4_alpha); 

	  alpha_pow.get_power(temp2, (d + 2*K) / 2);
	  if ((d + 2*K) & 1)
	    multiply(temp2, temp2, sqrt(alpha));
	  multiply(temp2, temp2, sqrt4_beta);
	  divide(b_K, temp2, temp);

	  // compute c_K
	  
	  multiply(sum1, p_mv[0], p_mv[d-K]);
	  for (i = 1; i < K; i++)
	    {
	      alpha_pow.get_power(temp_pow, i);

	      multiply(temp_mp1, p_mv[i], p_mv[d-K+i]);
	      multiply(temp_mp1, temp_pow, temp_mp1);
	      sum1 += temp_mp1;
	    }

	  square(sum1, sum1);
	  multiply(sum1, sqrt4_alpha, sum1);
	  sum2.assign_zero(); 
	  
	  for (i = 1; i <= K/2; i++)
	    {
	      if (n_over_k(d - K + 2*i, i))
		sum2 += p_mv[K - 2*i];
	    }
	  multiply(sum2, temp2, sum2);
	  add(c_K, sum1, sum2);
	  temp.invert();
	  multiply(c_K, temp, c_K);

	  // compute C
	  
	  square (temp, b_K);
	  temp.invert();
	  multiply(C, c_K, temp);

	  // compute trace 

	  C.trace_computation(Tr);
	  
	  // check Tr(C)

	  if (Tr.is_one())
	    {
	      lidia_error_handler("compute_sqrtP",
				  "Trace c_K/b_K^2 is one !!");
	      return;
	    }

	  if(!Tr.is_zero())
	      {
		k = Anz_bin - 1;
		
		while (!Tr.has_var_k(k))
		  k--;
		
		Tr.split_x_k(c_K, k);
		
                for (j = 1; j < K; j++)
		  {
		    if (p_mv[j].has_var_k(k))
		      substitute(p_mv[j], Tr, c_K, k);
		    
		    if (p_mv[d - 2*j - 1].has_var_k(k))
		      substitute(p_mv[d- 2*j - 1], Tr, c_K, k); 
		    
		    if (p_mv[d - 2*j - 2].has_var_k(k))
		      substitute(p_mv[d - 2*j - 2], Tr, c_K, k);
		  }
		substitute(C, Tr, c_K, k);
		
		if (k == Anz_bin - 1)
		  Anz_bin --;
	      }
	  
	  // now determine p_mv[K]

          const mv_list_entry *Term;
          mv_term P_mu;
          ff_element Cc_mu, Pc_mu;
          
          p_mv[K].assign_zero();
          
          for (Term = C.get_first(); Term != NULL ; 
               Term = Term->next)  
            {
              Cc_mu = Term->term.get_coeff();
              
              if (!solve_quadratic_with_table(Pc_mu, ff_element(1), Cc_mu,
                  true))
                {
		  lidia_error_handler("compute_sqrtP",
				      "Trace is null, but no roots found.");
		  return;
		}
	      
#ifdef DEBUG
              assert((Pc_mu*Pc_mu + Pc_mu + Cc_mu).is_zero());
#endif
              P_mu.assign_coeff(Pc_mu);
              P_mu.assign_var(Term->term.get_var());
              p_mv[K] += P_mu;

            } // end for
          
          // final value: p_mv[K] = b_K * X_{K-1} + b_K * p_mv[K]
	  
          P_mu.assign(b_K, bigint (1) << ((long) Anz_bin));
          multiply(p_mv[K], b_K, p_mv[K]);
          add(p_mv[K], P_mu, p_mv[K]); 
	  Anz_bin ++;

#ifdef DEBUG
          assert((p_mv[K] * p_mv[K] + b_K * p_mv[K] + c_K).is_zero());
#endif
	} // end else (velu_improvement was not possible)
      
      if (K >= d - 2*K - 1)
	break;
      
      // compute coefficient p_k[d - 2K - 1] via equation system ...
      
      eqnsys2n_1(p_mv, K, alpha_pow, d);

      if (K >= d - 2*K - 2)
	break;

      // compute coefficient p_k[d - 2K - 2] via equation system ...
	  
      eqnsys2n_2(p_mv, K+1, alpha_pow, d);            

      if (K >= d - 2*K - 3)
	break;

      // check length, if too large, call version that guesses some
      // variables
      
      if (p_mv[K].length() > length_for_guessed)
	{
	  compute_sqrtP_guess(K, d, p_mv, p_k, alpha_pow, 
			      sqrt4_alpha, sqrt4_beta, strategy);
	  mv_poly::delete_free_list();
	  delete[] p_mv;
	  return;
	}

      K++;
    }
  
  bigint  pi_values; 
  mv_poly * EQS;
  int number_eqs = Anz_bin + 3;


  EQS = new mv_poly[number_eqs];

  // build equations by using the remaining non-trivial equations 
  // from the linear equation system over GF(2^n)  

  switch (d % 3)
    {
    case 0: i = d/3 + 1; break;
    case 1: i = (d-1)/3 + 1; break;
    case 2: i = (d+1)/3 + 1; break;
    }
  
  for (j = 0; i <= (d-1)/2 && j < number_eqs; i++, j++)
    {
      build_eqn1(EQS[j], p_mv, i, alpha_pow, d);
      if (EQS[j].is_zero())
	j--;
    }

  switch (d % 3)
    {
    case 0: i = d/3 ; break;
    case 1: i = (d-1)/3; break;
    case 2: i = (d-2)/3; break;
    }
  
  for ( ; i <= d/2 && j < number_eqs; i++, j++)
    {
      build_eqn2(EQS[j], p_mv, i, alpha_pow, d);
      if (EQS[j].is_zero())
	j--;
    }

  // get values for all binary variables Pi by solving the equation system

  pi_values = solve_eqs(EQS, Anz_bin, j, 0, 0);
  if (pi_values == -1)
    {
      lidia_error_handler("compute_sqrt_P", "solve_eqs found no solution");
      return;
    }

  // evaluate all p_mv´s to get values in GF(2^n)

  for (i = 1; i < d - 2; i++)
    evaluate(p_k[i], p_mv[i], pi_values);
  
  delete[] EQS;
  delete[] p_mv; 
  mv_poly::delete_free_list();
}




//****************************************************************
// compute p_mv[d - 2k - 1] using the linear equation system    
// over GF(2^n) -- (22) in Lercier's paper
//
//   0 <= k <= floor((d-1)/2)
//
// p_mv[d - 2*k - 1]  = 
// (alpha^{2d -4k -1})^{-1} (p_mv[k]^4 + alpha^{2d - 4k - 1} \sum_{i=1}^k ...
//                                     + alpha^{2d - 4k} \sum_{i=0}^k ...

void eco_gf2n::eqnsys2n_1(mv_poly *& p_mv, lidia_size_t k, 
			  const power_tab<ff_element> & ap, lidia_size_t d)
{
  
  mv_poly temp;
  mv_poly sum1, sum2;
  ff_element pow;
  lidia_size_t i;

  sum1.assign_zero();
  sum2 = p_mv[d - 2*k];  

  for(i = 1; i <= k; i++)
    {
      if (n_over_k(d - 2*k + 2*i - 1, i)) 
       {
	 ap.get_power(pow, 2*i);
	 multiply(temp, pow, p_mv[d - 2*k - 1 + 2*i]);
	 sum1 += temp;
	}
      
      if (n_over_k(d - 2*k + 2*i, i))
	{
	  ap.get_power(pow, 2*i);
	  multiply(temp, pow, p_mv[d - 2*k + 2*i]);
	  sum2 += temp;
	}
    }

  ap.get_power(pow, 2*d - 4*k -1);
  multiply(sum1, pow, sum1);

  ap.get_power(pow, 2*d - 4*k);
  multiply(sum2, pow, sum2);

  square(temp, p_mv[k]);
  square(temp, temp);

  add(p_mv[d - 2*k - 1], sum1 + sum2, temp);

  ap.get_power(pow, 2*d - 4*k - 1);
  pow.invert();
  multiply(p_mv[d - 2*k - 1], pow, p_mv[d - 2*k - 1]);
}


//**************************************************************
// compute p_mv[d - 2k] using the linear equation system        
// over GF(2^n) -- (23) in Lerciers paper                             
//
//   1 <= k < floor(d/2)
//
// p_mv[d - 2*k] = p_mv[d - k]^4 + alpha *\sum_{i=0}^{k-1} ...
//                 + \sum_{i=1}^k ...


void eco_gf2n::eqnsys2n_2(mv_poly *& p_mv, lidia_size_t k, 
			  const power_tab<ff_element> & ap, lidia_size_t d)
{
  mv_poly sum1, sum2, temp;
  ff_element pow;
  lidia_size_t i;

  ap.get_power(pow, 1);
  multiply(sum1, pow, p_mv[d + 1 - 2*k]);
  sum2.assign_zero();

  for (i = 1; i < k; i++)
    {
      if (n_over_k (d + 1 -2*k + 2*i, i))
	{
	  ap.get_power(pow, 2*i + 1);
	  multiply(temp, pow, p_mv[d + 1 - 2*k + 2*i]);
	  sum1 += temp;
	}

      if (n_over_k(d - 2*k + 2*i, i))
	{
	  ap.get_power(pow, 2*i);
	  multiply(temp, pow, p_mv[d - 2*k + 2*i]);
	  sum2 += temp;
	}
    } 

  if (n_over_k(d, i))
    {
      ap.get_power(pow, 2*k);
      multiply(temp, pow, p_mv[d]);
      sum2 += temp;
    }

  square(temp, p_mv[d-k]);
  square(temp, temp);
  add(p_mv[d - 2*k], sum1 + sum2, temp);
}



/****************************************************************/
/* compute n over k modulo 2 and return it.                     */ 
/*                                                              */
/****************************************************************/

int eco_gf2n::n_over_k(lidia_size_t n, lidia_size_t k)
{
  if (!k)
    return 1;

  if (k == 1)
    if (n & 1)
      return 1;
    else
      return 0;

  lidia_size_t i, j;
  lidia_size_t anz_2 = 0;

  i = n-k+1;
  if (i & 1)
    i++;

  for(; i <= n; i += 2)
    {
      j = i >> 1;
      anz_2 ++;

      while (!(j & 1))
	{
	  j >>= 1;
	  anz_2 ++;
	}
    }

  for(i = 2; i <= k; i += 2)
    {
      j = i;
      while (!(j & 1))
	{
	  j >>= 1;
	  anz_2 --;
	}
    }
  if (anz_2)
    return 0;
  else
    return 1;
}



//*************************************************************************
// compute a non-trivial equation of multivariate polynomials  
// using the linear equation system over GF(2^n), stored in a fraction. 
//
// based on Lercier's formula (22)
//
// EQ = p_mv[k]^4 + alpha^{2d - 4k - 1} \sum_{i=0}^k ...
//                + alpha^{2d - 4k} \sum_{i=0}^k ...


void eco_gf2n::build_eqn1(mv_poly & EQ, mv_poly *& p_mv, 
			  lidia_size_t k, const power_tab<ff_element> & ap, 
			  lidia_size_t d)
{
  mv_poly temp;
  mv_poly sum1, sum2;
  ff_element pow;
  lidia_size_t i;

  sum1 = p_mv[d - 2*k - 1];
  sum2 = p_mv[d - 2*k];

  for(i = 1; i <= k; i++)
    {
      if (n_over_k(d - 2*k - 1 + 2*i, i))
	{
	  ap.get_power(pow, 2*i);
	  multiply(temp, pow, p_mv[d - 2*k - 1 + 2*i]);
	  sum1 += temp;
	}
      if (n_over_k(d - 2*k + 2*i, i))
	{
	  ap.get_power(pow, 2*i);
	  multiply(temp, pow, p_mv[d - 2*k + 2*i]);
	  sum2 += temp;
	}
    }

  ap.get_power(pow, 2*d - 4*k -1);
  multiply(sum1, pow, sum1);

  ap.get_power(pow, 2*d - 4*k);
  multiply(sum2, pow, sum2);

  square(temp, p_mv[k]);
  square(temp, temp);

  add(temp, sum1 + sum2, temp);
  EQ.assign(temp);
}


//***************************************************************************
// compute a non-trivial equation of multivariate polynomials                
// using the linear equation system over GF(2^n) (based on Lercier's
// formula (23))
//
// EQ = p_mv[d-k]^4 + alpha \sum_{i=0}^{k-1} ...
//                        + \sum_{i=0}^k ...

void eco_gf2n::build_eqn2(mv_poly & EQ, mv_poly *& p_mv, 
			  lidia_size_t k, const power_tab<ff_element> & ap, 
			  lidia_size_t d)
{
  mv_poly temp;
  mv_poly sum1, sum2;
  ff_element pow;
  lidia_size_t i;

  sum1.assign_zero();
  sum2.assign_zero();

  for (i = 0; i < k; i++)
    {
      if (n_over_k(d + 1 - 2*k + 2*i, i))
	{
	  ap.get_power(pow, 2*i+1);
	  multiply(temp, pow, p_mv[d + 1 - 2*k + 2*i]);
	  sum1 += temp;
	}
      if (n_over_k(d - 2*k + 2*i, i))
	{
	  ap.get_power(pow, 2*i);
	  multiply(temp, pow, p_mv[d - 2*k + 2*i]);
	  sum2 += temp;
	}
    }

  if (n_over_k(d, i))
    {
      ap.get_power(pow, 2*k);
      multiply(temp, pow, p_mv[d]);
      sum2 += temp;
    }
  
  square(temp, p_mv[d-k]);
  square(temp, temp);
  add(temp, sum1 + sum2, temp);
  EQ.assign(temp);
}


 


//****************************************************************************
//    Lerciers algorithm; Practical Improvement: V'elus Equation            
//    see Lercier-Paper (33)
//                                                                          
//   input:  p_mv[i] as multivariate polynomials for 
//           0 <= i < K and d - 2K <= i <= d,       
//           K odd, d >= 7, alpha-powers, sqrt(alpha), sqrt(beta)       
//   output: p_mv[K] without introduction of new binary variable       
//
//***************************************************************************

bool eco_gf2n::velu_improvement (mv_poly * & p_mv, lidia_size_t K, 
				 lidia_size_t d, 
				 const power_tab<ff_element> & ap, 
				 const ff_element & sqrt_alpha,  
				 const ff_element & sqrt_beta)
{
  lidia_size_t  i, r;
  ff_element h, h2;
  mv_poly temp, temp2;

  r = d - (K+1)/2;   // = d - ceil(K/2)

  ap.get_power(h, 2*r - d);
  multiply(h2, sqrt_alpha, sqrt_beta);
  multiply(h, h, h2);
  square(temp, p_mv[d - (K+1)/2]);
  multiply(temp, temp, h);

  ap.get_power(h, 1);
  square(temp2, p_mv[(K+1)/2]);
  multiply(temp2, temp2, h);
  add(temp, temp, temp2);

  if (d & 1) 
    {
      for (i = 0; i < K/2; i++)
	{
	  multiply(temp2,  p_mv[K-1-2*i], p_mv[2*i+1]);
	  add(temp, temp, temp2);
	}
    }
  else
    {
      for (i = 1; i <= (K-1)/2; i++)
	{
	  multiply(temp2,  p_mv[K - 2*i], p_mv[2*i]);
	  add(temp, temp, temp2);
	}
    }
  if (!p_mv[0].is_zero())
    h = p_mv[0].get_first()->term.get_coeff();
  else
    return false;

  if (h.is_zero())
    return false;
  h.invert();
  multiply(p_mv[K], temp, h);
  return true;
}



//---------------------------------------------------------------------
// Solve the equations that are stored in EQS.
// number_eq is number of equations, number_var is number of variables
// bit_mask is mask that stores the values for already guessed variables.
// number_of_guessed_var stores the number of such variables.
//
// Function returns -1 if no solution cound be found, otherwise a bigint
// that represents the solutions to the X_i


bigint eco_gf2n::solve_eqs(mv_poly* &EQS, lidia_size_t number_var, 
			   lidia_size_t number_eq, lidia_size_t bit_mask,
			   lidia_size_t number_of_guessed_var)
{  
  lidia_size_t l, i, j, max_len = 0;
  bigint pi_values;   
  bigint hbi, variables;
  mv_poly * EQS_d;
  lidia_size_t *eq_for_Xl;

  if (number_var > number_eq)
    {
      lidia_error_handler("solve_eqs","less equations than variables");
      return -1;
    }

  EQS_d = new mv_poly[number_eq];
  eq_for_Xl = new lidia_size_t[number_var];

  for (i = 0; i < number_var; i++)
    eq_for_Xl[i] = -1;

  // we assume that the EQS's are already reduced and have no variables
  // that we guessed and substituted

  for(i = number_eq - 1, l = number_var - 1 ; l >= 0 ; l--, i--)
    {
      max_len = 0;
      j = i;
      
      while (!EQS[j].has_var_k(l) && j > 0)
	j --;

      if (j == 0 && !EQS[0].has_var_k(l))
	{
	  i++;
	  continue;
	}

      if (j != i)
	swap(EQS[i], EQS[j]);

      eq_for_Xl[l] = i;      // now EQS[i] has variable l
      EQS[i].split_x_k (EQS_d[i], l);

      for (j = i-1; j >= 0; j--)
	{
	  substitute(EQS[j], EQS[i], EQS_d[i], l);
	  if (EQS[j].length() > max_len)
	    max_len = EQS[j].length();
	}

      if (max_len > length_for_guessed)
	{
	  mv_poly * EQS_copy;
	  int begin_try = i - 1, begin_var = l - 1;
	  bool error = false;
	  mv_poly gn, gd;
	  int number_of_local_guessed_var;
	  int num, var;

	  EQS_copy = new mv_poly [begin_try + 1];
	 
	  for (j = begin_try; j >= 0; j--) // save the already computed eqs
	    EQS_copy[j].assign(EQS[j]);

	  if (number_of_guessed_var > 0)
	    {
	      number_of_local_guessed_var = number_of_guessed_var + 1;
	      num = (1 << number_of_local_guessed_var) - 1;
	      variables.assign(num);
	      pi_values.assign(bit_mask);
	    }
	  else
	    {
	      number_of_local_guessed_var = 1;
	      num = 1;
	      variables.assign_one();
	      pi_values.assign_zero();
	    }
	  
	  gn.assign(EQS[begin_try]);  
	  evaluate(gn, pi_values, variables);
	  while(gn.length() > length_for_guessed 
		&& number_of_local_guessed_var <= max_number_of_guessed)    
	    {
	      number_of_local_guessed_var ++;
	      num = (num << 1) + 1;
	      variables.assign(num);
	      evaluate(gn, pi_values, variables);
	    }
	  
	  for (var = 0; var <= num; var++)
	    {
	      if (number_of_guessed_var > 0)
		{
		  int h;
		  
		  h = (var ^ bit_mask) & ((1 << number_of_guessed_var) - 1);
		  if (h != 0)
		    continue;
		}

	      if (error)
		for (j = begin_try; j >= 0; j--)
		  EQS[j].assign(EQS_copy[j]);

	      error = false;
	      variables.assign(num);
	      pi_values.assign(var);
	      
	      for (j = begin_try; j >= 0; j--)
		evaluate(EQS[j], pi_values, variables);

	      for (j = begin_var; j >= 0; j--)
		eq_for_Xl[j] = -1;

	      for(i = begin_try, l = begin_var; l >= 0 && !error; l--, i--)
		{
		  j = i;
		  
		  while (!EQS[j].has_var_k(l) && j > 0)
		    j --;
		  
		  if (j == 0 && !EQS[0].has_var_k(l))
		    {
		      i++;
		      continue;
		    }

		  if (j != i)
		    swap(EQS[i], EQS[j]);

		  eq_for_Xl[l] = i;  // EQS[i] stores X_l

		  EQS[i].split_x_k (EQS_d[i], l);

		  for (j = i - 1; j >= 0; j--)
		    {
		      substitute(EQS[j], EQS[i], EQS_d[i], l);
		      if (EQS[j].length() > max_len)
			max_len = EQS[j].length();
		    }

		  while (max_len > length_for_guessed && 
			 number_of_local_guessed_var <= max_number_of_guessed)
		    {
		      max_len = 0;
		      number_of_local_guessed_var ++;
		      num = (num << 1) + 1;
		      variables.assign(num);

		      for (j = i - 1; j >= 0; j--)
			{
			  evaluate(EQS[j], pi_values, variables);
			  if (EQS[j].length() > max_len)
			     max_len = EQS[j].length();
			}
		    }
		}

	      for (j = i; j >= 0 && !error; j--)  // check the remaining eqs
		{
		  if (!EQS[j].is_zero())
		    {
		      error = true;
		      break;
		    }
		}
	      if (error)
		continue;

	      for (l = variables.bit_length(); l < number_var; l++)
		{ 
		  i = eq_for_Xl[l];
		  if (i == -1)
		    {
		      variables.multiply_by_2();
		      inc(variables);
		      continue;
		    }
		  
		  evaluate(EQS[i], pi_values, variables);
		  shift_left(hbi, bigint(1), l);
		  
		  if (EQS[i].is_zero())  // X_l = 0
		    add(variables, variables, hbi);
		  else
		    {
		      evaluate(EQS_d[i], pi_values, variables);
		      if (EQS[i] == EQS_d[i])  // X_l = 1
			{
			  add(variables, variables, hbi);
			  add(pi_values, pi_values, hbi);
			}
		      else
			{
			  error = true;
			  break;
			}
		    }
		}
	      
	      if (!error)
		{
		  delete[] EQS_copy; delete[] EQS_d; delete[] eq_for_Xl;
		  return pi_values;
		}
	    }
	  delete[] EQS_copy; delete[] EQS_d; delete[] eq_for_Xl;
	  return bigint(-1);
	  } 
    }

  //---------------------------------------------------------------
  // no new variable was guessed, we can use the already reduced
  // equations to compute the values for X_i

  for (j = i; j >= 0; j--)
    {
      if (!EQS[j].is_zero())
	{
	  delete[] EQS_d;
	  delete[] eq_for_Xl;
	  return -1;
	}
    }
  
  variables.assign_one();

  l = 0;
  while (eq_for_Xl[l] == -1)
    {
      variables.multiply_by_2();
      inc(variables);
      l ++;
    }
  
  i = eq_for_Xl[l];

  if (EQS[i].is_zero())
    pi_values.assign_zero();
  else
    if (EQS[i] == EQS_d[i])
      shift_left(pi_values, bigint(1), l);
    else
      {
	delete[] EQS_d;
	delete[] eq_for_Xl;
	return -1;
      }
  
  for (l++; l < number_var; l++)
    { 
      i = eq_for_Xl[l];
      if (i == -1)
	{
	  variables.multiply_by_2();
	  inc(variables);
	  continue;
	}

      evaluate(EQS[i], pi_values, variables);
      shift_left(hbi, bigint(1), l);

      if (EQS[i].is_zero())  // X_l = 0
	add(variables, variables, hbi);
      else
	{
	  evaluate(EQS_d[i], pi_values, variables);
	  if (EQS[i] == EQS_d[i])  // X_l = 1
	    {
	      add(variables, variables, hbi);
	      add(pi_values, pi_values, hbi);
	    }
	  else
	    {
	      delete[] EQS_d;
	      delete[] eq_for_Xl;
	      return bigint(-1);
	    }
	}
    }
  
  if (number_of_guessed_var > 0)
    bitwise_or(pi_values, pi_values, 
       (bit_mask & ((1 << number_of_guessed_var) - 1)));

  delete[] EQS_d;
  delete[] eq_for_Xl;
  return pi_values;
}


//================================================================
//	      
	      
void eco_gf2n::compute_sqrtP_guess(int current_K, 
				   int d, mv_poly * & p_mv,
				   gf2n * & p_k,
				   const power_tab<gf2n> & alpha_pow, 
				   const gf2n & sqrt4_alpha, 
				   const gf2n & sqrt4_beta,
				   char strategy)
{
  mv_poly * p_mv_correct;
  mv_poly *EQS;
  mv_poly gn;
  int i;
  ff_element temp, temp2;
  ff_element b_K;             
  mv_poly    c_K;
  mv_poly    C;
  mv_poly    Tr;
  mv_poly    sum1, sum2, temp_mp1, temp_mp2;
  ff_element temp_pow;

  // save the already computed results

  p_mv_correct = new mv_poly[d+1];
  EQS = new mv_poly[d+1];

  p_mv_correct[d].assign(p_mv[d]);

  for (i = 0; i <= current_K; i++)
    {
      p_mv_correct[i].assign(p_mv[i]);
      p_mv_correct[d-2*i-1].assign(p_mv[d - 2*i - 1]);
      p_mv_correct[d-2*i-2].assign(p_mv[d - 2*i - 2]);
    }

  // find the smallest number of guessed variables

  int number_of_guessed_var = 1;
  int num = 1, K, j, k, Anz_bin;
  bigint variables(1), pi_values;
  bool error = false;
  int var;
  ff_element sqrt_alpha, sqrt_beta;
  
  square(sqrt_alpha, sqrt4_alpha);
  square(sqrt_beta, sqrt4_beta);

  gn.assign(p_mv[current_K]);
  evaluate(gn, bigint(0), variables);
  
  while(gn.length() > length_for_guessed 
	&& number_of_guessed_var <= max_number_of_guessed)    
    {
      number_of_guessed_var ++;
      num = (num << 1) + 1;
      variables.assign(num);
      evaluate(gn, bigint(0), variables);
    }

  for (var = 0; var <= num; var++)  // main loop
    {
      pi_values.assign(var);
      
      if (error)  // copy already computed values
	{
	  p_mv[d].assign(p_mv_correct[d]);

	  for (j = 0; j <= current_K; j++)
	    {
	      p_mv[j].assign(p_mv_correct[j]);
	      p_mv[d - 2*j - 1].assign(p_mv_correct[d - 2*j - 1]);
	      p_mv[d - 2*j - 2].assign(p_mv_correct[d - 2*j - 2]);
	    }
	  error = false;
	}
      
      evaluate(p_mv[d], pi_values, variables);
      Anz_bin = 0;

      for (j = 0; j <= current_K; j++)
	{
	  evaluate(p_mv[j], pi_values, variables);
	  evaluate(p_mv[d - 2*j - 1], pi_values, variables);
	  evaluate(p_mv[d - 2*j - 2], pi_values, variables);
	  Anz_bin = comparator<int>::max(Anz_bin, 
					 p_mv[j].max_index_variable());
	  Anz_bin = comparator<int>::max(Anz_bin, 
					 p_mv[d-2*j-1].max_index_variable());
	  Anz_bin = comparator<int>::max(Anz_bin, 
					 p_mv[d-2*j-2].max_index_variable());
	}
      Anz_bin ++;

      //-----------------------------------------------------------------
      // ready to start phase I
      //

      K = current_K + 1;

      while (!error) 
	{   // compute p_mv[K]: either with Velu or Lercier standard
	  
	  if (strategy == eco_gf2n::WITH_VELU && K >= 3 && d >= 7 && (K & 1))
	    {
	      // Velu Improvement is possible!
	      // compute p_mv[K] without introducing new binary variable 

	      if (!velu_improvement(p_mv, K, d, alpha_pow, sqrt_alpha, 
				    sqrt_beta))
		{
		  error = true;
		  break;
		}
	    }
	  else
	    {
	      // compute p_mv[K] with Lercier's base algorithm:
	      // compute b_K ..
	      
	      alpha_pow.get_power(temp, 2*K);
	      multiply(temp, temp, sqrt4_alpha); 
	      
	      alpha_pow.get_power(temp2, (d + 2*K) / 2);
	      if ((d + 2*K) & 1)
		multiply(temp2, temp2, sqrt_alpha);
	      multiply(temp2, temp2, sqrt4_beta);
	      divide(b_K, temp2, temp);

	      // compute c_K
	      
	      multiply(sum1, p_mv[0], p_mv[d-K]);
	      for (i = 1; i < K; i++)
		{
		  alpha_pow.get_power(temp_pow, i);
		  
		  multiply(temp_mp1, p_mv[i], p_mv[d-K+i]);
		  multiply(temp_mp1, temp_pow, temp_mp1);
		  sum1 += temp_mp1;
		}
	      
	      square(sum1, sum1);
	      multiply(sum1, sqrt4_alpha, sum1);
	      sum2.assign_zero(); 
	      
	      for (i = 1; i <= K/2; i++)
		{
		  if (n_over_k(d - K + 2*i, i))
		    sum2 += p_mv[K - 2*i];
		}
	      multiply(sum2, temp2, sum2);
	      add(c_K, sum1, sum2);
	      if (temp.is_zero())
		{
		  error = true;
		  break;
		}
	      else
		temp.invert();
	      multiply(c_K, temp, c_K);
	      
	      // compute C
	      
	      square (temp, b_K);
	      if (temp.is_zero())
		{
		  error = true;
		  break;
		}
	      else
		temp.invert();
	      multiply(C, c_K, temp);
	      
	      // compute trace 
	      
	      C.trace_computation(Tr);

	      // check Tr(C)
	      
	      if (Tr.is_one())
		{
		  error = true;
		  break;
		}

	      if(!Tr.is_zero()) 
		{
		  k = Anz_bin - 1;
		  
		  while (!Tr.has_var_k(k))
		    k--;
		  
		  Tr.split_x_k(c_K, k);
		  
		  for (j = 1; j < K; j++)
		    {
		      if (p_mv[j].has_var_k(k))
			substitute(p_mv[j], Tr, c_K, k);
		      
		      if (p_mv[d - 2*j - 1].has_var_k(k))
			substitute(p_mv[d- 2*j - 1], Tr, c_K, k); 
		      
                    if (p_mv[d - 2*j - 2].has_var_k(k))
                      substitute(p_mv[d - 2*j - 2], Tr, c_K, k);
		    }
		  substitute(C, Tr, c_K, k);
		  if (k == Anz_bin - 1)
		    Anz_bin --;
		}
  
	      // now iteratively determine p_mv[K]

	      evaluate(C, pi_values, variables);

	      const mv_list_entry *Term;
	      mv_term P_mu;
	      ff_element Cc_mu, Pc_mu;
	      
	      p_mv[K].assign_zero();
	      
	      for (Term = C.get_first(); Term != NULL ; 
		   Term = Term->next)  
		{
		  Cc_mu = Term->term.get_coeff();
		  
		  if (!solve_quadratic(Pc_mu, ff_element(1), Cc_mu))
		    {
		      error = true; break;
		    }
		  
		  P_mu.assign_coeff(Pc_mu);
		  P_mu.assign_var(Term->term.get_var());
		  p_mv[K] += P_mu;
		} // end for

	      if (error)
		break;
	      
	      // final value: p_mv[K] = b_K * X_{} + b_K * p_mv[K]

	      P_mu.assign(b_K, bigint (1) << ((long) Anz_bin));
	      multiply(p_mv[K], b_K, p_mv[K]);
	      add(p_mv[K], P_mu, p_mv[K]); 
	      Anz_bin ++;
	    } // end else (velu_improvement was not possible)

	  if (K >= d - 2*K - 1)
	    break;

	  // compute coefficient p_k[d - 2K - 1] via equation system ...
	  
	  eqnsys2n_1(p_mv, K, alpha_pow, d);

	  if (K >= d - 2*K - 2)
	    break;
	  
	  // compute coefficient p_k[d - 2K - 2] via equation system ...
          
	  eqnsys2n_2(p_mv, K+1, alpha_pow, d);   

	  if (K >= d - 2*K - 3)
	    break;

	  // check length, if length increases, add a guessed variable
	  
	  if (p_mv[K].length() > length_for_guessed)
	    {
	      number_of_guessed_var ++;
	      num = (num << 1) + 1;
	      variables.assign(num);
	      for (j = 0; j <= K; j++)
		{
		  evaluate(p_mv[j], pi_values, variables);
		  evaluate(p_mv[d - 2*j - 1], pi_values, variables);
		  evaluate(p_mv[d - 2*j - 2], pi_values, variables);
		}
	    }
	  K++;
	} // of while(1)

      if (error) 
	continue;

      // build equations by using the remaining non-trivial equations 
      // from the linear equation system over GF(2^n)  
      
      int number_eqs = Anz_bin + 3;

      switch (d % 3)
	{
	case 0: i = d/3 + 1; break;
	case 1: i = (d-1)/3 + 1; break;
	case 2: i = (d+1)/3 + 1; break;
	}
      
      for (j = 0; i <= (d-1)/2 && j < number_eqs; i++, j++)
	{
	  build_eqn1(EQS[j], p_mv, i, alpha_pow, d);
	  if (EQS[j].is_zero())
	    j--;
	}
      
      switch (d % 3)
	{
	case 0: i = d/3 ; break;
	case 1: i = (d-1)/3; break;
	case 2: i = (d-2)/3; break;
	}
      
      for ( ; i <= d/2 && j < number_eqs; i++, j++)
	{
	  build_eqn2(EQS[j], p_mv, i, alpha_pow, d);
	  if (EQS[j].is_zero())
	    j--;
	}
      number_eqs = j;

      pi_values = solve_eqs(EQS, Anz_bin, number_eqs, var, number_of_guessed_var);
      if (pi_values >= 0)
	{
	  bitwise_or(pi_values, pi_values, bigint(var));
	  break;
	}
      else
	error = true;
    }

  if (error)
    {
      lidia_error_handler("compute_sqrtP_guess","no solution found");
      return;
    }

  for (i = 1; i < d - 2; i++)
    evaluate(p_k[i], p_mv[i], pi_values);

  delete[] EQS;
  delete[] p_mv_correct;
}





#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
