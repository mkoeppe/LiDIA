//==============================================================================================
//
//      This file is part of LiDIA --- a library for computational number theory
//
//      Copyright (c) 1994--2001 the LiDIA Group.  All rights reserved.
//
//      See http://www.informatik.tu-darmstadt.de/TI/LiDIA/
//
//----------------------------------------------------------------------------------------------
//
//      $Id$
//
//      Author  :
//      Changes : See CVS log
//
//==============================================================================================


#ifdef HAVE_CONFIG_H
# include       "config.h"
#endif
#include        "LiDIA/lanczos.h"
#include        "LiDIA/random_generator.h"
#include        <ctime>
#include        <cstring>


#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif


lanczos::lanczos(const lanczos_sparse_matrix& matrix) : A(matrix)
{
  N = matrix.number_of_rows();
  M = matrix.number_of_columns();
  
  V_next = new lanczos_vector_block(M);
  V = new lanczos_vector_block(M);
  V_prev = new lanczos_vector_block(M);
  
  X = new lanczos_vector_block(M);
  dummy1_vec = new lanczos_vector_block(M > N ? M : N);
  dummy2_vec = new lanczos_vector_block(M > N ? M : N);
  
  for (size_type i = 0; i < M; i++)
    {
      random_generator rg;
      value_type random_value;
      
      rg >> random_value;
      X->put_row(i, random_value);
    }
}

  
  
  
lanczos::~lanczos()
{
  delete(V_next);
  
  // there exists an "unfindable bug" that V and V_prev have sometimes
  // the same address, double freeing leads to core dump.
  
  if (V != V_prev)  
    delete(V_prev);  
  delete(V); 
  
  delete(X);
  delete(dummy1_vec);
  delete(dummy2_vec);
}



size_t
lanczos::compute_tau2(const lanczos_small_matrix& mt,
                      value_type& S_entry,
                      lanczos_small_matrix& mw_inv) const
{
  lanczos_small_matrix In, h;
  size_type dim = 0, max = 0, i = 0, j = 0, k = 0;
  value_type temp = 0;
  size_type c[WordSize];
  size_type akt = 0, num = 0, lo = 0, hi = WordSize -1;
  
  for (i = 0; i < WordSize; i++)
    In.put_row(i, Bit_mask(i));
  h = mt;
  
  num = S_entry;
  i = WordSize;
  while (i-- > 0)
    if (num & Bit_mask(i))
      c[hi--] = i;
    else
      c[lo++] = i;
  S_entry = 0;
  
  dim = WordSize;
  
  for (j = 0; j < dim; j++) 
    {
      akt = c[j];
      k = j;
      
      while (k < dim) 
        {
          if (!((h.get_row(akt) >> (dim - akt - 1)) & 1)) 
            {
              if ((h.get_row(c[k]) >> (dim - akt - 1)) & 1) 
                {
                  temp = h.get_row(c[k]);
                  h.put_row(c[k], h.get_row(akt));
                  h.put_row(akt, temp);
                  temp = In.get_row(c[k]);
                  In.put_row(c[k], In.get_row(akt));
                  In.put_row(akt, temp);
                }
              k++;
            }
          else
            k = dim;
        }
      
      if ((h.get_row(akt) >> (dim - akt - 1)) & 1) 
        {
          S_entry += Bit_mask(akt);
          max++;
          for (k = 0; k < dim; k++) 
            {
              if (k != j) 
                {
                  // if Bit j in rows k == 1: eliminate it 
                  if ((h.get_row(c[k]) >> (dim - akt - 1)) & 1) 
                    {
                      h.put_row(c[k], h.get_row(c[k]) ^ h.get_row(akt));
                      In.put_row(c[k], In.get_row(c[k])^In.get_row(akt));
                    }
                }
            }
        }
      else 
        {
          // no pivot element in column j
          // Spalte j in T ist Linearkombination vorhergegangener Spalten in T.
          // Spalte wird verworfen
          k = j;
          while (k < dim) 
            {
              if (!((In.get_row(akt) >> (dim - akt - 1)) & 1)) 
                {
                  if ((In.get_row(c[k]) >> (dim - akt - 1)) & 1) 
                    {
                      // exchange row j and k in M
                      temp = h.get_row(c[k]);
                      h.put_row(c[k], h.get_row(akt));
                      h.put_row(akt, temp);
                      temp = In.get_row(c[k]);
                      In.put_row(c[k], In.get_row(akt));
                      In.put_row(akt, temp);
                    }
                  k++;
                }
              else k = dim;
            }
          
          // Reduction: eliminate rest of column j + dim
          for (k = 0; k < dim; k++) {
            if (k != j) {
              // if bit j in row k == 1: eliminate bit k
              if ((In.get_row (c[k]) >> (dim - akt - 1)) & 1) 
                {
                  h.put_row(c[k], h.get_row(c[k]) ^ h.get_row(akt));
                  In.put_row(c[k], In.get_row(c[k])^In.get_row(akt));
                }
            }
          }
          
          // eliminate row j
          h.put_row(akt, 0);
          In.put_row(akt, 0);
        }
    }
  
  mw_inv = In;
  
  return(max);
}



void
lanczos::compute_result()
{
  long j, i, k, h;
  value_type g;
  value_type maske1 = 0, maske2, muster1, muster2;
  lanczos_vector_block *ph1 = NULL, *phh1 = NULL;
  size_type rank = 0;
  long columns, max;
  bool fertig = false, pivot;

  dummy1_vec->clear();
  dummy2_vec->clear();

  //*****************************************************
  // determine Z as concatenation of (X-Y) and V_m
  // compute B*Z
  //*****************************************************
  
  i = 0;
  while (i < M)
    {
      j = 0;
      h = A.get_vector(i).get_number_of_entries();

      while (j < h) 
        {
          g = A.get_vector(i).get_entry(j);

          dummy2_vec->bitwise_xor(g, V->get_row(i));
          dummy1_vec->bitwise_xor(g, X->get_row(i));
          j++;
        }
      i++;
    }
  
  i = 0;
  while ((i < M) && (V->get_row(i++) == 0));
  
  if ((i != M) || (V->get_row(i-1) != 0)) {
    //*************************************************************
    // now : dummy1_vec = B * (X-Y), dummy2_vec = B * V_m
    // compute basis of zero space of BZ with Gauss' algorithm
    //*************************************************************
    
    columns = 2 * WordSize -1;
    i = M -1;
    j = columns;
    muster1 = One_mask;
    muster2 = One_mask;

    while ((j >= WordSize) && (!fertig)) {
      maske2 = Bit_mask(j-WordSize);
      
#ifdef DEBUG
      printf("Work on column %d \n", j);
      dummy1_vec->print();
      printf("       \n");
      dummy2_vec->print();
#endif
      pivot = false;
      
      while ((i >= 0) && (!(dummy1_vec->get_row(i) & muster1))
             && (!(dummy2_vec->get_row(i) & muster2)))
        i--;
      if (i == -1) {
        fertig = true;
        rank = j;
      }
      else {
        max = j;
        while ((max >= WordSize) && !pivot) {
          if (!(dummy2_vec->get_row(i) & Bit_mask(max - WordSize)))
            max --;
          else {
            pivot = true; // Pivotelement in Zeile i Spalte max
            ph1 = dummy2_vec;
            phh1 = V;
            maske1 = Bit_mask(max - WordSize);
          }
        }
        while ((max >= 0) && (!pivot)) {
          if (!(dummy1_vec->get_row(i) & Bit_mask(max)))
            max --;
          else {
            ph1 = dummy1_vec;
            phh1 = X;
            maske1 = Bit_mask(max);
            pivot = true; // Pivotelement in Zeile max
          }
        }
        if (max != j) {
          for (k = 0; k < M; k++)
            {
              if (dummy2_vec->get_row(k) & maske2) {
                if (!(ph1->get_row(k) & maske1)) {
                  ph1->put_row(k, ph1->get_row(k) | maske1);
                  dummy2_vec->put_row(k, dummy2_vec->
                                      get_row(k)^ maske2);
                }
              }
              else
                if (ph1->get_row(k) & maske1) {
                  ph1->put_row(k, ph1->get_row(k) ^ maske1);
                  dummy2_vec->put_row(k, dummy2_vec->
                                      get_row(k)|maske2);
                }
            }
          for (k = 0; k < M; k++) {
            if (V->get_row(k) & maske2) {
              if (!(phh1->get_row(k) & maske1)) {
                phh1->put_row (k, phh1->get_row(k) | maske1);
                V->put_row(k, V->get_row(k) ^ maske2);
              }
            }
            else
              if (phh1->get_row(k) & maske1) {
                phh1->put_row (k, phh1->get_row(k) ^ maske1);
                V->put_row(k, V->get_row(k) | maske2);
              }
          }
        }
        muster2 ^= maske2;

        if (dummy2_vec->get_row(i) & muster2)
          while (max >= WordSize) {
            if (dummy2_vec->get_row(i) &
                Bit_mask(max - WordSize)) {
              for (k = i; k >= 0; k--)
                if (dummy2_vec->get_row(k) &
                    Bit_mask(j - WordSize))
                  dummy2_vec->bitwise_xor(k, Bit_mask(max-WordSize));
              
              for (k = 0; k < M; k++)
                if (V->get_row(k) & Bit_mask(j - WordSize))
                  V->bitwise_xor(k, Bit_mask(max - WordSize));
            }
            max--;
          }
        if (dummy1_vec->get_row(i) != 0) {
          max = WordSize -1;
          while (max >= 0) {
            if (dummy1_vec->get_row(i) & Bit_mask(max)) {
              for (k = i; k >= 0; k--)
                if (dummy2_vec->get_row(k)
                    & Bit_mask(j - WordSize))
                  dummy1_vec->bitwise_xor(k, Bit_mask(max));
              
              for (k = 0; k < M; k++)
                if (V->get_row(k) & Bit_mask(j - WordSize))
                  X->bitwise_xor(k, Bit_mask(max));
            }
            max--;
          }
        }
      }
      i--;
      j--;
    }
    
    while ((j >= 0) && (!fertig)) {
      maske2 = Bit_mask(j);
      pivot = false;

      while ((i >= 0) &&
             (!(dummy1_vec->get_row(i) & muster1)) &&
             (!(dummy2_vec->get_row(i) & muster2)))
        i--;
      if (i == -1) {
        fertig = true;
        rank = j;
      }
      else {
        max = j;
        while ((max >= 0) && (!pivot)) {
          if (!(dummy1_vec->get_row(i) & Bit_mask(max)))
            max --;
          else {
            ph1 = dummy1_vec;
            phh1 = X;
            maske1 = Bit_mask(max);
            pivot = true; // Pivotelement in row max
          }
        }
        if (max != j) {
          for (k = 0; k < M; k++) {
            if (dummy1_vec->get_row(k) & maske2) {
              dummy1_vec->bitwise_or(k, maske1);
              dummy1_vec->bitwise_and(k, maske2);
            }
            else
              if (dummy1_vec->get_row(k) & maske1)
                dummy1_vec->bitwise_xor(k, maske1);
          }
          for (k = 0; k < M; k++) {
            if (X->get_row(k) & maske2) {
              X->bitwise_or(k, maske1);
              X->bitwise_xor(k, maske2);
            }
            else
              if (X->get_row(k) & maske1)
                X->bitwise_xor(k, maske1);
          }
        }
        muster1 ^= maske2;

        if (dummy1_vec->get_row(i) & muster1) {
          while (max >= 0) {
            if (dummy1_vec->get_row(i) & Bit_mask(max)) {
              for (k = 0; k < M; k++) {
                if (X->get_row(k) & Bit_mask(j)) {
                  if (X->get_row(k) & Bit_mask(max))
                    X->bitwise_xor(k, Bit_mask(max));
                  else
                    X->bitwise_or(k, Bit_mask(max));
                }
              }
              for (k = i; k >= 0; k--) {
                if (dummy2_vec->get_row(k) & Bit_mask(j)) {
                  if (dummy2_vec->get_row(k) & Bit_mask(max))
                    dummy2_vec->bitwise_xor(k, Bit_mask(max));
                  else
                    dummy2_vec->bitwise_or(k, Bit_mask(max));
                }
              }
            }
            max--;
          }
        }
      }
      i--;
      j--;
    }
  }
  else {
    rank = WordSize;
  }
  result_rank = rank;
}



void
lanczos::compute_D(lanczos_small_matrix& VAAV,
                   const lanczos_small_matrix& VAV,
                   const lanczos_small_matrix& W_inv,
                   const value_type S,
                   lanczos_small_matrix& eta,
                   lanczos_small_matrix& result) const
{
        lanczos_small_matrix dummy1_mat;
        //
        // computing D_(i+1)
        // Result will be stored in result
        // VAAV will be computed into VAAV SS^T
        // computing eta: = VAAVSS^T + VAV
        //

        //
        // Computing  V_0^T * A^2 * V_0 * S_0 S_0^T
        //

        if (S != One_mask)
                VAAV.eliminate(S);
        //
        // eta = V_0^T  A^2  V_0 * S_0 S_0^T + V_0^T A V_0
        //

        for (size_type i = 0; i < WordSize; i++)
                eta.put_row(i, VAAV.get_row(i) ^ VAV.get_row(i));

        //
        // dummy1_mat = winv * eta
        //
        W_inv.mult_right_to(eta, dummy1_mat);

        //
        // subtraction of I_N
        //
        for (size_type i = 0; i < WordSize; i++)
                result.put_row(i, Bit_mask(i) ^ dummy1_mat.get_row(i));
}



void
lanczos::compute_E(const lanczos_small_matrix& VAV,
                   const lanczos_small_matrix& W_inv_prev,
                   const value_type S,
                   lanczos_small_matrix& result) const
{
        //
        // computing E_(i+1)
        // Result will be stored in result
        // VAV will be computed into VAV SS^T
        //

        W_inv_prev.mult_right_to(VAV, result);

        if (S != One_mask) result.eliminate(S);
}



void
lanczos::compute_F(const lanczos_small_matrix& VAV_prev,
                   const lanczos_small_matrix& eta_prev,
                   const lanczos_small_matrix& W_inv_prev,
                   const lanczos_small_matrix& W_inv_pprev,
                   const value_type S,
                   lanczos_small_matrix& result) const
{
        lanczos_small_matrix dummy1_mat, dummy2_mat;
        //
        // computing F_(i+i)
        //

        VAV_prev.mult_right_to(W_inv_prev, dummy1_mat);

        //
        // dummy2 = I_N - dummy1
        //
        for (size_type i = 0; i < WordSize; i++)
                dummy2_mat.put_row(i, Bit_mask(i) ^
                                    dummy1_mat.get_row(i));

        //
        // dummy1 = W_(i-2)^(inv) dummy2
        //

        W_inv_pprev.mult_right_to(dummy2_mat, dummy1_mat);

        //
        //  delta_prev * S_i^T S_i
        //
        for (size_type i = 0; i < WordSize; i++)
                dummy2_mat.put_row(i, eta_prev.get_row(i) & S);

        //
        // F_(i+1) = dummy1 * dummy2_mat
        //

        dummy1_mat.mult_right_to(dummy2_mat, result);
}



void
lanczos::solve ()
{
  size_type counter, wur = 0, bound, i, comp_F;
  value_type S_entry = One_mask;
  bool end = false;
  lanczos_small_matrix *VAV, *VAAV, *VAV_prev,
    *eta_next, *eta, *eta_prev, *eta_pprev, *w_inv, *w_inv_prev,
    *w_inv_pprev, *delta, *delta_prev, *alpha, *beta, *gamma;
  
  lanczos_vector_block *J;
  //unsigned long  **W;
  

#if SIZEOF_LONG == 4
  bound = A.number_of_columns()  / (WordSize >> 1) + 10;
#elif SIZEOF_LONG == 8
  bound = A.number_of_columns()  / (WordSize >> 3) + 10;
#endif
   
  
  //**********************************************************************
  // Lege Struktur fuer kleine Matrizen an; rows = columns = 0;
  // Allokiere Speicherplatz (mit 0 initialisiert ) fuer 32 X 32 BIT)
  //**********************************************************************
  
  lanczos_small_matrix dummy1_mat;  
  
  VAV_prev = new lanczos_small_matrix();
  eta = new lanczos_small_matrix();
  eta_next = new lanczos_small_matrix();
  eta_pprev = new lanczos_small_matrix();
  eta_prev = new lanczos_small_matrix();
  alpha = new lanczos_small_matrix();
  beta = new lanczos_small_matrix();
  gamma = new lanczos_small_matrix();
  delta = new lanczos_small_matrix();
  w_inv = new lanczos_small_matrix();
  w_inv_prev = new lanczos_small_matrix();
  w_inv_pprev = new lanczos_small_matrix();
  delta_prev = new lanczos_small_matrix();
  VAV = new lanczos_small_matrix();
  VAAV = new lanczos_small_matrix();
  
  
  V_pprev = new lanczos_vector_block(A.number_of_columns());
  J = new lanczos_vector_block(A.number_of_columns());
  
  //*********************
  // V_0 = A * X
  //*********************
  
  A.mult_vectorblock_to(*X, *V_pprev);
  
  //*********************************
  // eta_pprev = V(_0)^T * V(_0)
  //*********************************
  
  V_pprev->mult_to(*V_pprev, *eta_pprev);
  
  //**********************************************************
  // J = A * V_i  (Einzige Anwendung der Matrix A)
  //**********************************************************
  A.mult_vectorblock_to(*V_pprev, *J);
  
  //**********************************
  // VAV = V_0^T * A * V_0
  //**********************************
  
  V_pprev->mult_to(*J, *VAV);
  
  // Testen, ob VAV Nullvektor
  
  if (VAV->is_zero()) {
    V = V_pprev;
    end = true;
  }
  
  if (!end) {
    //*********************************************************
    // VAAV = (A *V_0)^T * A * V_0 = V_0^T * A^2 * V_0
    //*********************************************************
    J->mult_to(*J, *VAAV);
    //***********************************
    // Berechnung von S_0,   W_0^inv
    //***********************************
    
    i = compute_tau2(*VAV, S_entry, *w_inv_pprev);
    wur += i;
    
    //**********************************************************************
    // Berechnung des aktuellen Teilsummanden der Partialsumme der Loesung;
    // X = \sum _{i=0}^{counter}V_i W_i^inv V_i^T V_0
    //***********************************************************************
    // X beinhaltet Y, so wird schon fUer spaeter X-Y berechnet
    w_inv_pprev->mult_right_to(*eta_pprev, dummy1_mat);
    
    V_pprev->mult_small_to(dummy1_mat, *dummy1_vec);
    X->add(*dummy1_vec);
    
    //
    // computing D_(1)
    //
    
    compute_D(*VAAV, *VAV, *w_inv_pprev, S_entry, *delta, *alpha);
    
    
    //********************************************************************
    // Bererchnung von V_(1)
    //********************************************************************

    for (i = 0; i < M; i++)
      V_prev->put_row(i, J->get_row(i) & S_entry);
    
    V_pprev->mult_small_to(*alpha, *dummy1_vec);

    for (i = 0; i < M; i++)
      V_prev->put_row(i, V_prev->get_row(i) ^ dummy1_vec->get_row(i));

    //*******************************
    // Berechnung von V_{1}^T V_0
    //*******************************
    V_prev->mult_to(*V_pprev, *eta_prev);
    
    // naechste Iteration
    
    //**********************************************************
    // J = A * V_1  (Einzige Anwendung der Matrix A)
    //**********************************************************
    A.mult_vectorblock_to(*V_prev, *J);
    
    
    //**********************************
    // VAV = V_1^T * A * V_1
    //**********************************
    V_prev->mult_to(*J, *VAV_prev);
    
    //*********************************************************
    // VAAV = (A *V_1)^T * A * V_1 = V_1^T * A^2 * V_1
    //*********************************************************
    J->mult_to(*J, *VAAV);
    
    if (VAV_prev->is_zero()) {
      V = V_prev;
      end = true;
    }
    if (!end) {
      //***********************************
      // Berechnung von S_1,   W_1^inv
      //***********************************
      //printf( "Iteration 1:\n");
      i = compute_tau2(*VAV_prev, S_entry, *w_inv_prev);
      //bin_out(S_entry);
      wur += i;
      //printf( "Dimension W_1 : %d, Groesse des Unterraums : %d\n", i, wur );
      
      
      //*************************************************************************
      // Berechnung des aktuellen Teilsummanden der Partialsumme der Loesung;
      // X = \sum _{i=0}^{counter}V_i W_i^inv V_i^T V_0
      //*************************************************************************
      
      w_inv_prev->mult_right_to(*eta_prev, dummy1_mat);
      V_prev->mult_small_to(dummy1_mat, *dummy1_vec);
      
      X->add(*dummy1_vec);
      
      //
      // computing D_(2)
      // computing E_(2)
      //
      
      compute_D(*VAAV, *VAV_prev, *w_inv_prev, S_entry, *delta_prev, *alpha);
      compute_E(*VAV_prev, *w_inv_pprev, S_entry, *beta);
      
      
      //******************************************************************
      // Bererchnung von V_2
      //******************************************************************
      
      for (i = 0; i < M; i++)
        V->put_row(i, J->get_row(i) & S_entry);
      V_prev->mult_small_to(*alpha, *dummy1_vec);
      for (i = 0; i < M; i++)
        V->put_row(i, V->get_row(i) ^ dummy1_vec->get_row(i));
      V_pprev->mult_small_to(*beta, *dummy1_vec);
      for (i = 0; i < M; i++)
        V->put_row(i, V->get_row(i) ^ dummy1_vec->get_row(i));
      
      //************************************************************
      //Berechnung von V_{2}^T V_0 durch Operationen mit kleinen Matrizen//
      //************************************************************
      V->mult_to(*V_pprev, *eta);
      counter = 2;
      
      
      while (!end) 
        {
          //**********************************************************
          // J = A * V_i  (Einzige Anwendung der Matrix A)
          //**********************************************************
          
          A.mult_vectorblock_to(*V, *J);
        
          
         //**********************************
         // VAV = V_i^T * A * V_i
         //**********************************
        
        V->mult_to(*J, *VAV);
                                //*********************************************************
                                // VAAV = (A *V_i)^T * A * V_i = V_i^T * A^2 * V_i
                                //*********************************************************
        
        J->mult_to(*J, *VAAV);
                                // vergl auf Nullvektor
        end = VAV->is_zero();
        if (!end) {
          //***********************************
          // Berechnung von S_i,   W_i^inv
          //***********************************
          comp_F = 1;
          //if (S_entry == One_mask)
          //compute_F = 0;
          // warum sollte fuer S_(i-1) gamme = 0 sein?!
          
          i = compute_tau2(*VAV, S_entry, *w_inv);
          //bin_out(S_entry);
          wur += i;
          
          //**********************************************************
          // Berechnung des aktuellen
          // Teilsummanden der Partialsumme der Loesung;
          // X = \sum _{i=0}^{counter}V_i W_i^inv V_i^T V_0
          //**********************************************************
          
          w_inv->mult_right_to(*eta, dummy1_mat);
          V->mult_small_to(dummy1_mat, *dummy1_vec);
          
          X->add(*dummy1_vec);
          
          //
          // computing D_(i+1)
          // computing E_(i+1)
          // computing F_(i+1)
          //
          
          compute_D(*VAAV, *VAV, *w_inv, S_entry, *delta, *alpha);
          compute_E(*VAV, *w_inv_prev, S_entry, *beta);
          compute_F(*VAV_prev, *delta_prev, *w_inv_prev, *w_inv_pprev,
                    S_entry, *gamma);
          
          //
          // computing V_(i+1)
          //
          
          if (S_entry != One_mask)
            J->eliminate(S_entry);
          V->mult_small_to(*alpha, *dummy1_vec);
          J->add(*dummy1_vec);
          V_prev->mult_small_to(*beta, *dummy1_vec);
          J->add(*dummy1_vec);
          
          V_pprev->mult_small_to(*gamma, *dummy1_vec);
          J->add(*dummy1_vec);
          
          //
          // computing V_{i+1}^T V_0
          //
          
          alpha->mult_right_transpose_to(*eta, *eta_next);
          beta->mult_right_transpose_to(*eta_prev, dummy1_mat);
          eta_next->add(dummy1_mat);
          gamma->mult_right_transpose_to(*eta_pprev, dummy1_mat);
          eta_next->add(dummy1_mat);
          
          counter++;
          if (counter > bound) {
            lidia_error_handler("lanczos", "Internal error: counter > bound");
          }
          else {
            lanczos_vector_block *hilf_vect;
            lanczos_small_matrix *hilf_matrix;
            
            hilf_vect = V_pprev;
            V_pprev = V_prev;
            V_prev = V;
            V = J;
            J = hilf_vect;
            
            hilf_matrix = eta_pprev;
            eta_pprev = eta_prev;
            eta_prev = eta;
            eta = eta_next;
            eta_next = hilf_matrix;

            hilf_matrix = VAV_prev;
            VAV_prev = VAV;
            VAV = hilf_matrix;

            hilf_matrix = delta_prev;
            delta_prev = delta;
            delta = hilf_matrix;

            hilf_matrix = w_inv_pprev;
            w_inv_pprev = w_inv_prev;
            w_inv_prev = w_inv;
            w_inv = hilf_matrix;
          }
        }
      }
    }
  }

  compute_result();
  post_process();
    
#ifdef DEBUG    
  if (X->is_zero()) {
    std::cout<<"\nResult was a NULL-VECTOR"<<std::flush;
  }
#endif
  
  delete(VAV_prev);
  delete(VAV);
  delete(VAAV);
  delete(eta);
  
  delete(eta_next);
  
  delete(eta_pprev);
  delete(eta_prev);
  delete(alpha);
  delete(beta);
  delete(gamma);
  delete(delta);
  delete(w_inv);
  delete(w_inv_prev);
  delete(w_inv_pprev);
  delete(delta_prev);
  delete(V_pprev);
  delete(J);
}



void
lanczos::post_process()
{
        // first of all, I think I sould finde solutions for B*x=0

        lanczos_vector_block test_vec(A.number_of_rows());

        size_type j, h, i;
        value_type nullmask;
        value_type g;

        i = 0;
        j = 0;

        while (i < A.number_of_columns()) {
                j = 0;
                h = A.get_vector(i).get_number_of_entries();
                // Wieviel Eintraeg gibt es in der i-ten Spalte ?
                // solange Spalteneintraege gehe alle Spalten einer Zeile durch
                while (j < h) {

                        g = A.get_vector(i).get_entry(j);
                        test_vec.put_row(g, X->get_row(i) ^ test_vec.get_row(g));
                        j++;
                }
                i++;
        }
        // Suche alle Nullvectoren
        nullmask = 0;

        for (i = 0; i < test_vec.get_length(); i++) {
                for (j = 0; j < WordSize; j++)
                        if (test_vec.get_row(i)&Bit_mask(j))
                                nullmask |= Bit_mask(j);

        }
        if (nullmask) {
                nullmask = nullmask ^ One_mask;
                for (i = 0; i < X->get_length(); i++)
                        X->bitwise_and(i, nullmask);
        }
}



//***********************************************************************
//* Aufruf : bin_out(a)                                                 *
//***********************************************************************
//* Parameter : long a : auszugebende Zahl                              *
//***********************************************************************
//* Rueckgabewert : %                                                   *
//***********************************************************************
//* Aufgabe : binaere Ausgabe einer long - Zahl                         *
//***********************************************************************
//* Literaturverzeichnis : D. E. Knuth : The art of computer            *
//*                                      programming Vol. 1             *
//***********************************************************************

void bin_out(const unsigned long a)
{
    /* for(size_t i=0; i < WordSize; i++) { */
    for(int i=WordSize-1; i >= 0; i--) {  /* so that LSB is on right */
        if(a & Bit_mask(i))
            printf("1 ");
        else
            printf("0 ");
    }
    printf("\n");
    fflush(stdout);
}

#ifdef LIDIA_NAMESPACE
}       // end of namespace LiDIA
#endif

