//==============================================================================================
//
// This file is part of LiDIA --- a library for computational number theory
//
// Copyright (c) 1994--2001 the LiDIA Group.  All rights reserved.
//
// See http://www.informatik.tu-darmstadt.de/TI/LiDIA/
//
//----------------------------------------------------------------------------------------------
//
// $Id: mpqs.cc,v 2.17 2006/03/06 12:08:36 lidiaadm Exp $
//
// Author : Volker Mueller (VM), based on an implementation
//                of Thomas Sosnowski
// Changes : See CVS log
//
// CAVEAT: This MPQS implementation is superseded by
//         single_factor<bigint>::mpqs(). In the near future, it will
//         delegate all calls to single_factor<bigint>.
//         rational_factorization will eventually be removed.
//
//==============================================================================================

// Description : see the diploma thesis of Thomas Sosnowski

#ifdef HAVE_CONFIG_H
#include	"config.h"
#endif
#include	"LiDIA/rational_factorization.h"
#include        "LiDIA/random_generator.h"
#include	"LiDIA/udigit.h"
#include	"LiDIA/lanczos.h"
#include	"LiDIA/lidia_signal.h"
#include	"LiDIA/mpqs_timing.h"
#include        "LiDIA/arith.inl"

#include	<cmath>
#include	<cstdlib>
#include	<cstdio>
#include        <fstream>


#ifdef LIDIA_NAMESPACE
namespace LiDIA
{
#endif



#ifdef DEBUG                    // for testing correctness of fulls, partials
  static bigint global_kN;     // and combined relations
  static int *global_FB;
#endif

  void concat (char *, char *, char *); // functions for file handling.
  void merge_unique (char *, char *, char *); // in future these functions
  char *get_name (const char *);     // should be collected in an own
  void sort_file_lp (FILE *, char *); // class for LP variations
  void sort_file_line (char *, char *);
  
  const unsigned int STLEN = 6000;
  const unsigned int CAND_NUMBER = 700;
  typedef unsigned char SIEBTYP;
  
  
//******************************************************
// our own signal handler
//******************************************************

#define lidia_error_handler_n(f, m)\
{ char *s; \
  s = get_name("RELATIONS"); std::remove(s);  \
  s = get_name("LANCZOS_FORMAT"); std::remove(s);  \
  s = get_name("NEW_PARTIALS"); std::remove(s);  \
  s = get_name("NEW_FULLS"); std::remove(s);  \
  s = get_name("LP_RELATIONS"); std::remove(s);  \
  s = get_name("TEMPORARY"); std::remove(s);  \
  s = get_name("SORT_LINE0"); std::remove(s);  \
  s = get_name("SORT_LINE1"); std::remove(s);  \
  s = get_name("SORT_LP0"); std::remove(s);  \
  s = get_name("SORT_LP1"); std::remove(s);  \
  lidia_error_handler(f, m); }



LIDIA_SIGNAL_FUNCTION (stop_rf_mpqs) 
{
  lidia_error_handler_n ("rational_factorization",
			 "mpqs::the factorization was interrupted");
}

  
  // the parameter table for mpqs has the following format:
  //  dd, approximation accuracy, size sieving interval, size FB,
  //  # prime factors in A, total # primes for determination of A,
  // starting sieving index
  //
  // these values are somewhat experimental and can be improved by
  // heavy testing.
  
  const float rational_factorization::qs_params[66][7] = {
    {15, 0.8, 1200, 35, 3, 10, 2},
    {16, 0.8, 1400, 35, 3, 10, 2},
    {17, 0.8, 3000, 40, 3, 10, 2},
    {18, 0.8, 3000, 60, 3, 10, 2},
    {19, 0.8, 3600, 80, 3, 10, 2},
    {20, 0.8, 4000, 100, 3, 10, 2},
    {21, 0.8, 4250, 100, 3, 10, 2},
    {22, 0.8, 4500, 120, 3, 10, 3},
    {23, 0.8, 4750, 140, 3, 10, 3},
    {24, 0.8, 5000, 160, 3, 12, 4},
    {25, 0.8, 5000, 180, 3, 12, 4},
    {26, 0.9, 6000, 200, 5, 8, 4},
    {27, 1.17, 6000, 220, 5, 8, 5},
    {28, 1.17, 6500, 240, 5, 8, 5},
    {29, 1.17, 6500, 260, 5, 8, 5},
    {30, 1.36, 7000, 325, 5, 8, 5},
    {31, 1.36, 7000, 355, 5, 8, 5},
    {32, 1.36, 7500, 375, 5, 8, 5},
    {33, 1.43, 7500, 400, 6, 8, 6},
    {34, 1.43, 7500, 425, 6, 8, 6},
    {35, 1.43, 7500, 550, 6, 10, 6},
    {36, 1.43, 8000, 650, 6, 10, 6},
    {37, 1.69, 9000, 750, 6, 10, 7},
    {38, 1.69, 10000, 850, 6, 10, 7},
    {39, 1.69, 11000, 950, 6, 10, 7},
    {40, 1.69, 14000, 1000, 6, 10, 7},
    {41, 1.69, 14000, 1150, 6, 10, 8},
    {42, 1.69, 15000, 1300, 6, 10, 8},
    {43, 1.69, 15000, 1600, 6, 10, 8},
    {44, 1.69, 15000, 1900, 7, 10, 9},
    {45, 1.69, 15000, 2200, 7, 10, 9},
    {46, 1.69, 20000, 2500, 7, 10, 9},
    {47, 1.69, 25000, 2500, 7, 10, 10},
    {48, 1.69, 27500, 2700, 7, 10, 10},
    {49, 1.69, 30000, 2800, 7, 10, 10},
    {50, 1.75, 35000, 2900, 7, 11, 10},
    {51, 1.75, 40000, 3000, 7, 11, 10},
    {52, 1.85, 50000, 3200, 7, 11, 11},
    {53, 1.85, 50000, 3500, 7, 11, 11},
    {54, 1.95, 80000, 3800, 7, 12, 11},
    {55, 1.95, 90000, 4100, 7, 12, 11},
    {56, 1.95, 100000, 4400, 7, 12, 11},
    {57, 1.95, 80000, 4700, 8, 12, 12},
    {58, 1.95, 80000, 5000, 8, 12, 12},
    {59, 2.15, 130000, 5500, 8, 12, 12},
    {60, 2.15, 140000, 5800, 8, 12, 12},
    {61, 2.15, 150000, 6100, 8, 13, 13},
    {62, 2.15, 160000, 6400, 8, 13, 13},
    {63, 2.35, 170000, 6700, 8, 13, 13},
    {64, 2.35, 180000, 7000, 8, 13, 13},
    {65, 2.35, 190000, 7300, 8, 13, 13},
    {66, 2.35, 200000, 7600, 8, 13, 13},
    {67, 2.4, 150000, 7900, 8, 13, 13},
    {68, 2.4, 150000, 8200, 8, 14, 13},
    {69, 2.4, 130000, 8600, 8, 14, 13},
    {70, 2.45, 130000, 8800, 8, 14, 13},
    {71, 2.45, 130000, 8800, 9, 14, 13},
    {72, 2.4, 260000, 9400, 9, 14, 13},
    {73, 2.4, 270000, 9700, 9, 14, 13},
    {74, 2.4, 280000, 9000, 9, 14, 13},
    {75, 2.6, 140000, 9000, 9, 14, 13},
    {76, 2.6, 160000, 9400, 9, 14, 13},
    {77, 2.6, 180000, 9600, 9, 14, 13},
    {78, 2.6, 200000, 9800, 9, 14, 13},
    {79, 2.65, 220000, 10000, 9, 14, 13},
    {80, 2.65, 250000, 10500, 9, 14, 13}
  };
  
  
  
  //************************************************************
  // compute the number of occurances of char c in the string s
  //************************************************************
  
inline unsigned int number_c (char *s, char c)
{
  unsigned int l = 0;
  
  while (*s++)
    if (*s == c)
      l++;
  
  return (l);
}

  
  
  //************************************************************
  // transform relations from internal format into format needed
  // in lanczos algorithm.
  //************************************************************


static void transform_relations ()
{
  char *s1;
  long val;
  
  s1 = get_name ("RELATIONS");   // file with relations in internal format
  FILE *fpin = fopen (s1, "r");
  
  s1 = get_name ("LANCZOS_FORMAT");   // file with relations in lanczos format
  FILE *fpout = fopen (s1, "w");
  
  char *p, *pp, zeile[STLEN], string[STLEN];
  int i = 0, j, number_entries;
  
  if (!fpin || !fpout) 
    lidia_error_handler_n ("rational_factorization",
			   "mpqs::transform_relations::can't open files");
  
  while (fgets (zeile, STLEN, fpin)) 
    {
      if ((pp = (char *) strchr (zeile, ':'))) 
	{
	  pp += 2;
	  number_entries = number_c (pp, ' ') + 2; // count # entries per row
	  string[0] = '\0';
	  
	  p = (char *) strtok (pp, " \n");
	  
	  while (p != NULL) 
	    {
	      j = std::strtol (p,NULL,10);    // handle only odd exponents
	      if (j & 1) 
		{
                  if(strlen(string) > STLEN - 3)
                    {
                      lidia_error_handler_n("rational_factorization",
                                            "mpqs::transform_relations::insufficent string length");
                    }
                  strcat(string, " 1");
		  p = (char *) strtok (NULL, " \n");
                  val = std::strtol (p,NULL,10);
                  if(val == LONG_MAX)
                  {
                    lidia_error_handler_n ("rational_factorization",
			   "mpqs::transform_relations::error in file content");
                    
                  }
                  else
                  {
                    char valbuffer[STLEN];
                    int vallength = snprintf(valbuffer, STLEN, " %ld", val);
                    if(STLEN <= strlen(string) + vallength)
                    {
                      lidia_error_handler_n("rational_factorization",
                                            "mpqs::transform_relations::insufficent string length");
                    }
                    strcat(string, valbuffer);
                  }
		  p = (char *) strtok (NULL, " \n");
		} 
	      else 
		{
		  number_entries -= 2;
		  p = (char *) strtok (NULL, " \n");
		  p = (char *) strtok (NULL, " \n");
		}
	    }
	  fprintf (fpout, "%d %d 0%s 0\n", ++i, number_entries, string);
	}
    }
  fflush(fpout);
  fclose (fpin);
  fclose (fpout);
}





//********************************************************************
// initialize mpqs parameters from precomputed table
//********************************************************************

void rational_factorization::
qs_read_par (unsigned int stellen, double &T, unsigned int &M,
	     unsigned int &size_FB, unsigned int &P_ONCE,
	     unsigned int &POLY, unsigned int &P_TOTAL,
	     unsigned int &smallstart)
{
  
  unsigned int i;
  
  if (stellen < 15)
    i = 0;
  else
    i = stellen - 15;
  
  T = qs_params[i][1];
  M = static_cast < unsigned int >(qs_params[i][2]); // length of sieve array is 2*M
  
  size_FB = static_cast < unsigned int >(qs_params[i][3]); // size of factor base
  
  P_ONCE = static_cast < unsigned int >(qs_params[i][4]); // number of prime factors of
  
  // coefficient A
  
  P_TOTAL = static_cast < unsigned int >(qs_params[i][5]); // total number of primes,
  
  // the P_ONCEs are taken from
  
  POLY = static_cast < unsigned int >(std::pow (2.0, static_cast < double >(P_ONCE - 1))); // number of coefficients B that
  
  // are used for
  smallstart = static_cast < unsigned int >(qs_params[i][6]); // begin sieving
  
  // with FB[smallstart]
  
}



//************************************************************
// assume that s is a relation of the form "... : FB_i exp_i ..."
// and exp is a vector of exponents.
// the function scans s and adds (sign == 1) or subtracts
// (sign == -1) s to the array exp
//
// -->used in combination of LP relations
//************************************************************

inline void add_sub_exp (char *s, int *exp, int sign)
{
  char *pp, *p, hilf[STLEN];
  long e;
  long val;
  
  strcpy (hilf, s);
  
  if ((pp = strchr (hilf, ':'))) {
    pp += 2;
    
    p = strtok (pp, " \n");
    while (p != NULL) {
      e = std::strtol (p,NULL,10);
      if (e==LONG_MAX) {
        lidia_error_handler_n ("rational_factorization",
                  "mpqs::add_sub_exp::error in input file (invalid exponent)");
      }
      else if(!e) {
        break;
      }
      else
      {
        p = strtok (NULL, " \n");
        val = std::strtol (p,NULL,10);
        if (val == LONG_MAX || val < 0) {
          lidia_error_handler_n ("rational_factorization",
                   "mpqs::add_sub_exp::error in input file (invalid prime)");
        }
        else {
          exp[val] += sign * e;
        }
        p = strtok (NULL, " \n");
      }
    }
  }
}



//******************************************************************
// count full relations which can be combined from partial relations
//******************************************************************

static unsigned long count_partials ()
{
  char *s1, *s2, *s3;
  FILE *fp;
  char line[STLEN];
  
  s1 = get_name ("NEW_PARTIALS");   // new found LP relations
  s2 = get_name ("TEMPORARY");  // temporary file
  s3 = get_name ("LP_RELATIONS");   // holds already processed LP relations
  
  if ((fp = fopen (s1, "r")) != NULL)  // if there are new LP relations
    {
      sort_file_lp (fp, s2); // sort file with partial relations with key LP
      fclose (fp);
      std::remove (s1);
      merge_unique (s3, s2, s1); // merge file to already processed LP
      std::remove (s3);
      rename (s1, s3);
    }
  
  fp = fopen (s3, "r");

  if (fp == NULL)  // no file --> no LP relations can be combined
    return 0; 

  // now the counting phase starts

  long counter = 0, large, large_old = -1, i = 1;
  
  while (fgets (line, STLEN, fp)) {
    // idea: count all LP relations with same LP
    large = std::strtol (line,NULL,10);
    
    if (large == large_old) 
      {
	counter += i;
	i++;
      } 
    else 
      {
	large_old = large;
	i = 1;
      }
  }
  fclose (fp);
  return (counter);
}



//**************************************************************
// combine LP relations to full relations
//**************************************************************

static void zusam (unsigned int size_exp, int * FB, const bigint & kN)
{
  char *s1, *s2, *s3;
  
  s1 = get_name ("LP_RELATIONS");   // file with LP relations
  s2 = get_name ("TEMPORARY");  // temporary file
  sort_file_line (s1, s2);

  FILE *fi = fopen (s2, "r");

  if (fi == NULL)
    return;
  
  FILE *fneu;
  char line[50][STLEN], faktor[STLEN], Hstring[STLEN], Hstring2[STLEN];
  bigint h, h_old, h_inv, neu;
  unsigned int large=0, large_old=0, i;
  int *exp = NULL;
  long pos, ii, j, ss = 0;
  bigint H1, H2, H3;
  
  s2 = get_name ("NEW_FULLS");   // file for new full relations
  fneu = fopen (s2, "w");
  if (fneu == NULL) 
    lidia_error_handler_n ("rational_factorization",
			   "zusam::can't open output file");
  
  exp = new int[size_exp];
  memset(exp,0,size_exp * sizeof(int));
  memset(line,0,sizeof(line));

  while (fgets (line[ss & 1], STLEN, fi)) 
    {
      sscanf (line[ss & 1], "%d @ %*s ", &large);
      ss++;
      
      if (large == large_old) // two LP relations with same large prime
	{                      // -->combination is done
	  pos = 2;
	  while (large == large_old && fgets (line[pos], STLEN, fi))
            {
            if (pos >= 50)
              {
              lidia_error_handler_n ("rational_factorization",
                                     "zusam::too many entries with same large");
              break;
              }
	    sscanf (line[pos], "%d @ %*s ", &large);
	    pos++;
	  }
	  pos--;
	  
	  for (ii = 0; ii < pos - 1; ii++)   
	    {
              Hstring2[0] = '\0';
	      sscanf (line[ii], "%d @ %s ", &large_old, Hstring2);
	      
	      for (j = ii + 1; j < pos; j++) 
		{
                  Hstring[0] = '\0';
		  sscanf (line[j], "%d @ %s ", &large, Hstring);
		  string_to_bigint(Hstring, H1);
		  string_to_bigint(Hstring2, H2);
		  H1 = (H1 * H1) % kN;
		  H2 = (H2 * H2) % kN;
		  
		  memset (exp, 0, size_exp * sizeof (int));
		  add_sub_exp (line[j], exp, 1);
		  add_sub_exp (line[ii], exp, -1);
		  memset (faktor, 0, STLEN);
		  
		  if (exp[1] != 0)
		    {
		      sprintf (faktor, "%s 1 1", faktor);
		      H2.negate();
		    }
		  
		  for (i = 2; i < size_exp; i++)
		    {
		      if (exp[i] != 0)
			{
			  sprintf (faktor, "%s %d %d", faktor, exp[i], i);
			  
			  if (exp[i] > 0)
			    {
			      power (H3, bigint (FB[i]), exp[i]);
			      H2 = (H2 * H3) % kN;
			    }
			  else
			    {
			      power (H3, bigint (FB[i]), -exp[i]);
			      H1 = (H1 * H3) % kN;
			    }
			}
		      
		      if (faktor[0] == (char) NULL)
			continue;
		    }
		  
		  if (H1 != H2 && H1 != kN - H2)
		    fprintf (fneu, "%s %s : %s 0\n", Hstring, Hstring2, 
			     faktor);
		
#ifdef DEBUG 
		  if ((H1 - H2) % global_kN != 0)
		    std::cout << "\n ERROR found: " << std::flush;
#endif
	      }
	    }
	  strcpy (line[(ss + 1) & 1], line[pos]);
	  sscanf (line[pos], "%d @ %*s ", &large_old);
	} 
      else
	large_old = large;
    }
  
  fclose (fi);
  fclose (fneu);  
  delete[]exp;

  s3 = get_name ("TEMPORARY");
  std::remove (s3);
  s1 = get_name ("RELATIONS");
  
  concat (s1, s2, s3);      // concat fulls and newly combined LP relations
  std::remove (s1);
  std::remove (s2);
  sort_file_line (s3, s1);  // automatically delete multiple relations
}



//**********************************************************
// determine approximate running time for mpqs on i digit
// number, if print = true, then print result on stdout
// seems to be very out-of-date !!
//**********************************************************

double rational_factorization::zeitqs (unsigned int dec_size, bool print)
{
  double sec;
  int day, st, min;
  double timings_celeron[] = {3.2, 3.4, 4.2, 4.6, 5.5, 8.0, 9.3, 11.5, 14.1,
  17.4, 17.8, 23.9, 29.5, 40.9, 61.7, 109.5, 142.6, 199.6, 264.2}; 

  // average timings for integers >= 39, <= 57 digits
  // on Celeron 

  if (dec_size <= 38)        // "const time" for small numbers
    sec = 2 * mpqs_machine_factor;  // mpqs_machine_factor is computed externally
  else 
    {
      if (dec_size <= 57)
	sec = timings_celeron[dec_size - 39] * mpqs_machine_factor;
      else
	sec = std::pow (2.0, LiDIA::log2(timings_celeron[57-39]) 
			+ (dec_size - 57)*0.53) * mpqs_machine_factor;
    }
  
  if (print == false)
    return sec;
  
  st = static_cast < int >(sec / 3600.0);
  min = static_cast < int >(sec / 60.0);
  
  if (st < 1) 
    {
      if (min < 1)
	std::cout << static_cast < int >(sec) << " seconds \n";
    else
      std::cout << min << " minutes " << 
	static_cast < int > (sec - min * 60.0) <<" seconds\n";
    } 
  else 
    {
      day = static_cast < int >(st / 24.0);
      min -= 60 * st;

      if (day < 1)
	std::cout << st << " hours " << min << " minutes \n";
      else 
	{
	  st -= 24 * day;
	  std::cout << day << " days " << st <<" hours ";
	  std::cout << min << " minutes";
	}
    }
  return (sec);
}



//*******************************************************************
// main sieving routine:
//
// FB is a pointer to an array which holds the factor basis
// LOGP is a pointer to an array which holds the approximations for
// the logarithms of the factor basis elements
// START1, START2 are arrays for starting points for different FB elements
// sieb points to a sieve array
// ende points to the end of the sieve array
// M is the size of the sieving interval
// CANDIDATE is an array  which is filled with candidates which might split
// smallstart marks the first FB element which is used for sieving
//
//******************************************************************

// FIXME: not portable!
#ifdef WORDS_BIGENDIAN
// big endian
# define INT32_OCTET_0_0x80	0x80000000U
# define INT32_OCTET_1_0x80	0x00800000U
# define INT32_OCTET_2_0x80	0x00008000U
# define INT32_OCTET_3_0x80	0x00000080U
#else
// little endian
# define INT32_OCTET_0_0x80	0x00000080U
# define INT32_OCTET_1_0x80	0x00008000U
# define INT32_OCTET_2_0x80	0x00800000U
# define INT32_OCTET_3_0x80	0x80000000U
#endif

static void
qs_sieve_interval (int *FB, SIEBTYP * LOGP, int *START1,
		   int *START2, SIEBTYP * sieb, SIEBTYP * ende,
		   unsigned int M, int *CANDIDATE, unsigned int smallstart)
{
  register int p, l, *fbp, *lsieb = (int *) sieb; // FIXME: SIEBTYP and int are incompatible types
  register SIEBTYP logp;
  register SIEBTYP *begin;
  
  register int x, counter = 0, M_2 = M << 1;
  register int oldstart1;
  
  memset (sieb, 0, (M_2) * sizeof (SIEBTYP));
  
  fbp = &FB[smallstart];
  l = smallstart;
  
  while ((p = *fbp++) != 0) {
    logp = LOGP[l];
    begin = sieb + START1[l];
    oldstart1 = START1[l];
    
    for (;;) {
      // sieving with FB[l] from START1[l]
      if (begin <= ende) {
	(*begin) += logp;
	begin += p;
      } else
	break;
    }
    
    if (oldstart1 != START2[l]) {
      begin = sieb + START2[l];
      for (;;) {
	if (begin <= ende) {
	  (*begin) += logp;
	  begin += p;
	} else
	  break;
      }
    }
    l++;
  }
  
  l = 0;
  while (l < M_2) {
    if ((p = (*lsieb)) & (0x80808080)) {
      // check whether at least one of
      // four consequent sieve entries
      // is a candidate
      
      x = l;
      
      if (p & (INT32_OCTET_0_0x80 | INT32_OCTET_1_0x80)) {
	if (p & (INT32_OCTET_0_0x80)) {
	  CANDIDATE[counter++] = x;
	  if (p & (INT32_OCTET_1_0x80))
	    CANDIDATE[counter++] = x + 1;
	} else
	  CANDIDATE[counter++] = x + 1;

	if (p & (INT32_OCTET_2_0x80))
	  CANDIDATE[counter++] = x + 2;
	
	if (p & (INT32_OCTET_3_0x80))
	  CANDIDATE[counter++] = x + 3;
      } else {
	if (p & (INT32_OCTET_2_0x80))
	  CANDIDATE[counter++] = x + 2;
	else
	  CANDIDATE[counter++] = x + 3;
      }
    }
    lsieb++;
    l += 4;
  }
  CANDIDATE[counter] = 0;
}



//******************************************************************
// compute optimal multiplier in the set cand of candidates and return
// this value
//******************************************************************

int rational_factorization::
compute_multiplier (const bigint & N, int bis, ecm_primes & prim)
{
  int cand[5] = {1, 3, 5, 7, 11};
  
  bigint kN;
  register unsigned long plauf;
  register int p, j, i, k = 1, nmod4;
  double wert, bestwert = 1, plus;
  
  nmod4 = static_cast < int >(N.least_significant_digit ()) & 0x3;
  
  for (j = 0; j <= 4; j++) 
    {
      if ((((p = cand[j]) * nmod4) & 0x00000003) != 1)
	continue;
      
      wert = -0.7 * LiDIA::log2 (static_cast < double >(p));
      
      multiply (kN, N, p);
      
      if ((kN.least_significant_digit () & 0x7) == 1)
	wert += 1.38629;
      
      plauf = prim.getprimes ();
      i = 0;
      
      while (i <= bis) 
	{
	  if (legendre
	      (static_cast <
	       int >(remainder (kN, static_cast < long >(plauf))),
	       static_cast < int >(plauf)) == 1)
	    {
	      i++;
	      plus =
		LiDIA::log2 (static_cast < double >(plauf)) / static_cast <
		double >(plauf);
	      if ((p % plauf) == 0)
		wert += plus;
	      else
		wert += 2 * plus;
	    }
	  
	  plauf = prim.getprimes ();
	}
      
      if (wert > bestwert) {
	bestwert = wert;
	k = p;
      }
      prim.resetprimes (1);
    }
  return (k);
}



//******************************************************************
// create the factor basis
//******************************************************************


int rational_factorization::
create_FB (unsigned int size, const bigint & kN, int **FB,
	   ecm_primes & prim)
{
  register unsigned int osize, p;
  register int *fbb;
  
  if (!(*FB = new int[size + 3])) 
    lidia_error_handler_n ("rational_factorization",
			   "mpqs::create_FB::can't allocate memory");
  
  fbb = FB[0];
  
  *fbb++ = static_cast < int >(size);
  
  *fbb++ = -1;
  
  osize = 0;
  prim.resetprimes (1);
  p = 2;
  
  while (osize < size) 
    {
      if (legendre
	  (static_cast < int >(remainder (kN, static_cast < long >(p))),
	   p) != -1)
      {
	// p in factor base!
	*fbb++ = p;
	osize++;
      }
      p = prim.getprimes ();
    }
  
  *fbb = 0;
  return (0);
}



//**************************************************************
// determine the number of ones in binary representation of bi
//**************************************************************

inline int count_ones (unsigned int bi)
{
  int s = 0;
  
  while (bi > 0) 
    {
      if (bi & 1)
	s++;
      bi >>= 1;
    }
  return s;
}



//*****************************************************************
// compute coefficients of sieving polynomial for self initializing
// variant. Coefficients A and B are returned abd several tables are
// updated -->see Thomas Sosnowskis diploma thesis.
//*****************************************************************

static void
compute_coeff (bigint & A, bigint & B, const bigint & kN, int *FB,
	       int *SQRTkN, int *START1, int *START2,
	       int P_ONCE, int P_TOTAL, int SIEBSTART,
	       int **vorb, int *Q_prime, int *Q_prime_glob,
	       bigint * BG, lidia_size_t index_i,
	       int start_fb, int *a_inv, bigint & A4_inverse,
	       unsigned int &bin_index)
{
  register int p, size_FB;
  int SIEBS, tmp, tmp1, tmp2;
  lidia_size_t j, nu_2, i;
  bigint reserve, TMP;
  
  if (index_i == 0) 
    {
      bin_index++;
      
      while (count_ones (bin_index) != P_ONCE)
	bin_index++;
      
      i = 0;
      for (j = 0; j < P_TOTAL; j++) // determine primes used for A
	{                      // in this iteration
	  if (bin_index & (1 << j)) {
	    Q_prime[i] = Q_prime_glob[j];
	    i++;
	  }
	}
      
      A.assign (Q_prime[0]); // compute coefficient A
      for (i = 1; i < P_ONCE; i++)
	multiply (A, A, Q_prime[i]);
      
      shift_left (A4_inverse, A, 2);
      
      // compute BG[0] to BG[P_ONCE-1]
      // each B is of the form
      // v_0*BG[0]+v_1*BG[1]+...+v_(P_ONCE-1)*BG[(P_ONCE-1)],
      // where each v_j is +1 or -1
      
      // this has to be done only once (for index_i==0) for the
      // coefficient A; if index_i > 0 then there is a linear
      // recursion for B
      
      for (i = 0; i < P_ONCE; i++) {
	p = Q_prime[i];
	divide (reserve, A, p);
	tmp =
	  static_cast < int >(remainder (reserve, static_cast < long >(p)));
	if (tmp < 0)
	  tmp += p;
	multiply (reserve, reserve, invert (tmp, p));
	tmp =
	  ressol (static_cast <
		  int >(remainder (kN, static_cast < long >(p))), p);
	if (tmp < 0)
	  tmp += p;
	multiply (reserve, reserve, tmp);
	remainder (BG[i], reserve, A);
      }
      
      B.assign (BG[0]);      // compute actual B coefficient
      for (i = 1; i < P_ONCE; i++)
	add (B, B, BG[i]);
      
      if (B.is_even ())      // assure   B = 1 mod 4
	add (B, (A.least_significant_digit () & 3) * A, B);
      
      size_FB = FB[0] + 1;   // a_inv[i] = 1/(2*A) mod p_i
      
      for (i = 2; i <= size_FB; i++)
	a_inv[i] =
	  invert (static_cast <
		  int >(remainder (A << 1, static_cast < long >(FB[i]))),
		  FB[i]);
      
      for (i = 0; i < P_ONCE; i++) {
	// vorb[i][j] = 1/A * B[i] mod p_j
	for (j = 2; j <= size_FB; j++) {
	  p = FB[j];
	  multiply (reserve, BG[i], a_inv[j] << 1);
	  if ((tmp =
	       static_cast <
	       int >(remainder (reserve, static_cast < long >(p)))) <0)
	    tmp += p;
	  
	  vorb[i][j] = tmp;
	}
      }
      
      for (j = 2; j <= size_FB; j++) {
	p = FB[j];
	SIEBS = SIEBSTART % p;
	
	tmp = static_cast < int >(remainder (-B, p));
	
	tmp1 = (tmp - SQRTkN[j]) % p;
	if (tmp1 < 0)
	  tmp1 += p;
	tmp = (tmp + SQRTkN[j]) % p;
	if (tmp < 0)
	  tmp += p;
	
	tmp2 =
	  static_cast <
	  int
	  >(multiply_mod
	    (static_cast < udigit > (tmp), udigit (a_inv[j]),
	     static_cast < udigit > (p)));
	
	tmp2 = (tmp2 + SIEBS) % p;
	if (tmp2 < 0)
	  START1[j] = tmp2 + p;
	else
	  START1[j] = tmp2;
	
	tmp2 =
	  static_cast <
	  int
	  >(multiply_mod
	    (static_cast < udigit > (tmp1), udigit (a_inv[j]),
	     static_cast < udigit > (p)));
	tmp2 = (tmp2 + SIEBS) % p;
	if (tmp2 < 0)
	  START2[j] = tmp2 + p;
	else
	  START2[j] = tmp2;
      }
      
      TMP = xgcd_left (A4_inverse, A4_inverse, kN); // determine 1/(4A) mod kN
    }
  
  else                      // no "real" computation -- use recursive formula
    {                         // first:   update of B, compute B[index_i], index_i > 0
      
      nu_2 = 0;              // nu_2 = nu_2(index_i)
      j = index_i;
      while ((j & 1) == 0) {
	nu_2++;
	j >>= 1;
      }
      
      shift_left (TMP, BG[nu_2], 1);
      
      if ((((j + 1) / 2) & 1) == 1) {
	i = -1;
	subtract (B, B, TMP);
      } else {
	i = 1;
	add (B, B, TMP);
      }
      
      size_FB = FB[0] + 1;   // determine new starting positions
      // for sieving
      
      if (i == -1) {
	for (j = 2; j <= size_FB; j++) {
	  p = FB[j];
	  START1[j] += vorb[nu_2][j];
	  if (START1[j] >= p)
	    START1[j] -= p;
	  START2[j] += vorb[nu_2][j];
	  if (START2[j] >= p)
	    START2[j] -= p;
	}
      } else {
	for (j = 2; j <= size_FB; j++) {
	  p = FB[j];
	  START1[j] -= vorb[nu_2][j];
	  if (START1[j] < 0)
	    START1[j] += p;
	  START2[j] -= vorb[nu_2][j];
	  if (START2[j] < 0)
	    START2[j] += p;
	}
      }
    }

  if (FB[2] == 2)           // note special situation for p = 2
    {
      START1[2] = 1;
      START2[2] = 1;
    }
  
  // now compute zeros of polynomials that have only one zero mod p
  // because p divides coefficient A

  square (reserve, B);      // compute coefficient -C
  subtract (reserve, kN, reserve);
  shift_left (TMP, A, 2);
  divide (TMP, reserve, TMP);
  
  for (j = 1; j <= P_TOTAL; j++)
    if (bin_index & (1 << (j - 1))) {
      p = FB[start_fb + j];
      tmp =
	invert (static_cast <
		int >(remainder (B, static_cast < long >(p))), p);
      if (tmp < 0)
	tmp += p;
      tmp2 =
	static_cast < int >(remainder (TMP, static_cast < long >(p)));
      if (tmp2 < 0)
	tmp2 += p;
      
      tmp = static_cast <int>(multiply_mod
			      (static_cast < udigit > (tmp2), 
			       static_cast < udigit > (tmp),
			       static_cast < udigit > (p)));
      START1[start_fb + j] = START2[start_fb + j] = 
	(tmp + SIEBSTART) % p;
    }
  
  /*
#ifdef DEBUG                    // check correctness of roots mod p
  if(FB[2] == 2)
    j = 3;
  else
    j = 2;
  for (; j <= FB[0] + 1; j++) 
    {
      p = FB[j];
      SIEBS = SIEBSTART % p;
      if (((A * (START1[j] - SIEBS) + B) * (START1[j] - SIEBS) +
	   (B * B - kN) / (4 * A)) % p != 0) 
	lidia_error_handler ("rational_factorization", "mpqs::compute_coeff::"
			     "found wrong polynomial in (1)");

      if (((A * (START2[j] - SIEBS) + B) * (START2[j] - SIEBS) +
	   (B * B - kN) / (4 * A)) % p != 0) 
	lidia_error_handler ("rational_factorization","mpqs::compute_coeff::"
			     "found wrong polynomial in (2)");
    }
#endif
  */
}



//-------------------------------------------------------------------
// insert ul into string *p, add blank and return pointer to first
// char after blank  --> much faster than sprintf

inline char *insert_at (char *p, unsigned long n)
{
  register int c, i, j, e;
  
  i = 0;
  do {
    p[i++] = (char) (n % 10 + '0');
  } while ((n /= 10) > 0);
  e = i;
  
  if (e > 1)
    for (i = 0, j = e - 1; i < j; i++, j--) {
      c = *(p + i);
      *(p + i) = *(p + j);
      *(p + j) = c;
    }
  p[e] = ' ';
  return (p + e + 1);
}



//***********************************************************************
// testing routine which filters correct full and LP relations out of
// candidates found in the sieving step.
//***********************************************************************

static int
teste (const bigint & kN, const bigint & A, int *FB, int *START1,
       int *START2, char *faktor, unsigned int M, double d_wurz,
       int *Q_prime, const bigint & B, unsigned int start_fb,
       unsigned int P_ONCE, unsigned int P_TOTAL, int *CANDIDATE,
       unsigned int smallstart, const bigint & A4_inverse,
       FILE * fpfull, FILE * fppart)
{
  int small_value = 0;
  unsigned int fak_i, ii, M_2, counter, upper_bound;
  bigint H, Qx, TMP;
  int x = 0, p, vorber, divides = 0, N1, N2;
  long rest;
  int geteilt = 0, ready, rest_i;
  double a, b;
  char Hstring[STLEN];
  char *faktorp;
  
  // compute the roots of the polynomial
  b = dbl (B);
  a = dbl (A);
  a *= 2;
  N1 = static_cast < int >((-b - d_wurz) / a);
  N2 = static_cast < int >((-b + d_wurz) / a);
  
  M_2 = M << 1;
  upper_bound = start_fb + P_TOTAL + 1;
  counter = 0;
  
  while ((x = CANDIDATE[counter++]) != 0) {
    // while there are candidates to test
    x -= M;
    multiply (Qx, A, x << 1);
    add (Qx, Qx, B);
    div_rem (TMP, H, Qx, kN); // needed for writing relation to file
    if (H.is_negative ())
      H.negate ();
    square (Qx, H);
    div_rem (TMP, Qx, Qx, kN);
    multiply (Qx, Qx, A4_inverse);
    div_rem (TMP, Qx, Qx, kN);
    
    faktor[0] = ' ';
    faktorp = faktor + 1;
    
    if (Qx.is_negative () && (x <= N1 || x >= N2))
      add (Qx, kN, Qx);
    else if (Qx.is_positive () && (N1 < x && x < N2)) {
      subtract (Qx, kN, Qx);
      faktorp = insert_at (faktorp, 1);
      faktorp = insert_at (faktorp, 1);
    }
    
    if (Qx.is_negative ()) {
      faktorp = insert_at (faktorp, 1);
      faktorp = insert_at (faktorp, 1);
      Qx.absolute_value ();
    }
    
    while (Qx.is_even ()) {
      divides++;
      Qx.divide_by_2 ();
    }
    
    if (divides) {
      faktorp = insert_at (faktorp, divides + 2); // '+2' because of 4*A
      faktorp = insert_at (faktorp, 2);
      divides = 0;
    } else {
      faktorp = insert_at (faktorp, 2);
      faktorp = insert_at (faktorp, 2);
    }
    
    fak_i = 2;
    divides = 0;
    ready = 0;
    
    while (fak_i++ < smallstart - 1) {
      p = FB[fak_i];
      
      vorber = (M + x) % p;
      
      if ((vorber == START1[fak_i]) || (vorber == START2[fak_i])) {
	do {
	  div_rem (TMP, rest, Qx, static_cast < long >(p));
	  
	  if (rest == 0) {
	    divides++;
	    Qx.assign (TMP);
	  }
	} while (rest == 0);
	
	faktorp = insert_at (faktorp, divides);
	faktorp = insert_at (faktorp, fak_i);
	divides = 0;
      }
    }
    
    while ((p = FB[fak_i]) != 0) {
      vorber = (M + x) % p;
      
      if ((fak_i <= upper_bound) && (fak_i > start_fb)) {
	ii = ready;
	while (ii < P_ONCE) {
	  if (p == Q_prime[ii]) {
	    geteilt = 1;
	    ready = ii + 1;
	    break;
	  }
	  ii++;
	}
      }
      
      if (geteilt == 1) {
	// p divides coefficient A ...
	if (vorber == START1[fak_i]) {
	  // ... and divides Qx
	  do {
	    div_rem (TMP, rest, Qx, static_cast < long >(p));
	    
	    if (rest == 0) {
	      divides++;
	      Qx.assign (TMP);
	    }
	  } while (rest == 0);
	}
	faktorp = insert_at (faktorp, divides + 1);
	faktorp = insert_at (faktorp, fak_i);
	
	geteilt = 0;
	divides = 0;
      } else {
	if ((vorber == START1[fak_i]) || (vorber == START2[fak_i])) {
	  do {
	    div_rem (TMP, rest, Qx, static_cast < long >(p));
	    
	    if (rest == 0) {
	      divides++;
	      Qx.assign (TMP);
	    }
	  } while (rest == 0);
	  
	  faktorp = insert_at (faktorp, divides);
	  faktorp = insert_at (faktorp, fak_i);
	  divides = 0;
	}
      }
      fak_i++;
    }
    
    
    if (Qx.is_one ()) {
      // full relation found
      small_value++;
      *(faktorp - 1) = '\0';
      
      bigint_to_string (H, Hstring);
      fprintf (fpfull, "%s 1 : %s 0\n", Hstring, faktor);
      
      /*
#ifdef DEBUG                    // test correctness of full relation
      char *p;
      bigint h1, h2;
      int e;
      
      h1.assign (1);
      p = (char *) strtok (faktor, " \n");
      while (p != NULL) {
	e = std::strtol (p,NULL,10);
	if (!e)
	  break;
	p = (char *) strtok (NULL, " \n");
	power (h2, bigint (FB[std::strtol (p,NULL,10)]), e);
	multiply (h1, h1, h2);
	p = (char *) strtok (NULL, " \n");
      }
      
      if ((H * H - h1) == 0) 
	cout<<"\nERROR: FULL RELATION IN Z !!";
	
      
      if ((H * H - h1) % kN != 0) 
	lidia_error_handler ("rational_factorization",
			     "mpqs::teste::found wrong full relation");
#endif
      */
    } 
    else 
      if ((Qx.intify (rest_i)) == 0) // LP relation found
	if (rest_i < 10000000) {
	  *(faktorp - 1) = '\0';
	  bigint_to_string (H, Hstring);
	  fprintf (fppart, "%9d @ %s : %s 0\n", rest_i, Hstring, faktor);
	
	  /*  
#ifdef DEBUG                    // test correctness of partial relation
	  char *p;
	  bigint h1, h2;
	  int e;

	  h1.assign (Qx);
	  p = (char *) strtok (faktor, " \n");
	  while (p != NULL) {
	    e = std::strtol (p,NULL,10);
	    if (!e)
	      break;
	    p = (char *) strtok (NULL, " \n");
	    power (h2, bigint (FB[std::strtol (p,NULL,10)]), e);
	    multiply (h1, h1, h2);
	    p = (char *) strtok (NULL, " \n");
	  }
	  
	  if ((H * H - h1) == 0) 
	    cout<<"\nERROR: PARTIAL RELATION IN Z !!";
	
	  
	  if ((H * H - h1) % kN != 0) 
	    lidia_error_handler ("rational_factoriztion",
				 "mpqs::teste::found wrong partial relation");
#endif
	  */
	}
  }
  return (small_value);
}



//*****************************************************************
// after LP relations have been combined, start linear system solver
// and use solutions of linear system to determine factors of kN
//*****************************************************************

bool rational_factorization::
qs_build_factors (const bigint & N,
		  const bigint & kN, unsigned int index, int *FB)
{
  rf_single_factor fact;
  FILE *faktorenmatrix;
  bigint X_quad, Y_quad, hilf, hilf2;
  int counter = 0, k, end_FB = FB[0] + 2, e;
  bool erg1, erg2;
  bool found = false;
  size_t i, zeilenindex = 0;
  char zeile[STLEN], Hstring[STLEN], Hstring2[STLEN], rest[STLEN];
  long *expo;
  long exp_i;
  char *p, *pp, *s1;
  
  if (!(expo = new long[end_FB + 1])) 
    lidia_error_handler_n ("rational_factorization",
			   "mpqs::qs_build_factors::can't allocate memory");

  memset (expo, 0, (end_FB + 1) * sizeof (long));
  transform_relations();

  if (info)
    std::cout << "\n\nLinear system built " << std::flush;
  
  // here starts the linear system solver
  
  s1 = get_name ("LANCZOS_FORMAT");
  lanczos_sparse_matrix solve_matrix(s1);

  preprocess pre;
  postprocess post;

  std::auto_ptr<index_list> correction_list = pre.process(solve_matrix);
  
  lanczos lan(solve_matrix);
  
  do {
    lan.solve();

#ifdef DEBUG
    if (lan.get_result_rank() <= 0)
      lidia_error_handler("rational_factorization", 
			  "mpqs::Lanczos found no solution!");
#endif
 
     } 
  while(lan.get_result_rank() <= 0);
  
  std::auto_ptr<lanczos_vector_block> solution = post.process(lan.get_result(),
                                                              *correction_list);

  lanczos_vector_block::result_vector_type lanczos_result = solution->result();

#ifdef DEBUG
  std::cout << "\nVerifying Lanczos solution ... " << std::flush;
  
  lanczos_sparse_matrix *matrix;
  lanczos_vector_block *vector;
  
  s1 = get_name ("LANCZOS_FORMAT");
  matrix = new lanczos_sparse_matrix (s1);
  if (!matrix) 
    lidia_error_handler("rational_factorization", 
			"mpqs::Error in Input Lanczos Matrix");

  vector = new lanczos_vector_block (matrix->number_of_columns ());
  vector->read(lanczos_result);
  
  unsigned long i2, j, h, g;
  lanczos_vector_block *result_vec;
  
  result_vec = new lanczos_vector_block (matrix->number_of_rows ());
  result_vec->clear ();
  i2 = 0;
  
  while (i2 < matrix->number_of_columns ()) 
    {
      j = 0;
      h = matrix->get_vector (i2).get_number_of_entries ();
      while (j < h) {
	g = matrix->get_vector (i2).get_entry (j);
	result_vec->put_row (g, vector->get_row (i2) ^
			     (result_vec->get_row (g)));
	j++;
      }
      i2++;
    }
  if (!result_vec->is_zero ())
    lidia_error_handler("rational_factorization", 
			"mpqs::Solution of Lanczos is wrong !!");
  else
    std::cout << "\nSolution of Lanczos is correct \n " << std::flush;
  
  delete (matrix);
  delete (vector);
  delete (result_vec);
#endif
  
  if (info)
    std::cout << "\nand solved with Block-Lanczos Algorithm\n" << std::flush;

  s1 = get_name ("RELATIONS");
  if (!(faktorenmatrix = fopen (s1, "r"))) 
    lidia_error_handler_n ("rational_factorization",
			   "mpqs::qs_build_factors::can't open relation"
			   " file");

  //--------------------------------------------------------------
  // now the combination of relations with these solutions starts

  int no_solution = 0;
  
  while (!found)            // use solutions to find X, Y
                            // with X^2 = Y^2 mod N
    {
      size_t lanczos_number_of_entries;

      X_quad.assign_one ();
      Y_quad.assign_one ();
      
      counter = 1;
      fgets (zeile, STLEN, faktorenmatrix); // skip .. relations

      if (lanczos_result[no_solution].empty())
	{
	  break;
	}	
      else
	lanczos_number_of_entries = lanczos_result[no_solution][0];

      for (i = 1; i <= lanczos_number_of_entries; i++)
	{
	  zeilenindex = lanczos_result[no_solution][i];
	  
	  while (counter <= zeilenindex) 
	    {
	      fgets (zeile, STLEN, faktorenmatrix);
	      counter++;
	    }
	  
	  sscanf (zeile, "%s %s : %s", Hstring, Hstring2, rest);
	  
	  string_to_bigint (Hstring, hilf);
	  string_to_bigint (Hstring2, hilf2);
	  
	  if (!hilf2.is_one ()) 
	    {
	      multiply (X_quad, X_quad, hilf2);
	      remainder (X_quad, X_quad, kN);
	    }
	  
	  multiply (Y_quad, Y_quad, hilf);
	  remainder (Y_quad, Y_quad, kN);
	  
	  // read single exponents and sum them up in array for later use
	  
	  if ((pp = (char *) strchr (zeile, ':'))) {
	    pp += 2;
	    p = (char *) strtok (pp, " \n");
	    while (p != NULL) {
	      e = std::strtol (p,NULL,10);
	      if (!e)
		break;
	      p = (char *) strtok (NULL, " \n");
	      expo[std::strtol (p,NULL,10)] += e;
	      p = (char *) strtok (NULL, " \n");
	    }
	  }
	}                      // solution read
      no_solution ++;

      
      for (i = 2; i <= end_FB; i++) 
	{
	  exp_i = expo[i];
	  
	  if (exp_i > 0) 
	    {
	      exp_i >>= 1;

	      while (exp_i != 0) 
		{
		  exp_i--;
		  multiply (X_quad, X_quad, FB[i]);
		}
	      if (X_quad.abs_compare (kN) >= 0)
		remainder (X_quad, X_quad, kN);
	    } 
	  else 
	    if (exp_i < 0) 
	      {
		exp_i = -((-exp_i) >> 1);
		
		while (exp_i != 0) 
		  {
		    exp_i++;
		    multiply (Y_quad, Y_quad, FB[i]);
		  }
		if (Y_quad.abs_compare (kN) >= 0)
		  remainder (Y_quad, Y_quad, kN);
	      }
	}

         add (hilf, X_quad, Y_quad);
         subtract (Y_quad, Y_quad, X_quad);

         X_quad = gcd (hilf, N);
         Y_quad = gcd (Y_quad, N);

         if (X_quad.is_one () || Y_quad.is_one ())  {
	     // no proper factor found
	     fseek (faktorenmatrix, 0, SEEK_SET);
	     memset (expo, 0, (end_FB + 1) * sizeof (long));
	     continue;
	 } 
	 else {
	     if((fact.single_exponent =
		 is_power(fact.single_base, X_quad)) == 0) {
		 fact.single_base.assign(X_quad);
		 fact.single_exponent = 1;
	     }
	     fact.factor_state = is_prime(fact.single_base) ?
		 rf_single_factor::prime : rf_single_factor::not_prime;
	     factors[index] = fact;
	     
	     if((fact.single_exponent =
		 is_power(fact.single_base, Y_quad)) == 0) {
		 fact.single_base.assign(Y_quad);
		 fact.single_exponent = 1;
	     }
	     fact.factor_state = is_prime(fact.single_base) ?
		 rf_single_factor::prime : rf_single_factor::not_prime;
	     factors[no_of_comp()] = fact;

	     found = true;
	     fseek (faktorenmatrix, 0, 0);
	     memset (expo, 0, (end_FB + 1) * sizeof (long));
	 }
    }
  fclose (faktorenmatrix);

  delete[] expo;

  return (!found);
}



//*****************************************************************
// main routine which factors one component of a rational_factorization
// with MPQS, calls all relevant functions in a loop until one factor was
// found. 
//******************************************************************

rational_factorization & rational_factorization::mpqs_impl(lidia_size_t index,
							   ecm_primes & prim) {
  bigint N (factors[index].single_base);
  bigint A, B, kN, A4_inverse;
  bigint *BG;
  long tmp, i, k, p;
  unsigned int stellenzahl;
  unsigned long counter_treff = 0;
  FILE *fpfull, *fppart;
  
  unsigned int M, size_FB, smalls = 0;
  int small_value;
  int *SQRTkN, *FB, *START1, *a_inv;
  int *CANDIDATE, *Q_prime, *Q_prime_glob;
  int *START2;
  int start_fb;
  unsigned int P_TOTAL, POLY, P_ONCE, smallstart, prozent, vergleich = 10;
  SIEBTYP *sieb = NULL, *ende = NULL, *LOGP = NULL;
  double d_wurz, T, LOGMUL;
  int **vorb;
  lidia_size_t index_i;
  lidia_size_t added_relations; 
  char *s1, *faktor;
  unsigned long last_cnt=0;
  
  i = is_power (A, N);      // perfect power handling
  if (i > 0) {
    rf_single_factor fact;
    
    if (is_prime (A, 5)) 
      {
	fact.single_base = A;
	fact.single_exponent =
	  factors[index].single_exponent * static_cast < int >(i);
	fact.factor_state = rf_single_factor::prime;
	factors[index] = fact;
      } 
    else 
      {
	rational_factorization A_fact;
	
	A_fact.assign (A);
	A_fact.factor ();
	
	int j, kk = A_fact.no_of_comp ();
	
	for (j = 0; j < kk; j++)
	  A_fact.factors[j].single_exponent *= static_cast < int >(i);
	
	divide (*this, *this, rational_factorization (N));
	multiply (*this, *this, A_fact);
      }
    return *this;
  }
  
  
  if (info) {
    std::cout << "\n\nQS                -- quadratic sieve\n == \n\n";
    std::
      cout << "number to factor: " << N << " (" << decimal_length (N) <<
      ")\n\n";
    std::cout.flush ();
  }
  
  k = compute_multiplier (N, 5, prim);
  multiply (kN, N, k);
  
  if (info)
    std::cout << "Multiplier: " << k << "\n" << std::flush;
  stellenzahl = decimal_length (kN);
  
  if (info) {
    std::cout << "\nestimated running time (user-time): ";
    zeitqs (stellenzahl, true);
    std::cout.flush ();
  }
  
  if ((stellenzahl > 80)) 
    {
      lidia_warning_handler ("rational_factorization",
			     "mpqs::Input Number too big to be factored"
			     " on one machine in reasonable time");
      stellenzahl = 80;
    }
  
  qs_read_par (stellenzahl, T, M, size_FB, P_ONCE, POLY, P_TOTAL,
	       smallstart);
  
  if (info) {
    std::cout << "\nSieve Interval: [-" << M << " , " << M << "]";
    std::cout << "\n# Factor Basis: " << size_FB << "\n";
    std::cout.flush ();
  }
  
  added_relations = static_cast < lidia_size_t >(size_FB * 0.05);
  if (added_relations < 50)  
    added_relations = 50;

  create_FB (size_FB, kN, &FB, prim);
  
#ifdef DEBUG
  global_kN = kN;
  global_FB = FB;
#endif
  
  BG = new bigint[P_ONCE + 1];
  sieb = new SIEBTYP[M << 1];
  faktor = new char[STLEN];
  
  LOGP = new SIEBTYP[size_FB + 2];
  SQRTkN = new int[size_FB + 2];
  START1 = new int[size_FB + 2];
  START2 = new int[size_FB + 2];
  a_inv = new int[size_FB + 2];
  CANDIDATE = new int[CAND_NUMBER];
  Q_prime_glob = new int[P_TOTAL];
  Q_prime = new int[P_ONCE];
  
  if (sieb == NULL || faktor == NULL || LOGP == NULL || SQRTkN == NULL ||
      START1 == NULL || START2 == NULL || a_inv == NULL || BG == NULL ||
      CANDIDATE == NULL || Q_prime_glob == NULL || Q_prime == NULL) 
    lidia_error_handler_n ("rational_factorization",
			   "mpqs::can't allocate memory");
  
  vorb = new int *[P_TOTAL];
  
  for (i = 0; i < static_cast < int >(P_TOTAL); i++) 
    {
      if (!(vorb[i] = new int[size_FB + 2])) 
	lidia_error_handler_n ("rational_factorization",
			       "mpqs::allocate vorb");
    }
  
  ende = sieb + (M << 1) - 1;
  
  // determine approximations for log(p) for FB elements p and sqrt(kN)
  // mod p (used for fast computation of sieving start points
  
  LOGMUL = (2 * static_cast < SIEBTYP > (0.5 * LiDIA::log2 (dbl (kN)) +
					 LiDIA::log2 (static_cast <
					       double >(M)) -T *
					 LiDIA::log2 (static_cast <
					       double >(FB[FB[0] + 1]))));
  LOGMUL = 127.0 / LOGMUL;
  
  tmp = size_FB + 2;
  for (i = 2; i < tmp; i++) {
    p = FB[i];
    LOGP[i] =
      static_cast < SIEBTYP >
      (LOGMUL * LiDIA::log2 (static_cast < double >(p)) * 2);

    if ((SQRTkN[i] =
	 ressol (static_cast <
		 int >(remainder (kN, static_cast < long >(p))),
		 static_cast < int >(p))) < 0)
      SQRTkN[i] += static_cast < int >(p); // compute sqrt(kN) modulo different moduli p
  }
  
  d_wurz = std::sqrt (dbl (kN));
  
  // the size of coefficient A should be approximately
  // sqrt(kN)/M, so the size of the primes p dividing
  // A should be approximately (sqrt(kN/M))^(1/P_ONCE)
  
  T = d_wurz / M;
  T = std::pow (2.0, LiDIA::log2 (T) / P_ONCE);
  
  if (T > FB[size_FB - 1]) 
    lidia_error_handler_n ("rational_factorization",
			   "mpqs: P_ONCE too small");
  
  i = 2;
  while (FB[i] < T)
    i++;
  start_fb = static_cast < int >(i); // P_TOTAL consecutive primes p[start_fb], ...,
  
  // are chosen from the factor basis
  
  if (i > 7)
    start_fb -= (P_ONCE >> 1);
  
  if (k >= FB[start_fb])    // Multiplier must not occur in factor basis
    while (k >= FB[start_fb])
      start_fb++;
  
  
  for (i = 0; i < static_cast < int >(P_TOTAL); i++) {
    Q_prime_glob[i] = FB[start_fb + i + 1]; // collect prime numbers which
                                            // will build the A coefficents
  }
  
  index_i = -1;
  
  lidia_signal sig1 (LIDIA_SIGTERM, stop_rf_mpqs), sig2 (LIDIA_SIGINT,
							 stop_rf_mpqs);
  lidia_signal sig3 (LIDIA_SIGHUP, stop_rf_mpqs), sig4 (LIDIA_SIGSEGV,
							stop_rf_mpqs);
  
  s1 = get_name ("LP_RELATIONS");
  fpfull = fopen (s1, "w");
  fprintf (fpfull, "\n");
  fclose (fpfull);
  
  s1 = get_name ("RELATIONS");
  if (!(fpfull = fopen (s1, "a"))) 
    lidia_error_handler_n ("rational_factorization",
			   "mpqs::can't open file RELATIONS");

  fprintf (fpfull, "\n");

  s1 = get_name ("NEW_PARTIALS");
  if (!(fppart = fopen (s1, "a"))) 
    lidia_error_handler_n ("rational_factorization",
			   "mpqs::can't open file NEW_PARTIALS");

  unsigned int bin_index = (1 << P_ONCE) - 1; // variable used for 
  //choosing the
  // correct A coeffs in compute_coeff
  
  while (true)              // central loop of
    // - computing polynomials and zeros
    // - sieving
    // - testing candidates of the sieve array
	{
	  if (index_i == static_cast < int >(POLY) - 1) // when all of the B's have already
            // been used, choose new A
            index_i = 0;
	  else
            index_i++;
	  
	  
	  compute_coeff (A, B, kN, FB, SQRTkN, START1, START2, P_ONCE, P_TOTAL,
			 M, vorb, Q_prime, Q_prime_glob, BG, index_i,
			 start_fb, a_inv, A4_inverse, bin_index);
	  
	  qs_sieve_interval (FB, LOGP, START1, START2, sieb, ende, M,
			     CANDIDATE, smallstart);
	  
	  small_value = teste (kN, A, FB, START1, START2, faktor, M, d_wurz, Q_prime,
			       B, start_fb, P_ONCE, P_TOTAL, CANDIDATE,
			       smallstart, A4_inverse, fpfull, fppart);
	  
	  smalls += small_value;
	  
	  prozent =
	    static_cast <
	    unsigned int >((static_cast < double >(smalls + counter_treff)
			    / static_cast < double >(size_FB)) *100);
	  
	  if (prozent >= vergleich) 
	    {
	      fflush (fppart);
	      fclose (fppart);
	      
	      counter_treff = count_partials ();
	      
	      prozent =
	        static_cast <
	        unsigned int >((static_cast < double >(smalls + counter_treff)
			        / static_cast < double >(size_FB)) *100);
	  
	      if (info && (last_cnt != (smalls + counter_treff)) ) 
		{
		  std::cout << smalls +
		    counter_treff << " (" << prozent << "%) relations ";
		  std::cout << "found. (" << smalls << " by fulls, " 
			    << counter_treff;
		  std::cout << " by partials)\n" << std::flush;
		}
	      
              last_cnt = smalls+counter_treff;
	      
              while (vergleich <= prozent)
                {
	        if (vergleich >= 90 || stellenzahl >= 75)
		  vergleich += 5;
	        else
		  vergleich += 10;
                }

	      s1 = get_name ("NEW_PARTIALS");
	      if (!(fppart = fopen (s1, "a"))) 
		lidia_error_handler_n ("rational_factorization",
				       "mpqs::can't open file NEW_PARTIALS");
	    }

	  if ((static_cast < long >(smalls) + counter_treff) > 
	      (static_cast < long >(size_FB) + added_relations) )
	    {
	      fflush (fppart);
	      fflush (fpfull);
	      fclose (fppart);
	      fclose (fpfull);
	      
              // need to include the partials that we have already.
              counter_treff = count_partials ();

	      if (info && (last_cnt != (smalls + counter_treff)) ) 
		{
	        prozent = static_cast < unsigned int >
                          ((static_cast < double >(smalls + counter_treff)
			          / static_cast < double >(size_FB)) *100);

		  std::cout << smalls +
		    counter_treff << " (" << prozent << "%) relations ";
		  std::cout << "found. (" << smalls << " by fulls, " 
			    << counter_treff;
		  std::cout << " by partials)\n" << std::flush;
		}
              last_cnt = smalls+counter_treff;
	      
	      zusam (size_FB + 2, FB, kN);

	      if (qs_build_factors (N, kN, index, FB) == false) {
		break;
	      }
	      else {
		  added_relations = static_cast <long> (static_cast < double > 
                                    (smalls + counter_treff) * 1.05) - size_FB;
		  if (info) {
		      std::cout << "\nNo non-trivial congruence found -->";
		      std::cout << " Restart sieving ...\n" << std::flush;
		  }

		  s1 = get_name ("RELATIONS");
		  if (!(fpfull = fopen (s1, "a"))) 
		    lidia_error_handler_n ("rational_factorization",
					   "mpqs::can't open file RELATIONS");
		  
		  s1 = get_name ("NEW_PARTIALS");
		  if (!(fppart = fopen (s1, "a"))) 
		    lidia_error_handler_n ("rational_factorization",
					   "mpqs::can't open file "
					   "NEW_PARTIALS");
		}
	    }
	}
      
  s1 = get_name ("RELATIONS");
  std::remove (s1);
  s1 = get_name ("LANCZOS_FORMAT");
  std::remove (s1);
  s1 = get_name ("NEW_PARTIALS");
  std::remove (s1);
  s1 = get_name ("NEW_FULLS");
  std::remove (s1);
  s1 = get_name ("LP_RELATIONS");
  std::remove (s1);
  s1 = get_name ("TEMPORARY");
  std::remove (s1);
  s1 = get_name ("SORT_LINE0");
  std::remove (s1);
  s1 = get_name ("SORT_LINE1");
  std::remove (s1);
  s1 = get_name ("SORT_LP0");
  std::remove (s1);
  s1 = get_name ("SORT_LP1");
  std::remove (s1);
  
  delete[]LOGP;
  delete[]SQRTkN;
  delete[]START1;
  delete[]START2;
  delete[]sieb;
  delete[]CANDIDATE;
  delete[]a_inv;
  delete[]FB;
  delete[]faktor;
  delete[]Q_prime_glob;
  delete[]Q_prime;
  delete[]BG;
  
  for (i = 0; i < static_cast < int >(P_TOTAL); i++)
    delete[]vorb[i];
  delete[]vorb;
  
  return *this;
}



//*********************************************************
// mpqs function for component index of a factorization
//*********************************************************


rational_factorization & rational_factorization::
mpqs_comp (lidia_size_t index) 
{
  if ((index < 0) || (index > no_of_comp ())) 
    lidia_error_handler_n ("rational_factorization",
			   "mpqs::index out of range");
  
  if (factors[index].factor_state == rf_single_factor::prime ||
      is_prime (factors[index].single_base, 8)) 
    {
      if (info)
	std::cout << "\nprime number : " << factors[index].
	  single_base << "\n";
      factors[index].factor_state = rf_single_factor::prime;
      return *this;
    }
  long B = 200000;
  
  if (decimal_length (factors[index].single_base) > 65)
    B = 500000;
  
  ecm_primes prim (2, B, 100000);
  
  trialdiv (index, 1, 180000, prim);
  prim.resetprimes (1);
  
  if (is_prime_factor (index)) // fact. found
    return *this;          // with TD
  
  mpqs (index, prim);
  compose ();
  return *this;
}



rational_factorization & rational_factorization::mpqs(lidia_size_t index,
						      ecm_primes & prim) {
    // Compute complete factorization 
    rf_single_factor& factor_ref = factors[index];
    rational_factorization f(factor_ref.single_base, 1);

    bool failed = false;
    while(!(failed || f.is_prime_factorization())) {
	failed = true;
	lidia_size_t f_len = f.no_of_comp();
	for(lidia_size_t i = 0; i < f_len; ++i) {
	    if(f.factors[i].factor_state == rf_single_factor::dont_know) {
		f.factors[i].factor_state =
		    is_prime(f.factors[i].single_base) ?
		    rf_single_factor::prime : rf_single_factor::not_prime;
	    }
	    if(f.factors[i].factor_state == rf_single_factor::not_prime) {
		prim.resetprimes(1);
		f.mpqs_impl(i, prim);
		f.refine();
		failed = false;
		break;
	    }
	}
    }

    lidia_size_t last_f_index = f.no_of_comp() - 1;
    for(lidia_size_t i = 0; i <= last_f_index; ++i) {
	f.factors[i].single_exponent *= factor_ref.single_exponent;
    }
    factor_ref = f.factors[last_f_index];
    f.factors.remove_from(last_f_index);
    factors.concat(factors, f.factors);
    
    return *this;
}

//********************************************************************
// mpqs function for rational integers
//********************************************************************

rational_factorization mpqs (const bigint & N)
{
  rational_factorization f (N);
  
  f.verbose (0);
  f.mpqs_comp (0);
  return f;
}



#ifdef LIDIA_NAMESPACE
}                               // end of namespace LiDIA
#endif
