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
//	Author	: Volker Mueller (VM)
//	Changes	: See CVS log
//
//==============================================================================================

#ifdef HAVE_CONFIG_H
# include	"config.h"
#endif
#include	"LiDIA/bigint.h"
#include	"LiDIA/gf_element.h"
#include	"LiDIA/point.h"
#include	"LiDIA/timer.h"

#ifdef DEBUG
#include        <cassert>
#endif


#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



// first we define ec_tables which are used in BG algorithm.

class ec_table_entry {
private:
  udigit h_value; // hash_value of x-coord of point
  int index; // index, -1 denotes empty
public:
  
  ec_table_entry()
  {
    index = -1;
  }
  ~ec_table_entry() { }
  
  ec_table_entry & operator = (const ec_table_entry & e)
  {
    h_value = e.h_value;
    index = e.index;
    return *this;
  }
  
  void set(udigit h, int i)
  {
    h_value = h;
    index = i;
  }
  
  int get_index()
  {
    return index;
  }
  udigit get_hash()
  {
    return h_value;
  }
};
  
  //---------------------------------------------------------------

//
// constructors / destructor
//
point< gf_element >
::point() : base_point< gf_element > ()
{
	debug_handler("point< gf_element >",
		      "point()");
	info = 0;
}



point< gf_element >
::point(const base_point< gf_element > & bp) : base_point< gf_element > (bp)
{
	debug_handler("point< gf_element >",
		      "point(const base_point< gf_element > &)");
	info = 0;
}



point< gf_element >
::point(const gf_element & xp,
	const gf_element & yp,
	const elliptic_curve< gf_element > & e) : base_point< gf_element > ()
{
	debug_handler("point< gf_element >",
		      "point(x, y, e)");

	this->assign(xp, yp, e);
}



point< gf_element >
::point(const gf_element & xp,
	const gf_element & yp,
	const gf_element & zp,
	const elliptic_curve< gf_element > & e) : base_point< gf_element > ()
{
	debug_handler("point< gf_element >",
		      "point(x, y, z, e)");

	this->assign(xp, yp, zp, e);
}



point< gf_element >
::point(const point< gf_element > & P) : base_point< gf_element > ()
{
	this->assign(P);
}



point< gf_element >
::point(const elliptic_curve< gf_element > &e) : base_point< gf_element > ()
{
	debug_handler("point< gf_element >",
		      "point(const elliptic_curve< gf_element > &");

	this->assign_zero(e);
}



point< gf_element >
::~point()
{
	debug_handler("point< gf_element >",
		      "~point()");
}



//
// assignment
//

point< gf_element > & point< gf_element >
::operator = (const point< gf_element > & P)
{
	this->assign(P);
	return *this;
}



//
// Computes a point (xx,y) in affine model and
// (xx,y,1) in projective model, if possible.
// Otherwise, the point is initialized by zero.
//
bool point< gf_element >::set_x_coordinate(const gf_element & xx)
{
	if (ec.get_model() == elliptic_curve_flags::PROJECTIVE) {
		z.assign_one(xx.get_field());
	}

	gf_element f1(xx), f0(xx), y(xx); // Polynomial Y^2 + f1*Y + f0

	add(f1, ec.get_a1()*xx, ec.get_a3());
	f0 = -(xx*(xx*(xx + ec.get_a2())+ec.get_a4())+ec.get_a6());

	if (!y.solve_quadratic(f1, f0)) {
		this->assign_zero(ec);
		return false;
	}
	else {
		this->assign(xx, y, ec);
		return true;
	}
}



//
// Frobenius map
//

void point< gf_element >
::frobenius(unsigned int d)
{
	bigint qd;
	gf_element nx, ny;

	power(qd, (this->get_x()).get_field().characteristic(),
	      ec.degree_of_definition() * d);

	power(nx, get_x(), qd);
	power(ny, get_y(), qd);

	if (ec.get_model() == elliptic_curve_flags::PROJECTIVE) {
		gf_element nz;
		power(nz, get_z(), qd);
		assign(nx, ny, nz);
	}
	else
		assign(nx, ny);
}



//---------------------------------------------------------------------

point< gf_element > frobenius(const point< gf_element > & P, unsigned int d)
{
	point< gf_element > Q(P);
	Q.frobenius(d);
	return Q;
}



//---------------------------------------------------------------------
// point order
//

rational_factorization rf_order_point(const point< gf_element > & P,
				      rational_factorization & mult_order)
{
   bigint m_order, h;
   rational_factorization Pord(1);
   unsigned int  j, i;
   point< gf_element > res(P);
   
   if (P.is_zero())
     return rational_factorization(1);
   
   m_order = evaluate_to_bigint(mult_order);

#ifdef DEBUG
   assert((m_order * P).is_zero());
#endif

   if (! mult_order.is_prime_factorization())
     mult_order.factor_completely();
   
   mult_order.refine();  // remove negative exponents if possible !!
   
   for (i = 0; i < static_cast<unsigned int>(mult_order.no_of_comp()); i++)
     {
       if (mult_order.exponent(i) <= 0)
	 {
	   lidia_error_handler("point_gf_element.c",
			       "factorization of order has a non-positive exponent");
	   return -1;
	 }
       power(h, mult_order.base(i), mult_order.exponent(i));
       divide(h, m_order, h);
       multiply(res, h, P);

       j = 0;
       while (j <= static_cast<unsigned int>(mult_order.exponent(i)) 
	      && !res.is_zero()) {
	 j ++;
	 multiply(Pord, Pord, mult_order.base(i));
	 multiply(res, mult_order.base(i), res);
       }
       if (j >= static_cast<unsigned int>(mult_order.exponent(i))
	   && !res.is_zero())
	 {
	   lidia_error_handler("point_gf_element.c",
			       "no multiple of point order");
	   return -1;
	 }
     }
   return (Pord);
}



//-----------------------------------------------------------------

rational_factorization rf_order_point(const point< gf_element > & P)
{
	rational_factorization f;
	point< gf_element > Q(P);
	elliptic_curve< gf_element > e(P.get_curve());

	f = e.rf_group_order(true); // true since we need prime factorization !

	return rf_order_point(P, f);
}



//-----------------------------------------------------------------

bigint order_point(const point< gf_element > & P,
		   rational_factorization & mult_order)
{
	return evaluate_to_bigint(rf_order_point(P, mult_order));
}



//-----------------------------------------------------------------

bigint order_point(const point< gf_element > & P)
{
	return evaluate_to_bigint(rf_order_point(P));
}



//-----------------------------------------------------------------
// discrete logarithm
//

//***************************************************************
// Babystep Giantstep Algorithms for elliptic curves:
//
// We have two versions. The first one is the standard version which 
//
// -->uses +/- trick for points
// -->uses hash values of x-coordinates
//
// Output:  return 0 <= lower <= x <= upper with x*P = Q
//          if such a x does not exist, then return -1.
//**************************************************************

const long MAX_BG_STEPS = 3000000;

bigint bg_algorithm(const point< gf_element > & PP,
		    const point< gf_element > & QQ,
                    const bigint & lower,
		    const bigint & upper,
		    bool info)
{
  if (PP.is_zero() && !QQ.is_zero())
    return (-1);

  if (lower.is_zero() && QQ.is_zero())
    return 0;

  if (PP.get_curve() != QQ.get_curve())
    lidia_error_handler("point_gf_element::bg_algorithm",
                        "Points P and Q on different curves");

  if (lower.is_negative() || upper < lower)
    lidia_error_handler("point_gf_element::bg_algorithm",
                        "lower bound > upper bound");

  point < gf_element > P(PP), Q(QQ);  // to make sure that points are affine
  bool was_projective = false;
 
  if (QQ.is_projective_point())
    {
      was_projective = true;
      P.make_affine_point(true);
      Q.make_affine_point(true);
    }

  point< gf_element > H(P.get_curve()), H2(P.get_curve()), 
    H3(P.get_curve());

  bigint h;
  long number_baby, size_table, number_giant;
  long i, j, l;
  udigit P_hash;

  if (info)
    std::cout<<"\nBabystep Giantstep algorithm: "<<std::flush;
  
  if (upper - lower < 30)    // for very small intervals
    {
      if (info)
        std::cout<<"\nTesting "<<(upper - lower) <<" possibilities ... "
                 <<std::flush; 
      H.assign(lower * P);
      h.assign(lower);
      if (H == Q)
	{
	  if (was_projective)
	    Q.make_projective_point(true);
	  return h;
	}      
      do
        {
          add(H, H, P);
          inc(h);
          if (H == Q)
	    {
	      if (was_projective)
		Q.make_projective_point(true);
	      return h;
	    }
        }
      while (h <= upper); 
      if (was_projective)
	Q.make_projective_point(true);
      return (-1);
    }
  
  //**** otherwise we use the Babystep Giantstep idea **************
  
  sqrt(h, (upper - lower) >> 1);  // compute number of babysteps
  h++;
  
  if (h > MAX_BG_STEPS)
    h = MAX_BG_STEPS;
  
  h.longify(number_baby);
  
  size_table = number_baby << 1;   // to use hashing, we double the table size
  ec_table_entry *HT;
  HT = new ec_table_entry[size_table];
  
  multiply(H2, lower, P);    //  Q - lower*P 
  subtract(H2, Q, H2);
  H.assign_zero();

  //****** Babysteps, store hash(x(i * P))  *********************
  
  if (info)
    std::cout << " (#Babysteps = " << number_baby << std::flush;
  
  for (i = 1; i <= number_baby; i++) 
    {
      add(H, H, P);

      if (H == H2)   // i * P = Q - lower* P
        {
	  if (info)
	    std::cout<<") "<<std::flush;

          delete[] HT;
#ifdef DEBUG
          assert((lower + bigint(i)) * P ==  Q);
#endif
	  if (was_projective)
	    Q.make_projective_point(true);
          return (lower + bigint(i));
        }

      if (H.is_negative_of(H2))
	{
	  // We have Q = (lower - i) * P.
	  if (info)
	    std::cout<<") "<<std::flush;
	  
          delete[] HT;
	  h = lower - i;   

	  // Compute x = (lower - i) + k * ord(P), k \in Z, 
	  // s.t. lower <= x <= upper.

	  // First compute ord(P)
	  bigint P_ord;
	  if (Q.is_zero())         
	  { // h is a multiple of order(P), use it
	    rational_factorization rf(h);
	    P_ord = order_point(P, rf);
	  }
	  else
	  {  
	    P_ord = order_point(P);
	  }

	  // Now move h into the requested interval.
	  bigint k = (upper - h) / P_ord;
	  h += k * P_ord;   // h <= upper < h + ord(P)
	  if (h < lower)
	  {
	    h = -1;
	  }
	  
	  if (was_projective)
	  {
	    Q.make_projective_point(true);
	  }
	  return (h);    // solution in [lower ... upper]
	}
	  
      if (!H.is_zero()) 
        {
          P_hash = hash(H.get_x());   // insert into hash table
          j = P_hash % size_table;
          while (HT[j].get_index() != -1)  // HT[j] is not empty
            {
              j++;
              if (j == size_table)
                j = 0;
            }
          HT[j].set(P_hash, i);
        }
    }

  //****** Giantsteps ***************************************************/
  
  if ((((upper - lower)/(number_baby << 1))
       + bigint(1)).longify(number_giant))
    lidia_error_handler("point_gf_element.c::bg_algorithm",
                        "#Giantsteps does not fit into a long -- "
			"consider using the ECO package");
  
  if (info)
    std::cout << ", #Giantsteps = " << number_giant << ") " << std::flush;
  
  bigint step_size = number_baby << 1; 
  multiply_by_2(H, H);

  for (j = 0; j <= number_giant; j++)
    {
      if (H2.is_zero()) 
        {
          delete[] HT;
	  h = lower + j * step_size;

#ifdef DEBUG
          assert(h*P == Q);
#endif
	  if (was_projective)
	    Q.make_projective_point(true);
	  if (h <= upper)  // solution in [lower, ... , upper]
	    return h;
	  else
	    return -1;
        }  
      
      P_hash = hash(H2.get_x());
      l = P_hash % size_table;
      
      while (HT[l].get_index() != -1)  
        {
          if (HT[l].get_hash() == P_hash) 
            {
              multiply(H3, bigint(HT[l].get_index()), P);
              if (H3 == H2) 
                {
                  h = lower + bigint(HT[l].get_index())
                    + j * step_size;
#ifdef DEBUG
                  assert(h * P ==  Q);
#endif
		  if (was_projective)
		    Q.make_projective_point(true);
                  delete [] HT;
		  if (h <= upper)  // solution in [lower, ... , upper]
		    return h;
		  else
		    return -1;
                }
	      if (H3.is_negative_of(H2)) 
                {
                  h = lower - bigint(HT[l].get_index()) 
                    + j * step_size;
#ifdef DEBUG
                  assert(h * P ==  Q);
#endif
		  if (was_projective)
		    Q.make_projective_point(true);
                  delete [] HT;
		  if (h <= upper)
		    return h;
		  else 
		    return -1;
                }
            }
          l++;
          if (l == size_table)
            l = 0;
        }
      subtract(H2, H2, H);
    }
  delete[] HT;

  if (was_projective)
    Q.make_projective_point(true);
  return -1;
}

//------------------------------------------------------------
// 2nd version of Babystep Giantstep, that does not use +/- trick,
// but assumes knowledge of some integers x, m such that
// solution = x mod m. This can be used to reduce the table size.
//

bigint bg_algorithm(const point< gf_element > & PP,
		    const point< gf_element > & QQ,
                    const bigint & lower,
		    const bigint & upper,
		    const bigint & xx,
		    const bigint & m,
		    bool info)
{
  if (PP.is_zero() && !QQ.is_zero())
    return (-1);

  if (PP.get_curve() != QQ.get_curve())
    lidia_error_handler("point_gf_element::bg_algorithm",
                        "Points P and Q on different curves");

  if (lower.is_negative() || upper < lower)
    lidia_error_handler("point_gf_element::bg_algorithm",
                        "lower bound > upper bound");

  if (m.is_negative() || xx.is_negative() || m.is_zero())
    lidia_error_handler("point_gf_element::bg_algorithm",
                        "Either value or modulus is non-positive");

  point < gf_element > P(PP), Q(QQ);  // to make sure that points are affine
  bool was_projective = false;
  bigint x (xx % m);
 
  if (QQ.is_projective_point())
    {
      was_projective = true;
      P.make_affine_point(true);
      Q.make_affine_point(true);
    }

  point< gf_element > H(P.get_curve()), H2(P.get_curve()), 
    H3(P.get_curve()), mP(P.get_curve()), xP(x * P);  

  if (xP == Q)  // test whether x already enough
    {
      if (was_projective)
	Q.make_projective_point(true);
      return x;
    }

  bigint h;
  long number_baby, size_table, number_giant;
  long i, j, l;
  udigit P_hash;

  multiply(mP, m, P);

  if (info)
    {
      std::cout<<"\nBabystep Giantstep algorithm: ";
      std::cout<<" (expect solution = "<<x;
      std::cout<<" mod "<<m<<") ";
      std::cout<<std::flush; 
    }
  
  if (upper - lower < 30 * m)    // for very small intervals
    {
      if (info)
        std::cout<<"\nTesting "<<(upper - lower)/m + 1<<" possibilities ... "
                 <<std::flush; 
      H.assign(lower * P + xP);
      h.assign(lower + x);
      if (H == Q)
	{
	  if (was_projective)
	    Q.make_projective_point(true);
	  return h;
	}      
      do
        {
          add(H, H, mP);
          add(h, h, m);
          if (H == Q)
	    {
	      if (was_projective)
		Q.make_projective_point(true);
	      return h;
	    }
        }
      while (h <= upper); 
      if (was_projective)
	Q.make_projective_point(true);
      return (-1);
    }
  
  //**** otherwise we use the Babystep Giantstep idea **************
  
  sqrt(h, (upper - lower) / m);  // compute number of babysteps
  h++;
  
  if (h > MAX_BG_STEPS)
    h = MAX_BG_STEPS;
  
  h.longify(number_baby);
  
  size_table = number_baby << 1;   // to use hashing, we double the table size
  ec_table_entry *HT;
  HT = new ec_table_entry[size_table];
  
  multiply(H2, lower, P);    //  Q - (lower + x)*P 
  add(H2, H2, xP);
  subtract(H2, Q, H2);
  H.assign_zero();

  //****** Babysteps, store hash(x(i * m * P))  *********************
  
  if (info)
    std::cout << " (#Babysteps = " << number_baby << std::flush;
  
  for (i = 1; i <= number_baby; i++) 
    {
      add(H, H, mP);

      if (H == H2)   // i * m * P = Q - (lower + x) * P
        {
	  if (info)
	    std::cout<<") "<<std::flush;

          delete[] HT;
#ifdef DEBUG
          assert((lower + bigint(i) * m + x) * P ==  Q);
#endif
	  if (was_projective)
	    Q.make_projective_point(true);
          return (lower + bigint(i) * m + x);
        }
      
      if (!H.is_zero()) 
        {
          P_hash = hash(H.get_x());   // insert into hash table
          j = P_hash % size_table;
          while (HT[j].get_index() != -1)  // HT[j] is not empty
            {
              j++;
              if (j == size_table)
                j = 0;
            }
          HT[j].set(P_hash, i);
        }
    }

  //****** Giantsteps ***************************************************/
  
  if ((((upper - lower)/ (number_baby * m))
       + bigint(1)).longify(number_giant))
    lidia_error_handler("point_gf_element.c::bg_algorithm",
                        "#Giantsteps does not fit into a long -- "
			"consider using the ECO package");
  
  if (info)
    std::cout << ", #Giantsteps = " << number_giant << ") " << std::flush;
  
  bigint step_size = m * bigint(number_baby);

  for (j = 0; j <= number_giant; j++)
    {
      if (H2.is_zero()) 
        {
          delete[] HT;
#ifdef DEBUG
          assert((lower + x + j * m * number_baby)*P == Q);
#endif
	  if (was_projective)
	    Q.make_projective_point(true);
          return(lower + x + j * m * number_baby);
        }  
      
      P_hash = hash(H2.get_x());
      l = P_hash % size_table;
      
      while (HT[l].get_index() != -1)  
        {
          if (HT[l].get_hash() == P_hash) 
            {
              multiply(H3, bigint(HT[l].get_index()), mP);
              if (H3 == H2) 
                {
                  h = lower + x + bigint(HT[l].get_index()) * m
                    + j * step_size;
#ifdef DEBUG
                  assert(h * P ==  Q);
#endif
		  if (was_projective)
		    Q.make_projective_point(true);
                  delete [] HT;
                  return h;
                }
            }
          l++;
          if (l == size_table)
            l = 0;
        }
      subtract(H2, H2, H);
    }
  delete[] HT;

  if (was_projective)
    Q.make_projective_point(true);
  return (-1);
}





//----------------------------------------------------------------
// Here we start all functions for discrete logarithms.
// First two helper functions for Pollard rho algorithm, 
// Teske's Variant, with tabel to store random points
//

const unsigned int long number_mults = 20; // number of mult's in Teske's walk
const unsigned long table_size = 50000;    // table of stored points
const double load_factor = 0.5;            // if 0.5 * size(table) is used,
                                           // we start deleting points

//-------------- function for correctness check of found point ----
// returns DL mod ord_P iff H (= a*P + b*Q) == ai*P + bi*Q


inline bigint check_res(const point < gf_element > & P,
			const point < gf_element > & Q,
			const point < gf_element > & H, 
			const bigint & a,
			const bigint & b, 
			const bigint & ai, 
			const bigint & bi,
			const bigint & ord_P)
{
  point <gf_element> G (Q.get_curve());
  bigint res;

  if (a == ai || b == bi)  // trivial match !!
    return -1;

  G= ai * P + bi * Q;

  if (H == G)
    {
      if (! xgcd_left(res, bi-b, ord_P).is_one())
        return -1;

      res = (res * (a - ai)) % ord_P;  // compute DL
      if (res.is_negative())
        add(res, res, ord_P);
      return res;
    }

  if (H.is_negative_of(G))
    {
      if (! xgcd_left(res, -bi-b, ord_P).is_one())
        return -1;
      
      res = (res * (a + ai)) % ord_P;
      if (res.is_negative())
        add(res, res, ord_P);
      return res;
    }
  return -1;
}

//---------------------------------------------------------------------
// Function to check for match of H = a*P + b*Q in table; if match found,
// DL is returned; otherwise -1 is returned

inline bigint search_table (gf_element * & table, bigint * & at, bigint * & bt,
			    const point <gf_element> & P,
			    const point<gf_element> & Q, 
			    const point<gf_element> & H,
			    const bigint & a, const bigint & b, 
			    const bigint & ord_P)
{
  int index;
  bigint res;

  if (H.is_zero())  // solution found: a * P + b * Q = O
    if (xgcd_left(res, -b, ord_P).is_one())
      {
	res = (res * a) % ord_P;
	if (res.is_negative())
	  add(res, res, ord_P);
	return res;
      }

  gf_element h = H.get_x();

  index = hash(h) % table_size;  
  while (!table[index].is_zero())
    {
      if (table[index] == h)
        {
          res = check_res(P, Q, H, a, b, at[index], 
                          bt[index], ord_P);
          if (res != -1)
            return res;
        }     
      index ++;
      if (index == table_size)
        index = 0;
    }
  return -1;
}

//--------------------------------------------------------------------
// Function to check for match of H (= a * P + b * Q) in table, if so, 
// return DL; otherwise insert (H.x, a, b) into tables. If the number of
// stored elements is larger than load_factor * size(table), a random
// element is deleted.

inline
bigint insert_and_delete_table(gf_element * & table, 
			       bigint * & at, bigint * & bt,
			       int & load, const point <gf_element> & P,
			       const point<gf_element> & Q, 
			       const point<gf_element> & H,
			       const bigint & a, const bigint & b, 
			       const bigint & ord_P)
{
  int index;
  bigint res;
  gf_element h;
  random_generator rg;
  
  if (H.is_zero())  // solution found: a * P + b * Q = O
    if (xgcd_left(res, -b, ord_P).is_one())
      {
	res = (res * a) % ord_P;
	if (res.is_negative())
	  add(res, res, ord_P);
	return res;
      }

  h = H.get_x();
  index = hash(h) % table_size;  // insert H.x into table and check
  while (!table[index].is_zero())
    {
      if (table[index] == h)
        {
          res = check_res(P, Q, H, a, b, at[index], 
                          bt[index], ord_P);
          if (res != -1)
            return res;
        }     
      index ++;
      if (index == table_size)
        index = 0;
    }

  table[index] = h;       // insert into table
  at[index] = a;
  bt[index] = b;
  load ++;
  
  if (load > table_size * load_factor)  // delete a random element
    {
      int next;
      do
        {
          rg >> index;
          index = index % table_size;
        }
      while (table[index].is_zero());
          
      next = index + 1;
      if (next == table_size)
        next = 0;
      
      while (hash(table[next]) == hash(table[index]))
        {
          table[index] = table[next];
	  at[index] = at[next];
	  bt[index] = bt[next];
	  index = next;
	  next++;
	  if (next == table_size)
	    next = 0;
        }
      table[index].assign_zero();
      at[index].assign_zero();
      bt[index].assign_zero();
      load --;
    }
  return -1;
}




bigint point< gf_element >::pollard_rho(point< gf_element > & QQ,
					 const bigint & ord_P,
					 bool info)
{
  point< gf_element > Pi[number_mults], H, H2, P(*this), Q(QQ);
  bigint ai[number_mults], bi[number_mults], a, b, a2, b2;
  bigint res;
  int i, iter, max_iter;
  udigit hsh;
  gf_element *table;
  bigint *at, * bt;
  int index, load = 0;
  bool projective_curve = false;

  if (info)
  {
    std::cout << "\nPollard Rho  DL algorithm: " << decimal_length(ord_P);
    std::cout << " dec. digits, " << std::flush;
  }

  if (Q.is_projective_point())  
    projective_curve = true;
  P.make_affine_point(true);
  Q.make_affine_point(true);

  // compute a upper bound before restart

  sqrt(res, 4*ord_P);   

  if (res.intify(max_iter))
    max_iter = 10000000;

  if (max_iter < 10000)
    max_iter = 10000;

  // allocate and initialize tables 

  table = new gf_element[table_size];
  at = new bigint [table_size];
  bt = new bigint [table_size];

  for (i = 0; i < table_size; i++)
    {
      table[i] = gf_element(Q.get_x().get_field());
      at[i] = bt[i] = 0;
    }


 restart:            // restart starts from here !!

  iter = 1;

  for (i = 0; i < number_mults; i++)  // initialize Teske's walks
    {
      Pi[i] = point < gf_element > (P.get_curve());       
      ai[i] = randomize(ord_P);
      bi[i] = randomize(ord_P);
      add(Pi[i], ai[i] * P, bi[i] * Q);
    }
  
  a = randomize(ord_P);  // initialize start points H_i and H_{2i}
  b = randomize(ord_P);
  add(H, a * P, b * Q);

  hsh = hash(H.get_x()) % number_mults;
  add(a2, a, ai[hsh]);
  if (a2 >= ord_P)
    subtract(a2, a2, ord_P);
  add(b2, b, bi[hsh]);
  if (b2 >= ord_P)
    subtract(b2, b2, ord_P);
  add(H2, H, Pi[hsh]);


  // main loop for walking

  while (1) 
    {
      res = search_table(table, at, bt, P, Q, H, a, b, ord_P);
      if (res != -1)  //solution found with match on H_i in table
        {
          if (info)
	    std::cout << iter << " iterations " << std::flush;
          break;
        }

      res = search_table(table, at, bt, P, Q, H2, a2, b2, ord_P);
      if (res != -1)  //solution found with match on H_{2i} in table
        {
          if (info)
	    std::cout << iter << " iterations " << std::flush;
          break;
        }
      
      iter ++;
      if (iter > max_iter)
         {
           if (info)
              std::cout<<"\nRestart random walk ... "<<std::flush;
           goto restart;
         }
      
      // compute next two points H_{i+1} and H_{2i+2}
      
      hsh = hash(H.get_x()) % number_mults;
      add(a, a, ai[hsh]);
      if (a >= ord_P)
        subtract(a, a, ord_P);
      add(b, b, bi[hsh]);
      if (b >= ord_P)
        subtract(b, b, ord_P);
      add(H, H, Pi[hsh]);
    
      hsh = hash(H2.get_x()) % number_mults;
      add(a2, a2, ai[hsh]);
      if (a2 >= ord_P)
        subtract(a2, a2, ord_P);
      add(b2, b2, bi[hsh]);
      if (b2 >= ord_P)
        subtract(b2, b2, ord_P);
      add(H2, H2, Pi[hsh]);

      res = insert_and_delete_table(table, at, bt, load, P, Q, 
				    H2, a2, b2, ord_P);
      if (res != -1)
	{
          if (info)
	      std::cout << iter << " iterations " << std::flush;
	  break;
        }
      
      hsh = hash(H2.get_x()) % number_mults;
      add(a2, a2, ai[hsh]);
      if (a2 >= ord_P)
        subtract(a2, a2, ord_P);
      add(b2, b2, bi[hsh]);
      if (b2 >= ord_P)
        subtract(b2, b2, ord_P);
      add(H2, H2, Pi[hsh]);

      res = insert_and_delete_table(table, at, bt, load, P, Q, 
				    H2, a2, b2, ord_P);
      if (res != -1)
	{
          if (info)
	     std::cout << iter << " iterations " << std::flush;
	  break;
	}
    }
  
  delete[] table;
  delete[] at;
  delete[] bt;

  if (projective_curve)        // as was the input !!
    QQ.make_projective_point(true);
  
  return res;
}






//****************************************************************
// Pohlig-Hellman DL Algorithm:
//
// determine 0 <= x < ord(P) such that x*P = Q,
// if Q \not \in< P >, return -1.
//****************************************************************

bigint discrete_logarithm(const point<gf_element> & P, 
                          const point<gf_element> & Q, 
                          bool info)
{
  rational_factorization ord_P;
  bigint ord_Q;

  if (P.ec != Q.ec)
    lidia_error_handler("DL", "points on different curves");
  
  if (Q.is_zero())
    return 0;
  
  if (P.is_zero())
    return -1;
  
  ord_P = rf_order_point(P);
  ord_Q = order_point(Q);

  return discrete_logarithm(P, ord_P, Q, ord_Q, info);
}


bigint discrete_logarithm(const point<gf_element> & P, 
                          const rational_factorization & ord_P,
                          const point<gf_element> & Q, 
                          const bigint & ord_Q,
                          bool info)
{
  if (P.ec != Q.ec)
    lidia_error_handler("DL","points on different curves");

  if (Q.is_zero())
    return 0;

  if (P.is_zero())
    return -1;

  bigint bi_ord_P = evaluate_to_bigint(ord_P);

#ifdef DEBUG
  assert((bi_ord_P*P).is_zero());
#endif

  if (!(bi_ord_P % ord_Q).is_zero())
    {
      if (info)
        std::cout << "\nQ no element of subgroup <P>\n\n" << std::flush;
      return (-1);       
   }

  point<gf_element> PP(P), QQ(Q), H(P.ec);

  bigint dl = 0, res, h, dl_mod_pn, dl_old, modulus_old;
  bigint modulus = 1, pn, help;
  lidia_size_t i, j;
  lidia_size_t unchanged_result = 0;  
         // negative results mean that -x res didn't change

  for (i = 0; i < ord_P.no_of_comp(); i++)
    {
      H.assign_zero(); 
      multiply(PP, bi_ord_P/ord_P.base(i), P);

      j = ord_P.exponent(i);
      pn.assign_one();
      dl_mod_pn.assign_zero();

      if (unchanged_result >= 1)  // already computed result sufficient ?
        {
          unchanged_result = 0;
          multiply(QQ, dl, P);
          if (QQ == Q)
            {
              if (info)
                std::cout<<"\nearly abort: x = "<<dl<<"\n"<<std::flush;
              return dl;
            }
        }
      if (unchanged_result <= -1)  // -(already computed result) sufficient ?
        {
          unchanged_result = 0;
          multiply(QQ, dl - modulus, P);
          if (QQ == Q)
            {
              if (info)
		std::cout<<"\nearly abort: x = "<<dl - modulus<<"\n"<<std::flush;
              return bi_ord_P + dl - modulus;
            }
        }

      while(j >= 1)
        {
          multiply(pn, pn, ord_P.base(i));
          divide(h, bi_ord_P, pn);
          multiply(QQ, h, Q-H);

          if (QQ.is_zero())   // dl = 0  mod new_modulus
            res.assign_zero();
          else
            {
              if (decimal_length(ord_P.base(i)) <= 8)
		{
		  res = bg_algorithm(PP, QQ, 1, ord_P.base(i), info);
#ifdef DEBUG
		  if (res != -1)
		    assert(res*PP == QQ);
#endif

		}
              else
		{ // sol = dl mod modulus
		  bigint bound = ord_P.base(i);
		  res = -1;
		  while (res == -1)
		    {
		      res = bg_algorithm(PP, QQ, 1, bound, dl, 
					 modulus, info);
#ifdef DEBUG
		      if (res != -1)
			assert(res*PP == QQ);
#endif
		      bound <<= 2;
		      if (info)
			std::cout<<"\nResizing search interval to "<<bound
				 <<" ... "<<std::flush;
		      if (decimal_length(bound / modulus) >= 14) 
			{
			  if (info)
			    std::cout<<"\nSwitching to Pollard Rho algorithm. "
				     << std::flush;
			  res = PP.pollard_rho(QQ, ord_P.base(i), info);
			}
		    }
		}
            }
              
          if (res == -1)  
            {
              if (info)
                std::cout<<"\nQ not element of subgroup <P>\n\n" << std::flush;
              return (-1);
            }
          else
            {
              divide(help, res*pn, ord_P.base(i));
              dl_mod_pn += help;
              H += help*P;
            }
          j--;
        }
      if (info)
        std::cout<<"\nx = "<<dl_mod_pn<<" mod "<<pn <<std::flush;

      dl_old = dl; modulus_old = modulus;
      dl = chinese_remainder(dl, modulus, dl_mod_pn, pn);
      multiply(modulus, modulus, pn);

      if (dl == dl_old)
        unchanged_result ++;
      else
        if (dl - modulus == dl_old - modulus_old)
          unchanged_result --;
    }
  if (info)
    std::cout <<"\n" <<std::flush;
  return dl;
}




#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
