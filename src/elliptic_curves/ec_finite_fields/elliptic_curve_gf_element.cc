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
//      Author  : Markus Maurer (MM), Thomas Pfahler (TPf)
//      Changes : See CVS log
//
//==============================================================================================

#ifdef HAVE_CONFIG_H
#include        "config.h"
#endif
#include        "LiDIA/bigint.h"
#include        "LiDIA/gf_element.h"
#include        "LiDIA/elliptic_curve.h"
#include        "LiDIA/elliptic_curves/elliptic_curve_rep.h"
#include        "LiDIA/point.h"
#include        "LiDIA/nmbrthry_functions.h"
#include        <cassert>

const int MAX_NUMBER_ELEMENTS_IN_FIELD = 200;

#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



//
// constructor / destructor
//

elliptic_curve< gf_element >
::elliptic_curve()
        : base_elliptic_curve< gf_element > ()

{
        debug_handler("elliptic_curve< gf_element >",
                      "elliptic_curve()");
        order_of_group = rational_factorization(-1);
}



elliptic_curve< gf_element >
::elliptic_curve(const elliptic_curve< gf_element > & E)
{
        debug_handler("elliptic_curve< gf_element >",
                      "elliptic_curve(const elliptic_curve< gf_element > &)");
        this->assign(E);
}



//------------------------------------------------------------------
// this is special, we forbid char p = 3 !!
//


elliptic_curve< gf_element >
::elliptic_curve(const gf_element & x4, const gf_element & x6,
                 elliptic_curve_flags::curve_model m)

        : base_elliptic_curve< gf_element > (x4, x6, m)

{
        debug_handler("elliptic_curve< gf_element >",
                      "elliptic_curve(const gf_element&x4, "
                      "const gf_element&x6, "
                      "elliptic_curve_flags::curve_model m)");

        order_of_group = rational_factorization(-1);
}



elliptic_curve< gf_element >
::elliptic_curve(const gf_element & x1, const gf_element & x2,
                 const gf_element & x3, const gf_element & x4,
                 const gf_element & x6,
                 elliptic_curve_flags::curve_model m)

        : base_elliptic_curve< gf_element > (x1, x2, x3, x4, x6, m)

{
        debug_handler("elliptic_curve< gf_element >",
                      "elliptic_curve(x1, x2, x3, x4, x6, "
                      "elliptic_curve_flags::curve_model m)");

        order_of_group = rational_factorization(-1);
}



elliptic_curve< gf_element >
::~elliptic_curve()
{
        debug_handler("elliptic_curve< gf_element >",
                      "~elliptic_curve()");

        // Cleaning done by ~base_elliptic_curve<T>.
}



//
// assignment
//

elliptic_curve< gf_element > & elliptic_curve< gf_element >
::operator = (const elliptic_curve< gf_element > & E)
{
        debug_handler ("elliptic_curve< gf_element >",
                       "operator = (const elliptic_curve< gf_element >");

        this->assign(E);
        return *this;
}



void elliptic_curve< gf_element >
::assign (const elliptic_curve< gf_element > & E)
{
  debug_handler ("elliptic_curve< gf_element >",
                 "assign(const elliptic_curve< gf_element >");
  
  if (this != &E) {
    if (E.e == NULL)
      lidia_error_handler("elliptic_curve< gf_element >::"
                          "assign(const elliptic_curve< gf_element > &)",
                          "e == NULL");
    else
      if (e == NULL) {
        e = E.e;
        e->inc_ref_counter();
      }
      else
        if (e != E.e) {
          if (e->get_ref_counter() == 1)
            delete e;
          else
            e->dec_ref_counter();
          
          e = E.e;
          e->inc_ref_counter();
        }
    order_of_group = E.order_of_group;
  }
}



//
// invariants
//
gf_element elliptic_curve< gf_element >
::j_invariant() const
{
  debug_handler ("elliptic_curve< gf_element >",
                 "j_invariant() const");
  
  if (e == NULL)
    lidia_error_handler("elliptic_curve< gf_element >::"
                        "j_invariant()const",
                        "e == NULL");
  return e->j_invariant();
}



//
// properties
//

bool elliptic_curve< gf_element >
::is_supersingular()
{
  debug_handler ("elliptic_curve< gf_element >",
                 "is_supersingular()");
  
  bigint trace_frobenius;
  
  trace_frobenius = this->discriminant().get_field().number_of_elements()
    + 1 - this->group_order();
  
  return ((trace_frobenius %
           this->discriminant().get_field().characteristic()) == 0);
}



unsigned int elliptic_curve< gf_element >
::degree_of_definition() const
{
  debug_handler ("elliptic_curve< gf_element >",
                 "degree_of_definition() const");
  
  if (e == NULL)
    lidia_error_handler("elliptic_curve< gf_element >::"
                        "degree_of_definition()const",
                        "e == NULL");
  return e->degree_of_definition();
}



bool elliptic_curve< gf_element >
::is_cyclic()
{
  debug_handler ("elliptic_curve< gf_element >",
                 "is_cyclic()");
  bigint n1, n2;
  isomorphism_type(n1, n2);
  return (n2.is_one());
}



//******************************************************************
// COMPUTING THE ISOMORPHISM TYPE
//
//  Algorithm:
//  1. first find a point P of large order by trying a few random points
//  2. find a relation: a) choose random point
//                      b) try to find minimal x such that x*Q = lambda *P
//
//
//  For b), we use a sophisticated strategy based on the follwoing ideas:
//
//   A) for all primes q dividing a multiple of x we compute the exact
//      power of q in x (idea: like in Shoup's equal-degree factorization)
//      -- more details see below
//
//   B) if we know already x[i] * Q[i] = lambda[i] * P i=1,2, then we can
//      determine a relation for Q[1] + Q[2]. The new scalar will be the
//      lcm(x[1], x[2]), which might be an improvement for free.



//------------------------------------------------------------------
// first we define a static function which solves a part of the problem.
// For given point P and Q, the function returns the minimal
// integer x such that x*Q = lambda * P. The idea is very similar to Shoup's
// algorithm for equal degree factorization. We construct multiples x[i] of
// x and compute x[i]*Q = lambda[i]*P. Then we can combine these values
// with a gcd algorithm, for free we determine lambda.
//
// cand is a multiple of x
// rf is the prime factorization of the elliptic curve order of P, Q
// --> returns x, lambda is dropped

static bigint find_x(const point< gf_element > & P,
                     const point< gf_element > & Q,
                     const bigint & cand,
                     const rational_factorization & rf)
{
  bigint x(1), h;
  int i, expo, j;
  point<gf_element> H (P.get_curve());
  
  for (i = 0; i < rf.no_of_comp(); i++)  // try all primes
    {
      if (cand % rf.base(i) != 0)
        continue;               // we just test divisors of the candidate
      
      h.assign(cand);
      expo = 0;
      
      do           // divide with maximal power of the current prime
        {
          expo ++;
          divide(h, h, rf.base(i));
        }
      while ((h % rf.base(i)).is_zero()); 
      
      multiply(H, h, Q);
      
      for (j = 0; j <= expo; j++)
        {
          if (LiDIA::discrete_logarithm(P, rf, H, cand/h, false) != -1)
            break;
          
          multiply(H, rf.base(i), H);
          multiply(x, x, rf.base(i));
        }
    }
  
#ifdef DEBUG  
  for (i = 0; i < rf.no_of_comp(); i++)
    if (x % rf.base(i) == 0)
      if (LiDIA::discrete_logarithm((x/rf.base(i))*Q, P, false) != -1)
        {          
          std::cout<<"\nERROR in find_x:\n    Q = "<<Q;
          std::cout<<"\n      P = "<<P;
          std::cout<<"\n      x = "<<x;
          std::cout<<"\n      x/base(i) = "<<x/rf.base(i);
          std::cout<<"\n      lambda = "<<discrete_logarithm((x/rf.base(i))*Q,
                                                             P, false);
          std::cout << std::flush;
          exit(1);
        }
#endif
  return x;
}



void elliptic_curve< gf_element >
::isomorphism_type(bigint & n1, bigint & n2,
                   point< gf_element > & P,
                   point< gf_element > & Q,
                   bool info)
{
  if (e == NULL)
    lidia_error_handler("elliptic_curve<gf_element>::"
                        "isomorphism_type(n1,n2,P1,P2,info)",
                        "e == NULL");

  bigint eorder, cand, h;
  rational_factorization rf_eord;
  
  rf_eord = rf_group_order(true, 0);
  eorder = evaluate_to_bigint(rf_eord);

  // determine maximal choice "cand" for second part: 
  // cand must divide q-1, cand^2 must divide group order.

  h = this->discriminant().get_field().number_of_elements() 
    - bigint(1);
  square(h, h);
  cand = gcd(h, eorder); 

  sqrt(h, cand);
  if (h*h != cand)
    {
      sort_vector<bigint> v(square_divides_n_divisors(cand));
      cand = v[v.size()-1];
   }
  else
    cand.assign(h);

  // We use the following invariant:
  // P stores point with maximally found order (P_ord), candidate for 1st part
  // Q stores current choice for 2nd part, Q_ord is its order in E
  // Q_old stores candidate for second part, x_old = order_<P> (Q_old)
  

  bigint P_ord(1), Q_ord, x_old(0), m_factor;
  point<gf_element> P_mfactor(*this), Q_old(*this);
  point<gf_element> H(*this);
  int tests;
  
  gf_element table[10];
  int next_table_index = 0;
  bool P_changed = false;
  rational_factorization rf_ord_P_mfactor;

  do 
    {
      tests = 1;
      do        // search for point with "maximal" order, store it in P
        {
          Q = random_point();
          tests++;
          if ((P_ord * Q).is_zero()) {// to avoid unnecessary order comp., this
                                      // Q does not give any improvement
            continue;
          }
          else
            Q_ord = order_point(Q);

          if (Q_ord > P_ord)  // Q is a "better" point than P, change P
            {
              P_changed = true;  // since it might decrease 'cand'
              P = Q; 
              P_ord.assign(Q_ord);
            }

          if (P_ord == eorder)  // curve is cyclic, we are done.
            {
              n1.assign(eorder); 
              n2.assign_one();
              Q = point<gf_element> (*this);  // = O

              return;
            }
        }
      while(tests < 5 || P_ord < eorder/cand);

      //----------------------------------------------------
      // second part, find relations:
      // If we find a new candidate for P, we change P. P always
      // stores the point with largest order P_ord
      //
      // P_ord is the minimal size for the first part
      // cand is the maximal size for the second part
     
      bigint x, order_Q, lambda;
     
      tests = 1;

      if (P_changed)  // check whether we can decrease 'cand'
        {
         if (info)
            std::cout<<"\nFound cyclic subgroup of order "<<P_ord<<std::flush;

          cand = gcd(cand, eorder/P_ord);  
          
          if (cand * P_ord < eorder)  // P_ord too small, search for new P
                                      // with bigger order (with first loop)
            {
              x_old.assign_zero();
              continue;
            }

          // determine m_factor as the biggest divisor of 
          // (group order / cand) which is coprime to cand 
          // --> must divide order of 1st part

          m_factor = eorder / cand;  
          x = gcd(m_factor, cand);

          while (x > 1)
            {
              m_factor /= x;
              x = gcd(m_factor, cand);
            }

          multiply(P_mfactor, m_factor, P);
          rf_ord_P_mfactor = rf_order_point(P_mfactor);
          P_changed = false;
        }

      rational_factorization cand_rf(eorder / m_factor);

      do 
        {                    
        start_loop2:
          tests ++;
          do {
              Q = random_point();

              // order_<P>(Q) must divide cand 

              multiply(Q, m_factor, Q);
              order_Q = order_point(Q, cand_rf);
              multiply (Q, order_Q / gcd(order_Q, cand), Q);
          } while(Q.is_zero());


          // First check whether this point can improve on a
          // previously computed point Q_old with order_<P>(Q_old) = x_old
	  // Improvement <=> x_old * Q \not \in <P>
	  // If there is no improvement, then we choose a new Q-candidate

          if (x_old > 0)    // there is already some candidate Q_old
            {
              int ll = 0;

              multiply(H, x_old, Q);
	      if (H.is_zero())          // no improvement !!
		goto start_loop2; 
              
	      // Practical improvement: x_old*Q often has
	      // the same value for different points Q, we store these in a
	      // small table to reduce the number of DL computations
	      // This table stores x-coordinates of points H' where we know
	      // that x_old * H' \in <P>

              while (ll < next_table_index) 
                {
                  if (H.get_x() == table[ll]) // no improvement !!
                    {
                      ll = next_table_index + 10;
                      break;
                    }
                  ll++;
                }
              if (ll == next_table_index + 10)  // iterate and choose new Q
                {
		  continue;
                }
              
              lambda = LiDIA::discrete_logarithm(P_mfactor, rf_ord_P_mfactor, 
						 H, order_point(H));

              if (lambda != -1)  // x_ord * Q \in <P>, no improvement
                {
                  if (next_table_index < 10 && !H.is_zero())
                    table[next_table_index++] = H.get_x();
                  continue;
                }
            }
          
	  // Now we have found a point Q which might give us an improvement
          // So find minimal x such that x*Q in <P>   <==>
	  // x*m_factor*Q in <m_factor*P> (and order(m_factor*P) <= ord_P)
	  
          x = find_x (P_mfactor, m_factor*Q, cand, rf_ord_P_mfactor); 

          if (info)
            std::cout<<"\nFound relation x * "<<Q<<" in < "<<P<<" > "
	      " for minimal x = "<<x<<std::flush;

          if (x_old.is_zero()) // initialize Q_old and x_old
            {
              x_old = x; Q_old = Q; 
	      table[next_table_index++] = Q.get_x();
            }
          else
            {
	      next_table_index = 0; // delete smaller order elements from table

	      if(x > x_old && (x % x_old).is_zero()) // x multiple of x_old
	        {
		  x_old = x;
		  Q_old = Q;
		}
	      else  
	        { // x(Q_old + Q) | lcm(x, x_old).
                  // Find the largest among x, x_old, x(Q_old + Q)

                  add(H, Q_old, Q);
                  n1 = find_x(P_mfactor, m_factor*H, cand, rf_ord_P_mfactor);

		  if (n1 > x_old && n1 > x)
		    {
		      x_old = n1;
                      Q_old = H;
		    }
                  else if (x > x_old)
                    {
                      x_old = x;
                      Q_old = Q;
                    }
                }
            }

          if (x_old == cand && (P_ord * cand) == eorder) // found a solution?
            {
              n1 = P_ord;
              n2 = cand;
              Q = Q_old;

              if (info)
                std::cout<<"\n"<<std::flush;
      

#ifdef DEBUG
              assert(order_point(P) == n1);
              assert(order_point(Q) == n2);
              assert(find_x(P, Q, eorder, rf_eord) == n2);
#endif
              return;
            }
        }
      while(tests < 20);

      // after 10 iterations no solution, start from the beginning.
      // We have to delete the information for different Q's, since P will
      // change !!

      P_changed = false;
      x_old.assign_zero();
      next_table_index = 0; // delete entries from table
    }
  while(1);
}



void elliptic_curve< gf_element >
::isomorphism_type(bigint & n1, bigint & n2, bool info)
{
	debug_handler ("elliptic_curve< gf_element >",
		       "isomorphism_type(bigint & n1, bigint & n2)");

	point< gf_element > P(*this), Q(*this);
	this->isomorphism_type(n1, n2, P, Q, info);
}



//
// random points in GF(p^d). It is expected that the curve
// is defined over GF(p^d), p is the characteristic
//

point< gf_element > elliptic_curve< gf_element >
::random_point(unsigned int d) const
{
  debug_handler ("elliptic_curve< gf_element >",
		 "random_point(unsigned int)");
  
  if (e == NULL)
    lidia_error_handler("elliptic_curve< gf_element >::"
			"random_point(d) const",
			"no curve defined");
  
  if ((d % this->degree_of_definition()) != 0)
    lidia_error_handler("elliptic_curve< gf_element >::"
			"random_point(d) const",
			"E(GF(p^d)) not defined");
  
  gf_element x, y;
  gf_element f1, f0;
  bool root_found;
  // Polynomial Y^2 + f1*Y + f0
  
  x.assign(this->discriminant());
  
  do {
    root_found = false;
    
    if (d == 0)
      x.randomize();
    else
      x.randomize(d);
    
    f1 = e->get_a1()*x + e->get_a3();
    f0 = -(x*(x*(x+e->get_a2())+e->get_a4())+e->get_a6());

    root_found = y.solve_quadratic(f1, f0);
    if (root_found && d != 0)  // check whether in subfield
      if (d % y.relative_degree() != 0)
	root_found = false;
  } while (!root_found);
  
  point< gf_element > P(x, y, *this);
  
  return P;
}



//
// group order: count points in GF(2^deg), deg very
//  small (no checking for the correctness of des is done
//  since it's a private function)
//

bigint elliptic_curve< gf_element >
::small_counting_points_gf2 (unsigned int deg) const
{
  debug_handler ("elliptic_curve< gf_element >",
		 "small_counting_points_gf2(unsigned int)");
  
  long qn, i, res;
  gf_element gen(this->discriminant().get_field()), x, h1, h0;
  bool found_prim_element;

  qn = 1 << deg;
  
  if (qn > MAX_NUMBER_ELEMENTS_IN_FIELD)
    lidia_error_handler("elliptic_curve< gf_element >",
			"small_counting_points_gf2",
			"field of definition has too much elements");
  
  if (qn == 2)
    gen.assign_one();
  else 
    {
      rational_factorization rf(qn - 1);
      rf.factor_completely();
      
      // now handle all elements in mult. group of GF(2^deg)
      // we first determine a primitive element for GF(2^deg)
      
      
      do 
	{
	  found_prim_element = true;
	  gen.randomize(deg);
	  for (i = 0; i < rf.no_of_comp(); i++) 
	    {
	      power(x, gen, (qn - 1)/rf.base(i));
	      if (x.is_one()) 
		{
		  found_prim_element = false;
		  break;
		}
	    }
	} 
      while (found_prim_element == false);
    }
  
  res = 2; // O and 2-TP with x(2-TP) = 0
  x.assign(gen);
  
  for (i = 1; i < qn; i++) 
    {
      h1 = e->get_a1()*x + e->get_a3();
      if (h1.is_zero())
	res++;
      else {
	h0 = (((x + e->get_a2())*x + e->get_a4())*x + e->get_a6());
	square(h1, h1);
	multiply(h0, h0, inverse(h1));
	h1 = h0;
	
        for (int i = 1; i < deg; i++)  // compute reduced trace
	  {
	    square(h1, h1);
	    add(h0, h0, h1);
	  }

	if (h0.is_zero())
	  res += 2;
      }
      multiply(x, x, gen);
    }
  return bigint(res);
}



//------------------------------------------------------------------
// same as above, but for p >= 3

bigint elliptic_curve< gf_element >
::small_counting_points_gfp (unsigned int deg) const
{
  debug_handler ("elliptic_curve< gf_element >",
		 "small_counting_points_gfp(unsigned int)");
  
  bigint r;
  long qn, i, res;
  gf_element gen(this->discriminant().get_field()), x, h;
  bool found_prim_element;
  
  power(r, this->discriminant().get_field().characteristic(), deg);
  r.longify(qn);
  
  if (qn > MAX_NUMBER_ELEMENTS_IN_FIELD)
    lidia_error_handler("elliptic_curve< gf_element >",
			"small_counting_points_gfp",
			"field of definition has too much elements");
  
  rational_factorization rf(qn - 1);
  rf.factor_completely();
  
  // now handle all elements in mult. group of GF(p^deg)
  // we first determine a primitive element for GF(p^deg)
  
  do {
    found_prim_element = true;
    do {
      gen.randomize(deg);
    } while (gen.is_zero());
    
    for (i = 0; i < rf.no_of_comp(); i++) {
      power(x, gen, (qn-1)/rf.base(i));
      if (x.is_one()) {
	found_prim_element = false;
	break;
      }
    }
  } while (found_prim_element == false);
  
  res = 1; // the point at infinity
  
  // test for point with x-coordinate 0
  
  square(h, e->get_a3());
  add(h, h, 4*e->get_a6());
  if (h.is_zero())  // exactly one point with X coord 0
    res ++;
  else {
    power(h, h, (qn - 1) >> 1);
    if (h.is_one())
      res += 2;
  }
  
  // now test all other x-coord's in mult. group
  
  x.assign(gen);
  
  for (i = 1; i < qn; i++) {
    h = e->get_a1()*x + e->get_a3();
    square(h, h);
    h += 4*(((x + e->get_a2())*x + e->get_a4())*x + e->get_a6());
    if (h.is_zero())  // exactly one point with X coord x
      res++;
    else {
      power(h, h, (qn - 1) >> 1);
      if (h.is_one())
	res += 2;
    }
    multiply(x, x, gen);
  }
  return bigint(res);
}



//-----------------------------------------------------------

bigint elliptic_curve< gf_element >
::group_order(int info)
{
  debug_handler ("elliptic_curve< gf_element >",
		 "group_order()");
  
  if (e == NULL)
    lidia_error_handler("elliptic_curve< gf_element >::"
			"group_order() const",
			"e == NULL");

  unsigned int def_degree = this->degree_of_definition();
  bigint res(-1), q, l, u, Porder;
  int i;
  
  res = evaluate_to_bigint(this->order_of_group);
  
  if (res != -1)
    return res;
  
  point< gf_element > P(*this), O(*this), max_P(*this);
  O.assign_zero();
  
  // we determine the group order over the field of definition
  // and lift it.
  
  bigint order_of_group_bi;
  bool found = false;
  
  power(q,
	this->discriminant().get_field().characteristic(),
	def_degree);

  sqrt(u, q << 2);
  l = q + bigint(1) - u;    // lower bound of Hasse interval
  if (l.is_negative())
    l.assign_one();
  u = q + bigint(1) + u;    // upper bound of Hasse interval  
  
  int counter = 0;
  bigint max_order(-1);
  
  do {
    if (q <= MAX_NUMBER_ELEMENTS_IN_FIELD)   
      // q is the order of elements in the field of definition
      {
	res = small_counting_points(def_degree);
	order_of_group = res;
	if (info)
	  std::cout << "\nGroup Order (over Field of Definition"
	    ") = " << res << std::flush;
	found = true;
      }
    else 
      {
	counter ++;
	P = random_point(def_degree);
	res = LiDIA::bg_algorithm(P, O, l, u, info);
	if (res == -1)
	  lidia_error_handler("elliptic_curve< gf_element >",
			      "group_order::no result found");
	
	// now we combine the new information with already computed infos.
	// order_of_group (and the bigint version in order_of_group_bi)
	// always stores a multiple of all found point orders 
	
	if (max_order == -1) 
	  {
	    order_of_group = rational_factorization(res);
	    order_of_group_bi = res;
	    max_order = order_point(P, order_of_group);
	    max_P.assign(P);
	  }
	else 
	  {
	    bigint h;

	    h = lcm(res, order_of_group_bi); 
	    if (h > order_of_group_bi)     // order_of_group grows
	      {
		multiply(order_of_group, order_of_group, h/order_of_group_bi);
		order_of_group_bi = h;
	      }
	    
	    Porder = order_point(P, order_of_group); 
            if (Porder > max_order)
	      {
		max_P.assign(P);
		max_order = Porder;
	      }
	
	    res = q + bigint(1) - ((q + bigint(1)) % max_order);
	    if (res < l)
	      res += max_order;
	    
	    if (res + max_order > u && res - max_order < l) 
	      {
		found = true; // group order uniquely determined 
		              // stored in res and order_of_group
		divide(order_of_group, order_of_group, order_of_group_bi/res);
		order_of_group.refine();  
		
		if (info)
		  std::cout << "\nGroup Order (over Field of Definition)"
		    " = " << res << std::flush;
	      }
	  }
      }
  } 
  while (!found && counter < 10);
    
  if (!found) 
    {
      divide(order_of_group, order_of_group, order_of_group_bi/max_order);
      order_of_group.refine();

      if (info)
	{
	  std::cout << "\nCurve seems to be of \"bad\" isomorphism type, ";
	  std::cout << " found only point "<<max_P<<" of";
	  std::cout << " order " << max_order << std::flush;
	}
      
      point< gf_element > Q(*this), H(*this);
      rational_factorization rf_max_P_ord, rf_Q_ord;
      bigint h, Q_ord;
      
      // now we know a point P_max of order max_order which has maximal
      // order among all random points tried so far

      rf_max_P_ord = rf_order_point (max_P, order_of_group); // its order
      
      do 
	{
	awal_loop:
	  Q = random_point (def_degree);

	  if (!((max_order * Q).is_zero())) // oops, order(Q) > order(P_max)
	    {
	      res = LiDIA::bg_algorithm(Q, O, l, u, info);
	      if (res == -1)
		lidia_error_handler("elliptic_curve< gf_element >",
				    "group_order::no result found");

	      res = lcm(res, max_order); // has to be > max_order
	      multiply(order_of_group, order_of_group, res / max_order);

	      bigint h;
	      h = order_point(Q, order_of_group);
	      if (h > max_order)
		{
		  max_P = Q;
		  max_order = h;
		  rf_max_P_ord = rf_order_point (max_P, order_of_group);
		}
	      else
		{
		  h = order_point(Q + max_P, order_of_group);
		   if (h > max_order)
		     {
		       max_P = Q + max_P;
		       max_order = h;
		       rf_max_P_ord = rf_order_point (max_P, order_of_group);
		     }
		}

	      res = q + bigint(1) - ((q + bigint(1)) % max_order);
	      if (res < l)
		res += max_order;
	      
	      if (res + max_order > u && res - max_order < l) 
		{
		  found = true; // group order uniquely determined 
		                // stored in res and order_of_group
		  order_of_group_bi = evaluate_to_bigint(order_of_group);
		  multiply(order_of_group, order_of_group, res);
		  divide(order_of_group, order_of_group, order_of_group_bi);
		  order_of_group.refine();  
		  
		  if (info)
		    std::cout << "\nGroup Order (over Field of Definition)"
		      " = " << res << std::flush;
		  goto before_lifting;
		}
	    }
	      
	  rf_Q_ord = rf_order_point (Q, order_of_group);//order(Q)|order(max_P)
	  Q_ord = evaluate_to_bigint (rf_Q_ord);

	  
	  if (LiDIA::discrete_logarithm(max_P, rf_max_P_ord, 
					Q, Q_ord, false) != -1)  // Q in <P>
	      goto awal_loop;

	  sort_vector< bigint > div = divisors(rf_Q_ord);
	  
	  for (i = 1; i < div.size(); i++) 
	    {
	      // test only divisors of ord(Q)
	      // a better strategy using only prime divisors
	      // is possible, but not yet implemented

	      if (i < div.size()-1)
		{
		  multiply(rf_max_P_ord, rf_max_P_ord, div[i-1]);
		  divide(rf_max_P_ord, rf_max_P_ord, gcd(max_order, div[i]));
		  rf_max_P_ord.refine(); // make sure no negative exponents

#ifdef DEBUG
		  assert(order_point(div[i]*max_P, rf_max_P_ord) == 
			 evaluate_to_bigint(rf_max_P_ord));
#endif

		  h = LiDIA::discrete_logarithm(div[i] * max_P,
						rf_max_P_ord,
						div[i] * Q, Q_ord / div[i], 
						false);
		}
	      
	      if (h != -1 || i == div.size()-1) 
		{
#ifdef DEBUG
		  if (i == div.size() - 1)
		    assert(LiDIA::discrete_logarithm(div[i] * max_P,
						     rf_max_P_ord,
						     div[i] * Q, 
						     Q_ord / div[i], 
						     false) != -1);
#endif
		  multiply(res, max_order, div[i]);
		  
		  h = q + bigint(1) - ((q + bigint(1)) % res);
		  if (h < l)
		    h += res;
		  
		  if (h + res > u && h - res < l) 
		    {
		      multiply(order_of_group, order_of_group, 
                               h/order_of_group);
		      order_of_group.refine();
		      res = h;
		      found = true;
		      i = div.size() + 1; // to break for loop
		    }
		}
	    }	      
	}
      while (! found);
    }

  // order_of_group stores the rational_factorization of the order
  // over the small field
  // res stores the group order over the small field

before_lifting:
  
  if (def_degree < this->discriminant().absolute_degree()) 
    {
      bigint c1, ci, ci_1(2), h, qpower;
      long number_of_liftings =
	this->discriminant().absolute_degree() / def_degree;
      
      c1 = q + bigint(1) - res; // trace over field of def.
      ci = c1; // trace over GF(q^i)
      qpower = q; // q^i
      
      for (i = 2; i <= number_of_liftings; i++) 
	{
	  h = ci_1; // now compute the new trace value
	  ci_1 = ci;
	  subtract(ci, c1 * ci, q * h);
	  multiply(qpower, qpower, q);
	}
      
      res = qpower + bigint(1) - ci;
      multiply(order_of_group, order_of_group, res / (q + 1 - c1));
      order_of_group.refine();
      
      if (info)
	std::cout << "\nGroup Order (over Extension Field"
	  ") = " << res << std::flush;
      
#ifdef DEBUG
      assert(probabilistic_test_of_group_order(res, 5, false) == true);
#endif
    }
  return res;
}



//------------------------------------------------------------------
rational_factorization  elliptic_curve< gf_element >
::rf_group_order(bool prime_factorization, int info)
{
	bigint h = evaluate_to_bigint(this->order_of_group);

	if (h == -1)
		h = group_order(info); // and set the internal variable

	if (prime_factorization == true)
		if (! order_of_group.is_prime_factorization())
			order_of_group.factor_completely();

	return this->order_of_group;
}



//------------------------------------------------------------------

void elliptic_curve< gf_element >
::set_group_order (const rational_factorization & f)
{
    debug_handler("elliptic_curve< gf_element >",
		  "set_group_order()");

    order_of_group.assign(f);

    if (!probabilistic_test_of_group_order(evaluate_to_bigint(f), 5, false))
      lidia_error_handler("elliptic_curve< gf_element >::"
			  "set_group_order(f)",
			  "f not the correct group order");
}



//-----------------------------------------------------------------------

bool elliptic_curve< gf_element >::probabilistic_test_of_group_order
(const bigint & res, lidia_size_t tests, bool use_internal_info)
{
	if (use_internal_info == true) {
		bigint h = evaluate_to_bigint(order_of_group);
		if (h != -1)
		  return (h == res);
	}

	point< gf_element > P(*this);
	bigint h_c, h_two_sqrt_p;

	//
	// check Hasse interval
	//
	bigint pn = this->get_a4().get_field().number_of_elements();

	h_c = pn + 1 - res;
	h_c.absolute_value();

	sqrt(h_two_sqrt_p, pn << 2);

	if (h_c > h_two_sqrt_p)
		return false;

	//
	// check with randomly chosen points
	//
	for (int i = 0; i < tests; i++) {
		P = this->random_point();
		multiply(P, res, P);
		if (! P.is_zero())
			return false;
	}
	return true;
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
