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
//      Author  : Thomas Pfahler (TPf)
//      Changes : See CVS log
//
//==============================================================================================


#include <memory>


#ifdef HAVE_CONFIG_H
# include       "config.h"
#endif
#include        "LiDIA/galois_field.h"
#include        "LiDIA/finite_fields/galois_field_rep.h"
#include        "LiDIA/gf_element.h"



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



//***********************************************************************
//*                     class galois_field_rep                          *
//***********************************************************************


  namespace { // anonymous

    inline
    gf_flags::gf_rep_enum
    set_flag(const bigint &p, int deg) {
      gf_flags::gf_rep_enum gf_rep;
      if (deg == 1)
        gf_rep = gf_flags::GFp;
      else if (deg > 1) {
        if (p == 2)
          gf_rep = gf_flags::GF2n;
        else
          gf_rep = gf_flags::GFpn;
      }
      else gf_rep = gf_flags::UNKNOWN;
      
      return gf_rep;
    }


    void
    update_factorization(rational_factorization &a,
                         const rational_factorization &b) {
      // refine factorization in a using factors in b
      for(int i = 0; i < b.no_of_comp(); i++)
        a.refine(b.base(i));
    }
        
  } // anonymous namespace


  galois_field_rep *galois_field_rep::head = 0;


  galois_field_rep::
  galois_field_rep(const bigint& characteristic,
                   unsigned int degree,
                   const Fp_polynomial& pol,
                   const rational_factorization& fact)
    : gf_flags(),
      I2(),
      poly_mod(),
      p(characteristic),
      n(degree),
      p_pow_n(),
      p_pow_n_minus_1(),
      gen(0),
      next(0) {
    // don't initialize the member function pointers, they are all
    // set in this->init_functions() below

    power(this->p_pow_n, characteristic, degree);

    this->gf_rep = set_flag(characteristic, degree);
    this->init_functions();
    if (!characteristic.is_zero()) {
      this->init_pol(pol);
      this->init_fact(fact);
    }
  }


  galois_field_rep::~galois_field_rep()
  {
    if(this->gen) {
      // not strictly necessary because `delete 0' is legal. But
      // let's be cautious...
      delete gen;
      gen = 0;
    }     
  
    if (this == galois_field_rep::head)
      galois_field_rep::head = next;
    else {
      galois_field_rep *K = galois_field_rep::head;
      while (K->next != this) {
        if (K->next == 0)
          lidia_error_handler("galois_field_rep",
                              "~galois_field_rep()::internal error");
        K = K->next;
      }
      K->next = next;
    }
  }



  galois_field_rep* galois_field_rep::search(galois_field_rep* h,
                                             const bigint &p, unsigned int deg)
    // search for field with p^deg elements in list starting with h
    // return 0 if no such field could be found
  {
    while (h != 0) {
      if (h->degree() == deg && h->characteristic() == p)
        return h;
      h = h->next;
    }
    return NULL;
  }



//
// static "constructors" for galois_field_rep
// (so we can keep a list of galois_field_rep)
//
  galois_field_rep*
  galois_field_rep::create(const bigint & characteristic, unsigned int degree,
                           const rational_factorization & fact) {
    galois_field_rep* K = search(head, characteristic, degree);
    if (K != 0)     // i.e. we have found a field with characteristic^degree elts.
    {
      update_factorization(K->p_pow_n_minus_1, fact);
      K->rc.inc_ref_counter(); // increase ref. counter
    }
    else {
      // need new galois_field_rep
      K = new galois_field_rep(characteristic, degree, Fp_polynomial(), fact);
      K->next = galois_field_rep::head;
      galois_field_rep::head = K; // insert into list of fields
    }
    return K;
  }

  galois_field_rep*
  galois_field_rep::create(const Fp_polynomial &pol,
                           const rational_factorization &fact)
  {
    debug_handler("galois_field_rep", "create(Fp_polynomial&, ...)");

    bigint p = pol.modulus();
    int deg = pol.degree();
    if (deg < 1 || p.is_zero())
      lidia_error_handler("galois_field_rep",
                          "create(Fp_polynomial, ...): bad polynomial");

    galois_field_rep* K = search(head, p, deg);
    while (K != 0) {
      // so far, we have found a field with p^deg elements, but we don't
      // know whether it has the right irred. polynomial
      if (K->irred_polynomial() == pol) {
        // yeah, we found the field in our list
        update_factorization(K->p_pow_n_minus_1, fact);
        K->rc.inc_ref_counter(); // increase ref. counter
        //output();
        return K;
      }
      K = search(K->next, p, deg);
    }

    // need new galois_field_rep
    K = new galois_field_rep(p, deg, pol, fact);
    K->next = galois_field_rep::head;
    galois_field_rep::head = K; // insert into list of fields

    return K;
  }



  void galois_field_rep::init_p_pow_n_minus_1() const
  {
    debug_handler("galois_field_rep", "init_p_pow_n_minus_1()");

    rational_factorization &F = p_pow_n_minus_1;

    if (!(F.is_prime_factorization())) {
      // p-1 divides p^n-1 for n>1  
      if (n > 1 && p != 2)  F.refine(p-1);

      F.factor_completely();
      if (!(F.is_prime_factorization()))
        lidia_error_handler("galois_field_rep",
                            "init_p_pow_n_minus_1: factorization failed");
    }
  }



  const rational_factorization&
  galois_field_rep::factorization_of_mult_order() const
  {
    init_p_pow_n_minus_1();
    return p_pow_n_minus_1;
  }



  void galois_field_rep::init_functions()
  {
    switch(gf_rep) {
      case GF2n:      construct = &galois_field_rep::construct_GF2n;
        destruct = &galois_field_rep::destruct_GF2n;
        copy = &galois_field_rep::copy_GF2n;
        as0 = &galois_field_rep::as0_GF2n;
        as1 = &galois_field_rep::as1_GF2n;
        iseq = &galois_field_rep::iseq_GF2n;
        is0 = &galois_field_rep::is0_GF2n;
        is1 = &galois_field_rep::is1_GF2n;
        neg = &galois_field_rep::neg_GF2n;
        add = &galois_field_rep::add_GF2n;
        sub = &galois_field_rep::sub_GF2n;
        mul = &galois_field_rep::mul_GF2n;
        inv = &galois_field_rep::inv_GF2n;
        sqr = &galois_field_rep::sqr_GF2n;
        hash = &galois_field_rep::hash_GF2n;
        get_pol_rep = &galois_field_rep::get_pol_rep_GF2n;
        set_pol_rep = &galois_field_rep::set_pol_rep_GF2n;
        break;
      case GFpn:      construct = &galois_field_rep::construct_GFpn;
        destruct = &galois_field_rep::destruct_GFpn;
        copy = &galois_field_rep::copy_GFpn;
        as0 = &galois_field_rep::as0_GFpn;
        as1 = &galois_field_rep::as1_GFpn;
        iseq = &galois_field_rep::iseq_GFpn;
        is0 = &galois_field_rep::is0_GFpn;
        is1 = &galois_field_rep::is1_GFpn;
        neg = &galois_field_rep::neg_GFpn;
        add = &galois_field_rep::add_GFpn;
        sub = &galois_field_rep::sub_GFpn;
        mul = &galois_field_rep::mul_GFpn;
        inv = &galois_field_rep::inv_GFpn;
        sqr = &galois_field_rep::sqr_GFpn;
        hash = &galois_field_rep::hash_GFpn;
        get_pol_rep = &galois_field_rep::get_pol_rep_GFpn;
        set_pol_rep = &galois_field_rep::set_pol_rep_GFpn;
        break;
      case UNKNOWN:
      case GFp:       construct = &galois_field_rep::construct_GFp;
        destruct = &galois_field_rep::destruct_GFp;
        copy = &galois_field_rep::copy_GFp;
        as0 = &galois_field_rep::as0_GFp;
        as1 = &galois_field_rep::as1_GFp;
        iseq = &galois_field_rep::iseq_GFp;
        is0 = &galois_field_rep::is0_GFp;
        is1 = &galois_field_rep::is1_GFp;
        neg = &galois_field_rep::neg_GFp;
        add = &galois_field_rep::add_GFp;
        sub = &galois_field_rep::sub_GFp;
        mul = &galois_field_rep::mul_GFp;
        inv = &galois_field_rep::inv_GFp;
        sqr = &galois_field_rep::sqr_GFp;
        hash = &galois_field_rep::hash_GFp;
        get_pol_rep = &galois_field_rep::get_pol_rep_GFp;
        set_pol_rep = &galois_field_rep::set_pol_rep_GFp;
        break;
      default:        lidia_error_handler("galois_field_rep",
                                          "init_functions: type not supported");
    }
  }



  void galois_field_rep::init_pol(const Fp_polynomial& pol)
  {
    Fp_polynomial tmp;
    char *s, *s2;
    int i;
    switch (gf_rep) {
      case GF2n:      if (pol.degree() < 0)
        I2.init(n);
      else {
        s = new char[n*10];
        s2 = new char[16];
        sprintf(s, "%d", n); // lead coeff
        for (i = n-1; i >= 0; i--)
          if (!pol[i].is_zero()) {
            sprintf(s2, " %d", i);
            strcat(s, s2);
          }
        I2.init(s, n);
        delete[] s;
        delete[] s2;
      }
        break;
      case GFpn:      if (pol.degree() < 0) {
        build_irred(tmp, p, n);
        poly_mod.build(tmp);
      }
      else
        poly_mod.build(pol);
        break;
      case GFp:       if (pol.degree() < 0) {
        tmp.assign_x(p);
        poly_mod.build(tmp);
      }
      else
        poly_mod.build(pol);
        break;
      default:        lidia_error_handler("galois_field_rep",
                                          "init_pol: type UNKNOWN not supported");
    }
  }



  void galois_field_rep::init_fact(const rational_factorization &fact)
  {
    if (fact.no_of_comp() == 0)//1 && fact.base(0).is_one() && fact.sign() == 1)
    {
      // no factorization was provided
      if (p == 2 && n > 30 && !(n & 1)) {
        bigint q;
        unsigned int i = n;
        while (!(i & 1) && i > 20)
          i >>= 1;
        shift_left(q, bigint(1), i);
        multiply(p_pow_n_minus_1, q-bigint(1), q+bigint(1));
        i <<= 1;
        while (i < n) {
          square(q, q);
          multiply(p_pow_n_minus_1, p_pow_n_minus_1, q+bigint(1));
          i <<= 1;
        }
      }
      else
        p_pow_n_minus_1.assign(p_pow_n - 1);
    }
    else {
      bigrational tmp = evaluate(fact);
      if (tmp+1 != p_pow_n)
        lidia_error_handler("galois_field_rep",
                            "galois_field_rep: bad factorization");
      p_pow_n_minus_1.assign(fact);
    }
  }



  Fp_polynomial galois_field_rep::irred_polynomial() const
  {
    Fp_polynomial pol;
    int i;
    switch (gf_rep) {
      case GFp:
      case GFpn:
        return poly_mod.modulus();
      case GF2n:
        pol.set_modulus(2);
        pol.set_coefficient(n);
        pol.set_coefficient(0);
        if (I2.invert_p == &info_gf2n::tri_invert)
          pol.set_coefficient(I2.exp1);
        else if (I2.invert_p == &info_gf2n::pent_invert) {
          pol.set_coefficient(I2.exp1);
          pol.set_coefficient(I2.exp2);
          pol.set_coefficient(I2.exp3);
        }
        else if (I2.invert_p == &info_gf2n::general_invert) {
          for (i = I2.anz_exponents-1; i >= 0; i--)
            pol.set_coefficient(I2.exponents[i]);
        }
        else
          lidia_error_handler("galois_field_rep",
                              "irred_polynomial(): internal error");
        break;
      default:
        lidia_error_handler("galois_field_rep",
                            "irred_polynomial(): no field specified");
    }
    return pol;
  }


  gf_element const&
  galois_field_rep::generator() const {
    if(!this->gen) {
      std::auto_ptr<gf_element> new_gen (new gf_element());
      new_gen->assign_primitive_element(galois_field(*this));
      this->gen = new_gen.release();
    }

    return *this->gen;
  }


#ifdef LIDIA_NAMESPACE
}       // end of namespace LiDIA
#endif
