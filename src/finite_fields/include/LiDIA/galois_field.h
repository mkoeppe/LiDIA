// -*- C++ -*-
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
//      Author  : Detlef Anton (DA), Thomas Pfahler (TPf)
//      Changes : See CVS log
//
//==============================================================================================


#ifndef LIDIA_GALOIS_FIELD_H_GUARD_
#define LIDIA_GALOIS_FIELD_H_GUARD_


#ifndef LIDIA_BIGINT_H_GUARD_
# include       "LiDIA/bigint.h"
#endif
#ifndef LIDIA_RATIONAL_FACTORIZATION_H_GUARD_
# include       "LiDIA/rational_factorization.h"
#endif
#ifndef LIDIA_FP_POLYNOMIAL_H_GUARD_
# include       "LiDIA/Fp_polynomial.h"
#endif
#ifndef LIDIA_GALOIS_FIELD_ITERATOR_H_GUARD_
# include       "LiDIA/galois_field_iterator.h"
#endif


#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif



  class galois_field_rep;
  class gf_element;


  class galois_field
  {
    friend class galois_field_rep;
    friend class gf_element;

    //
    // the C++ type we use to represent a Galois field
    //

    galois_field_rep const* rep;

    //
    // private constructor
    //

    galois_field(galois_field_rep const& internal_rep);

  public:

    //
    // constructors and destructor
    //

    galois_field();
    galois_field(const galois_field & K);
    galois_field(const bigint & characteristic, unsigned int degree = 1);
    galois_field(const bigint & characteristic, unsigned int degree,
                 const rational_factorization & fact);
    galois_field(const Fp_polynomial&);
    galois_field(const Fp_polynomial &, const rational_factorization &);
    ~galois_field();


    // iterator access
    galois_field_iterator begin() const;
    galois_field_iterator end() const;
        
    //
    // access functions
    //

    const bigint &                  characteristic()               const;
    unsigned int                    degree()                       const;
    const bigint &                  number_of_elements()           const;
    const rational_factorization &  factorization_of_mult_order()  const;
    Fp_polynomial                   irred_polynomial()             const;
    gf_element const& generator() const;


    //
    // assignment
    //

    void assign(const galois_field & K);
    galois_field& operator = (const galois_field&);

    friend void swap(galois_field&, galois_field&);

    //
    // comparisons
    //
    bool operator == (const galois_field& K) const
      {
        return rep == K.rep;
      }
    bool operator != (const galois_field& K) const
      {
        return rep != K.rep;
      }
    bool operator < (const galois_field&) const; // subfield
    bool operator > (const galois_field&) const; // superfield
  };

//
// input / output
//

  std::istream & operator >> (std::istream & in, galois_field & K);
  std::ostream & operator << (std::ostream & out, const galois_field & K);



#ifdef LIDIA_NAMESPACE
}       // end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif



#endif  // LIDIA_GALOIS_FIELD_H_GUARD_
