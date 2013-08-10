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


#ifndef LIDIA_BASE_POINT_CC_GUARD_
#define LIDIA_BASE_POINT_CC_GUARD_


#ifndef LIDIA_ELLIPTIC_CURVE_FLAGS_H_GUARD_
# include	"LiDIA/elliptic_curve_flags.h"
#endif
#ifndef LIDIA_BASE_POINT_H_GUARD_
# include	"LiDIA/elliptic_curves/base_point.h"
#endif



#ifdef LIDIA_NAMESPACE
# ifndef IN_NAMESPACE_LIDIA
namespace LiDIA {
# endif
#endif



//
// assignments
//

template< class T >
void
base_point< T >::set_verbose (int a)
{
	this->info = a;
}



template< class T >
void
base_point< T >::assign (const base_point< T > & P)
{
	debug_handler("base_point< T >", "assign(const base_point< T > &");

	if (&P != this) {
		this->is_0 = P.is_0;
		this->ec.assign(P.ec);
		this->info = P.info;

		this->x.assign(P.x);
		this->y.assign(P.y);
		this->z.assign(P.z);
	}
}



template< class T >
void
base_point< T >::assign_zero (const elliptic_curve< T > & e)
{
  debug_handler("base_point< T >", "assign_zero(const elliptic_curve< T > &");

  this->ec.assign(e);
  this->info = 0;
  this->is_0 = true;
  
  if (e.get_model() == elliptic_curve_flags::PROJECTIVE) {
    this->x.assign_zero();
    this->y.assign_one();
    this->z.assign_zero();
  }
}


template< class T >
void
base_point< T >::assign (const T & xx, const T & yy, 
			 const elliptic_curve< T > & e)
{
  debug_handler("base_point< T >", "assign_zero(const T &, const T &, "
		"const elliptic_curve< T > &");
  
  this->ec.assign(e);
  this->is_0 = false;
  this->x = xx;
  this->y = yy;
  
  this->z.assign(xx);
  this->z.assign_one();
  
  if (!on_curve())
    lidia_error_handler("base_point< T >::assign(x, y, e)",
			"Point not on curve.");
}



template< class T >
void
base_point< T >::assign (const T & xx, const T & yy, const T & zz,
			 const elliptic_curve< T > & e)
{
  this->ec.assign(e);
  this->x = xx;
  this->y = yy;
  this->z = zz;
  
  if (e.get_model() == elliptic_curve_flags::AFFINE || !on_curve())
    lidia_error_handler("base_point< T >::assign(x, y, z, e)",
			"Point not on curve or model affine.");
}



// for efficiency reasons

template< class T >
void
base_point< T >::assign_no_test (const T & xx, const T & yy,
				 const elliptic_curve< T > & e)
{
  this->ec.assign(e);
  this->is_0 = false;
  this->x = xx;
  this->y = yy;
  
  if (e.get_model() == elliptic_curve_flags::PROJECTIVE) 
    {
      this->z.assign(this->x);
      this->z.assign_one();
    }
}

template< class T >
void
base_point< T >::assign_no_test (const T & xx, const T & yy, const T & zz,
				 const elliptic_curve< T > & e)
{
  this->ec.assign(e);
  this->x = xx;
  this->y = yy;
  this->z = zz;
  if (e.get_model() == elliptic_curve_flags::AFFINE)
    lidia_error_handler("base_point< T >::assign(x, y, z)",
			"affine Model.");
}



template< class T >
void
base_point< T >::assign_zero ()
{
  this->is_0 = true;
  if (this->ec.get_model() == elliptic_curve_flags::PROJECTIVE) 
    {
      this->x.assign_zero();
      this->y.assign_one();
      this->z.assign_zero();
    }
}



template< class T >
void
base_point< T >::assign (const T & xx, const T & yy)
{
  this->is_0 = false;
  this->x = xx;
  this->y = yy;
  
  if (this->ec.get_model() == elliptic_curve_flags::PROJECTIVE) {
    this->z.assign(this->x);
    this->z.assign_one();
  }
  
  if (!on_curve())
    lidia_error_handler("base_point< T >::assign(x, y)",
			"Point not on curve.");
}



template< class T >
void
base_point< T >::assign (const T & xx, const T & yy, const T & zz)
{
  this->x = xx;
  this->y = yy;
  this->z = zz;
  
  if (!on_curve() || this->get_curve().get_model() == 
      elliptic_curve_flags::AFFINE)
    lidia_error_handler("base_point< T >::assign(x, y, z)",
			"Point not on curve or affine model.");
}



//
// swap
//

template< class T >
void
base_point< T >::swap (base_point< T > & P)
{
  LiDIA::swap(this->ec, P.ec);
  
  LiDIA::swap(this->x, P.x);
  LiDIA::swap(this->y, P.y);
  LiDIA::swap(this->z, P.z);
  
  bool b = P.is_0;
  P.is_0 = this->is_0;
  this->is_0 = b;
  
  int i = P.info;
  P.info = this->info;
  this->info = i;
}



//
// properties
//

template< class T >
bool
base_point< T >::on_curve() const
{
  if (this->is_0)
    return true;
  
  T l, r; 

  if (this->ec.get_model() == elliptic_curve_flags::PROJECTIVE) 
    {
      T h, h2; // experimental
      
      multiply(l, this->ec.get_a1(), this->x);
      square(h, this->z);
      multiply(r, h, this->ec.get_a3());
      add(l, l, r);
      multiply(l, l, this->z);
      add(l, l, this->y);
      multiply(l, l, this->y);
      
      multiply(r, this->ec.get_a2(), h);
      add(r, this->x, r);
      multiply(r, r, this->x);
      square(h2, h);
      multiply(h, h2, h);
      multiply(h2, this->ec.get_a4(), h2);
      add(r, r, h2);
      multiply(r, r, this->x);
      multiply(h, h, this->ec.get_a6());
      add(r, r, h);
    }
  else {
    multiply(l, this->ec.get_a1(), this->x);
    add(l, l, this->ec.get_a3());
    add(l, l, this->y);
    multiply(l, l, this->y);
    
    add(r, this->x, this->ec.get_a2());
    multiply(r, r, this->x);
    add(r, r, this->ec.get_a4());
    multiply(r, r, this->x);
    add(r, r, this->ec.get_a6());
  }
  
  subtract(l, l, r);
  return (l.is_zero());
}



//
// comparison
//

template< class T >
bool
base_point< T >::is_equal (const base_point< T > & P) const
{
  if (&P == this) {
    return true;
  }
  if (P.ec == this->ec) {
    if (P.ec.get_model() == elliptic_curve_flags::AFFINE)
      return ((P.is_0 == this->is_0) &&
	      (P.is_0 || (P.x == this->x && P.y == this->y)));
    
    else {
      // PROJECTIVE
      if (P.is_zero() && is_zero())
	return true;
      
      if (P.is_zero() || is_zero())
	return false;
      
      T x1(P.z), y1(P.z), x2(P.z), y2(P.z);

      square(x1, P.z); // x1 = (P.z)^2
      square(x2, this->z); // x2 = (z)^2
      multiply(y1, this->x, x1); // y1 = x * (P.z)^2
      multiply(y2, P.x, x2); // y2 = P.x * (z)^2
      
      if (y1 != y2)
	return false;
      
      multiply(x1, x1, P.z); // x1 = (P.z)^3
      multiply(x2, x2, this->z); // x2 = (z)^3
      multiply(y1, this->y, x1); // y1 = y * (P.z)^3
      multiply(y2, P.y, x2); // y2 = P.y * (z)^3
      
      if (y1 == y2)
	return true;
      else
	return false;
    }
  }
  return false;
}



//
// access
//

template< class T >
void
base_point< T >::make_affine_point (bool change_curve)
{
  if (this->ec.get_model() == elliptic_curve_flags::PROJECTIVE) 
    {
      if (is_zero()) 
	{
	  this->is_0 = true;
	  if (change_curve)
	    this->ec.set_model(elliptic_curve_flags::AFFINE);
	  return;
	}
    
      T h; 
      
      invert(this->z, this->z);
      square(h, this->z);
      multiply(this->x, this->x, h);
      multiply(h, h, this->z);
      multiply(this->y, this->y, h);
      this->z.assign_one();
      if (change_curve)
	this->ec.set_model(elliptic_curve_flags::AFFINE);
    }
}


template< class T >
void
base_point< T >::make_projective_point (bool change_curve)
{
  if (this->ec.get_model() == elliptic_curve_flags::AFFINE) {
    if (is_zero()) {
      this->x.assign_zero();
      this->y.assign_one();
      this->z.assign_zero();
    }
    else 
      this->z.assign_one();

    if (change_curve)
      this->ec.set_model(elliptic_curve_flags::PROJECTIVE);
  }
}



//
// arithmetic
//

//
// ::negate
//

template< class T >
void
negate (base_point< T > & r, const base_point< T > & x)
{
  if (x.is_zero())
    r = x;
  else {
    r.ec = x.ec;
    x.ec.get_point_operations().negate(r, x);
  }
}



//
// ::add
//

template< class T >
void
add (base_point< T > & R, const base_point< T > & P,
     const base_point< T > & Q)
{
  if (P.ec != Q.ec)
    lidia_error_handler("base_point< T >::add(R, P, Q)",
			"Different elliptic curves.");
  else
    R.ec = P.ec;

  if (Q.is_zero()) 
    {
      R.assign(P);
      return;
    }
  
  if (P.is_zero()) {
    R.assign(Q);
    return;
  }
  
  if (P.ec.get_model() == elliptic_curve_flags::PROJECTIVE)
    P.ec.get_point_operations().add(R, P, Q);
  else 
    {
      if (P.is_negative_of(Q)) {
	R.assign_zero();
	return;
      }
    
      if (P == Q)
	P.ec.get_point_operations().mult_by_2(R, P);
      else
	P.ec.get_point_operations().add(R, P, Q);
    }
}



//
// ::subtract
//

template< class T >
void
subtract (base_point< T > & R, const base_point< T > & P,
	  const base_point< T > & Q)
{
  if (P.ec != Q.ec)
    lidia_error_handler("base_point", "subtract::different elliptic curves");
  else
    R.ec = P.ec;
  
  if (P.ec.get_model() == elliptic_curve_flags::PROJECTIVE) {
    base_point< T > H;
    P.ec.get_point_operations().negate(H, Q);
    P.ec.get_point_operations().add(R, P, H);
    return;
  }
  
  if (P == Q) {
    R.assign_zero();
    return;
  }
  
  if (Q.is_zero()) {
    R.assign(P);
    return;
  }
  
  if (P.is_zero()) {
    P.ec.get_point_operations().negate(R, Q);
    return;
  }
  else {
    base_point< T > H;
    
    P.ec.get_point_operations().negate(H, Q);
    if (P == H)
      P.ec.get_point_operations().mult_by_2(R, P);
    else
      P.ec.get_point_operations().add(R, P, H);
  }
}



//
// ::multiply_by_2
//

template< class T >
void
multiply_by_2 (base_point< T > & R, const base_point< T > & P)
{
  if (P.is_zero()) {
    R.assign_zero();
    return;
  }
  P.ec.get_point_operations().mult_by_2(R, P);
}



//
// ::multiply
//


template< class T >
void
multiply (base_point< T > & R, const bigint & mm,
	  const base_point< T > & P)
{
  bigint m(mm), m3(3*m);
  int i;
  base_point< T > Q(P);
  
  if (m.is_zero()) {
    R.assign_zero();
    return;
  }
  
  if (m.is_negative()) {
    m.negate();
    m3.negate();
    negate(Q, Q);
  }
  
  if (m.is_one()) {
    R.assign(Q);
    return;
  }
  
  R.assign_zero(Q.ec);
  
  for (i = m3.bit_length() - 1; i > 0; i--) {
    multiply_by_2(R, R);
    if (m.bit(i) == 0 && m3.bit(i) == 1)
      add(R, R, Q);
    else
      if (m.bit(i) == 1 && m3.bit(i) == 0)
	subtract(R, R, Q);
  }
}



template< class T >
base_point< T >
base_point< T >::twice () const
{
  base_point< T > P(this->ec);
  multiply_by_2(P, (*this));
  return P;
}



//
// input / output
//

template< class T >
void
base_point< T >::write (std::ostream & out) const
{
  if (get_curve().get_model() == elliptic_curve_flags::PROJECTIVE)
    {
      if (is_zero())
	out << "(0 : 1 : 0)";
      else 
	out << "(" << this->x << " : " << this->y << " : " << this->z << ")";
    }
  else 
    {
     if (is_zero())
       out << "O";
     else
       out << "(" << this->x << " , " << this->y << ")";
    }
}



template< class T >
void
base_point< T >::read (std::istream & in)
{
  char c, c2, c3; // to eat commas and detect [;
  // `(' or '[' flags input format (x, y)
  // separators and terminators must then be present
  // (any nonnumeric will do after the first { or [)
  // else assumes x y separated by whitespace
  in >> c; // This eats any whitespace first, and c is the next char input
  
  // initialize field of coefficients
  this->x = this->ec.get_a4();
  this->y = this->x;
  
  switch (c) {
  case '[':
    in >> c;
    if (c != ']') {
      in.putback(c);
      if (this->get_curve().get_model() == elliptic_curve_flags::AFFINE) {
	in >> this->x >> c >> this->y >> c2;
	if (c != ',' || c2 != ']')
	  lidia_error_handler("base_point", "read::wrong input format");
	this->is_0 = false;
      }
      else {
	in >> this->x >> c >> this->y >> c2 >> this->z >> c3;
	if (c != ':' || c2 != ':' || c3 != ']')
	  lidia_error_handler("base_point", "read::wrong input format");
      }
    }
    else
      this->is_0 = true;
    break;
    
  case '(':
    in >> c;
    if (c != ')') {
      in.putback(c);
      if (this->get_curve().get_model() == elliptic_curve_flags::AFFINE) {
	in >> this->x >> c >> this->y >> c2;
	if (c != ',' || c2 != ')')
	  lidia_error_handler("base_point", "read::wrong input format");
	this->is_0 = false;
      }
      else {
	in >> this->x >> c >> this->y >> c2 >> this->z >> c3;
	if (c != ':' || c2 != ':' || c3 != ')')
	  lidia_error_handler("base_point", "read::wrong input format");
      }
    }
    else
      this->is_0 = true;
    break;
    
  case 'o':
  case 'O':
    if (this->get_curve().get_model() == elliptic_curve_flags::AFFINE)
      this->is_0 = true;
    else {
      this->x.assign_zero();
      this->y.assign_one();
      this->z.assign_zero();
    }
    break;
    
  default:
    in.putback(c);
    if (this->get_curve().get_model() == elliptic_curve_flags::AFFINE)
      in >> this->x >> this->y;
    else
      in >> this->x >> this->y >> this->z;
    
  }
  
  if (!on_curve())
    lidia_error_handler("base_point", "read::point not on curve");
}



#ifdef LIDIA_NAMESPACE
# ifndef IN_NAMESPACE_LIDIA
}	// end of namespace LiDIA
# endif
#endif



#endif	// LIDIA_BASE_POINT_CC_GUARD_
