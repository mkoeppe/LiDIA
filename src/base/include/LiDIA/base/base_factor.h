// -*- C++ -*-
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
//	Author	:
//	Changes	: See CVS log
//
//==============================================================================================


#ifndef LIDIA_BASE_FACTOR_H_GUARD_
#define LIDIA_BASE_FACTOR_H_GUARD_


#ifndef LIDIA_INTERFACE_LIB_H_GUARD_
# include	"LiDIA/base/interface_lib.h"
#endif



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif



//***************************************************************
//				class decomposable_object
//
// handles information about primality
// base class for
//	- base_factor< T >
//	- factorization< T >
//***************************************************************

class decomposable_object
{
public :
	enum decomp_state {
		not_prime = 0,
		unknown = 1,
		prime = 2
	};
	//if state = prime or state = not_prime, we know that this object is prime
	//or composite, respectively.
	//if state = unknown, we haven't tested it for primality, yet.

private:
	decomp_state S;

protected:

	decomposable_object() : S(unknown) {}
	decomposable_object(decomp_state s) : S(s) {}
	decomposable_object(const decomposable_object & dco) : S(dco.S) {}
	~decomposable_object() {}

	decomposable_object & operator = (const decomposable_object & dco)
	{
		S = dco.S;
		return *this;
	}

	void swap (decomposable_object & dco)
	{
		decomp_state	tmp = S;
		S = dco.S;
		dco.S = tmp;
	}

	void share_state(const decomposable_object &a)
	{
		//this function is called whenever two decomposable_object are discovered to be equal
		//if this->S is 'unknown', the value of a.S is copied
		//an error is raised if this->S = 'prime' and a.S = 'not_prime' or vice versa
		if (S == unknown)
			S = a.S;
		else if (S != a.S && a.S != unknown)
			lidia_error_handler("decomposable_object", "share_state(decomposable_object&)::conflicting states");
	}

	void concat_state(const decomposable_object &a, const decomposable_object &b)
	{
		//this->S = state of product of a and b
		if (b.S == not_prime)
			S = not_prime;
		else if (b.S == unknown && S == prime)
			S = unknown;
		else S = a.S;
	}


public:

#ifndef HEADBANGER
	void set_prime_flag(decomp_state new_S)
	{
		S = new_S;
	}

	decomp_state prime_flag() const
	{
		return S;
	}
#endif	// HEADBANGER

};



//****************************************************************
//				class base_factor< T >
//		only a base class for single_factor< T >
//****************************************************************

// forward declarations needed for template friends
template< class T >
class base_factor;

template<class T >
inline void multiply(base_factor< T > & c, const base_factor< T > & a,
		     const base_factor< T > & b);
template<class T >
inline void divide(base_factor< T > & c, const base_factor< T > & a,
		   const base_factor< T > & b);
template<class T >
base_factor< T > operator*(const base_factor< T > & a,
			   const base_factor< T > & b);
template<class T >
base_factor< T > operator/(const base_factor< T > & a,
			   const base_factor< T > & b);


template< class T >
class base_factor: public decomposable_object
{

protected:
	T rep;

protected:
	//
	// constructors, destructor
	//

	base_factor();
	base_factor(const T &);
	base_factor(const base_factor< T > &);
	~base_factor() {};

private:
	base_factor & operator = (const base_factor< T > &); // disallow
public:


	//
	// accessors
	//

	const T& base() const;
	T & base ();


	//
	// arithmetic operations
	//

#ifndef HEADBANGER
	friend
	void multiply< T >(base_factor< T > & c, const base_factor< T > & a,
				  const base_factor< T > & b);
	//c = a*b

	friend void divide< T >(base_factor< T > & c, const base_factor< T > & a,
				const base_factor< T > & b);
	//c = a/b
#endif	// HEADBANGER

#if 0
	friend void gcd< T >(base_factor< T > & c, const base_factor< T > & a,
			     const base_factor< T > & b);
		// c = gcd(a, b)
#endif



	void swap (base_factor & a);

	//
	// I/O
	//

	void read(std::istream & in);
	void write(std::ostream & out) const;


	//
	// c'tors and d'tor is protected, hence these must be declared friend
	//

	friend base_factor< T >
	operator * < T >(const base_factor< T > & a,
			 const base_factor< T > & b);

	friend base_factor< T >
	operator / < T >(const base_factor< T > & a,
			 const base_factor< T > & b);
};



//
// constructors, destructor
//

template< class T >
inline
base_factor< T >::base_factor()
	: decomposable_object(not_prime)
{
}



template< class T >
inline
base_factor< T >::base_factor(const T & x)
	: decomposable_object(unknown),
	  rep(x)
{
}



template< class T >
inline
base_factor< T >::base_factor(const base_factor< T > & x)
	: decomposable_object(x),
	  rep(x.rep)
{
}



//
// accessors
//


template< class T >
inline const T &
base_factor< T >::base() const
{
	return this->rep;
}



template< class T >
inline T &
base_factor< T >::base ()
{
	return this->rep;
}



//
//
//

template< class T >
inline void
multiply(base_factor< T > & c,
	 const base_factor< T > & a,
	 const base_factor< T > & b)
{
	multiply(c.rep, a.rep, b.rep);
	c.set_prime_flag(decomposable_object::unknown);
}



template< class T >
inline void
divide(base_factor< T > & c,
       const base_factor< T > & a,
       const base_factor< T > & b)
// ********************************************************************
// *WARNING : the software cannot check if 'b' is a divisor of 'a' !!!*
// ********************************************************************
{
	divide(c.rep, a.rep, b.rep);
	c.set_prime_flag(decomposable_object::unknown);
}



template< class T >
inline base_factor< T >
operator * (const base_factor< T > & a,
	    const base_factor< T > & b)
{
	base_factor< T > tmp;

	multiply(tmp, a, b);
	return tmp;
}



template< class T >
inline base_factor< T >
operator / (const base_factor< T > & a,
	    const base_factor< T > & b)
{
	base_factor< T > tmp;

	divide(tmp, a, b);
	return tmp;
}



template< class T >
inline bool
operator < (const base_factor < T > & a,
	   const base_factor< T > & b)
{
	return (a.base() < b.base());
}



template< class T >
inline bool
operator > (const base_factor< T > & a,
	    const base_factor< T > & b)
{
	return !(a.base() <= b.base());
}



template< class T >
inline bool
operator <= (const base_factor< T > & a,
	     const base_factor< T > & b)
{
	return (a.base() <= b.base());
}



template< class T >
inline bool
operator >= (const base_factor< T > & a,
	     const base_factor< T > & b)
{
	return !(a.base() < b.base());
}



template< class T >
inline bool
operator == (const base_factor< T > & a,
	     const base_factor< T > & b)
{
	return (a.base() == b.base());
}



template< class T >
inline bool
operator != (const base_factor< T > & a,
	     const base_factor< T > & b)
{
	return (a.base() != b.base());
}



template< class T >
inline void
base_factor< T >::read (std::istream & in)
{
	in >> this->rep;
	set_prime_flag(decomposable_object::unknown);
}



template< class T >
inline void
base_factor< T >::write (std::ostream & out) const
{
	out << this->rep;
}



template< class T >
inline void
base_factor< T >::swap (base_factor< T > & a)
{
	decomposable_object::swap(a);
	LiDIA::swap(this->rep, a.rep);
}



template< class T >
inline std::istream &
operator >> (std::istream &in, base_factor< T > &f)
{
	f.read(in);
	return in;
}



template< class T >
inline std::ostream &
operator << (std::ostream &out, const base_factor< T > &f)
{
	f.write(out);
	return out;
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif



#ifdef LIDIA_INCLUDE_CC
# include	"LiDIA/base/base_factor.cc"
#endif



#endif	// LIDIA_BASE_FACTOR_H_GUARD_
