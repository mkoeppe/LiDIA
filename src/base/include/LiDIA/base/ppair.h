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
//	Author	: Thomas Papanikolaou (TP), Thomas Pfahler (TPf)
//	Changes	: See CVS log
//
//==============================================================================================


#ifndef LIDIA_PPAIR_H_GUARD_
#define LIDIA_PPAIR_H_GUARD_



#ifndef LIDIA_LIDIA_H_GUARD_
#include	"LiDIA/LiDIA.h"
#endif



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif



#ifndef HEADBANGER

template< class T1, class T2 > class ppair
{
	T1 *l;
	T2 r;

public:

	ppair();
	ppair(T1 & pl, T2 & pr);
	ppair(const ppair< T1, T2 > & p);

	~ppair();


	void assign (const ppair< T1, T2 > & p);
	ppair< T1, T2 > & operator = (const ppair< T1, T2 > & p);


	const T1 & left() const;
	const T2 & right() const;
	T1 & left();
	T2 & right();


	void read(std::istream & in);
	void write(std::ostream & out) const;


	int compare(const ppair< T1, T2 > &p) const;

	bool operator == (const ppair< T1, T2 > & p) const;
	bool operator != (const ppair< T1, T2 > & p) const;
	bool operator > (const ppair< T1, T2 > & p) const;
	bool operator >= (const ppair< T1, T2 > & p) const;
	bool operator < (const ppair < T1, T2 > & p) const;
	bool operator <= (const ppair< T1, T2 > & p) const;

	void swap (ppair< T1, T2 > & p);
};



//
// c'tors and d'tor
//

template< class T1, class T2 >
inline ppair< T1, T2 >::ppair ()
{
	l = new T1;
	if (!l)
		lidia_error_handler("", "out of memory");
	// *l = 0;
	r = 0;
}



template< class T1, class T2 >
inline ppair< T1, T2 >::ppair (T1 &pl, T2 &pr)
{
	l = new T1;
	if (!l)
		lidia_error_handler("", "out of memory");
	*l = pl;
	r = pr;
}



template< class T1, class T2 >
inline ppair< T1, T2 >::ppair (const ppair< T1, T2 > &p)
{
	l = new T1;
	if (!l)
		lidia_error_handler("", "out of memory");
	*l = *p.l;
	r = p.r;
}



template< class T1, class T2 >
inline
ppair< T1, T2 >::~ppair ()
{
	delete l;
}



//
// accessors
//

template< class T1, class T2 >
inline const T1 &
ppair< T1, T2 >::left() const
{
	return *l;
}



template< class T1, class T2 >
inline const T2 &
ppair< T1, T2 >::right() const
{
	return r;
}



template< class T1, class T2 >
inline T1 &
ppair< T1, T2 >::left()
{
	return *l;
}



template< class T1, class T2 >
inline T2 &
ppair< T1, T2 >::right()
{
	return r;
}



//
// assigners
//

template< class T1, class T2 >
inline void
ppair< T1, T2 >::assign (const ppair< T1, T2 > & p)
{
	if (&p != this) {
		*l = *p.l;
		r = p.r;
	}
}



template< class T1, class T2 >
inline ppair< T1, T2 > &
ppair< T1, T2 >::operator = (const ppair< T1, T2 > & p)
{
	assign(p);
	return *this;
}



//
// comparators
//

template< class T1, class T2 >
inline bool
ppair< T1, T2 >::operator == (const ppair< T1, T2 > &p) const
{
  return ((&p == this) || ((*l == *p.l) && (r == p.r)));
}



template< class T1, class T2 >
inline bool
ppair< T1, T2 >::operator != (const ppair< T1, T2 > &p) const
{
	return ((*l != *p.l) || (r != p.r));
}



template< class T1, class T2 >
inline bool
ppair< T1, T2 >::operator > (const ppair< T1, T2 > &p) const
{
	return ((*l > *p.l) || ((*l == *p.l) && (r > p.r)));
}



template< class T1, class T2 >
inline bool
ppair< T1, T2 >::operator >= (const ppair< T1, T2 > &p) const
{
	return (*l >= *p.l);
}



template< class T1, class T2 >
inline bool
ppair< T1, T2 >::operator < (const ppair < T1, T2 > &p) const
{
	return ((*l < *p.l) || ((*l == *p.l) && (r < p.r)));
}



template< class T1, class T2 >
inline bool
ppair< T1, T2 >::operator <= (const ppair< T1, T2 > &p) const
{
	return (*l <= *p.l);
}



//
// I/O
//

template< class T1, class T2 >
inline std::istream &
operator >> (std::istream &in, ppair< T1, T2 > &p)
{
	p.read(in);
	return in;
}



template< class T1, class T2 >
inline std::ostream &
operator << (std::ostream &out, const ppair< T1, T2 > &p)
{
	p.write(out);
	return out;
}



template< class T1, class T2 >
inline int
compare (const ppair< T1, T2 > &p, const ppair< T1, T2 > &q)
{
	return p.compare(q);
}



template< class T1, class T2 >
inline void
swap (ppair< T1, T2 > &p, ppair< T1, T2 > &q)
{
	p.swap(q);
}



#endif	// HEADBANGER



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif



#ifdef LIDIA_INCLUDE_CC
# include	"LiDIA/ppair.cc"
#endif



#endif	// LIDIA_PPAIR_H_GUARD_
