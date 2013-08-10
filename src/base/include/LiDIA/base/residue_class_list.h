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
//	Author	: Oliver Morsch, Thomas Papanikolaou
//	Changes	: See CVS log
//
//==============================================================================================


#ifndef LIDIA_RESIDUE_CLASS_LIST_H_GUARD_
#define LIDIA_RESIDUE_CLASS_LIST_H_GUARD_



#ifndef LIDIA_LIDIA_H_GUARD_
# include	"LiDIA/LiDIA.h"
#endif



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif



template< class T >
class residue_class {
private:

	T mod; // a bigmods residue_class
	unsigned int rfc; // reference counter (how often is this class in use?)
	residue_class< T > * next; // pointer to the successor

public:

	//
	// c'tors and d'tor
	//

	residue_class();
	residue_class(const T & data);
	residue_class(const T & data, residue_class< T > * n);
	~residue_class();



private:

	// inhibit
	residue_class(const residue_class &);
	residue_class & operator = (const residue_class< T > &);



public:

	//
	// reference counter accessors and modifiers
	//

	unsigned int get_ref() const;
	unsigned int inc_ref();
	unsigned int dec_ref();



	//
	// accessors and modifiers
	//

	const T & get_mod() const;
	T & get_mod();
	void set_mod(const T & e);
	void set_succ(residue_class< T > * e);
	residue_class< T > * get_succ();

};



template< class T >
inline
residue_class< T >::residue_class()
	: mod(),
	  rfc(1),
	  next(NULL)
{
	debug_handler ("residue_class",
		       "residue_class()");
}



template< class T >
inline
residue_class< T >::residue_class(const T & data)
	: mod(data),
	  rfc(1),
	  next(NULL)
{
	debug_handler ("residue_class",
			       "residue_class(const T&)");
}



template< class T >
inline
residue_class< T >::residue_class(const T & data, residue_class< T > * n)
	: mod(data),
	  rfc(1),
	  next(n)
{
	debug_handler ("residue_class",
		       "residue_class(const T&, residue_class< T > *)");
}



template< class T >
inline
residue_class< T >::~residue_class()
{
	debug_handler ("residue_class",
		       "~residue_class()");
}



//
// reference counter accessors and modifiers
//

template< class T >
inline unsigned int
residue_class< T >::get_ref() const
{
	debug_handler ("residue_class",
		       "get_ref() const");
	return rfc;
}



template< class T >
inline unsigned int
residue_class< T >::inc_ref()
{
	debug_handler ("residue_class",
		       "inc_ref()");
	rfc++;
	return rfc;
}



template< class T >
inline unsigned int
residue_class< T >::dec_ref()
{
	debug_handler ("residue_class",
		       "dec_ref()");
	rfc--;
	return rfc;
}



//
// accessors and modifiers
//

template< class T >
inline const T &
residue_class< T >::get_mod() const
{
	debug_handler ("residue_class",
		       "const T& get_mod()");
	return mod;
}



template< class T >
inline T &
residue_class< T >::get_mod()
{
	debug_handler ("residue_class",
		       "T& get_mod()");
	return mod;
}



template< class T >
inline void
residue_class< T >::set_mod(const T & e)
{
	debug_handler ("residue_class",
		       "set_mod(const T&)");
	mod = e;
}



template< class T >
inline void
residue_class< T >::set_succ(residue_class< T > *e)
{
	debug_handler ("residue_class",
		       "set_succ(residue_class< T > *)");
	next = e;
}



template< class T >
inline residue_class< T > *
residue_class< T >::get_succ()
{
	debug_handler ("residue_class",
		       "succ()");
	return next;
}







template< class T >
class residue_class_list
{
protected:

	residue_class< T > * head;

public:

	//
	// c'tor and d'tor
	//

	residue_class_list();
	~residue_class_list();


private:
	// inhibit
	residue_class_list(const residue_class_list< T > &);
	residue_class_list< T > & operator = (const residue_class_list< T > &);

		

public:

	//
	// *** basic functions ***
	//

	const residue_class< T > * first() const;
	residue_class< T > * first();


	bool print(std::ostream & out = std::cout) const;


	residue_class< T > * insert(const T & data);


	residue_class< T > * set_to(residue_class< T > * a);


	int clear(residue_class< T > * elmnt);

};



//
// c'tor and d'tor
//

template< class T >
inline
residue_class_list< T >::residue_class_list ()
	: head(NULL)
{
	debug_handler ("residue_class_list",
		       "residue_class_list()");
}



template< class T >
inline const residue_class< T > *
residue_class_list< T >::first() const
{
	debug_handler ("residue_class_list",
		       "const residue_class< T > * first() const");
	return head;
}



template< class T >
inline residue_class< T > *
residue_class_list< T >::first()
{
	debug_handler ("residue_class_list",
		       "residue_class< T > * first()");
	return head;
}



template< class T >
inline residue_class< T > *
residue_class_list< T >::set_to(residue_class< T > * a)
{
	debug_handler ("residue_class_list",
		       "set_to(residue_class< T > *)");

	// WARNING: For efficiency reasons, the user has to take
	// the responsibility, that the pointer a really points to
	// an element of the list *this.

	a->inc_ref();
	return a;
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif



#ifdef LIDIA_INCLUDE_CC
# include	"LiDIA/base/residue_class_list.cc"
#endif



#endif  // LIDIA_RESIDUE_CLASS_LIST_H_GUARD_
