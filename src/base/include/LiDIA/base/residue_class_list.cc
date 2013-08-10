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


#ifndef LIDIA_RESIDUE_CLASS_LIST_CC_GUARD_
#define LIDIA_RESIDUE_CLASS_LIST_CC_GUARD_


#include	"LiDIA/base/residue_class_list.h"



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif



template< class T >
residue_class_list< T >::~residue_class_list()
{
	debug_handler ("residue_class",
		       "~residue_class_list()");

	residue_class< T > *next;

	while (head != NULL) {
		next = head->get_succ();
		delete head;
		head = next;
	}
}



template< class T >
residue_class< T > *
residue_class_list< T >::insert(const T & data)
{
	debug_handler ("residue_class_list",
		       "insert(const T&)");

	//
	// searches for an existing residue_class with mod = data
	//

	residue_class< T > *e = head;

	while (e != NULL) {
				// found ?
		if (e->get_mod() == data) {
			e->inc_ref();
			return e;
		}
		else
			e = e->get_succ();
	}

	//
	// residue_class does not exist yet =  > create it
	//

	e = new residue_class< T > (data, head);
	head = e;
	return e;
}



template< class T >
int
residue_class_list< T >::clear(residue_class< T > * elmnt)
{
	debug_handler ("residue_class_list",
		       "clear(residue_class< T > *)");

	// return value: 0  the residue_class to be deleted doesnt exist
	//               1  the rfc of the residue_class to be deleted is decreased
	//		2  the residue_class is deleted

	// WARNING: To avoid side effects, elmnt must point
	// to an element of the list *this or should be NULL, because,
	// if the rfc of elmnt is greater than one, the
	// rfc is decreased. In this case, it is not verified,
	// whether elmnt points to an elemnt of the list.

	int rc;

	if (head == NULL || elmnt == NULL)
		rc = 0;
	else {
				//
				// just decrease the ref. counter
				//

		if (elmnt->get_ref() > 1) {
			elmnt->dec_ref();
			rc = 1;
		}

				//
				// otherwise, remove the element from the list
				//

		else {
			residue_class< T > *e, *pred;
			rc = 2;

				//
				// elmnt is the head of the list
				//

			if (elmnt == head) {
				e = head->get_succ();
				delete head;
				head = e;
			}

				//
				// search for the predecessor of elmnt in the list
				//

			else {
				pred = head;
				e = head->get_succ();

				while (e != elmnt && e != NULL) {
					pred = e;
					e = e->get_succ();
				}

				// found ?

				if (e == elmnt) {
					pred->set_succ(e->get_succ());
					delete e;
				}
				else
					rc = 0;
			}
		}
	}

	return (rc);
}



template< class T >
bool
residue_class_list< T >::print(std::ostream & out) const
{
	debug_handler ("residue_class_list",
		       "print(std::ostream&) const");

	// return value: false:  list is empty -> nothing to print
	//               true :  all list elements have been printed

	if (head == NULL)
		return false;
	else {
		residue_class< T > *e = head;

		while (e != NULL) {
			out << "< " << e->get_mod() << ", " << e->get_ref() << " >";
			e = e->get_succ();
		}
		out << "\n";
		return true;
	}
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif



#endif	// LIDIA_RESIDUE_CLASS_LIST_CC_GUARD_
