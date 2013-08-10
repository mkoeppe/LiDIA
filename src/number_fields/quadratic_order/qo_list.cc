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
//	Author	: Michael Jacobson (MJ)
//	Changes	: See CVS log
//
//==============================================================================================


#ifdef HAVE_CONFIG_H
# include	"config.h"
#endif
#include	"LiDIA/number_fields/qo_list.h"
#include	"LiDIA/quadratic_order.h"



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



qo_node::qo_node()
{
	debug_handler ("qo_node", "qo_node()");
	rfc = 1;
	dynamic = false;
	next = NULL;
}



qo_node::qo_node(quadratic_order *QO2, bool is_dyn)
{
	debug_handler ("qo_node", "qo_node(quadratic_order*, bool)");
	QO = QO2;
	rfc = 1;
	dynamic = is_dyn;
	next = NULL;
}



qo_node::qo_node(quadratic_order *QO2, qo_node *n, bool is_dyn)
{
	debug_handler ("qo_node", "qo_node(quadratic_order*, qo_node*, bool)");
	QO = QO2;
	next = n;
	rfc = 1;
	dynamic = is_dyn;
}



qo_node::~qo_node()
{
	debug_handler ("qo_node", "~qo_node()");
	if (dynamic) {
	    quadratic_order* tmp = QO;
	    QO = 0;
	    delete tmp;
	}
}



unsigned int
qo_node::get_ref() const
{
	debug_handler ("qo_node", "get_ref()");
	return rfc;
}



unsigned int
qo_node::inc_ref()
{
	debug_handler ("qo_node", "inc_ref()");
	rfc++;
	return rfc;
}



unsigned int
qo_node::dec_ref()
{
	debug_handler ("qo_node", "dec_ref()");
	rfc--;
	return rfc;
}



const quadratic_order *
qo_node::get_qo() const
{
	debug_handler ("qo_node", "const quadratic_order* get_qo()");
	return QO;
}



quadratic_order *
qo_node::get_qo()
{
	debug_handler ("qo_node", "quadratic_order* get_qo()");
	return QO;
}



void
qo_node::set_qo(quadratic_order *QO2)
{
	debug_handler ("qo_node", "set_qo(const quadratic_order*)");
	QO = QO2;
}



void
qo_node::set_succ(qo_node *e)
{
	debug_handler ("qo_node", "set_succ(qo_node*)");
	next = e;
}



qo_node *
qo_node::succ()
{
	debug_handler ("qo_node", "succ()");
	return next;
}



void swap (qo_node* &x, qo_node* &y)
{
	debug_handler ("qo_node", "swap(qo_node*, qo_node*)");
	qo_node* z = x;
	x = y;
	y = z;
}



qo_list::qo_list()
{
	debug_handler ("qo_list", "entering qo_list()");
	head = NULL;
	last_qo = &quadratic_order::zero_QO;
	debug_handler ("qo_list", "leaving qo_list()");
}



qo_list::~qo_list()
{
	debug_handler("qo_list", "~qo_list()");
	qo_node *next;

	while (head != NULL) {
		next = head->succ();
		delete head;
		head = next;
	}
}



bool
qo_list::is_empty() const
{
	debug_handler ("qo_list", "bool is_empty()");
	return (head == NULL);
}



const qo_node *
qo_list::first() const
{
	debug_handler ("qo_list", "const qo_node* first()");
	return head;
}



qo_node *
qo_list::first()
{
	debug_handler ("qo_list", "qo_node* first()");
	return head;
}



const quadratic_order *
qo_list::last() const
{
	debug_handler ("qo_list", "const quadratic_order* last()");
	return last_qo;
}



quadratic_order *
qo_list::last()
{
	debug_handler ("qo_list", "quadratic_order* last()");
	return last_qo;
}



void
qo_list::set_last(quadratic_order *QO)
{
	debug_handler("qo_list", "set_last(quadratic_order*)");
	last_qo = QO;
}



void
qo_list::print(std::ostream & out) const
{
	debug_handler("qo_list", "print(std::ostream&)");

	if (head == NULL)
		std::cout << "qo_list::EMPTY LIST\n";
	else {
		qo_node *e = head;

		while (e != NULL) {
			if (e->get_qo())
				out << " < " << e->get_qo()->discriminant() << ", " << e->get_ref() << " >";
			else
				out << " < NULL, " << e->get_ref() << " >";
			e = e->succ();
		}
		out << "\n";
	}
}



//
// qo_list::nullify()
//
// Task:
//      if there is a node in the list with QO = QO2, set the pointer to the
//      quadratic order in this node to NULL.
//

void
qo_list::nullify(quadratic_order *QO2)
{
	debug_handler("qo_list", "nullify(quadratic_order*)");

	qo_node *e = head;

	while (e != NULL) {
		// found ?
		if (e->get_qo() == QO2) {
			if (clear(e) == 1)
				e->set_qo(&quadratic_order::zero_QO);
			break;
		}
		else
			e = e->succ();
	}

	if ((last_qo == QO2) || (!last_qo))
		if (head)
			last_qo = head->get_qo();
		else
			last_qo = &quadratic_order::zero_QO;
}



//
// qo_list::insert()
//
// Task:
//      if there is a node in the list with QO = QO2, then increase the
//      reference counter of this node by one and return a pointer to it.
//      Otherwise, create a new node in the list.
//

qo_node *
qo_list::insert(quadratic_order & QO2, bool is_dyn)
{
	debug_handler("qo_list", "insert(quadratic_order&, bool)");

	qo_node *e = head;

	while (e != NULL) {
		// found ?
		if (e->get_qo() == &QO2) {
			e->inc_ref();
			return e;
		}
		else
			e = e->succ();
	}

	// quadratic_order is not in the list yet => create it
	e = new qo_node(&QO2, head, is_dyn);
	head = e;
	return e;
}



//
// qo_list::clear()
//
// Task:
//      decrease the reference counter of elmt by one, and if the counter
//      is now zero, delete the node from the list.
//
// return value:
//      0  the qo_node to be deleted doesnt exist
//      1  the rfc of the qo_node to be deleted is decreased
//      2  the qo_node is deleted
//
// WARNING: To avoid side effects, elmnt must point
//      to an element of the list *this or should be NULL, because,
//      if the rfc of elmnt is greater than one, the
//      rfc is decreased. In this case, it is not verified,
//      whether elmnt points to an elemnt of the list.
//

int
qo_list::clear(qo_node * elmnt)
{
	debug_handler("qo_list", "clear(qo_node*)");

	int rc;

	if (head == NULL || elmnt == NULL)
		rc = 0;
	else {
		if (elmnt->get_ref() > 1) {
			// just decrease the ref. counter
			elmnt->dec_ref();
			rc = 1;
		}
		else {
			// otherwise, remove the element from the list
			qo_node *e, *pred;
			rc = 2;

			if (elmnt == head) {
				// elmnt is the head of the list
				e = head->succ();
				if (head->get_qo() == last_qo)
					last_qo = &quadratic_order::zero_QO;
				delete head;
				head = e;
			}
			else {
				// search for the predecessor of elmnt in the list
				pred = head;
				e = head->succ();
				while (e != elmnt && e != NULL) {
					pred = e;
					e = e->succ();
				}
				// found ?
				if (e == elmnt) {
					pred->set_succ(e->succ());
					if (e->get_qo() == last_qo)
						last_qo = &quadratic_order::zero_QO;
					delete e;
					e = NULL;
				}
				else
					rc = 0;
			}
		}
	}

	return (rc);
}



//
// qo_list::add_to_list()
//
// Task:
//      adds the new order to the list of quadratic orders at the given node.
//

qo_node *
qo_list::add_to_list(qo_node *old, quadratic_order & newQO)
{
	debug_handler("qo_list", "add_to_list");

	if (old)
		if (old->get_qo() == &newQO) {
			old->inc_ref();
			return old;
		}
		else
			clear(old);

	return insert(newQO, 0);
}



//
// qo_list::add_to_list()
//
// Author: Markus Maurer
//
// Task:
//    For efficiency reasons, this function works on pointers to
//    nodes. It adds the new order to the list of quadratic orders and
//    deletes the old one. It returns a pointer to the node of the new
//    order.
//
// Conditions:
//    old must be pointer to a node of the list or NULL.
//    new must be a pointer to a node of the list.
//
// WARNING:
//    new is not checked for NULL.
//

qo_node*
qo_list::add_to_list(qo_node *old_node, qo_node *new_node)
{
	debug_handler("qo_list", "add_to_list(qo_node*, qo_node*)");

	// Decrease reference for old_node and
	// increase reference for new_node.
	//
	if (old_node != new_node) {
		if (old_node)
			clear(old_node);
	}

	new_node->inc_ref();

	return new_node;
}



//
// qo_list::add_to_list(bigint)
//
// Task:
//      first, checks whether a quadratic order with discriminant Delta is
//      already in the list.  If so, a pointer to this quadratic order is
//      returned.  Otherwise, a quadratic order with discriminant Delta is
//      dynamically allocated, added to the list, and a pointer to it is
//      returned.
//

quadratic_order *
qo_list::add_to_list(const bigint & Delta)
{
	debug_handler("qo_list", "add_to_list(bigint)");

	qo_node *e = head;

	while (e != NULL) {
		// found ?
		if (e->get_qo()->discriminant() == Delta) {
			return e->get_qo();
		}
		else
			e = e->succ();
	}

	// quadratic_order is not in the list yet => create it
	quadratic_order *QO2 = new quadratic_order(Delta);

	e = new qo_node(QO2, head, 1);
	head = e;
	return QO2;
}



//
// qo_list::add_last_to_list()
//
// Task:
//      adds the last order to the list
//

qo_node *
qo_list::add_last_to_list(qo_node *old)
{
	debug_handler("qo_list", "add_last_to_list");

	if (!last_qo) {
		lidia_error_handler("qo_list", "add_last_to_list - last referenced quadratic order is NULL");
		return NULL;
	}

	if (old)
		clear(old);
	return insert(*last_qo, 0);
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
