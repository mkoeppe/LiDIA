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
//	Author	: Michael Jacobson (MJ)
//	Changes	: See CVS log
//
//==============================================================================================


#ifndef LIDIA_QO_LIST_H_GUARD_
#define LIDIA_QO_LIST_H_GUARD_


#ifndef HEADBANGER


#ifndef LIDIA_BIGINT_H_GUARD_
# include	"LiDIA/bigint.h"
#endif



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif



class quadratic_order;

//
// Class: qo_node
//
// This class represents an element in the list of currently defined
//    quadratic orders.  Each instance of qo_node contains the following:
//
//    *QO - pointer to a quadratic order
//    rfc - number of references to this quadratic order
//    dynamic - true if this quadratic order was created dynamically
//    *next - successor in the list
//

class qo_node
{
private:

	quadratic_order *QO; // pointer to a quadratic order
	unsigned int rfc; // reference counter
	bool dynamic; // true if the order was created dynamically
	qo_node *next; // pointer to the successor

public:

	//
	// constructors and destructor
	//

	qo_node();
	qo_node(quadratic_order *QO2, bool is_dyn);
	qo_node(quadratic_order *QO2, qo_node *n, bool is_dyn);
	~qo_node();



	//
	// rfc handling
	//

	unsigned int get_ref() const;
	unsigned int inc_ref();
	unsigned int dec_ref();



	//
	// basic functions
	//

	const quadratic_order *get_qo() const;
	quadratic_order *get_qo();
	void set_qo(quadratic_order *QO2);
	void set_succ(qo_node *e);
	qo_node *succ();
};




//
// Class: qo_list
//
// This class represents a list of all the currently defined quadratic
//    orders.  The idea is to avoid multiple definitions of the same order
//    by counting the number of references to a specific order.  The list
//    is implemented as a simple linked list of qo_node's.
//

class qo_list
{
protected:

	qo_node *head; // head of the list
	quadratic_order *last_qo; // last quadratic_order that was referenced



public:

	//
	// constructor and destructor
	//

	qo_list();
	~qo_list();



	//
	// basic functions
	//

	bool is_empty() const;

	const qo_node *first() const;
	qo_node *first();

	const quadratic_order *last() const;
	quadratic_order *last();

	void set_last(quadratic_order *QO);

	void print(std::ostream & out = std::cout) const;
	void nullify(quadratic_order *QO2);
	qo_node *insert(quadratic_order & QO2, bool is_dyn);
	int clear(qo_node * elmnt);

	qo_node *add_to_list(qo_node *old, quadratic_order & newQO);
	qo_node *add_last_to_list(qo_node *old);
	quadratic_order *add_to_list(const bigint & Delta);

	// < MM >
	qo_node *add_to_list(qo_node *old_node, qo_node *new_node);
	friend void swap (qo_node* &x, qo_node* &y);
	// < /MM >
};

void swap (qo_node* &x, qo_node* &y);


#endif	// HEADBANGER



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif



#endif	// LIDIA_QO_LIST_H_GUARD_
