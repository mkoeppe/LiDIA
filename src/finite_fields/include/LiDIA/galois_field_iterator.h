//==============================================================================================
//
//	This file is part of LiDIA --- a library for computational number theory
//
//	Copyright (c) 1994--2004 the LiDIA Group.  All rights reserved.
//
//	See http://www.informatik.tu-darmstadt.de/TI/LiDIA/
//
//----------------------------------------------------------------------------------------------
//
//	$Id$
//
//	Author	: Christoph Ludwig
//	Changes	: See CVS log
//
//==============================================================================================

#ifndef LIDIA_GALOIS_FIELD_ITERATOR_H_GUARD_
#define LIDIA_GALOIS_FIELD_ITERATOR_H_GUARD_

// ISO C++ headers
#include <iterator>
#include <memory>

// LiDIA headers
#include "LiDIA/LiDIA.h"


#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif

  class galois_field;
  class gf_element;

  class galois_field_iterator {
  private:
    bool past_the_end;
    std::auto_ptr<galois_field> enumerated_field;
    std::auto_ptr<gf_element> current_element;

    static std::auto_ptr<galois_field const> const dummy_field;

    void check_access(std::string const& caller) const;

    friend void swap(galois_field_iterator& lhs,
                     galois_field_iterator& rhs);

    friend bool operator==(galois_field_iterator const& lhs,
                           galois_field_iterator const& rhs);

  public:
    typedef gf_element const value_type;
    typedef std::bidirectional_iterator_tag iterator_category;
    typedef value_type* pointer;
    typedef value_type& reference;

    galois_field_iterator();
    explicit galois_field_iterator(galois_field const& gf);
    explicit galois_field_iterator(gf_element const& element);
    galois_field_iterator(galois_field_iterator const& iter);

    ~galois_field_iterator();
    
    galois_field_iterator& operator=(galois_field_iterator const& iter);

    galois_field const& get_field() const;
	
    value_type& operator*() const; 
    value_type* operator->() const;
    galois_field_iterator& operator++();
    galois_field_iterator operator++(int);
    galois_field_iterator& operator--();
    galois_field_iterator operator--(int);

    void move_past_the_end();
  };


  void swap(galois_field_iterator& lhs,
            galois_field_iterator& rhs);

  bool operator==(galois_field_iterator const& lhs,
                  galois_field_iterator const& rhs);
  bool operator!=(galois_field_iterator const& lhs,
                  galois_field_iterator const& rhs);


#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif


#endif // LIDIA_GALOIS_FIELD_ITERATOR_H_GUARD_


/**
 ** Local Variables:
 ** mode: C++
 ** c-file-offsets: ((case-label            . 2)
 **                  (statement-block-intro . +)
 **                  (knr-argdecl-intro     . 0)
 **                  (substatement-open     . 0)
 **                  (label                 . 0)
 **                  (statement-cont        . +))
 ** c-basic-offset: 2
 ** c-comment-only-line-offset: 0
 ** c-hanging-braces-alist: ((brace-list-open)
 **                          (brace-entry-open)
 **                          (substatement-open after)
 **                          (block-close . c-snug-do-while))
 ** c-cleanup-list: (brace-else-brace)
 ** indent-tabs-mode: nil
 ** End:
 **/
