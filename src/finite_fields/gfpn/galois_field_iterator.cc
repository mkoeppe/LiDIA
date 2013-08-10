// ISO C++ headers
#include <cassert>
#include <algorithm>

// LiDIA headers
#include "LiDIA/galois_field_iterator.h"
#include "LiDIA/galois_field.h"
#include "LiDIA/gf_element.h"
#include "LiDIA/precondition_error.h"


#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif

  std::auto_ptr<galois_field const> const
  galois_field_iterator::dummy_field(new galois_field());


  galois_field_iterator::
  galois_field_iterator()
    : past_the_end(true),
      enumerated_field(0),
      current_element(0) {
    
  }


  galois_field_iterator::
  galois_field_iterator(galois_field const& gf)
    : past_the_end(gf == *galois_field_iterator::dummy_field),
      enumerated_field(0),
      current_element(0) {
    if(!this->past_the_end) {
      this->enumerated_field.reset(new galois_field(gf));
      this->current_element.reset(new gf_element(gf));
    }
  }
    

  galois_field_iterator::
  galois_field_iterator(gf_element const& element)
    : past_the_end(false),
      enumerated_field(0),
      current_element(0) {
    this->enumerated_field.reset(new galois_field(element.get_field()));
    this->current_element.reset(new gf_element(element));
  }


  galois_field_iterator::
  galois_field_iterator(galois_field_iterator const& iter)
    : past_the_end(iter.past_the_end),
      enumerated_field(0),
      current_element(0) {
    if(iter.enumerated_field.get()) {
      this->enumerated_field.reset(new galois_field(*iter.enumerated_field));
      this->current_element.reset(new gf_element(*iter.current_element));
    }
  }

  
  galois_field_iterator::
  ~galois_field_iterator() {
    // this is not strictly necessary because the compiler will insert a call
    // to the destructors of enumerated_field and current_element anyway
    this->enumerated_field.release();
    this->current_element.release();
  }


  galois_field_iterator&
  galois_field_iterator::
  operator=(galois_field_iterator const& iter) {
    if(this != &iter) {
      galois_field_iterator tmp(iter);
      swap(*this, tmp);
    }
    return *this;
  }


  void
  galois_field_iterator::
  check_access(std::string const& caller) const {
    if(!this->enumerated_field.get()) {
      precondition_error ex(caller,
                            "galois_field_iterator",
                            "this->field() must not be a default constructed "
                            "dummy galois_field");
    }
    if(this->past_the_end) {
      precondition_error ex(caller,
                            "galois_field_iterator",
                            "*this != this->field().end() required");
      error_handler(ex);
    }
  }


  galois_field const&
  galois_field_iterator::
  get_field() const {
    if(this->enumerated_field.get()) {
      return *this->enumerated_field;
    }
    else {
      return *galois_field_iterator::dummy_field;
    }
  }


  galois_field_iterator::value_type&
  galois_field_iterator::
  operator*() const {
    this->check_access("operator*()");
    return *this->current_element;
  }
  

  galois_field_iterator::value_type*
  galois_field_iterator::
  operator->() const {
    this->check_access("operator->()");
    return this->current_element.get();
  }
  

  galois_field_iterator&
  galois_field_iterator::
  operator++() {
    this->check_access("operator++()");

    Fp_polynomial elem_rep = this->current_element->polynomial_rep();
    bigint const p_minus_one = elem_rep.modulus() - 1;
    lidia_size_t field_degree = this->enumerated_field->degree();

    bool finished = false;
    bigint coeff;
    for(lidia_size_t i = 0; i < field_degree && !finished; ++i) {
      elem_rep.get_coefficient(coeff, i);
      if(coeff == p_minus_one) {
        coeff.assign_zero();
      }
      else {
        ++coeff;
        finished = true;
      }
      elem_rep.set_coefficient(coeff, i);
    }

    if(!finished) {
      elem_rep.assign_zero();
      this->past_the_end = true;
    }
    this->current_element->set_polynomial_rep(elem_rep);
    return *this;
  }


  galois_field_iterator
  galois_field_iterator::
  operator++(int) {
    galois_field_iterator stored_value(*this);
    ++*this;
    return stored_value;
  }


  galois_field_iterator&
  galois_field_iterator::
  operator--() {
    if(!this->enumerated_field.get()) {
      precondition_error ex("operator--()",
                            "galois_field_iterator",
                            "this->field() must not be a default constructed "
                            "dummy galois_field");
    }
    if(this->current_element->is_zero() && !this->past_the_end) {
      precondition_error ex("operator--()",
                            "galois_field_iterator",
                            "*this != this->field().begin() required");
      error_handler(ex);
    }

    Fp_polynomial elem_rep = this->current_element->polynomial_rep();
    lidia_size_t field_degree = this->enumerated_field->degree();

    if(this->past_the_end) {

      bigint p_minus_one = elem_rep.modulus() - 1;
      elem_rep.set_max_degree(field_degree - 1);
      for(lidia_size_t i = 0; i < field_degree; ++i) {
        elem_rep.set_coefficient(p_minus_one, i);
      }
      this->past_the_end = false;

    }
    else {

      bool finished = false;
      bigint coeff;
      for(lidia_size_t i = 0; !finished && i < field_degree; ++i) {
        elem_rep.get_coefficient(coeff, i);
        if(coeff.is_zero()) {
          elem_rep.set_coefficient(elem_rep.modulus() - 1, i);
        }
        else {
          --coeff;
          elem_rep.set_coefficient(coeff, i);
          finished = true;
        }
      }
      assert(finished);

    }        

    this->current_element->set_polynomial_rep(elem_rep);
    return *this;
  }

    
  galois_field_iterator
  galois_field_iterator::
  operator--(int) {
    galois_field_iterator stored_value(*this);
    --*this;
    return stored_value;
  }


  void
  galois_field_iterator::
  move_past_the_end() {
    this->past_the_end = true;
  }


  namespace {

    template<typename T>
    void swap_auto_ptr(std::auto_ptr<T>& lhs, std::auto_ptr<T>& rhs) {
      std::auto_ptr<T> tmp = lhs;
      lhs = rhs;
      rhs = tmp;
    }

  }

  void swap(galois_field_iterator& lhs,
            galois_field_iterator& rhs) {
    std::swap(lhs.past_the_end, rhs.past_the_end);
    swap_auto_ptr(lhs.enumerated_field, rhs.enumerated_field);
    swap_auto_ptr(lhs.current_element, rhs.current_element);
  }


  bool
  operator==(galois_field_iterator const& lhs,
             galois_field_iterator const& rhs) {
    if(!lhs.enumerated_field.get() && !rhs.enumerated_field.get()) {
      // both iterators belong to the dummy field
      return true;
    }
    else if(*lhs.enumerated_field != *rhs.enumerated_field) {
      // the iterators cannot be equal since they belong to different fields
      return false;
    }
    else {
      // the iterators belong to the same field.
      // they are equal if they both point past the end or if they point
      // to the same element.
      return (lhs.past_the_end == rhs.past_the_end &&
              (lhs.past_the_end ||
               *lhs.current_element == *rhs.current_element));
    }
  }


  bool
  operator!=(galois_field_iterator const& lhs,
             galois_field_iterator const& rhs) {
    return !(lhs == rhs);
  }


#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif

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


