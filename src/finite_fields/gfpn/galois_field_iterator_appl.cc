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


// ISO C++ headers
#include <iostream>
#include <cstdlib>

// LiDIA headers
#include "LiDIA/galois_field.h"
#include "LiDIA/gf_element.h"
#include "LiDIA/galois_field_iterator.h"


#ifdef LIDIA_NAMESPACE
using namespace LiDIA;
#endif


bool test_case(galois_field const& gf) {
  bool result = true;
  galois_field_iterator begin_iter = gf.begin();
  galois_field_iterator end_iter = gf.end();

  long card_field;
  gf.number_of_elements().longify(card_field);
  if(card_field == 1) {
    // Unfortunately, if gf is the "dummy" field galois_field(), then
    // gf.number_of_elements() returns 0^0 = 1. Fix this here.
    card_field = 0;
  }

  long counter = 0;
  std::cout << "  All elements:";
  for(galois_field_iterator iter = begin_iter; iter != end_iter; ++iter) {
    std::cout << " " << *iter;
    ++counter;
  }
  std::cout << std::endl;

  if(counter == card_field) {
    std::cout << "  OK, we enumerated " << counter << " elements."
              << std::endl;
  }
  else {
    std::cerr << "  Error: We enumerated " << counter
              << " elements, but the field has " << gf.number_of_elements()
              << " elements!" << std::endl;
    result = false;
  }

  std::cout << "  All elements backwards:";
  galois_field_iterator iter = end_iter;
  while(iter != begin_iter) {
    --iter;
    std::cout << " " << *iter;
    --counter;
  }
  std::cout << std::endl;

  if(counter == 0) {
    std::cout << "  OK, we enumerated backwards the same number of elements."
              << std::endl;
  }
  else {
    std::string adv;
    if(counter > 0) {
      adv = "less";
    }
    else {
      adv = "more";
      counter *= -1;
    }
    std::cout << "  Error: We enumerated backwards " << counter
              << " elements " << adv << " than forward." << std::endl;
    result = false;
  }

  if(begin_iter != end_iter) {
    gf_element elem(gf);
    elem.randomize();

    std::cout << "  Starting from " << elem << "...   ";
    counter = 0;
    for(iter = galois_field_iterator(elem); iter != begin_iter; --iter) {
      ++counter;
    }
    for(iter = galois_field_iterator(elem); iter != end_iter; ++iter) {
      ++counter;
    }
    if(counter == card_field) {
      std::cout << "success" << std::endl;
    }
    else {
      std::cout << "failure" << std::endl;
      result = false;
    }

    iter = galois_field_iterator(elem);
    std::cout << "  Testing operator*()...   ";
    if(*iter == elem) {
      std::cout << "success" << std::endl;
    }
    else {
      std::cout << "failure" << std::endl;
      result = false;
    }
    
    iter = galois_field_iterator(elem);
    std::cout << "  Testing operator->()...   ";
    if(*(iter.operator->()) == elem) {
      std::cout << "success" << std::endl;
    }
    else {
      std::cout << "failure" << std::endl;
      result = false;
    }

    galois_field_iterator iter2 = iter;
    galois_field_iterator iter3 = iter;

    std::cout << "  Testing postincrement...   ";
    if(iter == iter2++) {
      std::cout << "success" << std::endl;
    }
    else {
      std::cout << "failure" << std::endl;
      result = false;
    }
    
    std::cout << "  Testing preincrement...   ";
    if(iter2 == ++iter3) {
      std::cout << "success" << std::endl;
    }
    else {
      std::cout << "failure" << std::endl;
      result = false;
    }
    
    std::cout << "  Testing postdecrement...   ";
    if(iter2 == iter3--) {
      std::cout << "success" << std::endl;
    }
    else {
      std::cout << "failure" << std::endl;
      result = false;
    }
    
    std::cout << "  Testing predecrement...   ";
    if(iter == --iter2) {
      std::cout << "success" << std::endl;
    }
    else {
      std::cout << "failure" << std::endl;
      result = false;
    }
  }

  return result;
}
  
  
int main_LiDIA(int argc, char** argv) {
  bigint prime = 7;

  galois_field dummy;
  galois_field prime_field(prime);
  galois_field deg2_field(prime, 2);
  galois_field deg3_field(prime, 3);

  std::cout << "Dummy field:" << std::endl;
  bool success = test_case(dummy);

  std::cout << "\n\nGF(" << prime << "):" << std::endl;
  success &= test_case(prime_field);

  std::cout << "\n\nGF(" << prime << "^2):" << std::endl;
  success &= test_case(deg2_field);

  std::cout << "\n\nGF(" << prime << "^3):" << std::endl;
  success &= test_case(deg3_field);

  std::cout << std::endl;
  if(success) {
    std::cout << "OK, all tests passed.\n" << std::endl;
    return EXIT_SUCCESS;
  }
  else {
    std::cout << "There were failures!\n" << std::endl;
    return EXIT_FAILURE;
  }
}



int main(int argc, char** argv) {

#if defined(LIDIA_EXCEPTIONS)
  try {
#endif

    return main_LiDIA(argc, argv);
	
#if defined(LIDIA_EXCEPTIONS)
  }
  catch(basic_error const& ex) {
    ex.traditional_error_handler();
    return 1;
  }
  catch(std::exception const& ex) {
    std::cerr << "unexpected exception: " << ex.what() << "\n";
    return 1;
  }
#endif
  
}

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
