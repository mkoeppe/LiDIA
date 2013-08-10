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
//	Author	: Thomas Papanikolaou (TP), Patrick Theobald (PT), 
//                Christoph Ludwig (CL)
//	Changes	: See CVS log
//
//==============================================================================================


#ifndef LIDIA_ERROR_H_GUARD_
#define LIDIA_ERROR_H_GUARD_



#ifndef LIDIA_LIDIA_H_GUARD_
# include "LiDIA/LiDIA.h"
#endif

#include  <stdexcept>   // defines class std::runtime_error
#include  <typeinfo>    // defines class std::bad_cast

#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif



/***********************************************************
 **
 ** The new (exception driven) error handling implementation
 **
 ************************************************************/


template<class ExceptionType>
inline
void error_handler(const ExceptionType& ex) {
  // Maybe we want a compile time assert to make
  // sure ExceptionType is derived from basic_error.
#ifdef LIDIA_EXCEPTIONS
  throw ex;
#else
  ex.traditional_error_handler();
#endif
}

/**
 ** class basic_error
 ** ================
 **
 ** basic_error is the root of LiDIA's exception hierarchy.
 ** It has a pure virtual member function, so it cannot be instantiated.
 **
 ** Note: We have to override the destructor in order to provide
 **       an empty exception specification! (That's an requirement
 **       inherited from std::exception. It holds for all derived
 **       classes, too.)
 **
 ** traditional_error_handler() calls lidia_error_handler() and does
 ** not return.
 **
 **/ 
class basic_error : public std::exception {
public:
  virtual ~basic_error() throw() = 0;
  virtual const char* what() const throw();
  virtual const std::string& offendingClass() const = 0;
  void traditional_error_handler() const;

protected:
  explicit basic_error(const std::string& what_msg);
  virtual  void traditional_error_handler_impl(std::string& classname,
					       std::string& msg) const;

private:
  std::string msg_;
};
  

/**
 ** class generic_error
 ** ==================
 **
 ** generic_error is the class actually thrown by LiDIA if, for some 
 ** particular error situation, no specialized exception
 ** class is available.
 **
 ** The two parameters of the second c'tor correspond to the parameters of
 ** lidia_error_handler. The first c'tor defaults classname to 
 ** "<generic_error>". offendingClass() returns the value of the
 ** c'tor's parameter classname.
 **/   
class generic_error : public basic_error {
public:
  explicit generic_error(const std::string& what_msg);
  generic_error(const std::string& classname,
	       const std::string& what_msg);
  virtual ~generic_error() throw();
  
  virtual const std::string& offendingClass() const;

private:
  std::string offendingClassName_;
};


/**
 ** index_out_of_bounds_error
 ** =====================
 **
 ** index_out_of_bounds_error is thrown if a sequence is subscripted with
 ** an invalid index. The index that caused the error is available through
 ** offendingIndex().
 **/
class index_out_of_bounds_error : public basic_error, std::out_of_range {
public:
  index_out_of_bounds_error(const std::string& classname, long n);
  virtual ~index_out_of_bounds_error() throw();
  virtual const char* what() const throw();
  
  virtual const std::string& offendingClass() const;
  long offendingIndex() const;
  
private:
  std::string offendingClassName_;
  long index_;
};


/**
 ** cast_error
 ** =========
 **
 ** If LiDIA has to perform a cast but cannot guarantee that
 ** the cast will succeed (for example in precondition_error.paramValue<T>),
 ** then it will throw a cast_error in case of failure.
 */
class cast_error : public basic_error, std::bad_cast {
public:
  cast_error(const std::string& classname,
	    const std::string& what_msg);
  virtual ~cast_error() throw();
  virtual const char* what() const throw();
  
  virtual const std::string& offendingClass() const;
  
private:
  std::string offendingClassName_;
};


/***********************************************************
 **
 ** The legacy error handling implementation
 **
 ************************************************************/


typedef void (*error_handler_ptr)(const char *, const char *);

extern error_handler_ptr lidia_error_handler_internal;
void default_error_handler(const char *, const char *);
error_handler_ptr set_error_handler(error_handler_ptr);


#define lidia_error_handler_c(f, m, code) { code; lidia_error_handler(f, m); }



inline void lidia_error_handler(const char *file, const char *msg)
{
	error_handler(generic_error(file, msg));
}



inline void lidia_error_handler(const char *proto, const char *file, const char *msg)
{
	std::cerr << "FUNCTION: " << std::endl << proto << std::endl << std::endl;
	std::cerr << "MESSAGE: " << std::endl;
	lidia_error_handler(file, msg);
}


#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif


#endif	// LIDIA_ERROR_H_GUARD_

