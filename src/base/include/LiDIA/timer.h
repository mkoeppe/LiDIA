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
//	Author	: Thomas Papanikolaou (TP)
//	Changes	: See CVS log
//
//==============================================================================================


#ifndef LIDIA_TIMER_H_GUARD_
#define LIDIA_TIMER_H_GUARD_



#ifndef LIDIA_LIDIA_H_GUARD_
# include	"LiDIA/LiDIA.h"
#endif



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif



enum {
	TIME_MODE,
	HMS_MODE
};



class timer {
private:

  static const double GRANULARITY;

  static long _clk_tcks;
  
  long user_t, sys_t;
  long t_user_t, t_sys_t;
  int print_mode;
  
  void print_hms(std::ostream & out = std::cout, long rt = 0) const;

public:


	static long clk_tcks ();

	timer();
	timer(const timer &);
	~timer();

	int set_print_mode(int x = 1);
	int get_print_mode() const;

	void start_timer();
	void stop_timer();
	void cont_timer();

	long user_time() const;
	long sys_time() const;
	long real_time() const;

	void print(std::ostream & out = std::cout) const;

	timer & operator = (const timer & t);

        timer & operator + (const timer &);
        timer & operator - (const timer &);
        timer & operator / (double);

	friend std::ostream & operator << (std::ostream & out, const timer & t);

};



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif



#endif	// LIDIA_TIMER_H_GUARD_
