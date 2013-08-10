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


#ifdef HAVE_CONFIG_H
# include	"config.h"
#endif
#include	"LiDIA/timer.h"
#include	<unistd.h>
#include	<sys/times.h>



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif


const double timer::GRANULARITY = 1000.0;
long timer::_clk_tcks = 0;



long
timer::clk_tcks ()
{
	if (_clk_tcks == 0) {
		_clk_tcks = sysconf(_SC_CLK_TCK);
	}
	return _clk_tcks;
}



void
timer::print_hms (std::ostream & out, long rt) const
{
  unsigned long d, h, m, s;

  if (rt) {
    d = static_cast< long > (rt / (24 * 60 * 60 * GRANULARITY));
    if (d) {
      out << d << " day ";
      rt -= static_cast< long >(d * (24 * 60 * 60 * GRANULARITY));
    }
    
    h = static_cast< long >(rt / (60 * 60 * GRANULARITY));
    if (h) {
      out << h << " hour ";
      rt -= static_cast< long >(h * (60 * 60 * GRANULARITY));
    }
    
    m = static_cast< long >(rt / (60 * GRANULARITY));
    if (m) {
      out << m << " min ";
      rt -= static_cast< long >(m * (60 * GRANULARITY));
    }
    
    s = static_cast< long >(rt / GRANULARITY);
    if (s) {
      out << s << " sec ";
      rt -= static_cast< long >(s * GRANULARITY);
    }
    
    if (rt)
      out << rt << " msec";
  }
  else 
    {
      out << "0 msec";
    }
}



timer::timer ()
{
  user_t = 0;
  sys_t = 0;
  t_user_t = 0;
  t_sys_t = 0;
  print_mode = TIME_MODE;
  if (_clk_tcks == 0) {
    _clk_tcks = sysconf(_SC_CLK_TCK);
  }
}



timer::timer (const timer & t)
{
  user_t = t.user_t;
  sys_t = t.sys_t;
  t_user_t = t.t_user_t;
  t_sys_t = t.t_sys_t;
  print_mode = t.print_mode;
}



timer::~timer ()
{
}



int
timer::set_print_mode (int m)
{
  int old_print_mode = print_mode;
  
  
  switch (m) {
  case TIME_MODE :
  case HMS_MODE :
    print_mode = m;
    break;
  default:
    warning_handler("timer", "unknown print mode. Set to HMS_MODE");
    print_mode = 1;
    break;
  }
  return old_print_mode;
}



int
timer::get_print_mode () const
{
  return print_mode;
}

void
timer::start_timer ()
{
#ifdef HAVE_POSIX_TIMES
  struct tms	buffer;
  
  
  times(&buffer);
  user_t = static_cast< long >(buffer.tms_utime *
			       (GRANULARITY / _clk_tcks));
  sys_t = static_cast< long >(buffer.tms_stime * (GRANULARITY / _clk_tcks));
#endif
  t_user_t = 0;
  t_sys_t = 0;
}


void
timer::stop_timer()
{
#ifdef HAVE_POSIX_TIMES
	struct tms	buffer;

	times(&buffer);
	t_user_t += static_cast< long >(buffer.tms_utime * 
					(GRANULARITY / _clk_tcks)) - user_t;
	t_sys_t += static_cast< long >(buffer.tms_stime * 
				       (GRANULARITY / _clk_tcks)) - sys_t;
#endif
}



void
timer::cont_timer ()
{
#ifdef HAVE_POSIX_TIMES
	struct tms	buffer;

	times(&buffer);
	user_t = static_cast< long >(buffer.tms_utime * 
				     (GRANULARITY / _clk_tcks));
	sys_t = static_cast< long >(buffer.tms_stime * 
				    (GRANULARITY / _clk_tcks));
#endif
}



long
timer::user_time () const
{
	return t_user_t;
}



long
timer::sys_time () const
{
	return t_sys_t;
}



long
timer::real_time() const
{
  return (t_user_t + t_sys_t);
}



void
timer::print(std::ostream & out) const
{
  switch (print_mode) {
  case TIME_MODE:
    out << real_time() << " real\t";
    out << user_time() << " user\t";
    out << sys_time() << " sys";
    break;
  case HMS_MODE:
    print_hms(out, real_time());
    break;
  default:
    warning_handler("timer", "unknown print mode. Use HMS_MODE");
    print_hms(out, real_time());
    break;
  }
}



timer &
timer::operator = (const timer & t)
{
	user_t = t.user_t;
	sys_t = t.sys_t;
	t_user_t = t.t_user_t;
	t_sys_t = t.t_sys_t;
	print_mode = t.print_mode;
	return *this;
}

timer & 
timer::operator + (const timer & t)
{
  user_t += t.user_t;
  sys_t += t.sys_t;
  t_user_t += t.t_user_t;
  t_sys_t += t.t_sys_t;
  return *this;
}

timer & 
timer::operator - (const timer & t)
{
  user_t -= t.user_t;
  if (user_t < 0)
    user_t = 0;
  sys_t -= t.sys_t;
  if (sys_t < 0)
    sys_t = 0;
  t_user_t -= t.t_user_t;
  if (t_user_t < 0)
    t_user_t = 0;
  t_sys_t -= t.t_sys_t;
  if (t_sys_t < 0)
    t_sys_t = 0;
  return *this;
}

timer & 
timer::operator / (double x)
{
  if (x <= 0.0)
    lidia_error_handler("timer","operator/::scalar non-positive");

  user_t = static_cast<long>(user_t / x);
  sys_t = static_cast<long>(sys_t / x);
  t_user_t = static_cast<long>(t_user_t / x);
  t_sys_t = static_cast<long>(t_sys_t / x);
  return *this;
}


std::ostream &
operator << (std::ostream & out, const timer & t)
{
	t.print(out);
	return out;
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
