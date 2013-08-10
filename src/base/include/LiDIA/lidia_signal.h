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
//	Author	: Markus Maurer (MM)
//	Changes	: See CVS log
//
//==============================================================================================


#ifndef LIDIA_LIDIA_SIGNAL_H_GUARD_
#define LIDIA_LIDIA_SIGNAL_H_GUARD_



#ifndef LIDIA_LIDIA_H_GUARD_
# include	"LiDIA/LiDIA.h"
#endif
#include	<signal.h>



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif



// defines and typedefs
// for signal handling

#ifndef SIGHUP
#define LIDIA_SIGHUP  1       // hangup
#else
#define LIDIA_SIGHUP SIGHUP
#endif

#ifndef SIGINT
#define LIDIA_SIGINT  2       // interrupt (rubout)
#else
#define LIDIA_SIGINT SIGINT
#endif

#ifndef SIGQUIT
#define LIDIA_SIGQUIT 3       // quit (ASCII FS)
#else
#define LIDIA_SIGQUIT SIGQUIT
#endif

#ifndef SIGILL
#define LIDIA_SIGILL  4       // illegal instruction (not reset when caught)
#else
#define LIDIA_SIGILL SIGILL
#endif

#ifndef SIGTRAP
#define LIDIA_SIGTRAP 5       // trace trap (not reset when caught)
#else
#define LIDIA_SIGTRAP SIGTRAP
#endif

#ifndef SIGIOT
#define LIDIA_SIGIOT  6       // IOT instruction
#else
#define LIDIA_SIGIOT SIGIOT
#endif

#ifndef SIGABRT
#define LIDIA_SIGABRT 6       // used by abort, replace SIGIOT in the future
#else
#define LIDIA_SIGABRT SIGABRT
#endif

#ifndef SIGEMT
#define LIDIA_SIGEMT  7       // EMT instruction
#else
#define LIDIA_SIGEMT SIGEMT
#endif

#ifndef SIGFPE
#define LIDIA_SIGFPE  8       // floating point exception
#else
#define LIDIA_SIGFPE SIGFPE
#endif

#ifndef SIGKILL
#define LIDIA_SIGKILL 9       // kill (cannot be caught or ignored)
#else
#define LIDIA_SIGKILL SIGKILL
#endif

#ifndef SIGBUS
#define LIDIA_SIGBUS  10      // bus error
#else
#define LIDIA_SIGBUS  SIGBUS
#endif

#ifndef SIGSEGV
#define LIDIA_SIGSEGV 11      // segmentation violation
#else
#define LIDIA_SIGSEGV SIGSEGV
#endif

#ifndef SIGSYS
#define LIDIA_SIGSYS  12      // bad argument to system call
#else
#define LIDIA_SIGSYS SIGSYS
#endif

#ifndef SIGPIPE
#define LIDIA_SIGPIPE 13      // write on a pipe with no one to read it
#else
#define LIDIA_SIGPIPE SIGPIPE
#endif

#ifndef SIGALRM
#define LIDIA_SIGALRM 14      // alarm clock
#else
#define LIDIA_SIGALRM SIGALRM
#endif

#ifndef SIGTERM
#define LIDIA_SIGTERM 15      // software termination signal from kill
#else
#define LIDIA_SIGTERM SIGTERM
#endif

#ifndef SIGUSR1
#define LIDIA_SIGUSR1 16      // user defined signal 1
#else
#define LIDIA_SIGUSR1 SIGUSR1
#endif

#ifndef SIGUSR2
#define LIDIA_SIGUSR2 17      // user defined signal 2
#else
#define LIDIA_SIGUSR2 SIGUSR2
#endif

#ifndef SIGCLD
#define LIDIA_SIGCLD  18      // child status change
#else
#define LIDIA_SIGCLD  SIGCLD
#endif

#ifndef SIGCHLD
#define LIDIA_SIGCHLD 18      // child status change alias (POSIX)
#else
#define LIDIA_SIGCHLD SIGCHLD
#endif

#ifndef SIGPWR
#define LIDIA_SIGPWR  19      // power-fail restart
#else
#define LIDIA_SIGPWR  SIGPWR
#endif

#ifndef SIGWINCH
#define LIDIA_SIGWINCH 20     // window size change
#else
#define LIDIA_SIGWINCH SIGWINCH
#endif

#ifndef SIGURG
#define LIDIA_SIGURG  21      // urgent socket condition
#else
#define LIDIA_SIGURG SIGURG
#endif

#ifndef SIGPOLL
#define LIDIA_SIGPOLL 22      // pollable event occured
#else
#define LIDIA_SIGPOLL SIGPOLL
#endif

#ifndef SIGIO
#define LIDIA_SIGIO   LIDIA_SIGPOLL // socket I/O possible (SIGPOLL alias)
#endif

#ifndef SIGSTOP
#define LIDIA_SIGSTOP 23      // stop (cannot be caught or ignored)
#else
#define LIDIA_SIGSTOP SIGSTOP
#endif

#ifndef SIGTSTP
#define LIDIA_SIGTSTP 24      // user stop requested from tty
#else
#define LIDIA_SIGTSTP SIGTSTP
#endif

#ifndef SIGCONT
#define LIDIA_SIGCONT 25      // stopped process has been continued
#else
#define LIDIA_SIGCONT SIGCONT
#endif

#ifndef SIGTTIN
#define LIDIA_SIGTTIN 26      // background tty read attempted
#else
#define LIDIA_SIGTTIN SIGTTIN
#endif

#ifndef SIGTTOU
#define LIDIA_SIGTTOU 27      // background tty write attempted
#else
#define LIDIA_SIGTTOU SIGTTOU
#endif

#ifndef SIGVTALRM
#define LIDIA_SIGVTALRM 28    // virtual timer expired
#else
#define LIDIA_SIGVTALRM SIGVTALRM
#endif

#ifndef SIGPROF
#define LIDIA_SIGPROF 29      // profiling timer expired
#else
#define LIDIA_SIGPROF SIGPROF
#endif

#ifndef SIGXCPU
#define LIDIA_SIGXCPU 30      // exceeded cpu limit
#else
#define LIDIA_SIGXCPU SIGXCPU
#endif

#ifndef SIGXFSZ
#define LIDIA_SIGXFSZ 31      // exceeded file size limit
#else
#define LIDIA_SIGXFSZ SIGXFSZ
#endif

#ifndef SIGWAITING
#define LIDIA_SIGWAITING 32   // process's lwps are blocked
#else
#define LIDIA_SIGWAITING SIGWAITING
#endif


// we allow only POSIX
#define LIDIA_SIGNAL_FUNCTION(f) extern "C" void f(int)
extern "C" { // signal handlers are C functions
    typedef void (*lidia_signal_handler_t)(int);
}


class lidia_signal {
private:

	int signal_num; // the signal, i.e. SIGINT
	// avoid name collision with function signal

public:

	typedef struct {
		int sig;
		int offset, width;
		int sp; // stack pointer
		int ss; // stack size
		lidia_signal_handler_t* stack;
	} signal_stack;

private:

	static signal_stack signals[];
	static const int    default_stack_size;


	static void allocate (lidia_signal_handler_t * &s,
			      int old_size,
			      int new_size);

	static void print_stack_index (int index);

	static void install_handler (int sig,
				     lidia_signal_handler_t handler);
	static void uninstall_handler (int sig);


public:

	lidia_signal(int sig, lidia_signal_handler_t handler);
	~lidia_signal();

	static void print_stack (int sig);

};



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif



#endif	// LIDIA_LIDIA_SIGNAL_H_GUARD_
