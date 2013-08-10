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


#ifdef HAVE_CONFIG_H
# include	"config.h"
#endif
#include	"LiDIA/lidia_signal.h"



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



const int lidia_signal::default_stack_size = 16;

lidia_signal::signal_stack lidia_signal::signals[] = {
#ifdef LIDIA_SIGHUP
	{ LIDIA_SIGHUP, 0, 2, -1, 0, NULL },
#endif
#ifdef LIDIA_SIGINT
	{ LIDIA_SIGINT, 2, 4, -1, 0, NULL },
#endif
#ifdef LIDIA_SIGQUIT
	{ LIDIA_SIGQUIT, 6, 2, -1, 0, NULL },
#endif
#ifdef LIDIA_SIGILL
	{ LIDIA_SIGILL, 8, 2, -1, 0, NULL },
#endif
#ifdef LIDIA_SIGTRAP
	{ LIDIA_SIGTRAP, 10, 2, -1, 0, NULL },
#endif
#ifdef LIDIA_SIGIOT
	{ LIDIA_SIGIOT, 12, 2, -1, 0, NULL },
#endif
#ifdef LIDIA_SIGABRT
	{ LIDIA_SIGABRT, 14, 2, -1, 0, NULL },
#endif
#ifdef LIDIA_SIGEMT
	{ LIDIA_SIGEMT, 16, 2, -1, 0, NULL },
#endif
#ifdef LIDIA_SIGTFPE
	{ LIDIA_SIGFPE, 18, 2, -1, 0, NULL },
#endif
#ifdef LIDIA_SIGKILL
	{ LIDIA_SIGKILL, 20, 2, -1, 0, NULL },
#endif
#ifdef LIDIA_SIGBUS
	{ LIDIA_SIGBUS, 22, 2, -1, 0, NULL },
#endif
#ifdef LIDIA_SIGSEGV
	{ LIDIA_SIGSEGV, 24, 4, -1, 0, NULL },
#endif
#ifdef LIDIA_SIGSYS
	{ LIDIA_SIGSYS, 28, 2, -1, 0, NULL },
#endif
#ifdef LIDIA_SIGPIPE
	{ LIDIA_SIGPIPE, 30, 2, -1, 0, NULL },
#endif
#ifdef LIDIA_SIGALRM
	{ LIDIA_SIGALRM, 32, 2, -1, 0, NULL },
#endif
#ifdef LIDIA_SIGTERM
	{ LIDIA_SIGTERM, 34, 2, -1, 0, NULL },
#endif
#ifdef LIDIA_SIGUSR1
	{ LIDIA_SIGUSR1, 36, 2, -1, 0, NULL },
#endif
#ifdef LIDIA_SIGUSR2
	{ LIDIA_SIGUSR2, 38, 2, -1, 0, NULL },
#endif
#ifdef LIDIA_SIGCLD
	{ LIDIA_SIGCLD, 40, 2, -1, 0, NULL },
#endif
#ifdef LIDIA_SIGCHLD
	{ LIDIA_SIGCHLD, 42, 2, -1, 0, NULL },
#endif
#ifdef LIDIA_SIGPWR
	{ LIDIA_SIGPWR, 44, 2, -1, 0, NULL },
#endif
#ifdef LIDIA_SIGWINCH
	{ LIDIA_SIGWINCH, 46, 2, -1, 0, NULL },
#endif
#ifdef LIDIA_SIGURG
	{ LIDIA_SIGURG, 48, 2, -1, 0, NULL },
#endif
#ifdef LIDIA_SIGPOLL
	{ LIDIA_SIGPOLL, 50, 2, -1, 0, NULL },
#endif
#ifdef LIDIA_SIGIO
	{ LIDIA_SIGIO, 52, 2, -1, 0, NULL },
#endif
#ifdef LIDIA_SIGSTOP
	{ LIDIA_SIGSTOP, 54, 2, -1, 0, NULL },
#endif
#ifdef LIDIA_SIGTSTP
	{ LIDIA_SIGTSTP, 56, 2, -1, 0, NULL },
#endif
#ifdef LIDIA_SIGCONT
	{ LIDIA_SIGCONT, 58, 2, -1, 0, NULL },
#endif
#ifdef LIDIA_SIGTTIN
	{ LIDIA_SIGTTIN, 60, 2, -1, 0, NULL },
#endif
#ifdef LIDIA_SIGTTOU
	{ LIDIA_SIGTTOU, 62, 2, -1, 0, NULL },
#endif
#ifdef LIDIA_SIGVTALRM
	{ LIDIA_SIGVTALRM, 64, 2, -1, 0, NULL },
#endif
#ifdef LIDIA_SIGPROF
	{ LIDIA_SIGPROF, 66, 2, -1, 0, NULL },
#endif
#ifdef LIDIA_SIGXCPU
	{ LIDIA_SIGXCPU, 68, 2, -1, 0, NULL },
#endif
#ifdef LIDIA_SIGXFSZ
	{ LIDIA_SIGXFSZ, 70, 2, -1, 0, NULL },
#endif
#ifdef LIDIA_SIGWAITING
	{ LIDIA_SIGWAITING, 72, 2, -1, 0, NULL },
#endif
	{ 0, 0, 0, 0, 0, NULL }
};



//
// Function: static lidia_signal::allocate
//
// Parameters and Method:
//
// s == NULL
//  allocates new_size elements for s
//
// old_size <= new_size:
//   enlarges s to new_size elements
//   with elements 0,...,old_size-1
//   unchanged
//
// old_size >  new_size:
//   reduces s to new_size elements
//   with deleting the elements
//   new_size,...,old_size-1
//

void
lidia_signal::allocate (lidia_signal_handler_t * &s,
			int old_size,
			int new_size)
{
#ifdef HAVE_POSIX_SIGNALS
	if (old_size > new_size) {
		old_size = new_size;
	}

	if (s == NULL) {
		s = new lidia_signal_handler_t[new_size];
	}
	else {
		lidia_signal_handler_t *t;
		int i;

		t = new lidia_signal_handler_t[new_size];

		for (i = 0; i < old_size; i++)
			t[i] = s[i];

		delete [] s;
		s = t;
	}
#endif	// HAVE_POSIX_SIGNALS
}


//
// Function: lidia_signal::print_stack_index
//
// Parameter: The index of a stack in the
//            stack array lidia_signal::signals
//
//            Must be a valid index. No check.
//
// Method: Prints the content of signals[index].
//

void
lidia_signal::print_stack_index (int index)
{
#ifdef HAVE_POSIX_SIGNALS
	std::cout << "signal        : " << signals[index].sig << "\n";
	std::cout << "stack pointer : " << signals[index].sp << "\n";
	std::cout << "stack size    : " << signals[index].ss << "\n";

#if 0
	int i;
	std::cout << "stack         : ";

	if (signals[index].sp < 0)
		std::cout << "empty\n";
	else {
		for (i = 0; i <= signals[index].sp; i++)
			std::cout << " (" << signals[index].stack[i] << ")";
		std::cout << "\n";
	}
#endif
	std::cout.flush();
#endif	// HAVE_POSIX_SIGNALS
}



//
// Function: lidia_signal::print_stack
//
// Parameter: A signal
//
// Method: Searches for the index of the stack of the signal.
//         If found, calls print_stack_index(index) to display
//         the content of the stack. Otherwise calls the
//	   lidia_error_handler.
//

void
lidia_signal::print_stack (int sig)
{
#ifdef HAVE_POSIX_SIGNALS
	int i;
	int sig_found;

	sig_found = 0;

	for (i = 0; signals[i].sig && !sig_found; ++i) {
		if (signals[i].sig == sig) {
			// mark found
			sig_found = 1;

			// display stack
			lidia_signal::print_stack_index(i);
		}
	}

	if (!sig_found) {
		lidia_error_handler("lidia_signal::print_stack",
				    "Signal not found in signal list.");
	}
#endif	// HAVE_POSIX_SIGNALS
}



//
// Function: static lidia_signal::install_handler
//
// Parameter: sig, the signal for which a handler
//                 should be installed
//            handler, the new handler for sig
//
// Method:
//
// Searches for signal sig in the list of signals
// and puts the handler on top of the stack.
// If the signal was found:
//
//  * First install for sig:
//     store current handler for sig in stack[0]
//     and set stack[1] = handler
//
//  * At least second call for sig:
//    verify, that the current handler for sig
//    is on top of the stack;
//    increase stack pointer sp and
//    store stack[sp] = handler;
//
// If the signal was not found:
//    call the lidia_error_handler
//

void
lidia_signal::install_handler (int sig, lidia_signal_handler_t handler)
{
#ifdef HAVE_POSIX_SIGNALS
	int i, sig_found;
	lidia_signal_handler_t h;

	struct sigaction       iact, oact;


	sig_found = 0;

	// Search for signal sig in the list of signals.

	for (i = 0; signals[i].sig && !sig_found; ++i) {
		if (signals[i].sig == sig) {
			// mark found
			sig_found = 1;

			// store current handler in h
			// and set new handler for sig
			// h = sigset(sig, handler);

			iact.sa_handler = handler;
			sigemptyset(&iact.sa_mask);
			iact.sa_flags = 0;
#ifdef SA_RESTART
			iact.sa_flags |= SA_RESTART;
#endif
			sigaction(sig, &iact, &oact);
			h = oact.sa_handler;

			// allocate space for stack
			// when called for the first time
			if (signals[i].ss == 0) {
				signals[i].ss = lidia_signal::default_stack_size;
				lidia_signal::allocate(signals[i].stack, 0, signals[i].ss);
			}

			// First install for sig:
			//  store current handler for sig in stack[0]
			//  and set stack[1] = handler
			if (signals[i].sp == -1) {
				signals[i].sp = 1;
				signals[i].stack[0] = h;
				signals[i].stack[1] = handler;
			}

			// At least second call for sig:
			//  verify, that the current handler for sig
			//  is on top of the stack;
			else if (h != signals[i].stack[signals[i].sp]) {
				lidia_error_handler("lidia_signal::install_handler",
						    "Inconsistent signal stack.");
			}
			else {
				// increase stack pointer
				++signals[i].sp;

				// double stack space, if end is reached
				if (signals[i].sp == signals[i].ss) {
					lidia_signal::allocate(signals[i].stack,
							       signals[i].ss,
							       2*signals[i].ss);
					signals[i].ss *= 2;
				}

				// store handler on top of stack
				signals[i].stack[signals[i].sp] = handler;
			}
			// end if (signals[i].sp == -1)
		}
		// end if (signals[i].sig == sig)
	}
	// end for (i = 0; signals[i].sig; ++i)


	if (!sig_found) {
		lidia_error_handler("lidia_signal::install_handler",
				    "signal not found in signal list.");
	}
#endif	// HAVE_POSIX_SIGNALS
}
// end install_handler




//
// Function: static lidia_signal::uninstall_handler
//
// Parameter: sig, the signal for which the current
//                 handler should be uninstalled
//
// Method:
//
// Removes the current handler for sig from the top
// of the stack, decreases the stack pointer and
// reinstalls the previous handler for sig which is
// now on top of the stack.
//
// If the stack does not contain a lidia signal handler
// anymore (stack pointer = 0), the user defined signal
// handler in stack[0] is remove from top of stack.
// If the user now installs a new handler for sig and
// calls lidia functions, which again install
// their own handlers for sig, the new user handler
// is stored in stack[0] by the lidia_signal::install_handler
// function (see above).
//

void
lidia_signal::uninstall_handler (int sig)
{
#ifdef HAVE_POSIX_SIGNALS
	int i, sig_found;
	lidia_signal_handler_t h;

	struct sigaction iact, oact;


	sig_found = 0;

	for (i = 0; signals[i].sig && !sig_found; ++i) {
		// search for signal
		if (signals[i].sig == sig) {
			// mark found
			sig_found = 1;

			// no handler installed for sig
			if (signals[i].sp == -1) {
				lidia_warning_handler("lidia_signal::uninstall_handler",
						      "uninstall without previous install.");
			}
			// found installed handler
			else {
				// store current handler in h
				// and install second handler from stack
				// h = sigset(sig, signals[i].stack[signals[i].sp-1]);

				iact.sa_handler = signals[i].stack[signals[i].sp-1];
				sigemptyset(&iact.sa_mask);
				iact.sa_flags = 0;
#ifdef SA_RESTART
				iact.sa_flags |= SA_RESTART;
#endif
				sigaction(sig, &iact, &oact);
				h = oact.sa_handler;

				// current handler is not on top of stack
				if (h != signals[i].stack[signals[i].sp]) {
					lidia_error_handler("lidia_signal::uninstall_handler",
							    "Inconsistent signal stack.");
				}

				// remove handler from top of the stack
				--signals[i].sp;

				// stack contains no lidia signal handler anymore:
				//  remove user handler from stack
				if (signals[i].sp == 0) {
					--signals[i].sp;
				}
			}
		}
	}
#endif	// HAVE_POSIX_SIGNALS
}



lidia_signal::lidia_signal (int sig, lidia_signal_handler_t handler)
{
	signal_num = sig;
	lidia_signal::install_handler(sig, handler);
}



lidia_signal::~lidia_signal ()
{
	lidia_signal::uninstall_handler(signal_num);
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
