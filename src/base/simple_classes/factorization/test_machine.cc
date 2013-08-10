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
//	Author	: Volker Mueller (VM)
//	Changes	: See CVS log
//
//==============================================================================================


#include	"LiDIA/timer.h"
#include        <fstream.h>


#ifdef LIDIA_NAMESPACE
using namespace LiDIA;
#endif



typedef unsigned char SIEBTYP;



int main()
{
	static unsigned int prime[54] = { 11, 13, 15, 17, 19, 21, 23, 29, 31, 33,
					  35, 37, 39, 41, 43, 47, 51, 53, 55, 57,
					  59, 61, 65, 67, 71, 73, 79, 83, 85, 89,
					  91, 93, 95, 97, 101, 103, 105, 107, 109,
					  113, 127, 131, 149, 151, 157, 163, 179,
					  239, 397, 401, 439, 479, 641};
	timer test;
	unsigned int p, treffer = 0;
	lidia_size_t x, j, jj;
	SIEBTYP *sieb;
	double mul;

	test.start_timer();
	sieb = new SIEBTYP[100000];

	for (jj = 0; jj < 7; jj++)
		for (j = 0; j < 50; j++) {
			p = prime[j];
			x = 0;
			while (x < 100000) {
				sieb[x] += 3;
				x += p;
			}
		}

	x = 0;
	while (x < 100000) {
		if (sieb[x] > 30)
			treffer++;
		x++;
	}
	delete[] sieb;
	test.stop_timer();
	mul = (test.user_time()+test.sys_time())/63.0; // 63 is the reference time
	//  on a SPARC ELC

	std::ofstream fp("mpqs_timing.h");
	fp << "const double mult_machine = " << mul << ";" << std::endl;
        fp.close();
	return 0;
}
