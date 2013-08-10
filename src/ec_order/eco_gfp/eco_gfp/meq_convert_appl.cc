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
//	Author	: Markus Maurer, Volker Mueller
//	Changes	: See CVS log
//
//==============================================================================================



#ifdef HAVE_CONFIG_H
# include	"config.h"
#endif
#include	"LiDIA/bigint.h"
#include	<cstdio>
#include	<cstdlib>



#ifdef LIDIA_NAMESPACE
using namespace LiDIA;
#endif



// --------------------------------------------------------------------
//  return true iff file with given input name exists

bool fexist(const char *fname)
{
	FILE *fp;

	if ((fp = fopen (fname, "r")) != NULL) {
		fclose (fp);
		return true;
	}
	else
		return false;
}



//----------------------------------------------------------------------
//  read ascii file 'name1' and write it to binary file 'name2'

void ascii_to_binary(const char * name1, const char * name2)
{
	FILE *ifp, *ofp;
	char c;
	bigint h;

	ifp = fopen(name1, "r");
	ofp = fopen(name2, "wb");

	if (ifp == NULL || ofp == NULL)
		lidia_error_handler("meq_convert", "can open files");

	while (true) {
		do {
			c = fgetc(ifp);
		}
		while (c == ' ' && c != EOF);
		if (c == EOF)
			break;
		if (c == '\n')
			putc(c, ofp);
		else ungetc(c, ifp);
		h.scan_from_file(ifp);
		h.write_to_file (ofp);
	}
	fclose(ifp); fclose(ofp);
}



//------------------------------------------------------------------
// read binary file 'name1' and write it to ascii file 'name2'


void binary_to_ascii(const char *name1, const char * name2)
{
	FILE *ifp, *ofp;
	char c;
	bigint h;

	ifp = fopen(name1, "rb");
	ofp = fopen(name2, "w");

	if (ifp == NULL || ofp == NULL)
		lidia_error_handler("meq_convert", "can open files");

	while (true) {
		do {
			c = fgetc(ifp);
		}
		while (c == ' ' && c != EOF);
		if (c == EOF)
			break;
		if (c == '\n')
			putc(c, ofp);
		else
			ungetc(c, ifp);
		h.read_from_file(ifp);
		h.print_to_file (ofp);
		putc(' ', ofp);
	}
	fclose(ifp); fclose(ofp);
}



//**********************************************************************
// main file which does the administration.

#if 0
char* formats[] = {"ASCII", "BINARY", "GZIPPED ASCII",
                   "GZIPPED BINARY"};
#endif



int main_LiDIA(int argc, char** argv)
{
	unsigned int l;
	unsigned int n = 0, m = 0;
	bool aexist = false, bexist = false;
	int type = 0;

	std::cout << "\nConversion Program for Modular Equations";
	std::cout << "\n========================================\n\n";

	std::cout << "\nPlease input prime l : "; std::cin >> l;


	char name[2][1024], bname[2][1024];
	char command[1024];


	if (n == 0) {
		if (fexist(name[1]))
			type = 1;

		if (m == 1 || m == 3) {
			if (!bexist)
				ascii_to_binary(name[type], bname[type]);
			if (m == 3) {
				sprintf(command, "gzip -9 %s", bname[type]);
				std::system(command);
				if (!bexist)
					std::remove(bname[type]);
			}
			std::cout << " done \n\n" << std::flush;
			return 0;
		}
		else {
			sprintf(command, "gzip -9 %s", name[type]);
			std::system(command);
			std::cout << " done \n\n" << std::flush;
			return 0;
		}
	}

	if (n == 1) {
		if (fexist(bname[1]))
			type = 1;

		if (m == 0 || m == 2) {
			if (!aexist)
				binary_to_ascii(bname[type], name[type]);
			if (m == 2) {
				sprintf(command, "gzip -9 %s", name[type]);
				std::system(command);
				if (!aexist)
					std::remove(name[type]);
			}
			std::cout << " done \n\n" << std::flush;
			return 0;
		}
		else {
			sprintf(command, "gzip -9 %s", bname[type]);
			std::system(command);
			std::cout << " done \n\n" << std::flush;
			return 0;
		}
	}
#if 0
	if (n == 2) {
		if (fexist(gname[1]))
			type = 1;

		if (!aexist) {
			sprintf(command, "gunzip %s", gname[type]);
			std::system(command);
		}
		if (m == 1 || m == 3) {
			if (!bexist)
				ascii_to_binary(name[type], bname[type]);

			if (m == 3) {
				sprintf(command, "gzip -9 %s", bname[type]);
				std::system(command);
				if (!bexist)
					std::remove(bname[type]);
			}
			std::cout << " done \n\n" << std::flush;
			return 0;
		}
		if (m != 0 && !aexist)
			std::remove(name[type]);
		std::cout << " done \n\n" << std::flush;
		return 0;
	}

	if (n == 3) {
		if (fexist(gbname[1]))
			type = 1;

		if (!bexist) {
			sprintf(command, "gunzip %s", gbname[type]);
			std::system(command);
		}
		if (m == 0 || m == 2) {
			if (!aexist)
				binary_to_ascii(bname[type], name[type]);

			if (m == 2) {
				sprintf(command, "gzip -9 %s", name[type]);
				std::system(command);
				if (!aexist)
					std::remove(name[type]);
			}
			std::cout << " done \n\n" << std::flush;
			return 0;
		}
		if (m != 1 && !bexist)
			std::remove(bname[type]);
		std::cout << " done \n\n" << std::flush;
		return 0;
	}
#endif
}


int main(int argc, char** argv) {

#if defined(LIDIA_EXCEPTIONS)
    try {
#endif

	main_LiDIA(argc, argv);
	
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
