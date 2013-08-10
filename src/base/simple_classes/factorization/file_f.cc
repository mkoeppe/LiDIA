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

// Description : this file contains several functions for name handling
//               and file handling. These are used in the QS
//               implementation.


#ifdef HAVE_CONFIG_H
# include	"config.h"
#endif
#include	"LiDIA/LiDIA.h"
#include	<cstdlib>
#include	<cstring>
#include	<cstdio>
#include	<unistd.h>
#include        <errno.h>


#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



//***********************************************************

static void
itoa (unsigned long n, char s[])
{
	int c, i, j;

	i = 0;
	do {
		s[i++] = static_cast<char>(n % 10 + '0');
	} while ((n /= 10) > 0);
	s[i] = '\0';

	for (i = 0, j = strlen(s) - 1; i < j; i++, j--) {
		c = *(s + i);
		*(s + i) = *(s + j);
		*(s + j) = c;
	}
}



//***********************************************************
// returns a  pointer to the name, that must not be deleted.
//

const int size_of_db = 20;

char*
get_name (const char *s)
{
  static char *string[size_of_db], *filename[size_of_db];
  static int initialized = 0;
  char *tmppath;
  int i;

  if (!initialized)
    {
      initialized = 1;
      for (i = 0; i < size_of_db; i++)
	string[i] = NULL;
    }
  
  for (i = 0; i < size_of_db; i++)
    {
      if (string[i] == NULL)
	break;

      if (strcmp(string[i], s) == 0)
	return filename[i];
    }
  
  if (i >= size_of_db)
    lidia_error_handler("rational_factorization","get_name::no new "
			"names available");
  
  string[i] = new char[strlen(s) + 1];
  strcpy(string[i], s);

  filename[i] = new char[1000];

  tmppath = getenv("TMPDIR");
  if(tmppath != NULL)
  {
    strcpy(filename[i], tmppath);
    strcat(filename[i],"/");
    strcat(filename[i], s);
  }
  else
  {
#ifdef P_tmpdir
    strcpy(filename[i], P_tmpdir);
    strcat(filename[i],"/");
    strcat(filename[i], s);
#else
    strcpy(filenmae[i], s);
#endif
  }

#ifdef HAVE_MKSTEMP
  strcat(filename[i], "XXXXXX");
  int fh = mkstemp(filename[i]);
  if (fh == -1)
    lidia_error_handler("mpqs", "get_name::can't generate temporary file");
  int errorcode;
  do { // loop in case close() is interrupted by a signal
    if (close(fh) != 0) {
      errorcode = errno;
    }
    else {
      errorcode = 0;
    }
  } while(errorcode == EINTR);
  if (errorcode != 0) {
    lidia_error_handler("mpqs", "get_name::can't close temporary file");
  }
#else
  if (tmpnam(filename[i]) == NULL)
    lidia_error_handler("mpqs", "get_name::can't generate temporary file");
#endif
  return filename[i];
}



//**************************************************************
// now we implement the file handling stuff.
//
// Note that BUFSIZE and STLEN influence the main memory needed
// in the sorting step. If you lower these values, then the routines
// are slightly slower, but need less memory.
//**************************************************************

const size_t BUFSIZE = 512*1024;
const int STLEN = 5000;
const int STLEN2 = 1000;


//*************************************************************
// concatenate the two input files with names f1, f2 and write
// result to f3. Note that the files have to be pairwise different.
//*************************************************************

void
concat (char* f1, char* f2, char* f3)
{
	FILE * fp1, *fp2, *fp3;
	char buffer[BUFSIZE+1];
	size_t i;

	fp3 = fopen(f3, "w");

	if (fp3 == NULL) {
		lidia_error_handler("rational_factorization",
				    "concat::can't open output file");
		return;
	}


	fp1 = fopen(f1, "r");

	if (fp1 != NULL) {
		while ((i = fread((char*) buffer, sizeof(char), BUFSIZE, fp1)) == BUFSIZE) {
			fwrite((char *) buffer, sizeof(char), BUFSIZE, fp3);
		}

		fwrite((char *) buffer, sizeof(char), i, fp3);
		fclose (fp1);
	}


	fp2 = fopen(f2, "r");

	if (fp2 != NULL) {
		while ((i = fread((char*) buffer, sizeof(char), BUFSIZE, fp2)) == BUFSIZE) {
			fwrite((char *) buffer, sizeof(char), BUFSIZE, fp3);
		}

		fwrite((char *) buffer, sizeof(char), i, fp3);
		fclose (fp2);
	}

	fclose(fp3);
}



//********************************************************************
// copy input file to output file (both given by file pointers)
//********************************************************************

inline void
copy_file (FILE *fpin, FILE *fpout)
{
	char buffer[BUFSIZE+1];
	size_t i;

	rewind(fpout);
	rewind(fpin);

	while ((i = fread((char*) buffer, sizeof(char), BUFSIZE, fpin)) == BUFSIZE) {
		fwrite((char *) buffer, sizeof(char), BUFSIZE, fpout);
	}

	fwrite((char *) buffer, sizeof(char), i, fpout);
	fflush(fpout);
}



//**************************************************************
// merge the two input files ifbs_{1, 2} together and write the result
// to the output file, delete multiple rows.
//
// Note: all file names have to be pairwise different.
//**************************************************************

void
merge_unique (char *ifbs1, char* ifbs2, char *ofbs)
{
	FILE *ifb1, *ifb2, *ofb;
	char buffer1[STLEN2];
	char buffer2[STLEN2];
	char buffer3[STLEN2];
	char *s1, *s2;
	int r;

	ofb = fopen(ofbs, "w");

	if (ofb == NULL) {
		lidia_error_handler("rational_factorization",
				    "merge_unique::can't open output file");
		return;
	}

	ifb1 = fopen(ifbs1, "r");
	ifb2 = fopen(ifbs2, "r");

	if (ifb1 == NULL && ifb2 == NULL) {
		fclose (ofb);
		return;
	}

	if (ifb1 == NULL) {
		copy_file(ifb2, ofb);
		fclose(ifb2); fclose(ofb);
		return;
	}

	if (ifb2 == NULL) {
		copy_file(ifb1, ofb);
		fclose(ifb1); fclose(ofb);
		return;
	}


	s1 = fgets(buffer1, STLEN2, ifb1);
	s2 = fgets(buffer2, STLEN2, ifb2);

	while (s1 != NULL && s2 != NULL)
	{
		r = strcmp(buffer1, buffer2);

		if (r < 0) {
			fputs(buffer1, ofb);
			strcpy(buffer3, buffer1);
			do {
				s1 = fgets(buffer1, STLEN2, ifb1);
			} while (s1 != NULL && strcmp(buffer1, buffer3) == 0);
		}
		else if (r > 0) {
			fputs(buffer2, ofb);
			strcpy(buffer3, buffer2);
			do {
				s2 = fgets(buffer2, STLEN2, ifb2);
			} while (s2 != NULL && strcmp(buffer2, buffer3) == 0);
		}
		else {
			s2 = fgets(buffer2, STLEN2, ifb2);
		}
	}

	if (s1 == NULL) {
		while (s2 != NULL) {
			fputs(buffer2, ofb);
			strcpy(buffer3, buffer2);
			do {
				s2 = fgets(buffer2, STLEN2, ifb2);
			} while (s2 != NULL && strcmp(buffer2, buffer3) == 0);
		}
	}
	else {
		while (s1 != NULL) {
			fputs(buffer1, ofb);
			strcpy(buffer3, buffer1);
			do {
				s1 = fgets(buffer1, STLEN2, ifb1);
			} while (s1 != NULL && strcmp(buffer1, buffer3) == 0);
		}
	}

	fclose(ifb1); fclose(ifb2); fclose(ofb);
}



//**********************************************************************
// we assume that array sbuf[0, ..., size-1] and input file are sorted,
// merge two structures uniquely (identical rows are deleted) and write
// result to output file
//**********************************************************************

inline void
merge_array_file (char ** sbuf, int size,
		  FILE *ifb, FILE *ofb)
{
	char buffer2[STLEN2];
	char buffer3[STLEN2];
	char *s2;
	int r, index = 0;

	rewind(ifb);
	rewind(ofb);
	s2 = fgets(buffer2, STLEN2, ifb);

	do {
		r = strcmp(sbuf[index], buffer2);

		if (r < 0) {
			fputs(sbuf[index], ofb);
			strcpy(buffer3, sbuf[index]);
			do {
				if (index < size-1)
					index ++;
			} while (index < size-1 && strcmp(sbuf[index], buffer3) == 0);

			if (index == size-1 && strcmp(sbuf[index], buffer3) == 0) {
				index ++;
			}
		}
		else {
			if (r > 0) {
				fputs(buffer2, ofb);
				strcpy(buffer3, buffer2);
				do {
					s2 = fgets(buffer2, STLEN2, ifb);
				} while (s2 != NULL && strcmp(buffer2, buffer3) == 0);
			}
			else {
				s2 = fgets(buffer2, STLEN2, ifb);
			}
		}
	} while (index < size && s2 != NULL);

	if (index == size) {
		// remaining rows in file
		do {
			fputs(buffer2, ofb);
			strcpy(buffer3, buffer2);
			do {
				s2 = fgets(buffer2, STLEN2, ifb);
			} while (s2 != NULL && strcmp(buffer2, buffer3) == 0);
		} while (s2 != NULL);
	}
	else {
		// remaining rows in array
		do {
			fputs(sbuf[index], ofb);
			strcpy(buffer3, sbuf[index]);
			do {
				if (index < size-1)
					index ++;
			} while (index < size-1 && strcmp(sbuf[index], buffer3) == 0);
			if (index == size-1 && strcmp(sbuf[index], buffer3) != 0) {
				fputs(sbuf[index], ofb);
			}
		} while (index < size-1);
	}
}



//***************************************************************
// given two rows, compare rows by large prime (first unsigned long
// in row)
//***************************************************************

extern "C" int
cmp_numeric (const void *sbuff1, const void *sbuff2)
{
	const char *s1, *s2;
	long l1, l2;

	s1 = *((char* const *) sbuff1);
	s2 = *((char* const *) sbuff2);
	l1 = atol(s1);
	l2 = atol(s2);

	if (l1 < l2)
		return -1;
	if (l1 > l2)
		return 1;
	return 0;
}



//***************************************************************
// given two rows, compare complete rows
//***************************************************************

extern "C" int
cmp_row (const void *sbuff1, const void *sbuff2)
{
	const char *s1, *s2;

	s1 = *((char* const *) sbuff1);
	s2 = *((char* const *) sbuff2);
	return (strcmp(s1, s2));
}



//***************************************************************
// sort input file by large primes, delete multiple lines
//***************************************************************

void
sort_file_lp (FILE* ifb1, char *ofbs)
{
	FILE *ofb;
	int i, j;
	int round = 0;
	FILE *tmp[2];

	// precheck 
	if(ifb1==NULL)
	{
	  lidia_error_handler("rational_factorization", "sort_file_numeric::can't open file fb1");
	  return;
	}


	ofb = fopen(ofbs, "w");
	if(ofb==NULL)
	{
	  lidia_error_handler("rational_factorization", "sort_file_numeric::can't open file ofbs");
	  return;
	}
	//char* name[2];
	//name[0] = new char[STLEN2];
	//name[1] = new char[STLEN2];

	//strcpy(name[0],"XXXXXX");
	//strcpy(name[1],"XXXXXX");

	//int fd = mkstemp(name[0]);
	//tmp[0] = fdopen(fd, "w+");
	//tmp[0] = tmpfile();
     tmp[0] = fopen(get_name("SORT_LP0"),"w+");
	 if(tmp[0] == NULL)
	 {
		fclose(ofb);
		lidia_error_handler("rational_factorization", "sort_file_numeric::can't open file SORT_LP0");
		return;
	 }

	//fd = mkstemp(name[1]);
	//tmp[1] = fdopen(fd, "w+");
	//tmp[1] = tmpfile();
        tmp[1] = fopen(get_name("SORT_LP1"),"w+");
	 if(tmp[1] == NULL)
	 {
		fclose(tmp[0]);
		fclose(ofb);
		lidia_error_handler("rational_factorization", "sort_file_numeric::can't open file SORT_LP1");
		return;
	 }

	char **sbuff, *c1;

	sbuff = new char*[STLEN];
	if (sbuff == NULL) {
		lidia_error_handler("rational_factorization", "sort_file_numeric::can't allocate memory");
		fclose(tmp[0]);
		fclose(tmp[1]);
		fclose(ofb);
		return;
	}

	for (i = 0; i < STLEN; i++) {
		sbuff[i] = new char[STLEN2];
		if (sbuff[i] == NULL) {
			lidia_error_handler("rational_factorization", "sort_file_numeric::can't allocate memory");
			fclose(tmp[0]);
			fclose(tmp[1]);
			fclose(ofb);
			delete [] sbuff;
			return;
		}
	}

	do {
		i = 0;
		do {
			c1 = fgets((char*) sbuff[i++], STLEN2, ifb1);
		} while (c1 != NULL && i < static_cast<int>(STLEN));

		if (i != STLEN)
			i--;

		qsort(static_cast<void *>(sbuff), i, sizeof(char **), cmp_numeric);

		if (round == 0) {
			for (j = 0; j < i; j++) {
				//look out: multiple rows !!
				fputs(sbuff[j], tmp[1]);
			}

			round ++;
			fflush(tmp[1]);
		}
		else {
			fflush(tmp[round%2]);
			merge_array_file(sbuff, i, tmp[round % 2], tmp[(round+1)%2]);
			round ++;
		}
	} while (c1 != NULL);

	for (i = 0; i < STLEN; i++) {
		// delete memory
		delete[] sbuff[i];
	}

	delete[] sbuff;

	copy_file(tmp[round%2], ofb);
	fclose(tmp[0]); fclose(tmp[1]);
        std::remove(get_name("SORT_LP0"));
        std::remove(get_name("SORT_LP1"));
	//std::remove(name[0]);
	//std::remove(name[1]);
	fclose(ofb);
	//delete[] name[0];
	//delete[] name[1];
}



//*************************************************************
// sort input file by rows, delete multiple lines
//*************************************************************

void
sort_file_line (char* ifbs1, char *ofbs)
{
	FILE *ifb1, *ofb;
	int i, j;
	int round = 0;
	FILE *tmp[2];

	ifb1 = fopen(ifbs1, "r");
	ofb = fopen(ofbs, "w");

	char* name[2];
	//name[0] = new char[STLEN2];
	//name[1] = new char[STLEN2];
	//strcpy(name[0],"XXXXXX");
	//strcpy(name[1],"XXXXXX");

	//int fd = mkstemp(name[0]);
	//tmp[0] = fdopen(fd, "w+");
	//tmp[0] = tmpfile();
        name[0] = get_name("SORT_LINE0");
        tmp[0] = fopen(name[0],"w+");

	//fd = mkstemp(name[1]);
	//tmp[1] = fdopen(fd, "w+");
	//tmp[1] = tmpfile();
        name[1] = get_name("SORT_LINE1");
        tmp[1] = fopen(name[1],"w+");

	if (ifb1 == NULL || ofb == NULL || tmp[1] == NULL) {
		lidia_error_handler("rational_factorization",
				    "sort_file_line::can't open file");
		return;
	}

	char **sbuff, *c1;

	sbuff = new char*[STLEN];

	if (sbuff == NULL) {
		lidia_error_handler("rational_factorization",
				    "sort_file_line::can't allocate memory");
		return;
	}

	for (i = 0; i < STLEN; i++) {
		sbuff[i] = new char[STLEN2];
		if (sbuff[i] == NULL) {
			lidia_error_handler("rational_factorization",
					    "sort_file_line::can't allocate memory");
			return;
		}
	}

	do {
		i = 0;
		do {
			c1 = fgets((char*) sbuff[i++], STLEN2, ifb1);
		} while (c1 != NULL && i < static_cast<int>(STLEN));

		if (i != STLEN)
			i--;

		qsort(static_cast<void *>(sbuff), i, sizeof(char **), cmp_row);

		if (round == 0) {
			fputs(sbuff[0], tmp[1]);
			for (j = 1; j < i; j++) {
				if (strcmp(sbuff[j], sbuff[j-1]) != 0) {
					fputs(sbuff[j], tmp[1]);
				}
			}

			round ++;
			fflush(tmp[1]);
		}
		else {
			//int fd;
			fclose(tmp[(round+1)&1]);
			std::remove (name[(round+1)&1]);
			//strcpy(name[(round+1) & 1], "XXXXXX");
			//fd = mkstemp(name[(round+1)&1]);
			//tmp[(round+1)&1] = fdopen(fd,"w+");
			//tmp[(round+1)&1] = tmpfile();
			tmp[(round+1)&1] = fopen(name[(round+1)&1],"w+");
			if (tmp[(round+1)&1] == NULL) {
				lidia_error_handler("rational_factorization",
						    "sort_file_line::can't open tmp file");
				return;
			}
			merge_array_file(sbuff, i, tmp[round &1], tmp[(round+1)&1]);
			fflush(tmp[(round+1)&1]);
			round ++;
		}
	} while (c1 != NULL); // first file done

	for (i = 0; i < STLEN; i++) {
		// delete memory
		delete[] sbuff[i];
	}

	delete[] sbuff;

	copy_file(tmp[round%2], ofb);

	fclose(tmp[0]);
	fclose(tmp[1]);
	std::remove(name[0]);
	std::remove(name[1]);
	fclose(ifb1);
	fclose(ofb);
	//delete[] name[0];
	//delete[] name[1];
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
