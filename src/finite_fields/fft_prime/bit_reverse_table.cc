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
//	Author	: Thorsten Rottschaefer (TR)
//                Victor Shoup (VS) and Thomas Pfahler (TPf)
//	Changes	: See CVS log
//
//==============================================================================================


#ifdef HAVE_CONFIG_H
# include	"config.h"
#endif
#include	"LiDIA/finite_fields/bit_reverse_table.h"



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



bit_reverse_table::~bit_reverse_table()
	// destructor deletes the table mem[]
{
	debug_handler("bit_reverse_table", "destructor");
	if (mem != 0) {
		for (lidia_size_t i = 0; i < allocated; i++)
			delete[] mem[i];
	}
	delete[] mem;
}



lidia_size_t bit_reverse_table::rev_inc(lidia_size_t a, lidia_size_t k)
	//increases 'a' in "bitreverse order" (size: 'k' bits)
{
	debug_handler("bit_reverse_table", "rev_inc(lidia_size_t, lidia_size_t)");
	lidia_size_t j, m;
	j = k; m = 1 << (k-1);
	while (j && (m & a)) {
		a ^= m;
		m >>= 1;
		j--;
    	}
	if (j)
		a ^= m;
        return a;
}



lidia_size_t* bit_reverse_table::table(lidia_size_t k)
	//returns mem[k], initializes if mem[k]==0 or k>allocated
{
	debug_handler("fft_table::bit_reverse_table", "table(lidia_size_t)");
	lidia_size_t *rev, i, j;

	if (k >= allocated) {
		//enlarge and copy old values
		lidia_size_t** new_mem = new lidia_size_t*[k+1];
		if (!k) lidia_error_handler("fft_table::bit_reverse_table", "table(lidia_size_t)::out of memory");

		lidia_size_t **op = mem;
		lidia_size_t **np = new_mem;

		for (i = allocated; i != 0; i--, op++, np++)
			*np = *op;

		for (i = allocated; i <= k; i++, np++)
			*np = 0; //new_mem[i] = 0;

		allocated = k+1;
		delete[] mem;
		mem = new_mem;
	}

	rev = mem[k];
	if (rev == 0) {
		//not initialized yet
		lidia_size_t n = 1 << k;
		rev = mem[k] = new lidia_size_t[n];
		if (!rev) lidia_error_handler("fft_table::bit_reverse_table", "table(lidia_size_t)::out of memory");

		lidia_size_t *p = rev;
		for (i = n, j = 0; i != 0; i--, p++, j = rev_inc(j, k))
			*p = j;
	}
	return rev;
}



void bit_reverse_table::copy(udigit* A, const udigit* a, lidia_size_t k)
{
	debug_handler("bit_reverse_table", "copy(udigit*, udigit*, lidia_size_t)");

	switch(k) {
	case 0: A[0] = a[0]; break;
	case 1: A[0] = a[0]; A[1] = a[1]; break;
	case 2: A[0] = a[0]; A[2] = a[1]; A[1] = a[2]; A[3] = a[3]; break;
	case 3: A[0] = a[0]; A[4] = a[1]; A[2] = a[2]; A[6] = a[3];
		A[1] = a[4]; A[5] = a[5]; A[3] = a[6]; A[7] = a[7]; break;
	case 4: A[0] = a[0]; A[8] = a[1]; A[4] = a[2]; A[12] = a[3];
		A[2] = a[4]; A[10] = a[5]; A[6] = a[6]; A[14] = a[7];
		A[1] = a[8]; A[9] = a[9]; A[5] = a[10]; A[13] = a[11];
		A[3] = a[12]; A[11] = a[13]; A[7] = a[14]; A[15] = a[15]; break;
	case 5: A[0] = a[0]; A[16] = a[1]; A[8] = a[2]; A[24] = a[3];
		A[4] = a[4]; A[20] = a[5]; A[12] = a[6]; A[28] = a[7];
		A[2] = a[8]; A[18] = a[9]; A[10] = a[10]; A[26] = a[11];
		A[6] = a[12]; A[22] = a[13]; A[14] = a[14]; A[30] = a[15];
		A[1] = a[16]; A[17] = a[17]; A[9] = a[18]; A[25] = a[19];
		A[5] = a[20]; A[21] = a[21]; A[13] = a[22]; A[29] = a[23];
		A[3] = a[24]; A[19] = a[25]; A[11] = a[26]; A[27] = a[27];
		A[7] = a[28]; A[23] = a[29]; A[15] = a[30]; A[31] = a[31]; break;
	case 6: A[0] = a[0]; A[32] = a[1]; A[16] = a[2]; A[48] = a[3];
		A[8] = a[4]; A[40] = a[5]; A[24] = a[6]; A[56] = a[7];
		A[4] = a[8]; A[36] = a[9]; A[20] = a[10]; A[52] = a[11];
		A[12] = a[12]; A[44] = a[13]; A[28] = a[14]; A[60] = a[15];
		A[2] = a[16]; A[34] = a[17]; A[18] = a[18]; A[50] = a[19];
		A[10] = a[20]; A[42] = a[21]; A[26] = a[22]; A[58] = a[23];
		A[6] = a[24]; A[38] = a[25]; A[22] = a[26]; A[54] = a[27];
		A[14] = a[28]; A[46] = a[29]; A[30] = a[30]; A[62] = a[31];
		A[1] = a[32]; A[33] = a[33]; A[17] = a[34]; A[49] = a[35];
		A[9] = a[36]; A[41] = a[37]; A[25] = a[38]; A[57] = a[39];
		A[5] = a[40]; A[37] = a[41]; A[21] = a[42]; A[53] = a[43];
		A[13] = a[44]; A[45] = a[45]; A[29] = a[46]; A[61] = a[47];
		A[3] = a[48]; A[35] = a[49]; A[19] = a[50]; A[51] = a[51];
		A[11] = a[52]; A[43] = a[53]; A[27] = a[54]; A[59] = a[55];
		A[7] = a[56]; A[39] = a[57]; A[23] = a[58]; A[55] = a[59];
		A[15] = a[60]; A[47] = a[61]; A[31] = a[62]; A[63] = a[63]; break;

	default :
		lidia_size_t n, i;
		const lidia_size_t *rev1, *rev2;
		const udigit *a2;
		n = (1 << (k-1));
		rev1 = table(k);
		rev2 = rev1 + n;
		a2 = a + n;

		for (i = n; i != 0; i--, a++, a2++, rev1++, rev2++) {
			A[*rev1] = *a;
			A[*rev2] = *a2;
		}
	}
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
