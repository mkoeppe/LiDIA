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
//	Author	: Thomas Pfahler (TPf), Volker Mueller (VM)
//	Changes	: See CVS log
//
//==============================================================================================


#ifdef HAVE_CONFIG_H
# include	"config.h"
#endif
#include	"LiDIA/udigit.h"



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



#if !INLINE_INTERFACE
# define inline
# include	"LiDIA/kernel/udigit_interface.h"
# undef inline
#endif



udigit
negate_mod (udigit a, udigit p)
	// a must be less than p, p > 1
{
	if (a == 0) {
		return 0;
	}
	else {
		return (p-a);
	}
}



// ============================================================
// COMPUTE a^e, WITHOUT TAKING CARE OF CARRYS               
// ============================================================

udigit
power (udigit a, unsigned int e)
{
	if (a == 1) {
		return a;
	}

	udigit n = 1;

	while (e > 0) {
		n *= a;
		e--;
	}
	return n;
}



// ============================================================
// INVERT a MOD m, BINARY GCD ALGORITHM, ONLY COEFFICIENT FOR 
// a IS COMPUTED                                              
// ============================================================

udigit
invert_mod (udigit a, udigit m)
{
	udigit	xu, yu, ui, ui1, v, u;
	int	single_two = 0;


	ui = a % m;

	if (m == 0 || ui == 0) {
		lidia_error_handler("udigit", "invert_mod::element not invertible");
	}

	ui1 = m;

	if ((!(ui & 1)) && (!(ui1 & 1))) {
		lidia_error_handler("udigit", "invert_mod::gcd(a, m) > 1");
	}

	while (true) {
		if (ui & 1) {
			break;
		}
		ui >>= 1;
		single_two++;
	}

	while (true) {
		if (ui1 & 1) {
			break;
		}
		ui1 >>= 1;
		single_two--;
	}

	v = ui1;
	u = ui;
	xu = 1;
	yu = v;

	int flagui = 1, flagui1 = 1;
	udigit vh = (v >> 1);


	// Invariants:    ui = flagui * ( xu * u - ?? * v)
	// ui1 = flagui1 * (yu * u - ?? * v)
	// u, v odd


	while (true) {
		if (ui == ui1 || ui == 1) {
			break;
		}

		if (ui1 == 1) {
			ui = 1;
			xu = yu;
			flagui = flagui1;
			break;
		}
		if (ui > ui1) {
			ui -= ui1;
			ui >>= 1;

			if (flagui == flagui1) {
				if (xu >= yu) {
					xu -= yu;
				}
				else {
					flagui = -flagui;
					xu = yu - xu;
				}

				if (xu & 1) {
					xu = (xu >> 1) + vh + 1;
				}
				else {
					xu >>= 1;
				}
			}
			else {
				// flagui != flagui1
				if ((xu & 1) != (yu & 1)) {
					// xu + yu odd
					if (xu >= (v-yu)) {
						// xu + yu >= v
						xu = (xu - (v-yu)) >> 1;
					}
					else {
						xu = ((xu + yu) >> 1) + vh + 1;
					}
				}
				else {
					xu = (xu >> 1) + (yu >> 1) + (yu & 1);
				}
			}

			while (true) {
				// drop 2's in ui
				if (ui & 1) {
					break;
				}
				ui >>= 1;

				if (!(xu & 1)) {
					xu >>= 1;
				}
				else {
					xu = (xu >> 1) + vh + 1;
				}
			}
		}
		else {
			// ui1 > ui
			ui1 -= ui;
			ui1 >>= 1;
			if (flagui == flagui1) {
				if (yu >= xu) {
					yu -= xu;
				}
				else {
					flagui1 = -flagui1;
					yu = xu - yu;
				}
				if (yu & 1) {
					yu = (yu >> 1) + vh + 1;
				}
				else {
					yu >>= 1;
				}
			}
			else {
				// flagui != flagui1
				if ((xu & 1) != (yu & 1)) {
					// xu + yu odd
					if (yu >= (v-xu)) {
						// xu + yu >= v
						yu = (yu - (v-xu)) >> 1;
					}
					else {
						yu = ((xu + yu) >> 1) + vh + 1;
					}
				}
				else {
					yu = (xu >> 1) + (yu >> 1) + (yu & 1);
				}
			}

			while (true) {
				// drop 2's in ui1
				if (ui1 & 1) {
					break;
				}
				ui1 >>= 1;

				if (!(yu & 1)) {
					yu >>= 1;
				}
				else {
					yu = (yu >> 1) + vh + 1;
				}
			}
		} // else "ui1 > ui"
	} // while (true)

	if (ui != 1) {
		lidia_error_handler("udigit", "invert_mod::gcd(a, m) > 1");
	}

	if (!single_two) {
		if (flagui == 1) {
			return xu;
		}
		else {
			return (m-xu);
		}
	}

	if (single_two > 0) {
		// 2's in a
		while (single_two > 0) {
			single_two --;
			if (!(xu & 1)) {
				xu >>= 1;
			}
			else {
				xu = (xu >> 1) + vh + 1;
			}
		}
		if (flagui == 1) {
			return xu;
		}
		else {
			return (m-xu);
		}
	}
	else {
		// 2's in m
		single_two ++;
		if (!(xu & 1)) {
			xu += v;
		}
		v <<= 1;
		vh = 3;

		while (single_two < 0) {
			single_two ++;

			if (((xu & vh) * (u & vh) - 1) & vh) {
				xu += v;
			}

			v <<= 1;
			vh = (vh << 1) + 1;
		}
		if (flagui == 1) {
			return xu;
		}
		else {
			return (m-xu);
		}
	}
}



// ============================================================
// BINARY GCD ALGORITHM                                       
// ============================================================

udigit
gcd (udigit a, udigit b)
{
	unsigned int	commontwo = 0;


	while (!(a & 1) && !(b & 1)) {
		commontwo++;
		a >>= 1;
		b >>= 1;
	}

	while (!(a & 1)) {
		a >>= 1;
	}

	while (!(b & 1)) {
		b >>= 1;
	}

	while (true) {
		if (a >= b) {
			a -= b;
			if (!a) {
				break;
			}
			while (!(a & 1)) {
				a >>= 1;
			}
		}
		else {
			b -= a;
			while (!(b & 1)) {
				b >>= 1;
			}
		}
	}
	return (b << commontwo);
}



// ============================================================
// EXTENDED BINARY GCD ALGORITHM, NOTE THAT WE TAKE CARE OF   
// OVERFLOWS                                                  
// ============================================================

udigit
xgcd (udigit & u, udigit & v, udigit a, udigit b)
{
	udigit		xu, yu, xv, yv, ui, ui1;
	unsigned int	common_two = 0;
	int		single_two = 0;


	if (b == 0) {
		v = 0;
		u = 1;
		return a;
	}

	if (a == 0) {
		u = 0;
		v = 1;
		return b;
	}

	ui = a;
	ui1 = b;

	while (true) {
		// drop common 2's in a and b
		if ((ui & 1) || (ui1 & 1)) {
			break;
		}
		common_two++;
		ui >>= 1;
		ui1 >>= 1;
	}

	while (true) {
		if (ui & 1) {
			break;
		}
		ui >>= 1;
		single_two++;
	}

	while (true) {
		if (ui1 & 1) {
			break;
		}
		ui1 >>= 1;
		single_two--;
	}

	u = ui;
	v = ui1;
	xu = 1; xv = 0;
	yu = v; yv = u-1;
	int flagui = 1, flagui1 = 1;

	//  Invariants:  ui   =  flagui * ( xu * u - xv * v )
	//  ui1 = flagui1 * (yu * u - yv * v)
	//  u, v  odd, 0 <= xu, yu <= v, 0 <= xv, yv <= u
	//  see (*) below

	udigit vh = (v >> 1);
	udigit uh = (u >> 1);


	while (true) {
		if (ui == ui1) {
			break;
		}

		if (ui > ui1) {
			ui -= ui1;
			ui >>= 1;

			if (flagui == flagui1) {
				if (xu >= yu && xv >= yv) {
					// (*) TPf, June 2000:
					// This used to be: if (xu >= yu) and produced incorrect results
					// now, the invariants above might not be fulfilled (but the
					// results appear to be correct nevertheless)
					xu -= yu;
					xv -= yv;
				}
				else {
					flagui = -flagui;
					xu = yu - xu;
					xv = yv - xv;
				}

				if (xu & 1) {
					xu = (xu >> 1) + vh + 1;
					xv = (xv >> 1) + uh + 1;
				}
				else {
					xu >>= 1;
					xv >>= 1;
				}
			}
			else {
				// flagui != flagui1
				if ((xu & 1) != (yu & 1)) {
					// xu + yu odd
					if (xu >= (v-yu)) {
						// xu + yu >= v
						xu = (xu - (v-yu)) >> 1;
						xv = (xv - (u-yv)) >> 1;
					}
					else {
						xu = ((xu + yu) >> 1) + vh + 1;
						xv = ((xv + yv) >> 1) + uh + 1;
					}
				}
				else {
					xu = (xu >> 1) + (yu >> 1) + (yu & 1);
					xv = (xv >> 1) + (yv >> 1) + (yv & 1);
				}
			}

			while (true) {
				// drop 2's in ui
				if (ui & 1) {
					break;
				}
				ui >>= 1;

				if (!(xu & 1)) {
					xu >>= 1;
					xv >>= 1;
				}
				else {
					xu = (xu >> 1) + vh + 1;
					xv = (xv >> 1) + uh + 1;
				}
			}
		}
		else {
			// ui1 > ui
			ui1 -= ui;
			ui1 >>= 1;
			if (flagui == flagui1) {
				if (yu >= xu) {
					// and yv >= xv !!
					yu -= xu;
					yv -= xv;
				}
				else {
					flagui1 = -flagui1;
					yu = xu - yu;
					yv = xv - yv;
				}

				if (yu & 1) {
					yu = (yu >> 1) + vh + 1;
					yv = (yv >> 1) + uh + 1;
				}
				else {
					yu >>= 1;
					yv >>= 1;
				}
			}
			else {
				// flagui != flagui1
				if ((xu & 1) != (yu & 1)) {
					// xu + yu odd
					if (yu >= (v-xu)) {
						// xu + yu >= v
						yu = (yu - (v-xu)) >> 1;
						yv = (yv - (u-xv)) >> 1;
					}
					else {
						yu = ((xu + yu) >> 1) + vh + 1;
						yv = ((xv + yv) >> 1) + uh + 1;
					}
				}
				else {
					yu = (xu >> 1) + (yu >> 1) + (yu & 1);
					yv = (xv >> 1) + (yv >> 1) + (yv & 1);
				}
			}

			while (true) {
				// drop 2's in ui1
				if (ui1 & 1) {
					break;
				}
				ui1 >>= 1;

				if (!(yu & 1)) {
					yu >>= 1;
					yv >>= 1;
				}
				else {
					yu = (yu >> 1) + vh + 1;
					yv = (yv >> 1) + uh + 1;
				}
			} // while-loop for dropping 2's in ui1
		} // else "ui1 > ui"
	} // while (true)

	if (flagui < 0) {
		xu = v - xu;
		xv = u - xv;
	}

	if (!single_two) {
		u = xu;
		v = xv;
		return (ui << common_two);
	}

	if (single_two > 0) {
		// single 2's in a
		while (single_two > 0) {
			single_two --;
			if (!(xu & 1)) {
				xu >>= 1;
			}
			else {
				xu = (xu >> 1) + vh + 1;
				xv += u;
			}
			u <<= 1;
		}
		u = xu;
		v = xv;
		return (ui << common_two);
	}
	else {
		while (single_two < 0) {
			single_two ++;
			if (!(xv & 1)) {
				xv >>= 1;
			}
			else {
				xv = (xv + u) >> 1;
				xu += v;
			}
			v <<= 1;
		}
		u = xu;
		v = xv;
		return (ui << common_two);
	}
}



// ============================================================
// COMPUTE base^exp MOD modul                               
// ============================================================

udigit
power_mod (udigit base, long  exp, udigit modul)
{
	udigit	erg, ergz;


	if (modul <= static_cast<udigit>(1)) {
		lidia_error_handler("power_mod", "modul <= 1");
	}

	if (exp == static_cast<udigit>(0)) {
		return 1;
	}

	if (exp < 0) {
		exp = -exp;
		ergz = invert_mod(base, modul);
	}
	else {
		ergz = base % modul;
		if (ergz == 0) {
			return 0;
		}
	}

	erg = static_cast<udigit>(1);

	while (true) {
		if (exp & 1) {
			erg = multiply_mod(erg, ergz, modul);
			if (exp == 1) {
				return erg;
			}
		}
		exp >>= 1;
		ergz = multiply_mod(ergz, ergz, modul);
	}
}



// *********************************************************************

udigit
power_mod (udigit base, udigit exp, udigit modul)
{
	udigit	erg, ergz;


	if (modul <= 1) {
		lidia_error_handler("power_mod", "modul <= 1");
	}

	if (exp == 0) {
		return 1;
	}

	ergz = base % modul;
	if (ergz == 0) {
		return 0;
	}

	erg = static_cast<udigit>(1);

	while (true) {
		if (exp & 1) {
			erg = multiply_mod(erg, ergz, modul);
			if (exp == 1) {
				return erg;
			}
		}
		exp >>= 1;
		ergz = multiply_mod(ergz, ergz, modul);
	}
}



// ***************************************************************

bool
is_prime (udigit n, unsigned int trials)
{
	static udigit	a[10] = {3, 5, 7, 11, 13, 17, 19, 23, 29, 31};

	udigit		b = 0;
	unsigned int	j, i, k = 0;


	if (n == 2) {
		return true;
	}

	if (n <= 1 || !(n & 1)) {
		return false;
	}

	if (n <= 36) {
		for (i = 0; i < 10; i++) {
			if (a[i] == n) {
				return true;
			}
		}
		return false;
	}

	udigit exp = n-1, odd_exp, erg; // Miller-Rabin Test
	odd_exp = exp;

	while (!(odd_exp & 1)) {
		odd_exp >>= 1;
		k++;
	}

	for (i = 1; i <= trials; i++) {
		if (i < 11) {
			b = a[i-1];
		}
		else {
			b = next_prime(b);
		}

		if (b == n) {
			return true;
		}
		erg = power_mod(b, odd_exp, n);

		j = k;
		while ((j > 0) && erg != 1 && erg != exp) {
			erg = multiply_mod(erg, erg, n);
			j--;
		}
		if (j == 0 && erg != 1) {
			return false;
		}
	}
	return true;
}



// ***************************************************************

int
jacobi (udigit a, udigit b)
{
	udigit	h;
	long	k = 1, v;

        // Setting up the table for table-lookup

	static int table[] = {0, 1, 0, -1, 0, -1, 0, 1};


        // test trivial cases

	if (!b) {
		if (a == 1) {
			return 1;
		}
		else {
			return 0;
		}
	}

	if (!(a & 1) && !(b & 1)) {
		return 0;
	}

	v = 0;
	while (!(b & 1)) {
		// b even
		v++;
		b >>= 1;
	}

	if (v & 1) {
		k *= table[a & 7];
	}

	while (a != 0) {
		// main-loop
		v = 0;
		while (!(a & 1)) {
			v++;
			a >>= 1;
		}
		if (v & 1) {
			k *= table[b & 7];
		}

		if (a < b) {
			// swap and correct intermediate result
			h = a;
			a = b;
			b = h;
			if (a & b & 2) {
				k = -k;
			}
		}
		a -= b;
	}

	if (b == 1) {
		return k;
	}
	else {
		return 0;
	}
}



// **********************************************************************

udigit
sqrt_mod (udigit a, udigit p)
	// p is assumed to be prime !!
{
	udigit		r, n, c, z, k;
	udigit		t;
	unsigned int	s;


	if (p == 2) {
		return a;
	}

	// p = 3 mod 4
	if ((p & 0x3) == 3) {
		return power_mod(a, (p+1) >> 2 , p);
	}

	// initialization:
	// compute k, s : p = 2^s (2k+1) + 1
	k = p - 1;
	s = 0;

	while (!(k & 1)) {
		++s;
		k >>= 1;
	}

	k = (k - 1) >> 1;

	r = power_mod(a, k, p);
	n = multiply_mod(r, r, p);
	n = multiply_mod(n, a, p);
	r = multiply_mod(r, a, p);

	if (n == 1) {
		return r;
	}

	// look for quadratic nonresidue
	z = 2;
	while (jacobi(z, p) == 1) {
		++z;
	}

	c = power_mod(z, (k << 1) + 1, p);

	// (5) iteration
	while (n > 1) {
		k = n;
		t = s;
		s = 0;

		while (k != 1) {
			k = multiply_mod(k, k, p);
			s++;
		}

		t = t - s;

		c = power_mod(c, udigit(1 << (t-1)), p);
		r = multiply_mod(r, c, p);
		c = multiply_mod(c, c, p);
		n = multiply_mod(n, c, p);
	}

	return r;
}



// ***********************************************************************

udigit
next_prime (udigit x)
{
	long	zmod3, zmod5, zmod7;
	udigit	prim;
	static udigit table[] = { 3, 5, 5, 7, 7, 11 };
	// table[i] = next prime greater than (i+2)


	if (x <= 1) {
		return static_cast<udigit>(2);
	}

	if (x <= 7) {
		return table[x-2];
	}

	// prim = x; make prim odd
	prim = x + 1;
	if (!(prim & 1)) {
		prim++;
	}

	// initialize modular counters
	zmod3 = prim % 3;
	zmod5 = prim % 5;
	zmod7 = prim % 7;

	// while not a prime number
	while (!is_prime(prim, 10)) {
		do {
			// increase prim by 2
			prim += 2;

			zmod3 = (zmod3 + 2) % 3;
			zmod5 = (zmod5 + 2) % 5;
			zmod7 = (zmod7 + 2) % 7;
		}
		// until it is not divisible by small primes
		while (!(zmod3 && zmod5 && zmod7));
	}
	return prim;
}



// ***********************************************************************

udigit
previous_prime (udigit x)
{
	long	zmod3, zmod5, zmod7;
	udigit	prim;
	static udigit table[] = { 0, 2, 3, 3, 5, 5, 7, 7, 7, 7 };
	// table[i] = greatest prime smaller than (i+2)


	if (x <= 2) {
		return (static_cast<udigit>(0));
	}

	if (x <= 11) {
		return table[x-2];
	}

	// prim = x; make prim odd
	prim = x - 1;
	if (!(prim & 1)) {
		prim--;
	}

	// initialize modular counters
	zmod3 = prim % 3;
	zmod5 = prim % 5;
	zmod7 = prim % 7;

	// while not a prime number
	while (!is_prime(prim, 10)) {
		do {
			// decrease prim by 2
			prim -= 2;

			zmod3 = (zmod3 + 1) % 3;
			zmod5 = (zmod5 + 3) % 5;
			zmod7 = (zmod7 + 5) % 7;
		}
		// until it is not divisible by small primes
		while (!(zmod3 && zmod5 && zmod7));
	}
	return prim;
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
