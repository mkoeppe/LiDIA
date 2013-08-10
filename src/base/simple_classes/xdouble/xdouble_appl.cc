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
//	Author	: Keith Briggs (KB)
//	Changes	: See CVS log
//
//==============================================================================================


#include	"LiDIA/xdouble.h"



#ifdef LIDIA_NAMESPACE
using namespace LiDIA;
#endif



int main_LiDIA(int argc, char** argv)
{
	std::cout.precision(17);
	xdouble one = 1.0, two = 2.0, x, y, z, e;
	int i;
	z = "1234567890.123456789012345678901234567890";
	std::cout << "Test 0: 1+1.0e-40 = ";
	std::cout << xdouble(1.0) + xdouble(1.0e-40) << std::endl;
	std::cout << "Test 1: input/output conversion\n";
	std::cout << "This should be 0.333333333: " << xdouble("0.333333333") << std::endl;
	std::cout << "This should be 0.0000333333333: " << xdouble("00000.0000333333333") << std::endl;
	std::cout << "This should be 333333333: " << xdouble("333333333") << std::endl;
	std::cout << "This should be 333333333e12: " << xdouble("333333333e12") << std::endl;
	std::cout << "This should be -333333333.3333: " << xdouble("-333333333.3333") << std::endl;
	std::cout << "This should be -333333333.33e-04: " << xdouble("-333333333.33e-04") << std::endl;
	std::cout << "The next output should be..\n";
	std::cout << "   1234567890.123456789012345678901234567890\n";
	std::cout << "z = " << z << "\n";
	std::cout << "Test 2: fmod(\", 10)\n";
	std::cout << "z mod 10: " << fmod(z, 10) << std::endl;
	z = floor(z);
	std::cout << "Test 3: floor(\"\")\n";
	std::cout << "Floor: " << (z) << std::endl;
	std::cout << "Test 4: exp(1)\n";
	std::cout << "Should be\n" <<
		" 2.7182818284590452353602874713526624977572470936999\n";
	std::cout.precision(32);
	z = exp(one);
	std::cout << z << "\n";
	e = exp(one);
	std::cout << "Test 5: log\n";
	z = log(recip(e * e));
	std::cout << "log(1/e/e): " << z << std::endl;
	z = log(recip(e));
	std::cout << "log(1/e): " << z << std::endl;
	z = log(e);
	std::cout << "log(e): " << z << std::endl;
	y = log(e * e);
	std::cout << "log(e^2): " << y << std::endl;
	y = log(e * e * e);
	std::cout << "log(e^3): " << y << std::endl;
	z = 0.159;
	std::cout << "log(0.159): " << log(z) << std::endl;
	std::cout << "Test 6: exp: next output should be 1.23456789\n";
	// exp(z) is exactly 1.23456789...
	z = "0.2107210222156525610500017104882905489049";
	z = exp(z);
	std::cout << z << "\n";
	std::cout << "Test 7: two reciprocals of previous result\n";
	z = one / z;
	std::cout << z << "\n";
	z = one / z;
	std::cout << z << "\n";
	std::cout << "Test 8: multiplication: next output should be 2\n";
	z = "1.4142135623730950488016887242096980785696718753769";
	std::cout << z * z << std::endl;
	std::cout << "Timing...\n";
	for (i = 1; i < 100000; i++)
		y = z * z;
	std::cout << " sqrt(2) should be " <<
		" 1.4142135623730950488016887242096980785696718753769\n";
	std::cout << "Test 9: Start with 2, sqrt five times and sqr five times\n";
	x = sqrt(two);
	std::cout << x << "\n";
	x = sqrt(x);
	std::cout << x << "\n";
	x = sqrt(x);
	std::cout << x << "\n";
	x = sqrt(x);
	std::cout << x << "\n";
	x = sqrt(x);
	std::cout << x << "\n";
	x = x * x;
	std::cout << x << "\n";
	x = x * x;
	std::cout << x << "\n";
	x = x * x;
	std::cout << x << "\n";
	x = x * x;
	std::cout << x << "\n";
	x = x * x;
	std::cout << x << "\n";
	std::cout << "Test 10: Comparison\n";
	x = "123456789.1234567890123456789";
	y = "123456789.1234567891";
	std::cout << "Next 3 should be true...\n";
	std::cout << (x < y) << std::endl;
	std::cout << (y >= x) << std::endl;
	std::cout << (y != x) << std::endl;
	std::cout << "Next should be false...\n";
	std::cout << (y == x) << std::endl;
	std::cout << (x > y) << std::endl;
	std::cout << "Test 11: sin\n";
	const xdouble exact_sin[] = {
		"0.000099999999833333333416666666646825396828152557318973",
		"0.0009999998333333416666664682539710097001513147348086",
		"0.0099998333341666646825424382690997290389643853601692",
		"0.099833416646828152306814198410622026989915388017982",
		"0.84147098480789650665250232163029899962256306079837",
		"-0.54402111088936981340474766185137728168364301291622",
		"-0.50636564110975879365655761045978543206503272129066",
		"0.82687954053200256025588742910921814121272496784779",
		"-0.30561438888825214136091003523250697423185004386181",
	};
	x = "0.0001";
	for (i = 0; i < 9; i++) {
		y = sin(x);
		std::cout << "sin(" << x << ") = " << y << "; abs error = ";
		std::cout << y - exact_sin[i] << std::endl;
		x = 10.0 * x;
	}
	std::cout << "Test 12: arctan\n";
	const xdouble exact_atan[] = {	// exact_atan[i] = atan(10^(i-4))
		"0.000099999999666666668666666652380952492063491154401162",
		"0.00099999966666686666652380963492054401162093455426801",
		"0.0099996666866652382063401162092795485613693525443766",
		"0.099668652491162027378446119878020590243278322504315",
		"0.78539816339744830961566084581987572104929234984378",
		"1.4711276743037345918528755717617308518553063771832",
		"1.5607966601082313810249815754304718935372153471432",
		"1.5697963271282297525647978820048308980869637651333",
		"1.5706963267952299525626550249873704896065212085332"
	};
	x = "0.0001";
	for (i = 0; i < 8; i++) {
		xdouble S = atan(x);
		std::cout << "arctan(" << x << ") = " << S << "; abs error = ";
		std::cout << S - exact_atan[i] << std::endl;
		x = 10.0 * x;
	}
	x = "0.1";
	std::cout << "Test 13: tan(arctan(0.1)) = " << sin(atan(x)) / cos(atan(x)) << std::endl;
	std::cout << "Test 14: input conversion\n";
	std::cout << "Please enter a number, terminate with< Enter > : ";
	std::cin >> x;
	std::cout << "\nYou entered " << x << std::endl;

	return 0;
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
