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
//	Author	: Thomas Papanikolaou (TP), Markus Maurer (MM)
//	Changes	: See CVS log
//
//==============================================================================================


#ifdef HAVE_CONFIG_H
# include	"config.h"
#endif
#include	"LiDIA/bigmod.h"



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



//
// class variables
//

residue_class< bigint > * bigmod::Mp = NULL;
bigint bigmod::M = bigint(0UL);



//
// class modifiers
//

void
bigmod::set_modulus (const bigint & m)
{
	if (m.is_zero ())
		lidia_error_handler ("bigmod", "set_modulus::zero modulus");

	if (base_bigmod::L == NULL)
		base_bigmod::L = new residue_class_list< bigint >;

	else if (bigmod::Mp != NULL)
		(base_bigmod::L)->clear (bigmod::Mp);

	if (m.is_positive())
		bigmod::M = m;
	else
		bigmod::M = -m;

	bigmod::Mp = (base_bigmod::L)->insert(bigmod::M);
}



//
// modifiers
//

void
bigmod::randomize()
{
	I.assign(LiDIA::randomize(bigmod::M));
	if (I.compare(bigmod::M) == 0)
		I.assign_zero();
}



//
// I/O
//

std::istream &
operator >> (std::istream & in, bigmod & a)
{
	in >> a.I;
	a.normalize();
	return in;
}



std::ostream &
operator << (std::ostream & out, const bigmod & a)
{
	out << a.I;
	return out;
}



int
string_to_bigmod (char *s, bigmod & a)
{
	int i = string_to_bigint(s, a.I);
	a.normalize();
	return i;
}



int
bigmod_to_string (const bigmod & a, char *s)
{
	return bigint_to_string(a.I, s);
}



#ifdef C_STDIO

//
// using fread/fwrite
//

void bigmod::read_from_file(FILE * fp)
{
	I.read_from_file(fp);
	this->normalize();
}



void bigmod::write_to_file (FILE * fp)
{
	I.write_to_file(fp);
}



//
// using fscanf/fprintf
//

void bigmod::scan_from_file (FILE * fp)
{
	I.scan_from_file(fp);
	this->normalize();
}



void bigmod::print_to_file(FILE * fp)
{
	I.print_to_file(fp);
}

#endif	// C_STDIO



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
