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
//	Author	: Nigel Smart (NS)
//	Changes	: See CVS log
//
//==============================================================================================


#ifdef HAVE_CONFIG_H
# include	"config.h"
#endif
#include	"LiDIA/curve_isomorphism.h"
#include	"LiDIA/curve_isomorphism.cc"
#include	"LiDIA/elliptic_curves/ec_arith.h"
#include	"LiDIA/bigint.h"
#include	"LiDIA/elliptic_curves/base_elliptic_curve.h"
#include	"LiDIA/elliptic_curves/base_elliptic_curve.cc"
#include	"LiDIA/point_bigint.h"



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



template class curve_isomorphism< bigrational, bigrational >;
template class curve_isomorphism< bigrational, bigint >;
template elliptic_curve< bigrational >
make_isomorphic_curve(const elliptic_curve< bigrational > &,
		      const bigrational&, const bigrational&,
		      const bigrational&, const bigrational&);

template elliptic_curve< bigrational >
make_isomorphic_curve(const elliptic_curve< bigint > &,
		      const bigrational&, const bigrational&,
		      const bigrational&, const bigrational&);


template void convert(bigrational &, const bigrational&);


//template std::ostream & operator << (std::ostream & , const
//				curve_isomorphism<bigrational,bigrational>&);

//template std::ostream & operator << (std::ostream & , const
//				curve_isomorphism<bigrational,bigint>&);
//template std::ostream & operator << (std::ostream & , const
//				curve_isomorphism<bigrational,bigrational>&);
//template std::ostream & operator << (std::ostream & , const
//				curve_isomorphism<bigrational,bigint>&);
template bool find_urst(const elliptic_curve< bigrational > &,
			const elliptic_curve< bigint > &, bigrational&,
			bigrational&, bigrational&, bigrational&);
template bool find_urst(const elliptic_curve< bigrational > &,
			const elliptic_curve< bigrational > &, bigrational&,
			bigrational&, bigrational&, bigrational&);




#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
