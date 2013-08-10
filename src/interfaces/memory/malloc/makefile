######################################################################
#
#  Name ............ LiDIA/src/base/interfaces/memory/malloc/makefile
#
#  Author       : Thomas Papanikolaou, Thomas Pfahler
#  Copyright    : 1994-1999 by the LiDIA Group
#                 Technische Universitaet Darmstadt
# $Id: makefile,v 1.1.1.1 2001/05/02 16:28:03 lidiaadm Exp $


#
# include the configuration file
#

include ../../../../../config/CONFIG

#
#
#

OBJ=new.o

#
# Rules
#

#
# Uncomment $(libp)LiDIA.a to compile new.o: this files includes
# oveloadings of the global new & delete operators. At the moment
# this feature has been only tested with g++ and may not work with
# other compilers.
#
 
all:	copyh # $(libp)LiDIA.a


#
# Dependencies
# 

$(libp)LiDIA.a:	$(OBJ)

copyh:
	-$(CP) gmm.h $(LIDIA_INC_DIR)

#
# Clean
#

clean:          _clean_
veryclean:      _veryclean_
distclean:      veryclean

#
# End Of Makefile
#
