# ===========================================================================
#  Makefile iS                                     Chun Shen Apr. 9, 2013
# ===========================================================================
##
##  Environments :	MAIN	= 	main sourcefile	
##
##  Usage : 	(g)make	[all]		compile the whole project		
##			install	make all and copy binary to $INSTPATH
##			clean		remove objectfiles in obj_$TYPE 
##			distclean	remove all objectsfiles and binaries
##  

CC := icpc
CFLAGS= -O3 -Wall -lz -std=c++11

RM		=	rm -f
O               =       .o
LDFLAGS         =       $(CFLAGS)
SYSTEMFILES     =       $(SRCGNU)

# --------------- Files involved ------------------

ifeq "$(MAIN)" ""
MAIN		=	hadronic_afterburner_tools.e
endif

SRC		=	main.cpp arsenal.cpp \
            ParameterReader.cpp particleSamples.cpp HBT_correlation.cpp \
            single_particleSpectra.cpp particle_yield_distribution.cpp \
			particle_decay.cpp BalanceFunction.cpp

INC		= 	parameters.h arsenal.h ParameterReader.h \
            Stopwatch.h particleSamples.h \
            HBT_correlation.h single_particleSpectra.h \
			particle_yield_distribution.h particle_decay.h \
			particle_info.h BalanceFunction.h

# -------------------------------------------------

OBJDIR		=	obj
SRCFILES 	= 	$(SRC) $(INC) GNUmakefile
OBJECTS		=	$(addprefix $(OBJDIR)/, $(addsuffix $O, \
			$(basename $(SRC))))
TARGET		=	$(MAIN)
INSTPATH	=	../

# --------------- Pattern rules -------------------

$(OBJDIR)/%.o: %.cpp
	$(CC) $(CFLAGS) -c $< -o $@

%.cpp:
	if [ -f $@ ] ; then touch $@ ; else false ; fi

# -------------------------------------------------

.PHONY:		all mkobjdir clean distclean install

all:		mkobjdir $(TARGET)

help:
		@grep '^##' GNUmakefile

mkobjdir:	
		-@mkdir -p $(OBJDIR)

$(TARGET):	$(OBJECTS)	
		$(CC) $(OBJECTS) -o $(TARGET) $(LDFLAGS) 
#		strip $(TARGET)

clean:		
		-rm $(OBJECTS)

distclean:	
		-rm $(TARGET)
		-rm -r obj

install:	$(TARGET)
		cp $(TARGET) $(INSTPATH)

# --------------- Dependencies -------------------
./main.cpp: arsenal.h parameters.h Stopwatch.h ParameterReader.h particleSamples.h HBT_correlation.h single_particleSpectra.h particle_yield_distribution.h particle_decay.h
./arsenal.cpp: arsenal.h
./ParameterReader.cpp: ParameterReader.h arsenal.h
./particleSamples.cpp: ParameterReader.h particleSamples.h 
./HBT_correlation.cpp: HBT_correlation.h ParameterReader.h particleSamples.h
./single_particleSpectra.cpp: single_particleSpectra.h ParameterReader.h particleSamples.h
./particle_yield_distribution.cpp: particle_yield_distribution.h ParameterReader.h particleSamples.h
./particle_decay.cpp: particle_decay.h particleSamples.h

