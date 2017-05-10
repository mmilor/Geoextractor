#-I $MGDODIR/Root -I $MGDODIR/Base -I $CLHEP_INCLUDE_DIR -I $GELATIODIR/Utilities -I $MGDODIR/Gerda -I $MGDODIR/GerdaTransforms -I $MGDODIR
#-L $MGDODIR/GerdaTransforms -L $MGDODIR/Majorana -L $MGDODIR/Root -L $MGDODIR/Transforms -L $MGDODIR/Base -L $MGDODIR/Gerda -L $MGDODIR/lib
#$(MGDODIR)/Root/MGTEvent.o

EXE	= Geoextractor

CC	=	g++

COPTS	=	-fPIC -DLINUX -Wall \
                $(shell root-config --cflags) 

FLAGS	=	-Wall

DEPLIBS	=   	-l m 

LIBS	=	$(shell root-config --libs) -lGeom

INCLUDEDIR =	-I$(ROOTSYS)/include

OBJS	=	$(EXE).o

INCLUDES =	

#########################################################################

.PHONY all	:	$(EXE)

.PHONY clean	:	
	rm -f $(OBJS) $(EXE)

$(EXE)	:	$(OBJS) 
		/bin/rm -f $(EXE)
		$(CC) $(FLAGS) -o $(EXE) $(OBJS) $(LIBS) $(DEPLIBS)

$(OBJS)	:	$(INCLUDES) 

%.o	:	%.cxx
		$(CC) $(COPTS) $(INCLUDEDIR) -c -o $@ $<
