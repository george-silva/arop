#############################################################################
# !make
#
# makefile name: makefile (for linux)
#
##!END
#############################################################################
CC = cc -m32

TARGET = ortho

INC = -I. -I$(MYINC)

LIB = -L$(MYLIB) -lgeotiff  -ltiff $(MYLIB)/geolib.a -lm


OBJ = lndortho_main.o lndortho_compute.o lndortho_io.o lndortho_util.o lndortho_gchk.o
INC_FILES = lndortho.h

all: $(TARGET)

$(OBJ) : $(INC_FILES)

# Make the process
$(TARGET) : $(OBJ)
	$(CC) $(OBJ) $(LIB)  -o $(TARGET)

clean:
	rm -f *.o $(TARGET)

install:
	cp $(TARGET) $(BIN)

#
# Rules
#
.c.o: $(INC_FILES)
	$(CC) $(INC) -c $< -o $@

#******************* End of make file *******************************


