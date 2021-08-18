DIR = ./
PROGRAMS = $(DIR)main
CC = g++
DIR       = OBJECTS/
DIREX     = RUN/


CFLAGS=-O3
LINKER= $(CC) $(CFLAGS)
LIBS=  

INCLUDES=



OBJ = $(DIR)main.o\


MAKEFILE = makefile

$(DIR)%.o: %.cpp $(MAKEFILE)
	$(CC) -c $(CFLAGS) -o $@ $(INCLUDES) $*.cpp

$(DIREX)deto2D:  $(OBJ) $(MAKEFILE)
	$(LINKER) $(OBJ)  -o $@ $(LIBS)

clean:  
	rm -f $(DIR)*.o	
