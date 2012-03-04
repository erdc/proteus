#
# Makefile for mesh development/testing.
#	Project: meshTest
#
###########################################
#
CC = gcc
CFLAGS = 
CURDIR = ~/

INCPATH =
#INCPATH = -I$(CURDIR)include -I/$(CURDIR)tests

#TESTPATH = -I/$(CURDIR)tests

LINK = gcc

LFLAGS =

LIB = 

EXE = mesh_013107 

HEADERS = mesh.h
SRC = meshTests.cpp
OBJ = meshTests.o

.SUFFIXES: .cpp

.cpp.o: 
	$(CC) -c $(CFLAGS) $(INCPATH) -o $@ $<

all : $(EXE)

$(EXE): $(OBJ)
	$(LINK) $(LFLAGS) -o $(EXE) $(OBJ) $(LIBS)

$(OBJ): $(SRC)

clean:
	-rm -f $(OBJ) $(EXE)
	-rm -f core
