CC = gcc
LFLAGS =  -lflint -lgmp -lm

%.o: %.c
	${CC} -c $<

test_AE.exe: AE.o test_AE.o
	${CC} ${LFLAGS} AE.o test_AE.o -o test_AE.exe

test_h90.exe: h90.o AE.o test_h90.o
	${CC} ${LFLAGS} h90.o AE.o test_h90.o -o test_h90.exe

test_minpoly.exe: minpoly.o test_minpoly.o
	${CC} ${LFLAGS} minpoly.o test_minpoly.o -o test_minpoly.exe

test: test_AE.exe test_h90.exe test_minpoly.exe
	./test_AE.exe
	./test_h90.exe
	./test_minpoly.exe

.PHONY: test
