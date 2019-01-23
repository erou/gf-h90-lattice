CC = gcc
LFLAGS = -lflint -lgmp -lm

%.o: %.c
	${CC} -c $<

test_AE.exe: AE.o test_AE.o
	${CC} AE.o test_AE.o -o test_AE.exe ${LFLAGS}

test_h90.exe: h90.o AE.o test_h90.o
	${CC} h90.o AE.o test_h90.o -o test_h90.exe ${LFLAGS}

test_minpoly.exe: minpoly.o test_minpoly.o
	${CC} minpoly.o test_minpoly.o -o test_minpoly.exe ${LFLAGS}

test: test_AE.exe test_h90.exe test_minpoly.exe
	./test_AE.exe
	./test_h90.exe
	./test_minpoly.exe

.PHONY: test
