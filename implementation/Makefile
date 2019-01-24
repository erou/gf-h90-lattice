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

test_tensor.exe: tensor.o test_tensor.o
	${CC} tensor.o test_tensor.o -o test_tensor.exe ${LFLAGS}

test_tensor.exe: tensor.o test_tensor_scalarmul.o
	${CC} tensor.o test_tensor_scalarmul.o -o test_tensor_scalarmul.exe ${LFLAGS}

test: test_AE.exe test_h90.exe test_minpoly.exe test_tensor.exe test_tensor_scalarmul.exe
	./test_AE.exe
	./test_h90.exe
	./test_minpoly.exe
	./test_tensor.exe
	./test_tensor_scalarmul.exe

.PHONY: test
