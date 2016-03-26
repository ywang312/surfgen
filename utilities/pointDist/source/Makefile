# Makefile for pointDist.x
OBJS	:= pointDist.o
CC	:= icc
COPT	:= -Ddebugging
# Executable
EXE	:= pointDist.x

%.o: %.c
	sh -c "$(CC) -c -g $(COPT) $< "

$(EXE):	$(OBJS)
	sh -c "$(CC) -o $@ $< $(BLAS_LIB)"

install: 
	@echo "Installing..."
	sh -c "mv pointDist.x ../../../bin/ "

clean:
	@echo "Cleaning..."
	rm -f ./*.o
