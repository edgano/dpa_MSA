CC    = cc
CFLAGS=  -O
OBJECTS= main.o faces.o bias.o primer.o ecalc.o fixer.o
msa : ${OBJECTS}
	${CC} ${CFLAGS} ${OBJECTS} -o msa

main.o : main.c
	${CC} ${CFLAGS} -c main.c 

${OBJECTS} : defs.h

clean:
	rm -f *.o
