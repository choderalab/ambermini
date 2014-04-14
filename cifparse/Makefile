include ../config.h

.SUFFIXES:
.SUFFIXES:  .c .o

.c.o:
	$(CC) -c $(CNOOPTFLAGS) $(CFLAGS) $(AMBERCFLAGS) $<

CIFOBS=	\
	lex.cif.o \
	cifp.tab.o \
	cifparse.o

libcifparse.a: $(CIFOBS)
	/bin/rm -f libcifparse.a
	$(AR) libcifparse.a $(CIFOBS)
	-ranlib libcifparse.a

cifp.tab.h: cifp.tab.c

cifp.tab.c:  cifparse.y lex.cif.c
	$(YACC) -v -d cifparse.y
	sed 's/yy/cifp/g' < y.tab.c > cifp.tab.c
	sed 's/yy/cifp/g' < y.tab.h > cifp.tab.h
	/bin/rm -f y.tab.c y.tab.h

cifp.tab.o: cifp.tab.c
	$(CC) -c -D$(LEX) $(CNOOPTFLAGS) $(CFLAGS) $(AMBERCFLAGS) cifp.tab.c

lex.cif.c: cifparse.l 
	$(LEX) -t cifparse.l | sed 's/yy/cifp/g' > lex.cif.c

install:  libcifparse.a
	mv libcifparse.a $(LIBDIR)

clean:
	/bin/rm -f *.o libcifparse.a y.output
	/bin/rm -f lex.cif.c cifp.tab.c cifp.tab.h

uninstall: clean
	/bin/rm -f $(LIBDIR)/libcifparse.a

lex.cif.o:	\
	cifp.tab.h
