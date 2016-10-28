include ../../../config.h

HEADERS = xpm.h

SRCS = data.c create.c misc.c rgb.c scan.c parse.c hashtable.c 	 XpmCrBufFrI.c XpmCrDataFrP.c XpmCrPFrBuf.c XpmRdFToI.c XpmWrFFrI.c 	 XpmCrBufFrP.c XpmCrIFrBuf.c XpmCrPFrData.c XpmRdFToP.c XpmWrFFrP.c 	 XpmCrDataFrI.c XpmCrIFrData.c XpmRdFToData.c XpmWrFFrData.c

OBJS = data.o create.o misc.o rgb.o scan.o parse.o hashtable.o 	 XpmCrBufFrI.o XpmCrDataFrP.o XpmCrPFrBuf.o XpmRdFToI.o XpmWrFFrI.o 	 XpmCrBufFrP.o XpmCrIFrBuf.o XpmCrPFrData.o XpmRdFToP.o XpmWrFFrP.o 	 XpmCrDataFrI.o XpmCrIFrData.o XpmRdFToData.o XpmWrFFrData.o

.c.o:
	$(CC) -c -I$(XHOME)/include $(CFLAGS) $(AMBERCFLAGS) -o $@ $<

libXpm.a:  $(OBJS)
	$(AR) $@ $(OBJS)
	$(RANLIB) $@

clean:
	/bin/rm -f *.o *.a
