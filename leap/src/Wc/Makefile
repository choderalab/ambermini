include ../../../config.h

HEADERS = COPY.h \
	COPY_X.h \
	MapAg.h \
	MapAgP.h \
	PortableC.h \
	WcCreate.h \
	WcCreateP.h \
	done.h

SRCS = 	MapAg.c \
	WcActCB.c \
	WcConvert.c \
	WcCreate.c \
	WcInvoke.c \
	WcLateBind.c \
	WcLoadRes.c \
	WcName.c \
	WcReg.c \
	WcRegXt.c \
	WcSetValue.c \
	WcTemplate.c \
	WcWarn.c \
	WcNameWidget.c

OBJS = 	MapAg.o \
	WcActCB.o \
	WcConvert.o \
	WcCreate.o \
	WcInvoke.o \
	WcLateBind.o \
	WcLoadRes.o \
	WcName.o \
	WcReg.o \
	WcRegXt.o \
	WcSetValue.o \
	WcTemplate.o \
	WcWarn.o \
	WcNameWidget.o

.c.o:
	$(CC) -c -I$(XHOME)/include $(CFLAGS) $(AMBERCFLAGS) -o $@ $<

libWcLeap.a:  $(OBJS)
	$(AR) $@ $(OBJS)
	$(RANLIB) $@

clean:
	rm -f *.o *.a
