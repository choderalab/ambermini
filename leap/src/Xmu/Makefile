include ../../../config.h

HEADERS = \
	Xmu.h \
	Atoms.h \
	CharSet.h \
	CloseHook.h \
	Converters.h \
	CurUtil.h \
	CvtCache.h \
	DisplayQue.h \
	Drawing.h \
	Editres.h \
	EditresP.h \
	Error.h \
	ExtAgent.h \
	Initer.h \
	Lookup.h \
	Misc.h \
	StdSel.h \
	StdCmap.h \
	SysUtil.h \
	WidgetNode.h \
	WinUtil.h \
	Xct.h


SRCS = \
	AllCmap.c \
	Atoms.c \
	ClientWin.c \
	Clip.c \
	CloseHook.c \
	CmapAlloc.c \
	CrCmap.c \
	CrPixFBit.c \
	CursorName.c \
	CvtCache.c \
	CvtStdSel.c \
	DefErrMsg.c \
	DelCmap.c \
	DisplayQue.c \
	Distinct.c \
	DrawLogo.c \
	DrRndRect.c \
	EditresCom.c \
	ExtAgent.c \
	FToCback.c \
	GetHost.c \
	GrayPixmap.c \
	Initer.c \
	LocBitmap.c \
	Lookup.c \
	LookupCmap.c \
	Lower.c \
	RdBitF.c \
	ScrOfWin.c \
	ShapeWidg.c \
	StdCmap.c \
	StrToBS.c \
	StrToBmap.c \
	StrToCurs.c \
	StrToGrav.c \
	StrToJust.c \
	StrToLong.c \
	StrToOrnt.c \
	StrToShap.c \
	StrToWidg.c \
	UpdMapHint.c \
	VisCmap.c \
	WidgetNode.c \
	Xct.c \
	sharedlib.c

OBJS = \
	AllCmap.o \
	Atoms.o \
	ClientWin.o \
	Clip.o \
	CloseHook.o \
	CmapAlloc.o \
	CrCmap.o \
	CrPixFBit.o \
	CursorName.o \
	CvtCache.o \
	CvtStdSel.o \
	DefErrMsg.o \
	DelCmap.o \
	DisplayQue.o \
	Distinct.o \
	DrawLogo.o \
	DrRndRect.o \
	EditresCom.o \
	ExtAgent.o \
	FToCback.o \
	GetHost.o \
	GrayPixmap.o \
	Initer.o \
	LocBitmap.o \
	Lookup.o \
	LookupCmap.o \
	Lower.o \
	RdBitF.o \
	ScrOfWin.o \
	ShapeWidg.o \
	StdCmap.o \
	StrToBS.o \
	StrToBmap.o \
	StrToCurs.o \
	StrToGrav.o \
	StrToJust.o \
	StrToLong.o \
	StrToOrnt.o \
	StrToShap.o \
	StrToWidg.o \
	UpdMapHint.o \
	VisCmap.o \
	WidgetNode.o \
	Xct.o

.c.o:
	$(CC) -c -I$(XHOME)/include $(CFLAGS) $(AMBERCFLAGS) -o $@ $<

libXmu.a:  $(OBJS)
	$(AR) $@ $(OBJS)
	$(RANLIB) $@

clean:
	/bin/rm -f *.o *.a
