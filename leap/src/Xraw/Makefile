include ../../../config.h

HEADERS = 	3d.h 	AllWidgets.h 	Arrow.h 	ArrowP.h  	AsciiSink.h 	AsciiSinkP.h 	AsciiSrc.h 	AsciiSrcP.h 	AsciiText.h 	AsciiTextP.h 	Box.h 	BoxP.h 	Cardinals.h 	Clock.h 	ClockP.h 	Command.h 	CommandP.h 	color.h 	Container.h 	ContainerP.h 	Dialog.h 	DialogP.h 	Form.h 	FormP.h 	Frame.h 	FrameP.h 	Grip.h 	GripP.h 	Label.h 	LabelP.h 	List.h 	ListP.h 	Logo.h 	LogoP.h 	MenuButton.h 	MenuButtoP.h 	Mailbox.h 	MailboxP.h 	Paned.h 	PanedP.h 	Panner.h 	PannerP.h 	Porthole.h 	PortholeP.h 	Repeater.h 	RepeaterP.h 	Reports.h 	Scrollbar.h 	ScrollbarP.h 	ScrolledTable.h 	ScrolledTableP.h 	Separator.h 	SeparatorP.h 	Simple.h 	SimpleP.h 	SimpleMenu.h 	SimpleMenP.h 	Sme.h 	SmeP.h 	SmeBSB.h 	SmeBSBP.h 	SmeLine.h 	SmeLineP.h 	StripChart.h 	StripCharP.h 	table.h 	Table2.h 	TableP.h 	Text.h 	TextP.h 	TextSink.h 	TextSinkP.h 	TextSrc.h 	TextSrcP.h 	Toggle.h 	ToggleP.h 	Tree.h 	TreeP.h 	Viewport.h 	ViewportP.h 	XawAll.h 	XawInit.h

SRCS = 	3d.c 	AllWidgets.c 	Arrow.c 	AsciiSink.c 	AsciiSrc.c 	AsciiText.c 	Box.c 	Clock.c 	Command.c 	color.c 	Container.c 	Dialog.c 	Form.c 	Frame.c 	Grip.c 	Label.c 	List.c 	Logo.c 	MenuButton.c 	Mailbox.c 	Paned.c 	Panner.c 	Porthole.c 	Repeater.c 	Scrollbar.c 	ScrolledTable.c 	Separator.c 	Simple.c 	SimpleMenu.c 	Sme.c 	SmeBSB.c 	SmeLine.c 	StripChart.c 	table.c 	Table2.c 	Text.c 	TextSink.c 	TextSrc.c 	TextAction.c 	TextPop.c 	TextTr.c 	Toggle.c 	Tree.c 	Viewport.c 	Vendor.c 	XawInit.c

OBJS = 	3d.o 	AllWidgets.o 	Arrow.o 	AsciiSink.o 	AsciiSrc.o 	AsciiText.o 	Box.o 	Clock.o 	Command.o 	color.o 	Container.o 	Dialog.o 	Form.o 	Frame.o 	Grip.o 	Label.o 	List.o 	Logo.o 	MenuButton.o 	Paned.o 	Panner.o 	Porthole.o 	Repeater.o 	Scrollbar.o  	ScrolledTable.o 	Separator.o 	Simple.o 	SimpleMenu.o 	Sme.o 	SmeBSB.o 	SmeLine.o 	StripChart.o 	table.o 	Table2.o 	Text.o 	TextSink.o 	TextSrc.o 	TextAction.o 	TextPop.o 	TextTr.o 	Toggle.o 	Tree.o 	Vendor.o 	Viewport.o 	XawInit.o

.c.o:
	$(CC) -c -I$(XHOME)/include $(CFLAGS) $(AMBERCFLAGS) -o $@ $<

libXaw.a:  $(OBJS)
	$(AR) $@ $(OBJS)
	$(RANLIB) $@

clean:
	/bin/rm -f *.o *.a
