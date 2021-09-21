#$Id: Makefile 525 2018-03-19 21:54:26Z mexas $

# This Makefile will build the main CASUP documentation pdf
# and also invoke a ROBOdoc processing script, target "robodoc".

DOC_SRC=	doc.1
DOC_PARTS=	intro.1 avail.1 build.1 ca.1 \
		cgca_m2lnklst.1 cgca_addhead.1 \
		cgca_addmiddle.1 cgca_halo.1 cgca_inithead.1 \
		cgca_lstdmp.1 cgca_rmhead.1 \
		cgca_rmmiddle.1 layout.1 cgca_payload.1 cgca_node.1 \
		cgca_gcupd.1 cgca_m3clvg.1 \
		messages.1 modules.1 mod2.1 parafem.1 minfo.1 mwarn.1 \
		merr.1 mdebug.1
DOC_FULL=	$(DOC_SRC:.1=.full)
DOC_REF=	$(DOC_SRC:.1=.ref)
DOC_PS=		$(DOC_SRC:.1=.ps)
DOC_PSSEL=	$(DOC_SRC:.1=.pssel)
DOC_PDF=	$(DOC_SRC:.1=.pdf)
CLEAN+=		$(DOC_FULL) $(DOC_REF) $(DOC_PS) $(DOC_PSSEL) $(DOC_PDF)

GV_SRC=		lnklst1.gv lnklst2.gv lnklst3.gv lnklst4.gv lnklst5.gv \
		lnklst6.gv lnklst7.gv
GV_FIG=		$(GV_SRC:.gv=.fig)
FIG_SRC=	$(GV_FIG) gen1.fig gen2.fig gen3.fig gen4.fig \
		gen5.fig gen6.fig centr1.fig centr2.fig cgca0.fig \
		cgca1.fig cgca2.fig cgca3.fig
FIG_PIC=	$(FIG_SRC:.fig=.pic)
CLEAN+=		$(GV_FIG) $(FIG_PIC)

ALL_FIGS=	$(FIG_PIC)

CLEAN+=		../robodoc

IND=		Ind
 		
GROFF_DIR=	/usr/local/bin
GROFF=		$(GROFF_DIR)/groff
PDFROFF=	$(GROFF_DIR)/pdfroff
REFER=		$(GROFF_DIR)/refer
SOELIM=		$(GROFF_DIR)/soelim
FIG2DEV=	fig2dev
PS2PDF=		ps2pdf
DOT=		dot
AWK=		awk

GV_CORR_AWK=	gv-bug.awk
 
PSSELECT=	psselect
PSSELECT_FLAGS= -p1,_2,_1,2-_3

FC=		gfortran7

PGPLOT_LIB=	-lpgplot
LIBS=		-lX11 -lm

FFLAGS+=	-fno-backslash
LDFLAGS+=	-Wl,-rpath=/usr/local/lib/gcc7

SRCF=		cleavage1.f90
EXEF=		$(SRCF:.f90=.x)

CLV_PROB_PS=	cleavage1.ps
CLV_PROB_PDF=	$(CLV_PROB_PS:.ps=.pdf)
CLV_PROB_CLEAN=	$(CLV_PROB_PS) $(CLV_PROB_PDF)

CLV_SRC=	cleavage-fc.gv
CLV_IMG_PS=	$(CLV_SRC:.gv=.ps)
CLV_IMG_PDF=	$(CLV_SRC:.gv=.pdf)
CLV_IMG_CLEAN=	$(CLV_IMG_PS) $(CLV_IMG_PDF)

JPG_SRC=	Rhombicuboctahedron.jpg
JPG_PNG=	$(JPG_SRC:.jpg=.png)
JPG_PNG_CLEAN=	$(JPG_PNG)

PNG_FIGS=	Elongated_square_gyrobicupola.png $(JPG_PNG)

.SUFFIXES:	.f90 .x .fig .gv .ps .pdf .jpg .png .pic

all:		robodoc $(GV_FIG) $(DOC_PDF)

# need to keep $(GV_FIG) in the list of targets to make.
# Otherwise linux make thinks this is an *intermediate*
# target which it is free to remove quietly.

robodoc:
	./build.sh

.gv.fig:
	$(DOT) -Tfig -o tmpfile $<
	$(AWK) -f $(GV_CORR_AWK) tmpfile > $@
	rm tmpfile

.fig.pic:
	$(FIG2DEV) -Lpic $< $@

$(DOC_FULL): $(DOC_SRC) $(DOC_PARTS) $(ALL_FIGS) $(EXT_TMAC_FILE)
	$(SOELIM) $(DOC_SRC) > $@

$(DOC_REF): $(DOC_FULL) $(IND)
	$(REFER) -p $(IND) -eSP $(DOC_FULL) > $@

$(DOC_PS): $(DOC_REF) $(DOC_COVER)
	$(GROFF) -ms -ept $(DOC_REF) > $@

$(DOC_PSSEL): $(DOC_PS)
	$(PSSELECT) $(PSSELECT_FLAGS) $> $@

$(DOC_PDF): $(DOC_PSSEL)
	$(PS2PDF) $> $@

.f90.x:
	${FC} ${FFLAGS} ${LDFLAGS} -o $@ $< ${PGPLOT_LIB} ${LIBS}

$(CLV_PROB_PS): $(EXEF)
	./$(EXEF)

$(CLV_PROB_PDF): $(CLV_PROB_PS)
	$(PS2PDF) $< $@

clean:
	rm -rf $(CLEAN)
