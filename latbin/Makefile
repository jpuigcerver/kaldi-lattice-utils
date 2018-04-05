all:

include ../kaldi.mk

BINFILES = \
	lattice-best-path2 \
	lattice-info \
	lattice-remove-ctc-blank \
	lattice-restrict-length \
	lattice-to-transcript-length-dist \
	lattice-to-word-frame-post \
	lattice-to-word-position-post

OBJFILES =

TESTFILES =

ADDLIBS = \
	$(KALDI_ROOT)/src/lat/kaldi-lat.a \
	$(KALDI_ROOT)/src/lm/kaldi-lm.a \
	$(KALDI_ROOT)/src/fstext/kaldi-fstext.a \
	$(KALDI_ROOT)/src/hmm/kaldi-hmm.a \
	$(KALDI_ROOT)/src/tree/kaldi-tree.a \
	$(KALDI_ROOT)/src/util/kaldi-util.a \
	$(KALDI_ROOT)/src/matrix/kaldi-matrix.a \
	$(KALDI_ROOT)/src/base/kaldi-base.a

include $(KALDI_ROOT)/src/makefiles/default_rules.mk

PREFIX=/usr/local
install: $(BINFILES)
	test -d $(PREFIX) || mkdir -p $(PREFIX)/bin
	install -m 0755 $(BINFILES) $(PREFIX)/bin
