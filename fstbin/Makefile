all:

include ../kaldi.mk

BINFILES = fst-compose fst-compose-sum fst-info fst-normalize

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
install: $(BINFILES) fst-determinizestar fst-minimize fst-project
	test -d $(PREFIX) || mkdir -p $(PREFIX)/bin
	install -m 0755 $(BINFILES) \
	fst-determinizestar fst-minimize fst-project \
	$(PREFIX)/bin
