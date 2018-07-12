all:

include ../kaldi.mk

TESTFILES = expand-subpaths-between-delimiters-test \
            expand-subpaths-labels-same-class-test \
            label-group-test \
            make-preceding-labels-same-class-test

ADDLIBS = \
	$(KALDI_ROOT)/src/fstext/kaldi-fstext.a \
	$(KALDI_ROOT)/src/util/kaldi-util.a \
	$(KALDI_ROOT)/src/matrix/kaldi-matrix.a \
	$(KALDI_ROOT)/src/base/kaldi-base.a

include $(KALDI_ROOT)/src/makefiles/default_rules.mk
