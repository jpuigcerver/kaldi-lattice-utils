ifndef KALDI_ROOT
  $(error KALDI_ROOT environment variable is undefined)
endif

include $(KALDI_ROOT)/src/kaldi.mk
EXTRA_CXXFLAGS += -Wno-sign-compare -Wno-unused-variable -I$(KALDI_ROOT)/src
