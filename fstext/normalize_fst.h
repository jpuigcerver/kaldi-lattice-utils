//
// Created by joapuipe on 28/03/18.
//

#ifndef KALDI_LATTICE_UTILS_NORMALIZE_FST_H
#define KALDI_LATTICE_UTILS_NORMALIZE_FST_H

#include <fst/fstlib.h>

template <typename Arc>
void NormalizeFst(fst::MutableFst<Arc>* fst) {
  typedef typename Arc::Weight Weight;
  fst::VectorFst<Arc> tmp(*fst);

  std::vector<Weight> costs;
  fst::ShortestDistance(tmp, &costs, true);

  if (costs[tmp.Start()] == Weight::Zero()) {
    fst->DeleteStates();
    return;
  }

  for (fst::StateIterator< fst::Fst<Arc> > si(tmp); !si.Done(); si.Next()) {
    const Weight final = tmp.Final(si.Value());
    if (final != Weight::Zero()) {
      tmp.SetFinal(si.Value(), fst::Divide(final, costs[tmp.Start()]));
    }
  }

  fst::Push<Arc, fst::REWEIGHT_TO_INITIAL>(tmp, fst, fst::kPushWeights);
}

#endif //KALDI_LATTICE_UTILS_NORMALIZE_FST_H
