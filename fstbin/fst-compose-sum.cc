// MIT License
//
// Copyright (c) 2018 Joan Puigcerver <joapuipe@gmail.com>
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

#include "base/kaldi-common.h"
#include "util/common-utils.h"
#include "fstext/fstext-lib.h"
#include "fstext/normalize_fst.h"


void PrepareFst(
    fst::VectorFst<fst::StdArc> *inp,  // Input FST (note: can be modified)
    fst::VectorFst<fst::LogArc> *out,  // Output FST
    float beam,  // Remove arcs from paths with cost > beam * best_path_cost
    bool normalize,  // If true, subtract total cost from the FST.
    bool project_input,
    bool ilabel_sort) {  // If true, sort arcs according to ilabels.
  if (project_input) {
    fst::Project(inp, fst::PROJECT_INPUT);
  } else {
    fst::Project(inp, fst::PROJECT_OUTPUT);
  }

  // If false, sort them according to olabels.
  // Prune FST using the given beam
  if (beam >= 0.0 && beam < std::numeric_limits<float>::infinity()) {
    fst::Prune(inp, fst::StdArc::Weight(beam));
  }

  // Convert StdArc to LogArc
  fst::ArcMap(*inp, out, fst::WeightConvertMapper<fst::StdArc, fst::LogArc>());

  // Subtract total cost from the FST.
  if (normalize) {
    NormalizeFst(out);
  }

  // Sort arcs
  if (ilabel_sort) {
    fst::ArcSort(out, fst::ILabelCompare<fst::LogArc>());
  } else {
    fst::ArcSort(out, fst::OLabelCompare<fst::LogArc>());
  }

  fst::Connect(out);
}

float ComputeTotalCostCompose(
    const fst::Fst<fst::LogArc> &fst1, const fst::Fst<fst::LogArc> &fst2) {
  // Compose FSTs.
  fst::VectorFst<fst::LogArc> comp_fst;
  fst::Compose(fst1, fst2, &comp_fst);
  // Compute total backward cost of each state.
  std::vector<fst::LogArc::Weight> costs;
  fst::ShortestDistance(comp_fst, &costs, true);
  // Return total sum (backward of initial state).
  return costs[comp_fst.Start()].Value();
}

int main(int argc, char *argv[]) {
  try {
    using namespace kaldi;
    using fst::VectorFst;
    using fst::StdArc;
    using fst::LogArc;

    float beam = std::numeric_limits<float>::infinity();
    bool use_inputs = false;
    bool normalize = true;

    const char *usage =
        "Compute the total cost of the composition of two FSTs.\n"
        "\n"
        "This can be used to compute: \\sum_{w} p(w | x_1) p(w | x_2). If the "
        "input FSTs represent p(w | x_1) and p(w | x_2), the result of the "
        "program will the the negative logarithm of this sum.\n"
        "\n"
        "Usage: fst-compose-sum [options] fsts1-rspecifier fsts2-rspecifier\n"
        "  e.g: fst-compose-sum --use-inputs=true ark:fst1.ark ark:fst2.ark\n";

    ParseOptions po(usage);
    po.Register("use-inputs", &use_inputs,
                "If true, use input symbols to compose the FSTs");
    po.Register("normalize", &normalize,
                "If true, normalize costs in each FST before the composition");
    po.Register("beam", &beam,
                "Prune the FSTs using this beam.");

    po.Read(argc, argv);
    if (po.NumArgs() != 2) {
      po.PrintUsage();
      exit(1);
    }

    const std::string fsts_rspecifier1 = po.GetArg(1);
    const std::string fsts_rspecifier2 = po.GetArg(1);

    SequentialTableReader<fst::VectorFstHolder> fst_reader1(fsts_rspecifier1);
    for (; !fst_reader1.Done(); fst_reader1.Next()) {
      const string &key1 = fst_reader1.Key();

      VectorFst<LogArc> fst1;
      PrepareFst(&fst_reader1.Value(), &fst1, beam, normalize, use_inputs,
                 false /* Sort w.r.t. olabel for fst::Compose */);
      fst_reader1.FreeCurrent();

      SequentialTableReader<fst::VectorFstHolder> fst_reader2(
          fsts_rspecifier2);
      for (; !fst_reader2.Done(); fst_reader2.Next()) {
        const string &key2 = fst_reader2.Key();

        VectorFst<LogArc> fst2;
        PrepareFst(&fst_reader2.Value(), &fst2, beam, normalize, use_inputs,
                   true /* Sort w.r.t. ilabel for fst::Compose */);
        fst_reader2.FreeCurrent();

        const float total_cost = ComputeTotalCostCompose(fst1, fst2);
        std::cout << key1 << " " << key2 << " "
                  << std::scientific << total_cost << std::endl;
      }
    }

    return 0;
  } catch (const std::exception &e) {
    std::cerr << e.what();
    return -1;
  }
}
