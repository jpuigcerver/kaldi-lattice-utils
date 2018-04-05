## lattice-remove-ctc-blank

Remove CTC blank symbols from Kaldi lattices, so that the output of the
lattices are sequences of characters decoded as done by CTC.

### Input lattices

Input lattices are expected to store the frame-posteriors over all symbols of a
neural-net-based recognizer, including the blank/no-character symbol used for
CTC. Thus, input lattices **must be acceptors and be acyclic**.

### Output lattices

Output lattices encode the same paths (alignments) as the input lattices, but the
labels in the output arcs are modified so that the sequence of output symbols in
the output lattice represents the character-level transcription of the utterance,
according to the CTC decoding rules.

### How does it work?

By a simple composition of the input lattice and a FST with the following structure:

![Composition FST](egs/lattice-remove-ctc-blank/C.png?raw=true)

The previous FST transduces all sequences of contiguous character symbols (`a` and
`b` in the example) to a single character emision, and replaces all the
blank/no-character output symbols (`$` in the image) by epsilon.

Notice that the character output symbols are emitted at the first symbol of a
sequence of equal input characters.

For instance, given the input sequence `$ $ a a a $ a $ b b $ $ $`, the previous
transducer will produce the output `a a b`, with the following alignment:

|  1  |  2  |  3  |  4  |  5  |  6  |  7  |  8  |  9  |  10 |  11 |  12 |  13 |
|-----|-----|-----|-----|-----|-----|-----|-----|-----|-----|-----|-----|-----|
| `$` | `$` | `a` | `a` | `a` | `$` | `a` | `$` | `b` | `b` | `$` | `$` | `$` |
|     |     | `a` |     |     |     | `a` |     | `b` |     |     |     |     |
