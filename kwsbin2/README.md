## lattice-word-index-utterance

This tool is used to create an inverted index of the given lattices, where the
score of each word is the probability that the word occurs in any of the
transcriptions of the utterance at least once.

For additional details, check the paper "Probabilistic interpretation and
improvements to the HMM-filler for handwritten keyword spotting", by
J. Puigcerver et al.

For instance, let's consider the following lattice, whose arcs represent
complete words in some pre-defined vocabulary:

![Composition FST](egs/word_lat.png?raw=true)

We can obtain the inverted index by simply calling:

```bash
$ lattice-word-index-utterance ark:egs/lattice.ark.txt ark,t:-
```

The previous command will output the following to the standard output:

```
lat1 2 0 ; 5 0 ; 6 0 ; 7 0 ; 8 0 ; 3 -0.2231435 ; 1 -1.609438 ; 4 -1.609438
```

First the key of the lattice (i.e. `lat1`) is shown. Then a sequence of
tuples (each tuple separated with the caracter `;`) is shown.
The first element in each tuple represents the word label, and the following
number is the log-probability that the word appears somewhere in the utterance.
The sequence of tuples is sorted by decreasing log-probability.

This information is more easily read in the following table:

| Word   | Probability |
|--------|-------------|
| the    | 1.0         |
| is     | 1.0         |
| best   | 1.0         |
| friend | 1.0         |
| dog    | 0.8         |
| a      | 0.2         |
| lizard | 0.2         |
