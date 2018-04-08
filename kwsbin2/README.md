## lattice-word-index-utterance

This tool is used to create an inverted index of the given lattices, where the
score of each word is the probability that the word occurs in any of the
transcriptions of the utterance at least once.

For additional details, check the paper "Probabilistic interpretation and
improvements to the HMM-filler for handwritten keyword spotting", by
J. Puigcerver et al. (2015).

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
tuples (each tuple separated with the character `;`) is presented.
The first element in each tuple represents the word label, and the following
number is the log-probability that the word appears somewhere in the utterance.
The sequence of tuples is sorted by decreasing log-probability.

This information is more easily read in the following table:

| Word   | Probability |
|--------|-------------|
| the    | 1.0         |
| is     | 1.0         |
| man's  | 1.0         |
| best   | 1.0         |
| friend | 1.0         |
| dog    | 0.8         |
| a      | 0.2         |
| lizard | 0.2         |


## lattice-word-index-segment

This tool is used to create a positional inverted index of the given lattices,
where the score of each word in a segment is the probability that the word
occurs in any of the transcriptions of the utterance at that specific time
segment.

For additional details, check the PhD thesis "A Probabilistic Formulation for
Keyword Spotting", by J. Puigcerver et al. (2018).

Suppose that associated time with each state in the previous Lattice is given
by the following table:

|       |     |     |     |     |     |     |     |     |     |     |
|-------|-----|-----|-----|-----|-----|-----|-----|-----|-----|-----|
| State |  0  |  1  |  2  |  3  |  4  |  5  |  6  |  7  |  8  |  9  |
| Time  |  0  |  2  |  9  |  4  |  8  | 12  | 16  | 22  | 27  | 33  |

We can obtain the time-positional inverted index with:

```bash
./lattice-word-index-segment ark:egs/lattice.ark.txt ark,t:-
```

The previous command will output the following to the standard output:

```
lat1 2 12 16 0 ; 6 16 22 0 ; 7 22 27 0 ; 8 27 33 0 ; 2 0 4 -0.2231435 ; 3 4 8 -0.2231435 ; 5 8 12 -0.2231435 ; 1 0 2 -1.609438 ; 4 2 9 -1.609438 ; 5 9 12 -1.609438
```

First the key of the lattice (i.e. `lat1`) is shown. Then a sequence of
tuples (each tuple separated with the character `;`) is presented.
The first element in each tuple represents the word label, and the next two
elements are the initial and final timesteps where the word is aligned
(i.e. the segment where the word appeared), and the finally number is the the
log-probability that the word appears in such segment.
The sequence of tuples is sorted by decreasing log-probability.

This information is more easily read in the following table:

| Word   | Segment | Probability |
|--------|---------|-------------|
| the    | 12--16  | 1.0         |
| man's  | 16--22  | 1.0         |
| best   | 22--27  | 1.0         |
| friend | 27--33  | 1.0         |
| the    | 0--4    | 0.8         |
| dog    | 4--8    | 0.8         |
| is     | 8--12   | 0.8         |
| a      | 0--2    | 0.2         |
| lizard | 2--9    | 0.2         |
| is     | 9--12   | 0.2         |
