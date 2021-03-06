## lattice-word-index-utterance

This tool is used to create an inverted index of the given lattices, where the
score of each word is the probability that the word occurs in any of the
transcriptions of the utterance at least once.

For additional details, check the paper "Probabilistic interpretation and
improvements to the HMM-filler for handwritten keyword spotting", by
J. Puigcerver et al. (2015).

For instance, let's consider the following lattice, whose arcs represent
complete words in some pre-defined vocabulary:

![Word FST](egs/word_lat.png?raw=true)

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

For additional details, check the PhD thesis "A Probabilistic Formulation of
Keyword Spotting", by J. Puigcerver et al. (2018).

For instance, let's consider the previous lattice and suppose that associated
time with each state in the previous Lattice is given by the following table:

|           |     |     |     |     |     |     |     |     |     |     |
|-----------|-----|-----|-----|-----|-----|-----|-----|-----|-----|-----|
| __State__ |  0  |  1  |  2  |  3  |  4  |  5  |  6  |  7  |  8  |  9  |
| __Time__  |  0  |  2  |  9  |  4  |  8  | 12  | 16  | 22  | 27  | 33  |

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
(i.e. the segment where the word appeared), and the final number is the
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


## lattice-word-index-position

This tool can be used to create a positional inverted index of the given
lattices, in the traditional meaning of "position" in the context of search
engines. That is, the probability that a word appears at some position within
the transcription, for all possible transcriptions of the utterance.

For additional details, check the PhD thesis "A Probabilistic Formulation of
Keyword Spotting", by J. Puigcerver et al. (2018).

For instance, using our previous example. We can obtain the positional inverted
index with:

```bash
./lattice-word-index-position ark:egs/lattice.ark.txt ark,t:-
```

The previous command will output the following to the standard output:

```
lat1 2 4 12 16 0 ; 5 3 8 12 0 ; 6 5 16 22 0 ; 7 6 22 27 0 ; 8 7 27 33 0 ; 2 1 0 4 -0.2231435 ; 3 2 4 8 -0.2231435 ; 1 1 0 2 -1.609438 ; 4 2 2 9 -1.609438
```

First the key of the lattice (i.e. `lat1`) is shown. Then a sequence of
tuples (each tuple separated with the character `;`) is presented.
The first element in each tuple represents the word label, the second element
is the position of the word within a word sequence, the third element is the
log-probability that the word appears in such position (in any sequence),
finally the last two elements represent the segment in which the keyword
was found (initial and final frames).
The sequence of tuples is sorted by decreasing log-probability.

This information is more easily read in the following table:

| Word   | Position | Segment | Probability |
|--------|----------|---------|-------------|
| the    | 4        | 12--16  | 1.0         |
| is     | 3        | 8--12   | 1.0         |
| man's  | 5        | 16--22  | 1.0         |
| best   | 6        | 22--27  | 1.0         |
| friend | 7        | 27--33  | 1.0         |
| the    | 1        | 0--4    | 0.8         |
| dog    | 2        | 4--8    | 0.8         |
| a      | 1        | 0--2    | 0.2         |
| lizard | 2        | 2--9    | 0.2         |

Notice that there is only a single entry for the word "is" in the last table,
whose probability is the sum of the two entries from the table in the previous
section, and the segmentation shows the most likely segmentation.

## lattice-char-index-segment

This tool is used to build a segment-level word index from character lattices.

Characters (arcs in the lattice) are grouped into groups, and subpaths from the
lattice are extracted so that all arcs in the subpath are part of the
same group. This creates a straightforward approach to obtain "words" from
the character lattices, by grouping the whitespace characters and the rest of
characters in two different groups.

Characters (labels) that should not be indexed should be part of the
"wspace-group" (e.g. different types of whitespaces). Other groups that should
be considered as isolated words can also be specified using the
"--other-groups" option. Labels within a group must be separated by spaces and
groups are separated by a semicolon (e.g. --other-groups="1 2 3 ; 4 5 6 7"
creates two groups with 3 and 4 elements respectively). Labels not present in
any specific group, are grouped together in a default group.

Let's see an example for a character-level version of the previous lattice
(see the [egs/char_lat.pdf](egs/char_lat.pdf) file for a better resolution).

![Char FST](egs/char_lat.png?raw=true)

```bash
./lattice-char-index-segment "28" ark:egs/lattice.char.ark.txt ark,t:-
```

The previous command will output the following to the standard output:

```
lat1 13_1_14_27_19 16 21 0 ; 20_8_5 12 15 0 ; 2_5_19_20 22 26 0 ; 6_18_9_5_14_4 27 33 0 ; 9_19 9 11 -2.524243e-05 ; 20_8_5 1 4 -0.2231432 ; 4_15_7 5 8 -0.2231432 ; 1 0 1 -1.609439 ; 4_9_26_1_18_4 2 8 -1.609439
```

First the key of the lattice (i.e. `lat1`) is shown. Then a sequence of
tuples (each tuple separated with the character `;`) is presented.
The first element in each tuple is a string containing the sequence of
characters (labels) representing the word (the labels are separated by `_`).
The next two elements represent the initial and final frames of the segment.
Finally, the last element is the log-probability of such segment.
The sequence of tuples is sorted by decreasing log-probability.

This information is more easily read in the following table:

| Word        | Segment | Probability |
|-------------|---------|-------------|
| m a n ' s   | 16--21  | 1.0         |
| t h e       | 12--15  | 1.0         |
| b e s t     | 22--26  | 1.0         |
| f r i e n d | 27--33  | 1.0         |
| i s         | 9--11   | 1.0         |
| t h e       | 1--4    | 0.8         |
| d o g       | 5--8    | 0.8         |
| a           | 0--1    | 0.2         |
| l i z a r d | 2--8    | 0.2         |

Notice that the index is not equivalent to the
[lattice-word-index-segment](#lattice-word-index-segment) tool, since the
whitespaces are not considered part of the word in this case, and are excluded
from the corresponding word segmentation.

## lattice-char-index-position

This tool is used to build a position-level word index from character lattices.

The character lattice is transformed in order to ensure that each input arc
to each state is part of the same group of labels, in the same fashion as the
previous tool does. Then, the lattice is transformed again to ensure that all
paths arriving to them have traversed the same number of "groups". The latter
gives the word-position of the state within the transcription.

Using this transformation, we can easily obtain all the subpaths that traverse
states of the same group, for each of the sentence positions.

```bash
./lattice-char-index-position "28" ark:egs/lattice.char.ark.txt ark,t:-
```

The previous command will output the following to the standard output:

```
lat1 13_1_14_27_19 5 16 21 0 ; 20_8_5 4 12 15 0 ; 2_5_19_20 6 22 26 0 ; 6_18_9_5_14_4 7 27 33 0 ; 9_19 3 9 11 -2.524243e-05 ; 20_8_5 1 1 4 -0.2231432 ; 4_15_7 2 5 8 -0.2231445 ; 1 1 0 1 -1.609439 ; 4_9_26_1_18_4 2 2 8 -1.609497
```


First the key of the lattice (i.e. `lat1`) is shown. Then a sequence of
tuples (each tuple separated with the character `;`) is presented.
The first element in each tuple is a string containing the sequence of
characters (labels) representing the word (the labels are separated by `_`).
The second element is the word position within the transcription.
The next two elements represent the initial and final frames of the most likely
segment for that word and position.
And finally, the last element is the log-probability of such segment.
The sequence of tuples is sorted by decreasing log-probability.

This information is more easily read in the following table:

| Word        | Position | Segment | Probability |
|-------------|----------|---------|-------------|
| m a n ' s   | 5        | 16--21  | 1.0         |
| t h e       | 4        | 12--15  | 1.0         |
| b e s t     | 6        | 22--26  | 1.0         |
| f r i e n d | 7        | 27--33  | 1.0         |
| i s         | 3        | 9--11   | 1.0         |
| t h e       | 1        | 1--4    | 0.8         |
| d o g       | 2        | 5--8    | 0.8         |
| a           | 1        | 0--1    | 0.2         |
| l i z a r d | 2        | 2--8    | 0.2         |