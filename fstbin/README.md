## fst-compose-sum

This tool computes the total sum (in the log-semiring) of the composition of
two FSTs.

By default, the composition is using the output language of the FSTs,
but the input labels can also be used, with `--use-input=true`.

Optionally, beam pruning can be used to reduce the size of the input FSTs so
that the composition is faster. For instance, with `--beam=20`.

The FSTs are normalized in the log-semiring before the composition, but this
can be disabled with `--normalize=false`.

For instance, let's consider an example FST:

```
example
0   1   1   1   0.35667494393873237891
0   1   2   2   1.20397280432593599262
1   2   1   1   0.22314355131420975576
1   2   2   2   -0.18232155679395462621
2
```

This FST has 4 path, with the following associated labels and costs:

| Path | Cost                         |
|------|------------------------------|
| 1 1  | -log(0.7 * 0.8) = 0.5798     |
| 1 2  | -log(0.7 * 1.2) = 0.1744     |
| 2 1  | -log(0.3 * 0.8) = 1.4271     |
| 2 2  | -log(0.3 * 1.2) = 1.0216     |

Let's compute the sum of the composition of this FST with itself.

```bash
$ fst-compose-sum ark:egs/example.fst.txt ark:egs/example.fst.txt
```

This produces the result:

```
example example 1.198654e+00
```

Notice that the composition of the normalized example FST with itself produces
a FST with the following paths and costs:

| Path | Cost                         |
|------|------------------------------|
| 1 1  | -log(0.7^2 * 0.4^2) = 2.5459 |
| 1 2  | -log(0.7^2 * 0.6^2) = 1.7350 |
| 2 1  | -log(0.3^2 * 0.4^2) = 4.2405 |
| 2 2  | -log(0.3^2 * 0.6^2) = 3.4296 |

Which the total sum (in the log-semiring) of the paths results in the previous
cost.


## fst-normalize

This tool "normalizes" the costs of a FST. Normalization usually means that
the sum of the costs of all paths is 0 (in the log-semiring).

For instance, let's consider an example FST:

```
example
0   1   1   1   0.35667494393873237891
0   1   2   2   1.20397280432593599262
1   2   1   1   0.22314355131420975576
1   2   2   2   -0.18232155679395462621
2
```

This FST has 4 path, with the following associated labels and costs:

| Path | Cost                         |
|------|------------------------------|
| 1 1  | -log(0.7 * 0.8) = 0.5798     |
| 1 2  | -log(0.7 * 1.2) = 0.1744     |
| 2 1  | -log(0.3 * 0.8) = 1.4271     |
| 2 2  | -log(0.3 * 1.2) = 1.0216     |

The sum of the costs of all paths (in the log-semiring) is roughly -0.6931
(i.e. -log(2.0)).


First, we normalize the FST so that the total sum (in the log-semiring) is 0.0:

```bash
$ fst-normalize ark:egs/example.fst.txt ark,t:-
```

This will produce the following normalized FST. Notice that weights are also
pushed to the initial state, but not labels.

```
example
0	1	1	1	0.356675
0	1	2	2	1.20397
1	2	1	1	0.916291
1	2	2	2	0.510826
2
```

The paths in this FST are:

| Path | Cost                       |
|------|----------------------------|
| 1 1  | -log(0.7 * 0.4) = 0.5798   |
| 1 2  | -log(0.7 * 0.6) = 0.1744   |
| 2 1  | -log(0.3 * 0.4) = 0.9163   |
| 2 2  | -log(0.3 * 0.6) = 0.5108   |


One can also perform a different kind of normalization, so that the path with
the smallest cost is 0 (i.e. normalization in the tropical-semiring).

In the previous example, notice that the path with the smallest cost
(i.e. `1 2`) has a cost of 0.1744 (i.e. -log(0.7 * 1.2)).

```bash
$ fst-normalize --use-log=false ark:egs/example.fst.txt ark,t:-
```

After the normalization in the tropical-semiring, the resulting FST is:

```
example
0	1	1	1
0	1	2	2	0.847298
1	2	1	1	0.405465
1	2	2	2
2
```

The paths in this FST are:

| Path | Cost                            |
|------|---------------------------------|
| 1 1  | -log(0.7 * 0.8 / 0.84) = 0.4055 |
| 1 2  | -log(0.7 * 1.2 / 0.84) = 0      |
| 2 1  | -log(0.3 * 0.8 / 0.84) = 1.2528 |
| 2 2  | -log(0.3 * 1.2 / 0.84) = 0.8473 |
