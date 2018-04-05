#!/bin/bash
set -e;
export LC_NUMERIC=C;

# Create output lattice representing the same language as in input.txt
# by removing the CTC blank output labels.
../../lattice-remove-ctc-blank 1 ark:input.txt ark,t:output.txt;

# Plot original lattice in PDF
lattice-to-fst --acoustic-scale=1 --lm-scale=1 ark:input.txt ark,t:- | \
    tail -n+2 | fstcompile | \
    fstdraw --portrait --isymbols=symbs.txt --osymbols=symbs.txt | \
    dot -Tpdf > input.pdf

# Plot output lattice in PDF
lattice-to-fst --acoustic-scale=1 --lm-scale=1 ark:output.txt ark,t:- | \
    tail -n+2 | fstcompile | \
    fstdraw --portrait --isymbols=symbs.txt --osymbols=symbs.txt | \
    dot -Tpdf > output.pdf

# Obtain the costs of all paths in the input lattice (there are 27 total
# paths). Note: in the input lattice all costs are written in the graph
# part of the cost.
lattice-to-nbest --n=27 ark:input.txt ark:- | \
    nbest-to-linear ark:- ark,t:- ark:/dev/null ark,t:- | \
    awk '{
      if (NF == 2) {
        COST[$1] = $2;
      } else {
        PATH[$1] = $2;
        for(i=3;i<=NF;++i) {
          PATH[$1] = PATH[$1]" "$i;
        }
      }
    }END{
      for (n in PATH) {
        print PATH[n], COST[n];
      }
    }' | sort -V > input_alignments.txt

# Obtain the costs of all paths in the output lattice (there are 27 total
# paths).
lattice-to-nbest --n=27 ark:output.txt ark:- | \
    nbest-to-linear ark:- ark,t:- ark:/dev/null ark,t:- | \
    awk '{
      if (NF == 2) {
        COST[$1] = $2;
      } else {
        PATH[$1] = $2;
        for(i=3;i<=NF;++i) {
          PATH[$1] = PATH[$1]" "$i;
        }
      }
    }END{
      for (n in PATH) {
        print PATH[n], COST[n];
      }
    }' | sort -V > output_alignments.txt

# Compare the costs of each of the alignments.
paste input_alignments.txt output_alignments.txt | awk '{
   d = ($4 - $8);
   r = sqrt(d * d) / sqrt($4 * $4);
   if (r > 1E-5) {
     print "Alignment costs do not match! "$0 > "/dev/stderr";
     error = 1;
   }
}END{
  if(error) { exit(1); }
}'

# Obtain the transcriptions from the input alignments.
sed -r 's|([0-9]+)( \1)+ |\1 |g' input_alignments.txt | \
    awk '{ for (i=1;i<NF;++i) { if ($i == 1) $i=""; } print $0; }' | \
    sed -r 's|[ ]+| |g;s|^[ ]+||g' | awk '{
      n = split($0, A, " ");
      if (n == 1) { COST["<eps> "] += A[n]; }
      else { $NF = ""; COST[$0] += A[n]; }
    }END{
      for(n in COST) print n, COST[n];
    }' | sort > input_transcripts.txt;

# Obtain the transcriptions directly from the output lattice.
lattice-to-nbest --n=27 ark:output.txt ark:- | \
    nbest-to-linear ark:- ark:/dev/null ark,t:- ark,t:- | \
    awk '{
      if (!($1 in PATH)) {
        PATH[$1] = $2;
        for(i=3;i<=NF;++i) {
          PATH[$1] = PATH[$1]" "$i;
        }
      } {
        COST[$1] = $2;
      }
    }END{
      for (n in PATH) {
        print PATH[n], COST[n];
      }
    }' | sed -r 's|[ ]+| |g;s|^[ ]+||g' | awk '{
      n = split($0, A, " ");
      if (n == 1) { COST["<eps> "] += A[n]; }
      else { $NF = ""; COST[$0] += A[n]; }
    }END{
      for(n in COST) print n, COST[n];
    }' | sort > output_transcripts.txt;

# Compare the costs of each of the transcriptions.
paste input_transcripts.txt output_transcripts.txt | awk '{
   a=$(NF / 2);
   b=$NF;
   d = (a - b);
   r = sqrt(d * d) / sqrt(a * a);
   if (r > 1E-5) {
     print "Transcription costs do not match! "$0 > "/dev/stderr";
     error = 1;
   }
}END{
  if(error) { exit(1); }
}'
