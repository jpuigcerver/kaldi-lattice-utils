#!/bin/bash
set -e;
export LC_NUMERIC=C;

# Directory where the prepare.sh script is placed.
SDIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)";
cd "$SDIR" || exit 1;

function get_lattice_nbest () {
  lattice-to-nbest --n="$1" ark:- ark:- 2> /dev/null |
  nbest-to-linear ark:- ark:/dev/null ark,t:- ark,t:- 2> /dev/null |
  awk '{
    if (NR % 2 == 0) {
      COST[$1] = $2;
    } else {
      PATH[$1] = $2;
      for(i=3;i<=NF;++i) {
        PATH[$1] = PATH[$1]" "$i;
      }
    }
  }END{
    for (n in PATH) {
      print COST[n], PATH[n];
    }
  }' | sort -g
}

function nbest_word_to_char () {
  awk -v MF="$1" '
BEGIN{
  while((getline < MF) > 0) { M[$2] = $1; }
}{
  for (i=2; i<=NF; ++i)  {
    if ($i in M) {
      $i = M[$i];
    } else {
      print "Symbol "$i" was not found in the table!" > "/dev/stderr";
      exit(1);
    }
  }
  print;
}'
}

# Determine which lattice-info command should be used
lattice_info_cmd="$(which lattice-info 2> /dev/null)" ||
lattice_info_cmd=lattice-info/lattice-info;

if [ ! -f "$lattice_info_cmd" ]; then
  echo "lattice-info was not found in your system, compiling!" >&2;
  [ -d lattice-info/.git ] ||
  git clone https://github.com/jpuigcerver/lattice-info.git || exit 1;
  ( cd lattice-info && make && cd .. ) || exit 1;
fi;

# Remove old symbol tables if they exist.
rm -f lattice.word.sym iam.word.sym;

# Convert character-level to word-level lattice
../../lattice-expand-subpaths \
  --symbol-table=lattice.word.sym \
  --symbol-table-text=true \
  --verbose=2 \
  3 ark:lattice.char.txt ark,t:lattice.word.txt;
echo "";

# Both lattices should have the same number of paths!
echo "CHECK NUMBER OF PATHS FOR THE TINY LATTICE...";
echo -n "Char lattice ";
"$lattice_info_cmd" ark:lattice.char.txt 2> /dev/null |
awk '$0 ~ /avg. of paths/{ print $1, $2, $3":", $4 }';
echo -n "Word lattice ";
"$lattice_info_cmd" ark:lattice.word.txt 2> /dev/null |
awk '$0 ~ /avg. of paths/{ print $1, $2, $3":", $4 }';
echo "";

# All paths should have the same cost!
echo "CHECK PATHS AND COSTS FOR THE TINY LATTICE...";
echo "Paths from char lattice:";
cat lattice.char.txt | get_lattice_nbest 3;
echo "Paths from word lattice:";
cat lattice.word.txt | get_lattice_nbest 3 | nbest_word_to_char lattice.word.sym;
echo "";


# Convert character to word-level lattices in IAM
../../lattice-expand-subpaths \
  --symbol-table=iam.word.sym \
  --symbol-table-text=true \
  --verbose=2 \
  "65 66 67 68 69 70 71 74 75 76 77 78 79" \
  "ark:zcat iam.char.ark.gz|" \
  "ark:|gzip -9 > iam.word.ark.gz";
echo "";

# Check number of paths in the IAM lattices
echo "CHECK NUMBER OF PATHS FOR IAM LATTICES...";
echo -n "Char lattice ";
"$lattice_info_cmd" "ark:zcat iam.char.ark.gz|" 2> /dev/null |
awk '$0 ~ /avg. of paths/{ print $1, $2, $3":", $4 }';
echo -n "Word lattice ";
"$lattice_info_cmd" "ark:zcat iam.word.ark.gz|" 2> /dev/null |
awk '$0 ~ /avg. of paths/{ print $1, $2, $3":", $4 }';
echo "";

# Check score of the best path in the IAM lattices
echo "CHECK BEST PATH SCORE FOR IAM LATTICES...";
echo "Best path score from char lattice:";
lattice-best-path "ark:zcat iam.char.ark.gz|"  2>&1 | awk '$0 ~ /best cost/{ print $(NF - 3);}'
echo "Best path score from word lattice:";
lattice-best-path "ark:zcat iam.word.ark.gz|"  2>&1 | awk '$0 ~ /best cost/{ print $(NF - 3);}'
echo "";

exit 0;
