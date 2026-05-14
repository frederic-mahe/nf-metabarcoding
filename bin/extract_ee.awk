#!/usr/bin/awk -f
#
# extract_ee.awk
#
# Convert paired FASTA header/sequence lines (as produced by
# `paste - -` over a vsearch `--eeout` fasta) into a three-column TSV:
#   <SHA1>  <ee>  <length>
#
# Input record layout (after paste):
#   >SHA1;ee=<float>\t<sequence>
#
# Field separator splits on '>', ';', '=', and TAB:
#   $1 = "" (before '>')
#   $2 = SHA1
#   $3 = "ee" (literal key)
#   $4 = ee value
#   $NF = sequence
#
# Output is space-separated to match the existing downstream
# `sort` / `uniq` invocations.

BEGIN { FS = "[>;=\t]" }
{ print $2, $4, length($NF) }
