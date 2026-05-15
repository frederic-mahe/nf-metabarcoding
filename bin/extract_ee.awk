#!/usr/bin/awk -f
#
# extract_ee.awk
#
# Read a vsearch fasta where headers carry both `ee=` (from --eeout)
# and `length=` (from --lengthout), and emit a three-column TSV:
#   <SHA1>  <ee>  <length>
#
# Header layout:
#   >SHA1;ee=<float>;length=<int>
#
# Field separator splits on '>', ';', and '=':
#   $1 = "" (before '>')
#   $2 = SHA1
#   $3 = "ee" (literal key)
#   $4 = ee value
#   $5 = "length" (literal key)
#   $6 = length value
#
# Output is space-separated to match the existing downstream
# `sort` / `uniq` invocations.

BEGIN { FS = "[>;=]" }
/^>/ { print $2, $4, $6 }
