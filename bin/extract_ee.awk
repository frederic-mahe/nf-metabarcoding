#!/usr/bin/awk -f
#
# extract_ee.awk
#
# Read a vsearch fasta where headers carry both `ee=` (from --eeout)
# and `length=` (from --lengthout), and emit a three-column TSV:
#   <hash>  <ee>  <length>
#
# <hash> is the amplicon name vsearch assigned via --relabel_sha1 /
# --relabel_md5 (a SHA1 or MD5 digest, see [S65]). This splitter is
# delimiter-based and never measures the name, so it is agnostic to
# the hash width.
#
# Header layout:
#   >hash;ee=<float>;length=<int>
#
# Field separator splits on '>', ';', and '=':
#   $1 = "" (before '>')
#   $2 = hash
#   $3 = "ee" (literal key)
#   $4 = ee value
#   $5 = "length" (literal key)
#   $6 = length value
#
# Output is space-separated to match the existing downstream
# `sort` / `uniq` invocations.

BEGIN { FS = "[>;=]" }
/^>/ { print $2, $4, $6 }
