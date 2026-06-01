#!/usr/bin/env bash
set -euo pipefail

infile=$1
extend_starts_bp=$2
extend_ends_bp=$3
starts_extended_file=$4
ends_extended_file=$5

# Use named pipes to avoid double decompression
mkfifo "${starts_extended_file}.fifo" "${ends_extended_file}.fifo" 2>/dev/null || true

# choose the right decompression command
if [[ $infile == *.gz ]]; then
  decompressor="zcat"
else
  decompressor="cat"
fi

# Single decompression, tee to both awk processes
$decompressor "$infile" | tee >(
  awk -v s=$extend_starts_bp '
    BEGIN { OFS="\t" }
    {
      new_start = $2 - s
      new_end = $2 + s
      printf("%s\t%d\t%d", $1, new_start, new_end)
      for(i=4; i<=NF; i++) printf("%s%s", OFS, $i)
      print ""
    }' > "$starts_extended_file"
) | awk -v e=$extend_ends_bp '
    BEGIN { OFS="\t" }
    {
      new_start = $3 - e
      new_end = $3 + e
      printf("%s\t%d\t%d", $1, new_start, new_end)
      for(i=4; i<=NF; i++) printf("%s%s", OFS, $i)
      print ""
    }' > "$ends_extended_file"

# Clean up FIFOs if they exist
rm -f "${starts_extended_file}.fifo" "${ends_extended_file}.fifo" 2>/dev/null || true

echo "Both jobs are done."
