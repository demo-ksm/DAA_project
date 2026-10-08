#!/bin/sh
# Re-run the (resumable) downloader until the final cache file exists.
cd "$(dirname "$0")/.."
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
n=0
while [ ! -f data/osm/chennai_roads_raw.rds ] && [ $n -lt 30 ]; do
  "$RS" scripts/02_download_osm.R
  n=$((n+1)); sleep 20
done
