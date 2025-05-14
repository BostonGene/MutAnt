#!/bin/bash
# Set to dbNSFP version to build
version="4.1a"

# define thread number for parallel processing where able
THREADS=7

do_prep() {
  extract_header
  custom_build_hg38
}

extract_header() {
  # grab header
  echo "Extracting header..."
  zcat dbNSFP${version}_variant.chr1.gz | head -n 1 | bgzip > header.gz

}

custom_build_hg38() {

  echo "Building hg38 version..."
  ### this section will produce data for hg38 capable pipelines
  ## hg38 version

  # Create a single file version
  # add header back into file
  cat header.gz dbNSFPv${version}_custom.gz > dbNSFPv${version}_custombuild.gz

  # NOTE: bgzip parameter -@ X represents number of threads
  cat dbNSFP${version}_variant.*.gz | zgrep -v '#chr' | bgzip -@ ${THREADS} >> dbNSFPv${version}_custombuild.gz

  # Create tabix index
  tabix -s 1 -b 2 -e 2 dbNSFPv${version}_custombuild.gz
}

do_prep
