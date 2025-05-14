#!/bin/bash
# Set to dbNSFP version to build
version="4.6a"

# define thread number for parallel processing where able
THREADS=7

do_prep() {
  extract_header
  custom_build_hg19
}

extract_header() {
  # grab header
  echo "Extracting header..."
  zcat dbNSFP${version}_variant.chr1.gz | head -n 1 | bgzip > header.gz

}

custom_build_hg19() {
  ### this section will produce data for hg19 capable pipelines
  ## hg19 version
  # for hg19 (coordinate data is located in columns 8 [chr] and 9 [position])
  # this takes output from above, filters out any variants with no hg19 coords and then sorts on hg19 chr and position, and then bgzips output
  # NOTE: bgzip parameter -@ X represents number of threads
  echo "Building hg19 version..."

  zcat dbNSFPv${version}_custombuild.gz | \
  awk '$8 != "."' | \
  awk 'BEGIN{FS=OFS="\t"} {$1=$8 && $2=$9; NF--}1' | \
  LC_ALL=C sort --parallel=${THREADS} -n -S 4G -T . -k 1,1 -k 2,2 --compress-program=gzip | \
  bgzip -@ ${THREADS} > dbNSFPv${version}.hg19.custombuild.gz
  # NOTE: removed target memory allocation

  # Create tabix index
  tabix -s 1 -b 2 -e 2 dbNSFPv${version}.hg19.custombuild.gz

}

do_prep
