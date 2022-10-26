
# make conda software available
. /path/to/Anaconda/3/2022.05/etc/profile.d/conda.sh

conda activate /path/to/splice/pangolin_env

pangolin inst/extdata/spliceai_output.vcf \
  /path/to/data/gatk_bundle/hg19/ucsc.hg19.fasta \
  /path/to/data/human/gencode/34/gencode.v34lift37.annotation.db \
  inst/extdata/spliceai_output.pangolin
