# This script downloads GRCh38 reference genome from NCBI

wget -P resources/ ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
bgzip -d resources/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
