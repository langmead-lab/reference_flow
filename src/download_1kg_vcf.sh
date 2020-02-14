# wget -P resources/1kg_vcf/ -r -l1 --no-parent -A "vcf.gz" http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL/
wget -P resources/1kg_vcf/ http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL/ALL.chr{1..22}.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz
wget -P resources/1kg_vcf/ http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL/ALL.chr{X}.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz
# http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL/ALL.chr*.vcf.gz
