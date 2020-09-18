METHOD="${1:-randflow_ld}"
wget -P resources/ https://genome-idx.s3.amazonaws.com/bt/flow/${METHOD}.tar.gz
tar -zxvf resources/randflow_ld.tar.gz --directory snakemake/
