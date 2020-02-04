'''
Index all genomes needed
'''
rule index_firstpass_genome:
    input:
        GENOME        
    output:
        expand(GENOME_IDX_PREFIX + '{idx}.bt2', idx = BT2_IDX_EXT)
    threads: THREADS
    params: GENOME_IDX_PREFIX
    shell:
        'bowtie2-build --threads {threads} {input} {params}'

rule index_secondpass_genome:
    input:
        os.path.join(DIR_SECONDPASS_GENOMES, '{SECONDPASS_GENOMES}.fa')
    output:
        os.path.join(DIR_SECONDPASS_GENOME_IDX + '{SECONDPASS_GENOMES}.1.bt2'),
        os.path.join(DIR_SECONDPASS_GENOME_IDX + '{SECONDPASS_GENOMES}.2.bt2'),
        os.path.join(DIR_SECONDPASS_GENOME_IDX + '{SECONDPASS_GENOMES}.3.bt2'),
        os.path.join(DIR_SECONDPASS_GENOME_IDX + '{SECONDPASS_GENOMES}.4.bt2'),
        os.path.join(DIR_SECONDPASS_GENOME_IDX + '{SECONDPASS_GENOMES}.rev.1.bt2'),
        os.path.join(DIR_SECONDPASS_GENOME_IDX + '{SECONDPASS_GENOMES}.rev.2.bt2')
    threads: THREADS
    params: DIR_SECONDPASS_GENOME_IDX + '{SECONDPASS_GENOMES}'
    shell:
        'bowtie2-build --threads {threads} {input} {params}'

rule check_index:
    input:
        expand(GENOME_IDX_PREFIX + '{idx}.bt2', idx = BT2_IDX_EXT),
        expand(os.path.join(DIR_SECONDPASS_GENOME_IDX + '{SECONDPASS_GENOMES}{idx}.bt2'), idx = BT2_IDX_EXT, SECONDPASS_GENOMES = SECONDPASS_GENOMES)
    output:
        touch(temp(os.path.join(DIR, 'index.done')))
