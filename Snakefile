import os

configfile: "files/config.yaml"
#TODO account for gene annotation files.
rule all:
    input:
        os.path.join(config['output']['dir'], 'L_var_w_mito_and_ribo.fa'),
        os.path.join(config['output']['dir'], 'L_var_w_mito_and_ribo.gff'),


rule format_rRNA:
    input:
        gff=config['input']['rRNA']
    output:
        gff=os.path.join(config['output']['dir'], 'rRNA.gff')
    script:
        "scripts/format_rRNA.py"

rule create_tmp_dir:
    output:
        config['output']['scratch']
    shell:
        "mkdir {output}"

rule gffs_to_bed:
    input:
        rRNA=config['input']['rRNA'],
        genome=config['input']['gff']
    output:
        rRNA=temp(os.path.join(config['output']['scratch'], 'rRNA.bed')),
        genome=temp(os.path.join(config['output']['scratch'], 'genome.bed'))
    shell:
        "gff2bed < {input.rRNA} > {output.rRNA}; "
        "gff2bed < {input.genome} > {output.genome}"

rule rRNA_overlap:
    input:
        rRNA=os.path.join(config['output']['scratch'], 'rRNA.bed'),
        genome=os.path.join(config['output']['scratch'], 'genome.bed')
    output:
        overlap=temp(os.path.join(config['output']['scratch'], 'overlap.bed'))
    shell:
        "bedtools intersect -a {input.genome} -b {input.rRNA} | "
        "awk '$8==\"gene\"' > {output.overlap}"

rule filter_overlaps:
    input:
        gff=config['input']['gff'],
        overlaps=os.path.join(config['output']['scratch'], 'overlap.bed')
    output:
        gff=temp(os.path.join(config['output']['scratch'], 'filtered.gff'))
    script:
        "scripts/filter_overlaps.py"

rule combine_annos:
    input:
        genome=os.path.join(config['output']['scratch'], 'filtered.gff'),
        rrna=os.path.join(config['output']['dir'], 'rRNA.gff'),
        mito=os.path.join(config['input']['mito_gff'])
    output:
        combined=os.path.join(config['output']['dir'],
                              'L_var_w_mito_and_ribo.gff')
    shell:
        "grep -v --no-filename '^#' {input.genome} {input.rrna} {input.mito} "
        "> {output.combined}"

rule combine_fastas:
    input:
        genome=os.path.join(config['input']['genome']),
        mito=os.path.join(config['input']['mito'])
    output:
        combined=os.path.join(config['output']['dir'],
                              'L_var_w_mito_and_ribo.fa')
    shell:
        "cat {input.genome} {input.mito} > {output.combined}"

# agat is not actually installed on the cluster, but can be accessed here:
# https://github.com/NBISweden/AGAT
rule sanitize_gff:
    input:
        gff=os.path.join(config['output']['dir'],
                         'L_var_w_mito_and_ribo.gff')
    output:
        gff=os.path.join(config['output']['dir'],
                         'L_var_w_mito_and_ribo_clean.gff')
    shell:
        "agat_convert_sp_gxf2gxf.pl -g {input.gff} -o {output.gff}"

rule gff2gtf:
    input:
        gff=os.path.join(config['output']['dir'],
                         'L_var_w_mito_and_ribo_clean.gff')
    output:
        gtf=os.path.join(config['output']['dir'],
                         'L_var_w_mito_and_ribo_clean.gtf')
    shell:
        "agat_convert_sp_gff2gtf.pl -gff {input.gff} -o {input.gtf}"