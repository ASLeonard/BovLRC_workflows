from pathlib import PurePath
import polars as pl

#workflow._singularity_args += f' -B /nfs/nas12.ethz.ch/fs1201/green_groups_tg_public/data -B $TMPDIR -B /cluster/work/pausch/inputs -B /cluster/spack/apps/'

minimap2_presets = {'PB':'map-hifi','ONT':'map-ont'}

df = pl.read_csv('long_read_metadata.csv')

def load_metadata_table():
    targets = []
    for row in df.to_dicts():

        targets.append(f'alignments/{row["SampleID"]}.{row["Technology"]}.mm2.stats')

        targets.append(f'clair3/{row["SampleID"]}_{row["Technology"]}/merge_output.gvcf.gz')
        targets.append(f'clair3/{row["SampleID"]}_{row["Technology"]}/merge_output.gvcf.gz.tbi')

        targets.append(f'sniffles2/{row["SampleID"]}.{row["Technology"]}.snf')

    return targets

rule all:
    input:
        load_metadata_table()

rule minimap2_index:
    input:
        reference = 'asset/genome_compact/ARS_UCD_v2.0.fa.gz'
    output:
        'asset/genome_compact/ARS_UCD_v2.0.{read}.mmi'
    params:
        preset = lambda wildcards: minimap2_presets[wildcards.read]
    threads: 1
    resources:
        mem_mb = 25000
    shell:
        '''
        minimap2 -x {params.preset} -d {output} {input}
        '''

rule minimap2_align:
    input:
        reference = 'asset/genome_compact/ARS_UCD_v2.0.fa.gz',
        bam = lambda wildcards: df.filter((pl.col('SampleID')==wildcards.sample)&(pl.col('Technology')==wildcards.read)).select('FASTQ_LR_Dir').item(),
        index = rules.minimap2_index.output
    output:
        #multiext('alignments/{sample}.{read}.mm2.bam','','.csi')
        multiext('alignments/{sample}.{read}.mm2.cram','','.crai')
    params:
        preset = lambda wildcards: minimap2_presets[wildcards.read],
        decoding_reference = lambda wildcards: df.filter((pl.col('SampleID')==wildcards.sample)&(pl.col('Technology')==wildcards.read)).select('Reference').item(),
        cram_fmt = lambda wildcards, input, output: f'--reference {input.reference} --output-fmt cram,version=3.1' if PurePath(output[0]).suffix == '.cram' else ''
    threads: lambda wildcards: 16
    resources:
        mem_mb = lambda wildcards: 5000 if wildcards.read == 'PB' else 7000,
        walltime = '24h',
        scratch = '50G'
    shell:
        '''
        samtools fastq --reference {params.decoding_reference} --threads {threads} {input.bam} |\
        minimap2 -t {threads} -ax {params.preset} -R '@RG\\tPL:{wildcards.read}\\tID:{wildcards.sample}\\tSM:{wildcards.sample}' --MD {input.index} - |\
        samtools sort - {params.cram_fmt} -m 3000M -@ {threads} -T $TMPDIR -o {output[0]} --write-index
        '''

rule cramino_stats:
    input:
        reference = 'asset/genome_compact/ARS_UCD_v2.0.fa.gz',
        bam = rules.minimap2_align.output
    output:
        'alignments/{sample}.{read}.mm2.stats'
    threads: 4
    resources:
        mem_mb = 2500,
        walltime = '4h'
    shell:
        '''
        cramino --reference {input.reference} -t {threads} -s {wildcards.sample}_{wildcards.read} {input.bam[0]} > {output}
        '''

rule clair3:
    input:
        reference = 'asset/genome_compact/ARS_UCD_v2.0.fa.gz',
        bam = rules.minimap2_align.output
    output:
        temp(expand('clair3/{{sample}}_{{read}}/{out}.vcf.gz{ext}',out=('full_alignment','merge_output','pileup'),ext=('','.tbi'))),
        temp(directory('clair3/{sample}_{read}/log')),
        temp('clair3/{sample}_{read}/run_clair3.log'),
        multiext('clair3/{sample}_{read}/merge_output.gvcf.gz','','.tbi')
    params:
        model = lambda wildcards: {'Sequel':'hifi_sequel2','Revio':'hifi_revio','r9_g3':'r941_prom_hac_g360+g422','r9_g4':'r941_prom_hac_g360+g422'}[df.filter((pl.col('SampleID')==wildcards.sample)&(pl.col('Technology')==wildcards.read)).select('Kit').item()],
        _dir = lambda wildcards, output: PurePath(output[0]).parent
    threads: 12
    resources:
        mem_mb = 6000,
        walltime = lambda wildcards: '24h' if wildcards.read == 'PB' else '72h',
        scratch = '200G'
    conda: 'clair3'
    shell:
        '''
        run_clair3.sh --bam_fn={input.bam[0]} \
                      --ref_fn={input.reference} \
                      --threads={threads} \
                      --platform=hifi \
                      --sample_name={wildcards.sample} \
                      --model_path=asset/models/{params.model} \
                      --include_all_ctgs \
                      --min_contig_size=40000000 \
                      --remove_intermediate_dir \
                      --gvcf \
                      --output={params._dir}
        '''

rule clair3_by_contig:
    input:
        reference = 'asset/genome_compact/ARS_UCD_v2.0.fa.gz',
        bam = rules.minimap2_align.output
    output:
        temp(expand('clair3/{{sample}}_{{read}}_{{contig}}/{out}.vcf.gz{ext}',out=('full_alignment','merge_output','pileup'),ext=('','.tbi'))),
        temp(directory('clair3/{sample}_{read}_{contig}/log')),
        temp('clair3/{sample}_{read}_{contig}/run_clair3.log'),
        multiext('clair3/{sample}_{read}_{contig}/merge_output.gvcf.gz','','.tbi')
    params:
        model = lambda wildcards: {'Sequel':'hifi_sequel2','Revio':'hifi_revio','r9_g3':'r941_prom_hac_g360+g422','r9_g4':'r941_prom_hac_g360+g422'}[df.filter((pl.col('SampleID')==wildcards.sample)&(pl.col('Technology')==wildcards.read)).select('Kit').item()],
        _dir = lambda wildcards, output: PurePath(output[0]).parent
    threads: 6
    resources:
        mem_mb = 5000,
        walltime = lambda wildcards: '24h' if wildcards.read == 'PB' else '24h',
        scratch = '20G'
    conda: 'clair3'
    shell:
        '''
        run_clair3.sh --bam_fn={input.bam[0]} \
                      --ref_fn={input.reference} \
                      --threads={threads} \
                      --ctg_name={wildcards.contig} \
                      --platform=hifi \
                      --sample_name={wildcards.sample} \
                      --model_path=asset/models/{params.model} \
                      --remove_intermediate_dir \
                      --gvcf \
                      --output={params._dir}
        '''

rule sniffles2:
    input:
        reference = 'asset/genome_compact/ARS_UCD_v2.0.fa.gz',
        bam = rules.minimap2_align.output
    output:
        vcf = temp(multiext('sniffles2/{sample}.{read}.vcf.gz','','.tbi')),
        snf = 'sniffles2/{sample}.{read}.snf'
    threads: 4
    resources:
        mem_mb = 3750
    conda: 'sniffles'
    shell:
        '''
        sniffles --input {input.bam[0]} \
         --sample-id {wildcards.sample} \
         --vcf {output.vcf[0]} \
         --snf {output.snf} \
         --threads {threads} \
         --minsvlen 50 \
         --mapq 20 \
         --reference {input.reference}
        '''
