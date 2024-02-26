from pathlib import PurePath
import pandas as pd

workflow._singularity_args += f' -B /nfs/nas12.ethz.ch/fs1201/green_groups_tg_public/data -B $TMPDIR -B /cluster/work/pausch/inputs'

minimap2_presets = {'PB':'map-hifi','ONT':'map-ont'}

def load_metadata_table():
    df = pd.read_csv('LR_metadata.csv')
    targets = []
    for _, row in df.iterrows():
        targets.append(f'clair3/{row["SampleID"]}_{row["Kit"]}/merge_output.gvcf.gz')
        targets.append(f'clair3/{row["SampleID"]}_{row["Kit"]}/merge_output.gvcf.gz.tbi')

        targets.append(f'sniffles2/{row["SampleID"]}.{row["Kit"]}.snf')
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
        mem_mb = 10000
    shell:
        '''
        minimap2 -x {params.preset} -d {output} {input}
        '''

rule minimap2_align:
    input:
        uBAM = lambda wildcards: rules.fibertools_predict_m6a.output[0] if wildcards.methylation == 'm6a' else 'alignments/uBAM/{sample}/{cell}.5mC.bam',
        rules.minimap2_index.output
    output:
        multiext('alignments/{sample}.read.mm2.cram','','.crai')
    params:
        preset = lambda wildcards: minimap2_presets[wildcards.read]
    threads: lambda wildcards: 24 if wildcards.cell[:2] == 'm8' else 16
    resources:
        mem_mb = 5000,
        walltime = '24h',
        scratch = '50G'
    shell:
        '''
        samtools fastq --threads {threads} {input.uBAM} |\
        minimap2 -t {threads} -ax {params.preset} -R '@RG\\tPL:PacBio\\tID:{wildcards.\\tSM:{wildcards.sample}' -y -Y {input.reference} - |\
        samtools sort - -m 3000M -@ {threads} -T $TMPDIR -o {output[0]} --write-index
        '''

rule clair3:
    input:
        reference = 'asset/genome_compact/ARS_UCD_v2.0.fa.gz',
        bam = '/nfs/nas12.ethz.ch/fs1201/green_groups_tg_public/data/BTA/bams_UCD2.0_eQTL_HiFi/{sample}.mm2.cram' #'bams/{sample}.cram'
    output:
        temp(expand('clair3/{{sample}}_{{read}}/{out}.vcf.gz{ext}',out=('full_alignment','merge_output','pileup'),ext=('','.tbi'))),
        temp(directory('clair3/{sample}_{read}/log')),
        temp('clair3/{sample}_{read}/run_clair3.log'),
        multiext('clair3/{sample}_{read}/merge_output.gvcf.gz','','.tbi')
    params:
        model = lambda wildcards: {'Sequel':'hifi_sequel2','Revio':'hifi_revio','r9_g3':'r941_prom_hac_g360+g422','r9_g4':'r941_prom_hac_g360+g422'}[wildcards.read],
        _dir = lambda wildcards, output: PurePath(output[0]).parent
    container: '/cluster/work/pausch/alex/software/images/clair3_v1.0.5.sif'
    threads: 12
    resources:
        mem_mb = 3000,
        walltime = '24h'
    shell:
        '''
        run_clair3.sh --bam_fn={input.bam} \
                      --ref_fn={input.reference} \
                      --threads={threads} \
                      --platform=hifi \
                      --sample_name={wildcards.sample} \
                      --model_path=/opt/models/{params.model} \
                      --include_all_ctgs \
                      --min_contig_size=40000000 \
                      --remove_intermediate_dir \
                      --gvcf \
                      --output={params._dir}
        '''

rule sniffles2:
    input:
        reference = 'asset/genome_compact/ARS_UCD_v2.0.fa.gz',
        bam = '/nfs/nas12.ethz.ch/fs1201/green_groups_tg_public/data/BTA/bams_UCD2.0_eQTL_HiFi/{sample}.mm2.cram' #'bams/{sample}.cram'
    output:
        vcf = 'sniffles2/{sample}.{read}.vcf.gz',
        snf = 'sniffles2/{sample}.{read}.snf'
    threads: 4
    resources:
        mem_mb = 3750
    conda: 'sniffles'
    shell:
        '''
        sniffles --input {input.bam} \
         --sample-id {wildcards.sample} \
         --vcf {output.vcf} \
         --snf {output.snf} \
         --threads {threads} \
         --minsvlen 50 \
         --mapq 20 \
         --reference {input.reference}
        '''
