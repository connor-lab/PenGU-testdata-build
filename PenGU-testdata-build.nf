#!/usr/bin/env nextflow

params.csv = 'genome_urls.csv'
params.refreads = 'SRR8062313'
params.fq_header = '@M04531:123:000000000-T3STP:1'
params.read_lengths = '150,200,225,250'
params.outdir = 'output'


fq_header = params.fq_header
read_lengths = params.read_lengths.split(',')
outdir = params.outdir


Channel.fromPath(params.csv).splitCsv(header: true).set{ assemblyUrlsCsv }
Channel.from( read_lengths ).set{ readLengths }
Channel.fromSRA(params.refreads).set{ referenceFastq }


process DOWNLOAD_ASSEMBLY {
    tag {cols.Name}

    publishDir "${outdir}/assembly_fasta", mode: 'copy', pattern: "${cols.Name}.fa"

    input:
    val cols from assemblyUrlsCsv

    output:
    set val("${cols.Name}"), val("${cols.targetCov}"), val("${cols.ReadPrefix}"), file("${cols.Name}.fa"), val("${cols.extract_amplicon}"), val("${cols.amplicon_region}") into choiceExtractRegion
    
    script:
    """
    curl -fssL '${cols.RefSeqFTPURL}' | gzip -d > ${cols.Name}.fa
    """
}

Channel.create().set{ extractRegion }
Channel.create().set{ simulateReadsNoRegion }

choiceExtractRegion.choice( extractRegion, simulateReadsNoRegion) { it[4].toLowerCase() ==~ /true/ ? 0 : 1}

process EXTRACT_REGION {
    tag { name }

    input:
    set name, coverage, prefix, file(assembly), extract, region from extractRegion

    output:
    set val(name), val(coverage), val(prefix), file("${assembly}") into simulateReadsRegion

    script:

    """
    samtools faidx ${assembly} ${region} > ${name}.${region}.fa
    mv ${name}.${region}.fa ${assembly}
    """
}

simulateReadsRegion.mix(
    simulateReadsNoRegion.map { 
        [ it[0], it[1], it[2], it[3] ] 
        }
    ).set{
        getFirstFastaHeader
    }

process GETFIRSTFASTAHEADER {
    tag { name }

    input:
    set name, coverage, prefix, file(assembly) from getFirstFastaHeader

    output:
    set val(name), val(coverage), val(prefix), file("${assembly}"), stdout into simulateReads

    script:
    """
    grep -m 1 ">" ${assembly} | awk '{print \$1}' | tr -d '\\n'
    """
}


process STAGE_AND_DECOMPRESS_REF_FASTQ {
    tag { name }

    stageInMode 'copy'

    input:
    set file(forward), file(reverse) from referenceFastq.map{ [ it[1][0], it[1][1] ] }

    output:
    set file("*_1.fastq"), file ("*_2.fastq") into decompressedReferenceFastq

    script:
    """
    gzip -d *.fastq.gz
    """
}


process SIMULATE_READS {
    tag { name }

    cpus 4

    maxForks 10

    input:
    set name, coverage, prefix, file(assembly), header, length, file(forward), file(reverse) from simulateReads.combine(readLengths).combine(decompressedReferenceFastq)

    output:
    set name, length, prefix, file("${name}.1.fastq"), file("${name}.2.fastq") into simulatedReads

    script:
    """
    afg -O ${name} -R ${assembly} -F1 ${forward} -F2 ${reverse} -RL ${length} -URQS true -CMP ${coverage} -TLM 500 -N 1000 -S \'${header}\'
    """

}

process CLEAN_SIMULATED_FASTQ {
    tag { name }

    publishDir "${outdir}/fastq", mode: 'copy', pattern: "*_R?_001.fastq.gz"

    input:
    set name, length, prefix, file(forward), file(reverse) from simulatedReads

    output:
    set name, file("*_R1_001.fastq.gz"), file("*_R2_001.fastq.gz") 

    script:
    """
    sed "s/^@HWI-ST745_0097:7/${fq_header}/g" ${forward} | sed 's/#0\\/1/ 1:N:0:10/g' | sed '/^+HWI/s/.*/+/' | gzip > ${prefix}-${name}_S${length}_L001_R1_001.fastq.gz
    sed "s/^@HWI-ST745_0097:7/${fq_header}/g" ${reverse} | sed 's/#0\\/2/ 2:N:0:10/g' | sed '/^+HWI/s/.*/+/' | gzip > ${prefix}-${name}_S${length}_L001_R2_001.fastq.gz
    """
}
