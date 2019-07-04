# PenGU-testdata-build
Simulate fastq files from a list of RefSeq URLs  


Given a .csv file of RefSeq asssembly URLs and read depths, and optionally a region defined in samtools faidx style, simulate fastq files at a range of read lengths.

`nextflow run PenGU-testdata-build.nf [OPTIONS]`



By default, submit to SLURM cluster. Change 'slurm' to 'local' in nextflow.config if you want to run locally.

```
--csv 		Input .csv file [default: 'genome_urls.csv']
--refreads 	SRR accession of a readset from a similar sequencer. Used for generating quality scores [default: 'SRR8062313']
--fq_header	String including "@" in fastq header [default: '@M04531:123:000000000-T3STP:1']
--read_lengths	Comma-delimited list of read lengths to simulate reads at [default: '125,150,175,200,225,250']
```

Requires an [NCBI_API_KEY](https://ncbiinsights.ncbi.nlm.nih.gov/2017/11/02/new-api-keys-for-the-e-utilities). 
Set `export NCBI_API_KEY=0123456789abcdef` somewhere in your environment (like `~/.bash_rc`).


This pipeline contains GPLv3 licensed code from [ArtificalFastqGenerator](https://sourceforge.net/projects/artfastqgen/).
