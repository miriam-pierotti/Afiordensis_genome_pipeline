# *Antipathella fiordensis* reference genome pipeline
This repository contains the code and documentation for the pseudo-chromosomal reference genome and RNA-informed annotation of the black coral *Antipathella fiordensis*, assembled with Oxford Nanopore long-read and cDNA sequencing, described in the manuscript:

CITATION HERE

___

## General Pipeline Description

```mermaid
%%{init: {"flowchart": {"defaultRenderer": "elk"}} }%%
flowchart TB
    subgraph GENOME ASSEMBLY
        direction TB
        ONTDNA[ONT DNA / Raw Reads] -->|.bam| SAM1(SAMTOOLS)
        SAM1 --> FQ1(((.fastq)))
        FQ1 --> DOR1(DORADO)
        FQ1 --> NP1(NanoPlot)
        DOR1 --> CHOP1(CHOPPER)
        CHOP1 --> |Short reads >8kbp| DORCOR(DORADO CORRECT)
        DORCOR ---> |Corrected filtered reads| HIFI(HIFIASM)
        DORCOR --> NP2(NANOPLOT)
        DORCOR --> SEQKIT(SEQKIT)
        CHOP1 --> |Ultralong reads >10kbp| HIFI(HIFIASM)
        HIFI --> HIFISTAT(HIFIASM STATS)
        HIFI --> |Assembled contigs| FASTA(((.fasta)))
        FASTA --> PURDUP(PURGEDUPS)
        FQ1 --> PURDUP
        PURDUP ---> MED(MEDAKA)
        PURDUP --> PURDUPSTATS(PURGEDUP STATS)
        FQ1 --> MED
        MED --> MEDSTATS(MEDAKA STATS)
        MED ---> RTTG(RAGTAG)

        subgraph COMPARATIVE GENOMICS
            DGEN(D-GENIES) --> RTTG
        end

        RTTG --> RESCAF(REORDERED SCAFFOLDS)
        RESCAF --> RESCAFSTATS(REORDERED SCAFFOLDS STATS)
        RESCAF --> FCSGX(FCS_GX)
        RESCAF --> BLOBTOOLS(BLOBTOOLS)

        subgraph  CONTAMINATION SCREENING
           FCSGX
           BLOBTOOLS 
        end
    end

    subgraph MITOGENOME ASSEMBLY
        direction TB
        RESCAF ---> MITHIFI(MitoHIFI)
        MITHIFI --> RAC(RACON)
        RAC --> MITFINDER(MITOFINDER)
    end

    subgraph GENOME ANNOTATION
        direction TB
        ONTRNA[ONT RNA / Raw Reads] -->|.bam| SAM2(SAMTOOLS)
        SAM2 --> FQ2(((.fastq)))
        FQ2 --> NP3(NANOPLOT)
        FQ2 --> DOR2(DORADO)
        DOR2 --> CHOP2(CHOPPER)
        CHOP2 --> MINMAP(MINIMAP)
        MINMAP --> RESCAF
        RESCAF --> BAMBU(BAMBU)
        
        subgraph STRUCTURAL ANNOTATION
            RESCAF --> REPMOD(REPEATMODELER)
            REPMOD --> BUS(BUSCO)
            REPMOD --> REPMAS(REPEATMASKER)
            REPMAS --> BRAK(BRAKER)
        end

        BRAK -->|.gtf file|BAMBU
        MINMAP -->|transcriptom aligned reads|BAMBU
        BAMBU -->|.gtf file|TRANS(TRANSDECODER)
        TRANS -->|.gff3 file|EGG(EGGNOG)
        TRANS -->|.gff3 file|INT(INTERPROSCAN)
    end
```

### Genome Assembly

The pseudo-chomosomal reference genome assembly is fully documented in [GENOME ASSEMBLY DOCUMENTATION](GENOME_ASSEMBLY.md).

### Annotation

The RNA-informed annotation is fully documented in the [ANNOTATION DOCUMENTATION](ANNOTATION.md).