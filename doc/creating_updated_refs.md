## Getting latest ENSEMBL fastas and GTFs for model organisms

2022-03-24

Setup. Define some functions and bash variables


    REF_DIR=/nfs/turbo/umms-brcfpipeline/references/ENSEMBL_genomes_105

    prep_refdir () {
        SPECIES_UC=${SPECIES^} #Uppercase first letter of species e.g. homo_sapiens -> Homo_sapiens
        cd $REF_DIR &&
        if [[ ! -d $SPECIES_UC ]]
            then mkdir $SPECIES_UC
        fi &&
        cd $SPECIES_UC &&
        mkdir $GENOME_BUILD &&
        cd $GENOME_BUILD
    }

    get_refs () {
        echo "Downloading fasta"
        wget ftp://ftp.ensembl.org/pub/release-${RELEASE}/fasta/${SPECIES}/dna/${SPECIES_UC}.${GENOME_BUILD}.${FASTATYPE}.fa.gz
        echo "Downloading GTF"
        wget ftp://ftp.ensembl.org/pub/release-${RELEASE}/gtf/${SPECIES}/${SPECIES_UC}.${GENOME_BUILD}.${RELEASE}.gtf.gz
        echo "Unzipping"
        gunzip *.gz
    }

Human GRCh38:

    SPECIES=homo_sapiens
    GENOME_BUILD=GRCh38
    RELEASE=105
    FASTATYPE="dna_sm.primary_assembly"
    prep_refdir
    get_refs

Human GRCh37:

    # remains unchanged - hardlinked from the previous ENSEMBL 98 refs

Mouse GRCm39:

    SPECIES=mus_musculus
    GENOME_BUILD=GRCm39
    RELEASE=105
    FASTATYPE="dna_sm.toplevel"
    prep_refdir
    get_refs

Mouse GRCm38:

    SPECIES=mus_musculus
    GENOME_BUILD=GRCm38
    RELEASE=102
    FASTATYPE="dna_sm.primary_assembly"
    prep_refdir
    get_refs

Mouse NCBIM37:

    # remains unchanged - hardlinked from the previous ENSEMBL 98 refs

Rat Rnor_6.0:

    SPECIES=rattus_norvegicus
    GENOME_BUILD=Rnor_6.0
    RELEASE=105
    FASTATYPE="dna_sm.toplevel"
    prep_refdir
    get_refs


Zebrafish

    SPECIES=danio_rerio
    GENOME_BUILD=GRCz11
    RELEASE=105
    FASTATYPE="dna_sm.primary_assembly"
    prep_refdir
    get_refs


C. elegans

    SPECIES=caenorhabditis_elegans
    GENOME_BUILD=WBcel235
    RELEASE=105
    FASTATYPE="dna_sm.toplevel"
    prep_refdir
    get_refs


Fruit fly

    SPECIES=drosophila_melanogaster
    GENOME_BUILD=BDGP6.32
    RELEASE=105
    FASTATYPE="dna_sm.toplevel"
    prep_refdir
    get_refs


## Creating annotation tables

You can use the Rscript included with Watermelon to gather current (or past) annotation information from biomaRt into a TSV file, e.g.:

    Rscript /path/to/Watermelon/scripts/ensembl_biomaRt_mapping.R --outfile Homo_sapiens/GRCh38/Homo_sapiens.GRCh38.105_annotation.tsv --dataset hsapiens_gene_ensembl


## Notes on getting ENSEMBL references for older genome builds

Can have slightly different steps depending on organism/build. These were ran somewhat manually. Note: these were not part of the update to ENSEMBL v105.

Human - GRCh37  (~ hg19)

    SPECIES=homo_sapiens
    SPECIES_UC=${SPECIES^}
    GENOME_BUILD=GRCh37
    RELEASE=75

    #Fasta difference - includes release num
    wget ftp://ftp.ensembl.org/pub/release-${RELEASE}/fasta/${SPECIES}/dna/${SPECIES_UC}.${GENOME_BUILD}.${RELEASE}.dna_sm.toplevel.fa.gz
    wget ftp://ftp.ensembl.org/pub/release-${RELEASE}/gtf/${SPECIES}/${SPECIES_UC}.${GENOME_BUILD}.${RELEASE}.gtf.gz
    gunzip *.gz

    #Also need to use archive host, and request slightly different attributes to work with the older biomart release
    Rscript ~/APtravis/Watermelon/scripts/ensembl_biomaRt_mapping.R --outfile Homo_sapiens.GRCh37.75_annotation.tsv --dataset hsapiens_gene_ensembl --host feb2014.archive.ensembl.org --attributes ensembl_gene_id,entrezgene,external_gene_id,description

Mouse - NCBIM37 (~ mm9)

    SPECIES=mus_musculus
    SPECIES_UC=${SPECIES^}
    GENOME_BUILD=NCBIM37
    RELEASE=67

    #Fasta difference - includes release num, and soft-masked not available
    wget ftp://ftp.ensembl.org/pub/release-${RELEASE}/fasta/${SPECIES}/dna/${SPECIES_UC}.${GENOME_BUILD}.${RELEASE}.dna.toplevel.fa.gz
    wget ftp://ftp.ensembl.org/pub/release-${RELEASE}/gtf/${SPECIES}/${SPECIES_UC}.${GENOME_BUILD}.${RELEASE}.gtf.gz
    gunzip *.gz

    #Also need to use archive host, and request slightly different attributes to work with the older biomart release
    Rscript ~/APtravis/Watermelon/scripts/ensembl_biomaRt_mapping.R --outfile Mus_musculus.NCBIM37.67_annotation.tsv --dataset mmusculus_gene_ensembl --host may2012.archive.ensembl.org --attributes ensembl_gene_id,entrezgene,external_gene_id,description
