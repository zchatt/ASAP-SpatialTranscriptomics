#!/bin/bash
#PBS -P FFbigdata
#PBS -N hybrid
#PBS -l select=1:ncpus=8:mem=16GB
#PBS -l walltime=4:00:00

#NOTE: modified from https://support.10xgenomics.com/single-cell-gene-expression/software/release-notes/build#GRCh38mm10_2020A

module load cellranger
source $CellRangerEnvFile

## input
cd /project/RDS-SMS-FFbigdata-RW/spatial_transcriptomics/xenograft/scrna_seq/reference

#################### SETUP ####################

human_genome="GRCh38"
rat_genome="mRatBN7"
version="2023-A"


build="GRCh38_and_mRatBN7-2023-A_build"
mkdir -p "$build"


# Download source files if they do not exist in reference_sources/ folder
source="reference_sources"
mkdir -p "$source"


human_fasta_url="http://ftp.ensembl.org/pub/release-98/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
human_fasta_in="${source}/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
human_gtf_url="http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/gencode.v32.primary_assembly.annotation.gtf.gz"
human_gtf_in="${source}/gencode.v32.primary_assembly.annotation.gtf"
human_fasta_cdna_url="http://ftp.ensembl.org/pub/release-98/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz"
human_fasta_cdna_in="${source}/Homo_sapiens.GRCh38.cdna.all.fa"

rat_fasta_url="ftp://ftp.ensembl.org/pub/release-110/fasta/rattus_norvegicus/dna/Rattus_norvegicus.mRatBN7.2.dna.toplevel.fa.gz"
rat_fasta_in="${source}/Rattus_norvegicus.mRatBN7.2.dna.toplevel.fa"
rat_gtf_url="https://ftp.ensembl.org/pub/release-110/gtf/rattus_norvegicus/Rattus_norvegicus.mRatBN7.2.110.gtf.gz"
rat_gtf_in="${source}/Rattus_norvegicus.mRatBN7.2.110.gtf"
rat_fasta_cdna_url="http://ftp.ensembl.org/pub/release-110/fasta/rattus_norvegicus/cdna/Rattus_norvegicus.mRatBN7.2.cdna.all.fa.gz"
rat_fasta_cdna_in="${source}/Rattus_norvegicus.mRatBN7.2.cdna.all.fa"

if [ ! -f "$human_fasta_in" ]; then
    curl -sS "$human_fasta_url" | zcat > "$human_fasta_in"
fi
if [ ! -f "$human_gtf_in" ]; then
    curl -sS "$human_gtf_url" | zcat > "$human_gtf_in"
fi
if [ ! -f "$human_fasta_cdna_in" ]; then
    curl -sS "$human_fasta_cdna_url" | zcat > "$human_fasta_cdna_in"
fi
if [ ! -f "$rat_fasta_in" ]; then
    curl -sS "$rat_fasta_url" | zcat > "$rat_fasta_in"
fi
if [ ! -f "$rat_gtf_in" ]; then
    curl -sS "$rat_gtf_url" | zcat > "$rat_gtf_in"
fi
if [ ! -f "$rat_fasta_cdna_in" ]; then
    curl -sS "$rat_fasta_cdna_url" | zcat > "$rat_fasta_cdna_in"
fi


# String patterns used for both genomes
ID="(ENS(MUS)?[GTE][0-9]+)\.([0-9]+)"


BIOTYPE_PATTERN=\
"(protein_coding|lncRNA|\
IG_C_gene|IG_D_gene|IG_J_gene|IG_LV_gene|IG_V_gene|\
IG_V_pseudogene|IG_J_pseudogene|IG_C_pseudogene|\
TR_C_gene|TR_D_gene|TR_J_gene|TR_V_gene|\
TR_V_pseudogene|TR_J_pseudogene)"
GENE_PATTERN="gene_type \"${BIOTYPE_PATTERN}\""
TX_PATTERN="transcript_type \"${BIOTYPE_PATTERN}\""
READTHROUGH_PATTERN="tag \"readthrough_transcript\""
PAR_PATTERN="tag \"PAR\""

GENE_PATTERN_ENSEMBL="gene_biotype \"${BIOTYPE_PATTERN}\"" # note that "gene_type" in gencode is "gene_biotype" in ensembl
TX_PATTERN_ENSEMBL="transcript_biotype \"${BIOTYPE_PATTERN}\"" # as above

#################### HUMAN ####################
# Please see the GRCh38-2020-A build documentation for details on these steps.


# Process FASTA -- translate chromosome names
human_fasta_modified="$build/$(basename "$human_fasta_in").modified"
cat "$human_fasta_in" \
    | sed -E 's/^>(\S+).*/>\1 \1/' \
    | sed -E 's/^>([0-9]+|[XY]) />chr\1 /' \
    | sed -E 's/^>MT />chrM /' \
    > "$human_fasta_modified"

# Process GTF -- split Ensembl IDs from version suffixes
human_gtf_modified="$build/$(basename "$human_gtf_in").modified"
cat "$human_gtf_in" \
    | sed -E 's/gene_id "'"$ID"'";/gene_id "\1"; gene_version "\3";/' \
    | sed -E 's/transcript_id "'"$ID"'";/transcript_id "\1"; transcript_version "\3";/' \
    | sed -E 's/exon_id "'"$ID"'";/exon_id "\1"; exon_version "\3";/' \
    > "$human_gtf_modified"


# Process GTF -- filter based on gene/transcript tags
cat "$human_gtf_modified" \
    | awk '$3 == "transcript"' \
    | grep -E "$GENE_PATTERN" \
    | grep -E "$TX_PATTERN" \
    | grep -Ev "$READTHROUGH_PATTERN" \
    | grep -Ev "$PAR_PATTERN" \
    | sed -E 's/.*(gene_id "[^"]+").*/\1/' \
    | sort \
    | uniq \
    > "${build}/gene_allowlist"


human_gtf_filtered="${build}/$(basename "$human_gtf_in").filtered"
grep -E "^#" "$human_gtf_modified" > "$human_gtf_filtered"
grep -Ff "${build}/gene_allowlist" "$human_gtf_modified" \
    >> "$human_gtf_filtered"


#################### MOUSE ####################
# Please see the mRatBN7-2023-A build documentation for details on these steps.

# Process FASTA -- translate chromosome names
# rat_fasta_modified="$build/$(basename "$rat_fasta_in").modified"
# cat "$rat_fasta_in" \
#     | sed -E 's/^>(\S+).*/>\1 \1/' \
#     | sed -E 's/^>([0-9]+|[XY]) />chr\1 /' \
#     | sed -E 's/^>MT />chrM /' \
#     > "$rat_fasta_modified"


# Process GTF -- split Ensembl IDs from version suffixes
rat_gtf_modified="$build/$(basename "$rat_gtf_in").modified"
cat "$rat_gtf_in" \
    | sed -E 's/gene_id "'"$ID"'";/gene_id "\1"; gene_version "\3";/' \
    | sed -E 's/transcript_id "'"$ID"'";/transcript_id "\1"; transcript_version "\3";/' \
    | sed -E 's/exon_id "'"$ID"'";/exon_id "\1"n; exon_version "\3";/' \
    > "$rat_gtf_modified"


# Process GTF -- filter based on gene/transcript tags
cat "$rat_gtf_modified" \
    | awk '$3 == "transcript"' \
    | grep -E "$GGENE_PATTERN_ENSEMBL" \
    | grep -E "$TX_PATTERN_ENSEMBL" \
    | grep -Ev "$READTHROUGH_PATTERN" \
    | sed -E 's/.*(gene_id "[^"]+").*/\1/' \
    | sort \
    | uniq \
    > "${build}/gene_allowlist"


rat_gtf_filtered="${build}/$(basename "$rat_gtf_in").filtered"
grep -E "^#" "$rat_gtf_modified" > "$rat_gtf_filtered"
grep -Ff "${build}/gene_allowlist" "$rat_gtf_modified" \
    >> "$rat_gtf_filtered"


#################### MKREF ####################

cellranger mkref --ref-version="$version" \
    --genome="$human_genome" --fasta="$human_fasta_modified" --genes="$human_gtf_filtered" \
    --genome="$rat_genome" --fasta="$rat_fasta_in" --genes="$rat_gtf_filtered"

#################### cat cdna together for kallisto indexing ####################

cat $human_fasta_cdna_in $rat_fasta_cdna_in > $build/GRCh38_and_mRatBN7-2023-A_cdna.fa
