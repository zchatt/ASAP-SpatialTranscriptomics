#!/bin/bash
# ### Kallisto for snRNAseq
# outlines here - https://github.com/pachterlab/kb_python 
# # NOTES: In order to generate the index file you will need to create a cDNA reference (transcriptome.fa) of each hybrid genome

# # pip install kb-python
export PATH=/home/zac/.local/bin/:$PATH
run_dir="/data/zac/asap_snrnaseq/snrnaseq/kallisto"
cd $run_dir

# # building references
# kb ref -i refdata-gex-GRCh38-and-mm10-2020-A.idx -g refdata-gex-GRCh38-and-mm10-2020-A_t2g.txt -f1 [transcriptome.fa] \
#  /data/zac/asap_snrnaseq/reference/refdata-gex-GRCh38-and-mm10-2020-A/fasta/genome.fa \
#  /data/zac/asap_snrnaseq/reference/refdata-gex-GRCh38-and-mm10-2020-A/genes/genes.gtf

# Build - screen -r 407996

# ## human-mouse
# kb ref -i refdata-gex-GRCh38-and-mm10-2020-A.idx \
#    -g refdata-gex-GRCh38-and-mm10-2020-A_t2g.txt \
#    -f1 transcriptome.fa \
#  /data/zac/asap_snrnaseq/reference/refdata-gex-GRCh38-and-mm10-2020-A/fasta/genome.fa \
#  /data/zac/asap_snrnaseq/reference/refdata-gex-GRCh38-and-mm10-2020-A/genes/genes.gtf

# ## human-rat
# kb ref -i refdata-gex-GRCh38_and_mRatBN7_281023.idx \
#    -g refdata-gex-GRCh38_and_mRatBN7_281023_t2g.txt \
#    -f1 transcriptome_GRCh38_and_mRatBN7_281023.fa \
#  /data/zac/asap_snrnaseq/reference/GRCh38_and_mRatBN7_281023/fasta/genome.fa \
#  /data/zac/asap_snrnaseq/reference/GRCh38_and_mRatBN7_281023/genes/genes.gtf.gz

# # # 12w_HJ377DRX2 - human/rat 
# r1="/data/zac/asap_snrnaseq/snrnaseq/proj-9520_10x_axoalguidance-1128.4.666/AGRF_CAGRF221012416_HJ377DRX2/12w/12w_S1_L001_R1_001.fastq.gz"
# r2="/data/zac/asap_snrnaseq/snrnaseq/proj-9520_10x_axoalguidance-1128.4.666/AGRF_CAGRF221012416_HJ377DRX2/12w/12w_S1_L001_R2_001.fastq.gz"
# ID="HJ377DRX2_12w_hr"
# INDEX="/data/zac/asap_snrnaseq/snrnaseq/kallisto/refdata-gex-GRCh38_and_mRatBN7_281023.idx"
# T2G="/data/zac/asap_snrnaseq/snrnaseq/kallisto/refdata-gex-GRCh38_and_mRatBN7_281023_t2g.txt"

# kb count -i $INDEX -g $T2G -o $ID -x 10xv3 $r1 $r2

# # # 6w_HJ377DRX2 - human/rat
# r1="/data/zac/asap_snrnaseq/snrnaseq/proj-9520_10x_axoalguidance-1128.4.666/AGRF_CAGRF221012416_HJ377DRX2/6w/6w_S1_L001_R1_001.fastq.gz"
# r2="/data/zac/asap_snrnaseq/snrnaseq/proj-9520_10x_axoalguidance-1128.4.666/AGRF_CAGRF221012416_HJ377DRX2/6w/6w_S1_L001_R2_001.fastq.gz"
# ID="HJ377DRX2_6w_hr"
# INDEX="/data/zac/asap_snrnaseq/snrnaseq/kallisto/refdata-gex-GRCh38_and_mRatBN7_281023.idx"
# T2G="/data/zac/asap_snrnaseq/snrnaseq/kallisto/refdata-gex-GRCh38_and_mRatBN7_281023_t2g.txt"

# kb count -i $INDEX -g $T2G -o $ID -x 10xv3 $r1 $r2

# # 6w_HLNGJDRX2 - human/rat with 100 cycle - try V2 chemistry before attempting technology string - https://github.com/pachterlab/kallisto/issues/287
r1="/data/zac/asap_snrnaseq/snrnaseq/proj-9520_10x_axoalguidance-1128.4.666/AGRF_CAGRF230213424_HLNGJDRX2/6w/6w_S1_L001_R1_001.fastq.gz"
r2="/data/zac/asap_snrnaseq/snrnaseq/proj-9520_10x_axoalguidance-1128.4.666/AGRF_CAGRF230213424_HLNGJDRX2/6w/6w_S1_L001_R2_001.fastq.gz"
r3="/data/zac/asap_snrnaseq/snrnaseq/proj-9520_10x_axoalguidance-1128.4.666/AGRF_CAGRF230213424_HLNGJDRX2/6w/6w_S1_L002_R1_001.fastq.gz"
r4="/data/zac/asap_snrnaseq/snrnaseq/proj-9520_10x_axoalguidance-1128.4.666/AGRF_CAGRF230213424_HLNGJDRX2/6w/6w_S1_L002_R2_001.fastq.gz"
ID="HLNGJDRX2_6w_hr"
INDEX="/data/zac/asap_snrnaseq/snrnaseq/kallisto/refdata-gex-GRCh38_and_mRatBN7_281023.idx"
T2G="/data/zac/asap_snrnaseq/snrnaseq/kallisto/refdata-gex-GRCh38_and_mRatBN7_281023_t2g.txt"

kb count --strand reverse -i $INDEX -g $T2G -o $ID -x 10xv2 $r1 $r2 $r3 $r4 


# # 1y_HLNGJDRX2 - human/mouse 
r1="/data/zac/asap_snrnaseq/snrnaseq/proj-9520_10x_axoalguidance-1128.4.666/AGRF_CAGRF230213424_HLNGJDRX2/1y/1y_S1_L001_R1_001.fastq.gz"
r2="/data/zac/asap_snrnaseq/snrnaseq/proj-9520_10x_axoalguidance-1128.4.666/AGRF_CAGRF230213424_HLNGJDRX2/1y/1y_S1_L001_R2_001.fastq.gz"
r3="/data/zac/asap_snrnaseq/snrnaseq/proj-9520_10x_axoalguidance-1128.4.666/AGRF_CAGRF230213424_HLNGJDRX2/1y/1y_S1_L002_R1_001.fastq.gz"
r4="/data/zac/asap_snrnaseq/snrnaseq/proj-9520_10x_axoalguidance-1128.4.666/AGRF_CAGRF230213424_HLNGJDRX2/1y/1y_S1_L002_R2_001.fastq.gz"
INDEX="/data/zac/asap_snrnaseq/snrnaseq/kallisto/refdata-gex-GRCh38-and-mm10-2020-A.idx"
T2G="/data/zac/asap_snrnaseq/snrnaseq/kallisto/refdata-gex-GRCh38-and-mm10-2020-A_t2g.txt"
ID="HLNGJDRX2_1y_hm"
kb count -i $INDEX -g $T2G -o $ID -x 10xv2 $r1 $r2 $r3 $r4
