# Commands for Illumina BaseSpace Sequence Hub command line interface (CLI)

cd ~

# authenticate account
$HOME/bs auth --api-server https://api.aps2.sh.basespace.illumina.com

# get config details
$HOME/bs load config 
#export BASESPACE_API_SERVER="https://api.aps2.sh.basespace.illumina.com"
#export BASESPACE_ACCESS_TOKEN="2ba3b51fca60457d878b974e49661cfd"

# export BASESPACE_API_SERVER="https://api.aps2.sh.basespace.illumina.com/"
# export BASESPACE_ACCESS_TOKEN="6c48fb1707d6456d8be5b1ae0986f1bd"

# find details of whoami
$HOME/bs whoami

# create project
#$HOME/bs create project --name=CHA12467 --exist-ok --api-server=https://api.aps2.sh.basespace.illumina.com/ --access-token=2ba3b51fca60457d878b974e49661cfd
$HOME/bs create project --name=geomx_281023_PB --exist-ok --api-server=https://api.aps2.sh.basespace.illumina.com/ --access-token=2ba3b51fca60457d878b974e49661cfd
$HOME/bs create project --name=geomx_281023_PF --exist-ok --api-server=https://api.aps2.sh.basespace.illumina.com/ --access-token=2ba3b51fca60457d878b974e49661cfd
$HOME/bs create project --name=geomx_281023_PEG --exist-ok --api-server=https://api.aps2.sh.basespace.illumina.com/ --access-token=2ba3b51fca60457d878b974e49661cfd
$HOME/bs create project --name=geomx_281023_PAB --exist-ok --api-server=https://api.aps2.sh.basespace.illumina.com/ --access-token=2ba3b51fca60457d878b974e49661cfd
$HOME/bs create project --name=geomx_281023_round1 --exist-ok --api-server=https://api.aps2.sh.basespace.illumina.com/ --access-token=2ba3b51fca60457d878b974e49661cfd
$HOME/bs create project --name=geomx_281023_PH --exist-ok --api-server=https://api.aps2.sh.basespace.illumina.com/ --access-token=2ba3b51fca60457d878b974e49661cfd

# see projects
$HOME/bs list projects --api-server=https://api.aps2.sh.basespace.illumina.com/ --access-token=2ba3b51fca60457d878b974e49661cfd

# # may need to combine technical replicates
# cd /Users/zacc/USyd/spatial_transcriptomics/analysis/geomx/geomx_oct2023
# mkdir -p merged_reps

# # combine technical replicates
# for i in $(awk 'NR>1 {print $1}' Annotations_remapped.txt | grep "649" | grep "11"); do
# 	echo $i
# 	du -sh */${i}*_R1_*.fastq.gz

# 	# cat files
# 	cat */${i}*_R1_*.fastq.gz > merged_reps/${i}_S0_L001_R1_001.fastq.gz
# 	cat */${i}*_R2_*.fastq.gz > merged_reps/${i}_S0_L001_R2_001.fastq.gz
# done

# # may need to rename files eg.
# cd /Users/zacc/CHA12602r2
# rename 's/_L(\d+)_R(\d+)_(\d+)\.fastq\.gz/|L$1|R$2|$3.fastq.gz/' *_L*_R*_*.fastq.gz
# rename 's/_/-/g' *.fastq.gz
# rename 's/\|/_/g' *.fastq.gz
# rename 's/-S/_S/g' *.fastq.gz

# upload datasets
$HOME/bs upload dataset --project=6206201 --api-server=https://api.aps2.sh.basespace.illumina.com/ --access-token=2ba3b51fca60457d878b974e49661cfd --concurrency=low --exclude '*' --include '*.fastq.gz' /Users/zacc/CHA12467

# OR upload datasets for large number of files 
# 1) move fastq to relevent directory
mkdir geomx_281023_PB
mv *1001660021814* geomx_281023_PB
mv *1001660021813* geomx_281023_PB
mv *1001660021812* geomx_281023_PB

mkdir geomx_281023_PF
mv *1001660021820* geomx_281023_PF
mv *1001660021821* geomx_281023_PF

mkdir geomx_281023_PEG
mv *1001660021819* geomx_281023_PEG

mkdir geomx_281023_PAB
mv *1001660021703* geomx_281023_PAB
mv *1001660021704* geomx_281023_PAB

mkdir geomx_281023_PH
mv *1001660021822* geomx_281023_PH

mkdir geomx_281023_round1
mv *1001660018649* geomx_281023_round1

#### upload

# geomx_281023_round1: upload Y, Run Y.
$HOME/bs whoami --api-server=https://api.aps2.sh.basespace.illumina.com/ --access-token=6c48fb1707d6456d8be5b1ae0986f1bd
$HOME/bs list projects --api-server=https://api.aps2.sh.basespace.illumina.com/ --access-token=6c48fb1707d6456d8be5b1ae0986f1bd
$HOME/bs create project --name=geomx_281023_round1 --exist-ok --api-server=https://api.aps2.sh.basespace.illumina.com/ --access-token=6c48fb1707d6456d8be5b1ae0986f1bd
cd ~/USyd/spatial_transcriptomics/analysis/geomx/geomx_oct2023/merged_reps/geomx_281023_round1
$HOME/bs upload dataset --project=6640635 --api-server=https://api.aps2.sh.basespace.illumina.com/ --access-token=6c48fb1707d6456d8be5b1ae0986f1bd --concurrency=low --exclude '*' --include '*.fastq.gz' .

# geomx_281023_PB: upload Y, Run Y.
$HOME/bs create project --name=geomx_281023_PB --exist-ok --api-server=https://api.aps2.sh.basespace.illumina.com/ --access-token=6c48fb1707d6456d8be5b1ae0986f1bd
cd ~/USyd/spatial_transcriptomics/analysis/geomx/geomx_oct2023/merged_reps/geomx_281023_PB
$HOME/bs upload dataset --project=6640636 --api-server=https://api.aps2.sh.basespace.illumina.com/ --access-token=6c48fb1707d6456d8be5b1ae0986f1bd --concurrency=low --exclude '*' --include '**1001660021812*' .
$HOME/bs upload dataset --project=6640636 --api-server=https://api.aps2.sh.basespace.illumina.com/ --access-token=6c48fb1707d6456d8be5b1ae0986f1bd --concurrency=low --exclude '*' --include '**1001660021813*' .
$HOME/bs upload dataset --project=6640636 --api-server=https://api.aps2.sh.basespace.illumina.com/ --access-token=6c48fb1707d6456d8be5b1ae0986f1bd --concurrency=low --exclude '*' --include '**1001660021814*' .

# geomx_281023_PF: upload Y, Run Y
$HOME/bs create project --name=geomx_281023_PF --exist-ok --api-server=https://api.aps2.sh.basespace.illumina.com/ --access-token=6c48fb1707d6456d8be5b1ae0986f1bd
cd ~/USyd/spatial_transcriptomics/analysis/geomx/geomx_oct2023/merged_reps/geomx_281023_PF
#$HOME/bs upload dataset --project=6645639 --api-server=https://api.aps2.sh.basespace.illumina.com/ --access-token=6c48fb1707d6456d8be5b1ae0986f1bd --concurrency=low --exclude '*' --include '**1001660021820*' .
$HOME/bs upload dataset --project=6645639 --api-server=https://api.aps2.sh.basespace.illumina.com/ --access-token=6c48fb1707d6456d8be5b1ae0986f1bd --concurrency=low --exclude '*' --include '**1001660021820-F-F11*' .
$HOME/bs upload dataset --project=6645639 --api-server=https://api.aps2.sh.basespace.illumina.com/ --access-token=6c48fb1707d6456d8be5b1ae0986f1bd --concurrency=low --exclude '*' --include '**1001660021820-F-G*' .
$HOME/bs upload dataset --project=6645639 --api-server=https://api.aps2.sh.basespace.illumina.com/ --access-token=6c48fb1707d6456d8be5b1ae0986f1bd --concurrency=low --exclude '*' --include '**1001660021820-F-H*' .
#$HOME/bs upload dataset --project=6645639 --api-server=https://api.aps2.sh.basespace.illumina.com/ --access-token=6c48fb1707d6456d8be5b1ae0986f1bd --concurrency=low --exclude '*' --include '**1001660021821*' .

# geomx_281023_PEG: upload Y, Run Y
$HOME/bs create project --name=geomx_281023_PEG --exist-ok --api-server=https://api.aps2.sh.basespace.illumina.com/ --access-token=6c48fb1707d6456d8be5b1ae0986f1bd
cd ~/USyd/spatial_transcriptomics/analysis/geomx/geomx_oct2023/merged_reps/geomx_281023_PEG
#$HOME/bs upload dataset --project=6645639 --api-server=https://api.aps2.sh.basespace.illumina.com/ --access-token=6c48fb1707d6456d8be5b1ae0986f1bd --concurrency=low --exclude '*' --include '**1001660021819*' .
$HOME/bs upload dataset --project=6645639 --api-server=https://api.aps2.sh.basespace.illumina.com/ --access-token=6c48fb1707d6456d8be5b1ae0986f1bd --concurrency=low --exclude '*' --include '**1001660021819-E-G10*' .
$HOME/bs upload dataset --project=6645639 --api-server=https://api.aps2.sh.basespace.illumina.com/ --access-token=6c48fb1707d6456d8be5b1ae0986f1bd --concurrency=low --exclude '*' --include '**1001660021819-E-G11*' .
$HOME/bs upload dataset --project=6645639 --api-server=https://api.aps2.sh.basespace.illumina.com/ --access-token=6c48fb1707d6456d8be5b1ae0986f1bd --concurrency=low --exclude '*' --include '**1001660021819-E-G12*' .
$HOME/bs upload dataset --project=6645639 --api-server=https://api.aps2.sh.basespace.illumina.com/ --access-token=6c48fb1707d6456d8be5b1ae0986f1bd --concurrency=low --exclude '*' --include '**1001660021819-E-H*' .

# geomx_281023_PAB: upload Y, Run Y
$HOME/bs create project --name=geomx_281023_PAB --exist-ok --api-server=https://api.aps2.sh.basespace.illumina.com/ --access-token=6c48fb1707d6456d8be5b1ae0986f1bd
cd ~/USyd/spatial_transcriptomics/analysis/geomx/geomx_oct2023/merged_reps/geomx_281023_PAB
#$HOME/bs upload dataset --project=6644640 --api-server=https://api.aps2.sh.basespace.illumina.com/ --access-token=6c48fb1707d6456d8be5b1ae0986f1bd --concurrency=low --exclude '*' --include '**1001660021703*' .
#$HOME/bs upload dataset --project=6644640 --api-server=https://api.aps2.sh.basespace.illumina.com/ --access-token=6c48fb1707d6456d8be5b1ae0986f1bd --concurrency=low --exclude '*' --include '**1001660021704*' .
$HOME/bs upload dataset --project=6644640 --api-server=https://api.aps2.sh.basespace.illumina.com/ --access-token=6c48fb1707d6456d8be5b1ae0986f1bd --concurrency=low --exclude '*' --include '**1001660021704-B-E*' .
$HOME/bs upload dataset --project=6644640 --api-server=https://api.aps2.sh.basespace.illumina.com/ --access-token=6c48fb1707d6456d8be5b1ae0986f1bd --concurrency=low --exclude '*' --include '**1001660021704-B-F*' .
$HOME/bs upload dataset --project=6644640 --api-server=https://api.aps2.sh.basespace.illumina.com/ --access-token=6c48fb1707d6456d8be5b1ae0986f1bd --concurrency=low --exclude '*' --include '**1001660021704-B-G*' .
$HOME/bs upload dataset --project=6644640 --api-server=https://api.aps2.sh.basespace.illumina.com/ --access-token=6c48fb1707d6456d8be5b1ae0986f1bd --concurrency=low --exclude '*' --include '**1001660021704-B-H*' .

# geomx_281023_PH: upload Y, Run Y
$HOME/bs create project --name=geomx_281023_PH --exist-ok --api-server=https://api.aps2.sh.basespace.illumina.com/ --access-token=6c48fb1707d6456d8be5b1ae0986f1bd
cd ~/USyd/spatial_transcriptomics/analysis/geomx/geomx_oct2023/merged_reps/geomx_281023_PH
$HOME/bs upload dataset --project=6657651 --api-server=https://api.aps2.sh.basespace.illumina.com/ --access-token=6c48fb1707d6456d8be5b1ae0986f1bd --concurrency=low --exclude '*' --include '**1001660021822*' .

# check projects and datasets
$HOME/bs list projects --api-server=https://api.aps2.sh.basespace.illumina.com/ --access-token=6c48fb1707d6456d8be5b1ae0986f1bd
$HOME/bs list dataset --api-server=https://api.aps2.sh.basespace.illumina.com/ --access-token=6c48fb1707d6456d8be5b1ae0986f1bd | grep "DSP-1001660021812-D-A06" | wc -l

## NOTE: upload .ini files online at https://aps2.sh.basespace.illumina.com/projects

# list workflows
$HOME/bs list applications --category=workflow
#| GeoMx® NGS Pipeline                            | 2163161 | 2.0.21 

# biosample/project/dataset IDs
#$HOME/bs list datasets --project-name=CHA12467
$HOME/bs list datasets --project=6402396


# launch application
$HOME/bs launch application -n 'GeoMx® NGS Pipeline' -i '2163161' --app-version 2.0.21 \
-o project-id:1232 -o sample-id:2323244

# record analysed projects
# geomx_281023_PB
# geomx_281023_PF
# geomx_281023_PEG
# geomx_281023_PAB
# geomx_281023_PH
# geomx_281023_round1


DSP-1001660021814-B-G10.dcc - error processing fastq - geomx_281023_PB
DSP-1001660021820-F-B10.dcc - error processing fastq - geomx_281023_PF

re-upload reads and re-run


# may need to combine technical replicates
cd /Users/zacc/USyd/spatial_transcriptomics/analysis/geomx/geomx_oct2023
mkdir -p merged_reps

# combine technical replicates
for i in $(awk 'NR>1 {print $1}' Annotations_remapped.txt | grep "DSP-1001660021820-F-B10" ); do
  echo $i
  du -sh */${i}*_R1_*.fastq.gz
  du -sh */${i}*_R2_*.fastq.gz

  # cat files
  cat */${i}*_R1_*.fastq.gz > merged_reps/${i}_S0_L001_R1_001.fastq.gz
  cat */${i}*_R2_*.fastq.gz > merged_reps/${i}_S0_L001_R2_001.fastq.gz
done

# geomx_281023_PEG: upload Y, Run Y
$HOME/bs create project --name=geomx_281023_rerun --exist-ok --api-server=https://api.aps2.sh.basespace.illumina.com/ --access-token=6c48fb1707d6456d8be5b1ae0986f1bd
$HOME/bs upload dataset --project=6698693 --api-server=https://api.aps2.sh.basespace.illumina.com/ --access-token=6c48fb1707d6456d8be5b1ae0986f1bd --concurrency=low --exclude '*' --include '*' .







