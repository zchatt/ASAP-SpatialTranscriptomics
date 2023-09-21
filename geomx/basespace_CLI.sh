# Commands for Illumina BaseSpace Sequence Hub command line interface (CLI)

cd ~

# authenticate account
$HOME/bs auth --api-server https://api.aps2.sh.basespace.illumina.com

# get config details
$HOME/bs load config 
#export BASESPACE_API_SERVER="https://api.aps2.sh.basespace.illumina.com"
#export BASESPACE_ACCESS_TOKEN="2ba3b51fca60457d878b974e49661cfd"

# find details of whoami
$HOME/bs whoami

# create project
$HOME/bs create project --name=CHA12467 --exist-ok --api-server=https://api.aps2.sh.basespace.illumina.com/ --access-token=2ba3b51fca60457d878b974e49661cfd

# see projects
$HOME/bs list projects --api-server=https://api.aps2.sh.basespace.illumina.com/ --access-token=2ba3b51fca60457d878b974e49661cfd

# # may need to rename files eg.
# cd /Users/zacc/CHA12602r2
# rename 's/_L(\d+)_R(\d+)_(\d+)\.fastq\.gz/|L$1|R$2|$3.fastq.gz/' *_L*_R*_*.fastq.gz
# rename 's/_/-/g' *.fastq.gz
# rename 's/\|/_/g' *.fastq.gz
# rename 's/-S/_S/g' *.fastq.gz

# upload datasets
$HOME/bs upload dataset --project=6206201 --api-server=https://api.aps2.sh.basespace.illumina.com/ --access-token=2ba3b51fca60457d878b974e49661cfd --concurrency=high --exclude '*' --include '*.fastq.gz' /Users/zacc/CHA12467

# list workflows
$HOME/bs list applications --category=workflow
#| GeoMx® NGS Pipeline                            | 2163161 | 2.0.21 

# biosample/project/dataset IDs
#$HOME/bs list datasets --project-name=CHA12467
$HOME/bs list datasets --project=6402396

# launch application
$HOME/bs launch application -n 'GeoMx® NGS Pipeline' -i '2163161' --app-version 2.0.21 \
-o project-id:1232 -o sample-id:2323244