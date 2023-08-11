# Commands for Illumina BaseSpace Sequence Hub command line interface (CLI)

cd ~

# authenticate account
./bs auth --api-server https://api.aps2.sh.basespace.illumina.com

# get config details
./bs load config 
#export BASESPACE_API_SERVER="https://api.aps2.sh.basespace.illumina.com"
#export BASESPACE_ACCESS_TOKEN="2ba3b51fca60457d878b974e49661cfd"

# find details of whoami
./bs whoami

# create project
./bs create project --name=CHA12467 --exist-ok --api-server=https://api.aps2.sh.basespace.illumina.com/ --access-token=2ba3b51fca60457d878b974e49661cfd

# see projects
./bs list projects --api-server=https://api.aps2.sh.basespace.illumina.com/ --access-token=2ba3b51fca60457d878b974e49661cfd

# upload datasets
./bs upload dataset --project=6206201 --api-server=https://api.aps2.sh.basespace.illumina.com/ --access-token=2ba3b51fca60457d878b974e49661cfd --concurrency=high --exclude '*' --include '*.fastq.gz' /Users/zacc/CHA12467

# list workflows
./bs list applications --category=workflow
#| GeoMx® NGS Pipeline                            | 2163161 | 2.0.21 

# biosample/project/dataset IDs
./bs list datasets --project-name=CHA12467

# launch application
./bs launch application -n 'GeoMx® NGS Pipeline' -i '2163161' --app-version 2.0.21 \
-o project-id:1232 -o sample-id:2323244