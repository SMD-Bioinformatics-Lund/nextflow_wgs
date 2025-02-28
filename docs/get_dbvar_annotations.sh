wget https://ftp.ncbi.nlm.nih.gov/pub/dbVar/sandbox/sv_datasets/nonredundant/deletions/GRCh38.nr_deletions.pathogenic.tsv.gz
wget https://ftp.ncbi.nlm.nih.gov/pub/dbVar/sandbox/sv_datasets/nonredundant/duplications/GRCh38.nr_duplications.pathogenic.tsv.gz
wget https://ftp.ncbi.nlm.nih.gov/pub/dbVar/sandbox/sv_datasets/nonredundant/insertions/GRCh38.nr_insertions.pathogenic.tsv.gz

zgrep -v ^# GRCh38.nr_deletions.pathogenic.tsv.gz | awk -F'\t' -v OFS='\t' '{split($12, arr, ";"); $13=arr[2]; $12=arr[1]; print $1, $2, $3, $10, $11, $12, $13}' > dbvar_del.bed
zgrep -v ^# GRCh38.nr_duplications.pathogenic.tsv.gz | awk -F'\t' -v OFS='\t' '{split($12, arr, ";"); $13=arr[2]; $12=arr[1]; print $1, $2, $3, $10, $11, $12, $13}' > dbvar_dup.bed
zgrep -v ^# GRCh38.nr_insertions.pathogenic.tsv.gz | awk -F'\t' -v OFS='\t' '{split($12, arr, ";"); $13=arr[2]; $12=arr[1]; print $1, $2, $3, $10, $11, $12, $13}' > dbvar_ins.bed
