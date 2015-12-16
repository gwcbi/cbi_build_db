#! /bin/bash

##########################################################################################
# This is an automated script for building NCBI databases
#
##########################################################################################
. dbbuild.sh

### Microbiome groups ####################################################################

micro="archaea bacteria fungi protozoa viral"

### Download assemblies for groups
for g in $micro; do download_group $g; done

### Count assemblies for groups
for g in $micro; do count_group $g; done
    
### Concatenate fasta files for groups
for g in $micro; do sequences_group $g; done

### Eukaryotic parasites #################################################################

### Download the assembly summary file for genbank
wget -O assembly_summary_genbank.txt ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/assembly_summary_genbank.txt

### Create assembly summary with selected assemblies
rm -f selected.assembly_summary.txt

# Tapeworms: class Cestoda - 6199
taxid_asmlist 6199 >> selected.assembly_summary.txt
# Flatworms: class Trematoda - 6178
taxid_asmlist 6178 >> selected.assembly_summary.txt
# Nematodes: class Chromadorea - 119089
taxid_asmlist 119089 >> selected.assembly_summary.txt
# Nematodes: class Enoplea - 119088
taxid_asmlist 119088 >> selected.assembly_summary.txt

### Download selected assemblies
download_asmsum selected.assembly_summary.txt

### Concatenate fasta files for selected assemblies
# Note: All these parasites are in the invertebrate group
python get_seqs_with_ti.py --sumfile assembly_summary_genbank.txt --db genbank invertebrate | gzip > sequences/invertebrate.fna.gz


### Build index ##########################################################################

### Including Multicellular eukaryotes
sbatch -N 1 -p defq -t 5760 <<'EOF'
#! /bin/bash

SN="build_index"

#--- Start the timer
t1=$(date +"%s")

#--- Concatenate all sequence files
module load pigz
rm -f /scratch/allseqs.fna
pigz -dc sequences/archaea.fna.gz >> /scratch/allseqs.fna
pigz -dc sequences/bacteria.fna.gz >> /scratch/allseqs.fna
pigz -dc sequences/fungi.fna.gz >> /scratch/allseqs.fna
pigz -dc sequences/protozoa.fna.gz >> /scratch/allseqs.fna
pigz -dc sequences/viral.fna.gz >> /scratch/allseqs.fna
pigz -dc sequences/invertebrate.fna.gz >> /scratch/allseqs.fna

#---Complete
t2=$(date +"%s")
diff=$(($t2-$t1))
echo "[---$SN---] ($(date)) $(($diff / 60)) minutes and $(($diff % 60)) seconds elapsed."

#--- Start the timer
t1=$(date +"%s")

#--- Build bowtie index
module load bowtie2

tdate=$(date "+%Y%m%d")
destdir="mbdb.$tdate"
mkdir -p $destdir
bowtie2-build --large-index --offrate 3 /scratch/allseqs.fna $destdir/index

#---Complete
t2=$(date +"%s")
diff=$(($t2-$t1))
echo "[---$SN---] ($(date)) $(($diff / 60)) minutes and $(($diff % 60)) seconds elapsed."

echo "[---$SN---] ($(date)) $SN COMPLETE."

EOF

### Excluding Multicellular eukaryotes
sbatch -N 1 -p defq -t 5760 <<'EOF'
#! /bin/bash

SN="build_index"

#--- Start the timer
t1=$(date +"%s")

#--- Concatenate all sequence files
module load pigz
rm -f /scratch/allseqs.fna
pigz -dc sequences/archaea.fna.gz >> /scratch/allseqs.fna
pigz -dc sequences/bacteria.fna.gz >> /scratch/allseqs.fna
pigz -dc sequences/fungi.fna.gz >> /scratch/allseqs.fna
pigz -dc sequences/protozoa.fna.gz >> /scratch/allseqs.fna
pigz -dc sequences/viral.fna.gz >> /scratch/allseqs.fna

#---Complete
t2=$(date +"%s")
diff=$(($t2-$t1))
echo "[---$SN---] ($(date)) $(($diff / 60)) minutes and $(($diff % 60)) seconds elapsed."

#--- Start the timer
t1=$(date +"%s")

#--- Build bowtie index
module load bowtie2

tdate=$(date "+%Y%m%d")
destdir="mbdb_noinv.$tdate"
mkdir -p $destdir
bowtie2-build --large-index --offrate 3 /scratch/allseqs.fna $destdir/index

#---Complete
t2=$(date +"%s")
diff=$(($t2-$t1))
echo "[---$SN---] ($(date)) $(($diff / 60)) minutes and $(($diff % 60)) seconds elapsed."

echo "[---$SN---] ($(date)) $SN COMPLETE."

EOF

