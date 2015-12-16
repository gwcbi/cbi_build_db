#! /bin/bash

### Download the latest refseq assemblies from group
download_group () {
    ### Set tgroup
    local tgroup=$1

    ### Download the assembly summary file for taxonomic group
    wget -O ${tgroup}.assembly_summary.txt ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/${tgroup}/assembly_summary.txt

    ### Calculate number of "latest" assemblies
    local expnum=$(awk -F "\t" '$11=="latest"{print $20}' $tgroup.assembly_summary.txt | wc -l)
    echo "Expecting $expnum assemblies"

    ### Find unique two-letter prefixes for all species
    rsync --list-only ftp.ncbi.nlm.nih.gov::genomes/refseq/${tgroup}/ | \
        perl -ne '/^d[rwxs-]{9}\s+\d+\s+[\d\/]{10}\s+\d+:\d+:\d+\s+(?!\.)(\S+)\n$/ && print "$1\n"' | \
        cut -c-2 | sort | uniq > prefixes.txt
    
    local dlnum=0
    local counter=0
    while [ $dlnum -lt $expnum ] && [ $counter -lt 5 ]; do
        ### Download using rsync (one-shot method)
        # rsync -avR ftp.ncbi.nlm.nih.gov::genomes/refseq/${tgroup}/*/latest_assembly_versions/*/*_genomic.{fna,gff,gbff}.gz ./

        ### Download using rsync -- each prefix
        cat prefixes.txt | while read px; do
            echo "[---------- Prefix $px ----------]"
            rsync -avR ftp.ncbi.nlm.nih.gov::genomes/refseq/${tgroup}/${px}*/latest_assembly_versions/*/*_genomic.{fna,gff,gbff}.gz ./
        done

        ### Number of downloaded assemblies
        dlnum=$(ls refseq/$tgroup/*/latest_assembly_versions/*/*_genomic.fna.gz | wc -l)
        echo "Expected $expnum, downloaded $dlnum"
        counter=$((counter+1))
    done

    if [ $dlnum != $expnum ] && [ $counter == 5 ]; then
        echo "WARNING: Counter maxed out. Continuing"
    fi
    
    rm -f prefixes.txt
}

### Count the number of assemblies in group
count_group () {
    ### Set tgroup
    local tgroup=$1
    ### Calculate number of "latest" assemblies
    local expnum=$(awk -F "\t" '$11=="latest"{print $20}' $tgroup.assembly_summary.txt | wc -l)
    ### Number of downloaded assemblies
    dlnum=$(ls refseq/$tgroup/*/latest_assembly_versions/*/*_genomic.fna.gz | wc -l)
    echo -e "$tgroup\t$expnum\t$dlnum"      
}

### Make merged sequence files for group
sequences_group() {
    ### Set tgroup
    local tgroup=$1
    
    ### Download the assembly summary file for taxonomic group
    if [[ ! -e ${tgroup}.assembly_summary.txt ]]; then
        wget -O ${tgroup}.assembly_summary.txt ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/${tgroup}/assembly_summary.txt
    fi

    ### Get the sequences
    mkdir -p sequences
    python get_seqs_with_ti.py ${tgroup} | gzip > sequences/${tgroup}.fna.gz
}

### Search assemblies for a taxonomy ID
taxid_asmlist() {
    local txid=$1
    ### Find assembly IDs matching taxonomy ID
    curl -s "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=assembly&term=txid${txid}\[Organism:exp\]" | \
        grep '^<Id>' | sed 's|</*Id>||g' | \
        while read aid; do
            curl -s "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=assembly&id=$aid" | \
            grep '<AssemblyAccession' | \
            perl -pe 's|.*<AssemblyAccession.*>(.*)</AssemblyAccession>.*|\1|' | \
            grep -f - assembly_summary_genbank.txt
        done
}

### Download assemblies from an assembly list
download_asmsum() {
    local asmsum=$1
    expnum=$(awk -F "\t" '$11=="latest"{print $20}' $asmsum | wc -l)
    echo "Expecting $expnum assemblies"

    local dlnum=0
    local counter=0
    while [ $dlnum -lt $expnum ] && [ $counter -lt 5 ]; do
        ### Download using rsync
        awk -F "\t" '$11=="latest"{print $8}' $asmsum | sed 's/ /_/g' | while read spname; do
            rsync -avR ftp.ncbi.nlm.nih.gov::genomes/genbank/*/$spname/latest_assembly_versions/*/*_genomic.{fna,gff,gbff}.gz ./
        done
        ### Number of downloaded assemblies
        dlnum=$(awk -F "\t" '$11=="latest"{print $8}' $asmsum | sed 's/ /_/g' | while read spname; do
                    echo $(ls genbank/*/$spname/latest_assembly_versions/*/*_genomic.fna.gz | wc -l)
                done | perl -nle '$sum += $_ } END { print $sum')
        echo "Expected $expnum, downloaded $dlnum"
        counter=$((counter+1))
    done

    if [ $dlnum != $expnum ] && [ $counter == 5 ]; then
        echo "WARNING: Counter maxed out. Continuing"
    fi
}

