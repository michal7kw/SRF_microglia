cd /beegfs/scratch/ric.cosr/moscato.giosue/Alberti_testseq
prj_dir=$(pwd | cut -d'/' -f6)
mkdir -p ./fastq/$prj_dir


# Definisci la directory contenente i file
cd /beegfs/scratch/ric.cosr/moscato.giosue/Alberti_testseq/fastq

# Loop attraverso tutti i file nella directory
for file in *.fastq.gz
do
    # Controlla se il file è un file regolare
    if [ -f "$file" ]; then
        # Rinomina il file sostituendo tutti i trattini con underscore
        nuovo_nome=$(echo "$file" | sed 's/-/_/g')
        mv "$file" "$nuovo_nome"
        echo "Il file $file è stato rinominato in $nuovo_nome"
    fi
done





for file in *_R1_001.fastq.gz
do
    echo $file
    folder=$(basename $file _R1_001.fastq.gz)
    echo $folder
    mkdir -p ./$prj_dir/$folder
    mv "$file" "./$prj_dir/$folder/${file%_R1_001.fastq.gz}_L001_R1_001.fastq.gz"

done

for file in *_R2_001.fastq.gz
do
    echo $file
    folder=$(basename $file _R2_001.fastq.gz)
    echo $folder
    #mkdir -p ./fastq/$prj_dir/$folder
    mv "$file" "./$prj_dir/$folder/${file%_R2_001.fastq.gz}_L001_R2_001.fastq.gz"

done


