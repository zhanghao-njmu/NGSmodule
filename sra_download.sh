#pysradb srp-to-srr --detailed --desc --expand --saveto ${SRP}.tsv ${SRP}

rawdata_dir=/data/database/SRR_collection/human/germ_cell_specification/rawdata/
SRPfile=$(find $rawdata_dir -maxdepth 1 -name "meta2.csv")
IFS=','
threads=1

for file in ${SRPfile[@]};
do
  while IFS=$IFS read line
  do
    srp=$(awk -F "$IFS" '{print $1}' <<<"$line")
    srr=$(awk -F "$IFS" '{print $2}' <<<"$line")
    
    if [[ $srr != "run_accession" ]];then
      {
      echo $srp $srr
      
      while [[ ! -e $rawdata_dir/$srp/$srr/$srr.sra ]] || [[ -e $rawdata_dir/$srp/$srr/$srr.sra.tmp ]] || [[ -e $rawdata_dir/$srp/$srr/$srr.sra.lock ]]
      do
        prefetch --output-directory ${srp} --max-size 1000000000000 $srr
      done
      
      if [[ -e $rawdata_dir/$srp/$srr/$srr.sra ]] && [[ ! -e $rawdata_dir/$srp/$srr/$srr.sra.tmp ]] && [[ ! -e $rawdata_dir/$srp/$srr/$srr.sra.lock ]];then
        cd $rawdata_dir/$srp/$srr
        rm -rf ./fasterq.tmp*
        fasterq-dump -f --threads $threads --split-3 ${srr}.sra -o $srr
        if [ -e ${srr} ]; then
          mv ${srr} ${srr}.fastq
        fi
        
        ls ./ | grep -E "(*.fastq$)|(*.fq$)" | xargs -i pigz -f --processes $threads {}
        
        if [ -e ${srr}_2.fastq.gz ] ;then 
          echo "$srp/$srr  pair-end"
          reformat.sh in1=${srr}_1.fastq.gz in2=${srr}_2.fastq.gz vpair >/dev/null 2>&1 
          [ $? -ne 0 ] && { echo -e "ERROR:${srr}_1.fastq.gz and ${srr}_2.fastq.gz have different numbers of reads or have different read names!\n"; continue; } 
        else
          echo "$srp/$srr  single-end"
        fi
      fi
      
      }&
    fi
  done < $file
done

wait 
echo "done"
