#pysradb srp-to-srr --detailed --desc --expand --saveto ${SRP}.tsv ${SRP}

rawdata_dir="$(pwd)/rawdata/"
SRPfile=$(find $rawdata_dir -maxdepth 1 -name "meta.csv")
IFS=','
threads=1
ntask_per_run=300


###### fifo ######
tempfifo=$$.fifo
trap "exec 1000>&-;exec 1000<&-;exit 0" 2
mkfifo $tempfifo;exec 1000<>$tempfifo;rm -rf $tempfifo
for ((i=1; i<=$ntask_per_run; i++));do
    echo >&1000
done

for file in ${SRPfile[@]};
do
  while IFS=$IFS read line
  do
    read -u1000
    {
    srp=$(awk -F "$IFS" '{print $1}' <<<"$line")
    srr=$(awk -F "$IFS" '{print $2}' <<<"$line")
    
    if [[ $srr != "run_accession" ]];then
      {
      echo $srp $srr
      
      while [[ ! -e $rawdata_dir/$srp/$srr/$srr.sra ]] || [[ -e $rawdata_dir/$srp/$srr/$srr.sra.tmp ]] || [[ -e $rawdata_dir/$srp/$srr/$srr.sra.lock ]]
      do
        echo "prefetch $srp/$srr"
        prefetch --output-directory ${srp} --max-size 1000000000000 $srr
      done
      
      if [[ -e $rawdata_dir/$srp/$srr/$srr.sra ]] && [[ ! -e $rawdata_dir/$srp/$srr/$srr.sra.tmp ]] && [[ ! -e $rawdata_dir/$srp/$srr/$srr.sra.lock ]];then
        
        cd $rawdata_dir/$srp/$srr
        
        if [[ ! -f $rawdata_dir/$srp/$srr/fasterq_dump.log ]] || [[ $(grep -i "error" $rawdata_dir/$srp/$srr/fasterq_dump.log) ]];then
          rm -rf ./fasterq.tmp*
          echo "fasterq-dump $srp/$srr"
          fasterq-dump -f --threads $threads --split-3 ${srr}.sra -o $srr 2>$rawdata_dir/$srp/$srr/fasterq_dump.log
          if [ -e ${srr} ]; then
            mv ${srr} ${srr}.fastq
          fi
        fi

        if [[ ! -f $rawdata_dir/$srp/$srr/pigz.log ]];then
          echo "pigz $srp/$srr"
          ls ./ | grep -E "(*.fastq$)|(*.fq$)" | xargs -i pigz -f --processes $threads {}
          echo -e "pigz finished">$rawdata_dir/$srp/$srr/pigz.log
        fi

        if [ -e ${srr}_2.fastq.gz ] ;then 
          echo "$srp/$srr pair-end"

          if [[ ! -f $rawdata_dir/$srp/$srr/reformat_vpair.log ]];then
            reformat.sh in1=${srr}_1.fastq.gz in2=${srr}_2.fastq.gz vpair allowidenticalnames=t 2>$rawdata_dir/$srp/$srr/reformat_vpair.log
            if [[ $? -ne 0 ]];then
              fq1_nlines=$(zcat ${srr}_1.fastq.gz |wc -l)
              fq2_nlines=$(zcat ${srr}_2.fastq.gz |wc -l)
              if [[ $fq1_nlines == $fq2_nlines ]];then
                echo -e "fq1_nlines:$fq1_nlines\nfq2_nlines:$fq2_nlines\nNames appear to be correctly paired." >>$rawdata_dir/$srp/$srr/reformat_vpair.log
              else
                echo -e "ERROR! R1 and R2 for $srp/$srr have different numbers of reads."
                echo -e "ERROR! R1 and R2 for $srp/$srr have different numbers of reads." >>$rawdata_dir/$srp/$srr/reformat_vpair.log
                echo -e "fq1_nlines:$fq1_nlines\nfq2_nlines:$fq2_nlines\n" >>$rawdata_dir/$srp/$srr/reformat_vpair.log
                continue
              fi
            fi
          elif [[ ! $(grep "Names appear to be correctly paired" $rawdata_dir/$srp/$srr/reformat_vpair.log) ]];then
            fq1_nlines=$(zcat ${srr}_1.fastq.gz |wc -l)
            fq2_nlines=$(zcat ${srr}_2.fastq.gz |wc -l)
            if [[ $fq1_nlines == $fq2_nlines ]];then
              echo -e "fq1_nlines:$fq1_nlines\nfq2_nlines:$fq2_nlines\nNames appear to be correctly paired." >>$rawdata_dir/$srp/$srr/reformat_vpair.log
            else
              echo -e "ERROR! R1 and R2 for $srp/$srr have different numbers of reads."
              echo -e "ERROR! R1 and R2 for $srp/$srr have different numbers of reads." >>$rawdata_dir/$srp/$srr/reformat_vpair.log
              echo -e "fq1_nlines:$fq1_nlines\nfq2_nlines:$fq2_nlines\n" >>$rawdata_dir/$srp/$srr/reformat_vpair.log
              continue
            fi
          elif [[ $(grep "ERROR! R1 and R2 for $srp/$srr have different numbers of reads." $rawdata_dir/$srp/$srr/reformat_vpair.log) ]];then
            echo -e "ERROR! R1 and R2 for $srp/$srr have different numbers of reads."
            continue  
          fi

        else
          echo "$srp/$srr single-end"
        fi
        
      fi
      
      }&
    fi

    echo >&1000
    }&

  done < $file
  wait

done

wait 
echo "done"
