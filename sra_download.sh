#!/usr/bin/env bash
trap 'trap - SIGTERM && kill -- -$$' SIGINT SIGTERM
#pysradb srp-to-srr --detailed --desc --expand --saveto ${SRP}.tsv ${SRP}

rawdata_dir="$(pwd)/rawdata/"
SRPfile="SRPmeta.tsv"
ifs='\t'
threads=3
ntask_per_run=120
force_extract="TRUE"

###### fifo ######
tempfifo=$$.fifo
trap "exec 1000>&-;exec 1000<&-;exit 0" 2
mkfifo $tempfifo
exec 1000<>$tempfifo
rm -rf $tempfifo
for ((i = 1; i <= $ntask_per_run; i++)); do
  echo >&1000
done

if [[ ! -d $rawdata_dir ]]; then
  mkdir -p $rawdata_dir
fi

while IFS=$ifs read line; do

  srp=$(awk -F "$ifs" '{print $1}' <<<"$line")
  srr=$(awk -F "$ifs" '{print $2}' <<<"$line")

  if [[ $srr != "run_accession" ]]; then

    while [[ ! -e $rawdata_dir/$srp/$srr/$srr.sra ]] || [[ -e $rawdata_dir/$srp/$srr/$srr.sra.tmp ]] || [[ -e $rawdata_dir/$srp/$srr/$srr.sra.lock ]]; do
      echo "prefetch $srp/$srr"
      cd $rawdata_dir
      prefetch --output-directory ${srp} --max-size 1000000000000 $srr &
      sleep 60
    done

    read -u1000
    {
      while :; do
        if [[ -e $rawdata_dir/$srp/$srr/$srr.sra ]] && [[ ! -e $rawdata_dir/$srp/$srr/$srr.sra.tmp ]] && [[ ! -e $rawdata_dir/$srp/$srr/$srr.sra.lock ]]; then
          echo $srp/$srr
          cd $rawdata_dir/$srp/$srr

          if [[ $force_extract == "TRUE" ]];then 
            rm -f $rawdata_dir/$srp/$srr/fasterq_dump.log $rawdata_dir/$srp/$srr/pigz.log $rawdata_dir/$srp/$srr/reformat_vpair.log
          fi
          
          if [[ ! -f $rawdata_dir/$srp/$srr/fasterq_dump.log ]] || [[ ! $(grep -i "error" $rawdata_dir/$srp/$srr/fasterq_dump.log) ]] || [[ $force_extract == "TRUE" ]]; then
            rm -rf ./fasterq.tmp*
            echo "fasterq-dump $srp/$srr"
            fasterq-dump -f --threads $threads --split-3 ${srr}.sra -o $srr 2>$rawdata_dir/$srp/$srr/fasterq_dump.log
            if [ -e ${srr} ]; then
              mv ${srr} ${srr}.fastq
            fi
          fi

          if [[ ! -f $rawdata_dir/$srp/$srr/pigz.log ]] || [[ $force_extract == "TRUE" ]]; then
            echo "pigz $srp/$srr"
            ls ./ | grep -E "(*.fastq$)|(*.fq$)" | xargs -i pigz -f --processes $threads {}
            echo -e "pigz finished" >$rawdata_dir/$srp/$srr/pigz.log
          fi

          if [ -e ${srr}_2.fastq.gz ]; then
            echo "$srp/$srr pair-end"

            if [[ ! -f $rawdata_dir/$srp/$srr/reformat_vpair.log ]] || [[ ! $(grep "Names appear to be correctly paired" $rawdata_dir/$srp/$srr/reformat_vpair.log) ]]; then
              reformat.sh in1=${srr}_1.fastq.gz in2=${srr}_2.fastq.gz vpair allowidenticalnames=t 2>$rawdata_dir/$srp/$srr/reformat_vpair.log
            fi

            if [[ ! $(grep "Names appear to be correctly paired" $rawdata_dir/$srp/$srr/reformat_vpair.log) ]]; then
              fq1_nlines=$(zcat ${srr}_1.fastq.gz | wc -l)
              fq2_nlines=$(zcat ${srr}_2.fastq.gz | wc -l)
              if [[ $fq1_nlines == $fq2_nlines ]]; then
                echo -e "fq1_nlines:$fq1_nlines\nfq2_nlines:$fq2_nlines\nNames appear to be correctly paired(custom)" >>$rawdata_dir/$srp/$srr/reformat_vpair.log
              else
                echo -e "fq1_nlines:$fq1_nlines\nfq2_nlines:$fq2_nlines\n" >>$rawdata_dir/$srp/$srr/reformat_vpair.log
                echo -e "ERROR! R1 and R2 for $srp/$srr have different numbers of reads."
                echo -e "ERROR! R1 and R2 for $srp/$srr have different numbers of reads." >>$rawdata_dir/$srp/$srr/reformat_vpair.log
                break
              fi
            fi

          else
            echo "$srp/$srr single-end"
          fi

          break

        else
          sleep 60
        fi

      done
      echo >&1000
    } &

  fi

done <"$SRPfile"

wait
echo "done"
