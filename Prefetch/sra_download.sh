#!/usr/bin/env bash
trap "trap - SIGTERM && kill -- -$$" SIGINT SIGTERM EXIT
#trap "exit" SIGINT SIGTERM
#trap "kill 0" EXIT

#### User can prepare the SRP meta file using the command 'pysradb srp-to-srr --detailed --desc --expand --saveto ${SRP}.tsv ${SRP}'
#### SRP meta file must include the fields "study_accession" "run_accession" "experiment_accession" "sample_accession" "run_total_spots"
#### Dependencies:
#### sra-tools (https://github.com/ncbi/sra-tools)
#### pigz (https://zlib.net/pigz)

#############################  Paramaters #############################
rawdata_dir="$(pwd)/rawdata/"
SRPfile="$(pwd)/SRP_meta_file.tsv"
ifs='\t'
threads=4
ntask_per_run=90
force_process="FALSE"

########################################################################

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

if [[ ! -f $SRPfile ]]; then
  echo "ERROR! No such file: $SRPfile"
  exit 1
fi

var_extract=$(awk -F "$ifs" '
  BEGIN {OFS=FS}
  NR==1 {
      for (i=1; i<=NF; i++) {
          f[$i] = i
      }
  }
  { print $(f["study_accession"]), $(f["run_accession"]), $(f["experiment_accession"]), $(f["sample_accession"]), $(f["run_total_spots"]) }
  ' "$SRPfile")

while IFS=$ifs read line; do
  srp=$(awk -F "$ifs" '{print $1}' <<<"$line")
  srr=$(awk -F "$ifs" '{print $2}' <<<"$line")
  if [[ "$srr" =~ SRR* ]]; then
    {
      while [[ ! -e $rawdata_dir/$srp/$srr/$srr.sra ]] || [[ -e $rawdata_dir/$srp/$srr/$srr.sra.tmp ]] || [[ -e $rawdata_dir/$srp/$srr/$srr.sra.lock ]]; do
        echo -e "+++++ $srp/$srr: Prefetching SRR +++++"
        cd $rawdata_dir
        prefetch --output-directory ${srp} --max-size 1000000000000 ${srr}
        sleep 60
      done
    } &
  fi
done <<<"$var_extract"

logfile="$rawdata_dir/process.log"
if [[ ! -f $logfile ]];then
  touch $logfile
fi

line_count=0
total_count=$(cat "$SRPfile" | wc -l)

while IFS=$ifs read line; do
  ((line_count++))
  read -u1000
  {
    srp=$(awk -F "$ifs" '{print $1}' <<<"$line")
    srr=$(awk -F "$ifs" '{print $2}' <<<"$line")
    srx=$(awk -F "$ifs" '{print $3}' <<<"$line")
    srs=$(awk -F "$ifs" '{print $4}' <<<"$line")
    nreads=$(awk -F "$ifs" '{print $5}' <<<"$line")
    nreads=$(echo "$nreads" | xargs)

    echo -e "########### $line_count/$total_count ###########"
    echo -e "RECORD INFO: $srp/$srr/$srx/$srs/$nreads"

    if [[ "$srr" =~ SRR* ]]; then

      force=${force_process}
      status="uncompleted"
      attempt=0

      while [[ $status == "uncompleted" ]] && (("$attempt" <= 1)); do
        if [[ -e $rawdata_dir/$srp/$srr/$srr.sra ]] && [[ ! -e $rawdata_dir/$srp/$srr/$srr.sra.tmp ]] && [[ ! -e $rawdata_dir/$srp/$srr/$srr.sra.lock ]]; then
          ((attempt++))
          if [[ $attempt != 1 ]]; then
            echo -e "+++++ $srp/$srr: Number of attempts: $attempt +++++"
          fi

          echo -e "+++++ $srp/$srr: Processing sra file +++++"
          cd $rawdata_dir/$srp/$srr

          if [[ $force == "TRUE" ]]; then
            rm -f $rawdata_dir/$srp/$srr/fasterq_dump.log $rawdata_dir/$srp/$srr/fasterq_dump_process.log $rawdata_dir/$srp/$srr/pigz.log
          fi

          rm -rf $rawdata_dir/$srp/$srr/fasterq.tmp*
          if [[ ! -f $rawdata_dir/$srp/$srr/fasterq_dump.log ]]; then
            fasterq-dump -f --threads $threads --split-3 ${srr}.sra -o $srr 2>$rawdata_dir/$srp/$srr/fasterq_dump_process.log
            if [ -e ${srr} ]; then
              mv ${srr} ${srr}.fastq
            fi
            echo -e "fasterq-dump finished" >$rawdata_dir/$srp/$srr/fasterq_dump.log
            echo -e "+++++ $srp/$srr: fasterq-dump done +++++"
          else
            echo -e "+++++ $srp/$srr: fasterq-dump skipped +++++"
          fi

          if ([[ ! $(ls $rawdata_dir/$srp/$srr/ | grep -E "(*.fastq.gz$)") ]] || [[ ! -f $rawdata_dir/$srp/$srr/pigz.log ]]) && [[ $(ls $rawdata_dir/$srp/$srr/ | grep -E "(*.fastq$)") ]]; then
            ls $rawdata_dir/$srp/$srr/ | grep -E "(*.fastq$)|(*.fq$)" | xargs -i pigz -f --processes $threads {}
            echo -e "pigz finished" >$rawdata_dir/$srp/$srr/pigz.log
            echo -e "+++++ $srp/$srr: pigz done +++++"
          fi

          if [[ -f ${srr}_1.fastq.gz ]] && [[ -f ${srr}_2.fastq.gz ]]; then
            pigz -t ${srr}_1.fastq.gz 2>/dev/null
            if [[ $? != 0 ]]; then
              echo -e "Warning! $srp/$srr: ${srr}_1.fastq.gz is not a completed .gz file."
              force="TRUE"
              continue
            fi
            pigz -t ${srr}_2.fastq.gz 2>/dev/null
            if [[ $? != 0 ]]; then
              echo -e "Warning! $srp/$srr: ${srr}_2.fastq.gz is not a completed .gz file."
              force="TRUE"
              continue
            fi
            echo -e "+++++ $srp/$srr: Integrity check passed +++++"
          elif [[ -f ${srr}.fastq.gz ]]; then
            pigz -t ${srr}.fastq.gz 2>/dev/null
            if [[ $? != 0 ]]; then
              echo -e "Warning! $srp/$srr: ${srr}.fastq.gz is not a completed .gz file."
              force="TRUE"
              continue
            fi
            echo -e "+++++ $srp/$srr: Integrity check passed +++++"
          else
            echo -e "Warning! $srp/$srr: Cannot find the correct fastq.gz file name."
            force="TRUE"
            continue
          fi

          if [[ -f ${srr}_1.fastq.gz ]] && [[ -f ${srr}_2.fastq.gz ]]; then

            fq1_nlines=$(unpigz -c "${srr}"_1.fastq.gz | wc -l)
            fq2_nlines=$(unpigz -c "${srr}"_2.fastq.gz | wc -l)
            echo -e "fq1_nlines:$fq1_nlines     fq1_nreads:$((fq1_nlines / 4))\nfq2_nlines:$fq2_nlines    fq2_nreads:$((fq2_nlines / 4))\n" >$rawdata_dir/$srp/$srr/fqcheck.log
            if [[ $fq1_nlines == "$fq2_nlines" ]]; then
              if [[ ! "$nreads" =~ ^[0-9]+$ ]]; then
                if [[ $((fq1_nlines % 4)) == 0 ]]; then
                  status="completed"
                  echo -e "+++++ $srp/$srr: Success! Processing completed. +++++"
                else
                  force="TRUE"
                  echo -e "Warning! $srp/$srr: Line count is not divisible by 4."
                fi
              else
                if [[ $fq1_nlines == $((nreads * 4)) ]]; then
                  status="completed"
                  echo -e "+++++ $srp/$srr: Success! Processing completed. +++++"
                else
                  force="TRUE"
                  echo -e "Warning! $srp/$srr has different numbers of reads with that SRP meta file recorded:\n         fq1=$((fq1_nlines / 4)) / Recorded=$nreads"
                fi
              fi

            else
              force="TRUE"
              echo -e "Warning! $srp/$srr has different numbers of reads between paired files:\n         fq1=$((fq1_nlines / 4))/ fq2=$((fq2_nlines / 4))"
            fi

          elif [[ -f ${srr}.fastq.gz ]]; then
            fq1_nlines=$(unpigz -c "${srr}".fastq.gz | wc -l)
            echo -e "fq1_nlines:$fq1_nlines   fq1_nreads:$((fq1_nlines / 4))\n" >$rawdata_dir/$srp/$srr/fqcheck.log

            if [[ ! "$nreads" =~ ^[0-9]+$ ]]; then
              if [[ $((fq1_nlines % 4)) == 0 ]]; then
                status="completed"
                echo -e "+++++ $srp/$srr: Success! Processing completed +++++"
              else
                force="TRUE"
                echo -e "Warning! $srp/$srr: Line count is not divisible by 4."
              fi
            else
              if [[ $fq1_nlines == $((nreads * 4)) ]]; then
                status="completed"
                echo -e "+++++ $srp/$srr: Success! Processing completed +++++"
              else
                force="TRUE"
                echo -e "Warning! $srp/$srr has different numbers of lines with that SRP meta file recorded:\n         fq1=$((fq1_nlines / 4)) / Recorded=$nreads"
              fi
            fi

          else
            force="TRUE"
            echo -e "Warning! Can not find any fastq.gz file! $srp/$srr has to restart the processing."
          fi

        else
          sleep 300
        fi
      done

      if [[ $status == "completed" ]]; then
        echo "Completed: $srp/$srr" >>"$logfile"
      else
        echo "Interrupted: $srp/$srr" >>"$logfile"  
        echo -e "\033[31mERROR! $srp/$srr interrupted. Please check the number of reads or re-download the SRA file.\033[0m"
      fi

    fi

    echo "\033[32m***** Completed:$(cat "$logfile" | grep "Completed" | uniq | wc -l) | Interrupted:$(cat "$logfile" | grep "Interrupted" | uniq | wc -l) | Total:$total_count *****\033[0m"

    echo >&1000
  } &
done <<<"$var_extract"
wait

ELAPSED="Elapsed: $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min $(($SECONDS % 60))sec"
echo -e "\n$ELAPSED"
echo -e "****************** SRA Prefetching Finished  ******************\n"
