#!/usr/bin/env bash
trap 'trap - SIGTERM && kill -- -$$' SIGINT SIGTERM
#### User can prepare the srp meta file using the command 'pysradb srp-to-srr --detailed --desc --expand --saveto ${SRP}.tsv ${SRP}'
#### Requirement:
#### sra-tools (https://github.com/ncbi/sra-tools)
#### pigz (https://zlib.net/pigz)
#### reformat (https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/reformat-guide)

#############################  Paramaters #############################
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

if [[ ! -f $SRPfile ]]; then
  echo "ERROR! No such file:$SRPfile"
  exit 1
fi

var_extract=$(awk -F $ifs '
  NR==1 {
      for (i=1; i<=NF; i++) {
          f[$i] = i
      }
  }
  { print $(f["study_accession"]) "\t" $(f["run_accession"]) "\t" $(f["experiment_accession"]) "\t" $(f["sample_accession"]) "\t" $(f["run_total_spots"]) "\t"  }
  ' "$SRPfile")

line_count=1
total_count=$(cat "$SRPfile" | wc -l)
while IFS=$ifs read line; do

  echo -e "########### $line_count/$total_count ###########"
  ((line_count++))

  srp=$(awk -F "$ifs" '{print $1}' <<<"$line")
  srr=$(awk -F "$ifs" '{print $2}' <<<"$line")
  srx=$(awk -F "$ifs" '{print $3}' <<<"$line")
  srs=$(awk -F "$ifs" '{print $4}' <<<"$line")
  nreads=$(awk -F "$ifs" '{print $5}' <<<"$line")

  if [[ "$srr" =~ SRR* ]]; then

    while [[ ! -e $rawdata_dir/$srp/$srr/$srr.sra ]] || [[ -e $rawdata_dir/$srp/$srr/$srr.sra.tmp ]] || [[ -e $rawdata_dir/$srp/$srr/$srr.sra.lock ]]; do
      echo -e "+++++ $srp/$srr: Prefetching SRR +++++"
      cd $rawdata_dir
      prefetch --output-directory ${srp} --max-size 1000000000000 $srr &
      sleep 60
    done

    read -u1000
    {

      force=${force_extract}
      status="uncompleted"

      while [[ $status == "uncompleted" ]]; do
        if [[ -e $rawdata_dir/$srp/$srr/$srr.sra ]] && [[ ! -e $rawdata_dir/$srp/$srr/$srr.sra.tmp ]] && [[ ! -e $rawdata_dir/$srp/$srr/$srr.sra.lock ]]; then
          echo -e "+++++ $srp/$srr: Processing sra file +++++"
          cd $rawdata_dir/$srp/$srr

          if [[ $force == "TRUE" ]]; then
            rm -f $rawdata_dir/$srp/$srr/fasterq_dump.log $rawdata_dir/$srp/$srr/fasterq_dump_process.log $rawdata_dir/$srp/$srr/pigz.log $rawdata_dir/$srp/$srr/reformat_vpair.log
          fi

          if [[ ! -f $rawdata_dir/$srp/$srr/fasterq_dump.log ]] || [[ $(grep -i "error" $rawdata_dir/$srp/$srr/fasterq_dump_process.log) ]]; then
            rm -rf ./fasterq.tmp*
            fasterq-dump -f --threads $threads --split-3 ${srr}.sra -o $srr 2>$rawdata_dir/$srp/$srr/fasterq_dump_process.log
            if [ -e ${srr} ]; then
              mv ${srr} ${srr}.fastq
            fi
            echo -e "fasterq-dump finished" >$rawdata_dir/$srp/$srr/fasterq_dump.log
            echo -e "+++++ $srp/$srr: fasterq-dump done +++++"
          else
            echo -e "+++++ $srp/$srr: fasterq-dump skipped +++++"
          fi

          if [[ ! -f $rawdata_dir/$srp/$srr/pigz.log ]]; then
            ls ./ | grep -E "(*.fastq$)|(*.fq$)" | xargs -i pigz -f --processes $threads {}
            echo -e "pigz finished" >$rawdata_dir/$srp/$srr/pigz.log
            echo -e "+++++ $srp/$srr: pigz done +++++"
          else
            echo -e "+++++ $srp/$srr: pigz skipped +++++"
          fi

          if [[ -f ${srr}_1.fastq.gz ]] && [[ -f ${srr}_2.fastq.gz ]]; then

            if [[ ! -f $rawdata_dir/$srp/$srr/reformat_vpair.log ]] || [[ ! $(grep "Names appear to be correctly paired" $rawdata_dir/$srp/$srr/reformat_vpair.log) ]]; then
              reformat.sh in1=${srr}_1.fastq.gz in2=${srr}_2.fastq.gz vpair allowidenticalnames=t 2>$rawdata_dir/$srp/$srr/reformat_vpair.log
            fi

            fq1_nlines=$(zcat ${srr}_1.fastq.gz | wc -l)

            if [[ ! $(grep "Names appear to be correctly paired" $rawdata_dir/$srp/$srr/reformat_vpair.log) ]]; then
              fq2_nlines=$(zcat ${srr}_2.fastq.gz | wc -l)
              if [[ $fq1_nlines == $(($nreads * 4)) ]] && [[ $fq1_nlines == $fq2_nlines ]]; then
                echo -e "fq1_nlines:$fq1_nlines\nfq2_nlines:$fq2_nlines\nNames appear to be correctly paired(custom)" >>$rawdata_dir/$srp/$srr/reformat_vpair.log
                status="completed"
                echo -e "+++++ $srp/$srr: Processing completed +++++"
              else
                echo -e "fq1_nlines:$fq1_nlines\nfq2_nlines:$fq2_nlines\n" >>$rawdata_dir/$srp/$srr/reformat_vpair.log
                echo -e "ERROR! R1 and R2 for $srp/$srr have different numbers of reads." >>$rawdata_dir/$srp/$srr/reformat_vpair.log
                force="TRUE"
                status="uncompleted"
                echo -e "ERROR! $srp/$srr has to restart the processing."
              fi
            else
              if [[ $fq1_nlines == $(($nreads * 4)) ]]; then
                status="completed"
                echo -e "+++++ $srp/$srr: Processing completed +++++"
              else
                force="TRUE"
                status="uncompleted"
                echo -e "ERROR! $srp/$srr has to restart the processing."
              fi
            fi

          elif [[ -f ${srr}.fastq.gz ]]; then
            fq1_nlines=$(zcat ${srr}.fastq.gz | wc -l)
            if [[ $fq1_nlines == $(($nreads * 4)) ]]; then
              status="completed"
              echo -e "+++++ $srp/$srr: Processing completed +++++"
            else
              force="TRUE"
              status="uncompleted"
              echo -e "ERROR! $srp/$srr has to restart the processing."
            fi
          else
            force="TRUE"
            status="uncompleted"
            echo -e "ERROR! $srp/$srr has to restart the processing."
          fi

        else
          sleep 60
        fi

      done
      echo >&1000
    } &

  fi

done <<<"$var_extract"

wait
echo "done"
