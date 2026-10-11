#!/bin/bash
#SBATCH --job-name=dsB
#SBATCH --account=a_obrien
#SBATCH --partition=general
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=700G
#SBATCH --time=48:00:00
#SBATCH --output=/home/uqnzhai/depth_study/logs/B_%j.out
set -uo pipefail
source $HOME/depth_study/config.sh
T=$SLURM_CPUS_PER_TASK
OUTDIR=$WORK/kraken; mkdir -p $OUTDIR $RDM_OUT/kraken

# Self-contained: reads are regenerated here with the same seqtk seed used in
# stage A, used once, and deleted. Nothing is left on scratch between datasets.
LIST=${LIST:?set LIST=/path/to/task.tsv}
TIMING=$OUTDIR/timing_$SLURM_JOB_ID.tsv
echo -e "id\tsubsample_s\tkraken_s\texit" > $TIMING
echo "task list: $LIST  ($(wc -l < $LIST) datasets)"

while IFS=$'\t' read -r SAMPLE TIER SEED FULL R1 R2; do
  [ -z "${SAMPLE:-}" ] && continue
  ID=${SAMPLE}__t${TIER}__s${SEED}
  if [ -s "$RDM_OUT/kraken/${ID}.kreport" ]; then
    echo "=== $ID already done, skip"; continue
  fi
  D=$WORK/kraken_tmp/$ID
  rm -rf $D; mkdir -p $D
  echo "=== $ID  $(date '+%F %T')"

  FRAC=$(python3 -c "print(min(1.0, $TIER/$FULL))")
  t0=$SECONDS
  ( run $SIF_SEQTK seqtk sample -s$SEED "$R1" $FRAC | gzip -1 > $D/${ID}_1.fq.gz ) &
  ( run $SIF_SEQTK seqtk sample -s$SEED "$R2" $FRAC | gzip -1 > $D/${ID}_2.fq.gz ) &
  wait
  t_sub=$((SECONDS-t0))

  t0=$SECONDS
  run $SIF_KRAKEN2 kraken2 --db $KRAKEN_DB --threads $T --paired --gzip-compressed \
    --report $OUTDIR/${ID}.kreport --output /dev/null \
    $D/${ID}_1.fq.gz $D/${ID}_2.fq.gz
  rc=$?
  t_kr=$((SECONDS-t0))
  echo -e "$ID\t$t_sub\t$t_kr\t$rc" >> $TIMING

  if [ $rc -eq 0 ]; then
    run $SIF_BRACKEN bracken -d $KRAKEN_DB -i $OUTDIR/${ID}.kreport \
      -o $OUTDIR/${ID}.bracken_S.tsv -w $OUTDIR/${ID}.bracken_S.kreport \
      -r $BRACKEN_LEN -l S -t 10
    cp $OUTDIR/${ID}.kreport $OUTDIR/${ID}.bracken_S.tsv $OUTDIR/${ID}.bracken_S.kreport $RDM_OUT/kraken/
    echo "=== $ID done (subsample ${t_sub}s, kraken ${t_kr}s)"
  else
    echo "=== $ID FAILED rc=$rc"
  fi
  rm -rf $D
done < $LIST

rm -rf $WORK/kraken_tmp
echo; column -t $TIMING
echo STAGEB_DONE
