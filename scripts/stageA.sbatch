#!/bin/bash
#SBATCH --job-name=dsA
#SBATCH --account=a_obrien
#SBATCH --partition=general
#SBATCH --nodes=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=120G
#SBATCH --time=16:00:00
#SBATCH --output=/home/uqnzhai/depth_study/logs/A_%A_%a.out
set -uo pipefail
source $HOME/depth_study/config.sh
export -f run run_argoap
T=$SLURM_CPUS_PER_TASK

LIST=${LIST:?set LIST=/path/to/task.tsv}
read SAMPLE TIER SEED FULL R1 R2 < <(sed -n "${SLURM_ARRAY_TASK_ID}p" $LIST)
[ -z "${SAMPLE:-}" ] && { echo "empty task line"; exit 1; }
FRAC=$(python3 -c "print(min(1.0, $TIER/$FULL))")
ID=${SAMPLE}__t${TIER}__s${SEED}
OUT=$WORK/$ID
rm -rf $OUT; mkdir -p $OUT/reads $OUT/rgi $OUT/idx $OUT/bins
cd $OUT
echo "id=$ID full=${FULL}Gbp frac=$FRAC cpus=$T node=$(hostname)"

TIMING=$OUT/timing.tsv
echo -e "step\tseconds\tmaxrss_kb\texit" > $TIMING
step() {
  local name=$1; shift
  echo "=== $name  $(date '+%F %T')"
  /usr/bin/time -f "%e\t%M" -o .t bash -c '"$@"' _ "$@"; local rc=$?
  echo -e "$name\t$(tail -1 .t)\t$rc" >> $TIMING
  echo "=== $name done rc=$rc"
}

step subsample bash -c "
  (run $SIF_SEQTK seqtk sample -s$SEED $R1 $FRAC | gzip -1 > reads/${ID}_1.fq.gz) &
  (run $SIF_SEQTK seqtk sample -s$SEED $R2 $FRAC | gzip -1 > reads/${ID}_2.fq.gz) &
  wait"
run $SIF_SEQTK seqtk size reads/${ID}_1.fq.gz | awk '{print "reads_R1="$1, "bases_R1="$2}' | tee size.txt
run $SIF_SEQTK seqtk size reads/${ID}_2.fq.gz | awk '{print "reads_R2="$1, "bases_R2="$2}' | tee -a size.txt

step argoap_s1 run_argoap args_oap stage_one -i reads -o argoap -f fq.gz -t $T
step argoap_s2 run_argoap args_oap stage_two -i argoap -t $T

step megahit run $SIF_MEGAHIT megahit -1 reads/${ID}_1.fq.gz -2 reads/${ID}_2.fq.gz \
  -o megahit -t $T --min-contig-len 1000
awk -v id=$ID '/^>/{split($1,a," "); print ">" id "__" substr(a[1],2); next}{print}' \
  megahit/final.contigs.fa > ${ID}.contigs.fa
run $SIF_SEQTK seqtk comp ${ID}.contigs.fa | awk '{n++; s+=$2} END{print "contigs="n, "assembly_bp="s}' | tee -a size.txt
rm -rf megahit

ln -sfn $CARD_DIR/localDB localDB
step rgi run $SIF_RGI rgi main -i ${ID}.contigs.fa -o rgi/${ID} -t contig -a DIAMOND -n $T --local --clean

step genomad run $SIF_GENOMAD genomad end-to-end --cleanup --splits 16 -t $T ${ID}.contigs.fa genomad $GENOMAD_DB

step bt2_build run $SIF_BOWTIE2 bowtie2-build --threads $T ${ID}.contigs.fa idx/${ID}
step bt2_map bash -c "
  run $SIF_BOWTIE2 bowtie2 -x idx/${ID} -1 reads/${ID}_1.fq.gz -2 reads/${ID}_2.fq.gz -p $T --no-unal 2> bt2.log \
  | run $SIF_SAMTOOLS samtools sort -@ 12 -m 2G -o ${ID}.bam - \
  && run $SIF_SAMTOOLS samtools index ${ID}.bam"

step jgi_depth run $SIF_METABAT2 jgi_summarize_bam_contig_depths --outputDepth depth.txt ${ID}.bam
step metabat2 run $SIF_METABAT2 metabat2 -i ${ID}.contigs.fa -a depth.txt -o bins/${ID}.bin -m 1500 -t 16 --seed 1
NB=$(ls bins/*.fa 2>/dev/null | wc -l); echo "bins=$NB" | tee -a size.txt

# ship results to RDM. RDM refuses chgrp/chmod, so never touch ownership or permissions.
mkdir -p $RDM_OUT/stageA/$ID $RDM_OUT/bins_all
rsync -rltD --no-perms --no-group --no-owner \
  --exclude reads --exclude idx --exclude '*.bam' --exclude '*.bai' \
  $OUT/ $RDM_OUT/stageA/$ID/ || echo "rsync exit $? (non-fatal)"
if [ "$NB" -gt 0 ]; then for b in bins/*.fa; do cp $b $RDM_OUT/bins_all/${ID}.$(basename $b); done; fi

echo; column -t $TIMING; echo; cat size.txt
echo "STAGEA_DONE $ID"

# always free scratch, whatever happened above
cd $WORK && rm -rf $OUT && echo "scratch cleaned: $ID"
