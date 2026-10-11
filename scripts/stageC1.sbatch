#!/bin/bash
#SBATCH --job-name=dsC1
#SBATCH --account=a_obrien
#SBATCH --partition=general
#SBATCH --nodes=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=200G
#SBATCH --time=24:00:00
#SBATCH --output=/home/uqnzhai/depth_study/logs/C1_%j.out
set -uo pipefail
source $HOME/depth_study/config.sh
T=$SLURM_CPUS_PER_TASK
W=$WORK/checkm2; rm -rf $W; mkdir -p $W/bins
echo "copying bins to scratch..."
cp $RDM_OUT/bins_all/*.fa $W/bins/
echo "bins: $(ls $W/bins/*.fa | wc -l)"
cd $W
run $SIF_CHECKM2 checkm2 predict -i bins -x fa -o checkm2 \
  --database_path $CHECKM2_DIR/uniref100.KO.1.dmnd -t $T --force
rc=$?
echo "checkm2 exit $rc"
if [ $rc -eq 0 ]; then
  mkdir -p $RDM_OUT/checkm2
  cp checkm2/quality_report.tsv $RDM_OUT/checkm2/
  echo "=== quality summary:"
  awk -F'\t' 'NR>1{n++; if($2>=50 && $3<=10) mq++; if($2>=90 && $3<=5) hq++} END{print "total",n," MQ(>=50,<=10)",mq," HQ(>=90,<=5)",hq}' checkm2/quality_report.tsv
fi
cd $WORK && rm -rf $W
echo STAGEC1_DONE
