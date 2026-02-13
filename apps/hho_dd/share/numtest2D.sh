set -x
for degree in 0 1 2 3; do
  for overlap in 2 4 8 16; do
    for lvl in 1 2 3; do
      for sds in 4 16; do
        OUTFILE=l${lvl}_s${sds}_k${degree}_o${overlap}.out
        sbatch --output ${OUTFILE} launchdd.sbatch $lvl $sds $degree $overlap
      done
    done
  done
done
