#!/bin/bash
#SBATCH --job-name=HHO-DD
#SBATCH --mem=10G
#SBATCH --time=1:00:00
#SBATCH --partition=short_cpu

lvl=$1
sds=$2
degree=$3
overlap=$4

module use /opt/share/sw2/modules/all/
module load silo gmsh Eigen MUMPS

PROGRAM=../../../build/apps/hho_dd/hhodd
OUTDIR=l${lvl}_s${sds}_k${degree}_o${overlap}
${PROGRAM} numtest_level${lvl}_${sds}sub.geo2s -d ${OUTDIR} -b bicgstab_eo.txt -s bicgstab_eo.silo
${PROGRAM} numtest_level${lvl}_${sds}sub.geo2s -d ${OUTDIR} -M -b bicgstab_mo.txt -s bicgstab_mo.silo
${PROGRAM} numtest_level${lvl}_${sds}sub.geo2s -d ${OUTDIR} -I -e error_eo.txt -s iterative_eo
${PROGRAM} numtest_level${lvl}_${sds}sub.geo2s -d ${OUTDIR} -I -M -e error_mo.txt -s iterative_mo

