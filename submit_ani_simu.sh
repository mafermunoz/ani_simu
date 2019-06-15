#!/bin/bash
#SBATCH --partition=rhel6-long
#SBATCH --ntasks=1
#SBATCH --mem=30G
#SBATCH --job-name=Ani


export HOME=/atlas/users/mmunozsa/

source /cvmfs/dampe.cern.ch/rhel6-64/etc/setup.sh
source ~/astro/bin/activate


python ani_simu_2.py
