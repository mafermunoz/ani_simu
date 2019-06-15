for i in {0..99}
do
  echo $i
  sbatch submit_ani_simu.sh  $i

done
