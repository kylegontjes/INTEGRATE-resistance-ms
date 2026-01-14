for i in $(ls)
 do echo $i
    sbatch $i
done