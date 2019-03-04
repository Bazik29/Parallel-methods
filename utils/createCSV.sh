num_for_repeat=1
max_threds=4
step_treads=1
toFile=1

dims=(100 300 500 800 1000)


echo "Threads,Id,Steps,Time,Integral,AbsErr,RungeErr" > out.csv

# for ((omp_threads=1; omp_threads <= $max_threds; omp_threads+=step_treads)) 
# do
# 	export OMP_NUM_THREADS=$omp_threads
# 	for dim in ${dims[@]}
# 	do
# 		for ((i=0; i < $num_for_repeat; i++))
# 		do
# 				./../main -n $dim -l >> out.csv
# 		done
# 	done
# done

for ((i=0; i < $num_for_repeat; i++))
do
    for ((omp_threads=1; omp_threads <= $max_threds; omp_threads+=step_treads)) 
    do
        export OMP_NUM_THREADS=$omp_threads
        for dim in ${dims[@]}
        do
            ./../main -n $dim -l >> out.csv
		done
	done
done