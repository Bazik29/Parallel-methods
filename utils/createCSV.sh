num_for_repeat=1000
max_threds=4
step_treads=1

dims=(10 500 1000 2000)
methods=(0 1 2 3 4 5)


echo "Threads,Id,Steps,Time,Integral,AbsErr,RungeErr" > out1.csv

for ((i=0; i < $num_for_repeat; i++))
do
    for ((omp_threads=1; omp_threads <= $max_threds; omp_threads+=step_treads)) 
    do
        export OMP_NUM_THREADS=$omp_threads
        for dim in ${dims[@]}
        do
            for m in ${methods[@]}
            do
                ./../test1 -n $dim -m $m >> out1.csv
            done
		done
	done
done