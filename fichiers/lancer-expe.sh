export OMP_NUM_THREADS
export GOMP_CPU_AFFINITY

ITE=$(seq 3) # nombre de mesures
THREADS=" 1 4 8 12 24" # nombre de threads
GOMP_CPU_AFFINITY=$(seq 0 23 ) # vérifier à l'aide de lstopo la bonne alternance des processeurs

PARAM="./prog -k sable -s 2048 -p alea -n -v " # parametres commun à toutes les executions 

execute (){
EXE="$PARAM $*"
OUTPUT="$(echo $* | tr -d ' ')"
for nb in $ITE; 
	do for OMP_NUM_THREADS in $THREADS ; 
		do echo -n "$OUTPUT $OMP_NUM_THREADS " >> ALL ; 
		$EXE 2>> ALL; 
	done; 
done
}

execute omp
execute omptiled
execute tasktiled
execute omp2tiled
execute task2tiled
