#!/bin/bash
c++ checkFiles.cpp -o CheckFiles
c++ Serial.cpp -o Serial
c++ GraphGenerator.cpp -o GraphGenerator
c++ Parallel_OMP.cpp -o Parallel_OMP -fopenmp
mpic++ Parallel_MPI.cpp -o Parallel_MPI

#NUM_PROC=3
#k=64
# shellcheck disable=SC2188
>Results.txt
rm -rf Outputs
rm -rf Inputs
mkdir -p "Outputs"
mkdir -p "Inputs"

#VERTICES=$((NUM_PROC * k))
allProcNo=(2 3 4 8)
allVertices=(96 192 384 768)

runGraphTest() {
  graph=$(./GraphGenerator "$1")

  serial=$(./Serial "$graph")
  sleep 1
  omp=$(OMP_NUM_THREADS=$2 ./Parallel_OMP "$graph")
  sleep 1
  mpi=$(mpirun -oversubscribe -np "$2" --mca opal_warn_on_missing_libcuda 0 ./Parallel_MPI "$graph")
  sleep 1
  ./CheckFiles "$serial" "$mpi" >>Results.txt
  ./CheckFiles "$serial" "$omp" >>Results.txt
  echo -e "\n" >>Results.txt
}

for i in "${allVertices[@]}"; do
  echo -e "Graph Size" "$i" >>Results.txt
  newName="${i}Results.txt"
  for j in "${allProcNo[@]}"; do
    echo -e "Number of Processors:" "$j" "\n" >>Results.txt
    runGraphTest "$i" "$j"
  done
  mv Results.txt "${newName}"
done

#Move into output directory
for f in *.txt; do
  if [[ "$f" == *"_Output"* ]]; then
    mv "$f" Outputs/
  elif [[ "$f" == *"Results"* ]]; then
    mv "$f" Results/
  elif [[ "$f" == *"_Graph"* ]]; then
    mv "$f" Inputs/
  fi
done

echo Success
