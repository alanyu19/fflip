#!/bin/bash

#if [ -e "jobs.id" ]; then
# rm jobs.id
#fi

current_dcd=$2
array_str=$2
COUNTER=1
# NUMBER in the next line is the total block number
while [ $COUNTER -lt $5 ]; do
    # the number in the next line is block size (in unit of number of trj files)
    # this should be provided by the python code, however this is fine for now...
    current_dcd=$(($current_dcd + 1))
    array_str="$array_str,$current_dcd"
    let COUNTER=COUNTER+1
done

sbatch --array $array_str --job-name=$4 submit.sh $1 $3

sleep 3

input="jobs.id"
index=0
while IFS= read -r line
do
  if [ $index -eq 0 ]; then
    idstr="$line"
  else
    idstr="$idstr:$line"
  fi
  let "index=index+1"
done < "$input"

sbatch --dependency=afterok:$idstr --job-name=fc-$4 finish.sh
