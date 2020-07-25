#!/bin/bash

if [ -e "jobs.id" ]; then
 rm jobs.id
fi

current_dcd=$2
array_str=$2
COUNTER=1
while [ $COUNTER -lt 8 ]; do
    #echo The counter is $COUNTER
    current_dcd=$(($current_dcd + 10))
    array_str="$array_str,$current_dcd"
    let COUNTER=COUNTER+1 
done

RES=$(sbatch --array $array_str --job-name=$4 submit.sh $1 $3) && echo ${RES##* } >> jobs.id 

sleep 5

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
