#!/usr/bin/env bash
if [ $# -lt 1 ];then
  echo "Usage:\n ./smore.sh model_name -train training_dataset -save embedding [model_options]"
  echo "Example:\n ./smore.sh hpe -train net.txt -save rep.txt"
  exit 1
fi
args=( "$@" )
for ((i=0; i < $#; i++)) ;do
  next_arg=$((i+1))
  if [ "${args[$i]}" == "-train" ] || [ "${args[$i]}" == "-save" ]; then
    args[${next_arg}]="data/"${args[${next_arg}]}
  fi
done
set "${args[@]}"
exec "./cli/$@"