#!/usr/bin/env bash
args=( "$@" )
for ((i=0; i < $#; i++)) ;do
  next_arg=$((i+1))
  if [ "${args[$i]}" == "-train" ] || [ "${args[$i]}" == "-save" ]; then
    args[${next_arg}]="data/"${args[${next_arg}]}
  fi
done
set "${args[@]}"
exec "./cli/$@"