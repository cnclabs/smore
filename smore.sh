#!/bin/sh
if [ $# -lt 1 ];then
  echo "Usage:\n ./smore.sh model_name -train training_dataset -save embedding [model_options]"
  echo "Example:\n ./smore.sh hpe -train net.txt -save rep.txt"
  exit 1
fi
IMAGE=josix/smore:latest
docker run -it --name smore --rm --user="$(id -u):$(id -g)" -v "$PWD":/usr/local/smore/data "$IMAGE" "$@"