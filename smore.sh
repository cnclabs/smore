#!/bin/sh
IMAGE=josix/smore:latest
docker run -it --name smore --rm --user="$(id -u):$(id -g)" -v "$PWD":/usr/local/smore/data "$IMAGE" "$@"