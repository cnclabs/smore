#!/bin/sh
IMAGE=smore:latest
docker run -it --name smore --rm -v "$PWD":/usr/local/smore/data "$IMAGE" "$@"