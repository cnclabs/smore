#!/bin/sh
IMAGE=josix/smore:latest
docker run -it --name smore --rm -v "$PWD":/usr/local/smore/data "$IMAGE" "$@"