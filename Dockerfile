FROM ubuntu:18.04

ENV DEBIAN_FRONTEND noninteractive
ENV BASE_DIR /usr/local
ENV APP_DIR $BASE_DIR/smore

RUN apt-get update \
  && apt-get -y install g++-7 make\
  && apt-get clean \
  && rm -rf /var/lib/apt/lists/*

COPY cli $APP_DIR/cli
COPY src $APP_DIR/src
COPY ./entrypoint.sh $APP_DIR
COPY Makefile $APP_DIR
WORKDIR $APP_DIR

RUN make CC=g++-7
ENTRYPOINT ["./entrypoint.sh"]
VOLUME ["$APP_DIR/data"]