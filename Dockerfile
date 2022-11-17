FROM ubuntu:latest

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends build-essential r-base \
     r-cran-randomforest python3.9 python3-pip python3-setuptools python3-dev python3-tk
RUN apt-get update && apt-get -y install cmake protobuf-compiler

WORKDIR /app

ADD . /app/
RUN pip install -r requirements.txt