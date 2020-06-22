FROM ubuntu:18.04
RUN apt-get update && apt-get upgrade -y \
    && apt-get install -y build-essential automake autoconf libtool zlib1g-dev libbz2-dev \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*
COPY . /build
WORKDIR /build
RUN ./autogen.sh
RUN ./configure --disable-shared
RUN make -j4
RUN make check
RUN make install
