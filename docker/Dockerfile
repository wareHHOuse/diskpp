FROM debian:bookworm-slim

RUN sed -i -e's/ main/ main contrib non-free/g' /etc/apt/sources.list.d/debian.sources
RUN apt update

run apt install -y apt-utils
RUN apt install -y build-essential
RUN apt install -y cmake
RUN apt install -y git
RUN apt install -y libeigen3-dev
RUN apt install -y libmumps-dev libmumps-seq-dev
RUN apt install -y liblua5.4-0 liblua5.4-dev
RUN apt install -y libsiloh5-0 libsilo-dev
RUN apt install -y gmsh libgmsh-dev
RUN apt install -y libmkl-dev

WORKDIR /diskpp
RUN git clone --recursive https://github.com/wareHHOuse/diskpp.git
WORKDIR /diskpp/diskpp/build
RUN cmake -DCMAKE_BUILD_TYPE=Release ..
RUN make -j $(nproc)

