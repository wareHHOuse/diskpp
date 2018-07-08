#!/bin/sh

EIGEN_URL="http://bitbucket.org/eigen/eigen/get/3.3.4.tar.bz2"
EIGEN_FILE=$(basename ${EIGEN_URL})

wget ${EIGEN_URL}
tar -jxvf ${EIGEN_FILE}

EIGEN_TAR_DIR=$(tar -tf ${EIGEN_FILE} | cut -f 1 -d '/' | uniq)

mv ${EIGEN_TAR_DIR} eigen

rm ${EIGEN_FILE}
