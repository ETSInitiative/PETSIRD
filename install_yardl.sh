#!/usr/bin/env bash

# Install yardl
# you can set version etc as follows
# env YARDL_VERSION:=0.6.7 DEST=yardl ./install_yardl.sh

# default values
: "${YARDL_VERSION:=0.6.7}"
: "${DEST:=yardl}"

set -e

# TODO will fail on Windows as no uname there
OS="$(uname -s)"
ARCH="$(uname -m)"
echo $OS
echo $ARCH
case "$OS" in
  Linux)
    YARDL_OS="linux"
    SUFFIX=".tar.gz"
    ;;
  Darwin)
    YARDL_OS="darwin"
    SUFFIX=".tar.gz"
    ;;
  Windows)
    YARDL_OS="windows"
    SUFFIX=".zip"
    ;;
esac
DL_SUFFIX="_${YARDL_OS}_${ARCH}${SUFFIX}"
FNAME="yardl_${YARDL_VERSION}${DL_SUFFIX}"

if command -v gh > /dev/null
then
  gh release download -R microsoft/yardl -p "*$DL_SUFFIX" v$YARDL_VERSION
elif command -v curl > /dev/null
then
  curl -L https://github.com/microsoft/yardl/releases/download/v${YARDL_VERSION}/${FNAME} -o ${FNAME}
elif command -v wget > /dev/null
then
  wget "https://github.com/microsoft/yardl/releases/download/v${YARDL_VERSION}/${FNAME}"
else
  echo "Need either gh or curl or wget"
  exit 1
fi

mkdir "${DEST}"
if test "${YARDL_OS}" = windows; then
  unzip "${FNAME}" -d "${DEST}"
else
  tar -xzf "${FNAME}" -C "${DEST}"
fi
rm "${FNAME}"

echo "Installed version in ${DEST}:"
ls -l "${DEST}"
