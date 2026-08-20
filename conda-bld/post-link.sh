#!/bin/bash

mkdir -p  $PREFIX/share/crpropa &> /dev/null  # eliminate possible warning outputs
cd $PREFIX/share/crpropa
#define CRPROPA_DATA_VERSION "2024-04-30"
CRPROPA_DATAFILE_VER=$(cat $PREFIX/include/crpropa/Version.h | grep --extended-regexp "#define +CRPROPA_DATA_VERSION +.+" | sed -E "s/#define +CRPROPA_DATA_VERSION +\"//g" | sed -E "s/\"//g")
USERPWD=Nag8KQLoCjpG5es

# download data
curl -u $USERPWD:$USERPWD -O https://ruhr-uni-bochum.sciebo.de/public.php/webdav/data-${CRPROPA_DATAFILE_VER}.tar.gz-CHECKSUM
curl -u $USERPWD:$USERPWD -O https://ruhr-uni-bochum.sciebo.de/public.php/webdav/data-${CRPROPA_DATAFILE_VER}.tar.gz

# check download hash, command depends on OS:
OS=$(uname)
if [ $OS = "Linux" ]
then
	if [ "$(cat data-${CRPROPA_DATAFILE_VER}.tar.gz-CHECKSUM)" = "$(md5sum data-${CRPROPA_DATAFILE_VER}.tar.gz)" ]
	then
		:
	else
		echo Checksums are not the same!
		exit 1
	fi
elif [ $OS = "Darwin" ]
then
	if [ $(cat data-${CRPROPA_DATAFILE_VER}.tar.gz-CHECKSUM | cut -d ' ' -f 1) = $(md5 data-${CRPROPA_DATAFILE_VER}.tar.gz | cut -d ' ' -f 4) ]
	then
		:
	else
		echo Checksums are not the same!
		exit 1
	fi
fi

tar xzf data-${CRPROPA_DATAFILE_VER}.tar.gz
mv data/* $PREFIX/share/crpropa/
# cleanup
rm -rf data data-${CRPROPA_DATAFILE_VER}.tar.gz data-${CRPROPA_DATAFILE_VER}.tar.gz-CHECKSUM
