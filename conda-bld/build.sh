#!/bin/bash

set -ex

PYTHON_INCLUDE_DIR=$(${PYTHON} -c "import sysconfig; print(sysconfig.get_paths()['include'])")
NUMPY_INCLUDE_DIR=$(${PYTHON} -c "import numpy; print(numpy.get_include())")

cd $SRC_DIR/crpropa
mkdir build && cd build
cmake .. -G Ninja \
	-DCMAKE_PREFIX_PATH="${PREFIX}" \
	-DPython_EXECUTABLE="${PYTHON}" \
	-DPython_NumPy_INCLUDE_DIR="${NUMPY_INCLUDE_DIR}" \
	-DPython_INCLUDE_DIR="${PYTHON_INCLUDE_DIR}" \
	-DPython_INSTALL_PACKAGE_DIR="${SP_DIR}" \
	-DCMAKE_INSTALL_PREFIX="${PREFIX}" \
	-DBUILD_DOC=OFF \
	-DDOWNLOAD_DATA=OFF \
	-DENABLE_COVERAGE=OFF \
	-DENABLE_GIT=ON \
	-DENABLE_HDF5=ON \
	-DENABLE_OPENMP=ON \
	-DENABLE_PYTHON=ON \
	-DENABLE_QUIMBY=OFF \
	-DENABLE_SWIG_BUILTIN=ON \
	-DENABLE_TESTING=ON \
	-DFAST_WAVES="${FAST_WAVES}" \
	-DOMP_SCHEDULE=dynamic \
	-DSIMD_EXTENSIONS="${SIMD_EXTENSIONS}" \
	-DUSE_ABSOLUTE_RPATH=ON \
	-DCRPROPA_TESTS_PATH="${PREFIX}/share/crpropa/test/"
cmake --build .
cmake --install .
$PREFIX/bin/pybind11-stubgen -o ${SP_DIR} crpropa
# copy tests to share folder so user can test crpropa:
mkdir $PREFIX/share/crpropa/test/
for file in test*
do
	cp $file $PREFIX/share/crpropa/test/
done
# copy lib binaries
cp -r libs/ $PREFIX/share/crpropa/test/
# copy ctest instructions
cp CTestTestfile.cmake $PREFIX/share/crpropa/test/
# generate executable test file
echo ctest --test-dir $PREFIX/share/crpropa/test/ --output-on-failure --repeat until-pass:3 > $PREFIX/bin/testCRPropa
chmod +x $PREFIX/bin/testCRPropa

# Copy the [de]activate scripts to $PREFIX/etc/conda/[de]activate.d.
# This will allow them to be run on environment activation.
for CHANGE in "activate" "deactivate"
do
    mkdir -p "${PREFIX}/etc/conda/${CHANGE}.d"
    cp "${RECIPE_DIR}/${CHANGE}.sh" "${PREFIX}/etc/conda/${CHANGE}.d/${PKG_NAME}_${CHANGE}.sh"
done