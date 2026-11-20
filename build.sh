cmake --preset release64 -DBUILD_DOCS=OFF
cmake --build --preset release64
ctest --preset release64
