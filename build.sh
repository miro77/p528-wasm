cmake -S cmake -B cmake_build
cmake --build cmake_build
cmake --install cmake_build --config Debug --prefix ./p528install
p528install/bin/p528GTest