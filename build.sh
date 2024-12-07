cmake -S cmake -B cmake_build
cmake --build cmake_build
cmake --install cmake_build --config Release --prefix ./p528install
p528install/bin/p528GTest