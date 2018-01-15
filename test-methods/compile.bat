@echo on
cd build
del main.exe
cmake -DCMAKE_BUILD_TYPE=Debug -G "MSYS Makefiles" .. && cmake --build . --config Debug -- -j6
make
cd ..