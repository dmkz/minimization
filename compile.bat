@echo off
PATH | findstr /I GnuWin32\bin >nul &&echo Check PATH to make.exe: OK ||echo Check path to make.exe: ERROR
PATH | findstr /I CMake\bin >nul &&echo Check PATH to cmake.exe: OK ||echo Check PATH to cmake.exe: ERROR
PATH | findstr /I MinGW\bin >nul &&echo Check PATH to MinGW compiler: OK||echo Check PATH to MinGW compiler: ERROR

@echo on
cmake -G "MinGW Makefiles"
make
main