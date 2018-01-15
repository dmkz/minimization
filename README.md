# Глобальная минимизация многомерной функции

Пространство поиска заполняется равномерной сеткой, а затем из каждого узла этой сетки запускаются три разных алгоритма поиска локального минимума. Все локальные минимумы сравниваются между собой и из них выбирается глобальный.

## Виды сеток:
* Sobol

## Алгоритмы локальной минимизации:
1. BFGS
2. Hessian Free
3. Nesterov

Статус сборки : [![Build Status](https://ci.worldfly.org/buildStatus/icon?job=Minimization-Container/Minimization-main)](https://ci.worldfly.org/job/Minimization-Container/job/Minimization-main/)

## Инструкция по возможной сборке проекта на Windows 7 и выше

В одном из вариантов сборки проекта на Windows предполагается наличие следующих утилит:
* Утилита cmake
* Утилита make
* Компилятор MinGW x64 (желательно последней версии) с поддержкой posix threads

Рассмотрим самый быстрый процесс установки всего необходимого для сборки и запуска проекта в данном варианте. Для этого понадобится только утилита msys2, которую можно скачать с официального сайта http://www.msys2.org/.

После установки msys2 необходимо открыть терминал cmd.exe и запустить msys2 командой:

    C:\msys64\msys2_shell.cmd -mingw64

Затем необходимо выполнить обновление баз данных пакетов:

    pacman -Syu

После успешного завершения данной команды может потребоваться закрытие окна терминала msys2 с последующим повторным открытием командой: 

    C:\msys64\msys2_shell.cmd -mingw64

Далее необходимо запустить полное обновление всех пакетов командой: 

    pacman -Su

И установить требуемые для сборки проекта утилиты командой:

    pacman -S mingw-w64-x86_64-gcc make mingw-w64-x86_64-cmake

Процесс установки и обновления всех пакетов завершен. Чтобы собрать проект и запустить, можно выполнить скрипт compile.bat в директории проекта.

## Инструкция по сборке на *nix

Для сборки требуется:
* cmake ( >= 3.1 )
* GNU Make ( >= 3.81 )
* gcc ( >= 7.2 )

Сборка проекта:

    cmake CMakeLists.txt
    make


Запуск проекта:

    ./main

## Интсрукция по сборке на Windows (вариант 2)

> Протестировано на Windows Server 2016.

Для сборки требуется:
* cmake ( >= 3.1 )
* Visual C++ Build Tools ( >= 2015 )
* Windows 8.1 или Windows 10 SDK

Сборка проекта в cmd:

    setlocal enabledelayedexpansion
    cmake CMakeLists.txt
    "C:\Program Files (x86)\Microsoft Visual C++ Build Tools\vcbuildtools_msbuild.bat"
    msbuild main.vcxproj
    
Запуск проекта в cmd:

    Debug\main.exe
    
Для запуска проекта может потребоваться ```ucrtbased.dll```. Проблему можно решить, скопировав файл ```copy "C:\Program Files (x86)\Windows Kits\10\bin\x86\ucrt\ucrtbased.dll"``` в ```Debug\ucrtbased.dll```.
