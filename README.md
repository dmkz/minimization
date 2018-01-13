# Глобальная минимизация многомерной функции

Пространство поиска заполняется равномерной сеткой, а затем из каждого узла этой сетки запускаются два разных алгоритма поиска локального минимума. Все локальные минимумы сравниваются между собой и из них выбирается глобальный.

## Виды сеток (в порядке эффективности):
* Sobol
* Faure
* Niederreiter

## Алгоритмы локальной минимизации:
1. BFGS
2. Hessian Free

Ограничения: размерность m от 32 до 64, количество узлов сетки 2^m.

Статус сборки для dmkz/minimization/develop: [![Build Status](https://ci.worldfly.org/buildStatus/icon?job=Minimization main)](https://ci.worldfly.org/job/Minimization%20main/)

Статус сборки для worldfly/minimization/develop: [![Build Status](https://ci.worldfly.org/buildStatus/icon?job=Minimization)](https://ci.worldfly.org/job/Minimization/)
