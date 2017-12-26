# Глобальная минимизация многомерной функции

Пространство поиска заполняется равномерной сеткой, а затем из каждого узла этой сетки запускаются два разных алгоритма поиска локального минимума. Все локальные минимумы сравниваются между собой и из них выбирается глобальный.

## Виды сеток (в порядке эффективности):
* Sobol
* Faure
* Niederreiter

## Алгоритмы локальной минимизации:
1. BFGS
2. Hessian Free
3. Nesterov

Ограничения: размерность m от 32 до 64, количество узлов сетки 2^m.
Статус сборки для worldfly/minimization/develop: ![alt text](https://travis-ci.org/worldfly/minimization.svg?branch=develop)
