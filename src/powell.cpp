#include "powell.hpp"

IterationData powell(Function func, Vector p, const StopCondition& stop_condition) {
// func - указатель на целевую функцию
// p - начальное приближение
// stop_condition - критерий остановы
// Результат работы метода будет лежать в структуре данных о последней итерации
    
    int i, j, ibig; //переменные для циклов
    Real del, fp, fptt, t; // del-дельта , fp и fptt - ф-ция в р и рtt
    
    int n = p.size(); // размерность вектора  p
    std::pair<Vector, Vector> ans;// Переменная, хранит ответ в виде пары векторов

    Matrix xi(n, Vector(n));//вектор, задающий направление поиска
    // Заполняем массив; xi - единичная матрица
    for (int i = 0; i < n; ++i) {
        xi[i][i] = 1.0;
    }

    Vector pt(n), ptt(n), xit(n);
    Real fret = func(p);// fret - хранит значение функции в точке p
    
    // Инициализируем начальной точкой структуру данных итерации:
    IterationData iter_data;
    iter_data.x_curr = p;
    iter_data.f_curr = fret;
    iter_data.iter_counter = 0;
    iter_data.method_title = "Powell";
    
    pt = p; //сохраняем начальную точку
    
    do { 
        fp = fret;
        ibig = 0; // номер направления
        del = 0.0; // разница в предыдущем иксе и нынешнем
        for (i = 0; i<n; i++) {
        // в каждой итерации, цикл по всем направлениям
            for (j = 0; j<n; j++) xit[j] = xi[j][i]; /* копируем направления*/
            fptt = fret;
            ans = linmin(p, xit, func); //минимизируем
            p = ans.first;
            xit = ans.second;
            fret = func(p);
            if (fptt - fret > del) {/* записываем del, если это наибольшее уменьшение*/
                del = fptt - fret;
                ibig = i + 1;
            }
        }
        
        for (j = 0; j < n; j++) {
            ptt[j] = 2.0*p[j] - pt[j]; /* меняем направление и  сохраняем старую точку*/
            xit[j] = p[j] - pt[j];
            pt[j] = p[j];
        }
        fptt = func(ptt); // значение функции в новой точке
        
        // Переходим к следующей итерации:
        iter_data.next(ptt, fptt);
        
        if (fptt < fp) {
            t = 2.0*(fp - 2.0*fret + fptt)*SQR(fp - fret - del) - del * SQR(fp - fptt);
            if (t < 0.0) {
                ans = linmin(p, xit, func); /* переходим к минимуму нового направления и сохраняем его*/
                p = ans.first;
                xit = ans.second;
                fret = func(p);
                for (j = 0; j<n; j++) {
                    xi[j][ibig - 1] = xi[j][n - 1];
                    xi[j][n - 1] = xit[j];
                }
            }
        }
    } while (!stop_condition(iter_data)); // Проверяем условие останова
    return iter_data;
}


Real f1dim(const Real x, int ncom, Vector *pcom_p, Vector *xicom_p, Function nrfunc) //Мы можем использовать одномерные функции минимизации, построив функтор F1DIM, 
    //который дает значение нашей функции (func), вдоль линии, проходящей через точку p в направлении xi.
    //LinmiN вызывает брент и инструктирует его минимизировать F1DIM.
{
    Vector xt(ncom);//задаем вектор определенной размерности
    Vector &pcom = *pcom_p, &xicom = *xicom_p; //передача значений по ссылке
    for (int j = 0; j<ncom; j++)
        xt[j] = pcom[j] + x * xicom[j];//заполнение вектора xt
    Real ans = nrfunc(xt); // Значение функции от полученного вектора xt
    return ans;
}

std::pair<Vector, Vector> linmin(Vector &pInit, Vector &xiInit, Function func)
{/*учитывая N мерную точку P и N мерное направление Xi, перемещает  P, где функция (f1dim) принимает минимум из точки P, вдоль направления Xi и заменяет значение Xi на фактическое перемещение вектора, которым была перемещена точка P . fret, как и (f1dim)возвращает  такое же значение - P. Это осуществляется путем вызова подпрограмм mnbrak и brent(перевёл как смог)*/
    int ncom;
    Vector *pcom_p, *xicom_p;

    int j;
    const Real TOL = 1.0e-8;
    Real xx, xmin, bx, ax;
    // Real fx, fb, fa; // не используется
    Vector p = pInit; Vector xi = xiInit;
    int n = p.size();
    ncom = n; //Задаем размерность
    pcom_p = new Vector(n); //Cоздаём начальную точку в пространстве размером n;
    xicom_p = new Vector(n); //Cоздаём направляющий вектор в пространстве размером n;

    Vector &pcom = *pcom_p, &xicom = *xicom_p;
    for (j = 0; j<n; j++) {
        pcom[j] = p[j];
        xicom[j] = xi[j];
    }
    ax = 0.0; //Инициализация
    xx = 1.0; //Инициализация
    std::pair <Vector, Vector> ans = mnbrak(ax, xx, f1dim,ncom,pcom_p,xicom_p,func);
    Vector temp = ans.first;
    ax = temp[0]; xx = temp[1]; bx = temp[2];
    temp = ans.second;
    // fa = temp[0]; fx = temp[1]; fb = temp[2];
    xmin = brent(ax, xx, bx, f1dim, TOL, ncom, pcom_p, xicom_p, func);
    for (j = 0; j<n; j++) {
        xi[j] *= xmin;
        p[j] += xi[j];
    }
    delete xicom_p; //освобождаем память
    delete pcom_p; //освобождаем память
    return { std::move(p), std::move(xi) };
}


/* Метод минимизации одномерных функций, с помоцью параболической интерполяции */
Real brent(const Real ax, const Real bx, const Real cx, Real f(const Real, int, Vector *, Vector *, Function),  
// минимум заключен в между a и b; x - это точка с наименьшим значением функции
// точка с наименьшим значением функции, найденным до сих пор
    const Real tol, int ncom, Vector *pcom_p, Vector *xicom_p, Function func)
{
    /*Для функции f, и набора абсцисс ax,bx,cx
    (такие ,что bx находистя между ax и cx, а f(bx)
    меньше, чем f(аx) и f(сx) вместе взятые, что помогает провести расчеты до дробногй точности tol. 
    По оси абсцисс
    минимум возвращается как x, а значение минимальной функции возвращается как brent,
    возвращаемое значение функции.*/
    const int ITMAX = 100; //заданный максимум итераций 
    const Real CGOLD = 0.3819660;
    const Real ZEPS = 1.0e-3;
    int iter;
    Real a, b, d = 0.0, etemp, fu, fv, fw, fx; /*
    w - точка со вторым наименьшим значением функции; v - предыдущее
        значение w; u - это точка, в которой функция была оценена на последней итерации.*/
    Real p, q, r, tol1, tol2, u, v, w, x, xm;
    Real e = 0.0; //Расстояние, на которое мы продвинулись

    a = (ax < cx ? ax : cx);
    b = (ax > cx ? ax : cx);
    x = w = v = bx;
    fw = fv = fx = f(x,ncom,pcom_p,xicom_p,func);
    for (iter = 0; iter<ITMAX; iter++) {
        xm = 0.5*(a + b);//cередина отрезка
        tol2 = 2.0*(tol1 = tol * std::abs(x) + ZEPS);
        if (std::abs(x - xm) <= (tol2 - 0.5*(b - a))) { // Критерий остановки !!! 
            return x;
        }
        if (std::abs(e) > tol1) { // условие для создания трех точек, для возможной параболичсекой интерполяции
            r = (x - w)*(fx - fv);
            q = (x - v)*(fx - fw);
            p = (x - v)*q - (x - w)*r;
            q = 2.0*(q - r);
            if (q > 0.0) p = -p;
            q = std::abs(q);
            etemp = e;
            e = d;
            if (std::abs(p) >= std::abs(0.5*q*etemp) || p <= q * (a - x) || p >= q * (b - x)) // условие, проверяющее возможность параболической интерполяции
                d = CGOLD * (e = (x >= xm ? a - x : b - x));
            else {
                d = p / q;
                u = x + d;
                if (u - a < tol2 || b - u < tol2)
                    d = SIGN(tol1, xm - x);
            }
        }
        else {
            d = CGOLD * (e = (x >= xm ? a - x : b - x));
        }
        u = (std::abs(d) >= tol1 ? x + d : x + SIGN(tol1, d));
        fu = f(u, ncom, pcom_p, xicom_p, func);
        if (fu <= fx) {
            if (u >= x) a = x; else b = x;
            //shft3(v, w, x, u);
            v = w;
            w = x;
            x = u;
            //shft3(fv, fw, fx, fu);
            fv = fw;
            fw = fx;
            fx = fu;
        }
        else {
            if (u < x) a = u; else b = u;
            if (fu <= fw || w == x) {
                v = w;
                w = u;
                fv = fw;
                fw = fu;
            }
            else if (fu <= fv || v == x || v == w) {
                v = u;
                fv = fu;
            }
        }
    }
    return x;
}


std::pair<Vector, Vector> mnbrak(Real axInit, Real bxInit, Real func(const Real, int, Vector *, Vector *, Function), int ncom, Vector *pcom_p, Vector *xicom_p, Function nrfunc)
{
    const Real GOLD = 1.618034, GLIMIT = 100.0, TINY = 1.0e-14;
    Real ulim, u, r, q, fu;
    Real fa, fb, fc, cx;
    Real ax = axInit;
    Real bx = bxInit;
    fa = func(ax, ncom, pcom_p, xicom_p, nrfunc); // значение функции в точках а и б
    fb = func(bx, ncom, pcom_p, xicom_p, nrfunc);
    if (fb > fa) {/* меняет значения ax и bx,fb и fa, чтобы мы могли идти вглубь по направлению от a до b*/
        SWAP(ax, bx);
        SWAP(fb, fa);
    }
    cx = bx + GOLD * (bx - ax);// Инициализация сx
    fc = func(cx, ncom, pcom_p, xicom_p, nrfunc); /* Вычисление значения функции в точке cx*/
    int iter = 0;
    while (fb > fc && iter++ < 100) {/* Вычисление значения переменной u с помощью параболической экстраполяции от ax,bx,cx. Переменная  TINY нужна для предотвращения возможного деления на ноль.*/
        r = (bx - ax)*(fb - fc);
        q = (bx - cx)*(fb - fa);
        u = bx - ((bx - cx)*q - (bx - ax)*r) /
            (2.0*SIGN(MAX(std::abs(q - r), TINY), q - r));
        ulim = bx + GLIMIT * (cx - bx); /* Дальше не идём. Теперь проверяем возможность  того что  u находится между bx и cx*/
        if ((bx - u)*(u - cx) > 0.0) {
            fu = func(u, ncom, pcom_p, xicom_p, nrfunc);
            if (fu < fc) {// получаем минимум между bx и cx
                ax = bx;
                bx = u;
                fa = fb;
                fb = fu;
                Vector x = { ax,bx,cx }; Vector f = { fa,fb,fc };
                return { x,f };
            }
            else if (fu > fb) {// получаем минимум между ax и u
                cx = u;
                fc = fu;
                Vector x = { ax,bx,cx }; Vector f = { fa,fb,fc };
                return { std::move(x), std::move(f) };
            }
            u = cx + GOLD * (cx - bx); /*параболическое приближение не требуется. используем формулу*/
            fu = func(u, ncom, pcom_p, xicom_p, nrfunc);
        }
        else if ((cx - u)*(u - ulim) > 0.0) {/*параболическое приближение на участке между сx и её допустимым пределом*/
            fu = func(u, ncom, pcom_p, xicom_p, nrfunc); /*значение функции в точке u*/
            if (fu < fc) {
                //shft3(bx, cx, u, u + GOLD * (u - cx));
                bx = cx;
                cx = u;
                u = u + GOLD * (u - cx);
                //shft3(fb, fc, fu, func(u, ncom, pcom_p, xicom_p, nrfunc));
                fb = fc;
                fc = fu;
                fu = func(u, ncom, pcom_p, xicom_p, nrfunc);
            }
        }
        else if ((u - ulim)*(ulim - cx) >= 0.0) { /*ограничить параболическое u до максимально допустимого значения*/
            u = ulim;
            fu = func(u, ncom, pcom_p, xicom_p, nrfunc);
        }
        else {//Используем формулу
            u = cx + GOLD * (cx - bx);
            fu = func(u, ncom, pcom_p, xicom_p, nrfunc);
        }
        //shft3(ax, bx, cx, u);
        ax = bx;
        bx = cx;
        cx = u;
        //shft3(fa, fb, fc, fu); //убрать последнюю точку и продолжить
        fa = fb;
        fb = fc;
        fc = fu;
    }
    Vector x = { ax,bx,cx }; Vector f = { fa,fb,fc };
    return { std::move(x), std::move(f) };
}