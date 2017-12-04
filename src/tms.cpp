#include <iostream>
#include <cmath>
#include <fstream>
#include <vector>
#include <algorithm>
#include <iomanip>

using namespace std;

/*
maxdim - макс размерность пространства
maxfig - максимальное количество чисел с основанием Q , которое мы можем обработать
nbits - число бит в слове (сейчас 31, для поддержания платформонезависимости)
c - хранит числа Нидеррайтера (I,J,R)
count - индекс текущего элемента в C
d - значения D(I,J)
nextq - номер следующего элемента
qpow - хранит степени Q
dimen - размерность генерируемой последовательности
nfigs - количество  чисел по основанию Q, которые мы используем по факту
recip - переменная, хранящая 1 / (Q ^ NFIGS)

maxq - максимальное значение Q
maxdeg - наибольшая степень полинома (многочлена)
deg - текущая степень многочлена
p, q - характеристики поля
add - результат сложения многочленов
mul - результат умножения
sub - результат деления
*/

// COMM
const int maxdim = 20;
const int maxfig = 20;
const int nbits = 31;

int c[maxdim][maxfig][maxfig/*!*/]; // /*!*/ - индексация как в си, т.е. C(MAXDIM, MAXFIG, 0:MAXFIG-1)
int count[maxfig/*!*/], d[maxdim][maxfig];
int nextq[maxdim], qpow[maxfig];
int dimen, nfigs;
double recip;

// COMM2
int cj[maxdim][nbits/*!*/];
int count2;
int dimen2;
int nextq2[maxdim];

// FIELD
const int maxq = 50;
const int maxdeg = 50, deg = 0;
int p, q, add[maxq/*!*/][maxq/*!*/], mul[maxq/*!*/][maxq/*!*/], sub[maxq/*!*/][maxq/*!*/];


vector<vector<double>> vResult;
vector<double> vector_;
int ResultCounter = 0;

struct point {
    int point;
    double value;
};

// эта подпрограмма инициализирует таблицы сложений, умножений, вычитаний
void SETFLD(int qin) {
    cout << "SETFLD" << endl;
    int i, j;
    if (qin <= 1 || qin > maxq) cout << "SETFLD: Bad value of Q";
    q = qin;
    p = q;
    if (p == 0) cout << "SETFLD: There is no field of order" << q;

    for (i = 0; i < q; ++i)
        for (j = 0; j < q; ++j) {
            add[i][j] = (i + j) % p;
            mul[i][j] = (i * j) % p;
        }

    for (i = 0; i < q; ++i)
        for (j = 0; j < q; ++j) sub[add[i][j]][i] = j;
}

//Программа используется для умножения полиномов
void PLYMUL(int *pa, int *pb, int *pc) {
    int i, j, dega, degb, degc, term;
    int pt[maxdeg + 2];
    dega = pa[deg];
    degb = pb[deg];
    if (dega == -1 || degb == -1) degc = -1;
    else degc = dega + degb;
    if (degc > maxdeg)
        cout << "PLYMUL: Degree of product exceeds MAXDEG" << endl;
    for (i = 0; i <= degc; ++i) {
        term = 0;
        for (j = max(0, i - dega); j <= min(degb, i); ++j)
            term = add[term][mul[pa[i - j + 1]][pb[j + 1]]];
        pt[i + 1] = term;
    }
    pc[deg] = degc;
    for (i = 1; i <= degc + 1; ++i) pc[i] = pt[i];
    for (i = degc + 2; i < maxdeg + 2; ++i) pc[i] = 0;
}

/*
  C помощью этой программы вычисляем значения констант V(J,R)

  px соответствующего неприводимого полинома для  текущего измерения
  На входе в подпрограмму, степень B  определяет параметр J, как deg(B) = Е*(J-1), где Е - степень px.
  MAXV дает размерность массива V.

  результат выполнения - B умноженный на px, его степень на данный момент Е*J.
  V содержит требуемыые значения.
  подпрограмма PLYMUL используется для умножения полиномов
*/

void CALCV(int *px, int *b, int *v, int maxv) {
    int h[maxdeg + 2];
    int bigm = 0, m = 0, kj, term;
    int arbit = 1, nonzer = 1;

    for (int i = 0; i < b[deg] + 2; ++i)
        h[i] = b[i];
    bigm = h[deg];
    PLYMUL(px, b, b);

    m = b[deg];

    kj = bigm;
    for (int r = 0; r < kj; ++r)
        v[r] = 0;
    v[kj] = 1;

    if (kj < bigm) {
        // weak branch
        term = sub[0][h[kj + 1]];
        for (int r = kj - 1; r < bigm - 1; ++r) {
            v[r] = arbit;
            term = sub[term][mul[h[r + 1]][v[r]]];
        }
        v[bigm] = add[nonzer][term];
        for (int r = bigm + 1; r < m - 1; ++r)
            v[r] = arbit;
    } else
        for (int r = kj + 1; r <= m - 1; ++r)
            v[r] = arbit;
    for (int r = 0; r <= maxv - m; ++r) {
        term = 0;
        for (int i = 0; i <= m - 1; ++i)
            term = sub[term][mul[b[i + 1]][v[r + i]]];
        v[r + m] = term;
    }
}

/*
C помощью этой программы вычисляем значения константы  C(I,J,R).
для каждого значения I, мы сначала рассчитываем все соответствующие значения  C : они хранятся в массиве CI.
Все  значения C  0 или 1. Далее мы помещаем значения в массив  СJ, таким образом,  CJ(I,R) имеет значения С
для указанных значений I и R, и для каждого J от 1 до NBITS
И когда значения CJ(I,R) расчитаны, мы возвращаемсяк этому массиву в вызываемой программе
*/

void CALCC2() {
    int maxe = 5, maxv = nbits + maxe;
    int px[maxdeg + 2], b[maxdeg + 2];
    int v[maxv + 1], ci[nbits][nbits];
    int e, i, j, r, u, term;
    int irred[maxdim][maxe + 2];
    {
        irred[0][0] = 1;
        irred[0][1] = 0;
        irred[0][2] = 1;

        irred[1][0] = 1;
        irred[1][1] = 1;
        irred[1][2] = 1;

        irred[2][0] = 2;
        irred[2][1] = 1;
        irred[2][2] = 1;
        irred[2][3] = 1;

        irred[3][0] = 3;
        irred[3][1] = 1;
        irred[3][2] = 1;
        irred[3][3] = 0;
        irred[3][4] = 1;

        irred[4][0] = 3;
        irred[4][1] = 1;
        irred[4][2] = 0;
        irred[4][3] = 1;
        irred[4][4] = 1;

        irred[5][0] = 4;
        irred[5][1] = 1;
        irred[5][2] = 1;
        irred[5][3] = 0;
        irred[5][4] = 0;
        irred[5][5] = 1;

        irred[6][0] = 4;
        irred[6][1] = 1;
        irred[6][2] = 0;
        irred[6][3] = 0;
        irred[6][4] = 1;
        irred[6][5] = 1;

        irred[7][0] = 4;
        irred[7][1] = 1;
        irred[7][2] = 1;
        irred[7][3] = 1;
        irred[7][4] = 1;
        irred[7][5] = 1;

        irred[8][0] = 5;
        irred[8][1] = 1;
        irred[8][2] = 0;
        irred[8][3] = 1;
        irred[8][4] = 0;
        irred[8][5] = 0;
        irred[8][6] = 1;

        irred[9][0] = 5;
        irred[9][1] = 1;
        irred[9][2] = 0;
        irred[9][3] = 0;
        irred[9][4] = 1;
        irred[9][5] = 0;
        irred[9][6] = 1;

        irred[10][0] = 5;
        irred[10][1] = 1;
        irred[10][2] = 1;
        irred[10][3] = 1;
        irred[10][4] = 1;
        irred[10][5] = 0;
        irred[10][6] = 1;

        irred[11][0] = 5;
        irred[11][1] = 1;
        irred[11][2] = 1;
        irred[11][3] = 1;
        irred[11][4] = 0;
        irred[11][5] = 1;
        irred[11][6] = 1;

        irred[12][0] = 5;
        irred[12][1] = 1;
        irred[12][2] = 1;
        irred[12][3] = 0;
        irred[12][4] = 1;
        irred[12][5] = 1;
        irred[12][6] = 1;

        irred[13][0] = 5;
        irred[13][1] = 1;
        irred[13][2] = 0;
        irred[13][3] = 1;
        irred[13][4] = 1;
        irred[13][5] = 1;
        irred[13][6] = 1;

        irred[14][0] = 6;
        irred[14][1] = 1;
        irred[14][2] = 1;
        irred[14][3] = 0;
        irred[14][4] = 0;
        irred[14][5] = 0;
        irred[14][6] = 0;
        irred[14][7] = 1;

        irred[15][0] = 6;
        irred[15][1] = 1;
        irred[15][2] = 0;
        irred[15][3] = 0;
        irred[15][4] = 1;
        irred[15][5] = 0;
        irred[15][6] = 0;
        irred[15][7] = 1;

        irred[16][0] = 6;
        irred[16][1] = 1;
        irred[16][2] = 1;
        irred[16][3] = 1;
        irred[16][4] = 0;
        irred[16][5] = 1;
        irred[16][6] = 0;
        irred[16][7] = 1;

        irred[17][0] = 6;
        irred[17][1] = 1;
        irred[17][2] = 1;
        irred[17][3] = 0;
        irred[17][4] = 1;
        irred[17][5] = 1;
        irred[17][6] = 0;
        irred[17][7] = 1;

        irred[18][0] = 6;
        irred[18][1] = 1;
        irred[18][2] = 0;
        irred[18][3] = 0;
        irred[18][4] = 0;
        irred[18][5] = 0;
        irred[18][6] = 1;
        irred[18][7] = 1;

        irred[19][0] = 6;
        irred[19][1] = 1;
        irred[19][2] = 1;
        irred[19][3] = 1;
        irred[19][4] = 0;
        irred[19][5] = 0;
        irred[19][6] = 1;
        irred[19][7] = 1;
    }
    SETFLD(2);

    for (i = 0; i < dimen; ++i) {
        e = irred[i][deg];
        b[deg] = 0;
        b[1] = 1;
        u = 0;
        for (j = 0; j < e + 2; ++j)
            px[j] = irred[i][j];

        for (j = 0; j <= nbits - 1; ++j) {
            if (u == 0) CALCV(px, b, v, maxv);
            for (r = 0; r <= nbits - 1; ++r)
                ci[j][r] = v[r + u];
            ++u;
            if (u == e) u = 0;
        }
        for (r = 0; r <= nbits - 1; ++r) {
            term = 0;
            for (j = 1; j <= nbits; ++j)
                term = (term << 1) + ci[j - 1][r];
            cj[i][r] = term;
        }
    }
}

/*
Эта подпрограмма вычисляет значения нидеррейтера C(I,J,R) , вызывая CALCC2 и выполняет другие необходимые инициализации перед вызовом GOLO2
dim - измерение последовательностей, которые должны быть созданы.
skip - количество значений, которые должны быть отброшены в начале последовательности
*/

void INLO2(int dim, int skip) {
    int r, gray;
    dimen2 = dim;
    if (dimen2 <= 0 || dimen2 > maxdim) {
        cout << "INLO2 : Bad dimension";
        return;
    }

    CALCC2();
    gray = skip xor skip >> 1;
    for (int i = 1; i <= dimen2; ++i) nextq2[i - 1] = 0;
    r = 0;
    while (gray != 0) {
        if (gray % 2 != 0)
            for (int i = 1; i <= dimen2; ++i)
                nextq2[i - 1] = nextq2[i - 1] xor cj[i - 1][r];
        gray >>= 1; //делим на 2
        r++;
    }
    count2 = skip;
}

/*
Эта функция создает новый квази-случайный вектор на каждый вызов программы
bounds - границы области, в которой генерируются точки.
*/

void GOLO2(double *quasi, vector<vector<double>> &bounds) {
    int r;
    recip = pow(2, -nbits);
    for (int i = 0; i < dimen2; ++i) {
        quasi[i] = nextq2[i] * recip;
        vResult[ResultCounter][i] = bounds[i][0] + (bounds[i][1] - bounds[i][0]) * quasi[i];
    }
    ++ResultCounter;

    r = 0;
    int i = count2;
    while (i % 2 != 0) {
        r++;
        i >>= 1;
    }
    if (r >= nbits) {
        cout << "GOLO2 : Too many calls";
        return;
    }
    for (int I = 0; I < dimen2; ++I)
        nextq2[I] = nextq2[I] xor cj[I][r];
    ++count2;
}

//
double func(vector<double> &p) {
    return p[0] * p[0] + p[1] * p[1];
}

//Функция сравнения точек для сортировки
int cmp(const void *a, const void *b) {
    struct point *x = (struct point *) a;
    struct point *y = (struct point *) b;
    return 2 * ((x->value - y->value) > 0) - 1;
}

/*
C помощью этой программы тестируем на точность численное интегрирование
C помощью малых отклонений двоичных последовательностей Нидеррайтера (1988), как это реализовано в INLO2, GOLO2.
Различные возможные проверки интегралов GENIN2 генерирует только  последовательности с основанием 2.
dime_ - размерность пространвства
sqlen_ - колличество точе
m_ - колличество минимальныйх точек для запуска функций минимизации.
bounds - границы области
*/

void
GENIN2(int dimen_, int seqlen_, int m_, double (*fun)(vector<double> &), vector<vector<double>> &bounds) // PROGRAM
{
    dimen = dimen_;
    int seqlen = seqlen_;
    int m = m_;

    // выделяем память для точек. количество точек: seqlen; координат на точку: dimen
    vResult.resize(seqlen);
    vector_.resize(seqlen);
    for (int i = 0; i < seqlen; ++i)
        vResult[i].resize(dimen);

    int skip = 0;

    INLO2(dimen, skip);
    cout << "GENIN2 :  Initialization complete" << endl;

    double quasi[maxdim];

    for (int i = 1; i <= seqlen; ++i)
        GOLO2(quasi, bounds);

    cout << "GENIN2:  iteration ends" << endl;


    //int k; //количество точек в сетке
    ofstream fout("out.txt");
    struct point *values = (struct point *) calloc(seqlen, sizeof(struct point));

    for (int i = 0; i < seqlen; ++i) {
        for (int j = 0; j < dimen; ++j) {
            vector_[j] = vResult[i][j];
            fout << fixed << setprecision(9) << vResult[i][j] << " ";
        }
        fout << endl;
        values[i].value = fun(vector_); //записываем в вектор значение функции в данной точке
        values[i].point = i; //записываем в вектор номер данной точки
    }
    fout.close();

    qsort(values, seqlen, sizeof(struct point),
          cmp); // сортируем вектор значений функции и запоминаем какой точке принадлежит каждое значение

    ofstream fout1("out_min.txt"); // файл, в котором хранятся все m точек, где функция минимальна

    for (int p = 0; p < m; ++p) {
        for (int t = 0; t < dimen; ++t)
            fout1 << vResult[values[p].point][t]
                  << " "; // записываем в файл координаты точек, в которых функция принимает минимальное значение
        fout1 << " " << values[p].value << endl;
    }

    fout1.close();
    cout << "Answer in out_min.txt" << endl;
    free(values);
}

/**
 * Calculate tsm net.
 * Read dimension and bounds from file.
 *
 * Params file example:
 *  5 - This is dimension. Dimension should be less than maxdim, which seems taken RANDOMLY.
 *  0 1 - This and below are bound.
 *  0 1
 *  0 1
 *  0 1
 *  0 1
 *  2^10 = 1024 or 2^15 = 32768 or 2^20 = 1048576 - This is sequence length.
 *                                                  There are no assumptions why this is so.
 *                                                  Length must be strictly positive.
 *  10 - This is dots count.
 *
 */
int main(int argc, char *argv[]) {
    if (argc != 2) {
        cout << "tsm <file>\n\n";
        cout << "Calculate tsm net.\n"
                "Read dimension and bounds from file.\n"
                "\n"
                "Params file example:\n"
                "  5 - This is dimension. Dimension should be less than maxdim.\n"
                "  0 1 - This and below are bound.\n"
                "  0 1\n"
                "  0 1\n"
                "  0 1\n"
                "  0 1\n"
                "  2^10 = 1024 or 2^15 = 32768 or 2^20 = 1048576 - This is sequence length.\n"
                "                                                  Length must be strictly positive.\n"
                "  10 - This is dots count.\n";
        return 1;
    }

    int m, dimension, s;
    vector<vector<double>> bounds;
    fstream inputFile(argv[1]);

    if (!inputFile) {
        cout << "Bad file";
        return 1;
    }

    inputFile >> dimension;
    if (dimension > maxdim) {
        cout << "Dimension may not exceed " << maxdim << endl;
        return 1;
    } else {
        cout << "Dimension: " << dimension << endl;
    }

    bounds.resize(dimension);
    for (int i = 0; i < dimension; ++i) {
        bounds[i].reserve(2);
        inputFile >> bounds[i][0] >> bounds[i][1];
        cout << "Bounds for " << i + 1 << " dimension: " << bounds[i][0] << " " << bounds[i][1] << endl;
    }

    inputFile >> s;
    if (s < 0) {
        cout << "Length must be strictly positive" << endl;
        return 1;
    } else {
        cout << "Sequence length: " << s << endl;
    }

    inputFile >> m;
    cout << "Dots count: " << m << endl;

    double (*f_pointer)(vector<double> &) = &func;

    GENIN2(dimension, s, m, f_pointer, bounds);

    return 0;
}
