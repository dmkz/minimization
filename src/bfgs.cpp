#include "bfgs.hpp"


// Авторы: Михнев Денис (реализация), Грошева Екатерина (теория)

IterationData bfgs(Function f, Vector start_point, const StopCondition& stop_condition) {
// f - указатель на целевую функцию
// start_point - начальное приближение
// stop_condition - критерий остановы
// Результат работы метода будет лежать в структуре данных о последней итерации

// Инициализируем начальной точкой структуру данных итерации:
    IterationData iter_data;
    iter_data.x_curr = start_point;
    iter_data.f_curr = f(start_point);
    iter_data.iter_counter = 0;
    iter_data.method_title = "BFGS";      
    
    int n = (int)start_point.size();
    Matrix B(n, Vector(n));
    Vector t = grad(f, start_point);
    for(int i = 0; i < n; ++i)
        B[i][i] = 1;
    Real alpha;
    for(int i = 0; i < n; ++i)
        t[i] = -t[i];
    alpha = search_alpha_bfgs(f, start_point, t, 100);
	if (alpha == -1) {
	    return iter_data;
    }
    Vector x_prv = start_point;
    Vector x_cur = start_point;
    for(int i = 0; i < n; ++i)
        x_cur[i] += alpha*t[i];
    B = hes_upd_bfgs(f, B, x_cur, x_prv);   
    Vector p(n);
    Vector cur_grad(n);
    cur_grad = grad(f, x_cur);
    do {
        Vector p(n);
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                p[i] -= B[i][j] * cur_grad[j];
            }
        }
        alpha = search_alpha_bfgs(f, x_cur, p, 100);
		if (alpha == -1) {
		    break;
        }
        x_prv = x_cur;
        x_cur += alpha * p;
        
        iter_data.next(x_cur, f(x_cur));
        
        B = hes_upd_bfgs(f, B, x_cur, x_prv);
        cur_grad = grad(f, x_cur);
        
    } while (!stop_condition(iter_data));
    return iter_data;
}

Matrix out_pr_bfgs(Vector& x, Vector& y) {
    int n = (int)x.size(); 
    Matrix res(n, Vector(n));
    for(int i = 0; i < n; ++i)
        for(int j = 0; j < n; ++j)
            res[i][j] = x[i]*y[j];
    return res;
}

Matrix hes_upd_bfgs(Function f, Matrix& B, Vector& x_cur, Vector& x_prv) {
    int n = (int)x_cur.size();
    Vector s(n);
    for(int i = 0; i < n; ++i) {
        s[i] = x_cur[i] - x_prv[i];
    }
    Vector y(n), y1(n), y2(n);
    y1 = grad(f, x_cur);
    y2 = grad(f, x_prv);
    for(int i = 0; i < n; ++i) {
        y[i] = y1[i] - y2[i];
    }
    Real ro = 1.0/dot(y, s);
    Matrix C;
    C = out_pr_bfgs(s, y);
    for(int i = 0; i < n; ++i) {
        for(int j = 0; j < n; ++j) {
            C[i][j] = -ro*C[i][j];
            if(i == j) C[i][j] += 1;
        }
    }
    Matrix A(n, Vector(n));
    for(int i = 0; i < n; ++i) {
        for(int j = 0; j < n; ++j) {
            for(int k = 0; k < n; ++k) {
                A[i][j] += C[i][k]*B[k][j];
            }
        }
    }
    C = out_pr_bfgs(y, s);
    for(int i = 0; i < n; ++i) {
        for(int j = 0; j < n; ++j) {
            C[i][j] = -ro*C[i][j];
            if(i == j) C[i][j] += 1;
        }
    }
    Matrix res(n, Vector(n));
    for(int i = 0; i < n; ++i) {
        for(int j = 0; j < n; ++j) {
            for(int k = 0; k < n; ++k) {
                res[i][j] += A[i][k]*C[k][j];
            }
        }
    }
    C = out_pr_bfgs(s, s);
    for(int i = 0; i < n; ++i) {
        for(int j = 0; j < n; ++j) {
            C[i][j] = ro*C[i][j];
            res[i][j] += C[i][j];
        }
    }
    return res;
}

Real search_alpha_bfgs(Function f, Vector& x, Vector& p, int iter_limit) {
    Real alpha0 = 1.0;
    Vector x_cur = x;
    for(size_t i = 0; i < x_cur.size(); ++i) {
        x_cur[i] += alpha0*p[i];
    }
    Real phi_a0 = f(x_cur);
    Real derphi0 = dot(grad(f, x), p);
    Real phi0 = f(x);
    if(phi_a0 < phi0 + 0.0001*alpha0*derphi0 + COMPARE_EPS) {
        return alpha0;
    }
    Real alpha1 = -(derphi0) * alpha0*alpha0 / 2.0 / (phi_a0 - phi0 - derphi0 * alpha0);
    for(size_t i = 0; i < x_cur.size(); ++i) {
        x_cur[i] = x[i] + alpha1*p[i];
    }
    Real phi_a1 = f(x_cur);
    if (phi_a1 <= phi0 + 0.0001*alpha1*derphi0 + COMPARE_EPS) {
        return alpha1;
    }
    for (int iter = 0; iter < iter_limit && alpha1 > 0; ++iter) {
        Real factor = alpha0 * alpha0 * alpha1 * alpha1 * (alpha1-alpha0);
        Real a = alpha0 * alpha0 * (phi_a1 - phi0 - derphi0*alpha1) - 
                   alpha1 * alpha1 * (phi_a0 - phi0 - derphi0*alpha0);
        a = a / factor;
        Real b = pow(-alpha0, 3) * (phi_a1 - phi0 - derphi0*alpha1) +
                   pow(alpha1, 3) * (phi_a0 - phi0 - derphi0*alpha0);
        b = b / factor;
        Real alpha2 = (-b + sqrt(fabs(b*b - 3 * a * derphi0))) / (3.0*a);
        for(size_t i = 0; i < x_cur.size(); ++i) {
            x_cur[i] = x[i] + alpha2 * p[i];
        }
        Real phi_a2 = f(x_cur);
        if (phi_a2 <= phi0 + 0.0001*alpha2*derphi0+COMPARE_EPS) {
            return alpha2;
        }
        if ((alpha1 - alpha2) + COMPARE_EPS > alpha1 / 2.0 || (1 - alpha2/alpha1) < 0.96 + COMPARE_EPS) {
            alpha2 = alpha1 / 2.0;
        }
        alpha0 = alpha1;
        alpha1 = alpha2;
        phi_a0 = phi_a1;
        phi_a1 = phi_a2;
    }
    return -1;
}