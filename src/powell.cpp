#include "powell.hpp"
#include "math.hpp"

// Авторы: Бабичев Денис (теория), Бессонов Трофим (тестирование, теория), Данилов Алексей (реализация), Киселев Николай (реализация)

std::pair<Vector, int> powell(Function func, Vector p, int iter_limit)
{
	// const int ITMAX = 200;
	const ld TINY = 1.0e-25;
	int i, j, ibig;
	ld del, fp, fptt, t;
	const ld ftol = 1.0e-6;
	int n = p.size();
	std::pair<Vector, Vector> ans;
    	
    std::vector<std::vector<ld>> xi(n, std::vector<ld>(n));
    for (int i = 0; i < n; ++i) {
        xi[i][i] = 1.0;
    }

	Vector pt(n), ptt(n), xit(n);
	ld fret = func(p);
	for (j = 0; j<n; j++) pt[j] = p[j];
	int iter;
    for (iter = 0; iter < iter_limit; ++iter) {
		fp = fret;
		ibig = 0;
		del = 0.0;
		for (i = 0; i<n; i++) {
			for (j = 0; j<n; j++) xit[j] = xi[j][i];
			fptt = fret;
			ans = linmin(p, xit, func);
			p = ans.first;
			xit = ans.second;
			fret = func(p);
			if (fptt - fret > del) {
				del = fptt - fret;
				ibig = i + 1;
			}
		}
		if (2.0*(fp - fret) <= ftol * (std::abs(fp) + std::abs(fret)) + TINY) {
			return {p, iter};
		}
		for (j = 0; j<n; j++) {
			ptt[j] = 2.0*p[j] - pt[j];
			xit[j] = p[j] - pt[j];
			pt[j] = p[j];
		}
		fptt = func(ptt);
		if (fptt < fp) {
			t = 2.0*(fp - 2.0*fret + fptt)*SQR(fp - fret - del) - del * SQR(fp - fptt);
			if (t < 0.0) {
				ans = linmin(p, xit, func);
				p = ans.first;
				xit = ans.second;
				fret = func(p);
				for (j = 0; j<n; j++) {
					xi[j][ibig - 1] = xi[j][n - 1];
					xi[j][n - 1] = xit[j];
				}
			}
		}
	}
    return {p, iter};
}


ld f1dim(const ld x, int ncom, Vector *pcom_p, Vector *xicom_p, Function nrfunc)
{
	Vector xt(ncom);
	Vector &pcom = *pcom_p, &xicom = *xicom_p;
	for (int j = 0; j<ncom; j++)
		xt[j] = pcom[j] + x * xicom[j];
	ld ans = nrfunc(xt);
	return ans;
}

std::pair<Vector, Vector> linmin(Vector &pInit, Vector &xiInit, Function func)
{
	int ncom;
	Vector *pcom_p, *xicom_p;

	int j;
	const ld TOL = 1.0e-8;
	ld xx, xmin, bx, ax;
    // ld fx, fb, fa; // не используется
	Vector p = pInit; Vector xi = xiInit;
	int n = p.size();
	ncom = n;
	pcom_p = new Vector(n);
	xicom_p = new Vector(n);
	
	Vector &pcom = *pcom_p, &xicom = *xicom_p;
	for (j = 0; j<n; j++) {
		pcom[j] = p[j];
		xicom[j] = xi[j];
	}
	ax = 0.0;
	xx = 1.0;
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
	delete xicom_p;
	delete pcom_p;
	return { p,xi };
}




ld brent(const ld ax, const ld bx, const ld cx, ld f(const ld, int, Vector *, Vector *, Function),
	const ld tol, int ncom, Vector *pcom_p, Vector *xicom_p, Function func)
{
	const int ITMAX = 100;
	const ld CGOLD = 0.3819660;
	const ld ZEPS = 1.0e-3;
	int iter;
	ld a, b, d = 0.0, etemp, fu, fv, fw, fx;
	ld p, q, r, tol1, tol2, u, v, w, x, xm;
	ld e = 0.0;

	a = (ax < cx ? ax : cx);
	b = (ax > cx ? ax : cx);
	x = w = v = bx;
	fw = fv = fx = f(x,ncom,pcom_p,xicom_p,func);
	for (iter = 0; iter<ITMAX; iter++) {
		xm = 0.5*(a + b);
		tol2 = 2.0*(tol1 = tol * std::abs(x) + ZEPS);
		if (std::abs(x - xm) <= (tol2 - 0.5*(b - a))) {
			return x;
		}
		if (std::abs(e) > tol1) {
			r = (x - w)*(fx - fv);
			q = (x - v)*(fx - fw);
			p = (x - v)*q - (x - w)*r;
			q = 2.0*(q - r);
			if (q > 0.0) p = -p;
			q = std::abs(q);
			etemp = e;
			e = d;
			if (std::abs(p) >= std::abs(0.5*q*etemp) || p <= q * (a - x) || p >= q * (b - x))
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


std::pair<Vector, Vector> mnbrak(ld axInit, ld bxInit, ld func(const ld, int, Vector *, Vector *, Function), int ncom, Vector *pcom_p, Vector *xicom_p, Function nrfunc)
{
	const ld GOLD = 1.618034, GLIMIT = 100.0, TINY = 1.0e-20;
	ld ulim, u, r, q, fu;
	ld fa, fb, fc, cx;
	ld ax = axInit;
	ld bx = bxInit;
	fa = func(ax, ncom, pcom_p, xicom_p, nrfunc);
	fb = func(bx, ncom, pcom_p, xicom_p, nrfunc);
	if (fb > fa) {
		SWAP(ax, bx);
		SWAP(fb, fa);
	}
	cx = bx + GOLD * (bx - ax);
	fc = func(cx, ncom, pcom_p, xicom_p, nrfunc);
	while (fb > fc) {
		r = (bx - ax)*(fb - fc);
		q = (bx - cx)*(fb - fa);
		u = bx - ((bx - cx)*q - (bx - ax)*r) /
			(2.0*SIGN(MAX(std::abs(q - r), TINY), q - r));
		ulim = bx + GLIMIT * (cx - bx);
		if ((bx - u)*(u - cx) > 0.0) {
			fu = func(u, ncom, pcom_p, xicom_p, nrfunc);
			if (fu < fc) {
				ax = bx;
				bx = u;
				fa = fb;
				fb = fu;
				Vector x = { ax,bx,cx }; Vector f = { fa,fb,fc };
				return { x,f };
			}
			else if (fu > fb) {
				cx = u;
				fc = fu;
				Vector x = { ax,bx,cx }; Vector f = { fa,fb,fc };
				return { x,f };
			}
			u = cx + GOLD * (cx - bx);
			fu = func(u, ncom, pcom_p, xicom_p, nrfunc);
		}
		else if ((cx - u)*(u - ulim) > 0.0) {
			fu = func(u, ncom, pcom_p, xicom_p, nrfunc);
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
		else if ((u - ulim)*(ulim - cx) >= 0.0) {
			u = ulim;
			fu = func(u, ncom, pcom_p, xicom_p, nrfunc);
		}
		else {
			u = cx + GOLD * (cx - bx);
			fu = func(u, ncom, pcom_p, xicom_p, nrfunc);
		}
		//shft3(ax, bx, cx, u);
		ax = bx;
		bx = cx;
		cx = u;
		//shft3(fa, fb, fc, fu);
		fa = fb;
		fb = fc;
		fc = fu;
	}
	Vector x = { ax,bx,cx }; Vector f = { fa,fb,fc };
	return { x,f };
}