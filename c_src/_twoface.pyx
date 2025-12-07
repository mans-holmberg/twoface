import numpy as np

cdef extern from "twoface.h":
	void integratetransit(int, int, int, double*, double*, double*, double, double*, double*, double*, double*, double*, double*, double*, double*)
	void circleangle(double*, double, double, int, double*)
	void integratetransit_asymmetric(int, int, int, double*, double*, double*, double, double, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*);
	void sky(int, double*, double, double, double, double, double, double, double*, double*, double*, double*);


cpdef _sky(int n, double[:] t_array, double tc, double per, double a, double inc, double ecc, double omega, double[:] z, double[:] x, double[:] y, double[:] psi):
	sky(n, &t_array[0], tc, per, a, inc, ecc, omega, &z[0], &x[0], &y[0], &psi[0])


cpdef _integratetransit(int m, int n, int k, double[:] planetx, double[:] planety, double[:] z, double p, double[:] r, double[:] f, double[:] spotx, double[:] spoty, double[:] spotradius, double[:] spotcontrast, double[:, :] planetangle, double[:] answer):
	integratetransit(m, n, k, &planetx[0], &planety[0], &z[0], p, &r[0], &f[0], &spotx[0], &spoty[0], &spotradius[0], &spotcontrast[0], &planetangle[0, 0], &answer[0])


cpdef _integratetransit_asymmetric(int m, int n, int k, double[:] r, double[:] f, double[:] planetx, double[:] planety, double[:] z, double[:] phi, double p1, double p2, double[:] spotx, double[:] spoty, double[:] spotradius, double[:] spotcontrast, double[:] spotcenterdistance, double[:] bounds, double[:] answer):
	integratetransit_asymmetric(m, n, k, &planetx[0], &planety[0], &z[0], p1, p2, &r[0], &f[0], &spotx[0], &spoty[0], &spotradius[0], &spotcontrast[0], &phi[0], &answer[0], &spotcenterdistance[0], &bounds[0])


cpdef _circleangle(double[:] r, double p, double[:] z , double[:, :] answer):
	cdef int n = r.shape[0]
	cdef int m = z.shape[0]
	cdef int i = 0
	for i in range(m):
		circleangle(&r[0], p, z[i], n, &answer[i][0])
