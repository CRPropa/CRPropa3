#ifndef MPC_EXPLICITRUNGEKUTTA_H_
#define MPC_EXPLICITRUNGEKUTTA_H_

#include <vector>
#include <stdexcept>
#include <iostream>
#include <sstream>

namespace mpc {

template<class Y>
class ExplicitRungeKutta {
public:
	std::vector<double> c, b, bs, a;
	unsigned int s;
	class F {
	public:
		virtual ~F() {
		}
		virtual Y operator()(double t, const Y &v) = 0;
	};
	void step(double t, const Y &initial, Y &result, double stepsize,
			F &f) const;
	void step(double t, const Y &initial, Y &result, Y &error, double stepsize,
			F &f) const;
	void loadCashKarp();
};

template<class Y>
inline void ExplicitRungeKutta<Y>::step(double t, const Y &y, Y &result,
		double stepsize, F &f) const {
	std::vector<Y> k;
	k.reserve(s);

	Y r = y;

	// calculate the sum of b_i * k_i
	for (unsigned int i = 0; i < s; i++) {
		// first parameter of f
		double t_n = t + c[i] * stepsize;

		// second parameter of f
		Y y_n = y;
		for (unsigned int j = 0; j < i; j++) {
			y_n += k[j] * a[i * s + j] * stepsize;
		}

		// update k_i
		k[i] = f(t_n, y_n);

		r += k[i] * b[i] * stepsize;
	}

	result = r;
}

template<class Y>
inline void ExplicitRungeKutta<Y>::step(double t, const Y &y, Y &result,
		Y &error, double stepsize, F &f) const {
	std::vector<Y> k;
	k.reserve(s);

	Y r = y;
	error = 0;

	// calculate the sum of b_i * k_i
	for (unsigned int i = 0; i < s; i++) {
		// first parameter of f
		double t_n = t + c[i] * stepsize;

		// second parameter of f
		Y y_n = y;
		for (unsigned int j = 0; j < i; j++) {
			y_n += k[j] * a[i * s + j] * stepsize;
		}

		// update k_i
		k[i] = f(t_n, y_n);

		r += k[i] * b[i] * stepsize;
		error += k[i] * (b[i] - bs[i]) * stepsize;
	}

	result = r;
}

// Cash-Karp coefficients
const double cash_karp_a[] = { 0., 0., 0., 0., 0., 0., 1. / 5., 0., 0., 0., 0.,
		0., 3. / 40., 9. / 40., 0., 0., 0., 0., 3. / 10., -9. / 10., 6. / 5.,
		0., 0., 0., -11. / 54., 5. / 2., -70. / 27., 35. / 27., 0., 0., 1631.
				/ 55296., 175. / 512., 575. / 13824., 44275. / 110592., 253.
				/ 4096., 0. };

const double cash_karp_b[] = { 37. / 378., 0, 250. / 621., 125. / 594., 0., 512.
		/ 1771. };
const double cash_karp_bs[] = { 2825. / 27648., 0., 18575. / 48384., 13525.
		/ 55296., 277. / 14336., 1. / 4. };
const double cash_karp_c[] = { 0., 1. / 5., 3. / 10., 3. / 5., 1., 7. / 8. };

template<class Y>
inline void ExplicitRungeKutta<Y>::loadCashKarp() {
	s = 6;
	a.assign(cash_karp_a, cash_karp_a + 36);
	b.assign(cash_karp_b, cash_karp_b + 6);
	bs.assign(cash_karp_bs, cash_karp_bs + 6);
	c.assign(cash_karp_c, cash_karp_c + 6);
}

} // namespace mpc

#endif /* MPC_EXPLICITRUNGEKUTTA_H_ */
