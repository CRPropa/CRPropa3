#ifndef _MPC_VECTOR3_H_
#define _MPC_VECTOR3_H_

#include <iostream>
#include <cmath>
#include <vector>
#include <limits>

namespace mpc {

template<typename T>
class Vector3 {
public:
	T x, y, z;

	Vector3() :
			x(0), y(0), z(0) {
	}

	template<typename U>
	Vector3(const Vector3<U> &v) :
			x(v.x), y(v.y), z(v.z) {
	}

	explicit Vector3(const double *v) :
			x(v[0]), y(v[1]), z(v[2]) {
	}

	explicit Vector3(const float *v) :
			x(v[0]), y(v[1]), z(v[2]) {
	}

	explicit Vector3(const T &X, const T &Y, const T &Z) :
			x(X), y(Y), z(Z) {
	}

	explicit Vector3(T t) :
			x(t), y(t), z(t) {
	}

	bool operator <(const Vector3<T> &v) const {
		if (x > v.x)
			return false;
		else if (x < v.x)
			return true;
		if (y > v.y)
			return false;
		else if (y < v.y)
			return true;
		if (z >= v.z)
			return false;
		else
			return true;
	}

	bool operator ==(const Vector3<T> &v) const {
		if (x != v.x)
			return false;
		if (y != v.y)
			return false;
		if (z != v.z)
			return false;
		return true;
	}

	Vector3<T> operator -(const Vector3<T> &v) const {
		return Vector3(x - v.x, y - v.y, z - v.z);
	}

	Vector3<T> operator +(const Vector3<T> &v) const {
		return Vector3(x + v.x, y + v.y, z + v.z);
	}

	Vector3<T> operator *(const Vector3<T> &v) const {
		return Vector3(x * v.x, y * v.y, z * v.z);
	}

	Vector3<T> operator *(const T &v) const {
		return Vector3(x * v, y * v, z * v);
	}

	Vector3<T> operator /(const T &f) const {
		return Vector3(x / f, y / f, z / f);
	}

	Vector3<T> &operator /=(const T &f) {
		x /= f;
		y /= f;
		z /= f;
		return *this;
	}

	Vector3<T> operator %(const T &f) const {
		return Vector3(x % f, y % f, z % f);
	}

	Vector3<T> &operator %=(const T &f) {
		x %= f;
		y %= f;
		z %= f;
		return *this;
	}

	Vector3<T> &operator *=(const T &f) {
		x *= f;
		y *= f;
		z *= f;
		return *this;
	}

	Vector3<T> &operator +=(const Vector3<T> &v) {
		x += v.x;
		y += v.y;
		z += v.z;
		return *this;
	}

	Vector3<T> &operator -=(const Vector3<T> &v) {
		x -= v.x;
		y -= v.y;
		z -= v.z;
		return *this;
	}

	Vector3<T> &operator =(const Vector3<T> &v) {
		x = v.x;
		y = v.y;
		z = v.z;
		return *this;
	}

	void clamp(T lower, T upper) {
		x = std::max(lower, std::min(x, upper));
		y = std::max(lower, std::min(y, upper));
		z = std::max(lower, std::min(z, upper));
	}

	T mag() const {
		return std::sqrt(x * x + y * y + z * z);
	}

	T mag2() {
		return (x * x + y * y + z * z);
	}

	T dot(const Vector3<T> & p) const {
		return x * p.x + y * p.y + z * p.z;
	}

	Vector3<T> cross(const Vector3<T> & p) const {
		return Vector3<T>(y * p.z - p.y * z, z * p.x - p.z * x,
				x * p.y - p.x * y);
	}

	Vector3<T> floor() const {
		return Vector3<T>(std::floor(x), std::floor(y), std::floor(z));
	}

	Vector3<T> ceil() const {
		return Vector3<T>(std::ceil(x), std::ceil(y), std::ceil(z));
	}

	T phi() const {
		T eps = std::numeric_limits<T>::min();
		if (fabs(x) < eps && fabs(y) < eps)
			return 0.0;
		else
			return std::atan2(y, x);
	}

	T theta() const {
		T eps = std::numeric_limits<T>::min();
		if (fabs(x) < eps && fabs(y) < eps && fabs(z) < eps)
			return 0.0;
		else
			return atan2((T)sqrt(x * x + y * y), z);
	}

	T angleTo(const Vector3<T> &v) const {
		T cosdistance = this->dot(v) / this->mag() / v.mag();
		// In some directions cosdistance is > 1 on some compilers
		// This ensures that the correct result is returned
		if (cosdistance >= 1.)
			return 0;
		if (cosdistance <= -1.)
			return M_PI;
		else
			return acos(cosdistance);
	}

	bool isParallelTo(const Vector3<T> &v, T maxAngle) const {
		T angle = this->angleTo(v);
		return angle < maxAngle;
	}

	T distanceTo(const Vector3<T> &point) const {
		Vector3<T> d = *this - point;
		return d.mag();
	}

	void set(const T x, const T y, const T z) {
		this->x = x;
		this->y = y;
		this->z = z;
	}

	void setLower(const Vector3<T> &v) {
		if (v.x < x)
			x = v.x;
		if (v.y < y)
			y = v.y;
		if (v.z < z)
			z = v.z;
	}

	void setUpper(const Vector3<T> &v) {
		if (v.x > x)
			x = v.x;
		if (v.y > y)
			y = v.y;
		if (v.z > z)
			z = v.z;
	}
};

template<>
inline Vector3<float> Vector3<float>::operator %(const float &f) const {
	return Vector3<float>(fmod(x, f), fmod(y, f), fmod(z, f));
}

template<>
inline Vector3<float> &Vector3<float>::operator %=(const float &f) {
	x = fmod(x, f);
	y = fmod(y, f);
	z = fmod(z, f);
	return *this;
}

template<>
inline Vector3<double> Vector3<double>::operator %(const double &f) const {
	return Vector3<double>(fmod(x, f), fmod(y, f), fmod(z, f));
}

template<>
inline Vector3<double> &Vector3<double>::operator %=(const double &f) {
	x = fmod(x, f);
	y = fmod(y, f);
	z = fmod(z, f);
	return *this;
}

template<typename T>
inline std::ostream &operator <<(std::ostream &out, const Vector3<T> &v) {
	out << v.x << " " << v.y << " " << v.z;
	return out;
}

template<typename T>
inline std::istream &operator >>(std::istream &in, Vector3<T> &v) {
	in >> v.x >> v.y >> v.z;
	return in;
}

template<typename T>
Vector3<T> operator *(T f, const Vector3<T> &v) {
	return Vector3<T>(v.x * f, v.y * f, v.z * f);
}

typedef Vector3<double> Vector3d;
typedef Vector3<float> Vector3f;

} // namespace mpc

#endif /* _MPC_VECTOR3_H_ */
