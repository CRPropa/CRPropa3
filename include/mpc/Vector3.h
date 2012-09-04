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

	// Provides implicit conversion
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

	void setX(const T X) {
		x = X;
	}

	void setY(const T Y) {
		y = Y;
	}

	void setZ(const T Z) {
		z = Z;
	}

	void setXYZ(const T X, const T Y, const T Z) {
		x = X;
		y = Y;
		z = Z;
	}

	void setRThetaPhi(const T r, const T theta, const T phi) {
		x = r * sin(theta) * cos(phi);
		y = r * sin(theta) * sin(phi);
		z = r * cos(theta);
	}

	T getX() const {
		return x;
	}

	T getY() const {
		return y;
	}

	T getZ() const {
		return z;
	}

	T getPhi() const {
		T eps = std::numeric_limits<T>::min();
		if ((fabs(x) < eps) and (fabs(y) < eps))
			return 0.0;
		else
			return std::atan2(y, x);
	}

	T getTheta() const {
		T eps = std::numeric_limits<T>::min();
		if ((fabs(x) < eps) and (fabs(y) < eps) and (fabs(z) < eps))
			return 0.0;
		else
			return atan2((T) sqrt(x * x + y * y), z);
	}

	T getAngleTo(const Vector3<T> &v) const {
		T cosdistance = this->dot(v) / this->getMag() / v.getMag();
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
		T angle = this->getAngleTo(v);
		return angle < maxAngle;
	}

	T getDistanceTo(const Vector3<T> &point) const {
		Vector3<T> d = *this - point;
		return d.getMag();
	}

	T getMag() const {
		return std::sqrt(x * x + y * y + z * z);
	}

	T getMag2() const {
		return x * x + y * y + z * z;
	}

	Vector3<T> getUnitVector() const {
		return Vector3<T>(x, y, z) / getMag();
	}

	void normalize() {
		*this /= getMag();
	}

	T dot(const Vector3<T> &v) const {
		return x * v.x + y * v.y + z * v.z;
	}

	Vector3<T> cross(const Vector3<T> &v) const {
		return Vector3<T>(y * v.z - v.y * z, z * v.x - v.z * x,
				x * v.y - v.x * y);
	}

	void rotate(const Vector3<T> &axis, T angle) {
		T ux = axis.getX() / axis.getMag();
		T uy = axis.getY() / axis.getMag();
		T uz = axis.getZ() / axis.getMag();
		T c = cos(angle);
		T s = sin(angle);
		Vector3<T> Rx(c + ux * ux * (1 - c), ux * uy * (1 - c) - uz * s,
				ux * uz * (1 - c) + uy * s);
		Vector3<T> Ry(uy * ux * (1 - c) + uz * s, c + uy * uy * (1 - c),
				uy * uz * (1 - c) - ux * s);
		Vector3<T> Rz(uz * ux * (1 - c) - uy * s, uz * uy * (1 - c) + ux * s,
				c + uz * uz * (1 - c));
		this->setXYZ(Rx.dot(*this), Ry.dot(*this), Rz.dot(*this));
	}

	void clamp(T lower, T upper) {
		x = std::max(lower, std::min(x, upper));
		y = std::max(lower, std::min(y, upper));
		z = std::max(lower, std::min(z, upper));
	}

	Vector3<T> abs() const {
		return Vector3<T>(std::abs(x), std::abs(y), std::abs(z));
	}

	Vector3<T> floor() const {
		return Vector3<T>(std::floor(x), std::floor(y), std::floor(z));
	}

	Vector3<T> ceil() const {
		return Vector3<T>(std::ceil(x), std::ceil(y), std::ceil(z));
	}

	T min() {
		return std::min(x, std::min(y, z));
	}

	T max() {
		return std::max(x, std::max(y, z));
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

	Vector3<T> operator +(const Vector3<T> &v) const {
		return Vector3(x + v.x, y + v.y, z + v.z);
	}

	Vector3<T> operator +(const T &f) const {
		return Vector3(x + f, y + f, z + f);
	}

	Vector3<T> operator -(const Vector3<T> &v) const {
		return Vector3(x - v.x, y - v.y, z - v.z);
	}

	Vector3<T> operator -(const T &f) const {
		return Vector3(x - f, y - f, z - f);
	}

	Vector3<T> operator *(const Vector3<T> &v) const {
		return Vector3(x * v.x, y * v.y, z * v.z);
	}

	Vector3<T> operator *(const T &v) const {
		return Vector3(x * v, y * v, z * v);
	}

	Vector3<T> operator /(const Vector3<T> &v) const {
		return Vector3(x / v.x, y / v.y, z / v.z);
	}

	Vector3<T> operator /(const T &f) const {
		return Vector3(x / f, y / f, z / f);
	}

	Vector3<T> operator %(const Vector3<T> &v) const {
		return Vector3(fmod(x, v.x), fmod(y, v.y), fmod(z, v.z));
	}

	Vector3<T> operator %(const T &f) const {
		return Vector3(fmod(x, f), fmod(y, f), fmod(z, f));
	}

	Vector3<T> &operator -=(const Vector3<T> &v) {
		x -= v.x;
		y -= v.y;
		z -= v.z;
		return *this;
	}

	Vector3<T> &operator -=(const T &f) {
		x -= f;
		y -= f;
		z -= f;
		return *this;
	}

	Vector3<T> &operator +=(const Vector3<T> &v) {
		x += v.x;
		y += v.y;
		z += v.z;
		return *this;
	}

	Vector3<T> &operator +=(const T &f) {
		x += f;
		y += f;
		z += f;
		return *this;
	}

	Vector3<T> &operator *=(const Vector3<T> &v) {
		x *= v.x;
		y *= v.y;
		z *= v.z;
		return *this;
	}

	Vector3<T> &operator *=(const T &f) {
		x *= f;
		y *= f;
		z *= f;
		return *this;
	}

	Vector3<T> &operator /=(const Vector3<T> &v) {
		x /= v.x;
		y /= v.y;
		z /= v.z;
		return *this;
	}

	Vector3<T> &operator /=(const T &f) {
		x /= f;
		y /= f;
		z /= f;
		return *this;
	}

	Vector3<T> &operator %=(const Vector3<T> &v) {
		x = fmod(x, v.x);
		y = fmod(y, v.y);
		z = fmod(z, v.z);
		return *this;
	}

	Vector3<T> &operator %=(const T &f) {
		x = fmod(x, f);
		y = fmod(y, f);
		z = fmod(z, f);
		return *this;
	}

	Vector3<T> &operator =(const Vector3<T> &v) {
		x = v.x;
		y = v.y;
		z = v.z;
		return *this;
	}
};

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
inline Vector3<T> operator *(T f, const Vector3<T> &v) {
	return Vector3<T>(v.x * f, v.y * f, v.z * f);
}

typedef Vector3<double> Vector3d;
typedef Vector3<float> Vector3f;

} // namespace mpc

#endif /* _MPC_VECTOR3_H_ */
