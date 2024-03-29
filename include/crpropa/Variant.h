//-------------------------------------------------------------
// Based on Variant.hh of the Physics eXtension Library (PXL) -
// http://vispa.physik.rwth-aachen.de/                        -
// Licensed under a LGPL-2 or later license                   -
//-------------------------------------------------------------

#ifndef CRPROPA_VARIANT_H
#define CRPROPA_VARIANT_H

#include <complex>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <limits>
#include <string>
#include <sstream>
#include <stdexcept>
#include <typeinfo>

#include "crpropa/Vector3.h"



// defines copy constructor X(const X&), isX(), asX(), fromX(), toX(), op(), op=, op==, op!=
#define VARIANT_ADD_TYPE_DECL_POD(NAME, TYPE, VALUE, FIELD) \
	Variant(const VALUE& v) { \
		data._t_##FIELD = v; \
		type = TYPE;  \
	} \
	bool is##NAME() const { \
		return type == TYPE; \
	} \
	VALUE& as##NAME() { \
		check(TYPE); \
		return data._t_##FIELD; \
	} \
	const VALUE& as##NAME() const { \
		check(TYPE); \
		return data._t_##FIELD; \
	} \
	static Variant from##NAME(const VALUE& v) { \
		return Variant(v); \
	} \
	VALUE to##NAME() const; \
	operator VALUE() const { \
		return to##NAME(); \
	} \
	Variant& operator=(const VALUE& v) { \
		clear(); \
		type = TYPE; \
		data._t_##FIELD = v; \
		return *this; \
	} \
	bool operator==(const VALUE& v) const { \
		check(TYPE); \
		return data._t_##FIELD == v; \
	} \
	bool operator!=(const VALUE& v) const { \
		check(TYPE); \
		return data._t_##FIELD != v; \
	} \

// defines isX(), asX(), fromX()
#define VARIANT_ADD_TYPE_DECL_PTR_BASE(NAME, TYPE, VALUE, FIELD) \
	bool is##NAME() const { \
		return type == TYPE; \
	} \
	VALUE& as##NAME() { \
		check(TYPE); \
		return *data._t_##FIELD; \
	} \
	const VALUE& as##NAME() const { \
		check(TYPE); \
		return *data._t_##FIELD; \
	} \
	static Variant from##NAME(const VALUE& v) { \
		return Variant(v); \
	} \

// defines isX(), asX(), fromX(), and copy constructor X(const X&), op=, op==, op!=
#define VARIANT_ADD_TYPE_DECL_PTR(NAME, TYPE, VALUE, FIELD) \
	VARIANT_ADD_TYPE_DECL_PTR_BASE(NAME, TYPE, VALUE, FIELD) \
	Variant(const VALUE& v) { \
		data._t_##FIELD = new VALUE(v); \
		type = TYPE; \
	} \
	Variant& operator=(const VALUE& v) { \
		if (type != TYPE) { \
			clear(); \
			data._t_##FIELD = new VALUE();\
		} \
		type = TYPE; \
		(*data._t_##FIELD) = v; \
		return *this; \
	} \
	bool operator==(const VALUE& v) const { \
		check(TYPE); \
		return *data._t_##FIELD == v; \
	} \
	bool operator!=(const VALUE& v) const { \
		return !(*this == v); \
	} \

#define VARIANT_ADD_ITER_DECL_PTR(NAME, TYPE, FIELD) \
	typedef FIELD##_t::iterator FIELD##_iterator; \
	typedef FIELD##_t::const_iterator FIELD##_const_iterator; \
	inline FIELD##_iterator begin##NAME() { \
		check(TYPE); \
		return data._t_##FIELD->begin(); \
	} \
	inline FIELD##_iterator end##NAME() { \
		check(TYPE); \
		return data._t_##FIELD->end(); \
	} \
	inline FIELD##_const_iterator begin##NAME() const { \
		check(TYPE); \
		return data._t_##FIELD->begin(); \
	} \
	inline FIELD##_const_iterator end##NAME() const { \
		check(TYPE); \
		return data._t_##FIELD->end(); \
	}                                                                                



namespace crpropa {

/**
 @class Variant
 @brief storage container for data types as e.g. int, float, string, etc.

 Allows storage of multiple data types in one base class. Used to construct a map of `arbitrary' data types.
 Note that most default C++ types allow default conversions from `Variant` to the corresponding type.
 Types that require an explicit call via `toTargetType()` are: complex (float and double), Vector3, and vector<Variant>.
 */
class Variant {
public:
	enum Type {
		TYPE_NONE = 0,
		TYPE_BOOL,
		TYPE_CHAR,
		TYPE_UCHAR,
		TYPE_INT16,
		TYPE_UINT16,
		TYPE_INT32,
		TYPE_UINT32,
		TYPE_INT64,
		TYPE_UINT64,
		TYPE_FLOAT,
		TYPE_DOUBLE,
		TYPE_LONGDOUBLE,
		TYPE_COMPLEXF,
		TYPE_COMPLEXD,
		TYPE_STRING,
		TYPE_VECTOR3F,
		TYPE_VECTOR3D,
		TYPE_VECTOR3C,
		TYPE_VECTOR
	};

	class bad_conversion: public std::exception {
		protected:
			std::string msg;
		public:
			const char* what() const throw () {
				return msg.c_str();
			}

			bad_conversion(Type f, Type t) {
				msg = "Variant: bad conversion from '";
				msg += Variant::getTypeName(f);
				msg += "' to '";
				msg += Variant::getTypeName(t);
				msg += "'";
			}

			~bad_conversion() throw () {
			}
	};

	typedef std::complex<float> complex_f;
	typedef std::complex<double> complex_d;
	typedef Vector3<std::complex<double>> Vector3c;
	typedef std::vector<Variant> vector_t;

protected:
	Type type;

	union {
		bool _t_bool;
		char _t_char;
		unsigned char _t_uchar;
		int16_t _t_int16;
		uint16_t _t_uint16;
		int32_t _t_int32;
		uint32_t _t_uint32;
		int64_t _t_int64;
		uint64_t _t_uint64;
		float _t_float;
		double _t_double;
		long double _t_ldouble;
		complex_f* _t_complex_f;
		complex_d* _t_complex_d;
		std::string* _t_string;
		Vector3f* _t_vector3f;
		Vector3d* _t_vector3d;
		Vector3c* _t_vector3c;
		vector_t* _t_vector;
	} data;


public:
	Variant();
	Variant(Type t);
	Variant(const Variant& v);
	Variant(const char* s);
	~Variant();
	inline Type getType() const {
		return type;
	}
	const char* getTypeName() const;
	static const char* getTypeName(Type t);
	const std::type_info& getTypeInfo() const;
	static Type toType(const std::string& name);		
	std::string toString(const std::string& delimiter = "\t") const;
	std::complex<float> toComplexFloat() const;
	std::complex<double> toComplexDouble() const;
	Vector3f toVector3f() const;
	Vector3d toVector3d() const;
	Vector3c toVector3c() const;
	vector_t toVector() const;
	static Variant fromString(const std::string& str, Type type);
	void clear(Type t = TYPE_NONE);
	bool isValid() const;
	size_t size() const;
	size_t getSizeOf() const;
	size_t getSize() const;
	void resize(size_t i);
	size_t copyToBuffer(void* buffer);
	operator std::string() const {
		return toString();
	}
	Variant& operator=(const Variant& v);
	bool operator==(const Variant& v) const;
	bool operator!=(const Variant& v) const;
	bool operator!=(const char* a) const;
	Variant& operator[](size_t i);
	inline Variant& operator[](int i) {
		return operator[]((size_t) i);
	}
	const Variant& operator[](size_t i) const;
	const Variant& operator[](int i) const {
		return operator[]((size_t) i);
	}
	operator vector_t&();
	operator const vector_t&() const;


	template<class T>
	T to() const {
		throw bad_conversion(type, TYPE_NONE);
	}

	// automatically-generated functions	
	VARIANT_ADD_TYPE_DECL_POD(Bool, TYPE_BOOL, bool, bool)
	VARIANT_ADD_TYPE_DECL_POD(Char, TYPE_CHAR, char, char)
	VARIANT_ADD_TYPE_DECL_POD(UChar, TYPE_UCHAR, unsigned char, uchar)
	VARIANT_ADD_TYPE_DECL_POD(Int16, TYPE_INT16, int16_t, int16)
	VARIANT_ADD_TYPE_DECL_POD(UInt16, TYPE_UINT16, uint16_t, uint16)
	VARIANT_ADD_TYPE_DECL_POD(Int32, TYPE_INT32, int32_t, int32)
	VARIANT_ADD_TYPE_DECL_POD(UInt32, TYPE_UINT32, uint32_t, uint32)
	VARIANT_ADD_TYPE_DECL_POD(Int64, TYPE_INT64, int64_t, int64)
	VARIANT_ADD_TYPE_DECL_POD(UInt64, TYPE_UINT64, uint64_t, uint64)
	VARIANT_ADD_TYPE_DECL_POD(Float, TYPE_FLOAT, float, float)
	VARIANT_ADD_TYPE_DECL_POD(Double, TYPE_DOUBLE, double, double)
	VARIANT_ADD_TYPE_DECL_POD(LongDouble, TYPE_LONGDOUBLE, long double, ldouble)
	VARIANT_ADD_TYPE_DECL_PTR(ComplexFloat, TYPE_COMPLEXF, std::complex<float>, complex_f)
	VARIANT_ADD_TYPE_DECL_PTR(ComplexDouble, TYPE_COMPLEXD, std::complex<double>, complex_d)
	VARIANT_ADD_TYPE_DECL_PTR(String, TYPE_STRING, std::string, string)
	VARIANT_ADD_TYPE_DECL_PTR(Vector3f, TYPE_VECTOR3F, Vector3f, vector3f)
	VARIANT_ADD_TYPE_DECL_PTR(Vector3d, TYPE_VECTOR3D, Vector3d, vector3d)
	VARIANT_ADD_TYPE_DECL_PTR(Vector3c, TYPE_VECTOR3C, Vector3c, vector3c)
	VARIANT_ADD_TYPE_DECL_PTR(Vector, TYPE_VECTOR, vector_t, vector)
	VARIANT_ADD_ITER_DECL_PTR(Vector, TYPE_VECTOR, vector)
		
private:
	void copy(const Variant& v);
	void check(const Type t) const;
	void check(const Type t);
};

#define VARIANT_TO_DECL(NAME, VALUE) \
	template<> inline VALUE Variant::to<VALUE>() const { \
		return to##NAME(); \
	} \

// declare type conversion functions
// not implemented for Vector3 and complex_*
VARIANT_TO_DECL(Bool, bool)
VARIANT_TO_DECL(Char, char)
VARIANT_TO_DECL(UChar, unsigned char)
VARIANT_TO_DECL(Int16, int16_t)
VARIANT_TO_DECL(UInt16, uint16_t)
VARIANT_TO_DECL(Int32, int32_t)
VARIANT_TO_DECL(UInt32, uint32_t)
VARIANT_TO_DECL(Int64, int64_t)
VARIANT_TO_DECL(UInt64, uint64_t)
VARIANT_TO_DECL(Float, float)
VARIANT_TO_DECL(Double, double)
VARIANT_TO_DECL(LongDouble, long double)
VARIANT_TO_DECL(String, std::string)
VARIANT_TO_DECL(Vector, Variant::vector_t)

std::ostream& operator <<(std::ostream& os, const Variant &v);


/**
 Taken from PXL
 https://git.rwth-aachen.de/3pia/pxl/pxl/-/blob/master/core/include/pxl/core/Functions.hh
 */
template <class T> 
inline void safeDelete(T*& p) {
	if (p) {
		delete p;
		p = 0;
	}
}



} // namespace crpropa 

#endif // CRPROPA_VARIANT
