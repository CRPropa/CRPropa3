//-------------------------------------------------------------
// Based on Variant.cc in the Physics eXtension Library (PXL) -
// http://vispa.physik.rwth-aachen.de/                        -
// Licensed under a LGPL-2 or later license                   -
//-------------------------------------------------------------

#include "crpropa/Variant.h"


namespace crpropa {


Variant::Variant() : type(TYPE_NONE) {
}

Variant::Variant(Type t) : type(TYPE_NONE) {
	check(t);
}

Variant::Variant(const Variant& v) : type(TYPE_NONE) {
	copy(v);
}

Variant::Variant(const char* s) {
	data._t_string = new std::string(s);
	type = TYPE_STRING;
}

Variant::~Variant() {
	clear(type);
	delete data._t_string;
}

const char* Variant::getTypeName() const {
	return getTypeName(type);
}

const char* Variant::getTypeName(Type t) {
	if (t == TYPE_NONE)
		return "none";
	else if (t == TYPE_BOOL)
		return "bool";
	else if (t == TYPE_CHAR)
		return "char";
	else if (t == TYPE_UCHAR)
		return "uchar";
	else if (t == TYPE_INT16)
		return "int16";
	else if (t == TYPE_UINT16)
		return "uint16";
	else if (t == TYPE_INT32)
		return "int32";
	else if (t == TYPE_UINT32)
		return "uint32";
	else if (t == TYPE_INT64)
		return "int64";
	else if (t == TYPE_UINT64)
		return "uint64";
	else if (t == TYPE_FLOAT)
		return "float";
	else if (t == TYPE_DOUBLE)
		return "double";
	else if (t == TYPE_LONGDOUBLE)
		return "ldouble";
	else if (t == TYPE_COMPLEXF)
		return "complex_f";
	else if (t == TYPE_COMPLEXD)
		return "complex_d";
	else if (t == TYPE_STRING)
		return "string";
	else if (t == TYPE_VECTOR3F)
		return "Vector3f";
	else if (t == TYPE_VECTOR3D)
		return "Vector3d";
	else if (t == TYPE_VECTOR3C)
		return "Vector3c";
	else if (t == TYPE_VECTOR)
		return "vector";
	else
		return "unknown";
}

const std::type_info& Variant::getTypeInfo() const {
	if (type == TYPE_BOOL) {
		const std::type_info& ti = typeid(data._t_bool);
		return ti;
	} else if (type == TYPE_CHAR) {
		const std::type_info& ti = typeid(data._t_char);
		return ti;
	} else if (type == TYPE_UCHAR) {
		const std::type_info& ti = typeid(data._t_uchar);
		return ti;
	} else if (type == TYPE_INT16) {
		const std::type_info& ti = typeid(data._t_int16);
		return ti;
	} else if (type == TYPE_UINT16) {
		const std::type_info& ti = typeid(data._t_uint16);
		return ti;
	} else if (type == TYPE_INT32) {
		const std::type_info& ti = typeid(data._t_int32);
		return ti;
	} else if (type == TYPE_UINT32) {
		const std::type_info& ti = typeid(data._t_uint32);
		return ti;
	} else if (type == TYPE_INT64) {
		const std::type_info& ti = typeid(data._t_int64);
		return ti;
	} else if (type == TYPE_UINT64) {
		const std::type_info& ti = typeid(data._t_uint64);
		return ti;
	} else if (type == TYPE_FLOAT) {
		const std::type_info& ti = typeid(data._t_float);
		return ti;
	} else if (type == TYPE_DOUBLE)	{
		const std::type_info& ti = typeid(data._t_double);
		return ti;
	} else if (type == TYPE_LONGDOUBLE)	{
		const std::type_info& ti = typeid(data._t_ldouble);
		return ti;
	} else if (type == TYPE_STRING) {
		const std::type_info& ti = typeid(*data._t_string);
		return ti;
	} else if (type == TYPE_COMPLEXF) { // pointer needed?
		const std::type_info& ti = typeid(*data._t_complex_f);
		return ti;
	} else if (type == TYPE_COMPLEXD) {
		const std::type_info& ti = typeid(*data._t_complex_d);
		return ti;
	} else if (type == TYPE_VECTOR3D) {
		const std::type_info& ti = typeid(data._t_vector3d);
		return ti;
	} else if (type == TYPE_VECTOR3F) {
		const std::type_info& ti = typeid(data._t_vector3f);
		return ti;
	} else if (type == TYPE_VECTOR) {
		const std::type_info& ti = typeid(*data._t_vector);
		return ti;
	} else {
		const std::type_info& ti = typeid(0);
		return ti;
	}
}

Variant::Type Variant::toType(const std::string& name) {
	if (name == "none")
		return TYPE_NONE;
	else if (name == "bool")
		return TYPE_BOOL;
	else if (name == "char") 
		return TYPE_CHAR;
	else if (name == "uchar")
		return TYPE_UCHAR;
	else if (name == "int16")
		return TYPE_INT16;
	else if (name == "uint16")
		return TYPE_UINT16;
	else if (name == "int32")
		return TYPE_INT32;
	else if (name == "uint32")
		return TYPE_UINT32;
	else if (name == "int64")
		return TYPE_INT64;
	else if (name == "uint64")
		return TYPE_UINT64;
	else if (name == "float")
		return TYPE_FLOAT;
	else if (name == "double")
		return TYPE_DOUBLE;
	else if (name == "long double")
		return TYPE_LONGDOUBLE;
	else if (name == "complex_f")
		return TYPE_COMPLEXF;
	else if (name == "complex_d")
		return TYPE_COMPLEXD;
	 else if (name == "string")
		return TYPE_STRING;
	else if (name == "Vector3f")
		return TYPE_VECTOR3F;
	else if (name == "Vector3d")
		return TYPE_VECTOR3D;
	else if (name == "Vector3c")
		return TYPE_VECTOR3C;
	else if (name == "vector")
		return TYPE_VECTOR;
	else
		return TYPE_NONE;
}

bool Variant::toBool() const {
	switch (type) {
		case TYPE_BOOL:
			return data._t_bool;
			break;
		case TYPE_CHAR:
			return data._t_char != 0;
			break;
		case TYPE_UCHAR:
			return data._t_uchar != 0;
			break;
		case TYPE_INT16:
			return data._t_int16 != 0;
			break;
		case TYPE_UINT16:
			return data._t_uint16 != 0;
			break;
		case TYPE_INT32:
			return data._t_int32 != 0;
			break;
		case TYPE_UINT32:
			return data._t_uint32 != 0;
			break;
		case TYPE_INT64:
			return data._t_int64 != 0;
			break;
		case TYPE_UINT64:
			return data._t_uint64 != 0;
			break;
		case TYPE_STRING: {
				std::string upperstr(*data._t_string);
				std::transform(upperstr.begin(), upperstr.end(), upperstr.begin(), (int(*) (int)) toupper);
				if (upperstr == "YES" || upperstr == "TRUE" || upperstr == "1")
					return true;
				else if (upperstr == "NO" || upperstr == "FALSE" || upperstr == "0")
					return false;
				else
					throw bad_conversion(type, TYPE_BOOL);
			}
			break;
		case TYPE_COMPLEXF: {
				if (data._t_complex_f->real() == 0 && data._t_complex_f->imag() == 0)
					return true;
				else
					return false;
			}
			break;
		case TYPE_COMPLEXD: {
				if (data._t_complex_d->real() == 0 && data._t_complex_d->imag() == 0)
					return true;
				else
					return false;
			}
			break;
		case TYPE_VECTOR:
			return data._t_vector->size() != 0;
			break;
		case TYPE_FLOAT:
		case TYPE_DOUBLE:
		case TYPE_LONGDOUBLE:
		default:
			throw bad_conversion(type, TYPE_BOOL);
			break;
	}

	return false;
}

float Variant::toFloat() const {
	switch (type) {
		case TYPE_BOOL:
			return data._t_bool ? (float) 1.0 : (float) 0.0; 
			break;
		case TYPE_CHAR:
			return static_cast<float>(data._t_char);
			break;
		case TYPE_UCHAR:
			return static_cast<float>(data._t_uchar);
			break;
		case TYPE_INT16:
			return static_cast<float>(data._t_int16);
			break;
		case TYPE_INT32:
			return static_cast<float>(data._t_int32);
			break;
		case TYPE_INT64:
			return static_cast<float>(data._t_int64);
			break;
		case TYPE_UINT16:
			return static_cast<float>(data._t_uint16);
			break;
		case TYPE_UINT32:
			return static_cast<float>(data._t_uint32);
			break;
		case TYPE_UINT64:
			return static_cast<float>(data._t_uint64);
			break;
		case TYPE_FLOAT:
			return static_cast<float>(data._t_float);
			break;
		case TYPE_DOUBLE:
			return static_cast<float>(data._t_double);
			break;
		case TYPE_LONGDOUBLE:
			return static_cast<float>(data._t_ldouble);
			break;
		case TYPE_STRING:
			return static_cast<float>(std::atof(data._t_string->c_str()));
			break;
		case TYPE_COMPLEXF: {
				if (data._t_complex_f->imag() == 0)
					return static_cast<float>(data._t_complex_f->real());
				else
					throw bad_conversion(type, TYPE_COMPLEXF);
			}
			break;
		case TYPE_COMPLEXD: {
				if (data._t_complex_d->imag() == 0)
					return static_cast<float>(data._t_complex_d->real());
				else
					throw bad_conversion(type, TYPE_COMPLEXD);
			}
			break;
		default:
			throw bad_conversion(type, TYPE_FLOAT);
			break;
	}

	return 0.;
}

double Variant::toDouble() const {
	switch (type) {
		case TYPE_BOOL:
			return data._t_bool ? (double) 1.0 : (double) 0.0; 
			break;
		case TYPE_CHAR:
			return static_cast<double>(data._t_char);
			break;
		case TYPE_UCHAR:
			return static_cast<double>(data._t_uchar);
			break;
		case TYPE_INT16:
			return static_cast<double>(data._t_int16);
			break;
		case TYPE_INT32:
			return static_cast<double>(data._t_int32);
			break;
		case TYPE_INT64:
			return static_cast<double>(data._t_int64);
			break;
		case TYPE_UINT16:
			return static_cast<double>(data._t_uint16);
			break;
		case TYPE_UINT32:
			return static_cast<double>(data._t_uint32);
			break;
		case TYPE_UINT64:
			return static_cast<double>(data._t_uint64);
			break;
		case TYPE_FLOAT:
			return static_cast<double>(data._t_float);
			break;
		case TYPE_DOUBLE:
			return static_cast<double>(data._t_double);
			break;
		case TYPE_LONGDOUBLE:
			return static_cast<double>(data._t_ldouble);
			break;
		case TYPE_COMPLEXF: {
				if (data._t_complex_f->imag() == 0)
					return static_cast<double>(data._t_complex_f->real());
				else
					throw bad_conversion(type, TYPE_COMPLEXF);
			}
			break;
		case TYPE_COMPLEXD: {
				if (data._t_complex_d->imag() == 0)
					return static_cast<double>(data._t_complex_d->real());
				else
					throw bad_conversion(type, TYPE_COMPLEXD);
			}
			break;
		case TYPE_STRING:
			return static_cast<double>(std::atof(data._t_string->c_str()));
			break;
		default:
			throw bad_conversion(type, TYPE_DOUBLE);
			break;
	}

	return 0.;
}

long double Variant::toLongDouble() const {
	switch (type) {
		case TYPE_BOOL:
			return data._t_bool ? (long double) 1.0 : (long double) 0.0; 
			break;
		case TYPE_CHAR:
			return static_cast<long double>(data._t_char);
			break;
		case TYPE_UCHAR:
			return static_cast<long double>(data._t_uchar);
			break;
		case TYPE_INT16:
			return static_cast<long double>(data._t_int16);
			break;
		case TYPE_INT32:
			return static_cast<long double>(data._t_int32);
			break;
		case TYPE_INT64:
			return static_cast<long double>(data._t_int64);
			break;
		case TYPE_UINT16:
			return static_cast<long double>(data._t_uint16);
			break;
		case TYPE_UINT32:
			return static_cast<long double>(data._t_uint32);
			break;
		case TYPE_UINT64:
			return static_cast<long double>(data._t_uint64);
			break;
		case TYPE_FLOAT:
			return static_cast<long double>(data._t_float);
			break;
		case TYPE_DOUBLE:
			return static_cast<long double>(data._t_double);
			break;
		case TYPE_LONGDOUBLE:
			return static_cast<long double>(data._t_ldouble);
			break;
		case TYPE_COMPLEXF: {
				if (data._t_complex_f->imag() == 0)
					return static_cast<long double>(data._t_complex_f->real());
				else
					throw bad_conversion(type, TYPE_COMPLEXF);
			}
			break;
		case TYPE_COMPLEXD: {
				if (data._t_complex_d->imag() == 0)
					return static_cast<long double>(data._t_complex_d->real());
				else
					throw bad_conversion(type, TYPE_COMPLEXD);
			}
			break;
		case TYPE_STRING:
			return static_cast<double>(std::atof(data._t_string->c_str()));
			break;
		default:
			throw bad_conversion(type, TYPE_LONGDOUBLE);
			break;
	}

	return 0.;
}

std::complex<float> Variant::toComplexFloat() const {
	switch (type) {
		case TYPE_COMPLEXF:
			return static_cast<std::complex<float>>(*data._t_complex_f);
			break;
		case TYPE_COMPLEXD:
			return static_cast<std::complex<float>>(*data._t_complex_d);
			break;
		default:
			throw bad_conversion(type, TYPE_COMPLEXF);
			break;
	}
}  

std::complex<double> Variant::toComplexDouble() const {
	switch (type) {
		case TYPE_COMPLEXF:
			return static_cast<std::complex<double>>(*data._t_complex_f);
			break;
		case TYPE_COMPLEXD:
			return static_cast<std::complex<double>>(*data._t_complex_d);
			break;
		default:
			throw bad_conversion(type, TYPE_COMPLEXD);
			break;
	}
}  

Vector3f Variant::toVector3f() const {
	switch (type) {
		case TYPE_VECTOR3F:
			return static_cast<Vector3f>(*data._t_vector3f);
			break;
		case TYPE_VECTOR3D:
			return static_cast<Vector3f>(*data._t_vector3d);
			break;
		default:
			throw bad_conversion(type, TYPE_VECTOR3F);
			break;
	}
}

Vector3d Variant::toVector3d() const {
	switch (type) {
		case TYPE_VECTOR3F:
			return static_cast<Vector3d>(*data._t_vector3f);
			break;
		case TYPE_VECTOR3D:
			return static_cast<Vector3d>(*data._t_vector3d);
			break;
		default:
			throw bad_conversion(type, TYPE_VECTOR3D);
			break;
	}
}

Vector3<std::complex<double>> Variant::toVector3c() const {
	switch (type) {
		case TYPE_VECTOR3C:
			return static_cast<Vector3c>(*data._t_vector3c);
			break;
		default:
			throw bad_conversion(type, TYPE_VECTOR3C);
			break;
	}
}

std::string Variant::toString(const std::string& delimiter) const {
	if (type == TYPE_STRING)
		return *data._t_string;

	std::stringstream ss;

	if (type == TYPE_BOOL) {
		ss << data._t_bool;
	} else if (type == TYPE_CHAR) {
		ss << data._t_char;
	} else if (type == TYPE_UCHAR) {
		ss << data._t_uchar;
	} else if (type == TYPE_INT16) {
		ss << data._t_int16;
	} else if (type == TYPE_UINT16) {
		ss << data._t_uint16;
	} else if (type == TYPE_INT32) {
		ss << data._t_int32;
	} else if (type == TYPE_UINT32) {
		ss << data._t_uint32;
	} else if (type == TYPE_INT64) {
		ss << data._t_int64;
	} else if (type == TYPE_UINT64) {
		ss << data._t_uint64;
	} else if (type == TYPE_FLOAT) {
		ss << std::scientific << data._t_float;
	} else if (type == TYPE_DOUBLE) {
		ss << std::scientific << data._t_double;
	} else if (type == TYPE_LONGDOUBLE) {
		ss << std::scientific << data._t_ldouble;
	} else if (type == TYPE_COMPLEXF) {
		ss << std::scientific << data._t_complex_f->real() << delimiter; 
		ss << std::scientific << data._t_complex_f->imag();
	} else if (type == TYPE_COMPLEXD) {
		ss << std::scientific << data._t_complex_d->real() << delimiter; 
		ss << std::scientific << data._t_complex_d->imag();
	} else if (type == TYPE_VECTOR3F) {
		ss << *data._t_vector3f;
	} else if (type == TYPE_VECTOR3D) {
		ss << *data._t_vector3d;
	} else if (type == TYPE_VECTOR) {
		ss << *data._t_vector;
	}

	return ss.str();
}

Variant::vector_t Variant::toVector() const {
	if (type == TYPE_VECTOR)
		return *data._t_vector;
	else
		throw bad_conversion(type, TYPE_VECTOR);
}

Variant Variant::fromString(const std::string& s, Type t) {
	std::stringstream ss(s);

	if (t == TYPE_BOOL) {
		std::string upperstr(s);
		std::transform(upperstr.begin(), upperstr.end(), upperstr.begin(), (int (*)(int)) toupper);
		if (upperstr == "YES" || upperstr == "TRUE" || upperstr == "1") 
			return Variant(true);
		else if (upperstr == "NO" || upperstr == "FALSE" || upperstr == "0")
			return Variant(false);
		throw bad_conversion(t, TYPE_BOOL);
	} else if (t == TYPE_CHAR) {
		char c;
		ss >> c;
		return Variant(c);
	} else if (t == TYPE_UCHAR) {
		unsigned char c;
		ss >> c;
		return Variant(c);
	} else if (t == TYPE_INT16) {
		int16_t c;
		ss >> c;
		return Variant(c);
	} else if (t == TYPE_INT32) {
		int32_t c;
		ss >> c;
		return Variant(c);
	} else if (t == TYPE_INT64) {
		int64_t c;
		ss >> c;
		return Variant(c);
	} else if (t == TYPE_UINT16) {
		uint16_t c;
		ss >> c;
		return Variant(c);
	} else if (t == TYPE_UINT32) {
		uint32_t c;
		ss >> c;
		return Variant(c);
	} else if (t == TYPE_UINT64) {
		uint64_t c;
		ss >> c;
		return Variant(c);
	} else if (t == TYPE_FLOAT) {
		float c;
		ss >> c;
		return Variant(c);
	} else if (t == TYPE_DOUBLE) {
		double c;
		ss >> c;
		return Variant(c);
	} else if (t == TYPE_LONGDOUBLE) {
		long double c;
		ss >> c;
		return Variant(c);
	} else if (t == TYPE_STRING) {
		return Variant(s);
	} else if (t == TYPE_COMPLEXF) {
		float _vr, _vi;
		ss >> _vr >> _vi;
		complex_f v(_vr, _vi);
		return Variant(v);
	} else if (t == TYPE_COMPLEXD) {
		double _vr, _vi;
		ss >> _vr >> _vi;
		complex_d v(_vr, _vi);
		return Variant(v);
	} else if (t == TYPE_VECTOR3F) {
		Vector3f v;
		float _val;
		ss >> _val;
		v.setX(_val);
		ss >> _val;
		v.setY(_val);
		ss >> _val;
		v.setZ(_val);
		return Variant(v);
	} else if (t == TYPE_VECTOR3D) {
		Vector3d v;
		double _val;
		ss >> _val;
		v.setX(_val);
		ss >> _val;
		v.setY(_val);
		ss >> _val;
		v.setZ(_val);
		return Variant(v);
	} else if (t == TYPE_VECTOR3C) {
		Vector3c v;
		std::complex<double> _val;
		ss >> _val;
		v.setX(_val);
		ss >> _val;
		v.setY(_val);
		ss >> _val;
		v.setZ(_val);
		return Variant(v);
	} else if (t == TYPE_VECTOR) {
		// std::regex useless("(|)|[|]| ");
		vector_t v;
		while (ss.good()) {
			std::string s;
			v.push_back(s);
		}
	} else {
		std::string msg; 
		msg += "fromString not implemented for type ";
        msg += getTypeName(t);
        throw std::runtime_error("Variant: " + msg);
	}
}

void Variant::clear(Type t) {
	if (t == TYPE_STRING) 
		safeDelete(data._t_string);
	else if (t == TYPE_VECTOR3F)
		safeDelete(data._t_vector3f);
	else if (t == TYPE_VECTOR3D)
		safeDelete(data._t_vector3d);
	else if (t == TYPE_VECTOR3C)
		safeDelete(data._t_vector3c);
	else if (t == TYPE_COMPLEXF)
		safeDelete(data._t_complex_f);
	else if (t == TYPE_COMPLEXD)
		safeDelete(data._t_complex_d);
	else if (t == TYPE_VECTOR)
		safeDelete(data._t_vector);

	// set the type to TYPE_NONE which will be used for checks
    type = TYPE_NONE;
    check(t);
}

void Variant::copy(const Variant& v) {
	Type t = v.type;

	if (t == TYPE_BOOL) 
		operator = (v.data._t_bool);
	else if (t == TYPE_CHAR)
		operator = (v.data._t_char);
	else if (t == TYPE_UCHAR)
		operator = (v.data._t_uchar);
	else if (t == TYPE_INT16)
		operator = (v.data._t_int16);
	else if (t == TYPE_UINT16)
		operator = (v.data._t_uint16);
	else if (t == TYPE_INT32)
		operator = (v.data._t_int32);
	else if (t == TYPE_UINT32)
		operator = (v.data._t_uint32);
	else if (t == TYPE_INT64)
		operator = (v.data._t_int64);
	else if (t == TYPE_UINT64)
		operator = (v.data._t_uint64);
	else if (t == TYPE_FLOAT)
		operator = (v.data._t_float);
	else if (t == TYPE_DOUBLE)
		operator = (v.data._t_double);
	else if (t == TYPE_LONGDOUBLE)
		operator = (v.data._t_ldouble);
	else if (t == TYPE_STRING)
		operator = (*v.data._t_string);
	else if (t == TYPE_COMPLEXF)
		operator = (*v.data._t_complex_f);
	else if (t == TYPE_COMPLEXD)
		operator = (*v.data._t_complex_d);
	else if (t == TYPE_VECTOR3F)
		operator = (*v.data._t_vector3f);
	else if (t == TYPE_VECTOR3D)
		operator = (*v.data._t_vector3d);
	else if (t == TYPE_VECTOR3C)
		operator = (*v.data._t_vector3c);
	else if (t == TYPE_VECTOR)
		operator = (*v.data._t_vector);
	else
		type = TYPE_NONE;
}

void Variant::check(const Type t) const {
	if (type != t) {
		throw bad_conversion(type, t);
	}
}

void Variant::check(const Type t) {
	if (type == TYPE_NONE) {
		memset(&data, 0, sizeof(data));
		if (t == TYPE_VECTOR3F)
			data._t_vector3f = new Vector3f(); 
		else if (t == TYPE_VECTOR3D)
			data._t_vector3d = new Vector3d(); 
		else if (t == TYPE_VECTOR3C)
			data._t_vector3c = new Vector3c(); 
		else if (t == TYPE_COMPLEXF)
			data._t_complex_f = new complex_f(); 
		else if (t == TYPE_COMPLEXD)
			data._t_complex_d = new complex_d(); 
		else if (t == TYPE_STRING)
			data._t_string = new std::string();
		else if (t == TYPE_VECTOR)
			data._t_vector = new vector_t();
		else
			type = t;
	} else if (type != t)  {
	throw bad_conversion(type, t);
	}
}

bool Variant::isValid() const {
	return (type != TYPE_NONE);
}

size_t Variant::size() const {
	if (type == TYPE_VECTOR) {
		return data._t_vector->size();
	} else {
		std::string msg; 
		msg += "size() not implemented for type ";
		msg += getTypeName(type);
		throw std::runtime_error("Variant: " + msg);
	}
}

size_t Variant::getSizeOf() const {
	switch (type) {
		case TYPE_BOOL:
			return sizeof(data._t_bool);
			break;
		case TYPE_CHAR:
			return sizeof(data._t_char);
			break;
		case TYPE_UCHAR:
			return sizeof(data._t_uchar);
			break;
		case TYPE_INT16:
			return sizeof(data._t_int16);
			break;
		case TYPE_UINT16:
			return sizeof(data._t_uint16);
			break;
		case TYPE_INT32:
			return sizeof(data._t_int32);
			break;
		case TYPE_UINT32:
			return sizeof(data._t_uint32);
			break;
		case TYPE_INT64:
			return sizeof(data._t_int64);
			break;
		case TYPE_UINT64:
			return sizeof(data._t_uint64);
			break;
		case TYPE_FLOAT:
			return sizeof(data._t_float);
			break;
		case TYPE_DOUBLE:
			return sizeof(data._t_double);
			break;
		case TYPE_LONGDOUBLE:
			return sizeof(data._t_ldouble);
			break;
		case TYPE_COMPLEXF:
			return sizeof(data._t_complex_f);
			break;
		case TYPE_COMPLEXD:
			return sizeof(data._t_complex_d);
			break;
		case TYPE_VECTOR3F:
			return sizeof(data._t_vector3f);
			break;
		case TYPE_VECTOR3D:
			return sizeof(data._t_vector3d);
			break;
		case TYPE_VECTOR3C:
			return sizeof(data._t_vector3c);
			break;
		case TYPE_STRING: {
				size_t len = strlen(data._t_string->c_str() + 1);
				return len;
			}
			break;
		case TYPE_NONE:
			return 0;
			break;
		default:
			throw std::runtime_error("Function getSize() cannot handle this type.");
	}
}

size_t Variant::getSize() const {
	return getSizeOf();
}

void Variant::resize(size_t i) {
	check(TYPE_VECTOR);
	return data._t_vector->resize(i);
}

size_t Variant::copyToBuffer(void* buffer) {
	if (type == TYPE_BOOL) {
		memcpy(buffer, &data._t_bool, sizeof(bool));
		return sizeof(data._t_bool);
	} else if (type == TYPE_CHAR) {
		memcpy(buffer, &data._t_char, sizeof(char));
		return sizeof(data._t_char);
	}  else if (type == TYPE_UCHAR) {
		memcpy(buffer, &data._t_uchar, sizeof(unsigned char));
		return sizeof(data._t_uchar);
	} else if (type == TYPE_INT16) {
		memcpy(buffer, &data._t_int16, sizeof(int16_t));
		return sizeof(data._t_int16);
	} else if (type == TYPE_UINT16) {
		memcpy(buffer, &data._t_uint16, sizeof(uint16_t));
		return sizeof(data._t_uint16);
	} else if (type == TYPE_INT32) {
		memcpy(buffer, &data._t_int32, sizeof(int32_t));
		return sizeof(data._t_int32);
	} else if (type == TYPE_UINT32) {
		memcpy(buffer, &data._t_uint32, sizeof(uint32_t));
		return sizeof(data._t_uint32);
	} else if (type == TYPE_INT64) {
		memcpy(buffer, &data._t_int64, sizeof(int64_t));
		return sizeof(data._t_int64);
	} else if (type == TYPE_UINT64) {
		memcpy(buffer, &data._t_uint64, sizeof(uint64_t));
		return sizeof(data._t_uint64);
	} else if (type == TYPE_FLOAT) {
		memcpy(buffer, &data._t_float, sizeof(float));
		return sizeof(data._t_float);
	} else if (type == TYPE_DOUBLE) {
		memcpy(buffer, &data._t_double, sizeof(double));
		return sizeof(data._t_double);
	} else if (type == TYPE_LONGDOUBLE) {
		memcpy(buffer, &data._t_ldouble, sizeof(long double));
		return sizeof(data._t_ldouble);
	}  else if (type == TYPE_STRING) {
		size_t len = data._t_string->size();
		memcpy(buffer, data._t_string->c_str(), len);
		return len;
	} else if (type == TYPE_NONE) {
		return 0;
	} else {
		throw std::runtime_error("This type cannot be handled by copyToBuffer().");
	}
}

Variant& Variant::operator=(const Variant &v) {
	copy(v);
	return *this;
}

bool Variant::operator==(const Variant& v) const {
	if (type != v.type)
		return false;

	if (type == TYPE_BOOL) {
		return (data._t_bool == v.data._t_bool);
	} else if (type == TYPE_CHAR) {
		return (data._t_char == v.data._t_char);
	} else if (type == TYPE_UCHAR) {
		return (data._t_uchar == v.data._t_uchar);
	} else if (type == TYPE_INT16) {
		return (data._t_int16 == v.data._t_int16);
	} else if (type == TYPE_UINT16) {
		return (data._t_uint16 == v.data._t_uint16);
	} else if (type == TYPE_INT32) {
		return (data._t_int32 == v.data._t_int32);
	} else if (type == TYPE_UINT32) {
		return (data._t_uint32 == v.data._t_uint32);
	} else if (type == TYPE_INT64) {
		return (data._t_int64 == v.data._t_int64);
	} else if (type == TYPE_UINT64) {
		return (data._t_uint64 == v.data._t_uint64);
	} else if (type == TYPE_FLOAT) {
		return (data._t_float == v.data._t_float);
	} else if (type == TYPE_DOUBLE) {
		return (data._t_double == v.data._t_double);
	} else if (type == TYPE_LONGDOUBLE) {
		return (data._t_ldouble == v.data._t_ldouble);
	} else if (type == TYPE_STRING) {
		return (*data._t_string == *v.data._t_string);
	} else if (type == TYPE_COMPLEXF) {
		return (*data._t_complex_f == *v.data._t_complex_f);
	} else if (type == TYPE_COMPLEXD) {
		return (*data._t_complex_d == *v.data._t_complex_d);
	} else if (type == TYPE_VECTOR3F) {
		return (*data._t_vector3f == *v.data._t_vector3f);
	} else if (type == TYPE_VECTOR3D) {
		return (*data._t_vector3d == *v.data._t_vector3d);
	} else if (type == TYPE_VECTOR3C) {
		return (*data._t_vector3c == *v.data._t_vector3c);
	} else if (type == TYPE_VECTOR) {
		return (*data._t_vector == *v.data._t_vector);
	} else {
		throw std::runtime_error("compare operator not implemented");
	}
}

bool Variant::operator!=(const Variant& v) const {
	if (type != v.type)
		return true;
	
	if (*this == v)
		return false;
	else
		return true;
}

bool Variant::operator!=(const char* v) const {
	check(TYPE_STRING);
	return data._t_string->compare(v) != 0;
}

Variant& Variant::operator[](size_t i) {
	check(TYPE_VECTOR);
	return (*data._t_vector)[i];
}

const Variant& Variant::operator[](size_t i) const {
	check(TYPE_VECTOR);
	return (*data._t_vector)[i];
}

Variant::operator std::vector<Variant>&() {
	check(TYPE_VECTOR);
	return *data._t_vector;
}

Variant::operator const std::vector<Variant>&() const {
	check(TYPE_VECTOR);
	return *data._t_vector;
}



#define INT_CASE(from_var, from_type, to_type, to) 													                   \
    case Variant::from_type: {                                                                                         \
		if (data._t_##from_var <std::numeric_limits<to>::min() || data._t_##from_var> std::numeric_limits<to>::max())  \
			throw bad_conversion(type, to_type);                                                                       \
		else                                                                                                           \
			return static_cast<to>(data._t_##from_var);                                                                \
		}                                                                                                              \
    	break;                                                                                                         \

#define INT_FUNCTION(to_type, fun, to)                                                                                 \
	to Variant::fun() const {                                                                                          \
		switch (type) {                                                                                                \
			case Variant::TYPE_BOOL:                                                                                   \
				return data._t_bool ? 1 : 0;                                                                           \
				break;                                                                                                 \
				INT_CASE(char, TYPE_CHAR, to_type, to)                                                                 \
				INT_CASE(uchar, TYPE_UCHAR, to_type, to)                                                               \
				INT_CASE(int16, TYPE_INT16, to_type, to)                                                               \
				INT_CASE(uint16, TYPE_UINT16, to_type, to)                                                             \
				INT_CASE(int32, TYPE_INT32, to_type, to)                                                               \
				INT_CASE(uint32, TYPE_UINT32, to_type, to)                                                             \
				INT_CASE(int64, TYPE_INT64, to_type, to)                                                               \
				INT_CASE(uint64, TYPE_UINT64, to_type, to)                                                             \
				INT_CASE(float, TYPE_FLOAT, to_type, to)                                                               \
				INT_CASE(double, TYPE_DOUBLE, to_type, to)                                                             \
				INT_CASE(ldouble, TYPE_LONGDOUBLE, to_type, to)                                                        \
			case Variant::TYPE_STRING: {                                                                               \
				long l = atol(data._t_string->c_str());                                                                \
				if (l <std::numeric_limits<to>::min() || l > std::numeric_limits<to>::max())                           \
					throw bad_conversion(type, to_type);                                                               \
				else                                                                                                   \
					return l;                                                                                          \
				}                                                                                                      \
				break;                                                                                                 \
			case Variant::TYPE_COMPLEXF:                                                                               \
			case Variant::TYPE_COMPLEXD:                                                                               \
			case Variant::TYPE_VECTOR3F:                                                                               \
			case Variant::TYPE_VECTOR3D:                                                                               \
			case Variant::TYPE_VECTOR3C:                                                                               \
			case Variant::TYPE_VECTOR:                                                                                 \
			case Variant::TYPE_NONE:                                                                                   \
				throw bad_conversion(type, to_type);                                                                   \
				break;                                                                                                 \
		}                                                                                                              \
		return 0;                                                                                                      \
	} 

INT_FUNCTION(TYPE_CHAR, toChar, char)
INT_FUNCTION(TYPE_UCHAR, toUChar, unsigned char)
INT_FUNCTION(TYPE_INT16, toInt16, int16_t)
INT_FUNCTION(TYPE_UINT16, toUInt16, uint16_t)
INT_FUNCTION(TYPE_INT32, toInt32, int32_t)
INT_FUNCTION(TYPE_UINT32, toUInt32, uint32_t)
INT_FUNCTION(TYPE_INT64, toInt64, int64_t)
INT_FUNCTION(TYPE_UINT64, toUInt64, uint64_t)



std::ostream& operator <<(std::ostream& os, const Variant& v) {
	switch (v.getType()) {
		case Variant::TYPE_BOOL:
			os << v.asBool();
			break;
		case Variant::TYPE_CHAR:
			os << v.asChar();
			break;
		case Variant::TYPE_UCHAR:
			os << v.asUChar();
			break;
		case Variant::TYPE_INT16:
			os << v.asInt16();
			break;
		case Variant::TYPE_UINT16:
			os << v.asUInt16();
			break;
		case Variant::TYPE_INT32:
			os << v.asInt32();
			break;
		case Variant::TYPE_UINT32:
			os << v.asUInt32();
			break;
		case Variant::TYPE_INT64:
			os << v.asInt64();
			break;
		case Variant::TYPE_UINT64:
			os << v.asUInt64();
			break;
		case Variant::TYPE_FLOAT:
			os << v.asFloat();
			break;
		case Variant::TYPE_DOUBLE:
			os << v.asDouble();
			break;
		case Variant::TYPE_LONGDOUBLE:
			os << v.asLongDouble();
			break;
		case Variant::TYPE_COMPLEXF:
			os << v.asComplexFloat();
			break;
		case Variant::TYPE_COMPLEXD:
			os << v.asComplexDouble();
			break;
		case Variant::TYPE_STRING:
			os << v.asString();
			break;
		case Variant::TYPE_VECTOR3F:
			os << v.asVector3f();
			break;
		case Variant::TYPE_VECTOR3D:
			os << v.asVector3d();
			break;
		case Variant::TYPE_VECTOR3C:
			os << v.asVector3c();
			break;
		default:
			break;
	}

	return os;
}



} // namespace crpropa
