//-------------------------------------------------------------
// Based on Variant.hh in the Physics eXtension Library (PXL) -
// http://vispa.physik.rwth-aachen.de/                        -
// Licensed under a LGPL-2 or later license                   -
//-------------------------------------------------------------

#ifndef VARIANT_HH
#define VARIANT_HH

#include <iostream>
#include <string>
#include <cstring>
#include <typeinfo>
#include <sstream>
#include <cstdlib>
#include <stdexcept>
#include <limits>
#include <stdint.h>

#define VARIANT_ADD_TYPE_DECL_POD(NAME, TYPE, VALUE) \
	bool is ## NAME() const { return (type == TYPE); } \
	operator VALUE () const { return to ## NAME(); } \
	VALUE &as ## NAME() { check(TYPE); return data._##NAME; } \
	const VALUE &as ## NAME() const	{ check(TYPE); return data._##NAME; } \
	static Variant from ## NAME(const VALUE &a) { return Variant(a); } \
	VALUE to ## NAME() const; \
	Variant &operator = (const VALUE &a) { clear(); type = TYPE; data._##NAME = a; return *this; } \
	bool operator != (const VALUE &a) const { check(TYPE); return data._##NAME != a; } \
	bool operator == (const VALUE &a) const { check(TYPE); return data._##NAME == a; } \
	Variant(const VALUE &a) { data._ ## NAME = a; type = TYPE; }

#define VARIANT_ADD_TYPE_DECL_PTR_BASE(NAME, TYPE, VALUE) \
	bool is ## NAME() const { return (type == TYPE); } \
	VALUE &as ## NAME() { check(TYPE); return *data._##NAME; } \
	const VALUE &as ## NAME() const	{ check(TYPE); return *data._##NAME; } \
	static Variant from ## NAME(const VALUE &a) { return Variant(a); } \

#define VARIANT_ADD_TYPE_DECL_PTR(NAME, TYPE, VALUE) \
	bool operator != (const VALUE &a) const { check(TYPE); return *data._##NAME != a; } \
	bool operator == (const VALUE &a) const { check(TYPE); return *data._##NAME == a; } \
	VARIANT_ADD_TYPE_DECL_PTR_BASE(NAME, TYPE, VALUE) \
	Variant &operator =(const VALUE &a) { if (type != TYPE) { clear(); data._##NAME = new VALUE; } type = TYPE; (*data._##NAME) = a; return *this; } \
	Variant(const VALUE &a) { data._ ## NAME = new VALUE(a); type = TYPE; }

namespace crpropa
{

/**
 @class Variant
 @brief storage container for data types as e.g. int, float, string, etc.

 Allows storage of multiple data types in one base class. used to construct a
 map of `arbitrarry' data types.

 */
class Variant
{
public:
	enum Type
	{
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
		TYPE_STRING
	};

	class bad_conversion: public std::exception
	{
		std::string msg;
	public:
		const char* what() const throw ()
		{
			return msg.c_str();
		}
		bad_conversion(Type f, Type t)
		{
			msg = "Variant: bad conversion from '";
			msg += Variant::getTypeName(f);
			msg += "' to '";
			msg += Variant::getTypeName(t);
			msg += "'";
		}
		~bad_conversion() throw ()
		{
		}
	};

	Variant();
	~Variant();

	Variant(const Variant& a);

	const std::type_info& getTypeInfo() const;

	const char * getTypeName() const
	{
		return getTypeName(type);
	}

	static Type toType(const std::string &name);

	static const char *getTypeName(Type type);

	template<class T>
	T to() const
	{
		throw bad_conversion(type, TYPE_NONE);
	}

	Type getType() const
	{
		return type;
	}

	bool operator ==(const Variant &a) const;

	bool operator !=(const Variant &a) const;

	Variant &operator =(const Variant &a)
	{
		copy(a);
		return *this;
	}

	bool isValid()
	{
		return (type != TYPE_NONE);
	}

	VARIANT_ADD_TYPE_DECL_POD(Bool, TYPE_BOOL, bool)

	VARIANT_ADD_TYPE_DECL_POD(Char, TYPE_CHAR, char)

	VARIANT_ADD_TYPE_DECL_POD(UChar, TYPE_UCHAR, unsigned char)

	VARIANT_ADD_TYPE_DECL_POD(Int16, TYPE_INT16, int16_t)

	VARIANT_ADD_TYPE_DECL_POD(UInt16, TYPE_UINT16, uint16_t)

	VARIANT_ADD_TYPE_DECL_POD(Int32, TYPE_INT32, int32_t)

	VARIANT_ADD_TYPE_DECL_POD(UInt32, TYPE_UINT32, uint32_t)

	VARIANT_ADD_TYPE_DECL_POD(Int64, TYPE_INT64, int64_t)

	VARIANT_ADD_TYPE_DECL_POD(UInt64, TYPE_UINT64, uint64_t)

	VARIANT_ADD_TYPE_DECL_POD(Float, TYPE_FLOAT, float)

	VARIANT_ADD_TYPE_DECL_POD(Double, TYPE_DOUBLE, double)

	VARIANT_ADD_TYPE_DECL_PTR(String, TYPE_STRING, std::string)
	Variant(const char *s);
	std::string toString() const;
	static Variant fromString(const std::string &str, Type type);
	operator std::string() const
	{
		return toString();
	}
	bool operator !=(const char *a) const
	{
		check(TYPE_STRING);
		return data._String->compare(a) != 0;
	}

	// io
	void clear();

protected:
	Type type;

	union
	{
		bool _Bool;
		char _Char;
		unsigned char _UChar;
		int16_t _Int16;
		uint16_t _UInt16;
		int32_t _Int32;
		uint32_t _UInt32;
		int64_t _Int64;
		uint64_t _UInt64;
		double _Double;
		float _Float;
		std::string *_String;
	} data;

private:
	void copy(const Variant &a);
	void check(const Type t) const;
	void check(const Type t);
};

#define VARIANT_TO_DECL(NAME, VALUE) \
	template<> inline VALUE Variant::to<VALUE>() const { return to ## NAME(); } \

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
VARIANT_TO_DECL(String, std::string)
VARIANT_TO_DECL(Double, double)

std::ostream& operator <<(std::ostream& os, const Variant &v);

} // namespace crpropa 

#endif // VARIANT_HH
