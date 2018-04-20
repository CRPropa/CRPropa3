#ifndef KISS_CONVERT_H
#define KISS_CONVERT_H

#include <iostream>
#include <sstream>
#include <string>
#include <typeinfo>
#include <stdexcept>

namespace kiss {

class bad_conversion: public std::runtime_error {
public:
	bad_conversion(std::string const& s) :
			std::runtime_error(s) {
	}
};

template<typename T>
inline std::string str(T const& x) {
	std::ostringstream o;
	o.imbue(std::locale("C"));
	o << x;
	if (!o) {
#ifdef DEBUG
		std::cerr << "bad_conversion: str(" << typeid(x).name() << ")";
#endif
		throw bad_conversion(std::string("str(") + typeid(x).name() + ")");
	}
	return o.str();
}

template<typename T>
inline void convert(const char *s, T& x) {
	if (s == 0) {
#ifdef DEBUG
		std::cerr << "bad_conversion: convert, null pointer";
#endif
		throw bad_conversion("null pointer");
	}
	std::istringstream i(s);
	i.imbue(std::locale("C"));
	i >> x;
	if (!i) {
#ifdef DEBUG
		std::cerr << "bad_conversion: convert (" << s << ")";
#endif
		throw bad_conversion(std::string(s) + " to " + typeid(x).name());
	}
}

template<typename T>
inline void convert(std::string const& s, T& x) {
	std::istringstream i(s);
	i.imbue(std::locale("C"));
	i >> x;
	if (!i) {
#ifdef DEBUG
		std::cerr << "bad_conversion: convert (" << s << ")";
#endif
		throw bad_conversion(s + " to " + typeid(x).name());
	}
}

template<typename T>
inline T convertTo(std::string const& s, bool failIfLeftoverChars = true) {
	T x;
	convert(s, x, failIfLeftoverChars);
	return x;
}

} // namespace kiss

#endif /* KISSCONVERT_H */
