#ifndef BASE64_H
#define BASE64_H

/// base64 encodig and decoding.
///
/// Based on the implementations by
/// Jouni Malinen <j@w1.fi> and contributors from wpa_supplicant and hostapd in
/// http://web.mit.edu/freebsd/head/contrib/wpa/ and
/// http://web.mit.edu/freebsd/head/contrib/wpa/src/utils/base64.c and
/// http://web.mit.edu/freebsd/head/contrib/wpa/src/utils/base64.h
///
/// Published under a 3-clause BSD license
///
#include <string>
#include <cstring>
#include <stdexcept>

namespace crpropa
{
	class Base64
	{
	public:
		static std::string encode(const unsigned char *src, size_t len);
		static std::string decode(const std::string &data);
	};
}


#endif // BASE64_H
