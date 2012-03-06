#include "kiss/uuid.h"

#ifdef _MSVC_VER
#define _CRT_RAND_S
#include <stdlib.h>
#else
#include <stdio.h>
#include <vector>
#endif

namespace kiss {

#define KISS_UUID_CREATE_BUFFER_SIZE 1000

uuid::uuid() {
	ints[0] = 0;
	ints[1] = 0;
	ints[2] = 0;
	ints[3] = 0;
}

uuid uuid::create() {
	unsigned int ints[4];
#ifdef _MSVC_VER
	for (int i = 0; i < 4; i++) {
		rand_s(ints[0]);
		rand_s(ints[1]);
		rand_s(ints[2]);
		rand_s(ints[3]);
	}
#elif _WIN32
#error "not implemented"
#else
	static std::vector<uuid> buffer;
	if (buffer.size() == 0) {
		buffer.resize(KISS_UUID_CREATE_BUFFER_SIZE);
		FILE *rd = fopen("/dev/urandom", "rb");
		size_t n = fread((void *) buffer.data(), sizeof(uuid),
				KISS_UUID_CREATE_BUFFER_SIZE, rd);
		fclose(rd);
		buffer.resize(n);
	}
	ints[0] = buffer.back().ints[0];
	ints[1] = buffer.back().ints[1];
	ints[2] = buffer.back().ints[2];
	ints[3] = buffer.back().ints[3];
	buffer.pop_back();
#endif

	unsigned char *bytes = (unsigned char *) &ints;
	/* set version 4 (random)*/
	bytes[6] = (4 << 4) + (bytes[6] & 0x0F);

	/* set variant (always DCE 1.1 only) */
	bytes[8] |= 0x80;
	bytes[8] &= ~0x40;

	return uuid(ints);
}

uuid::uuid(unsigned int *p) {
	ints[0] = p[0];
	ints[1] = p[1];
	ints[2] = p[2];
	ints[3] = p[3];
}

uuid::uuid(unsigned int a, unsigned int b, unsigned int c, unsigned int d) {
	ints[0] = a;
	ints[1] = b;
	ints[2] = c;
	ints[3] = d;
}

bool uuid::operator ==(const uuid& id) const {
	for (int i = 0; i < 4; i++) {
		if (ints[i] != id.ints[i])
			return false;
	}

	return true;
}

bool uuid::operator !=(const uuid& id) const {
	for (int i = 0; i < 4; i++) {
		if (ints[i] != id.ints[i])
			return true;
	}

	return false;
}

bool uuid::operator <(const uuid& op) const {
	for (int i = 0; i < 4; i++) {
		if (ints[i] < op.ints[i])
			return true;
		else if (ints[i] > op.ints[i])
			return false;
	}

	return false;
}

uuid uuid::parse(const char* id) {
	unsigned int ints[4];
	unsigned char *bytes = (unsigned char *) &ints;

	// reset
	for (size_t j = 0; j < 4; j++)
		ints[j] = 0;

	// read 32 char = 16 bytes
	unsigned char first = 0, second = 0;
	int byte = 0;
	const char* source = id;
	unsigned char* target = &first;

	while (*source != 0 && byte < 16) {
		// find next a valid character
		if (*source >= '0' && *source <= '9') {
			*target = *source - '0';
		} else if (*source >= 'a' && *source <= 'f') {
			*target = *source - 'a' + 10;
		} else if (*source >= 'A' && *source <= 'F') {
			*target = *source - 'A' + 10;
		} else {
			source++;
			continue;
		}

		//
		if (target == &first)
			target = &second;
		else {
			bytes[byte] = ((first << 4) | second);
			byte++;
			target = &first;
		}

		source++;
	}

	return uuid(ints);
}

uuid uuid::parse(const std::string& id) {
	unsigned char first, second;
	int byte = 0;
	unsigned char* target = &first;
	std::string::const_iterator source = id.begin();
	std::string::const_iterator end = id.end();
	unsigned int ints[4] = { 0, 0, 0, 0 };
	unsigned char *bytes = (unsigned char *) &ints;

	while (source != end && byte < 16) {
		if (*source >= '0' && *source <= '9') {
			*target = *source - '0';
		} else if (*source >= 'a' && *source <= 'f') {
			*target = *source - 'a' + 10;
		} else if (*source >= 'A' && *source <= 'F') {
			*target = *source - 'A' + 10;
		} else {
			source++;
			continue;
		}

		if (target == &first)
			target = &second;
		else {
			bytes[byte] = ((first << 4) | second);
			byte++;
			target = &first;
		}

		source++;
	}

	return uuid(ints);
}

void uuid::print(char *p) const {
	const char hex[] = { '0', '1', '2', '3', '4', '5', '6', '7', '8', '9', 'a',
			'b', 'c', 'd', 'e', 'f' };
	unsigned char *bytes = (unsigned char *) &ints;

	for (int i = 0; i < 4; i++) {
		*p++ = hex[bytes[i] >> 4];
		*p++ = hex[bytes[i] & 0x0F];
	}
	*p++ = '-';
	for (int i = 4; i < 6; i++) {
		*p++ += hex[bytes[i] >> 4];
		*p++ += hex[bytes[i] & 0x0F];
	}
	*p++ = '-';
	for (int i = 6; i < 8; i++) {
		*p++ += hex[bytes[i] >> 4];
		*p++ += hex[bytes[i] & 0x0F];
	}
	*p++ = '-';
	for (int i = 8; i < 10; i++) {
		*p++ += hex[bytes[i] >> 4];
		*p++ += hex[bytes[i] & 0x0F];
	}
	*p++ = '-';
	for (int i = 10; i < 16; i++) {
		*p++ += hex[bytes[i] >> 4];
		*p++ += hex[bytes[i] & 0x0F];
	}

	*p++ = 0;
}

bool uuid::valid() const {
	return !(ints[0] == 0 && ints[1] == 0 && ints[2] == 0 && ints[3] == 0);
}

unsigned int &uuid::operator[](size_t i) {
	return ints[i];
}

const unsigned int &uuid::operator[](size_t i) const {
	return ints[i];
}

const unsigned char *uuid::data() const {
	return (unsigned char *) &ints;
}

} // namespace kiss

std::ostream& operator <<(std::ostream& os, const kiss::uuid &id) {
	const char hex[] = { '0', '1', '2', '3', '4', '5', '6', '7', '8', '9', 'a',
			'b', 'c', 'd', 'e', 'f' };

	const unsigned char *bytes = id.data();

	for (int i = 0; i < 4; i++) {
		os << hex[bytes[i] >> 4];
		os << hex[bytes[i] & 0x0F];
	}
	os << '-';
	for (int i = 4; i < 6; i++) {
		os << hex[bytes[i] >> 4];
		os << hex[bytes[i] & 0x0F];
	}
	os << '-';
	for (int i = 6; i < 8; i++) {
		os << hex[bytes[i] >> 4];
		os << hex[bytes[i] & 0x0F];
	}
	os << '-';
	for (int i = 8; i < 10; i++) {
		os << hex[bytes[i] >> 4];
		os << hex[bytes[i] & 0x0F];
	}
	os << '-';
	for (int i = 10; i < 16; i++) {
		os << hex[bytes[i] >> 4];
		os << hex[bytes[i] & 0x0F];
	}

	return os;
}

unsigned int _char_to_byte(unsigned char c) {
	if (c >= '0' && c <= '9') {
		return c - '0';
	} else if (c >= 'a' && c <= 'f') {
		return c - 'a' + 10;
	} else if (c >= 'A' && c <= 'F') {
		return c - 'A' + 10;
	}
	return c;
}

std::istream& operator >>(std::istream& os, kiss::uuid &id) {
	unsigned int ints[4] = { 0, 0, 0, 0 };
#if 0
	unsigned char *bytes = id.ints;

	// read 32 char = 16 bytes
	unsigned char first = 0, second = 0;
	for (size_t i = 0; i < 4; i++) {
		os >> first;
		os >> second;
	}
	int byte = 0;
	const char* source = id;
	unsigned char* target = &first;

	while (*source != 0 && byte < 16) {
		// find next a valid character
		if (*source >= '0' && *source <= '9') {
			*target = *source - '0';
		} else if (*source >= 'a' && *source <= 'f') {
			*target = *source - 'a' + 10;
		} else if (*source >= 'A' && *source <= 'F') {
			*target = *source - 'A' + 10;
		} else {
			source++;
			continue;
		}

		//
		if (target == &first)
		target = &second;
		else {
			bytes[byte] = ((first << 4) | second);
			byte++;
			target = &first;
		}

		source++;
	}
#endif
	std::cerr << "not implemented" << std::endl;
	return os;
}
