#ifndef KISS_IO_H_
#define KISS_IO_H_

#include <iostream>

namespace kiss {

template<class T>
void write(std::ostream &out, const T &value) {
	out.write((char *) &value, sizeof(value));
}

inline
void write(std::ostream &out, const std::string &s) {
	write(out, s.size());
	out.write(s.data(), s.size());
}

template<class KEY, class VALUE>
inline void write(std::ostream &out, const std::map<KEY, VALUE> &m) {
	write(out, m.count());
	typename std::map<KEY, VALUE>::const_iterator i = m.begin();
	while (i != m.end()) {
		write(out, i->first);
		write(out, i->second);
		i++;
	}
}

template<class T>
inline bool read(std::istream &in, T &value) {
	in.read((char *) &value, sizeof(value));
	return (in);
}

template<class T>
inline T read(std::istream &in) {
	T t;
	in.read((char *) &t, sizeof(t));
	return t;
}

inline bool read(std::istream &in, std::string &s) {
	s.resize(read<size_t>(in));
	in.read((char *) &s[0], s.size());
	return (in);
}

template<class KEY, class VALUE>
inline bool read(std::istream &in, std::map<KEY, VALUE> &m) {
	size_t count = read<size_t>(in);
	for (size_t i = 0; i < count; i++) {
		VALUE &v = m[read<std::string>(in)];
		read(in, v);
		i++;
	}
	return (in);
}

} // namespace kiss

#endif /* KISS_IO_H_ */
