#ifndef KISS_IO_H
#define KISS_IO_H

#include <iostream>
#include <vector>
#include <map>

#include <stdint.h>

namespace kiss {

class Output {
protected:
	bool failed;
public:
	Output() :
			failed(false) {

	}
	virtual ~Output() {
	}
	virtual bool writeBytes(const char *data, size_t size) = 0;

	operator bool() {
		return (!failed);
	}
};

class StreamOutput: public Output {
	std::ostream &output;
public:
	StreamOutput(std::ostream &output) :
			output(output) {

	}

	bool writeBytes(const char *data, size_t size) {
		output.write(data, size);
		failed = output.fail();
		return (output);
	}

};

class CompressedOutput {
	Output &output;
	std::vector<char> buffer;
	size_t buffer_size;
	size_t buffer_position;
public:
	CompressedOutput(Output &output, size_t buffer_size = (1024 * 1024)) :
			output(output), buffer(buffer_size), buffer_size(buffer_size), buffer_position(
					0) {
	}

//	bool writeBytes(char *data, size_t size) {
//		data_left = size;
//		while (data_left) {
//
//			output.write(data, size);
//		}
//		return (output);
//	}
};

template<class T>
void write(Output &out, const T &value) {
	out.writeBytes((char *) &value, sizeof(value));
}

inline
void write(Output &out, const std::string &s) {
	write(out, s.size());
	out.writeBytes(s.data(), s.size());
}

template<class KEY, class VALUE>
inline void write(Output &out, const std::map<KEY, VALUE> &m) {
	write(out, (uint64_t) m.count());
	typename std::map<KEY, VALUE>::const_iterator i = m.begin();
	while (i != m.end()) {
		write(out, i->first);
		write(out, i->second);
		i++;
	}
}

template<class VALUE>
inline void write(Output &out, const std::vector<VALUE> &m) {
	write(out, (uint64_t) m.size());
	typename std::vector<VALUE>::const_iterator i = m.begin();
	while (i != m.end()) {
		write(out, *i);
		i++;
	}
}

class Input {
protected:
	bool failed;
public:
	virtual ~Input() {
	}
	virtual bool readBytes(char *data, size_t size) = 0;

	operator bool() {
		return (!failed);
	}
};

class StreamInput: public Input {
	std::istream &input;
public:
	StreamInput(std::istream &input) :
			input(input) {

	}

	bool readBytes(char *data, size_t size) {
		input.read(data, size);
		failed = input.fail();
		return (input);
	}

};

template<class T>
inline bool read(Input &in, T &value) {
	in.readBytes((char *) &value, sizeof(value));
	return (in);
}

template<class T>
inline T read(Input &in) {
	T t;
	in.readBytes((char *) &t, sizeof(t));
	return t;
}

inline bool read(Input &in, std::string &s) {
	s.resize(read<size_t>(in));
	bool result = in.readBytes((char *) &s[0], s.size());
	return result;
}

template<class KEY, class VALUE>
inline bool read(Input &in, std::map<KEY, VALUE> &m) {
	size_t count = read<uint64_t>(in);
	for (size_t i = 0; i < count; i++) {
		KEY key;
		if (!read(in, key))
			return false;
		VALUE &v = m[key];
		if (!read(in, v))
			return false;
	}
	return true;
}

template<class VALUE>
inline bool read(Input &in, std::vector<VALUE> &m) {
	size_t count = read<uint64_t>(in);
	m.resize(count);
	for (size_t i = 0; i < count; i++) {
		VALUE &v = m[i];
		read(in, v);
	}
	return (in);
}

} // namespace kiss

#endif /* KISS_IO_H */
