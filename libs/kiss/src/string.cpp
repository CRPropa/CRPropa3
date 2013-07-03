#include "kiss/string.h"

#include <iostream>
#include <string>
#include <vector>

namespace kiss {

std::string trim_right(const std::string &s, const std::string &t) {
	std::string::size_type i(s.find_last_not_of(t));

	if (i == std::string::npos)
		return "";
	else
		return std::string(s, 0, i);
}

std::string trim_left(const std::string &s, const std::string &t) {
	return std::string(s, s.find_first_not_of(t));
}

std::string trim(const std::string &s, const std::string &t) {
	std::string::size_type a = s.find_first_not_of(t);
	std::string::size_type b = s.find_last_not_of(t);

	if (a == std::string::npos || b == std::string::npos)
		return "";

	return std::string(s, a, b - a + 1);
}

void explode(const std::string &s, std::vector<std::string> &v,
		const bool trim_spaces, const std::string &t) {
	std::string::size_type a, b;

	a = s.find_first_not_of(t), b = s.find_first_of(t, a);

	while (a != std::string::npos) {
		if (trim_spaces)
			v.push_back(trim(s.substr(a, b - a)));
		else
			v.push_back(s.substr(a, b - a));

		a = s.find_first_not_of(t, b), b = s.find_first_of(t, a);
	}
}

std::string implode(const std::vector<std::string> &v, const std::string &t) {
	unsigned int i;
	std::string s;

	for (i = 0; i < (v.size() - 1); i++) {
		s.append(v[i]);
		s.append(t);
	}

	return s + v[i];
}

bool ends_with(const std::string &s, const std::string &w) {
	if (s.size() < w.size())
		return false;
	std::string::const_reverse_iterator si = s.rbegin();
	std::string::const_reverse_iterator wi = w.rbegin();
	while (wi != w.rend()) {
		if (*wi != *si)
			return false;
		wi++;
		si++;
	}

	return true;
}

bool starts_with(const std::string &s, const std::string &w) {
	if (s.size() < w.size())
		return false;
	std::string::const_iterator si = s.begin();
	std::string::const_iterator wi = w.begin();
	while (wi != w.end()) {
		if (*wi != *si)
			return false;
		wi++;
		si++;
	}
	return true;
}

} // namespace kiss
