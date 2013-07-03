#ifndef KISS_STRING_H
#define KISS_STRING_H

#include <string>
#include <map>
#include <vector>

namespace kiss {

#define SPACES " \t\r\n"

std::string trim_right(const std::string &s, const std::string &t);
std::string trim_left(const std::string &s, const std::string &t);

std::string trim(const std::string &s, const std::string &t = SPACES);
void explode(const std::string &s, std::vector<std::string> &v,
		const bool trim_spaces, const std::string &t);
std::string implode(const std::vector<std::string> &v, const std::string &t);
bool ends_with(const std::string &s, const std::string &w);
bool starts_with(const std::string &s, const std::string &w);

} // namespace kiss

#endif /* KISS_STRING_H */
