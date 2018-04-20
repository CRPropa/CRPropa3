#ifndef FILEUTILS_INCLUDED
#define FILEUTILS_INCLUDED

#include <string>
#include <vector>

struct mode {
	enum {
		none = 0, execute = 1, write = 2, read = 4, rw = 6, rx = 5, rwx = 7
	};
};

#ifdef _WIN32
const char path_seperator = '\\';
#else
const char path_seperator = '/';
#endif

bool create_directory(const std::string &path, size_t user_permission =
		mode::rwx, size_t group_permission = mode::rx, size_t other_permission =
		mode::none);

bool is_directory(const std::string &path);

bool list_directory(const std::string &directory,
		std::vector<std::string> &elements);

bool create_directory_recursive(const std::string &dir, size_t user_permission =
		mode::rwx, size_t group_permission = mode::rx, size_t other_permission =
		mode::none);
std::string concat_path(const std::string &a, const std::string &b);

std::string concat_path(const std::string &a, const std::string &b,
		const std::string &c);

void append_file(const std::string &target, const std::string &source, bool binary);

std::string executable_path();

#endif /* FILEUTILS_H */
