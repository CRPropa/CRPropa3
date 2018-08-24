#ifndef KISS_LOG_H
#define KISS_LOG_H

#include <iostream>

namespace kiss {

enum eLogLevel {
	LOG_LEVEL_ERROR, LOG_LEVEL_WARNING, LOG_LEVEL_INFO, LOG_LEVEL_DEBUG
};

class Logger {
	static std::ostream *stream;
	static eLogLevel level;
public:
	Logger(eLogLevel level);
	~Logger();
	static std::ostream &getLogStream();
	static void setLogStream(std::ostream *s);
	static void setLogStream(std::ostream &s);

	static void setLogLevel(eLogLevel level);
	static eLogLevel getLogLevel();

	static void loadEnvLogLevel();

	operator std::ostream &() {
		return getLogStream();
	}

	template<typename T> inline Logger& operator<<(T& data) {
		#pragma omp critical (KISS_LOGGER)
		{
		getLogStream() << data;
		}
		return *this;
	}

	inline Logger& operator<<(std::ostream& (*func)(std::ostream&)) {
		#pragma omp critical (KISS_LOGGER)
		{
		getLogStream() << func;
		}
		return *this;
	}
};

} // namespace kiss

#define KISS_LOG_ERROR if (kiss::Logger::getLogLevel() < kiss::LOG_LEVEL_ERROR) {} else kiss::Logger(kiss::LOG_LEVEL_ERROR)
#define KISS_LOG_WARNING if (kiss::Logger::getLogLevel() < kiss::LOG_LEVEL_WARNING) {} else kiss::Logger(kiss::LOG_LEVEL_WARNING)
#define KISS_LOG_INFO if (kiss::Logger::getLogLevel() < kiss::LOG_LEVEL_INFO) {} else kiss::Logger(kiss::LOG_LEVEL_INFO)
#define KISS_LOG_DEBUG if (kiss::Logger::getLogLevel() < kiss::LOG_LEVEL_DEBUG) {} else kiss::Logger(kiss::LOG_LEVEL_DEBUG)

#endif /* KISSLOG_H */
