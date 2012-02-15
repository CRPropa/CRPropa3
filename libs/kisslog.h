#ifndef KISSLOG_H_
#define KISSLOG_H_

#include <iostream>

namespace kiss {

enum eLogLevel {
	LOG_LEVEL_ERROR, LOG_LEVEL_WARNING, LOG_LEVEL_INFO, LOG_LEVEL_DEBUG
};

class Log {
	static std::ostream *stream;
	static eLogLevel level;
public:
	Log(eLogLevel level);
	~Log();
	static std::ostream &getLogStream();
	static void setLogStream(std::ostream *s);
	static void setLogStream(std::ostream &s);

	static void setLogLevel(eLogLevel level);
	static eLogLevel getLogLevel();

	static void loadEnvLogLevel();

	operator std::ostream &() {
		return getLogStream();
	}

	template<typename T> inline Log& operator<<(T& data) {
		getLogStream() << data;
		return *this;
	}

	inline Log& operator<<(std::ostream& (*func)(std::ostream&)) {
		getLogStream() << func;
		return *this;
	}
};

} // namespace kiss

#define KISS_LOG_ERROR if (kiss::Log::getLogLevel() < kiss::LOG_LEVEL_ERROR) {} else kiss::Log(kiss::LOG_LEVEL_ERROR)
#define KISS_LOG_WARING if (kiss::Log::getLogLevel() < kiss::LOG_LEVEL_WARNING) {} else kiss::Log(kiss::LOG_LEVEL_WARNING)
#define KISS_LOG_INFO if (kiss::Log::getLogLevel() < kiss::LOG_LEVEL_INFO) {} else kiss::Log(kiss::LOG_LEVEL_INFO)
#define KISS_LOG_DEBUG if (kiss::Log::getLogLevel() < kiss::LOG_LEVEL_DEBUG) {} else kiss::Log(kiss::LOG_LEVEL_DEBUG)

#endif /* KISSLOG_H_ */
