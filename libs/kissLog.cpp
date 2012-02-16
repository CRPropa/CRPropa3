#include "kissLog.h"

#include <stdlib.h>
#include <iostream>

namespace kiss {

std::ostream *Log::stream = &std::cerr;
eLogLevel Log::level = LOG_LEVEL_WARNING;
const char* sLogLevel[] = { "ERROR  ", "WARNING", "INFO  ", "DEBUG " };

class EnvLog {
public:
	EnvLog() {
		Log::loadEnvLogLevel();
	}
};
static EnvLog _env_log_;

Log::Log(eLogLevel level) {
	time_t rawtime;
	struct tm * timeinfo;
	char buffer[80];

	time(&rawtime);
	timeinfo = localtime(&rawtime);

	strftime(buffer, 80, "%Y-%m-%d %H:%M:%S ", timeinfo);
	*stream << buffer;
	*stream << "[" << sLogLevel[level] << "] ";
}

Log::~Log() {
	*stream << std::endl;
}

std::ostream &Log::getLogStream() {
	return (*stream);
}

void Log::setLogStream(std::ostream *s) {
	stream = s;
}

void Log::setLogStream(std::ostream &s) {
	stream = &s;
}

void Log::setLogLevel(eLogLevel l) {
	level = l;
}
eLogLevel Log::getLogLevel() {
	return (level);
}

void Log::loadEnvLogLevel() {
	if (::getenv("KISS_LOG_LEVEL")) {
		int level = atoi(::getenv("KISS_LOG_LEVEL"));
		switch (level) {
		case LOG_LEVEL_ERROR:
			Log::setLogLevel(LOG_LEVEL_ERROR);
			break;
		case LOG_LEVEL_WARNING:
			Log::setLogLevel(LOG_LEVEL_WARNING);
			break;
		case LOG_LEVEL_INFO:
			Log::setLogLevel(LOG_LEVEL_INFO);
			break;
		case LOG_LEVEL_DEBUG:
			Log::setLogLevel(LOG_LEVEL_DEBUG);
			break;
		default:
			std::cerr << "kiss::Log: unknown log level in KISS_LOG_LEVEL '"
					<< level << " values from 0-3 expected." << std::endl;
			break;
		}
	}
}

} // namespace kiss

