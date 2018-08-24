#include "kiss/logger.h"

#include <stdlib.h>
#include <iostream>

namespace kiss {

std::ostream *Logger::stream = &std::cerr;
eLogLevel Logger::level = LOG_LEVEL_WARNING;
const char* sLoggerLevel[] = { " ERROR ", "WARNING", " INFO  ", " DEBUG " };

class EnvLogger {
public:
	EnvLogger() {
		Logger::loadEnvLogLevel();
	}
};
static EnvLogger _env_log_;

Logger::Logger(eLogLevel level) {
	time_t rawtime;
	struct tm * timeinfo;
	char buffer[80];

	time(&rawtime);
	timeinfo = localtime(&rawtime);

	strftime(buffer, 80, "%Y-%m-%d %H:%M:%S ", timeinfo);
	*stream << buffer;
	*stream << "[" << sLoggerLevel[level] << "] ";
}

Logger::~Logger() {
	*stream << std::endl;
}

std::ostream &Logger::getLogStream() {
	return (*stream);
}

void Logger::setLogStream(std::ostream *s) {
	stream = s;
}

void Logger::setLogStream(std::ostream &s) {
	stream = &s;
}

void Logger::setLogLevel(eLogLevel l) {
	level = l;
}
eLogLevel Logger::getLogLevel() {
	return (level);
}




void Logger::loadEnvLogLevel() {
	if (::getenv("KISS_LOG_LEVEL")) {
		
		int level = atoi(::getenv("KISS_LOG_LEVEL"));
		switch (level) {
		case LOG_LEVEL_ERROR:
			Logger::setLogLevel(LOG_LEVEL_ERROR);
			break;
		case LOG_LEVEL_WARNING:
			Logger::setLogLevel(LOG_LEVEL_WARNING);
			break;
		case LOG_LEVEL_INFO:
			Logger::setLogLevel(LOG_LEVEL_INFO);
			break;
		case LOG_LEVEL_DEBUG:
			Logger::setLogLevel(LOG_LEVEL_DEBUG);
			break;
		default:
			std::cerr << "kiss::Logger: unknown log level in KISS_LOG_LEVEL '"
					<< level << " values from 0-3 expected." << std::endl;
			break;
		}
	}
}

} // namespace kiss

