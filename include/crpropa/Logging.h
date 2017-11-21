#ifndef LOGGING_H
#define LOGGING_H
#include <crpropa/Version.h>

#include <kiss/logger.h>

#include <fstream>

// make the kiss log functions available in python
void inline logError(const std::string &log)
{
  KISS_LOG_ERROR << log;
}

void inline logInfo(const std::string &log)
{
  KISS_LOG_INFO << log;
}

void inline logWarning(const std::string &log)
{
  KISS_LOG_WARNING << log;
}

void inline logDebug(const std::string &log)
{
  KISS_LOG_DEBUG << log;
}

void setLogStream(std::ostream &stream)
{
  kiss::Logger::setLogStream(stream);
}

void setLogLevel(int level)
{
  kiss::Logger::setLogLevel(static_cast<kiss::eLogLevel>(level));
}

#endif // LOGGING_H
