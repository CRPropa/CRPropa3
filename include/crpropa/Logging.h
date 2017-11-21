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

// A simple file logger
class LogFile
{
  std::ofstream of;
  public:
    LogFile(std::string _name)
    {
      of.open(_name);
      of << "CRPropa Version: " << g_GIT_REFSPEC << std::endl;
      of << "\n";
    }

    std::ostream& getStream()
    {
      return of;
    }
};

//void setLogStream(std::ofstream *stream)
//{
//  kiss::Logger::setLogStream(stream);
//}

//void setLogLevel(kiss::eLogLevel level)
//{
//  kiss::Logger::setLogLevel(level);
//}

#endif // LOGGING_H
