#ifndef LOGSTREAM_H
#define LOGSTREAM_H

#include <cstdio>
#include <sstream>
#include <iomanip>

// LOG output and file handler ///////
class LogStream
{

public:
    FILE *flog;
    const char *fname;
    int default_precision;
    int precision;
    void setprecision(int n);
    void resetprecision();
    LogStream(const char *fname);
    ~LogStream(void);
};

template <class T>
LogStream& operator<< (LogStream& st, T val)
{
  std::stringstream ss;

  ss << std::setprecision(st.precision) << std::fixed << val;
  fprintf(st.flog, "%s", ss.str().c_str());
  std::cout<<ss.str();
  return st;
}

#endif // LOGSTREAM_H
