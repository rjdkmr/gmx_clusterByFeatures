#include <cstdio>
#include <sstream>
#include <iostream>
#include <iomanip>

#include "gromacs/utility/futil.h"

#include "logstream.h"


LogStream::LogStream(const char *fname) {
    this->default_precision = std::cout.precision();
    this->precision = this->default_precision;
    this->flog = gmx_ffopen(fname, "w");
}

LogStream::~LogStream() {
     gmx_ffclose(this->flog);
}

void LogStream::setprecision(int n){
    this->precision = n;

}

void LogStream::resetprecision(){
    this->precision = this->default_precision;
}

