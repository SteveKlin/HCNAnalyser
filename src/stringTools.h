#ifndef STRINGTOOLS_H
#define STRINGTOOLS_H

#include <string>

typedef std::string streng;

streng formatString(streng format, ...);

streng intToStr(long i);
streng floatToStr(double f);

#endif // STRINGTOOLS_H
