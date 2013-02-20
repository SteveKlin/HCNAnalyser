#include "stringTools.h"

#include <stdarg.h>
#include <cstdio>

streng formatString(streng format, ...) {
    va_list args;
    va_start(args, format);

    int buffSize = 1024;
    streng result = "";
    int r = -1;
    while ((r < 0) && (buffSize < 65536)) {
        char* buffer = new char[buffSize];
        r = vsnprintf(buffer, buffSize, format.c_str(), args);
        if (r >= 0) { result = buffer; }
        else { buffSize *= 2; }
        delete [] buffer;
    }

    va_end (args);

    return result;
}

streng intToStr(long i) { return formatString("%li", i); }
streng floatToStr(double f) { return formatString("%f", f); }
