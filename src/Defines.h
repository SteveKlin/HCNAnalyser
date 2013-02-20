#ifndef DEFINES_H
#define DEFINES_H

#define ensure(a) { if (!(a)) qDebug("Ensure (%s) failed!", #a); }

template <class D> D min(D a, D b) { return (a < b) ? a : b; }
template <class D> D max(D a, D b) { return (a > b) ? a : b; }

#endif // DEFINES_H
