#include "mswrapper.h"

#include <QDebug>
/*
char **
cmatrix(int nsam, int len) {
int i;
char **m;

if( ! ( m = new char*[nsam] ) ) qDebug("alloc error in cmatrix") ;

for( i=0; i<nsam; i++) {
    if( ! ( m[i] = new char[len])) qDebug("alloc error in cmatric. 2");
}
return( m );
}

void free_cmatrix(char** list, int nsam) {
    unsigned int q;
    for (q = 0; q < nsam; ++q) {
        delete [] list[nsam];
    }
    delete [] list;
}
*/
