#include <vector>

#include "tools.h"

int search( vector<int> v, int x, int start, int end )
{
    if ( start==-1 && end==-1 ) {
        start = 0;
        end = v.size()-1;
    }
    if ( start > end ) return -1;
    
    int i = int( (start+end)/2 );
    if (v[i] == x) return i;
    else if (start == end) return -1;
    else if (x < v[i]) {
        i = search( v,x,start,i-1 );
        if (i>=0) return i;
        else return -1;
    }
    else {
        i = search( v,x,i+1,end );
        if (i>=0) return i;
        else return -1;
    }
}

