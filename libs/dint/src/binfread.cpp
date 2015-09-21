
#include "dint/binfread.h"

#include <netinet/in.h>

static int ToSwap() {
  static int ToSwap = -1 ;
  if ( ToSwap != -1 ) goto endToSwap ;
  else {
     if ( ntohs(1234) != ntohs(1234) ) ToSwap = 1;
     else ToSwap = 0;
  }

endToSwap:
  //cerr << "ToSwap = " << ToSwap << endl;
  return ToSwap ;
}

static int16_t Swap2(int16_t a) {
  return ( ((a >> 8)&0xff) | ((a << 8)&0xff00) ) ;
}

static int32_t Swap4(int32_t a) {
  return ( ((a >> 24)&0xff)      | ((a >> 8 )&0xff00 ) |
                   ((a << 8 )&0xff0000 ) | ((a << 24)&0xff000000) ) ;
}

static int64_t Swap8(int64_t a) {
  return ( ((a >> 56)&0xff)            | ((a >> 40 )&0xff00 ) |
           ((a >> 24)&0xff0000)        | ((a >> 8 )&0xff000000) |
           ((a << 8 )&(((int64_t)(0xff))<<32))   | 
           ((a << 24)&(((int64_t)(0xff))<<40))   |
           ((a << 40)&(((int64_t)(0xff))<<48))   |
           ((a << 56)&(((int64_t)(0xff))<<56)) ) ;
}

size_t binfread(void *ptr, size_t size, size_t nmemb, FILE *stream) {
  size_t i;
  size_t s = fread(ptr, size, nmemb, stream);

  if ( ToSwap() && size>1 ) {
    if ( size == 2 ) {
      int16_t *temp = (int16_t*)ptr;
      for ( i=0; i<nmemb; i++) 
        temp[i] = Swap2(temp[i]);
    }
    else if ( size == 4 ) {
      int32_t *temp = (int32_t*)ptr;
      for ( i=0; i<nmemb; i++) 
        temp[i] = Swap4(temp[i]);
    }
    else if ( size == 8 ) {
      int64_t *temp = (int64_t*)ptr;
      for ( i=0; i<nmemb; i++)
        temp[i] = Swap8(temp[i]);
    }
    else {
      fprintf(stderr,"DINT / binfread : Strange member size = %i. Exit", size);
      exit(1);
    }   
  }

  return s;
}
