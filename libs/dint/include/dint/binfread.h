#ifndef DINT__BINFREAD_H
#define DINT__BINFREAD_H

// T. Beau 2005

// Allow to deal with big / little endian machines.
// Data have been written on a Intel Machine...

#include <stdio.h>
#include <stdlib.h>
#include <iostream>

#if HAVE_ARPA_INET_H
#include <arpa/inet.h>
#endif

#if HAVE_NETINET_IN_H
#include <netinet/in.h>
#endif 

size_t binfread(void *, size_t, size_t, FILE *);

#endif

