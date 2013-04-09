#include <stdio.h>
#include <stdlib.h>
#include "dint/error.h"

FILE* SafeFOpen(const char* filename, const char* mode)
{
    FILE* file;
    if ((file = fopen(filename, mode)) == NULL)
    {
	printf("SafeFOpen: cannot open %s\n", filename);
	exit (IO_ERROR);
    }

    return file;
}
