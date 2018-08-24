/* name of the plugin: myPlugin*/
%module(directors="1", threads="1", allprotected="1") myPlugin

/* Exceptions required */
%include "exception.i"

/*  define headers to include into the wrapper. These are the plugin headers
 *  and the CRPRopa headers.
 */
%{
#include "CRPropa.h"
#include "myPlugin.h"
%}

/* import crpropa in wrapper */
%import (module="crpropa") "crpropa.i"

/* include plugin parts to generate wrappers for */
%include "myPlugin.h"





