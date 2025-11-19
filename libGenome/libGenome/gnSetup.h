/////////////////////////////////////////////////////////////////////////////
// File:        libGenome/gnSetup.h
// Purpose:     libGenome setup
// Description: Defines os/compiler specific constants, etc.
//              Included in libGenome/gnDefs.h.
// Rev:         A
// Author:      Aaron Darling 
// Modified by:      
// Copyright:   (c) Aaron Darling 
// Licenses:      
/////////////////////////////////////////////////////////////////////////////
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#ifndef _gnSetup_h_
#define _gnSetup_h_

#include <stdlib.h>

#ifdef GNMAKINGDLL   // build the libgenome dll
#define GNDLLEXPORT __declspec(dllexport)
#define GNDLLEXPORT_DATA(type) __declspec(dllexport) type
#elif defined(GNUSINGDLL)  
// the project uses a libGenome as a dll
#define GNDLLEXPORT __declspec(dllimport)
#define GNDLLEXPORT_DATA(type) __declspec(dllimport) type
#else	// static linking
#define GNDLLEXPORT
#define GNDLLEXPORT_DATA(type) type
#endif

// Auto-linking pragmas removed to prevent conflict with CMake's target_link_libraries.
// CMake generates 'libGenome.lib', but these pragmas were forcing a search for 'genome.lib'.

#endif
	//_gnSetup_h_
