/******************************************************************************
 *
 * MantaFlow fluid solver framework
 * Copyright 2011 Tobias Pfaff, Nils Thuerey 
 *
 * This program is free software, distributed under the terms of the
 * GNU General Public License (GPL) 
 * http://www.gnu.org/licenses
 *
 * Base class for particle systems
 *
 ******************************************************************************/

#ifndef _PYTHONINCLUDE_H
#define _PYTHONINCLUDE_H

#ifndef WIN32
#include <Python.h>
#else

// special handling for windows
// disable linking with debug version of python libs
#ifndef MANTA_WIN32_DEBUG
#include <Python.h>
#else

// note - we have to include these first!
#include <string>
#include <vector>
#include <iostream>

#undef _DEBUG
#define NDEBUG
#include <Python.h>
#define _DEBUG
#undef NDEBUG

#endif


#endif

#endif