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

#if defined(WIN32) || defined(_WIN32)

	// note - we have to include these first!
	#include <string>
	#include <vector>
	#include <iostream>

	#if defined(_DEBUG) && defined(PYTHON_DEBUG_AS_RELEASE)

		// special handling for windows
		// disable linking with debug version of python libs
		#undef _DEBUG
		#define NDEBUG
		#include <Python.h>
		#define _DEBUG
		#undef NDEBUG

	#else
		#include <Python.h>
	#endif
#else
	#include <Python.h>
#endif

#endif