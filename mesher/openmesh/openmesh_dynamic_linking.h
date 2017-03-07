/////////////////////////////////////////////////////////////////////////////
// linking_pragmas_static.h

// This file is specific to Micrsoft's compiler.

#pragma once

#if !defined(ON_DLL_IMPORTS)
#error This file contains DYNAMIC LIBRARY linking pragmas.
#endif

#if defined(WIN64)

// x64 (64 bit) dynamic libraries

#if defined(NDEBUG)

// Release x64 (64 bit) libs
#pragma message( " --- openmesh Release x64 (64 bit) build." )
#pragma comment(lib, "openmesh/lib/OpenMeshCore.lib")
#pragma comment(lib, "openmesh/lib/OpenMeshTools.lib")

#else // _DEBUG

// Debug x64 (64 bit) libs
#pragma message( " --- openmesh Debug x64 (64 bit) build." )
#pragma comment(lib, "openmesh/lib/OpenMeshCored.lib")
#pragma comment(lib, "openmesh/lib/OpenMeshToolsd.lib")

#endif // NDEBUG else _DEBUG

#else // WIN32

// x86 (32 bit) dynamic libraries

#if defined(NDEBUG)

// Release x86 (32 bit) libs
#pragma message( " --- openmesh Release x86 (32 bit) build." )
#pragma comment(lib, "openmesh/lib/OpenMeshCore.lib")
#pragma comment(lib, "openmesh/lib/OpenMeshTools.lib")

#else // _DEBUG

// Debug x86 (32 bit) libs
#pragma message( " --- openmesh Debug x86 (32 bit) build." )
#pragma comment(lib, "openmesh/lib/OpenMeshCored.lib")
#pragma comment(lib, "openmesh/lib/OpenMeshToolsd.lib")

#endif // NDEBUG else _DEBUG

#endif // WIN64 else WIN32