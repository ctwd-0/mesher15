/////////////////////////////////////////////////////////////////////////////
// linking_pragmas_static.h

// This file is specific to Micrsoft's compiler.

#pragma once

#if defined(ON_DLL_IMPORTS) || defined(ON_DLL_EXPORTS)
#error This file contains STATIC LIBRARY linking pragmas.
#endif

#if defined(WIN64)

// x64 (64 bit) static libraries

#if defined(NDEBUG)

// Release x64 (64 bit) libs
#pragma message( " --- openNURBS Release x64 (64 bit) build." )
#pragma comment(lib, "opennurbs/x64r/zlib.lib")
#pragma comment(lib, "opennurbs/x64r/opennurbs_staticlib.lib")

#else // _DEBUG

// Debug x64 (64 bit) libs
#pragma message( " --- openNURBS Debug x64 (64 bit) build." )
#pragma comment(lib, "opennurbs/x64d/zlib.lib")
#pragma comment(lib, "opennurbs/x64d/opennurbs_staticlib.lib")

#endif // NDEBUG else _DEBUG

#else // WIN32

// x86 (32 bit) static libraries

#if defined(NDEBUG)

// Release x86 (32 bit) libs
#pragma message( " --- openNURBS Release x86 (32 bit) build." )
#pragma comment(lib, "opennurbs/x86r/zlib.lib")
#pragma comment(lib, "opennurbs/x86r/opennurbs_staticlib.lib")

#else // _DEBUG

// Debug x86 (32 bit) libs
#pragma message( " --- openNURBS Debug x86 (32 bit) build." )
#pragma comment(lib, "opennurbs/x86d/zlib.lib")
#pragma comment(lib, "opennurbs/x86d/opennurbs_staticlib.lib")

#endif // NDEBUG else _DEBUG

#endif // WIN64 else WIN32