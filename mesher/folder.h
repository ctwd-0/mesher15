#pragma once

#ifdef _WIN32
#include <direct.h>
#include <io.h>
#elif _LINUX
#include <stdarg.h>
#include <sys/stat.h>
#endif

#ifdef _WIN32
#define ACCESS _access
#define MKDIR(a) _mkdir((a))
#elif _LINUX
#define ACCESS access
#define MKDIR(a) mkdir((a), 0755)
#endif

inline bool make_dir(const char* dir) {
	int ret = ACCESS(dir, 0);
	if (ret != 0) {
		ret = MKDIR(dir);
		if (ret != 0) {
			return false;
		}
	}
	return true;
}