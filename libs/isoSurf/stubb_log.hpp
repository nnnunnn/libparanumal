// declare stubbs for missing MPI logging routines.

#pragma once

#ifndef _MSC_VER

#include <string>
#include <cstdio>

void  umLOG(int n, const char* format_str, ...);
void  umLOG(const std::string& msg, int n);
char* umSTR(const char* fmt, ...);

#define umMSG umLOG
#define umTRC umLOG

// NBN: global g_procid = comm.rank;
extern int g_procid;

#endif
