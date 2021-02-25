// Useful c-style macros
// Copyright (C) 2018  David Zuñiga-Noël <dzuniga at uma.es>

#ifndef MACROS_H_
#define MACROS_H_

// STL
#include <cstdio>
#include <iostream>
#include <stdexcept>
#include <string>

#define STRINGIFY(x) #x
#define TOSTRING(x) STRINGIFY(x)

#define __AT__ __FILE__ ":" TOSTRING(__LINE__)
#define RUNTIME_ASSERT(x) { if (!(x)) throw std::runtime_error("Runtime assertion failed (at " __AT__ "):\n" TOSTRING(x)); }

#define DO_PRAGMA(x) _Pragma (#x)
#define TODO(x) DO_PRAGMA(message ("TODO - " x))

#define USE(x) (void)(x);

#define COUT_LOG(x) { std::cout << std::string("[ ") + __PRETTY_FUNCTION__ + std::string(" ] ") + TOSTRING(x) + std::string(": ") << (x) << std::endl; }

#endif // MACROS_H_
