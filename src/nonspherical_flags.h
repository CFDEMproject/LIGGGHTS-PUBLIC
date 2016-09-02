#ifndef LMP_NONSPERICAL_FLAGS_H
#define LMP_NONSPERICAL_FLAGS_H

// comment in our out
//#define SUPERQUADRIC_ACTIVE_FLAG

// comment in our out
//#define NONSPHERICAL_ACTIVE_FLAG

#ifdef SUPERQUADRIC_ACTIVE_FLAG
    #ifndef NONSPHERICAL_ACTIVE_FLAG
        #define NONSPHERICAL_ACTIVE_FLAG
    #endif
#include "superquadric.h"
#endif

#endif

