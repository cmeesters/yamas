// $Rev: 622 $
// $LastChangedDate: 2010-12-18 15:42:06 +0100 (Sa, 18 Dez 2010) $
// $Author: meesters $

#include "ya_interface.hpp"

#ifndef YA_GLOBALS
#define YA_GLOBALS

#ifdef MAIN
    #define EXTERN
#else
    #define EXTERN extern
#endif

// yamas user data as kept in the UsageData struct
EXTERN yamas::UsageData ud;

#endif // YA_GLOBALS
