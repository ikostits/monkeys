/* 
 * File:   errors.cc
 * Author: irina
 * 
 * Created on May 7, 2014, 4:41 PM
 */

#include "errors.h"

std::string error_to_str(error_num err) {
    switch (err) {
        case ERR_IS_NOT_ZERO:
            return "ERR_IS_NOT_ZERO";
        case ERR_IS_ZERO:
            return "ERR_IS_ZERO";
        case ERR_LESS_THAN_ZERO:
            return "ERR_LESS_THAN_ZERO";
        default:
            return "UNKNOWN";
    }
}
