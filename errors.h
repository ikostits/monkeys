/* 
 * File:   errors.h
 * Author: irina
 *
 * Created on May 7, 2014, 4:41 PM
 */

#ifndef ERRORS_H
#define	ERRORS_H

#include <boost/exception/all.hpp>
#include <string>

enum error_num {
    ERR_IS_NOT_ZERO,
    ERR_IS_ZERO,
    ERR_LESS_THAN_ZERO
};

struct tag_errno_code { };
struct tag_code_line { };
struct tag_err_description { };
typedef boost::error_info<tag_errno_code,error_num> errno_code;
typedef boost::error_info<tag_code_line,int> code_line;
typedef boost::error_info<tag_err_description,std::string> err_description;

struct exception_base: virtual std::exception, virtual boost::exception { };
struct integral_error: virtual exception_base { };
struct geometry_error: virtual exception_base { };

std::string error_to_str(error_num err);

#endif	/* ERRORS_H */

