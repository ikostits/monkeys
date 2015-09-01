//
//  precision.cc
//  monkeys
//
//  Created by Irina Kostitsyna on 6/11/14.
//
//

#include "precision.h"

mpfr_class logarithm(const mpfr_class& vt) {
    mpfr_class result;
    mpfr_log(result.var_, vt.var_, MPFR_RNDNA);

    return result;
}

mpfr_class power(const mpfr_class& vt, long i) {
    mpfr_class result;
    mpfr_pow_si(result.var_, vt.var_, i, MPFR_RNDNA);

    return result;
}

mpfr_class operator+(const mpfr_class& lhs, const mpfr_class& rhs) {
    mpfr_class result;
    mpfr_add(result.var_, lhs.var_, rhs.var_, MPFR_RNDNA);

    return result;
}

mpfr_class operator-(const mpfr_class& lhs, const mpfr_class& rhs) {
    mpfr_class result;
    mpfr_sub(result.var_, lhs.var_, rhs.var_, MPFR_RNDNA);

    return result;
}

mpfr_class operator*(const mpfr_class& lhs, const mpfr_class& rhs) {
    mpfr_class result;
    mpfr_mul(result.var_, lhs.var_, rhs.var_, MPFR_RNDNA);

    return result;
}

mpfr_class operator*(double lhs, const mpfr_class& rhs) {
    mpfr_class result;
    mpfr_mul_d(result.var_, rhs.var_, lhs, MPFR_RNDNA);

    return result;
}

mpfr_class operator*(const mpfr_class& lhs, double rhs) {
    mpfr_class result;
    mpfr_mul_d(result.var_, lhs.var_, rhs, MPFR_RNDNA);

    return result;
}

mpfr_class operator/(const mpfr_class& lhs, const mpfr_class& rhs) {
    mpfr_class result;
    mpfr_div(result.var_, lhs.var_, rhs.var_, MPFR_RNDNA);

    return result;
}

mpfr_class operator/(const mpfr_class& lhs, double rhs) {
    mpfr_class result;
    mpfr_div_d(result.var_, lhs.var_, rhs, MPFR_RNDNA);

    return result;
}

bool operator<(const mpfr_class& lhs, const mpfr_class& rhs) {
    return mpfr_cmp(lhs.var_, rhs.var_) < 0;
}

bool operator<(const mpfr_class& lhs, double rhs) {
    return mpfr_cmp_d(lhs.var_, rhs) < 0;
}

bool operator>(const mpfr_class& lhs, const mpfr_class& rhs) {
    return mpfr_cmp(lhs.var_, rhs.var_) > 0;
}

bool operator>(const mpfr_class& lhs, double rhs) {
    return mpfr_cmp_d(lhs.var_, rhs) > 0;
}

mpfr_class abs(const mpfr_class& val) {
    mpfr_class result;
    mpfr_abs(result.var_, val.var_, MPFR_RNDNA);

    return result;
}

std::ostream& operator<<(std::ostream& os, const mpfr_class& val) {
//    char* s;
//    mp_exp_t  e;
//    s = mpfr_get_str (NULL, &e, 0, 0, val.var_, MPFR_RNDNA);
//    mpfr_get_d(val.var_, MPFR_RNDNA);
//    os << s;
//    mpfr_free_str(s);
    os << mpfr_get_d(val.var_, MPFR_RNDNA);
    return os;
}
/*
std::istream& operator>>(std::istream &is, mpfr_class& val) {
    std::string tmp;
    std::getline(is, tmp);
    mpfr_set_str(val.var_, tmp.c_str(), 0, MPFR_RNDNA);
    return is;
}
*/
mpfr_class mpfr_pi() {
    mpfr_class result;
    mpfr_const_pi(result.var_, MPFR_RNDNA);
    return result;
}
