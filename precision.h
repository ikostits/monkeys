//
//  precision.h
//  monkeys
//
//  Created by Irina Kostitsyna on 6/11/14.
//
//

#ifndef PRECISION_H
#define PRECISION_H

#include <mpfr.h>
#include <CGAL/Gmpq.h>
#include <iostream>

#define MPFR_PRECISION 512
#define PRECISION 1e-10
#define SQUARED_PRECISION 1e-20

class mpfr_class {
    friend mpfr_class operator+(const mpfr_class& lhs, const mpfr_class& rhs);
    friend mpfr_class operator-(const mpfr_class& lhs, const mpfr_class& rhs);
    friend mpfr_class operator*(const mpfr_class& lhs, const mpfr_class& rhs);
    friend mpfr_class operator*(double lhs, const mpfr_class& rhs);
    friend mpfr_class operator*(const mpfr_class& lhs, double rhs);
    friend mpfr_class operator/(const mpfr_class& lhs, const mpfr_class& rhs);
    friend mpfr_class operator/(const mpfr_class& lhs, double rhs);
    friend bool operator<(const mpfr_class& lhs, const mpfr_class& rhs);
    friend bool operator<(const mpfr_class& lhs, double rhs);
    friend bool operator>(const mpfr_class& lhs, const mpfr_class& rhs);
    friend bool operator>(const mpfr_class& lhs, double rhs);
    friend mpfr_class logarithm(const mpfr_class& vt);
    friend mpfr_class power(const mpfr_class& vt, long i);
    friend mpfr_class abs(const mpfr_class& val);
    friend std::ostream& operator<<(std::ostream& os, const mpfr_class& val);
//    friend std::istream& operator>>(std::istream& is, mpfr_class& val);
    friend mpfr_class mpfr_pi();
public:
    mpfr_class() {
        mpfr_init2(var_, MPFR_PRECISION);
    };
    explicit mpfr_class(const mpfr_t& rhs) {
        mpfr_init2(var_, mpfr_get_prec(rhs));
        mpfr_set(var_, rhs, MPFR_RNDNA);
    };
    mpfr_class(const mpfr_class& rhs) {
        mpfr_init2(var_, mpfr_get_prec(rhs.var_));
        mpfr_set(var_, rhs.var_, MPFR_RNDNA);
    };
    explicit mpfr_class(double rhs) {
        mpfr_init2(var_, MPFR_PRECISION);
        mpfr_set_d(var_, rhs, MPFR_RNDNA);
    };
    explicit mpfr_class(int rhs) {
        mpfr_init2(var_, MPFR_PRECISION);
        mpfr_set_si(var_, rhs, MPFR_RNDNA);
    };
    template <class NT>
    explicit mpfr_class(const NT& rhs) {
        mpfr_init2(var_, MPFR_PRECISION);
        mpfr_set_q(var_, rhs.exact().mpq(), MPFR_RNDNA);
    };
    virtual ~mpfr_class() {
        mpfr_clear(var_);
    };

    mpfr_class& operator=(mpfr_class arg) {
        mpfr_swap(var_, arg.var_);
        return *this;
    }
    bool is_zero() const {
//        return mpfr_zero_p(var_) != 0;
        return abs(*this) < PRECISION;
    };
    bool is_not_zero() const {
        return mpfr_zero_p(var_) == 0;
    };
    mpfr_class operator-() const {
        mpfr_class result;
        mpfr_neg(result.var_, var_, MPFR_RNDNA);

        return result;
    };
    bool operator==(double rhs) const {
        if (rhs == 0)
            return this->is_zero();
        else
            return mpfr_cmp_d(var_, rhs) == 0;
    };
    bool operator!=(double rhs) const {
        if (rhs == 0)
            return this->is_not_zero();
        else
            return mpfr_cmp_d(var_, rhs) != 0;
    };
    bool operator==(const mpfr_class& rhs) const {
        return (mpfr_equal_p(var_, rhs.var_) != 0);
    };
    mpfr_class& operator+=(const mpfr_class& rhs) {
        mpfr_add(var_, var_, rhs.var_, MPFR_RNDNA);
        return *this;
    };
    mpfr_class& operator-=(const mpfr_class& rhs) {
        mpfr_sub(var_, var_, rhs.var_, MPFR_RNDNA);
        return *this;
    };
    double get_d() {
        return mpfr_get_d(var_, MPFR_RNDNA);
    }
private:
    mpfr_t var_;
};

mpfr_class mpfr_pi();

#endif
