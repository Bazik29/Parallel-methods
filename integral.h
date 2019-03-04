#pragma once
#include <functional>

double runge(double I_h, double I_2h, int k);

// Методы одномерного интегрирования
double rectangle_l(std::function<double(double)> f, double a, double b, size_t n);

double rectangle_r(std::function<double(double)> f, double a, double b, size_t n);

double rectangle_m(std::function<double(double)> f, double a, double b, size_t n);

double abs_err_rect_rl(std::function<double(double)> f_1, double a, double b, size_t n);

double abs_err_rect_m(std::function<double(double)> f_2, double a, double b, size_t n);

double trapezoidal(std::function<double(double)> f, double a, double b, size_t n);

double abs_err_trap(std::function<double(double)> f_2, double a, double b, size_t n);

double simpson(std::function<double(double)> f, double a, double b, size_t n);

double abs_err_simps(std::function<double(double)> f_4, double a, double b, size_t n);

double newton_38(std::function<double(double)> f, double a, double b, size_t n);

double abs_err_newton_38(std::function<double(double)> f_4, double a, double b, size_t n);
