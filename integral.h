#pragma once
#include <functional>

typedef double(*func_x)(double);
typedef double(*func_xy)(double, double);

double runge(double I_h, double I_2h, int k);

// Методы одномерного интегрирования
double rectangle_l(func_x f, double a, double b, size_t n);

double rectangle_r(func_x f, double a, double b, size_t n);

double rectangle_m(func_x f, double a, double b, size_t n);

double abs_err_rect_rl(func_x f_1, double a, double b, size_t n);

double abs_err_rect_m(func_x f_2, double a, double b, size_t n);

double trapezoidal(func_x f, double a, double b, size_t n);

double abs_err_trap(func_x f_2, double a, double b, size_t n);

double simpson(func_x f, double a, double b, size_t n);

double abs_err_simps(func_x f_4, double a, double b, size_t n);

double newton_38(func_x f, double a, double b, size_t n);

double abs_err_newton_38(func_x f_4, double a, double b, size_t n);

double monte_carlo_1d(func_x f, double a, double b, size_t n);

double abs_err_monte_carlo_1d(func_x f, double a, double b, size_t n);

// Методы двумерного интегрирования
double newton_38_2d(func_xy f, double a, double b, double c, double d, size_t n);

double abs_err_newton_38_2d(func_xy f, double a, double b, double c, double d, size_t n);

double rectangle_r_2d(func_xy f, double a, double b, double c, double d, size_t n);

double rectangle_l_2d(func_xy f, double a, double b, double c, double d, size_t n);

double abs_err_rect_rl_2d(func_xy f, double a, double b, double c, double d, size_t n);

double rectangle_m_2d(func_xy f, double a, double b, double c, double d, size_t n);

double abs_err_rect_m_2d(func_xy f, double a, double b, double c, double d, size_t n);

double trapezoidal_2d(func_xy f, double a, double b, double c, double d, size_t n);

double abs_err_trap_2d(func_xy f, double a, double b, double c, double d, size_t n);

double simpson_2d(func_xy f, double a, double b, double c, double d, size_t n);

double abs_err_simps_2d(func_xy f, double a, double b, double c, double d, size_t n);

double monte_carlo_2d(func_xy f, double a, double b, double c, double d, size_t n);

double abs_err_monte_carlo_2d(func_xy f, double a, double b, double c, double d, size_t n);

