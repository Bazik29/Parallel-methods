#pragma once
#include <functional>

// using func_x = double(double);
// using func_xy = double(double, double);

using func_x = std::function<double(double)>;
using func_xy = std::function<double(double, double)>;

double runge(double I_h, double I_2h, int k);

// Методы одномерного интегрирования
double rectangle_l(func_x f, double a, double b, size_t n);

double rectangle_r(func_x f, double a, double b, size_t n);

double rectangle_m(func_x f, double a, double b, size_t n);

double trapezoidal(func_x f, double a, double b, size_t n);

double simpson(func_x f, double a, double b, size_t n);

double newton_38(func_x f, double a, double b, size_t n);

// Погрешности

double abs_err_rect_rl(double max_1, double a, double b, size_t n);

double abs_err_rect_m(double max_2, double a, double b, size_t n);

double abs_err_trap(double max_2, double a, double b, size_t n);

double abs_err_simps(double max_4, double a, double b, size_t n);

double abs_err_newton_38(double max_4, double a, double b, size_t n);

// Методы двумерного интегрирования
double newton_38_2d(func_xy f, double a, double b, double c, double d, size_t n);

double abs_err_newton_38_2d(func_xy f, double a, double b, double c, double d, size_t n);

double rectangle_r_2d(func_xy f, double a, double b, double c, double d, size_t n);

double rectangle_l_2d(func_xy f, double a, double b, double c, double d, size_t n);

double abs_err_rect_rl_2d(func_xy f, double a, double b, double c, double d, size_t n);

double rectangle_m_2d(func_xy f, double a, double b, double c, double d, size_t n);

double abs_err_rect_m_2d(func_xy f_2xy, func_xy f_1x, func_xy f_1y, double a, double b, double c, double d, size_t n);

double trapezoidal_2d(func_xy f, double a, double b, double c, double d, size_t n);

double abs_err_trap_2d(func_xy f, double a, double b, double c, double d, size_t n);

double simpson_2d(func_xy f, double a, double b, double c, double d, size_t n);

double abs_err_simps_2d(func_xy f, double a, double b, double c, double d, size_t n);

double monte_carlo_2d(func_xy f, double a, double b, double c, double d, size_t n);

double abs_err_monte_carlo_2d(func_xy f, double a, double b, double c, double d, size_t n);
