#include "integral.h"
#include <omp.h>
#include <cmath>
#include <random>
#include <vector>

double runge(double I_h, double I_2h, int k)
{
    return std::abs(I_h - I_2h) / (std::pow(2, k) - 1);
}

double rectangle_l(func_x f, double a, double b, size_t n)
{
    double h = (b - a) / n;
    double result = 0;
#pragma omp parallel for schedule(static) reduction(+ \
                                                    : result)
    for (size_t i = 0; i < n - 1; i++)
    {
        //double x = a + h * i;
        result += f(a + h * i);
    }
    return h * result;
}

double rectangle_r(func_x f, double a, double b, size_t n)
{
    double h = (b - a) / n;
    double result = 0;
#pragma omp parallel for schedule(static) reduction(+ \
                                                    : result)
    for (size_t i = 1; i < n; i++)
    {
        double x = a + h * i;
        result += f(x);
    }
    return h * result;
}

double rectangle_m(func_x f, double a, double b, size_t n)
{
    double h = (b - a) / n;
    double result = 0;
#pragma omp parallel for reduction(+ \
                                   : result)
    for (size_t i = 0; i < n - 1; i++)
    {
        double x = a + h * i + h / 2;
        result += f(x);
    }
    return h * result;
}

double trapezoidal(func_x f, double a, double b, size_t n)
{
    double h = (b - a) / n;
    double result = 0;

#pragma omp parallel for schedule(static) reduction(+ \
                                                    : result)
    for (size_t i = 0; i < n; i++)
    {
        double x1 = a + h * i;
        double x2 = a + h * (i + 1);
        result += f(x1) + f(x2);
    }
    return h / 2 * result;
}

double simpson(func_x f, double a, double b, size_t n)
{
    double h = (b - a) / n;
    double result = 0;

#pragma omp parallel for schedule(static) reduction(+ \
                                                    : result)
    for (size_t i = 0; i < n - 1; i += 2)
    {
        double x0 = a + h * i;
        double x1 = a + h * (i + 1);
        double x2 = a + h * (i + 2);
        result += f(x0) + 4 * f(x1) + f(x2);
    }
    return h / 3 * result;
}

double newton_38(func_x f, double a, double b, size_t n)
{
    double h = (b - a) / n;
    double result = 0;

#pragma omp parallel for schedule(static) reduction(+ \
                                                    : result)
    for (size_t i = 0; i < n - 1; i += 3)
    {
        double x0 = a + h * i;
        double x1 = a + h * (i + 1);
        double x2 = a + h * (i + 2);
        double x3 = a + h * (i + 3);
        result += f(x0) + 3 * (f(x1) + f(x2)) + f(x3);
    }
    return 3. / 8 * h * result;
}

double abs_err_rect_rl(double max_1, double a, double b, size_t n)
{
    return max_1 * std::pow(b - a, 2) / (2 * n);
}

double abs_err_rect_m(double max_2, double a, double b, size_t n)
{
    return max_2 * std::pow(b - a, 3) / (24 * n * n);
}

double abs_err_trap(double max_2, double a, double b, size_t n)
{
    return max_2 * std::pow(b - a, 3) / (12 * n * n);
}

double abs_err_simps(double max_4, double a, double b, size_t n)
{
    return max_4 * std::pow(b - a, 5) / (2880 * std::pow(n, 4));
}

double abs_err_newton_38(double max_4, double a, double b, size_t n)
{
    return max_4 * std::pow(b - a, 5) / (80 * std::pow(n, 4));
}

// -------------------------------------------------------------------------
// -------------------------------------------------------------------------
// -------------------------------------------------------------------------

double newton_38_2d(func_xy f, double a, double b, double c, double d, size_t n)
{
    double h_x = (b - a) / n;
    double h_y = (d - c) / n;
    double sum = 0;

#pragma omp parallel for schedule(static) reduction(+ \
                                                    : sum)
    for (size_t i = 0; i < n - 1; i += 3)
        for (size_t j = 0; j < n - 1; j += 3)
        {
            double x0 = a + h_x * i;
            double x1 = a + h_x * (i + 1);
            double x2 = a + h_x * (i + 2);
            double x3 = a + h_x * (i + 3);

            double y0 = c + h_y * j;
            double y1 = c + h_y * (j + 1);
            double y2 = c + h_y * (j + 2);
            double y3 = c + h_y * (j + 3);

            sum += f(x0, y0) + 3 * f(x1, y0) + 3 * f(x2, y0) + f(x3, y0) +
                   3 * (f(x0, y1) + 3 * f(x1, y1) + 3 * f(x2, y1) + f(x3, y1)) +
                   3 * (f(x0, y2) + 3 * f(x1, y2) + 3 * f(x2, y2) + f(x3, y2)) +
                   f(x0, y3) + 3 * f(x1, y3) + 3 * f(x2, y3) + f(x3, y3);
        }
    return 3. / 8 * 3. / 8 * h_x * h_y * sum;
}

double abs_err_newton_38_2d(func_xy f, double a, double b, double c, double d, size_t n)
{
    // TO DO
    return -1;
}

double rectangle_r_2d(func_xy f, double a, double b, double c, double d, size_t n)
{
    double h_x = (b - a) / n;
    double h_y = (d - c) / n;
    double sum = 0;
#pragma omp parallel for schedule(static) reduction(+ \
                                                    : sum)
    for (size_t i = 1; i < n; i++)
        for (size_t j = 1; j < n; j++)
        {
            double x = a + h_x * i;
            double y = c + h_y * j;
            sum += f(x, y);
        }
    return h_x * h_y * sum;
}

double rectangle_l_2d(func_xy f, double a, double b, double c, double d, size_t n)
{
    double h_x = (b - a) / n;
    double h_y = (d - c) / n;
    double sum = 0;
#pragma omp parallel for schedule(static) reduction(+ \
                                                    : sum)
    for (size_t i = 0; i < n - 1; i++)
        for (size_t j = 0; j < n - 1; j++)
        {
            double x = a + h_x * i;
            double y = c + h_y * j;
            sum += f(x, y);
        }
    return h_x * h_y * sum;
}

double abs_err_rect_rl_2d(func_xy f, double a, double b, double c, double d, size_t n)
{
    // TO DO
    return -1;
}

double rectangle_m_2d(func_xy f, double a, double b, double c, double d, size_t n)
{
    double h_x = (b - a) / n;
    double h_y = (d - c) / n;
    double sum = 0;
#pragma omp parallel for schedule(static) reduction(+ \
                                                    : sum)
    for (size_t i = 0; i < n - 1; i++)
        for (size_t j = 0; j < n - 1; j++)
        {
            double x = a + h_x * i + h_x / 2;
            double y = c + h_y * j + h_y / 2;
            sum += f(x, y);
        }
    return h_x * h_y * sum;
}

double abs_err_rect_m_2d(func_xy f_xy, func_xy f_x, func_xy f_y, double a, double b, double c, double d, size_t n)
{
    double res = 0;
    double h_x = (b - a) / n;
    double h_y = (c - d) / n;
    double max_x = 0;
    double max_y = 0;
    double max_xy = 0;

#pragma omp parallel for schedule(static)
    for (size_t i = 0; i < n; i++)
    {
        double x = a + h_x * i;
        max_x = std::max(max_x, std::fabs(f_x(x, (c + d) / 2)));
    }
#pragma omp parallel for schedule(static)
    for (size_t i = 0; i < n; i++)
    {
        double y = c + h_y * i;
        max_y = std::max(max_y, std::fabs(f_y((a + b) / 2, y)));
    }
#pragma omp parallel for schedule(static)
    for (size_t i = 0; i < n - 1; i++)
        for (size_t j = 0; j < n - 1; j++)
        {
            double x = a + h_x * i + h_x / 2;
            double y = c + h_y * j + h_y / 2;
            max_xy = std::max(max_xy, std::fabs(f_xy(x, y)));
        }

    res = (d - c) * (d - c) / 2 * (b - a) * (b - a) / 2 * max_xy - (d - c) * (d - c) / 2 * max_y - (d - c) * (b - a) * (b - a) / 2 * max_x;
    return res;
}

double trapezoidal_2d(func_xy f, double a, double b, double c, double d, size_t n)
{
    double h_x = (b - a) / n;
    double h_y = (d - c) / n;
    double sum = 0;

#pragma omp parallel for schedule(static) reduction(+ \
                                                    : sum)
    for (size_t i = 0; i < n - 1; i++)
        for (size_t j = 0; j < n - 1; j++)
        {
            double x0 = a + h_x * i;
            double x1 = a + h_x * (i + 1);

            double y0 = c + h_y * j;
            double y1 = c + h_y * (j + 1);

            sum += f(x0, y0) + f(x1, y1);
        }
    return h_x * h_y * sum / 2;
}

double abs_err_trap_2d(func_xy f, double a, double b, double c, double d, size_t n)
{
    // TO DO
    return -1;
}

double simpson_2d(func_xy f, double a, double b, double c, double d, size_t n)
{
    double h_x = (b - a) / n;
    double h_y = (d - c) / n;
    double sum = 0;

#pragma omp parallel for schedule(static) reduction(+ \
                                                    : sum)
    for (size_t i = 0; i < n - 1; i += 2)
        for (size_t j = 0; j < n - 1; j += 2)
        {
            double x0 = a + h_x * i;
            double x1 = a + h_x * (i + 1);
            double x2 = a + h_x * (i + 2);

            double y0 = c + h_y * j;
            double y1 = c + h_y * (j + 1);
            double y2 = c + h_y * (j + 2);

            sum += f(x0, y0) + 4 * f(x1, y0) + f(x2, y0) +
                   4 * (f(x0, y1) + 4 * f(x1, y1) + f(x2, y1)) +
                   f(x0, y2) + 4 * f(x1, y2) + f(x2, y2);
        }
    return h_x * h_y * sum / 9;
}

double abs_err_simps_2d(func_xy f, double a, double b, double c, double d, size_t n)
{
    // TO DO
    return -1;
}

double monte_carlo_2d(func_xy f, double a, double b, double c, double d, size_t n)
{
    // TO DO
    return -1;
}

double abs_err_monte_carlo_2d(func_xy f, double a, double b, double c, double d, size_t n)
{
    // TO DO
    return -1;
}
