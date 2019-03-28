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
#pragma omp parallel for reduction(+ \
                                   : result)
    for (size_t i = 0; i < n - 1; i++)
    {
        double x = a + h * i;
        result += f(x);
    }
    return h * result;
}

double rectangle_r(func_x f, double a, double b, size_t n)
{
    double h = (b - a) / n;
    double result = 0;
#pragma omp parallel for reduction(+ \
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

double abs_err_rect_rl(func_x f_1, double a, double b, size_t n)
{
    double h = (b - a) / n;
    double maximum = 0;

#pragma omp parallel for
    for (size_t i = 0; i < n; i++)
    {
        double x = a + h * i;
        maximum = std::max(maximum, std::fabs(f_1(x)));
    }
    return maximum * std::pow(b - a, 2) / (2 * n);
}

double abs_err_rect_m(func_x f_2, double a, double b, size_t n)
{
    double h = (b - a) / n;
    double maximum = 0;

#pragma omp parallel for
    for (size_t i = 0; i < n; i++)
    {
        double x = a + h * i;
        maximum = std::max(maximum, std::fabs(f_2(x)));
    }
    return maximum * std::pow(b - a, 3) / (24 * n * n);
}

double trapezoidal(func_x f, double a, double b, size_t n)
{
    double h = (b - a) / n;
    double result = 0;

#pragma omp parallel for reduction(+ \
                                   : result)
    for (size_t i = 0; i < n; i++)
    {
        double x1 = a + h * i;
        double x2 = a + h * (i + 1);
        result += f(x1) + f(x2);
    }
    return h / 2 * result;
}

double abs_err_trap(func_x f_2, double a, double b, size_t n)
{
    double h = (b - a) / n;
    double maximum = 0;

#pragma omp parallel for
    for (size_t i = 0; i < n; i++)
    {
        double x = a + h * i;
        maximum = std::max(maximum, std::fabs(f_2(x)));
    }
    return maximum * std::pow(b - a, 3) / (12 * n * n);
}

double simpson(func_x f, double a, double b, size_t n)
{
    double h = (b - a) / n;
    double result = 0;

#pragma omp parallel for reduction(+ \
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

double abs_err_simps(func_x f_4, double a, double b, size_t n)
{
    double h = (b - a) / n;
    double maximum = 0;

#pragma omp parallel for
    for (size_t i = 0; i < n; i++)
    {
        double x = a + h * i;
        maximum = std::max(maximum, std::fabs(f_4(x)));
    }
    return maximum * std::pow(b - a, 5) / (2880 * std::pow(n, 4));
}

double newton_38(func_x f, double a, double b, size_t n)
{
    double h = (b - a) / n;
    double result = 0;

#pragma omp parallel for reduction(+ \
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

double abs_err_newton_38(func_x f_4, double a, double b, size_t n)
{
    double h = (b - a) / n;
    double maximum = 0;

#pragma omp parallel for
    for (size_t i = 0; i < n; i++)
    {
        double x = a + h * i;
        maximum = std::max(maximum, std::fabs(f_4(x)));
    }
    return maximum * std::pow(b - a, 5) / (80 * std::pow(n, 4));
}

double monte_carlo_1d(func_x f, double a, double b, size_t n)
{
    std::random_device randD;
    std::mt19937 randMT(randD());
    std::uniform_real_distribution<> rand_x(a, b);

    double result = 0;
#pragma omp parallel for reduction(+ \
                                   : result)
    for (size_t i = 0; i < n; i++)
    {
        result += f(rand_x(randMT));
    }
    return (b - a) / n * result;
}

double abs_err_monte_carlo_1d(func_x f, double a, double b, size_t n)
{
    // TO DO
    return -1;
}

double newton_38_2d(func_xy f, double a, double b, double c, double d, size_t n)
{
    double h_x = (b - a) / n;
    double h_y = (d - c) / n;
    double sum = 0;

#pragma omp parallel for reduction(+ \
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
#pragma omp parallel for reduction(+ \
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
#pragma omp parallel for reduction(+ \
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
#pragma omp parallel for reduction(+ \
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

double abs_err_rect_m_2d(func_xy f, double a, double b, double c, double d, size_t n)
{
    // TO DO
    return -1;
}

double trapezoidal_2d(func_xy f, double a, double b, double c, double d, size_t n)
{
    double h_x = (b - a) / n;
    double h_y = (d - c) / n;
    double sum = 0;

#pragma omp parallel for reduction(+ \
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

#pragma omp parallel for reduction(+ \
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
