#include "integral.h"
#include <omp.h>
#include <cmath>
#include <random>
#include <vector>

double runge(double I_h, double I_2h, int k)
{
    return std::abs(I_h - I_2h) / (std::pow(2, k) - 1);
}

double rectangle_l(std::function<double(double)> f, double a, double b, size_t n)
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

double rectangle_r(std::function<double(double)> f, double a, double b, size_t n)
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

double rectangle_m(std::function<double(double)> f, double a, double b, size_t n)
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

double abs_err_rect_rl(std::function<double(double)> f_1, double a, double b, size_t n)
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

double abs_err_rect_m(std::function<double(double)> f_2, double a, double b, size_t n)
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

double trapezoidal(std::function<double(double)> f, double a, double b, size_t n)
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

double abs_err_trap(std::function<double(double)> f_2, double a, double b, size_t n)
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

double simpson(std::function<double(double)> f, double a, double b, size_t n)
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

double abs_err_simps(std::function<double(double)> f_4, double a, double b, size_t n)
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

double newton_38(std::function<double(double)> f, double a, double b, size_t n)
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

double abs_err_newton_38(std::function<double(double)> f_4, double a, double b, size_t n)
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

double monte_carlo_1d(std::function<double(double)> f, double a, double b, size_t n)
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
