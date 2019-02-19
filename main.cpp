#include <omp.h>

#include <iostream>
#include <iomanip>
#include <functional>
#include <cmath>

#include <chrono>

const double PI = 3.141592653589793238463;

// 1.3737
double foo_a = 1;
double foo_b = 2;
double foo(double x)
{
    //x^2*arcctg(x/2)/x
    //x^2*(pi/2 - arctg(x/2))/x
    return x * x * (PI / 2 - std::atan(x / 2)) / x;
}

double foo_1(double x)
{
    //arcctg(x/2) - 2x/(x^2+4)
    //pi/2 - arctg(x/2) - 2x/(x^2+4)
    return PI / 2 - std::atan(x / 2) - 2 * x / (x * x + 4);
}

double foo_2(double x)
{
    //-16/(x^2 + 4)^2
    return -16 / ((x * x + 4) * (x * x + 4));
}

double foo_4(double x)
{
    //-16*((24*x^2/(x^2+4)^4) - 4/(x^2+4)^2)
    return -16 * (24 * x * x / std::pow((x * x + 4), 4) - 4 / ((x * x + 4) * (x * x + 4)));
}

double rectangle_l(std::function<double(double)> f, double a, double h, size_t n)
{
    double result = 0;
#pragma omp parallel for reduction(+ : result)
    for (size_t i = 0; i < n - 1; i++)
    {
        double x = a + h * i;
        result += f(x);
    }
    return h * result;
}

double rectangle_r(std::function<double(double)> f, double a, double h, size_t n)
{
    double result = 0;
#pragma omp parallel for reduction(+ : result)
    for (size_t i = 1; i < n; i++)
    {
        double x = a + h * i;
        result += f(x);
    }
    return h * result;
}

double rectangle_m(std::function<double(double)> f, double a, double h, size_t n)
{
    double result = 0;
#pragma omp parallel for reduction(+ : result)
    for (size_t i = 0; i < n - 1; i++)
    {
        double x = a + h * i + h / 2;
        result += f(x);
    }
    return h * result;
}

double abs_err_rect_m(std::function<double(double)> f_2, double a, double b, double h, size_t n)
{
    double maximum = 0;

#pragma omp parallel for
    for (size_t i = 0; i < n; i++)
    {
        double x = a + h * i;
        maximum = std::max(maximum, std::fabs(f_2(x)));
    }
    return maximum * std::pow(b - a, 3) / (24 * n * n);
}

double abs_err_rect_rl(std::function<double(double)> f_1, double a, double b, double h, size_t n)
{
    double maximum = 0;

#pragma omp parallel for
    for (size_t i = 0; i < n; i++)
    {
        double x = a + h * i;
        maximum = std::max(maximum, std::fabs(f_1(x)));
    }
    return maximum * std::pow(b - a, 2) / (2 * n);
}

double trapezoidal(std::function<double(double)> f, double a, double h, size_t n)
{
    double result = 0;

#pragma omp parallel for reduction(+ : result)
    for (size_t i = 0; i < n; i++)
    {
        double x1 = a + h * i;
        double x2 = a + h * (i + 1);
        result += f(x1) + f(x2);
    }
    return h / 2 * result;
}

double abs_err_trap(std::function<double(double)> f_2, double a, double b, double h, size_t n)
{
    double maximum = 0;

#pragma omp parallel for
    for (size_t i = 0; i < n; i++)
    {
        double x = a + h * i;
        maximum = std::max(maximum, std::fabs(f_2(x)));
    }
    return maximum * std::pow(b - a, 3) / (12 * n * n);
}

double simpson(std::function<double(double)> f, double a, double h, size_t n)
{
    double result = 0;

#pragma omp parallel for reduction(+ : result)
    for (size_t i = 1; i < n - 1; i += 2)
    {
        double x1 = a + h * (i - 1);
        double x2 = a + h * i;
        double x3 = a + h * (i + 1);
        result += f(x1) + 4 * f(x2) + f(x3);
    }
    return h / 3 * result;
}

double abs_err_simps(std::function<double(double)> f_4, double a, double b, double h, size_t n)
{
    double maximum = 0;

#pragma omp parallel for
    for (size_t i = 0; i < n; i++)
    {
        double x = a + h * i;
        maximum = std::max(maximum, std::fabs(f_4(x)));
    }
    return maximum * std::pow(b - a, 5) / (2880 * std::pow(n, 4));
}

double newton_38(std::function<double(double)> f, double a, double h, size_t n)
{
    // double result = 0;

// #pragma omp parallel for reduction(+ : result)
//     for (size_t i = 1; i < n - 1; i += 2)
//     {
//         double x1 = a + h * (i - 1);
//         double x2 = a + h * i;
//         double x3 = a + h * (i + 1);
//         result += f(x1) + 4 * f(x2) + f(x3);
//     }
    return 0;
}

double abs_err_newton_38(std::function<double(double)> f_4, double a, double b, double h, size_t n)
{
    // double maximum = 0;

// #pragma omp parallel for
//     for (size_t i = 0; i < n; i++)
//     {
//         double x = a + h * i;
//         maximum = std::max(maximum, std::fabs(f_4(x)));
//     }
    return 0;
}

double runge(double I_h, double I_2h, int k)
{
    return std::abs(I_h - I_2h) / (std::pow(2, k) - 1);
}

using std::cout;
using std::endl;
using std::chrono::steady_clock;
using std::chrono::duration;
using std::chrono::duration_cast;

int main(int argc, char const *argv[])
{
    cout.setf(cout.fixed);
    cout.precision(15);

    size_t n = 100;
    if (argc > 1)
        n = std::atoi(argv[1]);

    double h = (foo_b - foo_a) / n;

    // параметры для оценки погрешности по правилу Рунге
    size_t n_2 = n / 2;
    double h_2 = (foo_b - foo_a) / n_2;

    double result, result_test, rung, abs_err;

    int num_procs = omp_get_num_procs();

    cout << "Число процессоров: " << num_procs << endl;
    cout << "------Методы одномерного численного интегрирования------\n";
    cout << "Число разбиений:           " << n << endl;
    cout << "Размер шага:               " << h << endl;
    cout << "Для правила Рунге\n";
    cout << "Число разбиений:           " << n_2 << endl;
    cout << "Размер шага:               " << h_2 << endl;

    // Формула левых прямоугольников
    auto begTime = steady_clock::now();
    result = rectangle_l(foo, foo_a, h, n);
    auto rect_lTime = duration_cast<duration<double>>(steady_clock::now() - begTime);

    result_test = rectangle_l(foo, foo_a, h_2, n_2);
    rung = runge(result_test, result, 1);

    cout << "\nФормула левых прямоугольников\n";
    cout << "Время:                     " << rect_lTime.count() << endl;
    cout << "Результат:                 " << result << endl;
    cout << "Ошибка по правилу Рунге:   " << rung << endl;

    // Формула правых прямоугольников
    begTime = steady_clock::now();
    result = rectangle_r(foo, foo_a, h, n);
    auto rect_rTime = duration_cast<duration<double>>(steady_clock::now() - begTime);

    result_test = rectangle_r(foo, foo_a, h_2, n_2);
    rung = runge(result_test, result, 1);

    cout << "\nФормула правых прямоугольников\n";
    cout << "Время:                     " << rect_rTime.count() << endl;
    cout << "Результат:                 " << result << endl;
    cout << "Ошибка по правилу Рунге:   " << rung << endl;

    abs_err = abs_err_rect_rl(foo_1, foo_a, foo_b, h, n);
    cout << "Абсолютная погрешность формулы левых и правых прямоугольников: " << abs_err << endl;

    // Формула средних прямоугольников
    begTime = steady_clock::now();
    result = rectangle_m(foo, foo_a, h, n);
    auto rect_mTime = duration_cast<duration<double>>(steady_clock::now() - begTime);

    result_test = rectangle_m(foo, foo_a, h_2, n_2);
    rung = runge(result_test, result, 2);
    abs_err = abs_err_rect_m(foo_2, foo_a, foo_b, h, n);

    cout << "\nФормула правых прямоугольников\n";
    cout << "Время:                     " << rect_mTime.count() << endl;
    cout << "Результат:                 " << result << endl;
    cout << "Ошибка по правилу Рунге:   " << rung << endl;
    cout << "Абсолютная погрешность:    " << abs_err << endl;

    // Формула трапеций
    begTime = steady_clock::now();
    result = trapezoidal(foo, foo_a, h, n);
    auto trapez_Time = duration_cast<duration<double>>(steady_clock::now() - begTime);

    result_test = trapezoidal(foo, foo_a, h_2, n_2);
    rung = runge(result_test, result, 2);
    abs_err = abs_err_trap(foo_2, foo_a, foo_b, h, n);

    cout << "\nФормула трапеций\n";
    cout << "Время:                     " << trapez_Time.count() << endl;
    cout << "Результат:                 " << result << endl;
    cout << "Ошибка по правилу Рунге:   " << rung << endl;
    cout << "Абсолютная погрешность:    " << abs_err << endl;

    // Формула Симпсона
    begTime = steady_clock::now();
    result = simpson(foo, foo_a, h, n);
    auto simps_Time = duration_cast<duration<double>>(steady_clock::now() - begTime);

    result_test = simpson(foo, foo_a, h_2, n_2);
    rung = runge(result_test, result, 4);
    abs_err = abs_err_simps(foo_4, foo_a, foo_b, h, n);

    cout << "\nФормула Симпсона\n";
    cout << "Время:                     " << simps_Time.count() << endl;
    cout << "Результат:                 " << result << endl;
    cout << "Ошибка по правилу Рунге:   " << rung << endl;
    cout << "Абсолютная погрешность:    " << abs_err << endl;

    // Формула Ньютона (правило трех восьмых)
    begTime = steady_clock::now();
    result = newton_38(foo, foo_a, h, n);
    auto newton_Time = duration_cast<duration<double>>(steady_clock::now() - begTime);

    result_test = newton_38(foo, foo_a, h_2, n_2);
    rung = runge(result_test, result, 4);
    abs_err = abs_err_newton_38(foo_2, foo_a, foo_b, h, n);

    cout << "\nФормула Ньютона (правило трех восьмых)\n";
    cout << "Время:                     " << newton_Time.count() << endl;
    cout << "Результат:                 " << result << endl;
    cout << "Ошибка по правилу Рунге:   " << rung << endl;
    cout << "Абсолютная погрешность:    " << abs_err << endl;

    return 0;
}
