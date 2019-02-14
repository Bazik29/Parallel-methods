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

double runge(double I_h, double I_2h, int k)
{
    return std::abs(I_h - I_2h) / (std::pow(2, k) - 1);
}

int main(int argc, char const *argv[])
{
    std::cout.setf(std::cout.fixed);
    std::cout.precision(15);
    size_t n = 100;
    if (argc > 1)
        n = std::atoi(argv[1]);

    double h = (foo_b - foo_a) / n;

    // параметры для оценки погрешности по правилу Рунге
    size_t n_2 = n / 2;
    double h_2 = (foo_b - foo_a) / n_2;

    double result, result_test, rung, abs_err;

    // auto startTime = std::chrono::steady_clock::now();
    // auto trapezoidalTime = std::chrono::steady_clock::now();

    // Метод левых прямоугольников
    result = rectangle_l(foo, foo_a, h, n);
    result_test = rectangle_l(foo, foo_a, h_2, n_2);
    rung = runge(result_test, result, 1);

    std::cout << "Метод левых прямоугольников\n";
    // std::cout << "Время:     " << runtimeDuration.count() << std::endl;
    std::cout << "Результат:                " << result << std::endl;
    std::cout << "Ошибка по правилу Рунге:  " << rung << std::endl;

    // Метод правых прямоугольников
    result = rectangle_r(foo, foo_a, h, n);
    result_test = rectangle_r(foo, foo_a, h_2, n_2);
    rung = runge(result_test, result, 1);

    std::cout << "\nМетод правых прямоугольников.\n";
    // std::cout << "Время:     " << runtimeDuration.count() << std::endl;
    std::cout << "Результат:                " << result << std::endl;
    std::cout << "Ошибка по правилу Рунге:  " << rung << std::endl;

    abs_err = abs_err_rect_rl(foo_1, foo_a, foo_b, h, n);
    std::cout << "Абсолютная погрешность методов левых и правых прямоугольников:\n";
    std::cout << abs_err << std::endl;

    // Метод средних прямоугольников
    result = rectangle_m(foo, foo_a, h, n);
    result_test = rectangle_m(foo, foo_a, h_2, n_2);
    rung = runge(result_test, result, 1);
    abs_err = abs_err_rect_m(foo_2, foo_a, foo_b, h, n);

    std::cout << "\nМетод правых прямоугольников.\n";
    // std::cout << "Время:     " << runtimeDuration.count() << std::endl;
    std::cout << "Результат:                " << result << std::endl;
    std::cout << "Ошибка по правилу Рунге:  " << rung << std::endl;
    std::cout << "Абсолютная погрешность:   " << abs_err << std::endl;

    // Метод трапеций
    result = trapezoidal(foo, foo_a, h, n);
    result_test = trapezoidal(foo, foo_a, h_2, n_2);
    rung = runge(result_test, result, 2);
    abs_err = abs_err_trap(foo_2, foo_a, foo_b, h, n);

    std::cout << "\nМетод трапеций.\n";
    // std::cout << "Время:     " << runtimeDuration.count() << std::endl;
    std::cout << "Результат:                " << result << std::endl;
    std::cout << "Ошибка по правилу Рунге:  " << rung << std::endl;
    std::cout << "Абсолютная погрешность:   " << abs_err << std::endl;

    // auto finishTime = std::chrono::steady_clock::now();
    // auto runtimeDuration = std::chrono::duration_cast<std::chrono::duration<double>>(finishTime - startTime);


    // std::cout << "Time:     " << runtimeDuration.count() << std::endl;
    // std::cout << "Result:   " << result << std::endl;
    // std::cout << "Eps:      " << eps << std::endl;
    // std::cout << "Step:     " << h << std::endl;
    // std::cout << "N steps:  " << n << std::endl;
    // std::cout << "Num of attempts: " << count << std::endl;

    return 0;
}
