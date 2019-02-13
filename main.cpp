#include <omp.h>

#include <iostream>
#include <iomanip>
#include <functional>
#include <cmath>

#include <chrono>

// 0.21548
double foo_a = 0.4;
double foo_b = 1;
double foo(double x)
{
    return std::exp(-x * x + 0.38) / (2 + std::sin(1 / (1.5 + x * x)));
}

double trapezoidal(std::function<double(double)> f, double a, double h, size_t n)
{
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

double runge(double I_h, double I_2h, int k)
{
    return std::abs(I_h - I_2h) / (std::pow(2, k) - 1);
}

int main(int argc, char const *argv[])
{
    for (int i = 1; i < 10; i++)
    {
    auto startTime = std::chrono::steady_clock::now();

    double eps = std::pow(10,-i);
    // if (argc > 1)
    //     eps = std::atof(argv[1]);

    double h = std::pow(eps, 1. / 4.);
    size_t n = std::ceil((foo_b - foo_a) / h);
    //if (!(n & 1)) n+=1;

    double result, result_test, rung;
    int count = 1;
    // std::cout << " n: ";
    for (;;)
    {
        // std::cout << n << " ";
        result = trapezoidal(foo, foo_a, h, n);

        double h_2 = h * 2;
        size_t n_2 = std::ceil((foo_b - foo_a) / h_2);

        result_test = trapezoidal(foo, foo_a, h_2, n_2);

        rung = runge(result_test, result, 1);

        if (rung <= eps)
            break;

        h = h / 2.;
        n = std::ceil((foo_b - foo_a) / h);

        count += 1;
    }
    // std::cout << "\n";

    auto finishTime = std::chrono::steady_clock::now();

    auto runtimeDuration = std::chrono::duration_cast<std::chrono::duration<double>>(finishTime - startTime);

    std::cout.setf(std::cout.fixed);
    std::cout.precision(15);
    //std::cout << "Time:     " << runtimeDuration.count() << std::endl;
    //std::cout << "Result:   " << result << std::endl;
    //std::cout << "Eps:      " << eps << std::endl;
    //std::cout << "Step:     " << h << std::endl;
    //std::cout << "N steps:  " << n << std::endl;
    //std::cout << "Num of attempts: " << count << std::endl;

    std::cout << result << " eps: " << eps << " h: " << h << " n: " << n << "\trunge: " << rung << " " << count <<std::endl;
}
    return 0;
}
