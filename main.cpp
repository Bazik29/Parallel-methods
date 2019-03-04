#include <iostream>
#include <iomanip>
#include <cmath>
#include <chrono>
#include <fstream>

#include <omp.h>

#include "integral.h"
#include "utils.h"

const double PI = 3.141592653589793238463;

// 1.373720859104567
double foo_a = 1;
double foo_b = 2;
double foo(double x)
{
    //x^2*arcctg(x/2)/x
    //x^2*(pi/2 - arctg(x/2))/x
    //x*(pi/2 - arctg(x/2))
    return x * (PI / 2 - std::atan(x / 2));
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

using std::cout;
using std::endl;
using std::chrono::duration;
using std::chrono::duration_cast;
using std::chrono::steady_clock;

int main(int argc, char const *argv[])
{
    cout.setf(cout.fixed);
    cout.precision(15);
    /*
    0 - Формула левых прямоугольников
    1 - Формула правых прямоугольников
    2 - Формула средних прямоугольников
    3 - Формула трапеций
    4 - Формула Симпсона
    5 - "Правило трех восьмых"
    */

    // число разбиений (для каждого метода)
    size_t aN[6] = {100, 100, 100, 100, 100, 99};

    if (argc > 1)
    {
        int n = std::atoi(argv[1]);

        for (int i = 0; i < 6; i++)
            aN[i] = n;
        // число разбиений кратно 2 для формулы Симпсона
        if (aN[4] % 2 != 0)
            aN[4] -= 1;
        // число разбиений кратно 3 для "правила трех восьмых"
        while (aN[5] % 3 != 0)
            aN[5] -= 1;
    }

    // число разбиений (для каждого метода) для оценки погрешности по Рунге
    size_t aN_runge[6];
    for (int i = 0; i < 6; i++)
        aN_runge[i] = aN[i] / 2;
    // число разбиений кратно 2 для формулы Симпсона
    if (aN_runge[4] % 2 != 0)
        aN_runge[4] -= 1;
    // число разбиений кратно 3 для "правила трех восьмых"
    while (aN_runge[5] % 3 != 0)
        aN_runge[5] -= 1;

    std::ofstream output;
    output.open("out.csv", std::ios_base::trunc);
    output.setf(cout.fixed);
    output.precision(15);

    double result,
        result_test, rung, abs_err, calc_time;

    int num_procs = omp_get_num_procs();

    cout << "Число процессоров: " << num_procs << endl;
    cout << "------Методы одномерного численного интегрирования------\n";
    cout << "Верный результат:             1.373720859104567\n";
    cout << "--------------------------------------------------------\n";

    // Формула левых прямоугольников
    auto begTime = steady_clock::now();
    result = rectangle_l(foo, foo_a, foo_b, aN[0]);
    calc_time = duration_cast<duration<double>>(steady_clock::now() - begTime).count();

    abs_err = abs_err_rect_rl(foo_1, foo_a, foo_b, aN[0]);

    result_test = rectangle_l(foo, foo_a, foo_b, aN_runge[0]);
    rung = runge(result_test, result, 1);

    print_log(0, aN[0], result, abs_err, aN_runge[0], rung, calc_time);
    file_CSV(output, num_procs, 0, aN[0], calc_time, result, abs_err, rung);

    // Формула правых прямоугольников
    begTime = steady_clock::now();
    result = rectangle_r(foo, foo_a, foo_b, aN[1]);
    calc_time = duration_cast<duration<double>>(steady_clock::now() - begTime).count();

    result_test = rectangle_r(foo, foo_a, foo_b, aN_runge[1]);
    rung = runge(result_test, result, 1);

    print_log(1, aN[1], result, abs_err, aN_runge[1], rung, calc_time);
    file_CSV(output, num_procs, 1, aN[1], calc_time, result, abs_err, rung);

    // Формула средних прямоугольников
    begTime = steady_clock::now();
    result = rectangle_m(foo, foo_a, foo_b, aN[2]);
    calc_time = duration_cast<duration<double>>(steady_clock::now() - begTime).count();

    abs_err = abs_err_rect_m(foo_2, foo_a, foo_b, aN[2]);

    result_test = rectangle_m(foo, foo_a, foo_b, aN_runge[2]);
    rung = runge(result_test, result, 2);

    print_log(2, aN[2], result, abs_err, aN_runge[2], rung, calc_time);
    file_CSV(output, num_procs, 2, aN[2], calc_time, result, abs_err, rung);

    // Формула трапеций
    begTime = steady_clock::now();
    result = trapezoidal(foo, foo_a, foo_b, aN[3]);
    calc_time = duration_cast<duration<double>>(steady_clock::now() - begTime).count();

    abs_err = abs_err_trap(foo_2, foo_a, foo_b, aN[3]);

    result_test = trapezoidal(foo, foo_a, foo_b, aN_runge[3]);
    rung = runge(result_test, result, 2);

    print_log(3, aN[3], result, abs_err, aN_runge[3], rung, calc_time);
    file_CSV(output, num_procs, 3, aN[3], calc_time, result, abs_err, rung);

    // Формула Симпсона
    begTime = steady_clock::now();
    result = simpson(foo, foo_a, foo_b, aN[4]);
    calc_time = duration_cast<duration<double>>(steady_clock::now() - begTime).count();

    abs_err = abs_err_simps(foo_4, foo_a, foo_b, aN[4]);

    result_test = simpson(foo, foo_a, foo_b, aN_runge[4]);
    rung = runge(result_test, result, 4);

    print_log(4, aN[4], result, abs_err, aN_runge[4], rung, calc_time);
    file_CSV(output, num_procs, 4, aN[4], calc_time, result, abs_err, rung);

    // Формула Ньютона (правило трех восьмых)
    begTime = steady_clock::now();
    result = newton_38(foo, foo_a, foo_b, aN[5]);
    calc_time = duration_cast<duration<double>>(steady_clock::now() - begTime).count();

    abs_err = abs_err_newton_38(foo_2, foo_a, foo_b, aN[5]);

    result_test = newton_38(foo, foo_a, foo_b, aN_runge[5]);
    rung = runge(result_test, result, 4);

    print_log(5, aN[5], result, abs_err, aN_runge[5], rung, calc_time);
    file_CSV(output, num_procs, 5, aN[5], calc_time, result, abs_err, rung);

    output.close();

    return 0;
}
