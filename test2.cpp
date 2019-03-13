#include <iostream>
#include <iomanip>
#include <chrono>
#include <unistd.h>
#include <omp.h>

#include "integral.h"
#include "utils.h"

#include "function.h"

using std::cout;
using std::endl;
using std::chrono::duration;
using std::chrono::duration_cast;
using std::chrono::steady_clock;

int main(int argc, char *argv[])
{
    cout.setf(cout.fixed);
    cout.precision(15);
    /*
    0 - Формула левых прямоугольников
    1 - Формула правых прямоугольников
    2 - Формула средних прямоугольников
    3 - Формула трапеций
    4 - Формула Симпсона
    5 - Формула трех восьмых
    6 - Метод Монте-Карло
    */

    const int N_MET = 7;
    // число разбиений (для каждого метода)
    size_t aN[N_MET] = {100, 100, 100, 100, 100, 99};
    int m = 0;

    bool withError = false;

    int opt;
    while ((opt = getopt(argc, argv, "n:m:e?")) != -1)
    {
        switch (opt)
        {
        case 'n':
        {
            int n_steps = std::atoi(optarg);
            for (int i = 0; i < N_MET; i++)
                aN[i] = n_steps;
            // число разбиений кратно 2 для формулы Симпсона
            if (aN[4] % 2 != 0)
                aN[4] -= 1;
            // число разбиений кратно 3 для "правила трех восьмых"
            while (aN[5] % 3 != 0)
                aN[5] -= 1;
            break;
        }
        case 'm':
        {
            m = std::atoi(optarg);
            break;
        }
        case 'e':
        {
            withError = true;
            break;
        }
        }
    }

    // число разбиений (для каждого метода) для оценки погрешности по Рунге
    size_t aN_runge[N_MET];
    for (int i = 0; i < N_MET; i++)
        aN_runge[i] = aN[i] / 2;
    // число разбиений кратно 2 для формулы Симпсона
    if (aN_runge[4] % 2 != 0)
        aN_runge[4] -= 1;
    // число разбиений кратно 3 для "правила трех восьмых"
    while (aN_runge[5] % 3 != 0)
        aN_runge[5] -= 1;

    double result = 0,
           result_test = 0, rung = -1, abs_err = -1, calc_time = 0;

    int num_procs = omp_get_max_threads(); //omp_get_num_procs();
    auto begTime = steady_clock::now(); // power of "auto"=)))

    switch (m)
    {
    case 0:
    {
        // Формула левых прямоугольников
        begTime = steady_clock::now();
        result = rectangle_l_2d(My::bar, My::bar_a, My::bar_b, My::bar_c, My::bar_d, aN[0]);
        calc_time = duration_cast<duration<double>>(steady_clock::now() - begTime).count();

        if (withError)
        {
            abs_err = abs_err_rect_rl_2d(My::bar, My::bar_a, My::bar_b, My::bar_c, My::bar_d, aN[0]);
            result_test = rectangle_l_2d(My::bar, My::bar_a, My::bar_b, My::bar_c, My::bar_d, aN_runge[0]);
            rung = runge(result_test, result, 1);
        }

        print_CSV(num_procs, 0, aN[0], calc_time, result, abs_err, rung);
        break;
    }
    case 1:
    {
        // Формула правых прямоугольников
        begTime = steady_clock::now();
        result = rectangle_r_2d(My::bar, My::bar_a, My::bar_b, My::bar_c, My::bar_d, aN[1]);
        calc_time = duration_cast<duration<double>>(steady_clock::now() - begTime).count();

        if (withError)
        {
            abs_err = abs_err_rect_rl_2d(My::bar, My::bar_a, My::bar_b, My::bar_c, My::bar_d, aN[0]);
            result_test = rectangle_r_2d(My::bar, My::bar_a, My::bar_b, My::bar_c, My::bar_d, aN_runge[1]);
            rung = runge(result_test, result, 1);
        }

        print_CSV(num_procs, 1, aN[1], calc_time, result, abs_err, rung);
        break;
    }
    case 2:
    {
        // Формула средних прямоугольников
        begTime = steady_clock::now();
        result = rectangle_m_2d(My::bar, My::bar_a, My::bar_b, My::bar_c, My::bar_d, aN[2]);
        calc_time = duration_cast<duration<double>>(steady_clock::now() - begTime).count();

        if (withError)
        {
            abs_err = abs_err_rect_m_2d(My::bar, My::bar_a, My::bar_b, My::bar_c, My::bar_d, aN[2]);
            result_test = rectangle_m_2d(My::bar, My::bar_a, My::bar_b, My::bar_c, My::bar_d, aN_runge[2]);
            rung = runge(result_test, result, 2);
        }

        print_CSV(num_procs, 2, aN[2], calc_time, result, abs_err, rung);
        break;
    }
    case 3:
    {
        // Формула трапеций
        begTime = steady_clock::now();
        result = trapezoidal_2d(My::bar, My::bar_a, My::bar_b, My::bar_c, My::bar_d, aN[3]);
        calc_time = duration_cast<duration<double>>(steady_clock::now() - begTime).count();

        if (withError)
        {
            abs_err = abs_err_trap_2d(My::bar, My::bar_a, My::bar_b, My::bar_c, My::bar_d, aN[3]);
            result_test = trapezoidal_2d(My::bar, My::bar_a, My::bar_b, My::bar_c, My::bar_d, aN_runge[3]);
            rung = runge(result_test, result, 2);
        }

        print_CSV(num_procs, 3, aN[3], calc_time, result, abs_err, rung);
        break;
    }
    case 4:
    {
        // Формула Симпсона
        begTime = steady_clock::now();
        result = simpson_2d(My::bar, My::bar_a, My::bar_b, My::bar_c, My::bar_d, aN[4]);
        calc_time = duration_cast<duration<double>>(steady_clock::now() - begTime).count();

        if (withError)
        {
            abs_err = abs_err_simps_2d(My::bar, My::bar_a, My::bar_b, My::bar_c, My::bar_d, aN[4]);
            result_test = simpson_2d(My::bar, My::bar_a, My::bar_b, My::bar_c, My::bar_d, aN_runge[4]);
            rung = runge(result_test, result, 4);
        }

        print_CSV(num_procs, 4, aN[4], calc_time, result, abs_err, rung);
        break;
    }
    case 5:
    {
        // Формула Ньютона-Котеса (правило трех восьмых)
        begTime = steady_clock::now();
        result = newton_38_2d(My::bar, My::bar_a, My::bar_b, My::bar_c, My::bar_d, aN[5]);
        calc_time = duration_cast<duration<double>>(steady_clock::now() - begTime).count();

        if (withError)
        {
            abs_err = abs_err_newton_38_2d(My::bar, My::bar_a, My::bar_b, My::bar_c, My::bar_d, aN[5]);
            result_test = newton_38_2d(My::bar, My::bar_a, My::bar_b, My::bar_c, My::bar_d, aN[5]);
            rung = runge(result_test, result, 4);
        }

        print_CSV(num_procs, 5, aN[5], calc_time, result, abs_err, rung);
        break;
    }
    case 6:
    {
        // Метод Монте-Карло
        begTime = steady_clock::now();
        result = monte_carlo_2d(My::bar, My::bar_a, My::bar_b, My::bar_c, My::bar_d, aN[6]);
        calc_time = duration_cast<duration<double>>(steady_clock::now() - begTime).count();

        if (withError)
        {
            abs_err = abs_err_monte_carlo_2d(My::bar, My::bar_a, My::bar_b, My::bar_c, My::bar_d, aN[6]);
            result_test = monte_carlo_2d(My::bar, My::bar_a, My::bar_b, My::bar_c, My::bar_d, aN_runge[6]);
            rung = runge(result_test, result, 4);
        }

        print_CSV(num_procs, 6, aN[6], calc_time, result, abs_err, rung);
        break;
    }
    default:
        break;
    }

    return 0;
}
