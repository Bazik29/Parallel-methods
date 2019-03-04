#include "utils.h"
#include <string>
#include <iostream>

void print_log(int id, int N, double result, double err, int N_runge, double err_runge, double calc_time)
{
    std::string names[] = {
        "Формула левых прямоугольников",
        "Формула правых прямоугольников",
        "Формула средних прямоугольников",
        "Формула трапеций",
        "Формула Симпсона",
        "\"Правило трех восьмых\""};

    std::cout << names[id] << std::endl;
    std::cout << "Число разбиений:              " << N << std::endl;
    std::cout << "Результат:                    " << result << std::endl;
    std::cout << "Погрешность:                  " << err << std::endl;
    std::cout << "Число разбиений для Рунге:    " << N_runge << std::endl;
    std::cout << "Погрешность по Рунге:         " << err_runge << std::endl;
    std::cout << "Время вычисления:             " << calc_time << std::endl;
    std::cout << std::endl;
}

void print_CSV(int threads, int id, int N, double calc_time, double result, double err, double err_runge)
{
    std::cout << threads << "," << id << "," << N << "," << calc_time << "," << result << "," << err << "," << err_runge << std::endl;
}

void file_CSV(std::ofstream &file, int threads, int id, int N, double calc_time, double result, double err, double err_runge)
{
    file << threads << "," << id << "," << N << "," << calc_time << "," << result << "," << err << "," << err_runge << std::endl;
}
