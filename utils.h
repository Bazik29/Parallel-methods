#pragma once
#include <string>
#include <fstream>

void print_log(std::string name, int N, double result, double err, int N_runge, double err_runge, double calc_time);

void print_CSV(int threads, int id, int N, double calc_time, double result, double err, double err_runge);

void file_CSV(std::ofstream &file, int threads, int id, int N, double calc_time, double result, double err, double err_runge);