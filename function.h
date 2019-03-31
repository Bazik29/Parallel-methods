#pragma once
#define _USE_MATH_DEFINES
#include <cmath>

namespace My
{

double foo_a = 1;
double foo_b = 2;

double foo_max_1 = 1.0000059608609866;
double foo_max_2 = 1.0000059608609866;
double foo_max_3 = 1.0000059608609866;
double foo_max_4 = 1.5491926629263049;

double foo_res = 1.3737208591045669;

double foo(double x)
{
    return x*(-std::atan((1.0/2.0)*x) + M_PI_2);
}

double foo_1(double x)
{
   return -1.0/2.0*x/((1.0/4.0)*std::pow(x, 2) + 1) - std::atan((1.0/2.0)*x) + M_PI_2;
}

double foo_2(double x)
{
    return 4*(std::pow(x, 2)/(std::pow(x, 2) + 4) - 1)/(std::pow(x, 2) + 4);
}

double foo_3(double x)
{
    return 16*(6*std::pow(x, 4)/std::pow(std::pow(x, 2) + 4, 2) - 7*std::pow(x, 2)/(std::pow(x, 2) + 4) + 1)/std::pow(std::pow(x, 2) + 4, 2);
}

double foo_4(double x)
{
    return 16*x*(-std::pow(x, 2)/(std::pow(x, 2) + 4) + 1)/std::pow(std::pow(x, 2) + 4, 2);
}

double bar_a = 1;
double bar_b = 3;
double bar_c = 0;
double bar_d = 3;

double bar_res = 13.36;

double bar(double x, double y)
{
    return ((y < x) && (y > x / 2)) ? std::pow(x + y, 2)/x : 0;
}

double bar_x(double x, double y)
{
    return ((y < x) && (y > x / 2)) ? (2*x + 2*y)/x - std::pow(x + y, 2)/std::pow(x, 2) : 0;
}

double bar_y(double x, double y)
{
    return ((y < x) && (y > x / 2)) ? (2*x + 2*y)/x : 0;
}

double bar_xx(double x, double y)
{
    return ((y < x) && (y > x / 2)) ? 2*(1 - 2*(x + y)/x + std::pow(x + y, 2)/std::pow(x, 2))/x : 0;
}

double bar_yy(double x, double y)
{
    return ((y < x) && (y > x / 2)) ? (2*x + 2*y)/x : 0;
}

double bar_xy(double x, double y)
{
    return ((y < x) && (y > x / 2)) ? 2*(1 - (x + y)/x)/x : 0;
}

} // namespace My