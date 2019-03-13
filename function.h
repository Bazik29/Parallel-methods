#pragma once

#include <cmath>

namespace My
{

const double PI = 3.141592653589793238463;

// foo: 1.373720859104567
// bar: ?

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

double bar_a = 1;
double bar_b = 3;
double bar_c = 0;
double bar_d = 3;
double bar(double x, double y)
{
    return ((y < x) && (y > x / 2)) ? std::pow(x + y, 2) / x : 0;
}

double bar_1(double x, double y)
{
    // TO DO
    return 0;
}

double bar_2(double x, double y)
{
    // TO DO
    return 0;
}

double bar_4(double x, double y)
{
    // TO DO
    return 0;
}

} // namespace My
