#include <iomanip>
#include "Color.h"

Color::Color() {
    this->r = 0;
    this->g = 0;
    this->b = 0;
}

Color::Color(double r, double g, double b)
{
    this->r = r;
    this->g = g;
    this->b = b;
}

Color::Color(const Color &other)
{
    this->r = other.r;
    this->g = other.g;
    this->b = other.b;
}

