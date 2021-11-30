#include <sstream>
#include <iomanip>
#include <cmath>
#include "includes/CodeUtils/numbers.hpp"

bool codeutils::isNumberIntegral(double number)
{
    double difference = std::fabs(number - (std::round(number)));
    return ((difference < 0.00001) || (difference > 0.99999));
}
