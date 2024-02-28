
#include <iostream>
#include <math.h>

using namespace std;

double function_1(double x) {//y= |_frac_{{1};{({x^{5}+x^{4}+1})}}-log({x})
    return 1 / (1 + (pow(x, 5) + pow(x, 4) + 1)) - log(x) / log(10);
}

double function_2(double x) {//y=x^{5}+0.1x^{4}+0.4abs({x})-1.2
    return pow(x, 5) + 0.1 * pow(x, 4) + 0.4 * abs(x) - 1.2;
}

int main() {


}