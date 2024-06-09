//
// Created by Иван on 5/24/2024.
//

#include "IntegralProblem.h"

IntegralProblem::IntegralProblem(const double &x0_init, const double &X_init, const double &hx_init) {
    x0 = x0_init;
    X = X_init;
    hx = hx_init;
    if(hx == 0){
        std::cout << "LOG[ERROR] hx or hy is set to 0!!!" << std::endl;
    }
    num_x_steps = static_cast<int>((X - x0)/hx);
}
