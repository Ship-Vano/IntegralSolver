//
// Created by Иван on 5/24/2024.
//

#ifndef CODE_INTEGRALPROBLEM_H
#define CODE_INTEGRALPROBLEM_H

#include <cmath>
#include <iostream>
#include <functional>
#include "Libs/algebra.h"

class IntegralProblem {
public:

    /// Параметры задачи

    double lambda = 1.;

    // Функция правой части уравнения Фредгольма f(x)
    std::function<double(const double&)> f;
    bool  f_isSet = false;

    // Ядро оператора K(x,s)
    std::function<double(const double&, const double&)> K;
    bool  K_isSet = false;

    // Параметры расчёта
    double x0;  // начало отсчёта по x
    double X;   // координата правого конца по x
    double hx;  // шаг по x-координате
    double EPS; // точность расчётов
    bool EPS_is_set = false;
    int num_x_steps; // количество шагов по пространству x

    IntegralProblem(const double &x0_init, const double &X_init, const double &hx_init);

};

#endif //CODE_INTEGRALPROBLEM_H


