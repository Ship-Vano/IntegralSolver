//
// Created by Иван on 5/24/2024.
//

#ifndef CODE_INTEGRALPROBLEMSOLVER_H
#define CODE_INTEGRALPROBLEMSOLVER_H

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <cmath>
#include <functional>
#include "FileIO.h"
#include "IntegralProblem.h"
#include "SLAEsolver.h"

bool QuadratureScheme(const IntegralProblem &problem, const string &filename="UntitledTest");

bool IterativeScheme(const IntegralProblem &problem, const string &filename="UntitledTest");

bool DegenerateCoreScheme(const IntegralProblem &problem, const string &filename="UntitledTest");
#endif //CODE_INTEGRALPROBLEMSOLVER_H
