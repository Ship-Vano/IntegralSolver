#include <iostream>
#include"IntegralProblemSolver.h"

void QuadtratureTest1(){
    IntegralProblem problem(0., 1., 0.0526316);
    problem.lambda = 0.5;
    problem.K = ([] (double x, double s){ return 1. - x * std::cos(x * s); });
    problem.K_isSet = true;
    problem.f = ([](double x) { return x * x + std::sqrt(x); });
    problem.f_isSet = true;
    problem.EPS = 1e-4;
    problem.EPS_is_set = true;
    QuadratureScheme(problem, "QuadtratureTest1.txt");
}

int main() {
    std::cout<<"hello" << std::endl;
    QuadtratureTest();
    return 0;
}
