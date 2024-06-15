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

void IterativeTest1(){
    IntegralProblem problem(0., 1., 0.0526316);
    problem.lambda = 0.5;
    problem.K = ([] (double x, double s){ return 1. - x * std::cos(x * s); });
    problem.K_isSet = true;
    problem.f = ([](double x) { return x * x + std::sqrt(x); });
    problem.f_isSet = true;
    problem.EPS = 1e-4;
    problem.EPS_is_set = true;
    IterativeScheme(problem, "IterativeTest1.txt");
}

void DegenerateTest1() {

    // параметры уравнения
    IntegralProblem problem(0., 0.5*M_PI, 0.002);
    problem.lambda = 2.;
    problem.f = ([](double x) { return x; });
    problem.f_isSet = true;
    problem.EPS = 1e-9;
    problem.EPS_is_set = true;

    // факторы ядра (множители)
    int amount_of_core_funcs = 11;
    std::vector<std::function<double(double)>> phi;
    phi.resize(amount_of_core_funcs);
    std::vector<std::function<double(double)>> psi;
    psi.resize(amount_of_core_funcs);
    double factorial_placeholder = 1.;
    bool is_negtv = false;
    bool is_cos = false;
    for(int i = 0; i < amount_of_core_funcs; ++i){
        is_cos = false;
        is_negtv = false;
        if(i/2 % 2 == 1){
            is_negtv = true;
        }
        if(i%2 == 0){
            is_cos = true;
        }
        //std::cout << (is_negtv ? "-" : "+") << "1 / " << factorial_placeholder << " " << (is_cos ? "cos" : "sin") << "(s) x^" << i << "\n";
        psi[i] = [=](double s) { return (is_cos ? std::cos(s) : std::sin(s)); };
        phi[i] = [=](double x) { return (is_negtv ? -1 : 1) * std::pow(x, i) / factorial_placeholder; };
        factorial_placeholder *= (i + 1);
    }

    problem.phi = phi;
    problem.phi_is_set = true;
    problem.psi = psi;
    problem.psi_is_set = true;

    DegenerateCoreScheme(problem, "DegenerateTest1.txt");
}


int main() {
    //std::setlocale(LC_ALL, "rus");
    std::cout<<"hello" << std::endl;
    QuadtratureTest1();
    IterativeTest1();
    //DegenerateTest1();
    return 0;
}
