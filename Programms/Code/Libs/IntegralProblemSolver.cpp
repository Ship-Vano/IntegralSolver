//
// Created by Иван on 5/24/2024.
//
#include "IntegralProblemSolver.h"

/*ЛР 5:*/

bool QuadratureScheme(const IntegralProblem &problem, const string &filename) {

    // Создание файла
    std::string path = "./OutputData/" + filename;
    std::ofstream fpoints(path);
    std::cout << "log[INFO]: Starting QuadratureScheme" << std::endl;
    std::cout << "log[INFO]: Opening a file \"" << filename << "\" to write..." << std::endl;
    if (fpoints.is_open()) {

        double left_edge = problem.x0;
        double right_edge = problem.X;
        double hx = problem.hx;
        int num_steps = problem.num_x_steps;

        std::function<double(double, double)> Core;
        if(problem.K_isSet) {
             Core = problem.K;
        }
        else{
            std::cout << "Core function K is not set!" << std::endl;
            return false;
        }

        double lambda = problem.lambda;

        std::function<double(double)> f;
        if(problem.f_isSet) {
            f = problem.f;
        }
        else{
            std::cout << "Right function f is not set!" << std::endl;
            return false;
        }

        double EPS = 1e-4;
        if(problem.EPS_is_set){
            EPS = problem.EPS;
        }

        std::vector<std::vector<double>> A(num_steps+1, std::vector<double>(num_steps+1, 0.));
        std::vector<double> b(num_steps+1, 0.);
        std::vector<double> sol(num_steps+1, 0.);

        double x_i = left_edge;
        std::vector<double> net(num_steps+1, 0.);
        for(int i = 0; i < num_steps+1; ++i){
            net[i] = x_i;
            x_i += hx;
            for(int j = 0; j < num_steps+1; ++j){
                A[i][j] = -lambda * Core(x_i, left_edge + hx*j) * hx;
            }
            A[i][0] *= 0.5;
            A[i][num_steps] *= 0.5;
            A[i][i] += 1.;
            b[i] = f(x_i);
        }
        System<double> GaussSys(A, b);
        if(GaussSys.GaussianPartChoiceSolve(1e-3)) {
            sol = GaussSys.SolutionX;
        }
        else {
            std::cout << "LOG[WARN] The solution is Empty." << std::endl;
            return false;
        }
        writeVectorToFile(fpoints, net);
        writeVectorToFile(fpoints, sol);
        fpoints.close();
    }
    else {
        std::cout << "log[ERROR]: Couldn't open or create a file" << std::endl;
        return false;
    }
    return true;
}



bool IterativeScheme(const IntegralProblem &problem, const string &filename) {

    // Создание файла
    std::string path = "./OutputData/" + filename;
    std::ofstream fpoints(path);
    std::cout << "log[INFO]: Starting IterativeScheme" << std::endl;
    std::cout << "log[INFO]: Opening a file \"" << filename << "\" to write..." << std::endl;
    if (fpoints.is_open()) {

        double left_edge = problem.x0;
        double right_edge = problem.X;
        double hx = problem.hx;
        int num_steps = problem.num_x_steps;

        std::function<double(double, double)> Core;
        if(problem.K_isSet) {
            Core = problem.K;
        }
        else{
            std::cout << "Core function K is not set!" << std::endl;
            return false;
        }

        double lambda = problem.lambda;

        std::function<double(double)> f;
        if(problem.f_isSet) {
            f = problem.f;
        }
        else{
            std::cout << "Right function f is not set!" << std::endl;
            return false;
        }

        double EPS = 1e-4;
        if(problem.EPS_is_set){
            EPS = problem.EPS;
        }

        std::vector<std::vector<double>> A(num_steps+1, std::vector<double>(num_steps+1, 0.));
        std::vector<double> b(num_steps+1, 0.);
        std::vector<double> sol_k(num_steps+1, 0.);

        double x_i = left_edge;
        std::vector<double> net(num_steps+1, 0.);

        for(int i = 0; i < num_steps+1; ++i){
            net[i] = x_i;
            sol_k[i] = f(x_i);
            x_i += hx;
        }
        std::vector<double> sol_kp(sol_k);

        int its = 0;
        int max_allowed_its = 2000.;
        double temp_int_sum = 0.;
        double fL = 0.;
        double fR = 0.;
        double s = 0.;
        do{
            sol_kp.swap(sol_k);
            x_i = left_edge;
            for (int i = 0; i <= num_steps; ++i) {
                x_i += hx;
                temp_int_sum = 0.;
                fL = 0.,
                        fR = Core(x_i, left_edge) * sol_k[0];
                for (int j = 1; j <= num_steps; ++j) {
                    std::swap(fL, fR);
                    s = left_edge + hx * j;
                    fR = Core(x_i, s) * sol_k[j];
                    temp_int_sum += (fL + fR);
                }
                sol_kp[i] = f(x_i) + lambda * 0.5 * hx* temp_int_sum;
            }
            ++its;
        }while(its < max_allowed_its && vec_norm_inf(sol_k + (-1)*sol_kp) > EPS);

        writeVectorToFile(fpoints, net);
        writeVectorToFile(fpoints, sol_kp);

        fpoints.close();
    }
    else {
        std::cout << "log[ERROR]: Couldn't open or create a file" << std::endl;
        return false;
    }
    return true;
}
