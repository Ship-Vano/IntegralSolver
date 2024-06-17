//
// Created by Иван on 5/24/2024.
//
#include "IntegralProblemSolver.h"

/*ЛР 5:*/

// Метод квадратур
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

// Метод простой итерации
bool IterativeScheme(const IntegralProblem &problem, const string &filename, const int& max_allowed_its) {

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
       // int max_allowed_its = 2000.;
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
                fL = 0.;
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
        fpoints << its << std::endl;
        fpoints.close();
    }
    else {
        std::cout << "log[ERROR]: Couldn't open or create a file" << std::endl;
        return false;
    }
    return true;
}

// Метод вырожденного ядра
bool DegenerateCoreScheme(const IntegralProblem &problem, const string &filename) {

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

        double x_i = left_edge;
        std::vector<double> net(num_steps+1, 0.);

        // настройка метода

        // Функции \phi_i
        std::vector<std::function<double(double)>> phi;

        // Функции \psi_i
        std::vector<std::function<double(double)>> psi;
        int amount_of_core_funcs = 0;
        if(problem.phi_is_set && problem.psi_is_set){
            phi = problem.phi;
            psi = problem.psi;
            amount_of_core_funcs = std::min(phi.size(), psi.size());
        }
        else{
            std::cout << "Phis and psis were not set! Can't determine a degenerate core." << std::endl;
            return false;
        }

        std::vector<std::vector<double>> A(amount_of_core_funcs, std::vector<double>(amount_of_core_funcs, 0.));
        std::vector<double> b(amount_of_core_funcs, 0.);
        std::vector<double> solution(num_steps+1, 0.);

        // основной блок решения
        double temp_int_sum = 0.;
        double fL = 0.;
        double fR = 0.;
        double s = 0.;
        for (int i = 0; i < amount_of_core_funcs; ++i) {
            for (int j = 0; j < amount_of_core_funcs; ++j) {
                temp_int_sum = 0.;
                fL = 0.;
                fR = psi[i](left_edge) * phi[i](left_edge);
                for(int k = 1; k <= num_steps; ++k){
                    std::swap(fL, fR);
                    s = left_edge + hx * k;
                    fR = psi[i](s) * phi[j](s);
                    //fL = psi[i](s - hx) * phi[j](s - hx);
                    temp_int_sum += (fL + fR);
                }
                A[i][j] = - lambda * 0.5 * hx * temp_int_sum;
            }
            temp_int_sum = 0.;
            fL = 0.;
            fR = psi[i](left_edge) * phi[i](left_edge);
            for(int k = 1; k <= num_steps; ++k){
                std::swap(fL, fR);
                s = left_edge + hx * k;
                fR = psi[i](s) * f(s);
                //fL = psi[i](s - hx) * f(s - hx);
                temp_int_sum += (fL + fR);
            }
            b[i] = 0.5 * hx * temp_int_sum;
            A[i][i] += 1.;
        }
//        for(int i = 0; i <= num_steps; ++i){
//                out(A[i]);
//        }
        System<double> GaussSys(A, b);
//        for(int i = 0; i < amount_of_core_funcs; ++i){
//            out(GaussSys.MatrixA[i]);
//        }

        GaussSys.GaussianPartChoiceSolve(EPS);
        std::vector<double> G_sol = GaussSys.SolutionX;
        for(int i = 0; i < num_steps+1; ++i){
            net[i] = x_i;
            temp_int_sum = 0.;
            for(int j = 0; j < amount_of_core_funcs; ++j){
                temp_int_sum += G_sol[j] * phi[j](x_i);
            }
            solution[i] = f(x_i) + lambda * temp_int_sum;
            x_i += hx;
        }

        writeVectorToFile(fpoints, net);
        writeVectorToFile(fpoints, solution);

        fpoints.close();
    }
    else {
        std::cout << "log[ERROR]: Couldn't open or create a file" << std::endl;
        return false;
    }
    return true;
}

// Метод для ингулярных уравнений
bool SingularScheme(const IntegralProblem &problem, const string &filename) {

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

        double x_i = left_edge;
        std::vector<double> net(num_steps+1, 0.);

        // настройка метода
        int N = num_steps + 1;
        std::vector<std::vector<double>> A(N+1, std::vector<double>(N+1, 0.));
        std::vector<double> b(N+1, 0.);
        std::vector<double> sol(N, 0.);
        double angle = 0.;
        // основной блок
        for(int i = 0; i < N; ++i){
            for(int j = 0; j < N; ++j){
                angle = (M_PI + 2. * M_PI * i - 2. * M_PI * j) / (N*1.0);
                A[i][j] = std::sin(angle) / (2. * N * (-1. + std::cos(angle)));
            }
            A[i][N] = 1.;
            b[i] = f(2.*M_PI * (i+0.5) / (N*1.0));
        }
        for(int j = 0; j < N; ++j){
            A[N][j] = 1.;
        }

        System<double> GaussSys(A, b);
        GaussSys.GaussianPartChoiceSolve(EPS);

        double R = GaussSys.SolutionX.back();
        for(int i = 0; i < N; ++i)
        {
            sol[i] = GaussSys.SolutionX[i];
            net[i] = 2 * M_PI * (double)i / (double)N;
        }
        std::cout << "R = " << R << std::endl;

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
