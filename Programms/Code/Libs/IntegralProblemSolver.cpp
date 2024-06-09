//
// Created by Иван on 5/24/2024.
//
#include "IntegralProblemSolver.h"

/*ЛР 5:*/

bool QuadratureSсheme(const IntegralProblem &problem, const string &filename) {

    // Создание файла
    std::string path = "./OutputData/" + filename;
    std::ofstream fpoints(path);
    std::cout << "log[INFO]: Starting QuadratureSсheme" << std::endl;
    std::cout << "log[INFO]: Opening a file \"" << filename << "\" to write..." << std::endl;
    if (fpoints.is_open()) {



        fpoints.close();
    }
    else {
        std::cout << "log[ERROR]: Couldn't open or create a file" << std::endl;
        return false;
    }
    return true;
}

