#include <iostream>
#include"IntegralProblemSolver.h"
#include <typeinfo>


/* ЧАСТЬ 1 ТЕСТОВ ИЗ МЕТОДИЧКИ: КВАДРАТУРЫ И ИТЕРАЦИИ*/

// Тест1
void Test1(){
    IntegralProblem problem(0., 1., 0.05);
    problem.lambda = 0.5;
    problem.K = ([] (double x, double s){ return 1. - x * std::cos(x * s); });
    problem.K_isSet = true;
    problem.f = ([](double x) { return 0.5 * (1. + std::sin(x)); });
    problem.f_isSet = true;
    problem.EPS = 1e-7;
    problem.EPS_is_set = true;
    QuadratureScheme(problem, "QuadtratureTest1.txt");
    IterativeScheme(problem, "IterativeTest1.txt");

    /* Анал решение u(x) = 1 */
}

// Тест2
void Test2(){
    IntegralProblem problem(0., 1., 0.05);
    problem.lambda = 0.5;
    problem.K = ([] (double x, double s){ return 1. - x * std::cos(x * s); });
    problem.K_isSet = true;
    problem.f = ([](double x) { return x * x + std::sqrt(x);  });
    problem.f_isSet = true;
    problem.EPS = 1e-7;
    problem.EPS_is_set = true;
    QuadratureScheme(problem, "QuadtratureTest2.txt");
    IterativeScheme(problem, "IterativeTest2.txt");
}

// Тест3
void Test3(){
    IntegralProblem problem(0.1, 1., 0.05);
    problem.lambda = 0.5;
    problem.K = ([] (double x, double s){ return 1. - x * std::cos(x * s); });
    problem.K_isSet = true;
    problem.f = ([](double x) { return 0.5 * (1. + std::sin(x)); });
    problem.f_isSet = true;
    problem.EPS = 1e-7;
    problem.EPS_is_set = true;
    QuadratureScheme(problem, "QuadtratureTest3.txt");
    IterativeScheme(problem, "IterativeTest3.txt");
}

// Тест4
void Test4(){
    IntegralProblem problem(0.1, 1., 0.05);
    problem.lambda = 0.5;
    problem.K = ([] (double x, double s){ return 1. - x * std::cos(x * s); });
    problem.K_isSet = true;
    problem.f = ([](double x) { return x * x + std::sqrt(x);  });
    problem.f_isSet = true;
    problem.EPS = 1e-7;
    problem.EPS_is_set = true;
    QuadratureScheme(problem, "QuadtratureTest4.txt");
    IterativeScheme(problem, "IterativeTest4.txt");
}

/* ЧАСТЬ 2 ТЕСТОВ ИЗ МЕТОДИЧКИ: СИНГУЛЯРНОСТИ*/

void SingularTest1() {
    // параметры уравнения (любые переменные в конструктор --- они ни на что не влияют)
    IntegralProblem problem(0, 480, 10);
    double N = 5.;
    problem.f = ([=](double phi) { return std::sin((N + 1) * 0.5 * phi);});
    problem.f_isSet = true;
    problem.EPS = 1e-7;
    problem.EPS_is_set = true;
    problem.num_x_steps = 48; //число контрольных точек
    SingularScheme(problem, "SingularTest1.txt");
}

void SingularTest2() {
    // параметры уравнения (любые переменные в конструктор --- они ни на что не влияют)
    IntegralProblem problem(0, 480, 10);
    double N = 16.;
    problem.f = ([=](double phi) { return std::cos(N * 0.5 * phi);});
    problem.f_isSet = true;
    problem.EPS = 1e-7;
    problem.EPS_is_set = true;
    problem.num_x_steps = 48; //число контроьных точек
    SingularScheme(problem, "SingularTest2.txt");
}

void SingularTest3() {
    // параметры уравнения (любые переменные в конструктор --- они ни на что не влияют)
    IntegralProblem problem(0, 480, 10);
    //double N = 16.;
    problem.f = ([=](double phi) { return std::sin( 7 * phi);});
    problem.f_isSet = true;
    problem.EPS = 1e-7;
    problem.EPS_is_set = true;
    problem.num_x_steps = 5; //число контроьных точек
    SingularScheme(problem, "SingularTest3_1.txt");

    problem.num_x_steps = 10; //число контроьных точек
    SingularScheme(problem, "SingularTest3_2.txt");

    problem.num_x_steps = 15; //число контроьных точек
    SingularScheme(problem, "SingularTest3_3.txt");

    problem.num_x_steps = 20; //число контроьных точек
    SingularScheme(problem, "SingularTest3_4.txt");

    problem.num_x_steps = 25; //число контроьных точек
    SingularScheme(problem, "SingularTest3_5.txt");

    problem.num_x_steps = 30; //число контроьных точек
    SingularScheme(problem, "SingularTest3_6.txt");

}

/* КАСТОМНЫЕ ТЕСТЫ*/

// Проверка точности квадратурной формулы (вопрос 8 отчёта)
void QuadCheck(){
    double q = 0.5;
    double h = (std::exp(1.0) -1.)/30.;
    for(int i = 1; i < 4; ++i) {
        IntegralProblem problem(1., std::exp(1.0), h * std::pow(q,i-1));
        problem.lambda = std::exp(1.0);
        problem.K = ([](double x, double s) { return std::log(s) / x; });
        problem.K_isSet = true;
        problem.f = ([](double x) { return std::log(x); });
        problem.f_isSet = true;
        problem.EPS = 1e-17;
        problem.EPS_is_set = true;
        QuadratureScheme(problem, "QuadCheck" + to_string(i) + ".txt");
    }
}


// Проверка точности квадратурной формулы (вопрос 8 отчёта)
void StopCriterionCheck(){
    double eps = 1e-2;
    double its = 1;
    for(int i = 1; i < 13; ++i) {
        IntegralProblem problem(0., M_PI, M_PI/50.);
        problem.lambda = 1./(2*M_PI);
        problem.K = ([](double x, double s) { return s * std::sin(x); });
        problem.K_isSet = true;
        problem.f = ([](double x) { return std::cos(x); });
        problem.f_isSet = true;
        problem.EPS = 1e-17;
        problem.EPS_is_set = true;
        IterativeScheme(problem, "StopCriterionCheck" + to_string(i) + ".txt", its);
        ++its;
    }
}
// Вырожденные ядра
void DegenerateTest1(){
    // параметры уравнения
    IntegralProblem problem(0., 1., 0.025);
    problem.lambda = 1.;
    problem.f = ([](double x) { return 5.*x; });
    problem.f_isSet = true;
    problem.EPS = 1e-8;
    problem.EPS_is_set = true;

    // факторы ядра (множители)
    int amount_of_core_funcs = 3;
    std::vector<std::function<double(double)>> phi;
    phi.resize(amount_of_core_funcs);
    std::vector<std::function<double(double)>> psi;
    psi.resize(amount_of_core_funcs);
    phi[0] = ([](double x) { return 2. * std::sin(2.*M_PI * x); });
    psi[0] = ([](double s) { return std::cos(2.*M_PI * s); });
    phi[1] = ([](double x) { return -2.*std::cos(2.*M_PI * x); });
    psi[1] = ([](double s) { return std::sin(2.*M_PI * s); });
    phi[2] = ([](double x) { return -4.; });
    psi[2] = ([](double s) { return 1.; });

    problem.phi = phi;
    problem.phi_is_set = true;
    problem.psi = psi;
    problem.psi_is_set = true;

    DegenerateCoreScheme(problem, "DegenerateTest1.txt");

    problem.K = ([](double x, double s) { return 2. * (std::sin(2.*M_PI*(x-s))-2); });
    problem.K_isSet = true;
    QuadratureScheme(problem, "DegenerateTest1_quad.txt");
    IterativeScheme(problem, "DegenerateTest1_iter.txt");
}
/*void DegenerateTest0(){

    // параметры уравнения
    IntegralProblem problem(-M_PI, M_PI, 0.02*M_PI);
    problem.lambda = 1.;
    problem.f = ([](double x) { return std::sin(2.*x); });
    problem.f_isSet = true;
    problem.EPS = 1e-17;
    problem.EPS_is_set = true;

    // факторы ядра (множители)
    int amount_of_core_funcs = 2;
    std::vector<std::function<double(double)>> phi;
    phi.resize(amount_of_core_funcs);
    std::vector<std::function<double(double)>> psi;
    psi.resize(amount_of_core_funcs);
    phi[0] = ([](double x) { return 1./M_PI * std::sin(x); });
    psi[0] = ([](double s) { return std::sin(s); });
    phi[1] = ([](double x) { return 1.; });
    psi[1] = ([](double s) { return s; });

    problem.phi = phi;
    problem.phi_is_set = true;
    problem.psi = psi;
    problem.psi_is_set = true;

    DegenerateCoreScheme(problem, "DegenerateTest0.txt");

    problem.K = ([](double x, double s) { return 1./M_PI * std::sin(x) * std::sin(s) + s; });
    problem.K_isSet = true;
    QuadratureScheme(problem, "DegenerateTest0_quad.txt");
    IterativeScheme(problem, "DegenerateTest0_iter.txt");
}*/

void DegenerateTest2() {

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

    DegenerateCoreScheme(problem, "DegenerateTest2.txt");
}


/* Файлы для табличек эиткена для метода квадратур на основе теста 1*/
void make_data_ffor_tables_1() {

    double h0 = 0.1;
    for (int i = 0; i < 6; i++) {
        IntegralProblem problem(0., 1., h0);
        problem.lambda = 0.5;
        problem.K = ([](double x, double s) { return 1. - x * std::cos(x * s); });
        problem.K_isSet = true;
        problem.f = ([](double x) { return 0.5 * (1. + std::sin(x)); });
        problem.f_isSet = true;
        problem.EPS = 1e-7;
        problem.EPS_is_set = true;
        QuadratureScheme(problem, "/data_for_tables/QuadtratureTest_" + to_string(i) + ".txt");

        h0 = h0 / 2.;
    }
}



int main() {
    //std::setlocale(LC_ALL, "rus");
    std::cout<<"hello" << std::endl;
    //std::cout << typeid( 1.0/5).name() << std::endl;
    //Test1();
    //Test2();
    //Test3();
    //Test4();
    //SingularTest1();
    //SingularTest2();
    //SingularTest3();
    //DegenerateTest1();
    //QuadCheck();
    StopCriterionCheck();


    make_data_ffor_tables_1();

    return 0;
}
