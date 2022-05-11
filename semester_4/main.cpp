#include <bits/stdc++.h>
#include "Linalg.cpp"
#include "algs.cpp"
using namespace std;

double eps = 1e-6;

vector<pair<double,long long>> sim_iter;
vector<pair<double,long long>> seidel_iter;

vector<double> gauss_iter;
vector<double> reflection_iter;

//функция для вывода результат работы итерационных методов 
void calculate_iter_test(const Matrix & A, const Vector &b, const Vector &x, const vector<double> &eps_range, bool print = true){
    if(print)
        cout << "|===========МПИ==========|=========Зейдель=========|\n";
    for(auto&e:eps_range){
        if(print)
            cout << "e = " << e << '\n';
        pair<Vector,long long> simp_iter =  simple_iteration(A,b,e);
        pair<Vector,long long> seidel =  seidel_method(A,b,e);
        Vector d1 = simp_iter.first - x;
        Vector d2 = seidel.first - x;
        sim_iter.push_back(make_pair(d1.get_norm2(), simp_iter.second));
        seidel_iter.push_back(make_pair(d2.get_norm1(), seidel.second));

        if(print){
            printf("x:\n");
            for(int i = 0; i < seidel.first.size; ++i){
                printf("\t%.7f\t\t%.7f\n",simp_iter.first.vec[i],seidel.first.vec[i]);
            }
            printf("Δ:\n");
            for(int i = 0; i < seidel.first.size; ++i){
                printf("\t%.7f\t\t%.7f\n",abs(simp_iter.first.vec[i] - x.vec[i]),abs(seidel.first.vec[i] - x.vec[i]));
            }
            printf("k:\n");
            printf("\t   %lld\t\t\t%lld\n", simp_iter.second, seidel.second);
        }
    }
}

//функция для вывода результатов работы точных методов
void calculate_accurate_test(const Matrix & A, const Vector &b, const Vector &x, bool print = true){
    if(print)
        cout << "|==========Гаусс=========|========Отражений========|\n";
    Vector gauss =  gauss_method(A,b);
    Vector reflection =  reflection_method(A,b);
    Vector d1 = gauss - x;
    Vector d2 = reflection - x;
    gauss_iter.push_back(d1.get_norm2());
    reflection_iter.push_back(d2.get_norm2());

    if(print){
        printf("x:\n");
        for(int i = 0; i < gauss.size; ++i){
            printf("\t%.7f\t\t%.7f\n",gauss.vec[i],reflection.vec[i]);
        }
        printf("Δ:\n");
        for(int i = 0; i < gauss.size; ++i){
            printf("\t%.7f\t\t%.7f\n",abs(gauss.vec[i] - x.vec[i]),abs(reflection.vec[i] - x.vec[i]));
        }
    }
}

//функция для вывода таблицы
void calculate_table(const vector<pair<Matrix,Vector>> &tests, const vector<Vector> &solutions, const vector<double> &iter_e_range = {1e-2,1e-3,1e-4,1e-5,1e-6}, bool print = true){
    for(int i = 0; i < tests.size(); ++i){
        if(print){
            cout << "Тест " << i << ":\n";
            cout << "Точное решение: " << solutions[i];
        }
        //Чтобы таблица не была громоздской, в качестве последнего аргумента можно передать список точностей, 
        //с которыми будут вычисляться итерационные методы
        calculate_iter_test(tests[i].first, tests[i].second, solutions[i], iter_e_range, print);
        calculate_accurate_test(tests[i].first, tests[i].second, solutions[i], print);
    }
}

//генерация плохо обусловленной матрицы для 5го теста 
Matrix get_bad_matrix(int size, double eps){
    Matrix A(size);
    Matrix A2(size);

    for(int i = 0;i < size; ++i){
        for(int j = 0; j < size; ++j){
            if(j == i)A.matrix[i][j] = 1;
            else if(j > i)A.matrix[i][j] = -1;
            else A2.matrix[i][j] = 1;
        }
    }

    A2 = A2 + A;
    A = A + A2*eps*size;
    return A;
}

//генерация правой части СЛАУ для 5го теста
Vector get_bad_vector(int size){
    Vector b(size);
    for(int i = 0;i < size; ++i)
        b.vec[i] = -1;
    b.vec[size-1] = 1;
    return b;
}

//вывод результата работы алгоритмов для теста 5
void test5(const vector<double> eps_arr = {1e-3,1e-6}, double print = true){
    cout << "Тест 5:\n";
    for(auto &eps:eps_arr)
        for(int size = 3; size < 11; size++){
            //генерация ответа
            Vector x(size);
            x.vec[size-1] = 1;
            // генерация данных из задания
            Matrix A = get_bad_matrix(size, eps); 
            Vector b = get_bad_vector(size);

            cout << "N = " << size << '\n';
            cout << "Точное решение: " << x;
            cout << "eps = " << eps << '\n';
            calculate_iter_test(A,b,x,{1e-6}, print);
            calculate_accurate_test(A,b,x, print);
        }
}

//инициализация тестов и ответов для тестов 0..4
void init_tests(vector<pair<Matrix,Vector>> & tests, vector<Vector> & solutions){
    double n = 9;
    //Точные решения получены с помощью стороннего сервиса: https://math.semestr.ru/gauss/gauss.php
    //Тест 0
    tests.push_back(make_pair(
        Matrix({{0,2,3}, {1,2,4}, {4,5,6}}),
        Vector({13,17,32})));
    solutions.push_back(Vector({1,2,3}));
    //Тест 1
    tests.push_back(make_pair(
        Matrix({{n+2,1,1}, {1,n+4,1}, {1,1,n+6}}),
        Vector({n+4,n+6,n+8})));
    solutions.push_back(Vector({1, 1, 1}));
    //Тест 2  
    tests.push_back(make_pair(
        Matrix({{-n-2,1,1}, {1,-n-4,1}, {1,1,-n-6}}),
        Vector({-n-4,-n-6,-n-8})));
    solutions.push_back(Vector({1.42585551, 1.36501901, 1.31939163}));
    //Тест 3  
    tests.push_back(make_pair(
        Matrix({{-n-2,n+1,n+1}, {n+1,-n-4,n+1}, {n+1,n+1,-n-6}}),
        Vector({n+4,n+6,n+8})));
    solutions.push_back(Vector({2.36884154, 2.0758988, 1.8298269}));
    //Тест 4
    tests.push_back(make_pair(
        Matrix({{n+2,n+1,n+1}, {n+1,n+4,n+1}, {n+1,n+1,n+6}}),
        Vector({n+4,n+6,n+8})));
    solutions.push_back(Vector({-0.10204082, 0.63265306, 0.77959184}));
}

void print_norm_arr(){
    cout << "SIM\t\tZEIDEL\t\tGAUSS\t\tREFLECT\n";
    for(int i = 0; i < gauss_iter.size(); ++i){
        printf("%0.7f\t\t%0.7f\t\t%0.7f\t\t%0.7f\n",sim_iter[i].first, seidel_iter[i].first, gauss_iter[i], reflection_iter[i]);
    }   
}

void print_iter_arr(){
    cout << "SIM\t\tZEIDEL\n";
    for(int i = 0; i < gauss_iter.size(); ++i){
        printf("%d\t\t%d\n",sim_iter[i].second, seidel_iter[i].second);
    }   
}

int main(){
    vector<pair<Matrix,Vector>> tests;
    vector<Vector> solutions;
    init_tests(tests, solutions);
    calculate_table(tests, solutions);
    //Последний аргумент - массив точностей, с которыми будут считаться итерационные методы
    //Если в качестве третьего аргумента функции ничего не передавать, то итерационные методы будут вычисляться с 
    //c точностями {1e-2,1e-3,1e-4,1e-5,1e-6}
    test5({1e-6});
    //В test5() можно передать массив констант eps, которая фигурирует в 5-ом тесте. По умолчанию это {1e-3,1e-6}
    return 0;   
}