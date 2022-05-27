#pragma once
#include <bits/stdc++.h>
#include "Linalg.cpp"
using namespace std;

static const int N = 9;

long double func(Vector vec){
    long double x = vec.vec[0], y = vec.vec[1], z = vec.vec[2];
    return 2*x*x + (3 + 0.1*N)*y*y + (4 + 0.1*N)*z*z + x*y - y*z + x*z + x - 2*y + 3*z + N; 
}

pair<Vector,int> method_of_gradiend_steepest_descent(const Matrix &A, const Vector &b, const long double &eps = 1e-6){
    Vector x_pred(A.size);
    cout << "\tНачальное приближение: " << x_pred;
    int num_of_iter = 0;
    while(true){
        Vector q = A*x_pred + b;
        long double mu = (-1.0)*q.get_norm2()*q.get_norm2()/((A*q)*q);
        Vector x_curr = x_pred + q*mu;
        Vector dx = x_pred - x_curr;
        if(dx.get_norm2() < eps){
            x_pred = x_curr;
            break;
        }
        x_pred = x_curr;
        num_of_iter++;
    }
    return {x_pred, num_of_iter};
}

pair<Vector,int> method_of_coord_steepest_descent(const Matrix &A, const Vector &b, const long double &eps = 1e-6){
    Vector x_pred(A.size);
    cout << "\tНачальное приближение: " << x_pred;
    int num_of_iter = 0;
    while(true){
        Vector x_curr = x_pred;
        long double min_value = 1e9;

        for(int i = 0;i < A.size; ++i){
            Vector q = A*x_pred + b;
            Vector e_i = get_basis_vector(A.size, i);
            long double mu = (-1.0)*(e_i*q)/((A*e_i)*e_i);
            Vector x_ = x_pred + e_i*mu;
            if(min_value > func(x_)){
                min_value = func(x_);
                x_curr = x_;
            }
        }

        Vector dx = x_pred - x_curr;
        if(dx.get_norm2() < eps){
            x_pred = x_curr;
            break;
        }

        x_pred = x_curr;
        num_of_iter++;
    }
    return {x_pred, num_of_iter};
}

int main(){
    Matrix A({{4, 1, 1}, {1, 7.8, -1}, {1, -1, 9.8}});
    Vector b({1,-2,3});
    Vector answer({-0.25035440884604479727814006237596, 0.25588318684434363481712503544088, -0.25446555146016444570456478593706});
    //Matrix A({{2,1}, {1,2}});
    //Vector b({-2, 1});
    //Vector answer({5./3, -4./3});
    cout << "N: " << N << '\n';
    cout << "Вектор b: " << b << '\n';
    cout << "Матрица A:\n" << A << '\n';
    cout << "Градиентный метод:\n";
    auto ans_1 = method_of_gradiend_steepest_descent(A, b, 1e-6);
    cout << "\tТочка минимума: " << ans_1.first << "\tЗначение функции: " << func(ans_1.first) << '\n';
    cout << "\tАбсолютная погрешность: " << ans_1.first - answer;
    cout << "\tКоличество итераций: " << ans_1.second << '\n'; 
    cout << "Покоординатный метод:\n";
    auto ans_2 = method_of_coord_steepest_descent(A, b, 1e-6);
    cout << "\tТочка минимума: " << ans_2.first << "\tЗначение функции: " << func(ans_2.first) << '\n';
    cout << "\tАбсолютная погрешность: " << ans_2.first - answer;
    cout << "\tКоличество итераций: " << ans_2.second << '\n';
    
    //answer: ()
    return 0;
}