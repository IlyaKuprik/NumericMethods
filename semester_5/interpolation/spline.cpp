#define _CRT_SECURE_NO_DEPRECATE
#pragma once
#include "roots.cpp"
#include <bits/stdc++.h>
#include "../../Linalg.cpp"
#include "../../algs.cpp"
using namespace std;

vector<vector<long double>> spline10_coeffs(const vector<long double> roots, long double func(long double)){

    vector<vector<long double>> coefficents;
    for(int i = 0; i < roots.size()-1; ++i){
        long double x_i0 = roots[i], x_i1 = roots[i+1];
        long double y_i0 = func(x_i0), y_i1 = func(x_i1);
        long double a_i1 = (y_i1 - y_i0)/(x_i1 - x_i0);
        long double a_i0 = y_i0 - a_i1*x_i0;
        coefficents.push_back({a_i0, a_i1});
    }
    return coefficents;
}

vector<vector<long double>> spline32_coeffs(const vector<long double> roots, long double func(long double)){
    Vector gamma(roots.size()-2);
    Matrix H(roots.size()-2);
    for(int i = 0; i < roots.size()-2; ++i){
        long double x1 = roots[i], x2 = roots[i+1], x3 = roots[i+2];
        long double h1 = x2 - x1, h2 =  x3 - x2;
        H.matrix[i][i] = 2*(h1 + h2);
        if(i < roots.size()-3){
            H.matrix[i+1][i] = h2;
            H.matrix[i][i+1] = h2;
        }

        gamma.vec[i] = 6*((func(x3) - func(x2))/h2 - (func(x2) - func(x1))/h1);
    }

    vector<long double> y_dev2 = reflection_method(H, gamma).vec;
    y_dev2.push_back(0);
    y_dev2.insert(y_dev2.begin(),0);

    vector<long double> y_dev1(y_dev2.size());
    for(int i = 0; i < y_dev1.size() -1; ++i){
        long double h = roots[i+1] - roots[i];
        y_dev1[i] = (func(roots[i+1]) - func(roots[i]))/h - h*y_dev2[i+1]/6 - h*y_dev2[i]/3;
    }

    vector<vector<long double>> coefficents;
    for(int i = 0; i < y_dev2.size(); ++i)
        coefficents.push_back({y_dev1[i], y_dev2[i]});
    return coefficents;
}

vector<vector<long double>> spline21_coeffs(const vector<long double> roots, long double func(long double)){
    int matrix_size = roots.size()*3 - 3;
    Matrix A(matrix_size);
    Vector b(matrix_size);
    
    for(int i = 0; i < matrix_size; i+=3){
        long double curr_x = roots[i/3], next_x = roots[i/3+1];
        long double curr_y = func(curr_x), next_y = func(next_x);

        b.vec[i] = curr_y;
        b.vec[i + 1] = next_y;

        A.matrix[i][i] = curr_x*curr_x;
        A.matrix[i][i+1] = curr_x;
        A.matrix[i][i+2] = 1;

        A.matrix[i+1][i] = next_x*next_x;
        A.matrix[i+1][i+1] = next_x;
        A.matrix[i+1][i+2] = 1;

        A.matrix[i+2][i] = 2*next_x;
        A.matrix[i+2][i+1] = 1;
        if(i + 3 < matrix_size){
            A.matrix[i+2][i+3] = -2*next_x;
            A.matrix[i+2][i+4] = -1;
        }
    }

    vector<long double> a = reflection_method(A, b).vec;

    vector<vector<long double>> coefficents;
    for(int i = 0; i < matrix_size; i+=3) coefficents.push_back({a[i], a[i+1], a[i + 2]});
    return coefficents;
}


vector<vector<long double>> spline20_coeffs(const vector<long double> roots, long double func(long double)){
    vector<vector<long double>> coefficents;
    for(int k = 0; k < roots.size(); k+=2){
        Matrix A(3);
        Vector b(3);
        if(k+2 >= roots.size()) break;
        for(int i = 0; i < 3; ++i){
            b.vec[i] = func(roots[k+i]); 
            for(int j = 0; j < 3; ++j){
                A.matrix[i][j] = pow(roots[k+i], 2-j);
            }
        }
        vector<long double> c = reflection_method(A, b).vec;
        coefficents.push_back(c);
    }
    return coefficents;
}

long double spline_at_point(long double x, vector<vector<long double>> coefficents, const vector<long double> roots, long double func(long double), string method = "linear"){
    if(method == "spline_20"){
        for(int k = 0; k < roots.size(); k+=2){
            if((roots[k] < x || abs(x - roots[k]) < 1e-8) && (x < roots[k+2] || abs(x - roots[k+2]) < 1e-8))
                return coefficents[k/2][0]*x*x + coefficents[k/2][1]*x + coefficents[k/2][2];
        }
    }
    for(int i = 0;i < roots.size()-1; ++i){
        if((roots[i] < x || abs(x - roots[i]) < 1e-8) && (x < roots[i+1] || abs(x - roots[i+1]) < 1e-8)){
            if(method == "spline_10")
                return coefficents[i][0] + coefficents[i][1]*x;
            if(method == "spline_32"){
                return func(roots[i]) + coefficents[i][0]*(x - roots[i]) + coefficents[i][1]*pow(x - roots[i],2)/2 +
                        (coefficents[i+1][1] - coefficents[i][1])*pow(x - roots[i],3)/(6*(roots[i+1] - roots[i]));
            }
            if(method == "spline_21"){
                return coefficents[i][0]*x*x + coefficents[i][1]*x + coefficents[i][2];
            }
            else return 0;
        }
    }
    return 0;
}

pair<vector<long double>, vector<long double>> spline_on_interval(int n_dots, long double a, long double b,
                            int n_roots, long double func(long double), vector<long double> roots, string method){
   
    vector<long double> dots(n_dots);
    vector<long double> values(n_dots);

    vector<vector<long double>> coefficents;
    if(method == "spline_10")
        coefficents = spline10_coeffs(roots, func);
    else if(method == "spline_32")
        coefficents = spline32_coeffs(roots, func);
    else if(method == "spline_20")
        coefficents = spline20_coeffs(roots, func);
    else if(method == "spline_21")
        coefficents = spline21_coeffs(roots, func);
    a = roots[0];
    b = roots[roots.size()-1];
    for(int i = 0; i < n_dots; ++i){
        long double x = a + i*((b-a)/(n_dots - 1));
        dots[i] = x;

        values[i] = spline_at_point(x, coefficents, roots, func, method);
        cout << x << " " << values[i] << '\n';
    }
    return make_pair(dots, values);
}