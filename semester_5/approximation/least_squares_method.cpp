#include <bits/stdc++.h>
#include "../../Linalg.cpp"
#include "../../algs.cpp"
using namespace std;

long double sum(vector<pair<long double, long double>> dots, int deg){
    long double su = 0;
    for(int i = 0; i < dots.size(); ++i){
        su += pow(dots[i].first, deg);
    }
    return su;
}

vector<long double> normal_equation(int n, vector<pair<long double, long double>> dots){
    Matrix X(n+1);
    Vector f(n+1);
    int m = dots.size();
    for(int i = 0;i < n+1; ++i){
        for(int k = 0; k < dots.size(); ++k)
            f.vec[i] += dots[k].second*pow(dots[k].first, i);
        for(int j = 0; j < n+1; ++j){
            X.matrix[i][j] = sum(dots, i + j);
        }
    }
    
    vector<long double> a = reflection_method(X, f).vec;
    return a;
}

long double alpha(int j_1, vector<pair<long double, long double>> dots, vector<vector<long double>> q_matrix){
    int j = j_1 - 1;
    long double numerator = 0;
    long double denominator = 0;
    for(int i = 0; i < dots.size(); ++i){
        long double tmp = pow(q_matrix[j][i],2);
        numerator += dots[i].first*tmp;
        denominator += tmp;
    }
    return numerator/denominator;
}

long double betta(int j, vector<pair<long double, long double>> dots, vector<vector<long double>> q_matrix){
    long double numerator = 0;
    long double denominator = 0;
    for(int i = 0; i < dots.size(); ++i){
        long double tmp = q_matrix[j-1][i];
        numerator += dots[i].first*q_matrix[j][i]*tmp;
        denominator += tmp*tmp;
    }
    return numerator/denominator;
}

long double q(int j,long double x, vector<pair<long double, long double>> dots, vector<vector<long double>> q_matrix, int for_matrix){
    if(j == 0) return 1;
    if(j == 1) return x - (1./dots.size())*sum(dots, 1);
    return x*q_matrix[j-1][for_matrix] - alpha(j, dots, q_matrix)*q_matrix[j-1][for_matrix] - betta(j-1, dots, q_matrix)*q_matrix[j-2][for_matrix];
}

long double a_k(int k, vector<pair<long double, long double>> dots, vector<vector<long double>> q_matrix){
    long double numerator = 0;
    long double denominator = 0;
    for(int i = 0; i < dots.size(); ++i){
        //5.2.38
        long double tmp = q_matrix[k][i];
        numerator += tmp*dots[i].second;
        denominator += tmp*tmp;
    }
    return numerator/denominator;
}

long double ortoghonal_polynoms_at_point(long double x, int n, vector<pair<long double, long double>> dots, vector<vector<long double>> q_matrix,
                                        vector<long double> alpha_list, vector<long double> betta_list){
    vector<long double> q_x_list(n+1);
    q_x_list[0] = 1;
    q_x_list[1] = x - (1./dots.size())*sum(dots, 1);

    for(int i = 2; i < n+1; ++i){;
        q_x_list[i] = x*q_x_list[i-1] - alpha_list[i]*q_x_list[i-1] - betta_list[i-1]*q_x_list[i-2];
    }

    vector<long double> a_list(n+1);
    for(int i = 0; i < n+1; ++i)
        a_list[i] = a_k(i, dots, q_matrix);

    long double su = 0;

    for(int i = 0; i < n+1; ++i)
        su += a_list[i]*q_x_list[i];
    return su;
}

long double least_squares_method_at_point(long double x, vector<long double> coefficents){
    long double value = coefficents[0];
    for(int i = 1; i < coefficents.size(); ++i){
        value += pow(x, i)*coefficents[i];
    }
    return value;
}

long double least_squares_method_on_interval(int n, vector<pair<long double, long double>> dots, string method = "normal", bool print = true){
    vector<long double> coefficents;

    vector<vector<long double>> q_matrix(n+1, vector<long double>(dots.size()));
    vector<long double> alpha_list(n+1);
    vector<long double> betta_list(n+1);
    if(method == "normal")
        coefficents = normal_equation(n, dots);
    else{
        for(int i = 0;i < n+1; ++i)
            for(int j = 0;j < dots.size(); ++j)
                q_matrix[i][j] = q(i, dots[j].first, dots, q_matrix, j);

        for(int i = 1;i < n+1; ++i){
            alpha_list[i] = alpha(i+1, dots, q_matrix);
            betta_list[i] = betta(i, dots, q_matrix);
    }
    }
    long double error = 0;
    for(int i = 0; i < dots.size(); i++){
        long double x = dots[i].first;

        long double val;
        if(method == "normal")
            val = least_squares_method_at_point(x, coefficents);
        else
            val = ortoghonal_polynoms_at_point(x, n, dots, q_matrix, alpha_list, betta_list);
        error += pow(val - dots[i].second, 2);
        if(print)
            cout << x << " " << val << endl;
    }
    return error;
}
