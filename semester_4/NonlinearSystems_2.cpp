#include <iostream>
#include <cmath>
#include "Linalg.cpp"
#include "algs.cpp"
using namespace std;

Matrix Jacobi_0(double xk, double yk, double k = 1){
    return Matrix({{k*cos(xk + 0.5), -1}, {1, -k*sin(yk-2)}});
}

Vector F_0(double xk, double yk, double k = 1){
    return Vector({k*sin(xk + 0.5) - yk - 1, xk + k*cos(yk-2) - 2});
}

pair<int, Vector> nonlinear_newthon_method(Vector apr_x, Matrix Jacobi(double, double, double), Vector F(double, double, double), double k = 1, double eps = 1e-4, bool continue_iter = true){
    double xk = apr_x.vec[0], yk = apr_x.vec[1];
    int iter_n = 0;
    do{
        iter_n++;
        Matrix dF = Jacobi(xk, yk, k);
        Vector mF = F(xk, yk, k);
        Vector dx = reflection_method(dF, mF*(-1));
        double gn = dx.vec[0], hn = dx.vec[1];
        double xk1 = xk + gn;
        double yk1 = yk + hn;

        if(abs(yk - yk1) < eps && abs(xk - xk1) < eps){
            yk = yk1;
            xk = xk1;
            break;
        }

        yk = yk1;
        xk = xk1;
    }while(continue_iter);

    return make_pair(iter_n, Vector({xk, yk}));
}

Vector get_note_appr(int N){
    printf("N = %d\n",N);
    double x = 2,y = -1;
    for(int k = 1; k <= N; ++k){//когда строим нужно одну итерацию 
        pair<int, Vector> temp = nonlinear_newthon_method(Vector({x, y}),Jacobi_0, F_0, (1.0*k)/N, 1e-6, false);
        Vector xk = temp.second;
        x = xk.vec[0];
        y = xk.vec[1];
    }
    return Vector({x, y});
}

int main(){
    // для того чтобы получить верный ответ, использовал сервис https://www.kontrolnaya-rabota.ru/s/equal-many/system-any/
    //x = 2.981636760368604
    //y = -1.333528673687926

    //графический метод (график можно найти в папке с заднием)
    pair<int, Vector> temp = nonlinear_newthon_method(Vector({3, -1}), Jacobi_0, F_0, 1.0, 1e-6);
    Vector x = temp.second;
    printf("С помощью приближения было получено решение (%.7f, %.7f), за %d итераций\n", x.vec[0], x.vec[1], temp.first);
    
    // В этой функции используется замечание 1.4
    Vector appr = get_note_appr(10);
    temp = nonlinear_newthon_method(Vector({appr.vec[0], appr.vec[1]}), Jacobi_0, F_0, 1.0, 1e-6);
    x = temp.second;
    printf("Решение с использованием замечания 1.4: (%.7f, %.7f), за %d итераций\n \n", x.vec[0], x.vec[1],temp.first );
    return 0;
}