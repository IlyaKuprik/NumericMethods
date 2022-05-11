#pragma once
#include <bits/stdc++.h>
#include "Linalg.cpp"
using namespace std;


pair<Vector,long long> simple_iteration(Matrix A, Vector b, double eps = 1e-6){

    //Инициализация еденичной матрицы
    Matrix E(A.size);
    for(int i = 0; i < A.size; ++i){
        E.matrix[i][i] = 1;
    }

    double u = 1/(A.get_norm2());
    Matrix B = E - A*u;
    double k;
    //Проверка нормы матрицы B
    if(B.get_norm1() < 1){
        k = B.get_norm1()/(1 - B.get_norm1());
    }else if(B.get_norm2() < 1){
        k = B.get_norm2()/(1 - B.get_norm2());
    }else if(B.get_norm_inf() < 1){
        k = B.get_norm_inf()/(1 - B.get_norm_inf());
    }else{
        //Приведение матрицы к ПО виду
        Matrix At = A.transpose();
        A = At*A;
        b = At*b;
        u = 1/(A.get_norm2());
        B = E - A*u;
        k = B.get_norm1()/(1 - B.get_norm1());
    }
    Vector c = b*u;
    Vector x = c;
    long long num_of_iter = 0;
    while(true){
        Vector x_1 = B*x + c;
        Vector dx = A*x_1 - b;
        x = x_1;
        num_of_iter++;
        //Использую более общий критерий остановки(т.к не получилось реализовать спектральную норму)
        if(dx.get_norm2() <= eps)
            break;
    }
    
    return make_pair(x,num_of_iter);
}

pair<Vector,long long> seidel_method(Matrix A, Vector b, double eps = 1e-6){
    int size = A.size;
    bool is_diagonal_dominated = true;
    //проверка на диагональное преобладание
    for(int i = 0;i < size; ++i){
        double su = 0;
        for(int j = 0;j < size; ++j){
            if(j == i)continue;
            su += abs(A.matrix[i][j]);
        }
        if(abs(A.matrix[i][i]) - su < 0) is_diagonal_dominated = false;
    }
    //если матрица не имеет диагонального преобладания, то делаем её положительно определенной
    if(!is_diagonal_dominated){
        Matrix At = A.transpose();
        A = At*A;
        b = At*b;
    }

    
    Matrix C(size);
    Vector d(size);
    //Инициализация матрицы C
    for(int i = 0;i < size; ++i){
        d.vec[i] = b.vec[i]/A.matrix[i][i];
        for(int j = 0;j < size; ++j){
            if(i == j)
                C.matrix[i][j] = 0;
            else
                C.matrix[i][j] = (-1)*A.matrix[i][j]/A.matrix[i][i];
        }
    }

    Vector x0 = d;
    // Цикл для построения решения
    long long num_of_iter = 0;
    while(true){
        Vector x_k(size);
        for(int j = 0; j < size; ++j){
            x_k.vec[0] += C.matrix[0][j]*x0.vec[j];
        }
        x_k.vec[0] += d.vec[0];
        for(int i = 1;i < size; ++i){
            for(int j = 0; j < i; ++j){
                x_k.vec[i] += x_k.vec[j]*C.matrix[i][j];
            }
            for(int j = i+1; j < size; ++j){
                x_k.vec[i] += x0.vec[j]*C.matrix[i][j];
            }
            x_k.vec[i] += d.vec[i];
        }
        x0  = x_k;
        Vector dx = A*x_k - b;
        num_of_iter++;
        if(dx.get_norm1() <= eps)
            break;
    }

    return make_pair(x0, num_of_iter);
}



Vector gauss_method(Matrix A, Vector b){
    //перестановка строк матрицы, чтобы в i-ой строке i-ый элемент не был равен нулю.
    A.sort(b);
    //cout << A;
    int size = A.size;
    //для уобства представляем матрицу как список векторов
    vector<Vector> a;
    for(int i = 0; i < size; ++i){
        a.push_back(Vector(size));
        for(int j = 0; j < size; ++j){
            a[i].vec[j] = A.matrix[i][j];
        }
    }
    //прямой ход

    for(int i = 0;i < size;++i){
        //Разделим i-ое уравнение на aii
        b.vec[i] = b.vec[i] * (1/a[i].vec[i]);;
        a[i] = a[i] * (1/a[i].vec[i]);
        
        for(int j = i+1; j < size; j++){
            //Умножим i-ое уравнение на aji
            double bj = b.vec[i] * a[j].vec[i];
            Vector aj = a[i] * a[j].vec[i];
            //вычтем из j первое
            b.vec[j] = b.vec[j] - bj;
            a[j]= a[j] - aj;
        }
    }
    //обратный ход
    Vector x(size);
    for(int i = 0;i < size;++i)x.vec[i]=b.vec[i];
    for(int i = size -1; i >= 0;--i){
        for(int j = size-2; j >= i; j--){
            x.vec[i] = x.vec[i] - x.vec[j+1]*a[i].vec[j+1];
        }   
    }
    
    return x;
}

Vector reflection_method(Matrix A, Vector b){
    int size = A.size;

    Matrix Q(size);
    for(int i = 0; i < A.size; ++i){
        Q.matrix[i][i] = 1;
    }

    Matrix R = A;
    //Основной цикл для построения матриц R и Q
    for(int i = 0;i < size-1; ++i){

        Matrix E(size-i);
        //Еденичная матрица размера size-i
        for(int k = 0; k < size-i; ++k){
            E.matrix[k][k] = 1;
        }

        Vector y(size-i); //вектор y, (a[j][i],..,a[size][i]) j = i..size
        Vector z(size-i);//вектор z вида (1, 0, .. 0)
        z.vec[0] = 1;
        for(int j = i; j < size; ++j){
            y.vec[j-i] = R.matrix[j][i];
        }
        
        long double alpha = y.get_norm2();

        Vector p = y - z*alpha;
        long double r = p.get_norm2();

        Vector w = (y - z*alpha)*(1/r);

        Matrix W(size-i); // W = w*w^T

        for(int l = 0; l < size-i; ++l){
            for(int j = 0; j < size-i; ++j){
                W.matrix[l][j] = w.vec[l]*w.vec[j];
            }
        }

        Matrix Ri(size-i);//часть основной матрицы R
        Matrix Qi(size-i);
        for(int l = i; l < size;++l){
            for(int j = i; j < size; ++j){
                Ri.matrix[l-i][j-i] = R.matrix[l][j];
            }
        }

        Qi = E - W*2;//вычисление следующего множителя для матрицы Q.
        Ri = Qi*Ri;

        
        Matrix Q_(size);// Qi имеет размерность, отличную от Q. Поэтому Q_ - матрица размера size, содержащая в себе Qi
        for(int l = 0; l < size; ++l)
            Q_.matrix[l][l] = 1;

        for(int l = i; l < size;++l){
            for(int j = i; j < size; ++j){
                Q_.matrix[l][j] = Qi.matrix[l-i][j-i];
            }
        }

        Q = Q*Q_; // обновление матрицы Q

        for(int l = i; l < size;++l){
            for(int j = i; j < size; ++j){
                R.matrix[l][j] = Ri.matrix[l-i][j-i]; // обновление матрицы R
            }
        }
    }
    Vector y = Q.transpose()*b;

    //обратный ход
    Vector x = y;
    for(int i = size -1; i >= 0;--i){
        for(int j = size-2; j >= i; j--){
            x.vec[i] = x.vec[i] - x.vec[j+1]*R.matrix[i][j+1];
        }
        x.vec[i] /= R.matrix[i][i];   
    }

    return x;
}