//Ilnaz Magizov DSAI-02
//i.magizov@innopolis.university

#include <iostream>
#include <vector>
#include <bits/stdc++.h>
#include <random>
#define GNUPLOT_NAME "C:\\gnuplot\\bin\\gnuplot -persist"

using namespace std;

//epsilon for comparisons of numbers
double eps = 0.0000001;

//class column vector with size and vector of variables
class ColumnVector{
public:
    int size;
    vector <double> mas;

    //constructor
    ColumnVector(int size){
        this->size = size;
        mas.resize(size, 0);
    }
    //function for resize column vector
    void resize(int size){
        this->size = size;
        mas.resize(size, 0);
    }
    //function to finding norm
    double norm(){
        double ans = 0;
        for(int i = 0; i < size; ++i){
            ans += mas[i] * mas[i];
        }
        return sqrt(ans);
    }
    //overloading operator of summary
    ColumnVector operator + (ColumnVector cv){
        ColumnVector ans(this->size);
        for(int i = 0; i < this->size; ++i){
            ans.mas[i] = this->mas[i] + cv.mas[i];
        }
        return ans;
    }
    //overloading operator of equating
    ColumnVector& operator = (ColumnVector cv){
        for(int i = 0; i < this->size; ++i){
            this->mas[i] = cv.mas[i];
        }
        return *this;
    }
};

//class for matrices with number of rows, number of columns and matrix
class Matrix{
public:
    int sizeR;
    int sizeC;
    vector <vector <double>> mas;

    //constructor
    Matrix(int sizeR, int sizeC){
        this->sizeC = sizeC;
        this->sizeR = sizeR;
        mas.resize(sizeR, vector<double>(sizeC, 0));
    }
    //function for transposition matrix
    Matrix transposition(){
        Matrix ans(sizeC, sizeR);
        for(int i = 0; i < sizeR; ++i){
            for(int j = 0; j < sizeC; ++j){
                ans.mas[j][i] = mas[i][j];
            }
        }
        return ans;
    }
    //overloading operator for summation
    Matrix operator + (Matrix m2){
        Matrix ans(sizeR, sizeC);
        for(int i = 0; i < sizeR; ++i){
            for(int j = 0; j < sizeC; ++j){
                ans.mas[i][j] = mas[i][j] + m2.mas[i][j];
            }
        }
        return ans;
    }
    //overloading operator for difference
    Matrix operator - (Matrix m2){
        Matrix ans(sizeR, sizeC);
        for(int i = 0; i < sizeR; ++i){
            for(int j = 0; j < sizeC; ++j){
                ans.mas[i][j] = mas[i][j] - m2.mas[i][j];
            }
        }
        return ans;
    }
    //overloading operator for equation
    virtual Matrix& operator = (Matrix m1){
        for(int i = 0; i < sizeR; ++i){
            for(int j = 0; j < sizeC; ++j){
                mas[i][j] = m1.mas[i][j];
            }
        }
        return *this;
    }
    //overloading operator for multiplication with matrix
    Matrix operator * (Matrix m2){
        Matrix ans(sizeR, m2.sizeC);
        for(int i = 0; i < sizeR; ++i){
            for(int j = 0; j < m2.sizeC; ++j){
                for(int k = 0; k < sizeC; ++k){
                    ans.mas[i][j] += mas[i][k] * m2.mas[k][j];
                }
            }
        }
        return ans;
    }
    //overloading operator for multiplication with column vector
    ColumnVector operator * (ColumnVector m){
        ColumnVector ans(this->sizeR);
        for(int i = 0; i < this->sizeR; ++i){
            for(int j = 0; j < this->sizeC; ++j){
                ans.mas[i] += this->mas[i][j] * m.mas[j];
            }
        }
        return ans;
    }
};

//class for elimination matrix with same functions as in classic matrix
class EliminationMatrix: public Matrix{
public:
    EliminationMatrix(Matrix m, int indexI, int indexJ): Matrix(m.sizeR, m.sizeR){
        indexI--; indexJ--;
        for(int i = 0; i < m.sizeR; ++i){
            this->mas[i][i] = 1;
        }
        int newIndexJ = 0;
        for(int i = 0; i < m.sizeR; ++i){
            if(abs(m.mas[indexJ][i]) > eps){
                newIndexJ = i;
                break;
            }
        }
        this->mas[indexI][newIndexJ] = -1.0 * (m.mas[indexI][newIndexJ]) / (m.mas[indexJ][newIndexJ]);
    }

    EliminationMatrix& operator = (Matrix m) override{
        for(int i = 0; i < sizeR; ++i){
            for(int j = 0; j < sizeC; ++j){
                mas[i][j] = m.mas[i][j];
            }
        }
        return *this;
    }
};

//class for permutation matrix with same functions as in classic matrix
class PermutationMatrix: public Matrix{
public:
    PermutationMatrix(Matrix m, int indexI, int indexJ): Matrix(m.sizeR, m.sizeR){
        indexI--; indexJ--;
        this->mas[indexI][indexJ] = 1;
        this->mas[indexJ][indexI] = 1;
        for(int i = 0; i < m.sizeR; ++i){
            if(i == indexI || i == indexJ) continue;
            this->mas[i][i] = 1;
        }
    }
    PermutationMatrix& operator = (Matrix m) override{
        for(int i = 0; i < sizeR; ++i){
            for(int j = 0; j < sizeC; ++j){
                mas[i][j] = m.mas[i][j];
            }
        }
        return *this;
    }
};

//class for augmented matrix with same functions as in classic matrix
class AugmentedMatrix: public Matrix{
public:
    AugmentedMatrix(Matrix m): Matrix(m.sizeR, 2 * m.sizeC){
        for(int i = 0; i < m.sizeR; ++i){
            this->mas[i][m.sizeC + i] = 1;
        }
        for(int i = 0; i < m.sizeR; ++i){
            for(int j = 0; j < m.sizeC; ++j){
                this->mas[i][j] = m.mas[i][j];
            }
        }
    }
    AugmentedMatrix& operator = (Matrix m) override{
        for(int i = 0; i < m.sizeR; ++i){
            for(int j = 0; j < m.sizeC; ++j){
                this->mas[i][j] = m.mas[i][j];
            }
        }
        return *this;
    }
};

//overloading input of Matrix
istream& operator >> (istream& stream, Matrix& object){
    for(int i = 0; i < object.sizeR; i++){
        for(int j = 0; j < object.sizeC; j++){
            stream >> object.mas[i][j];
        }
    }
    return stream;
}

//overloading input of column vector
istream& operator >> (istream& stream, ColumnVector& object){
    for (int i = 0; i < object.size; ++i) {
        stream >> object.mas[i];
    }
    return stream;
}

//overloading output of Matrix
ostream& operator << (ostream& stream, Matrix& object){
    for(int i = 0; i < object.sizeR; i++){
        for(int j = 0; j < object.sizeC; j++){
            if(abs(object.mas[i][j]) <= eps){
                stream << 0.0 << ' ';
            }
            else {
                stream << object.mas[i][j] << ' ';
            }
        }
        stream << '\n';
    }
    return stream;
}

//overloading output of column vector
ostream& operator << (ostream& stream, ColumnVector& object){
    for (int i = 0; i < object.size; ++i) {
        if(abs(object.mas[i]) <= eps){
            stream << 0.0 << '\n';
        }
        else {
            stream << object.mas[i] << '\n';
        }
    }
    return stream;
}

//function for gaussian elimination
void GaussianElimination(Matrix& m){
    for(int i = 0; i < m.sizeR - 1; ++i){
        double maxValue = 0;
        int coordinateMaxValue = 0;
        //find maximum value in column
        for(int j = i; j < m.sizeR; ++j){
            if(abs(maxValue) < abs(m.mas[j][i]) + eps && maxValue != m.mas[j][i]){
                maxValue = m.mas[j][i];
                coordinateMaxValue = j;
            }
        }
        //check if all variables in column 0
        if(abs(maxValue) <= eps) continue;
        //check if algorithm need swap places of rows
        if(coordinateMaxValue != i) {
            PermutationMatrix pm(m, coordinateMaxValue + 1, i + 1);
            m = pm * m;
        }
        //elimination remaining lines
        for(int j = i + 1; j < m.sizeR; ++j){
            if(abs(m.mas[j][i]) <= eps) continue;
            EliminationMatrix em(m, j + 1, i + 1);
            m = em * m;
        }
    }
}

//function for zeroing out all without main diagonal
void WayBack(Matrix& m){
    for(int i = m.sizeR - 1; i >= 0; --i){
        int coordinate = 0;
        //function for find element which not equal to 0
        for(int j = m.sizeR - 1; j >= 0; --j){
            if(abs(m.mas[i][j]) > eps){
                coordinate = j;
                break;
            }
        }
        //elimination remaining lines
        for(int j = coordinate - 1; j >= 0; --j){
            if(abs(m.mas[j][i]) <= eps) continue;
            EliminationMatrix em(m, j + 1, i + 1);
            m = em * m;
        }
    }
}

//function for diagonal normalization
void DiagonalNormalization(Matrix& m){
    for(int i = 0; i < m.sizeR; ++i){
        double x = m.mas[i][i];
        if(abs(x) <= eps || (x >= 1 - eps && x <= 1 + eps)) continue;
        for(int j = i; j < m.sizeC; ++j){
            m.mas[i][j] /= x;
        }
    }
}

Matrix InverseMatrixInAugmented(Matrix m){
    Matrix ans(m.sizeR, m.sizeR);
    for(int i = 0; i < m.sizeR; ++i){
        for (int j = 0; j < m.sizeR; ++j) {
            ans.mas[i][j] = m.mas[i][m.sizeR + j];
        }
    }
    return ans;
}

Matrix inverse(const Matrix& m){
    AugmentedMatrix augmentedMatrix(m);
    GaussianElimination(augmentedMatrix);
    WayBack(augmentedMatrix);
    DiagonalNormalization(augmentedMatrix);
    return InverseMatrixInAugmented(augmentedMatrix);
}

int main() {
    FILE* pipe = _popen(GNUPLOT_NAME, "w");
    cout << fixed << setprecision(4);
    int m, n;
    m = 5;
    vector <double> mas(m);
    ColumnVector b(m);
    for(int i = 0; i < m; ++i){
        mas[i] = rand() % 50;
        b.mas[i] = rand() % 50;
        cout << mas[i] << ' ' << b.mas[i] << '\n';
    }
    n = 3;
    Matrix A(m, n + 1);
    for(int i = 0; i < m; ++i){
        for(int j = 0; j <= n; ++j){
            A.mas[i][j] = pow(mas[i], j);
        }
    }
    Matrix AT = A.transposition();
    Matrix ATbyA = AT * A;
    Matrix ATbyA_inverse = inverse(ATbyA);
    ColumnVector ATbyB = AT * b;
    ColumnVector x = ATbyA_inverse * ATbyB;
    cout << "A:\n" << A << "A_T*A:\n" << ATbyA << "(A_T*A)^-1:\n" << ATbyA_inverse << "A_T*b:\n" << ATbyB << "x~:\n" << x;
    fprintf(pipe, "plot [0 : 50] [0 : 50] %lf*x**3 + %lf*x**2 + %lf*x**1 + %lf*x**0 , '-' using 1:2 with points\n", x.mas[3], x.mas[2], x.mas[1], x.mas[0]);
    for(int i = 0; i < m; ++i){
        fprintf(pipe, "%f\t%f\n", mas[i], b.mas[i]);
    }
    fprintf(pipe, "e\n");
    fflush(pipe);
    _pclose(pipe);
    return 0;
}