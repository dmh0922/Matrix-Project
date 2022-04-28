#include <iostream>
#include "Matrix.h"
#define EPSILON 1e-7

Matrix::Matrix(const int N, const int M) : LEN_Y(N), LEN_X(M)
{
    this->Content = new double* [LEN_Y];
    for (int y=0; y<LEN_Y; y++) { this->Content[y] = new double [LEN_X]; }
}

Matrix::Matrix(double** Data, const int N, const int M) : LEN_Y(N), LEN_X(M)
{
    this->Content = new double* [LEN_Y];
    for (int y=0; y<LEN_Y; y++) { this->Content[y] = new double [LEN_X]; }

    for (int y=0; y<LEN_Y; y++) {
        for (int x=0; x<LEN_X; x++) {
            Content[y][x] = Data[y][x];
        }
    }
}

Matrix::Matrix(const Matrix& Origin) : LEN_Y(Origin.LEN_Y) , LEN_X(Origin.LEN_X)
{
    this->Content = new double* [LEN_Y];
    for (int y=0; y<LEN_Y; y++) { this->Content[y] = new double [LEN_X]; }

    for (int y=0; y<LEN_Y; y++) {
        for (int x=0; x<LEN_X; x++) {
            Content[y][x] = Origin.Content[y][x];
        }
    }
}

Matrix::~Matrix() { for (int y=0; y<LEN_Y; y++) { delete[] this->Content[y]; } }

// ----------------------------------

// Can be used to change its Content
void Matrix::set(double** Data)
{
    for (int y=0; y<this->LEN_Y; y++) {
        for (int x=0; x<this->LEN_X; x++) {
            this->Content[y][x] = Data[y][x];
        }
    }
}

// Multiplies itself with the given SCALAR
void Matrix::mult(const double SCALAR)
{
    for (int y=0; y<LEN_Y; y++) {
        for (int x=0; x<LEN_X; x++) { this->Content[y][x] *= SCALAR; }
    }
    return;
}

// Shows its Content
void Matrix::Show(void) const
{
    std::cout << "------------------\n";
    for (int y=0; y<LEN_Y; y++) {
        for (int x=0; x<LEN_X; x++) {
            std::cout << Content[y][x] << ' ';
        }
        std::cout << '\n';
    }
    std::cout << "------------------\n";
    return;   
}

// Returns its size : [ LEN_Y , LEN_X ]
const std::pair<int,int> Matrix::getSize(void) const { return std::pair(this->LEN_X, this->LEN_Y); }

// Returns the number at Content[ _y ][ _x ]
const double Matrix::getValue(const int _y, const int _x) const { return this->Content[_y-1][_x-1]; }

// --------------------------------

// Checking wether this Matrix is same to the Origin Matrix or not
const bool Matrix::equals(const Matrix& Origin) const
{
    if ( this->getSize() != Origin.getSize() ) return false;

    for (int y=0; y<LEN_Y; y++) {
        for (int x=0; x<LEN_X; x++) {
            if ((std::abs(this->Content[y][x] - Origin.Content[y][x])) > EPSILON )
                return false;
        }
    }
    return true;
}

const bool Matrix::isSquare(void) const { return this->LEN_X == this->LEN_Y; }
const bool Matrix::isSymmetric(void) const { return this->equals(this->transpose()); }
const bool Matrix::isOrthogonal(void) const { return this->equals(this->transpose().inverse()); }

// --------------------------------

// Returns the determinant of this Matrix
const double Matrix::det(void) const
{
    if ( !this->isSquare() ) return 0;
    if (this->LEN_X == 1) { return Content[0][0]; }
    if (this->LEN_X == 2) { return Content[0][0]*Content[1][1] - Content[0][1]*Content[1][0] ; }

    double ret = 0;
    for (int x=0; x<LEN_X; x++) { ret += cofactor(x,0) * this->Content[0][x]; }
    return ret;
};

// Returns the cofactor of the number at : Content[ _y ][ _x ]
const double Matrix::cofactor(const int _y, const int _x) const
{
    if ( !this->isSquare() ) return 0;

    int N = this->LEN_X - 1;
    int M = this->LEN_Y - 1;

    int idx_x = 0, idx_y = 0;
    double** Data = new double* [M];
    for (int iter=0; iter<M; iter++) { Data[iter] = new double[N]; }

    for (int y=0; y<LEN_Y; y++) {
        if (y == _y) continue;
        for (int x=0; x<LEN_X; x++) {
            if (x == _x) continue;
            Data[idx_y][idx_x++] = this->Content[y][x];
        }
        idx_y++;
        idx_x = 0;
    }

    Matrix* temp = new Matrix(Data,N, M);
    double ret = temp->det();
    for (int iter=0; iter<M; iter++) { delete[] Data[iter]; }
    delete temp;

    return ( (_x + _y) % 2 == 0 ) ? ret : -ret;
}

// --------------------------------

// Returns its transpose Matrix
const Matrix Matrix::transpose(void) const
{
    double** Data = new double* [LEN_X];
    for (int iter=0; iter<LEN_X; iter++) { Data[iter] = new double [LEN_Y]; }

    for (int x=0; x<LEN_X; x++) {
        for (int y=0; y<LEN_Y; y++) {
            Data[y][x] = this->Content[x][y];
        }
    }

    Matrix ret = Matrix(Data,LEN_Y,LEN_X);
    for (int iter=0; iter<LEN_X; iter++) { delete[] Data[iter]; }

    return ret;
}

// Returns the Matrix with its cofactors
const Matrix Matrix::cofact_matrix(void) const
{
    double** Data = new double* [LEN_Y];
    for (int iter=0; iter<LEN_X; iter++) { Data[iter] = new double [LEN_X]; }

    for (int y=0; y<LEN_Y; y++) {
        for (int x=0; x<LEN_X; x++) {
            Data[y][x] = this->cofactor(x,y);
        }
    }

    Matrix ret = Matrix(Data,LEN_X,LEN_Y);
    for (int iter=0; iter<LEN_X; iter++) { delete[] Data[iter]; }

    return ret;
}

// Returns its Inverse Matrix
const Matrix Matrix::inverse(void) const
{
    const double& DET = this->det();
    if (DET == 0) { std::cout << "This Matrix Can't Have an inverse : det = 0.\n"; return Matrix(); }
    if ( !this->isSquare() ) { std::cout << "This Matrix Can't Have an inverse : Not a Square.\n"; return Matrix(); }

    Matrix ret = this->cofact_matrix().transpose();
    ret.mult(1/DET);
    return ret;
}

//Returns the product with Matrix Origin
const Matrix Matrix::product(const Matrix& Origin) const
{
    if (this->LEN_X != Origin.LEN_Y) return Matrix();

    double** Data = new double* [this->LEN_Y];
    for (int y=0; y<this->LEN_Y; y++) { Data[y] = new double [Origin.LEN_X]; }

    for (int y=0; y<this->LEN_Y; y++) {
        for (int x=0; x<Origin.LEN_X; x++) {
            double value = 0;
            for (int iter=0; iter<this->LEN_X; iter++) {
                value += this->Content[y][iter] * Origin.Content[iter][x];
            }
            Data[y][x] = value;
        }
    }

    Matrix ret = Matrix(Data,this->LEN_Y, Origin.LEN_X);
    for (int iter=0; iter<this->LEN_Y; iter++) { delete[] Data[iter]; }
    return ret;
}
