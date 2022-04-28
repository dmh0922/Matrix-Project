
class Matrix
{
private:
    double** Content;
    const int LEN_Y;
    const int LEN_X;

public:
    Matrix() : LEN_Y(-1) , LEN_X(-1) {};
    Matrix(const int N, const int M);
    Matrix(double** Data, const int N, const int M);
    Matrix(const Matrix& Origin);
    ~Matrix();

    void set(double** Data);
    void mult(const double SCALAR);

    void Show(void) const;
    const std::pair<int,int> getSize(void) const;
    const double getValue(const int _x, const int _y) const;

    const bool equals(const Matrix& Origin) const;
    const bool isSquare(void) const;
    const bool isSymmetric(void) const;
    const bool isOrthogonal(void) const;

    const double det(void) const;
    const double cofactor(const int _x, const int _y) const;

    const Matrix transpose(void) const;
    const Matrix cofact_matrix(void) const;
    const Matrix inverse(void) const;
    const Matrix product(const Matrix& M) const;
};