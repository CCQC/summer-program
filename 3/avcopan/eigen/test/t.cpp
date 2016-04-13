#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>
#include <iostream>
#include <cmath>
using namespace Eigen;
using namespace std;

int main()
{
  const double pi = acos(-1.0);
  MatrixXd A(3,3);
  A << cos(1), -sin(1), 0,
       sin(1),  cos(1), 0,
           0 ,      0 , 1;

  MatrixXd B(3,3);
  B = A.pow(pi/4);
  std::cout << "The matrix A is:\n" << A << "\n\n"
               "The matrix power A^(pi/4) is:\n" << B << std::endl;
  return 0;
}
