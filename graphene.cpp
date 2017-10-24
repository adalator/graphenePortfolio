#include <cmath>
#include <complex>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include <Eigen/Dense>

typedef std::complex<double> complexd; //Convenient shorthand

bool isNeighbour(int boundRight, int boundBottom, int x, int y);

//Generate square lattice hamiltonian
void HGen(std::vector<std::vector<int> >& unitCell, Eigen::MatrixXcd& hamiltonian,
          const int Nx, const int Ny, const double k, const double a, const complexd t)
{
  //!!!!!!!!!!!!!! Eigen::Matrix is COLUMN MAJOR !!!!!!!!!!!!!!!!!!!!!
  //Normal C++: matrix[column][row] || Eigen: matrix[row][column]
  for (int ocnt = 0; ocnt < Ny; ++ocnt) {
    for (int icnt = 0; icnt < Nx; ++icnt) {
      //Set on-site energy for i=j elements
      hamiltonian(unitCell[icnt][ocnt]-1, unitCell[icnt][ocnt]-1) = complexd(0.0, 0.0);
      //Check if right neighbour in cell for unitCell matrix
      if (isNeighbour(Nx, Ny, icnt+1, ocnt)) {
        hamiltonian(unitCell[icnt][ocnt]-1, unitCell[icnt+1][ocnt]-1) += -t;
      }
      else {
        hamiltonian(unitCell[icnt][ocnt]-1, unitCell[0][ocnt]-1) += -t*complexd(cos(k*a), sin(k*a));
      }
      //Check if left neighbour in cell for unitCell matrix
      if (isNeighbour(Nx, Ny, icnt-1, ocnt)) {
        hamiltonian(unitCell[icnt][ocnt]-1, unitCell[icnt-1][ocnt]-1) += -t;
      }
      else {
        hamiltonian(unitCell[icnt][ocnt]-1, unitCell[Nx-1][ocnt]-1) += -t*complexd(cos(k*a), -sin(k*a));
      }
      //Check if down neighbour in cell for unitCell matrix, if not ignore
      if (isNeighbour(Nx, Ny, icnt, ocnt+1)) {
        hamiltonian(unitCell[icnt][ocnt]-1, unitCell[icnt][ocnt+1]-1) += -t;
      }
      if (isNeighbour(Nx, Ny, icnt, ocnt-1)) {
        hamiltonian(unitCell[icnt][ocnt]-1, unitCell[icnt][ocnt-1]-1) += -t;
      }
    }
  }
}

//Generate armchair hamiltonian
void HGenArmchair(Eigen::MatrixXcd& hamiltonian, const int Nx, const int Ny, const double k,
                 const double a, const complexd t)
{
  for (int ocnt = 0; ocnt < Ny; ++ocnt) {
    for (int icnt = 0; icnt < Nx; ++icnt) {
      //Sort internal hopping
      //Left exception
      if (icnt%6 == 0) {
        hamiltonian(icnt+6*ocnt, (icnt+6*ocnt)+1) += -t;
        hamiltonian(icnt+6*ocnt, (icnt+6*ocnt)+5) += -t;
      }
      //Right exception
      else if (icnt%6 == 5) {
        hamiltonian(icnt+6*ocnt, (icnt+6*ocnt)-1) += -t;
        hamiltonian(icnt+6*ocnt, (icnt+6*ocnt)-5) += -t;
      }
      //General case
      else {
        hamiltonian(icnt+6*ocnt, (icnt+6*ocnt)+1) += -t;
        hamiltonian(icnt+6*ocnt, (icnt+6*ocnt)-1) += -t;
      }
    }
    //Internal interhexagon
    if (ocnt < (Ny-1)) {
      hamiltonian(2+6*ocnt, 11+6*ocnt) += -t;
      hamiltonian(11+6*ocnt, 2+6*ocnt) += -t;
    }
  }
  //Intercell periodic boundary conditions
  hamiltonian(5, 2+6*(Ny-1)) += -t*complexd(cos(k*a), -sin(k*a));
  hamiltonian(2+6*(Ny-1), 5) += -t*complexd(cos(k*a), sin(k*a));
}

//Generate zigzag hamiltonian
void HGenZigzag(Eigen::MatrixXcd& hamiltonian, const double k, const double a,
                  const complexd t)
{
  hamiltonian(0, 1) += -t;
  hamiltonian(0, 2) += -t;
  hamiltonian(0, 1) += -t*complexd(cos(k*a), -sin(k*a));

  hamiltonian(1, 0) += -t;
  hamiltonian(1, 0) += -t*complexd(cos(k*a), sin(k*a));

  hamiltonian(2, 0) += -t;
  hamiltonian(2, 3) += -t;
  hamiltonian(2, 3) += -t*complexd(cos(k*a), -sin(k*a));

  hamiltonian(3, 2) += -t;
  hamiltonian(3, 2) += -t*complexd(cos(k*a), sin(k*a));

}


bool isNeighbour(int boundRight, int boundBottom, int x, int y)
{
  if (x < 0) {
    return false;
  }
  else if (x >= boundRight) {
    return false;
  }
  else if (y < 0) {
    return false;
  }
  else if (y >= boundBottom) {
    return false;
  }
  else {
    return true;
  }
}

//Prints a 2-D vector's components to the terminal
void printArray(const std::vector<std::vector<int> >& array, const int Nx, const int Ny)
{
  for (int outercount = 0; outercount < Ny; ++outercount) { 
    for (int innercount = 0; innercount< Nx; ++innercount) {
      std::cout << array[innercount][outercount] << " ";
    }
    std::cout << "\n";
  }
}

//Set values of hamiltonian to 0
void resetH(Eigen::MatrixXcd& hamiltonian, const int Nx, const int Ny)
{
  for (int ocnt = 0; ocnt < Nx*Ny; ++ocnt) {
    for (int icnt= 0; icnt < Nx*Ny; ++icnt) {
      hamiltonian(icnt, ocnt) = complexd(0.0, 0.0);
    }
  }
}

int main()
{
  int Nx, Ny, iterations, cellType;
  const double pi = 3.14159265359;
  std::string filename;

  std::cout << "Select cell type: 1 - Square Lattice, 2 - Armchair, 3 - Zigzag\n";
  std::cin >> cellType;
  if (cellType == 1) {
    std::cout << "Please enter a height for the unit cell\n";
    std::cin >> Ny;
    std::cout << "Please enter a width for the unit cell\n";
    std::cin >> Nx;
  }
  else if (cellType == 2) {
    Nx = 6;
    std::cout << "Please enter number of hexagons in unit cell:\n";
    std::cin >> Ny;
  }
  else if (cellType == 3) {
    Ny = 2;
    Nx = 2;
  }
  std::cout << "Please enter a number of k-iterations\n";
  std::cin >> iterations;
  std::cout << "Please enter a filename\n";
  std::cin >> filename;

  std::cout << " Ny: " << Ny << " Nx: " << Nx << " Iterations: " << iterations << "\n";
  std::cout << "Filename: " << filename << "\n";

  //Create unit cell vector initialized to 0
  std::vector<std::vector<int> > unitCell;
  unitCell.resize(Nx, std::vector<int>(Ny, 0));

  //Fill unit cell with numbers from 1 to Nx*Ny
  {
  int count = 1;
  for (int outercount = 0; outercount < Ny; ++outercount) {
    for (int innercount = 0; innercount < Nx; ++innercount) {
      unitCell[innercount][outercount] = count;
      ++count; 
    }
  }
  }

  Eigen::MatrixXcd hamiltonian(Nx*Ny, Nx*Ny);
  double k = 0.0;
  double a = 1.0;
  complexd t(2.7, 0.0);
  std::ofstream outf(filename);

  outf << "k ";
  for (int count = 1; count <= Nx*Ny; ++count) {
    outf << "E" << count << " ";
  }
  outf << "\n";

  //Generate hamiltonian for various k and output eigenenergies
  for (int count = 0; count <= iterations; ++count) {
    if (cellType == 1) {
      HGen(unitCell, hamiltonian, Nx, Ny, k, a, t);
    }
    else if (cellType == 2){
      HGenArmchair(hamiltonian, Nx, Ny, k, a, t);
    }
    else if (cellType == 3) {
      HGenZigzag(hamiltonian, k, a, t);
    }
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> solver;
    solver.compute(hamiltonian);
    //std::cout << k << "\n" << solver.eigenvalues() << "\n\n";
    outf << k << " ";
    for (int cnt = 0; cnt < Nx*Ny; ++cnt) {
      outf << solver.eigenvalues()[cnt] << " ";
    }
    outf << "\n";
    k += pi/double(iterations);
    resetH(hamiltonian, Nx, Ny);
  }

  return 0;
}
