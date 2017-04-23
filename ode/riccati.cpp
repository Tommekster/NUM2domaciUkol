/*
 * File:   riccati.cpp
 * Author: oberhuber
 *
 * Created on February 25, 2016, 10:41 AM
 */

#include <cstdlib>
#include <fstream>
#include "RiccatiProblem.h"
#include "Euler.h"
#include "hw02RungeKutta.h"
#include "Merson.h"
#include "ODESolver.h"
#include "ODESolution.h"

using namespace std;

typedef RiccatiProblem Problem; // musi dat pocet stupnu volnosti a funci F (vycislit pravou stranu)
typedef Euler< Problem > EulerIntegrator; /* řešič */
typedef RungeKutta< Problem > RKIntegrator; /* řešič */
typedef Merson< Problem > MersonIntegrator; /* řešič */
// interval ve kterem to chceme resit (a,b) =
/* (a = */const double initialTime( -9.0 );/* , */
const double finalTime( -1.0 );/* = b ) */
const double timeStep( 1.0e-3 ); // v techto casovych krocich ukladam do souboru pro vykresleni
const double integrationTimeStep( 1.0e-3 ); // v techto casovych krocich to resim ... presnost;
const double constParam( 1.0 );

// (a,b) = (-9.031,0.473)
// predchozi (0.12,0.15)

double error1l(
  Problem &problem,
  ODESolution &solution,
  double &timeStep){
  const double timeStepsCount(std::ceil( std::max( 0.0, finalTime - initialTime ) / timeStep ));
  double norm = 0;

  for( int k = 0; k <= timeStepsCount; k++ ){
    double time( initialTime + k * timeStep );
    norm += std::abs(problem.getExactSolution(time) - solution.getElement(k, 0)) * timeStep;
  }

  return norm;
}

double error2l(
  Problem &problem,
  ODESolution &solution,
  double &timeStep){
  const double timeStepsCount(std::ceil( std::max( 0.0, finalTime - initialTime ) / timeStep ));
  double norm = 0;

  for( int k = 0; k <= timeStepsCount; k++ ){
    double time( initialTime + k * timeStep );
    norm += std::pow(problem.getExactSolution(time) - solution.getElement(k, 0), 2) * timeStep;
  }

  return norm;
}

double errorInf(
  Problem &problem,
  ODESolution &solution,
  double &timeStep){
  const double timeStepsCount(std::ceil( std::max( 0.0, finalTime - initialTime ) / timeStep ));
  double norm = 0;

  for( int k = 0; k <= timeStepsCount; k++ ){
    double time( initialTime + k * timeStep );
    double diff = std::abs(problem.getExactSolution(time) - solution.getElement(k, 0)) * timeStep;
    if(diff > norm) norm = diff;
  }

  return norm;
}

double eol(double error1, double tau1, double error2, double tau2){
  return std::log(error1/error2)/std::log(tau1/tau2);
}

void countEOLs(const double *error, double frac, double *eols){
  for(int i = 0; i < 9; i++) eols[i] = eol(error[i],error[i+9],frac,1.0);
}

template<typename Integrator>
void method(Problem &problem, ODESolution &solution, const double &timeStep){
  Integrator integrator(problem);

  integrator.setIntegrationTimeStep( timeStep );
  ODESolver< Problem, Integrator > solver( problem, integrator );

  double initialCondition[ 1 ];
  initialCondition[ 0 ] = problem.getExactSolution( initialTime );
  solver.setInitialCondition( initialCondition );

  solver.solve( solution, initialTime, finalTime, timeStep );
}

void saveErrors(const char *fileName, const double *error, int len, const char **labels = 0){
  std::fstream file;
  file.open( fileName, std::ios::out );
  if(labels == 0) for(int i = 0; i < len; i++) file << error[i] << std::endl;
  else for(int i = 0; i < len; i++) file << labels[i] << error[i] << std::endl;
}

void solve(Problem &problem, ODESolution &solution, double *error, double timeStep, const char *prefix){
  char buff[100];
  //double error[3*3];
  // Euler
  method< EulerIntegrator >(problem, solution, timeStep);
  strcpy(buff,prefix);
  strcat(buff,"riccati-euler.txt");
  solution.write( buff, initialTime, timeStep );
  error[0] = error1l(problem, solution, timeStep);
  error[1] = error2l(problem, solution, timeStep);
  error[2] = errorInf(problem, solution, timeStep);

  // Runge-Kutta
  method< RKIntegrator >(problem, solution, timeStep);
  strcpy(buff,prefix);
  strcat(buff,"riccati-rungekutta.txt");
  solution.write( buff, initialTime, timeStep );
  error[3] = error1l(problem, solution, timeStep);
  error[4] = error2l(problem, solution, timeStep);
  error[5] = errorInf(problem, solution, timeStep);

  // Merson
  method< MersonIntegrator >(problem, solution, timeStep);
  strcpy(buff,prefix);
  strcat(buff,"riccati-merson.txt");
  solution.write( buff, initialTime, timeStep );
  error[6] = error1l(problem, solution, timeStep);
  error[7] = error2l(problem, solution, timeStep);
  error[8] = errorInf(problem, solution, timeStep);

  strcpy(buff,prefix);
  strcat(buff,"riccati-exact.txt");
  problem.writeExactSolution( buff, initialTime, finalTime, timeStep, constParam );

  strcpy(buff,prefix);
  strcat(buff,"riccati-errors.txt");
  saveErrors(buff,error,9);
}

int main( int argc, char** argv )
{
    Problem problem;
    ODESolution solution;

    double error[3*3*6];
    double eols[3*3*5];
    solve(problem, solution, &error[0], 1.0e-3, "01/");
    solve(problem, solution, &error[9], 1.0e-3/2.0, "02/");
    countEOLs(&error[0],2.0,&eols[0]);
    solve(problem, solution, &error[18], 1.0e-3/4.0, "04/");
    countEOLs(&error[9],2.0,&eols[9]);
    solve(problem, solution, &error[27], 1.0e-3/8.0, "08/");
    countEOLs(&error[18],2.0,&eols[18]);
    solve(problem, solution, &error[36], 1.0e-3/16.0, "16/");
    countEOLs(&error[27],2.0,&eols[27]);
    solve(problem, solution, &error[45], 1.0e-3/32.0, "32/");
    countEOLs(&error[36],2.0,&eols[36]);

    //problem.writeExactSolution("riccati-exact.txt", initialTime, finalTime,
    //  1.0e-3/32.0, constParam);

    const char *labels[] = {
      "timeStep: 1e-3\nEuler\nE_1l = ",
      "E_2l = ",
      "E_In = ",
      "Runge-Kutta\nE_1l = ",
      "E_2l = ",
      "E_In = ",
      "Merson\nE_1l = ",
      "E_2l = ",
      "E_In = ",

      "\ntimeStep: 1e-3/2\nEuler\nE_1l = ",
      "E_2l = ",
      "E_In = ",
      "Runge-Kutta\nE_1l = ",
      "E_2l = ",
      "E_In = ",
      "Merson\nE_1l = ",
      "E_2l = ",
      "E_In = ",

      "\ntimeStep: 1e-3/4\nEuler\nE_1l = ",
      "E_2l = ",
      "E_In = ",
      "Runge-Kutta\nE_1l = ",
      "E_2l = ",
      "E_In = ",
      "Merson\nE_1l = ",
      "E_2l = ",
      "E_In = ",

      "\ntimeStep: 1e-3/8\nEuler\nE_1l = ",
      "E_2l = ",
      "E_In = ",
      "Runge-Kutta\nE_1l = ",
      "E_2l = ",
      "E_In = ",
      "Merson\nE_1l = ",
      "E_2l = ",
      "E_In = ",

      "\ntimeStep: 1e-3/16\nEuler\nE_1l = ",
      "E_2l = ",
      "E_In = ",
      "Runge-Kutta\nE_1l = ",
      "E_2l = ",
      "E_In = ",
      "Merson\nE_1l = ",
      "E_2l = ",
      "E_In = ",

      "\ntimeStep: 1e-3/32\nEuler\nE_1l = ",
      "E_2l = ",
      "E_In = ",
      "Runge-Kutta\nE_1l = ",
      "E_2l = ",
      "E_In = ",
      "Merson\nE_1l = ",
      "E_2l = ",
      "E_In = ",
    };
    saveErrors("riccati-errors.txt",error,3*3*6,labels);

    const char *labels2[] = {
      "timeStep: 1e-3 -> 1e-3/2\nEuler\nE_1l = ",
      "E_2l = ",
      "E_In = ",
      "Runge-Kutta\nE_1l = ",
      "E_2l = ",
      "E_In = ",
      "Merson\nE_1l = ",
      "E_2l = ",
      "E_In = ",

      "\ntimeStep: 1e-3/2 -> 1e-3/4\nEuler\nE_1l = ",
      "E_2l = ",
      "E_In = ",
      "Runge-Kutta\nE_1l = ",
      "E_2l = ",
      "E_In = ",
      "Merson\nE_1l = ",
      "E_2l = ",
      "E_In = ",

      "\ntimeStep: 1e-3/4 -> 1e-3/8\nEuler\nE_1l = ",
      "E_2l = ",
      "E_In = ",
      "Runge-Kutta\nE_1l = ",
      "E_2l = ",
      "E_In = ",
      "Merson\nE_1l = ",
      "E_2l = ",
      "E_In = ",

      "\ntimeStep: 1e-3/8 -> 1e-3/16\nEuler\nE_1l = ",
      "E_2l = ",
      "E_In = ",
      "Runge-Kutta\nE_1l = ",
      "E_2l = ",
      "E_In = ",
      "Merson\nE_1l = ",
      "E_2l = ",
      "E_In = ",

      "\ntimeStep: 1e-3/16 -> 1e-3/32\nEuler\nE_1l = ",
      "E_2l = ",
      "E_In = ",
      "Runge-Kutta\nE_1l = ",
      "E_2l = ",
      "E_In = ",
      "Merson\nE_1l = ",
      "E_2l = ",
      "E_In = ",
    };
    saveErrors("riccati-eol.txt",error,3*3*5,labels2);

    return EXIT_SUCCESS;
}
