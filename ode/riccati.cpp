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
/* (a = */const double initialTime( 0.12 );/* , */
const double finalTime( 0.15 );/* = b ) */
const double timeStep( 1.0e-3 ); // v techto casovych krocich ukladam do souboru pro vykresleni
const double integrationTimeStep( 1.0e-3 ); // v techto casovych krocich to resim ... presnost;
const double constParam( 1.0 );


class Error{
public:
  Error(Problem &riccati, ODESolution &u,
              const double& _initialTime,
              const double& _finalTime,
              const double& _timeStep,
              const double& _constParam )
              :problem(riccati),
              solution(u),
              initialTime(_initialTime),
              finalTime(_finalTime),
              timeStep(_timeStep),
              constParam(_constParam),
              timeStepsCount(std::ceil( std::max( 0.0, finalTime - initialTime ) / timeStep ))
              {}

  virtual double getError() = 0;

protected:
  Problem &problem;
  ODESolution &solution;
  const double initialTime;
  const double finalTime;
  const double timeStep;
  const double constParam;
  const int timeStepsCount;
};

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

template<typename Integrator>
void method(Problem &problem, ODESolution &solution, const double &integrationTimeStep){
  Integrator integrator(problem);

  integrator.setIntegrationTimeStep( integrationTimeStep );
  ODESolver< Problem, Integrator > solver( problem, integrator );

  double initialCondition[ 1 ];
  initialCondition[ 0 ] = problem.getExactSolution( initialTime );
  solver.setInitialCondition( initialCondition );

  solver.solve( solution, initialTime, finalTime, timeStep );
}

void saveErrors(const char *fileName, const double *error){
  std::fstream file;
  file.open( fileName, std::ios::out );
  for(int i = 0; i < 3; i++) file << error[i] << std::endl;
}

void solve(Problem &problem, ODESolution &solution, double timeStep, const char *prefix){
  char buff[100];
  double error[3*3];
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
  solution.write( "riccati-rungekutta.txt", initialTime, timeStep );
  error[3] = error1l(problem, solution, timeStep);
  error[4] = error2l(problem, solution, timeStep);
  error[5] = errorInf(problem, solution, timeStep);

  // Merson
  method< MersonIntegrator >(problem, solution, timeStep);
  strcpy(buff,prefix);
  strcat(buff,"riccati-merson.txt");
  solution.write( "riccati-merson.txt", initialTime, timeStep );
  error[6] = error1l(problem, solution, timeStep);
  error[7] = error2l(problem, solution, timeStep);
  error[8] = errorInf(problem, solution, timeStep);

  problem.writeExactSolution( "riccati-exact.txt", initialTime, finalTime, timeStep, constParam );

  strcpy(buff,prefix);
  strcat(buff,"riccati-errors.txt");
  saveErrors(buff,error);
}

int main( int argc, char** argv )
{
    Problem problem;
    ODESolution solution;

    solve(problem, solution, 1.0e-3, "01");
    solve(problem, solution, 1.0e-3/2.0, "02");
    solve(problem, solution, 1.0e-3/4.0, "04");
    solve(problem, solution, 1.0e-3/8.0, "08");
    solve(problem, solution, 1.0e-3/16.0, "16");
    solve(problem, solution, 1.0e-3/32.0, "32");

    return EXIT_SUCCESS;
}
