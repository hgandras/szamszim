#include <cmath>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
using namespace std;

#include "vector.hpp" // vectors with components of type double
#include "odeint.hpp" // ODE integration routines, Runge-Kutta ...
using namespace cpl;

const double pi = 4 * atan(1.0);

const double g = 9.8; // acceleration of gravity

double L1 = 2.0; // length of pendulum
double L2 = 1.0;
double m1 = 1.0;
double m2 = 1.0;
double q = 0.5;             // damping coefficient
double Omega_D = 2.0 / 3.0; // frequency of driving force
double F_D = 0.9;           // amplitude of driving force
bool nonlinear;             // linear if false

// Euler stepper
void Euler(Vector &x, double tau, Vector derivs(const Vector &))
{
    Vector v = tau * derivs(x);
    x += v;
}

// Euler-Cromer stepper
void EulerCromer(Vector &x, double tau, Vector derivs(const Vector &))
{
    Vector v = tau * derivs(x);
    double x_t = x[1];
    x += v;
    x[1] = x_t + x[2] * tau;
}

Vector f(const Vector &x)
{ // extended derivative vector
    double t = x[0];
    double theta = x[1];
    double omega = x[2];
    Vector f(3); // Vector with 3 components
    f[0] = 1;
    f[1] = omega;
    if (nonlinear)
        f[2] = -(g / L1) * sin(theta) - q * omega + F_D * sin(Omega_D * t);
    else
        f[2] = -(g / L1) * theta - q * omega + F_D * sin(Omega_D * t);
    return f;
}

Vector dp(const Vector &x)
{

    double t = x[0];
    double theta_1 = x[1];
    double omega_1 = x[2];
    double theta_2 = x[3];
    double omega_2 = x[4];

    double c = cos(theta_1 - theta_2);
    double s = sin(theta_1 - theta_2);

    Vector f(5);
    f[0] = 1;
    f[1] = omega_1;
    f[3] = omega_2;

    f[2] = (m2 * g * sin(theta_2) * c - m2 * s * (L1 * pow(omega_1, 2) * c + L2 * pow(omega_2, 2)) -
            (m1 + m2) * g * sin(theta_1)) /
           L1 / (m1 + m2 * pow(s, 2));

    f[4] = ((m1 + m2) * (L1 * pow(omega_1, 2) * s - g * sin(theta_2) + g * sin(theta_1) * c) +
            m2 * L2 * pow(omega_2, 2) * s * c) /
           L2 / (m1 + m2 * pow(s, 2));

    return f;
}

int main()
{
    cout << " Nonlinear damped driven pendulum\n"
         << " --------------------------------\n";
    cout << " Enter theta_1(0),omega_1(0),theta_2(0),omega_2(0):";
    double theta_1, omega_1, theta_2, omega_2, tMax;
    cin >> theta_1 >> omega_1 >> theta_2 >> omega_2;
    cout << " Enter integration time t_max: ";
    cin >> tMax;

    double dt = 0.05;
    double accuracy = 1e-6;
    ofstream dataFile("pendulum.data");

    double t = 0;
    Vector x(5);
    x[0] = t;
    x[1] = theta_1;
    x[2] = omega_1;
    x[3] = theta_2;
    x[4] = omega_2;

    while (t < tMax)
    {
        RK4Step(x, dt, dp);
        t = x[0], theta_1 = x[1], omega_1 = x[2], theta_2 = x[3], omega_2 = x[4];

        dataFile << t << '\t' << theta_1 << '\t' << omega_1 << '\t' << theta_2 << '\t' << omega_2 << '\n';
    }

    cout << " Output data to file pendulum.data" << endl;
    dataFile.close();
}
