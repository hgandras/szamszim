#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <vector>
using namespace std;

#include "vector.hpp"
#include "odeint.hpp"
using namespace cpl;

const double pi = 4 * atan(1.0);
const double G = 1;
double m1, m2, m3;
bool switch_t_with_y = false; //  to interpolate to y = 0

double Cubed(std::vector<double> r1, std::vector<double> r2)
{
    double r_squared = pow(r1[0] - r2[0], 2) + pow(r1[1] - r2[1], 2);
    double r_cubed = r_squared * sqrt(r_squared);
    return r_cubed;
}

//  Derivative vector for Newton's law of gravitation
Vector f(const Vector &x)
{
    double t = x[0];
    std::vector<double> r1 = {x[1], x[2]};
    std::vector<double> r2 = {x[3], x[4]};
    std::vector<double> r3 = {x[5], x[6]};

    double r12_cubed = Cubed(r1, r2);
    double r13_cubed = Cubed(r1, r3);
    double r23_cubed = Cubed(r2, r3);

    Vector f(13);
    f[0] = 1;
    for (int i = 1; i < 7; i++)
    {
        f[i] = x[i + 6];
    }

    for (int i = 7; i < 13; i++) // Sebességek léptetése
    {
        if (i == 7 or i == 8)
        {
            f[i] = -G * m2 * (x[i - 6] - x[i + 2 - 6]) / r12_cubed - G * m3 * (x[i - 6] - x[i + 4 - 6]) / r13_cubed;
        }
        if (i == 9 or i == 10)
        {
            f[i] = -G * m1 * (x[i - 6] - x[i - 2 - 6]) / r12_cubed - G * m3 * (x[i - 6] - x[i + 2 - 6]) / r23_cubed;
        }
        if (i == 11 or i == 12)
        {
            f[i] = -G * m1 * (x[i - 6] - x[i - 4 - 6]) / r13_cubed - G * m2 * (x[i - 6] - x[i - 2 - 6]) / r23_cubed;
        }
    }

    return f;
}

//  Change independent variable from t to y and step back to y = 0
void interpolate_crossing(Vector x, int &crossing)
{
    ++crossing;
    switch_t_with_y = true;
    RK4Step(x, -x[2], f);
    cout << " crossing " << crossing << "\t t = " << x[0]
         << "\t x = " << x[1] << endl;
    switch_t_with_y = false;
}

int main()
{
    cout << " 3 body problem\n"
         << " -----------------------------------------------------\n"
         << " Enter initial positions (r1_x,r1_y,r2_x,r2_y,r3_x,r3_y):";
    double a1;
    Vector x0; // Vektor ami tartalmazza a pozíciókat és a sebességeket
    x0[0] = 0; // t_0 kezdőpillanat, ezt 0-ra rakom
    for (int i = 0; i < 6; i++)
    {
        cin >> a1;
        x0.push_back(a1);
    }
    cout << "Enter initial velocities (v1_x,v1_y,v2_x,v2_y,v3_x,v3_y):";
    for (int i = 0; i < 6; i++)
    {
        cin >> a1;
        x0.push_back(a1);
    }

    cout << " Enter masses(m1,m2,m3):";

    cin >> m1 >> m2 >> m3;

    cout << " Enter step size, integrating time, and adaptive accuracy: ";
    double dt, t_max, accuracy;
    cin >> dt >> t_max >> accuracy;

    ofstream dataFile("3body_fixed.data");
    Vector x = x0;
    cout << "\n Integrating with fixed step size" << endl;
    do
    {
        for (int i = 0; i < 13; i++)
            dataFile << x[i] << '\t';
        dataFile << '\n';
        RK4Step(x, dt, f);
    } while (x[0] < t_max);
    cout << " data in file fixed.data" << endl;
    dataFile.close();

    dataFile.open("3body_adaptive.data");
    x = x0;
    double dt_max = 0, dt_min = 100;
    cout << "\n Integrating with adaptive step size" << endl;
    do
    {
        for (int i = 0; i < 13; i++)
            dataFile << x[i] << '\t';
        dataFile << '\n';
        adaptiveRK4Step(x, dt, accuracy, f);
    } while (x[0] < t_max);
    cout << " step size: min = " << dt_min << "  max = " << dt_max << endl;
    cout << " data in file adaptive.data" << endl;
    dataFile.close();
}
