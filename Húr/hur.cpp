#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <vector>
using namespace std;

const int L = 1;
double density;
int num_points;
double c, dt;

vector<double> stepper(vector<double> &x0, vector<double> &xm1)
{
    vector<double> x1;
    for (int i = 0; i <= num_points; i++)
    {
        x1[i] = 2 * x0[i] * xm1[i] + pow(c / (density / dt), 2) * (x0[i + 1] + x0[i - 1] - 2 * x0[i]);
    }
    x0 = xm1;

    return x1;
}

int main()
{

    cout << "Húr szimulációja \n";
    cout << "------------------- \n";
    cout << "Kezdeti feltétel, hogy melyik pontnál emeljük meg a húrt, ez egy [0;1] tartományba eső szám legyen:";
    double x0;
    cin >> x0;
    cout << "Pontok száma,terjedési sebessség:";
    cin >> num_points >> c;
    density = L / num_points;
    cout << "Lépéshossz,integrálási idő";
    double tmax;
    cin >> dt >> tmax;
    vector<double> x(num_points);
    vector<double> x_m1(num_points, 0.0);

    for (int i = 0; i <= num_points * x0; i++)
    {
        x[i] = density * i / x0;
    }

    for (int i = num_points * x0; i < num_points; i++)
    {
        x[i] = (L - density * i) / (L - x0);
    }

    return 0;
}
