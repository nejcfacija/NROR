#define _USE_MATH_DEFINES
#include <vector>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>
#include <chrono>
#include <omp.h>
#include <ctime>
#include <unistd.h>

using namespace std;

int main() {

    vector<vector<double>> A;
  
    vector<double> b;

    std::string filename = "./datoteka_A_b.txt";

    std::ifstream infile;
    infile.open(filename);

    std::string string_first_line;
    std::getline(infile, string_first_line);

    std::replace(string_first_line.begin(), string_first_line.end(), '=', ' ');

    std::istringstream iss(string_first_line);
    std::string nepomemben_del1;
    std::string nepomemben_del2;
    int n;

    iss >> nepomemben_del1;
    iss >> nepomemben_del2;
    iss >> n;

    std::cout << "Velikost matrike A: " << n << "x" << n << std::endl;;

    for (int iiA = 0; iiA < n; iiA++)
    {
        std::string line;
        std::getline(infile, line);
        std::replace(line.begin(), line.end(), ';', ' ');

        std::istringstream iss_column(line);

        vector<double> row;

        for (int column = 0; column < n; column++)
        {
            double element_a = 0;
            iss_column >> element_a;
            row.push_back(element_a);
        }

        A.push_back(row);
    }

    std::string empty_line;
    std::getline(infile, empty_line);

    std::string string_line_b;
    std::getline(infile, string_line_b);

    std::replace(string_line_b.begin(), string_line_b.end(), '>', ' ');
    std::istringstream iss_b(string_line_b);
    int n_b;

    iss_b >> nepomemben_del1;
    iss_b >> nepomemben_del2;
    iss_b >> n_b;

    std::cout << "Velikost vektorja b: " << n_b << std::endl;;

    for (int iib = 0; iib < n_b; iib++)
    {
        std::string line_b_element;
        std::getline(infile, line_b_element);
        std::istringstream iss_b_element(line_b_element);

        double b_element = 0;
        iss_b_element >> b_element;

        b.push_back(b_element);
    }

    vector<double> T;
    for (int iiT = 0; iiT < n_b; iiT++)
    {
        T.push_back(100);
    }

    auto start_time = std::chrono::high_resolution_clock::now();
    for (int ii = 0; ii < 2000; ii++)
    {
        vector<double> T_old = T;

        for (int i = 0; i < n; i++) {
            double sum = 0.0;

            for (int j = 0; j < n; j++) {
                if (i != j) {
                    sum += A[i][j] * T_old[j];
                }
            }

            T[i] = (b[i] - sum) / A[i][i];
        }

    }

    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> time_duration = end_time - start_time;
    std::cout << "Time of Gauss-Seidel: " << time_duration.count() << " seconds" << std::endl;

    double max_T = 0;
    for (int iiT = 0; iiT < n_b; iiT++)
    {
        cout << T[iiT] << endl;

        if (T[iiT] > max_T) {
            max_T = T[iiT];
        }
    }
    std::cout << "Serial Time: " << time_duration.count() << " seconds" << std::endl;
    cout << "Max. temperature: " << max_T << " degree C." << endl;



    return 0;
}
