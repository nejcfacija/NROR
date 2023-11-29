#include <iostream>
#include <fstream>
#include <vector>

double naprejRazlika(std::vector<std::pair<double, double>>& podatki, int indeks, double dx) {
    return (-3 * podatki[indeks].second + 4 * podatki[indeks + 1].second - podatki[indeks + 2].second) / (2 * dx);
}

double nazajRazlika(std::vector<std::pair<double, double>>& podatki, int indeks, double dx) {
    return (3 * podatki[indeks].second - 4 * podatki[indeks - 1].second + podatki[indeks - 2].second) / (2 * dx);
}

double centralnaRazlika(std::vector<std::pair<double, double>>& podatki, int indeks, double dx) {
    return (podatki[indeks + 1].second - podatki[indeks - 1].second) / (2 * dx);
}

int main() {
    std::ifstream vhodnaDatoteka("poly.txt");
    if (!vhodnaDatoteka.is_open()) {
        std::cerr << "Ni mogoče odpreti datoteke poly.txt\n";
        return 1;
    }

    int steviloVrednosti;
    vhodnaDatoteka >> steviloVrednosti;

    std::vector<std::pair<double, double>> podatki;
    double x, fx;
    for (int i = 0; i < steviloVrednosti; ++i) {
        vhodnaDatoteka >> x >> fx;
        podatki.emplace_back(x, fx);
    }

    vhodnaDatoteka.close();

    double dx = 0.0294117647058822; // ∆x
    std::ofstream izhodnaDatoteka("diff_poly.txt");
    if (!izhodnaDatoteka.is_open()) {
        std::cerr << "Ni mogoče ustvariti datoteke diff_poly.txt\n";
        return 1;
    }

    for (int i = 0; i < steviloVrednosti; ++i) {
        double odvod;
        if (i == 0) {
            odvod = naprejRazlika(podatki, i, dx);
        } else if (i == steviloVrednosti - 1) {
            odvod = nazajRazlika(podatki, i, dx);
        } else {
            odvod = centralnaRazlika(podatki, i, dx);
        }
        izhodnaDatoteka << "Odvod f'(" << podatki[i].first << ") = " << odvod << std::endl;
    }

    izhodnaDatoteka.close();

    std::cout << "Diferencialni rezultati shranjeni v diff_poly.txt\n";

    return 0;
}
