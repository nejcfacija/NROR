#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <algorithm>
#include <chrono>
#include <nlohmann/json.hpp>
using json = nlohmann::json;




int main() {
    auto start = std::chrono::high_resolution_clock::now();

    // Potrebovali bomo sledeče parametre
    std::vector<double> X;
    std::vector<double> Y;
    std::vector <std::vector<int>> celice;
    std::vector<std::string> tipi_robnih_pogojev_str;
    std::vector<int> tipi_robnih_pogojev;
    std::vector<double> vrednosti_robnih_pogojev;
    std::vector<double> vrednosti_prestopa_toplote;
    std::vector<std::vector<int>> vozlisca_robnih_pogojev;

    // Definiramo pot do datoteke z mrežo in branje datoteke
    std::string filename = "primer4mreza.txt";

    std::ifstream file;

    file.open(filename);

    std::string str_tocke;
    file >> str_tocke;
    std::string str_points;
    file >> str_points;
    int st_neznank = std::stoi(str_points);

    std::string skrita_vrstica1;
    std::getline(file, skrita_vrstica1);

    for (size_t i = 0; i < st_neznank; ++i) {

        std::string s;
        std::getline(file, s);
        std::replace(s.begin(), s.end(), ';', ' ');
        std::replace(s.begin(), s.end(), ',', ' ');
        std::istringstream iss(s);
        int node_id;
        double x;
        double y;
        iss >> node_id >> x >> y;

        X.push_back(x);
        Y.push_back(y);

    }

    std::string prazna_vrstica1;
    std::getline(file, prazna_vrstica1);

    std::string str_celice;
    file >> str_celice;
    std::string str_cells;
    file >> str_cells;
    int st_celic = std::stoi(str_cells);

    std::string skrita_vrstica2;
    std::getline(file, skrita_vrstica2);

    for (size_t i = 0; i < st_celic; ++i) {

        std::string s;
        std::getline(file, s);
        std::replace(s.begin(), s.end(), ';', ' ');
        std::replace(s.begin(), s.end(), ',', ' ');
        std::istringstream iss(s);
        std::vector<int> cell;
        int cell_id;
        int node0_id;
        int node1_id;
        int node2_id;
        int node3_id;
        iss >> cell_id >> node0_id >> node1_id >> node2_id >> node3_id;

        cell = { node0_id, node1_id, node2_id, node3_id };
        celice.push_back(cell);

    }

    std::string prazna_vrstica2;
    std::getline(file, prazna_vrstica2);



    std::string str_robni;
    file >> str_robni;
    std::string str_pogoji;
    file >> str_pogoji;
    std::string str_rp;
    file >> str_rp;
    int st_pogojev = std::stoi(str_rp);


    std::string skrita_vrstica3;
    std::getline(file, skrita_vrstica3);


    for (size_t i = 0; i < st_pogojev; ++i) {

        std::string s;
        std::getline(file, s);
        std::istringstream iss(s);
        std::string pogoj;
        std::string ena;
        std::string tip_pogoja;

        iss >> pogoj >> ena >> tip_pogoja;
        tipi_robnih_pogojev_str.push_back(tip_pogoja);

        if (tip_pogoja == "temperatura") {

            tipi_robnih_pogojev.push_back(0);

            std::string s;
            std::getline(file, s);
            std::istringstream iss(s);
            std::string temperatura;
            double vrednost_temperature;
            iss >> temperatura >> vrednost_temperature;

            vrednosti_robnih_pogojev.push_back(vrednost_temperature);
            vrednosti_prestopa_toplote.push_back(-1.);

        }

        else if (tip_pogoja == "toplotni") {

            tipi_robnih_pogojev.push_back(1);

            std::string s;
            std::getline(file, s);
            std::istringstream iss(s);
            std::string toplotni;
            std::string tok;
            double toplotni_tok;
            iss >> toplotni >> tok >> toplotni_tok;

            vrednosti_robnih_pogojev.push_back(toplotni_tok);
            vrednosti_prestopa_toplote.push_back(-1.);

        }

        else if (tip_pogoja == "prestop") {

            tipi_robnih_pogojev.push_back(2);

            std::string s;
            std::getline(file, s);
            std::istringstream iss1(s);
            std::string temperatura;
            double vrednost_temperature;
            iss1 >> temperatura >> vrednost_temperature;

            std::string t;
            std::getline(file, t);
            std::istringstream iss2(t);
            std::string dum1;
            std::string dum2;
            double vrednost_prestopa;
            iss2 >> dum1 >> dum2 >> vrednost_prestopa;

            vrednosti_robnih_pogojev.push_back(vrednost_temperature);
            vrednosti_prestopa_toplote.push_back(vrednost_prestopa);

        }

        std::getline(file, s);
        std::string st_vozl_v_rp_str = s;
        int st_vozl_v_rp = std::stoi(st_vozl_v_rp_str);
        std::vector<int> vozl_v_rp;

        for (size_t j = 0; j < st_vozl_v_rp; ++j) {

            std::string s;
            std::getline(file, s);
            std::istringstream iss(s);
            int id_vozlisca;
            iss >> id_vozlisca;
            vozl_v_rp.push_back(id_vozlisca);

        }

        vozlisca_robnih_pogojev.push_back(vozl_v_rp);
        std::string prazna_vrstica3;
        std::getline(file, prazna_vrstica3);

    }

    file.close();


    // Sosednja vozlisca
    std::vector<std::vector<int>> sosednja_vozlisca;


    for (size_t i = 0; i < st_neznank; ++i) {

        std::vector<int> sosedje_i = { -1, -1, -1, -1 };

        for (size_t j = 0; j < st_celic; ++j) {

            std::vector<int> trenutna_celica = celice[j];
            int vozlisce0 = trenutna_celica[0];
            int vozlisce1 = trenutna_celica[1];
            int vozlisce2 = trenutna_celica[2];
            int vozlisce3 = trenutna_celica[3];

            if (i == vozlisce0 || i == vozlisce1 || i == vozlisce2 || i == vozlisce3) {

                for (size_t k = 0; k < 4; ++k) {

                    int sosednje_vozlisce = trenutna_celica[k];
                    if (sosednje_vozlisce != i) {
                        int pozicija = -1;
                        int x_obravnavano_vozl = X[i];
                        int y_obravnavano_vozl = Y[i];
                        int x_sosed = X[sosednje_vozlisce];
                        int y_sosed = Y[sosednje_vozlisce];
                        if (x_obravnavano_vozl - x_sosed < 0.000000001 && x_obravnavano_vozl - x_sosed > -0.000000001) {

                            if (y_obravnavano_vozl - y_sosed > 0) {
                                pozicija = 1;
                            }
                            else {
                                pozicija = 3;
                            }

                        }
                        else if (y_obravnavano_vozl - y_sosed < 0.000000001 && y_obravnavano_vozl - y_sosed > -0.000000001) {

                            if (x_obravnavano_vozl - x_sosed > 0) {
                                pozicija = 0;
                            }
                            else {
                                pozicija = 2;
                            }

                        }
                        else {
                            pozicija = -1;
                        }
                        if (pozicija != -1) {
                            sosedje_i[pozicija] = sosednje_vozlisce;
                        }

                    }

                }

            }

        }
        sosednja_vozlisca.push_back(sosedje_i);

    }

    // Gradnja matrike
    std::vector<std::vector<int>> A(st_neznank, std::vector<int>(st_neznank, 0));
    std::vector<int> b(st_neznank, 0);
    int deltaX = 1;
    int deltaY = 1;
    int k = 1;


    for (size_t i = 0; i < st_neznank; ++i) {

        std::vector<int> sosedi = sosednja_vozlisca[i];
        int levi_sosed = sosedi[0];
        int spodnji_sosed = sosedi[1];
        int desni_sosed = sosedi[2];
        int zgornji_sosed = sosedi[3];
        if (levi_sosed != -1 && spodnji_sosed != -1 && desni_sosed != -1 && zgornji_sosed != -1) {

            A[i][levi_sosed] = 1;
            A[i][spodnji_sosed] = 1;
            A[i][desni_sosed] = 1;
            A[i][zgornji_sosed] = 1;
            A[i][i] = -4;

        }
        else {

            int tip_robnega_pogoja = 0;
            int vrednost = 0;
            int vrednost_prestopa = 0;
            for (size_t j = 0; j < 4; ++j) {

                std::vector<int> vozlisca_v_trenutnem_rp = vozlisca_robnih_pogojev[j];
                for (size_t k = 0; k < vozlisca_v_trenutnem_rp.size(); ++k) {

                    int vozlisce_v_trenutnem_rp = vozlisca_v_trenutnem_rp[k];
                    if (i == vozlisce_v_trenutnem_rp) {

                        tip_robnega_pogoja = tipi_robnih_pogojev[j];
                        vrednost = vrednosti_robnih_pogojev[j];
                        vrednost_prestopa = vrednosti_prestopa_toplote[j];

                    }

                }

            }
            if (tip_robnega_pogoja == 0) {

                A[i][i] = 1;
                b[i] = vrednost;

            }
            else if (tip_robnega_pogoja == 1) {

                int stevilo_sosedov = 0;
                for (size_t l = 0; l < 4; ++l) {

                    if (sosedi[l] != -1) {
                        stevilo_sosedov = stevilo_sosedov + 1;
                    }

                }
                if (stevilo_sosedov == 3) {

                    if (levi_sosed == -1) {

                        A[i][i] = A[i][i] - 4;
                        A[i][desni_sosed] = A[i][desni_sosed] + 2;
                        A[i][spodnji_sosed] = A[i][spodnji_sosed] + 1;
                        A[i][zgornji_sosed] = A[i][zgornji_sosed] + 1;
                        b[i] = -2 * (vrednost * deltaX / k);

                    }
                    if (desni_sosed == -1) {

                        A[i][i] = A[i][i] - 4;
                        A[i][levi_sosed] = A[i][levi_sosed] + 2;
                        A[i][spodnji_sosed] = A[i][spodnji_sosed] + 1;
                        A[i][zgornji_sosed] = A[i][zgornji_sosed] + 1;
                        b[i] = -2 * (vrednost * deltaX / k);

                    }
                    if (spodnji_sosed == -1) {

                        A[i][i] = A[i][i] - 4;
                        A[i][zgornji_sosed] = A[i][zgornji_sosed] + 2;
                        A[i][levi_sosed] = A[i][levi_sosed] + 1;
                        A[i][desni_sosed] = A[i][desni_sosed] + 1;
                        b[i] = -2 * (vrednost * deltaX / k);

                    }
                    if (zgornji_sosed == -1) {

                        A[i][i] = A[i][i] - 4;
                        A[i][spodnji_sosed] = A[i][spodnji_sosed] + 2;
                        A[i][levi_sosed] = A[i][levi_sosed] + 1;
                        A[i][desni_sosed] = A[i][desni_sosed] + 1;
                        b[i] = -2 * (vrednost * deltaX / k);

                    }

                }

            }
            else if (tip_robnega_pogoja == 2) {

                int stevilo_sosedov = 0;
                for (size_t l = 0; l < 4; ++l) {

                    if (sosedi[l] != -1) {
                        stevilo_sosedov = stevilo_sosedov + 1;
                    }

                }
                if (stevilo_sosedov == 3) {

                    if (levi_sosed == -1) {

                        A[i][i] = A[i][i] - 2 * (vrednost_prestopa * deltaX / k + 2);
                        A[i][desni_sosed] = A[i][desni_sosed] + 2;
                        A[i][spodnji_sosed] = A[i][spodnji_sosed] + 1;
                        A[i][zgornji_sosed] = A[i][zgornji_sosed] + 1;
                        b[i] = b[i] - 2 * vrednost_prestopa * deltaX * vrednost / k;

                    }
                    if (desni_sosed == -1) {

                        A[i][i] = A[i][i] - 2 * (vrednost_prestopa * deltaX / k + 2);
                        A[i][levi_sosed] = A[i][levi_sosed] + 2;
                        A[i][spodnji_sosed] = A[i][spodnji_sosed] + 1;
                        A[i][zgornji_sosed] = A[i][zgornji_sosed] + 1;
                        b[i] = b[i] - 2 * vrednost_prestopa * deltaX * vrednost / k;

                    }
                    if (spodnji_sosed == -1) {

                        A[i][i] = A[i][i] - 2 * (vrednost_prestopa * deltaX / k + 2);
                        A[i][levi_sosed] = A[i][levi_sosed] + 1;
                        A[i][desni_sosed] = A[i][desni_sosed] + 1;
                        A[i][zgornji_sosed] = A[i][zgornji_sosed] + 2;
                        b[i] = -2 * vrednost_prestopa * deltaX * vrednost / k;

                    }
                    if (zgornji_sosed == -1) {

                        A[i][i] = A[i][i] - 2 * (vrednost_prestopa * deltaX / k + 2);
                        A[i][levi_sosed] = A[i][levi_sosed] + 1;
                        A[i][desni_sosed] = A[i][desni_sosed] + 1;
                        A[i][spodnji_sosed] = A[i][spodnji_sosed] + 2;
                        b[i] = -2 * vrednost_prestopa * deltaX * vrednost / k;

                    }

                }

            }

        }

    }


    //Resitev sistema enacb
    std::vector<double> T(st_neznank, 100.0);
    int iteracije = 1000;
    for (int iitt = 0; iitt < iteracije; ++iitt) {

        for (int jj = 0; jj < st_neznank; ++jj) {

            double d = b[jj];
            for (int ii = 0; ii < st_neznank; ++ii) {

                if (jj != ii) {

                    d = d - A[jj][ii] * T[ii];

                }
                T[jj] = d / A[jj][jj];

            }

        }

    }

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    std::cout << "Program se je izvajal: " << duration.count() / 1000000.00 << " sekund" << std::endl;

    //VTKfile

    const std::string vtkFileName = "resitev4mreza.vtk";

    std::ofstream vtkFile(vtkFileName);

    if (!vtkFile.is_open()) {
        std::cerr << "Error: Unable to open the VTK file for writing." << std::endl;
        return 1;
    }

    vtkFile << "# vtk DataFile Version 3.0\n";
    vtkFile << "Mreza1\n";
    vtkFile << "ASCII\n";
    vtkFile << "DATASET UNSTRUCTURED_GRID\n";

    vtkFile << "POINTS " << st_neznank << " float\n";
    for (int i = 0; i < st_neznank; ++i) {
        vtkFile << X[i] << " " << Y[i] << " " << 0 << "\n";
    }
    vtkFile << "\n";
    vtkFile << "CELLS " << st_celic << " " << st_celic * 5 << "\n";
    for (int i = 0; i < st_celic; ++i) {
        vtkFile << celice[i].size() << " " << celice[i][0] << " " << celice[i][1] << " " << celice[i][2] << " " << celice[i][3] << "\n";
    }
    vtkFile << "\n";
    vtkFile << "CELL_TYPES " << st_celic << "\n";
    for (int i = 0; i < st_celic; ++i) {
        vtkFile << 9 << "\n";
    }
    vtkFile << "\n";
    vtkFile << "POINT_DATA" << " " << st_neznank << "\n";
    vtkFile << "SCALARS Temperature float 1" << "\n";
    vtkFile << "LOOKUP_TABLE default" << "\n";
    for (int i = 0; i < st_neznank; ++i) {
        vtkFile << T[i] << "\n";
    }

    vtkFile.close();

    std::cout << vtkFileName << " se nahaja v mapi" << std::endl;

    std::cout << "" << std::endl;
    std::cout << "Vrednosti T" << std::endl;
    for (int i = 0; i < st_neznank; ++i) {
        std::cout << "T[" << i << "] = " << T[i] << std::endl;
    }

    // Ustvarimo JSON datoteko
    json jsonData;
    jsonData["Array"] = A;
    jsonData["Vector"] = b;

    // Pišemo v datoteko
    std::ofstream outputFile("output.json");
    outputFile << std::setw(4) << jsonData << std::endl;
    outputFile.close();
    std::cout << "JSON podatki so zapisani v datoteki output.json" << std::endl;

    return 0;
}