#include <vector>
#include <utility>
#include <string>
#include <complex>
#include <tuple>
#include <iostream>
#include <fstream>

#include "utils.h"

using std::ofstream;

bool RepetitionCheck(std::vector<int>, int);

std::vector<int> FindNeighbours(int[], std::vector<int>, int, int, std::ostream&);

std::pair<int, int> Next(std::vector<int>, int[], int);

double TfromTDC(double t1, double t2, double L);
double AttenuationFactor(double d, int planeID);
double EfromADC(double adc1, double adc2, double d1, double d2, int planeID);

//Preclustering
void Preclustering(std::string input)
{
    const char* finname = input.c_str();
    TFile f(finname, "READ");

    TTree* t = (TTree*)f.Get("tDigit");
    Int_t TotHits = t->GetEntries();
    // Access Calorimeter leaves
    TLeaf* cell_id = (TLeaf*)t->GetLeaf("dg_cell.id");
    
    TString filename = input;
    filename.ReplaceAll(".root", ".preclustered.txt");
    ofstream out;
    out.open(filename);
    //loop over events
    for (int i = 0; i <2; i++) {
        t->GetEntry(i);
        int lentry = cell_id->GetLen();
        int digits[lentry];
        std::vector<int> checked_array;
        int n_cluster = 0, cluster[lentry];
        // Loop over event digits
        for (int j = 0; j < lentry; j++) {
            digits[j] = cell_id->GetValue(j);
        }
        out<<"999"<<endl;
//        cout << "Evento: " << i + 1 << " Digits: " << lentry << endl;
        for (int k = 0; k < lentry; k++) {
            if (k == 0) {
                cluster[k] = digits[0];
            }
            else {
                cluster[k] = Next(checked_array, digits, lentry).first;
                n_cluster = Next(checked_array, digits, lentry).second;
                if (cluster[k] == 0) {
                    break;
                }
            }
//            cout<< cluster[k] << " ("<<n_cluster<<") ";
            out << n_cluster << " ";

            checked_array.push_back(cluster[k]);
            checked_array = FindNeighbours(digits, checked_array, lentry, cluster[k], out);
 //           cout << endl;
            out << "888"<<endl;
        }
        //for (std::vector<int>::const_iterator i = checked_array.begin(); i != checked_array.end(); ++i)
//            std::cout << *i << ' ';
//        cout << endl;
//        cout << "------------" << endl;
 //       out << "-" << endl;
    }
}


// Clustering
int Clustering(std::string input)
{
    const char* finname = input.c_str();
    TFile f(finname, "READ");
    TTree* t = (TTree*)f.Get("tDigit");
    TLeaf* cell_id = (TLeaf*)t->GetLeaf("dg_cell.id");
    TLeaf* cell_adc1 = (TLeaf*)t->GetLeaf("dg_cell.adc1");
    TLeaf* cell_x = (TLeaf*)t->GetLeaf("dg_cell.x");
    TLeaf* cell_y = (TLeaf*)t->GetLeaf("dg_cell.y");
    TLeaf* cell_z = (TLeaf*)t->GetLeaf("dg_cell.z");
    TLeaf* cell_tdc1 = (TLeaf*)t->GetLeaf("dg_cell.tdc1");
    TLeaf* cell_adc2 = (TLeaf*)t->GetLeaf("dg_cell.adc2");
    TLeaf* cell_tdc2 = (TLeaf*)t->GetLeaf("dg_cell.tdc2");
    TString filename = input;
    filename.ReplaceAll(".root", ".preclustered.txt");
    std::ifstream LUT(filename);
    int n_event = 0;
    int n_line = 0;
    int clusteradc12 = 0;
    Double_t x_weighted, y_weighted, z_weighted, Etot;
    x_weighted = 0;
    y_weighted = 0;
    z_weighted = 0;
    for (std::string line; getline(LUT, line); )
    {
        n_line++;
        std::istringstream in(line);
        int value;
        while (in >> value) {
            if (value == 888) {
                cout << "Cluster Energy: " << clusteradc12 << endl;
                cout << "X_weigh: " << x_weighted / Etot << " Y_weigh: "<<y_weighted/Etot<<" Z_weigh: "<<z_weighted/Etot<<endl;
                cout << "------" << endl;
                x_weighted = 0;
                y_weighted = 0;
                z_weighted = 0;
                Etot = 0;
                clusteradc12 = 0;
            }
            else if (value == 999) {
//               cout << "Evento numero: " << n_event << " inizia a riga: "<<n_line<<endl;
                t->GetEntry(n_event);
                //cout << cell_id->GetLen() << endl;
                n_event++;
                cout << "*******************" << endl;
            }           
            else 
            {
                //cout << cell_adc1->GetValue(value) + cell_adc2->GetValue(value) <<" ";
                x_weighted = x_weighted + (cell_x->GetValue(value) * (cell_adc1->GetValue(value) + cell_adc2->GetValue(value)));
                y_weighted = y_weighted + (cell_y->GetValue(value) * (cell_adc1->GetValue(value) + cell_adc2->GetValue(value)));
                z_weighted = z_weighted + (cell_z->GetValue(value) * (cell_adc1->GetValue(value) + cell_adc2->GetValue(value)));
                Etot = Etot + cell_adc1->GetValue(value) + cell_adc2->GetValue(value);
                clusteradc12 = clusteradc12 + cell_adc1->GetValue(value) + cell_adc2->GetValue(value);
                //cout << cell_adc1->GetValue(value) + cell_adc2->GetValue(value) << " ";
                cout << " Z: " << cell_z->GetValue(value) << " X: "<< cell_x->GetValue(value)<<" Y: "<< cell_y->GetValue(value)<<" -- ";
            }
        }
    }
    return 0;
}


std::vector<int> FindNeighbours(int digits[], std::vector<int> already_checked, int size, int check, std::ostream& out) {
    for (int i = 0; i < size; i++) {
        int id_check = digits[i];
        if (id_check == check) {
            continue;
        }
        else if (id_check == check + 101 || id_check == check + 100 || id_check == check + 99 || id_check == check + 1 || id_check == check - 1 || id_check == check - 99 || id_check == check - 100 || id_check == check - 101) {
            bool RepCheck = RepetitionCheck(already_checked, id_check);
            if (RepCheck == true) {
                continue;
            } 
            else {
//                cout << id_check <<"  ("<<i<<") ";
                out <<i<<" ";
                already_checked.push_back(id_check);
                already_checked=FindNeighbours(digits, already_checked, size, id_check, out);
            }
        }
        else if (id_check == check - 889 || id_check == check - 989 || id_check == check - 1089) {
            bool RepCheck = RepetitionCheck(already_checked, id_check);
            if (RepCheck == true) {
                continue;
            } 
            else {
//                cout << id_check << " ("<<i<<") ";
                out << i << " ";
                already_checked.push_back(id_check);
                already_checked=FindNeighbours(digits, already_checked, size, id_check, out);
            }
        }
        else if ( id_check == check + 23111 ||id_check == check + 23011 || id_check == check + 22911) {
            bool RepCheck = RepetitionCheck(already_checked, id_check);
            if (RepCheck == true) {
                continue;
            }
            else {
//                cout << id_check << " ("<<i<<") ";
                out << i << " ";
                already_checked.push_back(id_check);
                already_checked = FindNeighbours(digits, already_checked, size, id_check, out);
            }
        }
        else if (id_check == check + 1089 || id_check == check + 989 || id_check == check + 911) {
            bool RepCheck = RepetitionCheck(already_checked, id_check);
            if (RepCheck == true) {
                continue;
            }
            else {
//                cout << id_check << " (" << i << ") ";
                out << i << " ";
                already_checked.push_back(id_check);
                already_checked = FindNeighbours(digits, already_checked, size, id_check, out);
            }
        }
        else if (id_check == check -23011 || id_check == check -22911 || id_check == check -23111) {
            bool RepCheck = RepetitionCheck(already_checked, id_check);
            if (RepCheck == true) {
                continue;
            }
            else {
//                cout << id_check << " ("<<i<<") ";
                out << i << " ";
                already_checked.push_back(id_check);
                already_checked = FindNeighbours(digits, already_checked, size, id_check, out);
            }
        }
    }
    return already_checked;
}

bool RepetitionCheck(std::vector<int> v, int check) {
    if (std::find(v.begin(), v.end(), check) != v.end()) {
        return true;
    }
    else {
        return false;
    }
}

std::pair<int,int> Next(std::vector<int> already_checked, int digits[], int size) {
    int next=0;
    int entry = 0;
    for (int i = 0; i < size; i++) {
        bool RepCheck = RepetitionCheck(already_checked, digits[i]);
        if (RepCheck == true) {
            continue;
        }
        else {
            next = digits[i];
            entry = i;
            break;
        }
    }
    return std::make_pair(next, entry);
}

double kloe_simu::AttenuationFactor(double d, int planeID)
{
    /*
         dE/dx attenuation - Ea=p1*exp(-d/atl1)+(1.-p1)*exp(-d/atl2)
           d    distance from photocatode - 2 cells/cell; d1 and d2
          atl1  50. cm
          atl2  430 cm planes 1-2    innermost plane is 1
                380 cm plane 3
                330 cm planes 4-5
           p1   0.35
    */
    double atl2 = 0.0;

    switch (planeID) {
    case 0:
    case 1:
        atl2 = kloe_simu::atl2_01;
        break;

    case 2:
        atl2 = kloe_simu::atl2_2;
        break;

    case 3:
    case 4:
        atl2 = kloe_simu::atl2_34;
        break;

    default:
        // std::cout << "planeID out if range" << std::endl;
        atl2 = -999.0;
        break;
    }

    return kloe_simu::p1 * TMath::Exp(-d / kloe_simu::atl1) +
        (1. - kloe_simu::p1) * TMath::Exp(-d / atl2);
}

// reconstruct t of the hit from tdc1 and tdc2
double kloe_simu::TfromTDC(double t1, double t2, double L)
{
    return 0.5 * (t1 + t2 - kloe_simu::vlfb * L / kloe_simu::m_to_mm);
}

// energy deposit of the hit from adc1 and adc2 and
// reconstructed longidutinal coordinate
double kloe_simu::EfromADC(double adc1, double adc2, double d1, double d2,
    int planeID)
{
    double f1 = AttenuationFactor(d1, planeID);
    double f2 = AttenuationFactor(d2, planeID);

    return 0.5 * (adc1 / f1 + adc2 / f2) * kloe_simu::adc2MeV;
}
