//
//  Mesh.h
//  Mes Project
//
//  Created by Krystian Zapart 
//

#ifndef grid_h
#define grid_h

#include <iostream>
#include <cmath>
#include "Define.h"

struct Mesh
{

    Node* nodes;
    Element* elements;
    double** mat_N;//macierz funkcji kształtu
    double** global_H;
    double** global_C;
    double* global_P;
    double* T0;// wektor temp pocz
    double* final_P;

    int nN = nH * nB;
    int nE = (nH - 1) * (nB - 1);
    Mesh()
    {
        nodes = new Node[nN];
        elements = new Element[nE];

       
        double dx = b / (nB - 1), dy = h / (nH - 1);
        int k = 0;
        for (int i = 0; i < nB; ++i)
            for (int j = 0; j < nH; ++j)
            {
                nodes[k] = Node(k + 1, dx * i, dy * j);
                k++;
            }

        int ID = 1;
        for (int i = 0; i < nE; ++i)
        {
            if (i != 0)
                if (i % (nB - 1) == 0)
                    ID++;
            elements[i] = Element(i + 1, ID, nH, nodes, nN);
            ID = elements[i].nodes[3].id;
        }

        global_P = new double[nN];
        final_P = new double[nN];
        global_H = new double* [nN];
        global_C = new double* [nN];
        T0 = new double[nN];
        for (int i = 0; i < nN; ++i)
        {
            global_P[i] = 0.0;
            final_P[i] = 0.0;
            T0[i] = 100.0;
            global_H[i] = new double[nN];
            global_C[i] = new double[nN];
            for (int j = 0; j < nN; ++j)
            {
                global_H[i][j] = 0.0;
                global_C[i][j] = 0.0;
            }
        }
    }

    Mesh(double wys, double szer, int nodePion, int nodePoziom, int integrationPoints, int tPocz, int krokCzasowy, int tOtocz, int wspAlfa, int cieploWlasciwe, int cond, int gestosc)
    {
        //przypisanie wartosci
        nN = nH * nB;
        nE = (nH - 1) * (nB - 1);

        nodes = new Node[nN];
        elements = new Element[nE];
        //stworzenie node'ow
        double dx = b / (nodePoziom - 1);
        double dy = h / (nodePion - 1);

        int g = 0;
        for (int i = 0; i < nB; i++) {
            for (int j = 0; j < nH; j++) {
                //nodes[g].uzupelnijNode(g + 1, dx * i, dy * j, h, b);
                nodes[g] = Node(g + 1, dx * i, dy * j);
                g++;
            }
        }
        //stworzenie elementow
        
        int ID = 1;
        for (int i = 0; i < nE; i++)
        {
            if (i != 0) {
                if (i % (nH - 1) == 0)
                {
                    ID++;
                }
            }
            //elements[i].uzupelnijElement(i + 1, ID, nH);
            elements[i] = Element(i + 1, ID, nH, nodes, nN);
            ID = elements[i].nodes[3].id;
        }
        //Tworzenie globalnych wektorów i macierzy
        global_H = new double*[nN];
        global_C = new double*[nN];

        for (int i = 0; i < nN; i++)
        {
            global_C[i] = new double[nN];
            global_H[i] = new double[nN];
        }

        global_P = new double[nN];
        T0 = new double[nN];
        final_P = new double[nN];

        //Ustawienie wartości początkowych dla stworzonych wektorów i macierzy
        for (int i = 0; i < nN; i++)
        {
            global_P[i] = 0.0;
            T0[i] = tPocz;
            final_P[i] = 0.0;

            for (int j = 0; j < nN; j++)
            {
                global_C[i][j] = 0.0;
                global_H[i][j] = 0.0;
            }
        }

    }

//    void show_Elementy()
//    {
//        std::cout << "\nElementy:\n";
//        for (int i = 0; i < nE; i++)
//        {
//            //elements[i].showElement();
//        }
//        std::cout << "\n";
//    }

//    void show_Nodes()
//    {
//        std::cout << "\nNode'y\n";
//        for (int i = 0; i < nN; i++)
//        {
//            nodes[i].showNodeNumer();
//        }
//        std::cout << "\n";
//    }

//    void show_Mesh()
//    {
//        std::cout << "\nmesh:\n";
//        for (int i = nH; i > 0; i--)
//        {
//            int k = i - 1;
//            for (int j = 0; j < nB; j++)
//            {
//                nodes[k].showNode();
//                k += nH;
//            }
//            std::cout << "\n";
//        }
//        std::cout << "\n";
//
//
//        std::cout<< "Elementy: \n";
//        for(int i = 0; i < (nB - 1)*(nH-1); i++)
//        {
////            elements[i].show_Element();
//        }
//        std::cout << "\n";
//
//        std::cout << "Node'y:\n";
//        for (int i = 0; i < nN; i++)
//        {
//            nodes[i].showNodeNumer();
//        }
//        std::cout << "\n";
//
//    }

    void fJakobian(double** dKsi, double** dEta)
    {

        double** Jakobian = new double* [2];
        double** Jakobian_inv = new double* [2];
        double** dx = new double* [pc * pc];
        double** dy = new double* [pc * pc];


        Jakobian[0] = new double[2];
        Jakobian[1] = new double[2];
        Jakobian_inv[0] = new double[2];
        Jakobian_inv[1] = new double[2];

        for (int i = 0; i < pc * pc; i++)
        {
            dx[i] = new double[4];
            dy[i] = new double[4];
        }

        for (int i = 0; i < nE; i++)
        {
            double w1 = 5.0 / 9.0, w2 = 8.0 / 9.0, w3 = w1, w4 = w1;
            
            for (int j = 0; j < pc * pc; j++)
            {
                //Macierz H
                count_Jakobian(i, j, dEta, dKsi, Jakobian);
                count_inv_Jakobian(Jakobian, Jakobian_inv);
                count_DxDy(dx, dy, dEta, dKsi, Jakobian_inv, j);
                count_H(elements[i].H, dx, dy, Jakobian, j, w1 * w4);

                //Macierz C
                elements[i].count_C(j, mat_N, Jakobian, w1 * w4);
    
                w1 = w2;
                w2 = w3;
                w3 = w1;

                if (j == 7)
                {
                    //środkowy element
                    w1 = 8.0 / 9.0;
                    w4 = w1;
                }
            }

//            final_H(elements[i].H, Jakobian, conductivity);
//            elements[i].printH();
        }
    }

    void count_Jakobian(int i, int j, double** dEta, double** dKsi, double** Jakobian)
    {
        for (int i = 0; i < 2; i++)
        {
            for (int j = 0; j < 2; j++)
            {
                Jakobian[i][j] = 0;
            }
        }

        for (int k = 0; k < 4; k++)
        {
            Jakobian[0][0] += dKsi[j][k] * elements[i].nodes[k].y;
            //std::cout << Jakobian[0][0] << "\n";
            Jakobian[0][1] += dEta[j][k] * elements[i].nodes[k].y;
            //std::cout << Jakobian[0][1] << "\n";
            Jakobian[1][0] += dKsi[j][k] * elements[i].nodes[k].x;
            //std::cout << Jakobian[1][0] << "\n";
            Jakobian[1][1] += dEta[j][k] * elements[i].nodes[k].x;
            //std::cout << Jakobian[1][1] << "\n";
        }


    }

    void count_inv_Jakobian(double** Jakobian, double** Jakobian_rev)
    {
        double det = (Jakobian[0][0] * Jakobian[1][1] - Jakobian[0][1] * Jakobian[1][0]);
        //std::cout << "det: " << det << "\n";
        double det_rev = 1 / det;
        //std::cout << "det _rev: " << det_rev << "\n";
        Jakobian_rev[0][0] = det_rev * Jakobian[0][0];
        Jakobian_rev[0][1] = det_rev * Jakobian[0][1];
        Jakobian_rev[1][0] = det_rev * Jakobian[1][0];
        Jakobian_rev[1][1] = det_rev * Jakobian[1][1];

        /*
        std::cout << "TEST JAKOBIAN _rev:\n\n";
        std::cout << Jakobian_rev[0][0] << "\n";
        std::cout << Jakobian_rev[0][1] << "\n";
        std::cout << Jakobian_rev[1][0] << "\n";
        std::cout << Jakobian_rev[1][1] << "\n";
        std::cout << "\n\n";
        */
    }

    void count_DxDy(double** dx, double** dy, double** dEta, double** dKsi, double** Jakobian_rev, int i)
    {
        for (int j = 0; j < 4; j++)
        {
            dx[i][j] = Jakobian_rev[0][0] * dEta[i][j] + Jakobian_rev[0][1] * dKsi[i][j];
            dy[i][j] = Jakobian_rev[1][0] * dEta[i][j] + Jakobian_rev[1][1] * dKsi[i][j];
            //std::cout << "\n\nIteracja: (" << i << "," << j << ")\n dX = " << dX[i][j] << "\t" << "dY = " << dY[i][j] << "\n";
        }
    }

    void count_H(double** H, double** dx, double** dy, double** Jakobian, int a, double w)
    {
        double det = Jakobian[0][0] * Jakobian[1][1] - Jakobian[0][1] * Jakobian[1][0];
        if (pc == 2)
        {
//            std::cout<<"________________H_________________"<<std::endl;
            for (int i = 0; i < 4; i++)
            {
                for (int j = 0; j < 4; j++)
                {
                    H[i][j] += (dx[a][i] * dx[a][j]) *det *conductivity;
                    H[i][j] += (dy[a][i] * dy[a][j]) *det *conductivity;
//                    std::cout << h[i][j] << " "<<std::endl;
                }
            }
        }
        else if (pc == 3)
        {
            
            for (int i = 0; i < 4; i++)
            {
                for (int j = 0; j < 4; j++)
                {
                    H[i][j] += (dx[a][i] * dx[a][j]) * w *det *conductivity;
                    H[i][j] += (dy[a][i] * dy[a][j]) * w *det *conductivity;
                }
            }
        }
    }

    void agreg_H()
    {
        double** agreg_tab = new double* [nN];
        for (int i = 0; i < nN; i++)
        {
            agreg_tab[i] = new double[nN];
        }

        for (int i = 0; i < nE; i++)
        {
            for (int i = 0; i < nN; i++)
            {
                for (int j = 0; j < nN; j++)
                {
                    agreg_tab[i][j] = 0;
                }
            }
            elements[i].agreg(agreg_tab);
            add_To_Agregation(agreg_tab);
        }


        //wyświetlGlobalneH();
    }

//    void show_H_global()
//    {
//        std::cout << "Global H: \n";
//        for (int i = 0; i < nN; i++)
//        {
//            for (int j = 0; j < nN; j++)
//            {
//                std::cout << global_H[i][j] << "\t";
//            }
//            std::cout << "\n";
//        }
//        std::cout << "\n";
//    }

//    void show_H_final()
//    {
//        // [H] = [H] + [C]/dT
//        std::cout << "final H: \n";
//        for (int i = 0; i < nN; i++)
//        {
//            for (int j = 0; j < nN; j++)
//            {
//                std::cout << global_H[i][j] << "\t";
//            }
//            std::cout << "\n";
//        }
//        std::cout << "\n";
//    }

    void add_To_Agregation(double** tab)
    {
        for (int i = 0; i < nN; i++)
        {
            for (int j = 0; j < nN; j++)
            {
                global_H[i][j] += tab[i][j];
            }
        }
    }

    void add_To_Global_vector(double* vec)
    {
        for (int i = 0; i < nN; i++)
        {
            global_P[i] += vec[i];
        }
    }

    void agreg_C()
    {
        double** tmp_agreg_tab = new double* [nN];
        for (int i = 0; i < nN; i++)
        {
            tmp_agreg_tab[i] = new double[nN];
        }

        for (int i = 0; i < nE; i++)
        {
            for (int i = 0; i < nN; i++)
            {
                for (int j = 0; j < nN; ++j)
                {
                    tmp_agreg_tab[i][j] = 0;
                }
            }
            elements[i].agreg_C(tmp_agreg_tab);
            add_to_C_glob(tmp_agreg_tab);
        }
    }

//    void show_C()
//    {
//        std::cout << "Global C: " << "\n";
//        for (int i = 0; i < nN; i++)
//        {
//            for (int j = 0; j < nN; j++)
//            {
//                std::cout << global_C[i][j] << "\t";
//            }
//            std::cout << "\n";
//        }
//        std::cout << "\n";
//    }

    void add_to_C_glob(double** tab)
    {
        for (int i = 0; i < nN; i++)
        {
            for (int j = 0; j < nN; j++)
            {
                global_C[i][j] += tab[i][j];
           }
        }
    }

    void sum_H_C()
    {
        for (int i = 0; i < nN; i++)
        {
            for (int j = 0; j < nN; j++)
            {
                global_H[i][j] += global_C[i][j] / sim_step;
            }
        }
        //wyświetlKońcoweH();
    }

//    void show_P_Global()
//    {
//        std::cout << "Global P: " << "\n";
//        for (int i = 0; i < nN; i++)
//        {
//            std::cout << global_P[i] << "\t";
//        }
//        std::cout << "\n\n";
//    }

    void count_final_P()
    {
        //{P} + {[C]/dT * {T0}}
        //std::cout << "--------------------------------\n\n";
        double suma;
        for (int i = 0; i < nN; i++)
        {
            suma = 0.0;
            for (int j = 0; j < nN; j++)
            {
                suma += (global_C[i][j] / sim_step) * T0[j];
            }
            final_P[i] = global_P[i] + suma;

        }
        //std::cout << "--------------------------------\n\n";
    }
};

#endif /* grid_h */
