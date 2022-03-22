//
//  element_Univ.h
//  Mes Project
//
//  Created by Krystian Zapart
//

#ifndef element_Univ_h
#define element_Univ_h

#include <iostream>
#include "wall.h"

struct Element_Univ
{
    double** dKsi;
    double** dEta;
    double** nc;

    wall Lwall, Rwall, Uwall, Dwall;
    Mesh* mesh;

    Element_Univ(Mesh* mesh)
    {
        this->mesh = mesh;

        dKsi = new double* [pc * pc];
        dEta = new double* [pc * pc];
        nc = new double* [pc * pc];

        for (int i = 0; i < pc * pc; i++)
        {
            dKsi[i] = new double[4];
            dEta[i] = new double[4];
            nc[i] = new double[4];
        }

        //Uzupełnianie funkcji kształtu
        if (pc == 2)
        {
            calc_N(nc,-(1.0/sqrt(3.0)),-(1.0/sqrt(3.0)));
        }
        if (pc == 3)
        {
            calc_N(nc, sqrt(3.0 / 5.0), 1.0);
        }
        mesh->mat_N = nc;
        
        //Przypisanie punktów do ścian
        if (pc == 2)
        {
            Lwall = wall(nH, -1.0, (1.0/sqrt(3.0)), -1.0, -(1.0/sqrt(3.0)));
            Rwall = wall(nH, 1.0, -(1.0/sqrt(3.0)), 1.0, (1.0/sqrt(3.0)));
            Uwall = wall(nB, -(1.0/sqrt(3.0)), 1.0, (1.0/sqrt(3.0)), 1.0);
            Dwall = wall(nB, -(1.0/sqrt(3.0)), -1.0, (1.0/sqrt(3.0)), -1.0);
        }

        if (pc == 3)
        {
            Lwall = wall(nH, -1.0, (sqrt(3.0/5.0)), -1.0, 0, -1.0, -(sqrt(3.0/5.0)));
            Rwall = wall(nH, 1.0, -(sqrt(3.0/5.0)), 1.0, 0, 1.0, (sqrt(3.0/5.0)));
            Uwall = wall(nB, (sqrt(3.0/5.0)), 1.0, 0, 1.0, -(sqrt(3.0/5.0)), 1.0);
            Dwall = wall(nB, -(sqrt(3.0/5.0)), -1.0, 0, -1.0, (sqrt(3.0/5.0)), -1.0);
        }
        calc_walls(mesh->elements);
    }
//uzupełnienei macierzy funkcji kształtu
    void calc_N(double** tab, double ksi, double eta)
    {

        if (pc == 2)
        {
            tab[0][0] = 0.25 * (1.0 - ksi) * (1 - eta);
            tab[0][1] = 0.25 * (1 + ksi) * (1 - eta);
            tab[0][2] = 0.25 * (1 + ksi) * (1 + eta);
            tab[0][3] = 0.25 * (1 - ksi) * (1 + eta);

            ksi = -ksi;
            tab[1][0] = 0.25 * (1.0 - ksi) * (1 - eta);
            tab[1][1] = 0.25 * (1 + ksi) * (1 - eta);
            tab[1][2] = 0.25 * (1 + ksi) * (1 + eta);
            tab[1][3] = 0.25 * (1 - ksi) * (1 + eta);

            eta = -eta;
            tab[2][0] = 0.25 * (1.0 - ksi) * (1 - eta);
            tab[2][1] = 0.25 * (1 + ksi) * (1 - eta);
            tab[2][2] = 0.25 * (1 + ksi) * (1 + eta);
            tab[2][3] = 0.25 * (1 - ksi) * (1 + eta);

            ksi = -ksi;
            tab[3][0] = 0.25 * (1.0 - ksi) * (1 - eta);
            tab[3][1] = 0.25 * (1 + ksi) * (1 - eta);
            tab[3][2] = 0.25 * (1 + ksi) * (1 + eta);
            tab[3][3] = 0.25 * (1 - ksi) * (1 + eta);
        }
        
        if (pc == 3)
        {
            
            double a = ksi;

            ksi = -a;
            eta = -a;
            tab[0][0] = 0.25 * (1 - ksi) * (1 - eta);
            tab[0][1] = 0.25 * (1 + ksi) * (1 - eta);
            tab[0][2] = 0.25 * (1 + ksi) * (1 + eta);
            tab[0][3] = 0.25 * (1 - ksi) * (1 + eta);
//
//            std::cout<<"kształt\n";
//            for(int i = 0;i<4;i++)
//            {
//                std::cout<<tab[0][i]<<" ";
//            }
//            std::cout<<std::endl;
            
            ksi = 0;
            tab[1][0] = 0.25 * (1 - ksi) * (1 - eta);
            tab[1][1] = 0.25 * (1 + ksi) * (1 - eta);
            tab[1][2] = 0.25 * (1 + ksi) * (1 + eta);
            tab[1][3] = 0.25 * (1 - ksi) * (1 + eta);
//
//            for(int i = 0;i<4;i++)
//            {
//                std::cout<<tab[1][i]<<" ";
//            }
//            std::cout<<std::endl;
//
            ksi = a;
            tab[2][0] = 0.25 * (1 - ksi) * (1 - eta);
            tab[2][1] = 0.25 * (1 + ksi) * (1 - eta);
            tab[2][2] = 0.25 * (1 + ksi) * (1 + eta);
            tab[2][3] = 0.25 * (1 - ksi) * (1 + eta);
            
//            for(int i = 0;i<4;i++)
//            {
//                std::cout<<tab[2][i]<<" ";
//            }
//            std::cout<<std::endl;
            
            eta = 0;
            tab[3][0] = 0.25 * (1 - ksi) * (1 - eta);
            tab[3][1] = 0.25 * (1 + ksi) * (1 - eta);
            tab[3][2] = 0.25 * (1 + ksi) * (1 + eta);
            tab[3][3] = 0.25 * (1 - ksi) * (1 + eta);
            
//            for(int i = 0;i<4;i++)
//            {
//                std::cout<<tab[3][i]<<" ";
//            }
//            std::cout<<std::endl;
            
            eta = a;
            tab[4][0] = 0.25 * (1 - ksi) * (1 - eta);
            tab[4][1] = 0.25 * (1 + ksi) * (1 - eta);
            tab[4][2] = 0.25 * (1 + ksi) * (1 + eta);
            tab[4][3] = 0.25 * (1 - ksi) * (1 + eta);
            
//            for(int i = 0;i<4;i++)
//            {
//                std::cout<<tab[4][i]<<" ";
//            }
//            std::cout<<std::endl;
            
            ksi = 0;
            tab[5][0] = 0.25 * (1 - ksi) * (1 - eta);
            tab[5][1] = 0.25 * (1 + ksi) * (1 - eta);
            tab[5][2] = 0.25 * (1 + ksi) * (1 + eta);
            tab[5][3] = 0.25 * (1 - ksi) * (1 + eta);
            
//            for(int i = 0;i<4;i++)
//            {
//                std::cout<<tab[5][i]<<" ";
//            }
//            std::cout<<std::endl;
            
            ksi = -a;
            tab[6][0] = 0.25 * (1 - ksi) * (1 - eta);
            tab[6][1] = 0.25 * (1 + ksi) * (1 - eta);
            tab[6][2] = 0.25 * (1 + ksi) * (1 + eta);
            tab[6][3] = 0.25 * (1 - ksi) * (1 + eta);
            
//            for(int i = 0;i<4;i++)
//            {
//                std::cout<<tab[6][i]<<" ";
//            }
//            std::cout<<std::endl;
            
            eta = 0;
            tab[7][0] = 0.25 * (1 - ksi) * (1 - eta);
            tab[7][1] = 0.25 * (1 + ksi) * (1 - eta);
            tab[7][2] = 0.25 * (1 + ksi) * (1 + eta);
            tab[7][3] = 0.25 * (1 - ksi) * (1 + eta);
            
//            for(int i = 0;i<4;i++)
//            {
//                std::cout<<tab[7][i]<<" ";
//            }
//            std::cout<<std::endl;
//
            ksi = 0;
            tab[8][0] = 0.25 * (1 - ksi) * (1 - eta);
            tab[8][1] = 0.25 * (1 + ksi) * (1 - eta);
            tab[8][2] = 0.25 * (1 + ksi) * (1 + eta);
            tab[8][3] = 0.25 * (1 - ksi) * (1 + eta);
            
//            for(int i = 0;i<4;i++)
//            {
//                std::cout<<tab[8][i]<<" ";
//            }
//            std::cout<<std::endl;
        }

    }

    void calc_walls(Element* el)
    {
        //zmienne pomocnicze
        int L = 0, R = 0, U = 0, D = 0;
        for (int i = 0; i < mesh->nE; i++)
        {
            if ((el[i].nodes[1].bc == 1) && (el[i].nodes[2].bc == 1))
            {
                Rwall.w[R] = el[i];
                R++;
            }
            if ((el[i].nodes[3].bc == 1) && (el[i].nodes[0].bc == 1))
            {
                Lwall.w[L] = el[i];
                L++;
            }
            if ((el[i].nodes[0].bc == 1) && (el[i].nodes[1].bc == 1))
            {
                Dwall.w[D] = el[i];
                D++;
            }
            if ((el[i].nodes[2].bc == 1) && (el[i].nodes[3].bc == 1))
            {
                Uwall.w[U] = el[i];
                U++;
            }
        }
    }

//    void show_dksi_deta()
//    {
//
//        if (pc == 2)
//        {
//
//            std::cout << "\ndN/dKsi:\n";
//            for (int i = 0; i < 4; i++)
//            {
//                for (int j = 0; j < 4; j++)
//                {
//                    std::cout << dKsi[i][j] << "\t";
//                }
//                std::cout << "\n";
//            }
//
//            std::cout << "\ndN/dEta:\n";
//            for (int i = 0; i < 4; i++)
//            {
//                for (int j = 0; j < 4; j++)
//                {
//                    std::cout << dEta[i][j] << "\t";
//                }
//                std::cout << "\n";
//            }
//        }
//
//        else if (pc == 3) {
//
//            std::cout << "\ndN/dKsi:\n";
//            for (int i = 0; i < 9; i++)
//            {
//                for (int j = 0; j < 4; j++)
//                {
//                    std::cout << dKsi[i][j] << "\t";
//                }
//                std::cout << "\n";
//            }
//
//            std::cout << "\ndN/dEta:\n";
//            for (int i = 0; i < 9; i++)
//            {
//                for (int j = 0; j < 4; j++)
//                {
//                    std::cout << dEta[i][j] << "\t";
//                }
//                std::cout << "\n";
//            }
//        }
//    }


    double f1(double x)
    {
        return (1.0 / 4.0) * (1 - x);
    }

    double f2(double x)
    {
        return (1.0 / 4.0) * (1 + x);
    }

    double f3(double x)
    {
        return (1.0 / 4.0) * (x - 1);
    }

    void count_dksi_deta()
    {
        if (pc == 2)
        {
            double x = 1.0 / sqrt(3.0);

            double tabXY[4][2] = { {-x,-x},
                                    {x,-x},
                                    {x,x},
                                    {-x, x} };


            for (int i = 0; i < 4; i++)
            {
                dKsi[i][0] = f3(tabXY[i][0]);
                dKsi[i][1] = -f2(tabXY[i][0]);
                dKsi[i][2] = f2(tabXY[i][0]);
                dKsi[i][3] = f1(tabXY[i][0]);

                dEta[i][0] = f3(tabXY[i][1]);
                dEta[i][1] = f1(tabXY[i][1]);
                dEta[i][2] = f2(tabXY[i][1]);
                dEta[i][3] = -f2(tabXY[i][1]);
            }
        }
        else if (pc == 3)
        {
            double x = sqrt(3.0 / 5.0), x1 = 0.0;
            double tabXY[9][2] = { {-x,-x},
                                    {x1,-x},
                                    {x,-x},
                                    {x, x1},
                                    {x,x},
                                    {x1,x},
                                    {-x,x},
                                    {-x,x1},
                                    {x1,x1} };


            for (int i = 0; i < 9; i++)
            {
                dKsi[i][0] = f3(tabXY[i][0]);
                dKsi[i][1] = -f2(tabXY[i][0]);
                dKsi[i][2] = f2(tabXY[i][0]);
                dKsi[i][3] = f1(tabXY[i][0]);

                dEta[i][0] = f3(tabXY[i][1]);
                dEta[i][1] = f1(tabXY[i][1]);
                dEta[i][2] = f2(tabXY[i][1]);
                dEta[i][3] = -f2(tabXY[i][1]);
            }
        }
    }

    void count_Hbc()
    {
        for (int i = 0; i < nB - 1; i++)
        {
            calc_Hbc(Dwall.w[i], Dwall, 0);
            calc_Hbc(Uwall.w[i], Uwall, 2);
        }
        for (int i = 0; i < nH - 1; i++)
        {
            calc_Hbc(Lwall.w[i], Lwall, 3);
            calc_Hbc(Rwall.w[i], Rwall, 1);
        }
    }

    void calc_Hbc(Element el, wall wall, int sId)
    {
        double det = el.length(sId);
        if (pc == 2) {
            for (int j = 0; j < 4; j++)
            {
                for (int k = 0; k < 4; k++)
                {
                    el.Hbc[j][k] += alfa * det / 2 * wall.n1[j] * wall.n1[k];
//                    std::cout<<e.Hbc[j][k]<<" ";
                    el.Hbc[j][k] += alfa * det / 2 * wall.n2[j] * wall.n2[k];
//                    std::cout<<e.Hbc[j][k]<<" ";
                }
            }
        }
        else if (pc == 3)
        {
            double w1 = 5.0 / 9.0, w2 = 8.0 / 9.0;
            for (int j = 0; j < 4; j++)
            {
                for (int k = 0; k < 4; k++)
                {
                    el.Hbc[j][k] += alfa * det / 2 * w1 * wall.n1[j] * wall.n1[k];
                    el.Hbc[j][k] += alfa * det / 2 * w2 * wall.n2[j] * wall.n2[k];
                    el.Hbc[j][k] += alfa * det / 2 * w1 * wall.n3[j] * wall.n3[k];
                }
            }
        }
    }

   
    void count_vec_P()
    {
        for (int i = 0; i < nB - 1; i++)
        {
            Dwall.w[i].count_P(Dwall.n1, Dwall.n2, Dwall.n3, ambient_Temp, 0);
            Uwall.w[i].count_P(Uwall.n1, Uwall.n2, Uwall.n3, ambient_Temp, 2);
        }
        for (int i = 0; i < nH - 1; i++)
        {
            Lwall.w[i].count_P(Lwall.n1, Lwall.n2, Lwall.n3, ambient_Temp, 3);
            Rwall.w[i].count_P(Rwall.n1, Rwall.n2, Rwall.n3, ambient_Temp, 1);
        }
    }
    
    void agreg_vec_P()
    {

        count_vec_P();
        double* wektorTemp = new double[mesh->nN];
        for (int i = 0; i < mesh->nE; i++) {
            for (int i = 0; i < mesh->nN; i++)
            {
                wektorTemp[i] = 0;
            }
            mesh->elements[i].agreg_P(wektorTemp);
            mesh->add_To_Global_vector(wektorTemp);
        }
        //mesh->wyświetlGlobalnyP();
    }

    void sim_result()
    {
        double* result = solve_eq(mesh->global_H, mesh->final_P, mesh->nN);
        //std::cout << "Test case : \t" << *sol << "\n";
        double min = result[0], max = 0.0;
        for (int i = 0; i < mesh->nN; i++)
        {
            mesh->T0[i] = result[i];
            mesh->final_P[i] = 0;
            //std::cout << sol[i] << "\n";
            if (result[i] < min) {
                min = result[i];
            }
            else if (result[i] > max)
            {
                max = result[i];
            }
        }
        std::cout << "MaxTemp: " << max << " MinTemp: " << min << "\n";
        std::cout << "\n";

    }

    double* solve_eq(double** gh, double* ep, int num)
    {
        double* n = new double[num];
        double* x1 = new double[num];
        double* x2 = new double[num];
        double** M = new double* [num];

        for (int i = 0; i < num; i++)
        {
            x1[i] = 0.0;
            x2[i] = 0.0;
            n[i] = 0.0;
        }

        for (int i = 0; i < num; i++)
        {
            M[i] = new double[num];
            for (int j = 0; j < num; ++j)
            {
                M[i][j] = 0.0;
            }
        }

        for (int i = 0; i < num; i++)
        {
            n[i] = 1 / gh[i][i];
        }
        // Calculate M = -D^-1 (L + U)
        for (int i = 0; i < num; i++)
        {
            for (int j = 0; j < num; j++)
            {
                if (i == j)
                {
                    M[i][j] = 0.0;
                }
                else
                {
                    M[i][j] = -(gh[i][j] * n[i]);
                }
            }
        }

        for (int k = 0; k < 100; k++)
        {
            for (int i = 0; i < num; i++)
            {
                x2[i] = n[i] * ep[i];
                for (int j = 0; j < num; j++)
                {
                    x2[i] += M[i][j] * x1[j];
                }
            }
            for (int i = 0; i < num; i++) {
                x1[i] = x2[i];
            }
        }
        return x1;
    }
};


#endif /* element_Univ_h */
