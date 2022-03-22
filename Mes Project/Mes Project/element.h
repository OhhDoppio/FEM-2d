//
//  element.h
//  Mes Project
//
//  Created by Krystian Zapart
//

#ifndef element_h
#define element_h

#include <iostream>
#include <cmath>

#include "Node.h"

struct Element
{
    int value;

    double** H;
    double** Hbc;
    double** C;
    double* P;
    Node* nodes;
    
    Element() {}
    
    Element(int value,int ID, int c, Node* arr, int nn)
    {
        this->value = value;
        
    
        nodes = new Node[4];
        nodes[0] = arr[ID - 1];
        nodes[1] = arr[ID + c - 1];
        nodes[2] = arr[nodes[1].id];
        nodes[3] = arr[ID];

        H = new double* [4];
        Hbc = new double* [4];
        C = new double* [4];
        P = new double[4];
        
        //Ustawienie wartości początkowych dla y
        for (int i = 0; i < 4; i++)
        {
            H[i] = new double[4];
            Hbc[i] = new double[4];
            C[i] = new double[4];
            P[i] = 0;
            for (int j = 0; j < 4; j++)
            {
                H[i][j] = 0;
                Hbc[i][j] = 0;
                C[i][j] = 0;
            }
        }

   
    }
    
//    void show_Element()
//    {
//        std::cout << "value: " << value << "\n";
//        std::cout << "nodes: " << nodes[0].id << ", " << nodes[1].id << ", " << nodes[2].id << ", " << nodes[3].id << "\n";
//    }
//
//    void show_H()
//    {
//        for (int i = 0; i < 4; i++)
//        {
//            for (int j = 0; j < 4; j++)
//            {
//                std::cout << H[i][j] << "\t";
//            }
//            std::cout << "\n";
//        }
//        std::cout << "\n";
//    }
//
//    void show_HBC()
//    {
//        for (int i = 0; i < 4; i++)
//        {
//            for (int j = 0; j < 4; j++)
//            {
//                std::cout << Hbc[i][j] << "\t";
//            }
//            std::cout << "\n";
//        }
//        std::cout << "\n";
//    }
     
    void count_C(int i, double** n, double** J, double waga)
    {
        double det = J[0][0] * J[1][1] - J[0][1] * J[1][0];

        if (pc == 2)
        {
            for (int j = 0; j < 4; j++)
            {
                for (int k = 0; k < 4; k++)
                {
                    C[j][k] += specific_heat * ro * det * (n[i][j] * n[i][k]);
                    //std::cout << C[j][k] << "\n";
                }
            }
        }
        else if (pc == 3)
        {
            for (int j = 0; j < 4; j++)
            {
                for (int k = 0; k < 4; k++)
                {
                    C[j][k] += specific_heat * ro * det * (n[i][j] * n[i][k]) * waga;
                    //std::cout << C[j][k] << "\n";
                }
            }
        }
    }
    
    double length(int sId)
    {
        int x1 = 0, x2 = 0;
        
        if (sId == 0)
        {
            x1 = 0;
            x2 = 1;
        }
        
        if (sId == 1)
        {
            x1 = 1;
            x2 = 2;
        }
        
        if (sId == 2)
        {
            x1 = 2;
            x2 = 3;
        }

        if (sId == 3)
        {
            x1 = 3;
            x2 = 0;
        }
        
        double xDiff = nodes[x1].x - nodes[x2].x;
        double yDiff = nodes[x1].y - nodes[x2].y;
        return sqrt(xDiff * xDiff + yDiff * yDiff);
    }

    void agreg(double** tab)
    {
        for (int i = 0; i < 4; i++)
        {
            for (int j = 0; j < 4; j++)
            {
                Hbc[i][j] += H[i][j];
            }
        }
        
        for (int i = 0; i < 4; i++)
        {
            for (int j = 0; j < 4; j++)
            {
                int x = nodes[i].id - 1;
                int y = nodes[j].id - 1;
                tab[y][x] = Hbc[i][j];
            }
        }
        
    }

    void agreg_P(double* vec)
    {
        for (int i = 0; i < 4; i++)
        {
            vec[nodes[i].id - 1] = P[i];
        }
    }

    void count_P(double* n1, double* n2, double* n3, double t, int sId)
    {
        double det = length(sId) / 2;
        if (pc == 2)
        {
            double w1 = 1, w2 = 1;
            for (int i = 0; i < 4; i++)
            {
                P[i] += 300.0 * ((w1 * n1[i]) * t + (w2 * n2[i]) * t) * det;
            }
        }
        else if (pc == 3)
        {
            double w1 = 5.0 / 9.0, w2 = 8.0 / 9.0;
            for (int i = 0; i < 4; i++)
            {
                P[i] += 300.0 * ((w1 * n1[i]) * t + (w2 * n2[i]) * t + (w1 * n3[i]) * t) * det;
            }
        }
    }

    void agreg_C(double** tab)
    {
        for (int i = 0; i < 4; i++)
        {
            for (int j = 0; j < 4; j++)
            {
                int x = nodes[i].id - 1;
                int y = nodes[j].id - 1;
                tab[y][x] = C[i][j];
            }
        }
    }
};
#endif /* element_h */
