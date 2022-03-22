//
//  wall.h
//  Mes Project
//
//  Created by Krystian Zapart
//

#ifndef wall_h
#define wall_h

#include <iostream>
#include <cmath>
#include "element.h"

struct wall {
    Element* w;
    double* n1;
    double* n2;
    double* n3;
    int size;

    wall() {};
    //pc == 2
    wall(int size, double ksi1, double eta1, double ksi2, double eta2)
    {
        this->size = size;
        w = new Element[size];
        n1 = new double[4];
        n2 = new double[4];
        n3 = new double[4];

        n1[0] = 0.25 * (1 - ksi1) * (1 - eta1);
        n1[1] = 0.25 * (1 + ksi1) * (1 - eta1);
        n1[2] = 0.25 * (1 + ksi1) * (1 + eta1);
        n1[3] = 0.25 * (1 - ksi1) * (1 + eta1);

        n2[0] = 0.25 * (1 - ksi2) * (1 - eta2);
        n2[1] = 0.25 * (1 + ksi2) * (1 - eta2);
        n2[2] = 0.25 * (1 + ksi2) * (1 + eta2);
        n2[3] = 0.25 * (1 - ksi2) * (1 + eta2);

        n3[0] = 0;
        n3[1] = 0;
        n3[2] = 0;
        n3[3] = 0;
    }

    //pc == 3
    wall(int size, double ksi1, double eta1, double ksi2, double eta2, double ksi3, double eta3)
    {
        this->size = size;
        w = new Element[size];
        n1 = new double[4];
        n2 = new double[4];
        n3 = new double[4];

        n1[0] = 0.25 * (1 - ksi1) * (1 - eta1);
        n1[1] = 0.25 * (1 + ksi1) * (1 - eta1);
        n1[2] = 0.25 * (1 + ksi1) * (1 + eta1);
        n1[3] = 0.25 * (1 - ksi1) * (1 + eta1);

        n2[0] = 0.25 * (1 - ksi2) * (1 - eta2);
        n2[1] = 0.25 * (1 + ksi2) * (1 - eta2);
        n2[2] = 0.25 * (1 + ksi2) * (1 + eta2);
        n2[3] = 0.25 * (1 - ksi2) * (1 + eta2);

        n3[0] = 0.25 * (1 - ksi3) * (1 - eta3);
        n3[1] = 0.25 * (1 + ksi3) * (1 - eta3);
        n3[2] = 0.25 * (1 + ksi3) * (1 + eta3);
        n3[3] = 0.25 * (1 - ksi3) * (1 + eta3);
    }
};

#endif /* wall_h */
