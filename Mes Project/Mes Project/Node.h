//
//  Node.h
//  Mes Project
//
//  Created by Krystian Zapart
//

#ifndef Node_h
#define Node_h

#include <iostream>
#include <cmath>

#include "Define.h"

struct Node
{
    int id;
    double x, y;
    int bc = 0;
 
    Node()
    {
        this->id = 0;
        this->x = 0;
        this->y = 0;
        this->bc = 0;
    }

    Node(int ID, double X, double Y)
    {
        id = ID;
        x = X;
        y = Y;

        if (x == 0) {
            bc = 1;
        }
        else if (y == 0) {
            bc = 1;
        }
        else if (x == b) {
            bc = 1;
        }
        else if (y == h) {
            bc = 1;
        }
    }

//    void showNode()
//    {
//        std::cout << "[" << x << "," << y << "]\t"<<"BC: "<<bc<<"\t";
//    }
//
//    void showNodeNumer()
//    {
//        std::cout << "Node num: " << id << "[" << x << "," << y << "]\n";
//    }

};
#endif /* Node_h */
