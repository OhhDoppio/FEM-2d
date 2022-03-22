//
//  main.cpp
//  Mes Project
//
//  Created by Krystian Zapart
//

#include <iostream>
#include <algorithm>

#include "Define.h"

#include "Node.h"
#include "element.h"
#include "Mesh.h"
#include "wall.h"
#include "element_univ.h"


int main(int argc, const char * argv[])
{
             
        Mesh mesh(h,b,nH, nB,pc,initial_Temperature,sim_step,ambient_Temp,alfa,specific_heat,conductivity,ro);
       
//        mesh.show_Mesh();
    
        Element_Univ el(&mesh);

        el.count_dksi_deta();
        
//        el.show_dksi_deta();
    
        mesh.fJakobian(el.dKsi, el.dEta);
       
        el.count_Hbc();
        
        mesh.agreg_H();
    
//        mesh.show_H_global();
    
        std::cout<<"\n";
       
        el.agreg_vec_P();
    
        mesh.agreg_C();

        mesh.sum_H_C();
    

//    mesh.show_Elementy();
    
//    mesh.show_C();

//    mesh.show_H_global();
    
       for (int i = 0; i < sim_time / sim_step; i++)
       {
           std::cout << " "<< i+1 << "\n";
           mesh.count_final_P();
           el.sim_result();
       }
       
   
    
    return 0;
}
