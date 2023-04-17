/*
###############################################################################
# If you use PhysiCell in your project, please cite PhysiCell and the version #
# number, such as below:                                                      #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1].    #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# See VERSION.txt or call get_PhysiCell_version() to get the current version  #
#     x.y.z. Call display_citations() to get detailed information on all cite-#
#     able software used in your PhysiCell application.                       #
#                                                                             #
# Because PhysiCell extensively uses BioFVM, we suggest you also cite BioFVM  #
#     as below:                                                               #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1],    #
# with BioFVM [2] to solve the transport equations.                           #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# [2] A Ghaffarizadeh, SH Friedman, and P Macklin, BioFVM: an efficient para- #
#     llelized diffusive transport solver for 3-D biological simulations,     #
#     Bioinformatics 32(8): 1256-8, 2016. DOI: 10.1093/bioinformatics/btv730  #
#                                                                             #
###############################################################################
#                                                                             #
# BSD 3-Clause License (see https://opensource.org/licenses/BSD-3-Clause)     #
#                                                                             #
# Copyright (c) 2015-2022, Paul Macklin and the PhysiCell Project             #
# All rights reserved.                                                        #
#                                                                             #
# Redistribution and use in source and binary forms, with or without          #
# modification, are permitted provided that the following conditions are met: #
#                                                                             #
# 1. Redistributions of source code must retain the above copyright notice,   #
# this list of conditions and the following disclaimer.                       #
#                                                                             #
# 2. Redistributions in binary form must reproduce the above copyright        #
# notice, this list of conditions and the following disclaimer in the         #
# documentation and/or other materials provided with the distribution.        #
#                                                                             #
# 3. Neither the name of the copyright holder nor the names of its            #
# contributors may be used to endorse or promote products derived from this   #
# software without specific prior written permission.                         #
#                                                                             #
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" #
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE   #
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE  #
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE   #
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR         #
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF        #
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS    #
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN     #
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)     #
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE  #
# POSSIBILITY OF SUCH DAMAGE.                                                 #
#                                                                             #
###############################################################################
*/

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <ctime>
#include <cmath>
#include <omp.h>
#include <fstream>


// put custom code modules here! 
#include "./addons/keras/src/model.h"

#include <chrono>
using namespace std::chrono;


int main( int argc, char* argv[] )
{
    auto start = high_resolution_clock::now();

    
    auto model = keras2cpp::Model::load("WT_DNN_NO_INTRACELLULAR.model"); //model input
    
/*     double glc_ub = 0.200;
    double gln_ub = 0.003;
    double lac_i_ub = 14.4;
    double glc_i_ub = 1.08;
    double gln_i_ub = 0.39;
 */
    float glc_ub = 0.200;
    float gln_ub = 0.003;
    float lac_i_ub = 14.4;
    float glc_i_ub = 1.08;
    float gln_i_ub = 0.39;
    
    std::vector <double> glc_bds;
    std::vector <double> gln_bds;
    std::vector <double> glc_i_bds;
    std::vector <double> gln_i_bds;
    std::vector <double> lac_i_bds;
    
 
    keras2cpp::Tensor in{2};
    keras2cpp::Tensor out;
    std::vector<double> result;
    in.data_ = {glc_ub,gln_ub};
    out = model(in);
    out.print();
    
    result = out.result_vector();
    double biomass_creation_flux = result[0]/100;
    std::cout << "Biomass = " << biomass_creation_flux << std::endl;
    
    /* double number_of_evaluations = 5;
    
    for (double i = 0; i < number_of_evaluations+1; ++i)
        glc_bds.push_back(i/number_of_evaluations * glc_ub);
    
    for (double i = 0; i < number_of_evaluations+1; ++i)
        gln_bds.push_back(i/number_of_evaluations * gln_ub);
    
    for (double i = 0; i < number_of_evaluations+1; ++i)
        glc_i_bds.push_back(i/number_of_evaluations * glc_i_ub);

    for (double i = 0; i < number_of_evaluations+1; ++i)
        gln_i_bds.push_back(i/number_of_evaluations * gln_i_ub);
    
    for (double i = 0; i < number_of_evaluations+1; ++i)
        lac_i_bds.push_back(i/number_of_evaluations * lac_i_ub);    


    keras2cpp::Tensor in{5};
    keras2cpp::Tensor out;
    std::vector<double> result;

    std::ofstream i_myfile;
    std::ofstream o_myfile;
    o_myfile.open ("outputs.csv");
    i_myfile.open ("inputs.csv");

//    #pragma omp parallel for 
    for (int i = 0; i < glc_bds.size(); i++)
    {
        for (int j = 0; j < gln_bds.size(); j++)
        {
            for (int k = 0; k < glc_i_bds.size(); k++)
            {
                for (int l = 0; l < gln_i_bds.size(); l++)
                {
                    for (int m = 0; m < lac_i_bds.size(); m++)
                    {
                        float glc_bd = glc_bds[i];
                        float gln_bd = gln_bds[j];
                        float lac_i_bd = lac_i_bds[m];
                        float glc_i_bd = glc_i_bds[k];
                        float gln_i_bd = gln_i_bds[l];
                        in.data_ = {glc_bd,gln_bd,lac_i_bd,glc_i_bd,gln_i_bd};
                        out = model(in);
                        //out.print();
                        result = out.result_vector();
                        i_myfile << glc_bd << "," << gln_bd << "," << lac_i_bd << "," << glc_i_bd << "," << gln_i_bd << '\n';
                        o_myfile << result[0] << "," << result[1] << "," << result[2] << "," << result[3] <<'\n';
                    }
                }
            }
        }
    }

    o_myfile.close();
    i_myfile.close(); */


    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop - start);
 
    std::cout << "Time taken by function: " << duration.count()/1000 << " miliseconds" << std::endl;

    return 0;
    

}
