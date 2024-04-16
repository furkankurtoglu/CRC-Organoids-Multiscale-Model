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

#include "./core/PhysiCell.h"
#include "./modules/PhysiCell_standard_modules.h" 

// put custom code modules here! 
#include "./custom_modules/custom.h" 

// #include "./addons/keras/src/model.h"
	
using namespace BioFVM;
using namespace PhysiCell;

int main( int argc, char* argv[] )
{
	// load and parse settings file(s)
	bool XML_status = false; 
	char copy_command [1024]; 
    
    time_t t; // t passed as argument in function time()
    struct tm * tt; // decalring variable for localtime()
    
    
	if( argc > 1 )
	{
		XML_status = load_PhysiCell_config_file( argv[1] ); 
		sprintf( copy_command , "cp %s %s" , argv[1] , PhysiCell_settings.folder.c_str() ); 
	}
	else
	{
		XML_status = load_PhysiCell_config_file( "./config/PhysiCell_settings.xml" );
		sprintf( copy_command , "cp ./config/PhysiCell_settings.xml %s" , PhysiCell_settings.folder.c_str() ); 
	}
	if( !XML_status )
	{ exit(-1); }
	
	// copy config file to output directry 
	system( copy_command ); 
	
	// OpenMP setup
	omp_set_num_threads(PhysiCell_settings.omp_num_threads);
	
	// time setup 
	std::string time_units = "min"; 

	/* Microenvironment setup */ 
	
	setup_microenvironment(); // modify this in the custom code 
	std::cout << "setup microvironment is done" << std::endl;
	bool intracellular_simulation = parameters.bools( "intracellular_simulation" );
	bool whole_well =  parameters.bools( "whole_well_simulation" );
 

    // Create Coarse Microenvironment
    Microenvironment coarse_well;
    Microenvironment* PCoarse_well = &coarse_well;
    create_coarse_microenvironment(PCoarse_well);

    // Create Transfer Region
    Microenvironment transfer_region;
    Microenvironment* PTransfer = &transfer_region;
    create_transfer_region(PTransfer);
    			
	std::cout << "whole well is created" << std::endl;
	
	/* PhysiCell setup */ 
 	
	// set mechanics voxel size, and match the data structure to BioFVM
	double mechanics_voxel_size = 32; 
	Cell_Container* cell_container = create_cell_container_for_microenvironment( microenvironment, mechanics_voxel_size );
	
	/* Users typically start modifying here. START USERMODS */ 
	create_cell_types();
	std::cout << "create cell types is done" << std::endl;
	
	setup_tissue();

	std::cout << "setup tissue is done" << std::endl;

    double DNN_intracellular_dt = 0.01;
    double last_intracellular_time  = 0.0; 
    double intracellular_dt_tolerance = 0.001 * DNN_intracellular_dt; 
    double next_intracellular_update = DNN_intracellular_dt; 
	/* Users typically stop modifying here. */ 
	
	// set MultiCellDS save options 

	set_save_biofvm_mesh_as_matlab( true ); 
	set_save_biofvm_data_as_matlab( true ); 
	set_save_biofvm_cell_data( true ); 
	set_save_biofvm_cell_data_as_custom_matlab( true );
	
	// save a simulation snapshot 
	
	char filename[1024];
	sprintf( filename , "%s/initial" , PhysiCell_settings.folder.c_str() ); 
	save_PhysiCell_to_MultiCellDS_v2( filename , microenvironment , PhysiCell_globals.current_time ); 
	
	// save a quick SVG cross section through z = 0, after setting its 
	// length bar to 200 microns 

	PhysiCell_SVG_options.length_bar = 200; 

	// for simplicity, set a pathology coloring function 
	
	std::vector<std::string> (*cell_coloring_function)(Cell*) = my_coloring_function; 
	
	sprintf( filename , "%s/initial.svg" , PhysiCell_settings.folder.c_str() ); 
	SVG_plot( filename , microenvironment, 0.0 , PhysiCell_globals.current_time, cell_coloring_function );
	
	sprintf( filename , "%s/legend.svg" , PhysiCell_settings.folder.c_str() ); 
	create_plot_legend( filename , cell_coloring_function ); 
	
	display_citations(); 
	
	// set the performance timers 

	BioFVM::RUNTIME_TIC();
	BioFVM::TIC();
	
	std::ofstream report_file;
	if( PhysiCell_settings.enable_legacy_saves == true )
	{	
		sprintf( filename , "%s/simulation_report.txt" , PhysiCell_settings.folder.c_str() ); 
		
		report_file.open(filename); 	// create the data log file 
		report_file<<"simulated time\tnum cells\tnum division\tnum death\twall time"<<std::endl;
	}
	
	// main loop 
	
	try 
	{		
		while( PhysiCell_globals.current_time < PhysiCell_settings.max_time + 0.1*diffusion_dt )
		{
			// save data if it's time. 
			if( fabs( PhysiCell_globals.current_time - PhysiCell_globals.next_full_save_time ) < 0.01 * diffusion_dt )
			{
				display_simulation_status( std::cout ); 
				if( PhysiCell_settings.enable_legacy_saves == true )
				{	
					log_output( PhysiCell_globals.current_time , PhysiCell_globals.full_output_index, microenvironment, report_file);
				}
				
				if( PhysiCell_settings.enable_full_saves == true )
				{	
					sprintf( filename , "%s/output%08u" , PhysiCell_settings.folder.c_str(),  PhysiCell_globals.full_output_index ); 
					
					save_PhysiCell_to_MultiCellDS_v2( filename , microenvironment , PhysiCell_globals.current_time ); 

                    sprintf( filename , "%s/output%08u_microenvironment1.mat" , PhysiCell_settings.folder.c_str(),  PhysiCell_globals.full_output_index );      
                    coarse_well.write_to_matlab(filename);
                    
                    sprintf( filename , "%s/output%08u_microenvironment2.mat" , PhysiCell_settings.folder.c_str(),  PhysiCell_globals.full_output_index );      
                    transfer_region.write_to_matlab(filename);
                    

                time (&t); //passing argument to time()
                tt = localtime(&t);
                std::cout << "Current Day, Date and Time is = "<< asctime(tt) << std::endl;
				}
				
				PhysiCell_globals.full_output_index++; 
				PhysiCell_globals.next_full_save_time += PhysiCell_settings.full_save_interval;
			}
			
			// save SVG plot if it's time
			if( fabs( PhysiCell_globals.current_time - PhysiCell_globals.next_SVG_save_time  ) < 0.01 * diffusion_dt )
			{
				if( PhysiCell_settings.enable_SVG_saves == true )
				{	
					sprintf( filename , "%s/snapshot%08u.svg" , PhysiCell_settings.folder.c_str() , PhysiCell_globals.SVG_output_index ); 
					SVG_plot( filename , microenvironment, 0.0 , PhysiCell_globals.current_time, cell_coloring_function );
					
					PhysiCell_globals.SVG_output_index++; 
					PhysiCell_globals.next_SVG_save_time  += PhysiCell_settings.SVG_save_interval;
				}
			}

			// update the microenvironment
			microenvironment.simulate_diffusion_decay( diffusion_dt );

            coarse_well.simulate_diffusion_decay(diffusion_dt);
 
            //simulate_non_regular_mesh(Microenvironment* Pcoarse_well, Microenvironment* PTransfer);
        
             //std::cout << "Coarse Well diffusion is done" << std::endl;
            // Obtain coarse well concentrations
            std::vector<double> v1 = {0, 0, 0}; 
            std::vector<double> v2 = {0, 0, 0};
            for ( int m = 0; m < coarse_well.mesh.voxels.size() ; m++)
            {
                // Coarse microenvironment is 1D domain in Y and Z directions. There are multiple voxels in X directions.
                double mic_cen_x = coarse_well.mesh.voxels[m].center[0]; // Matching one big bread slice in a X-axis
                if (abs(mic_cen_x - 528) < 0.1 * 32 ) // 32 = coarse_well.mesh.dx???? 
                { 
                    v1[0]=coarse_well(m)[0]; //glucose    
                    v1[1]=coarse_well(m)[1]; //glutmaine
                    v1[2]=coarse_well(m)[2]; //lactate
                }
                //std::cout << "Obtaining Coarse well Concentrations" 
                //std::cout << v1 << std::endl;
            }
        
            // Copy coarse well concentrations into "coarse" side of transfer region
            for ( int m = 0; m < transfer_region.mesh.voxels.size() ; m++)
            {
                // Transfer region is 3D domain. The top voxels in X dimension corresponds the lowest X voxel in coarse microenvironment. 
                
                double mic_cen_x = transfer_region.mesh.voxels[m].center[0];
                if (abs(mic_cen_x - 528) < 0.1 * 32 )
                { 
                    transfer_region(m)[0]=v1[0]; // glucose
                    transfer_region(m)[1]=v1[1]; // glutamine
                    transfer_region(m)[2]=v1[2]; // lactate
                }
                //std::cout << v1 << std::endl;
            }

            // Obtain fine microenvironment concentrations
            for ( int m = 0; m < microenvironment.mesh.voxels.size() ; m++)
            {  
                double mic_cen_x = microenvironment.mesh.voxels[m].center[0];
                
                if (abs(mic_cen_x - 496) < 0.1 * 32 )
                {
                    int tr_id = transfer_region.nearest_voxel_index(microenvironment.mesh.voxels[m].center);
                    //std::cout << microenvironment.mesh.voxels[m].center << std::endl;
                    transfer_region(tr_id)[0]=microenvironment(m)[0]; //glucose
                    transfer_region(tr_id)[1]=microenvironment(m)[1]; //glutamine
                    transfer_region(tr_id)[2]=microenvironment(m)[2]; //lactate

                }
            }
            
            // Simulate Transfer Region
            transfer_region.simulate_diffusion_decay(diffusion_dt);
            
            // Updating Coarse and Fine Microenvironments based on the transfer region
            std::vector<double> v3 = {0, 0, 0};
            double voxel_counter_for_transfer_region = 0;
            
            // Coarsen "coarse" side of transfer region
            for ( int m = 0; m < transfer_region.mesh.voxels.size() ; m++)
            {
                double mic_cen_x = transfer_region.mesh.voxels[m].center[0];
                if (abs(mic_cen_x - 528) < 0.1 * 32 ) 
                { 
                    v3[0] += transfer_region(m)[0]; // Glucose
                    v3[1] += transfer_region(m)[1]; // Glutamine
                    v3[2] += transfer_region(m)[2]; // Lactate
                    voxel_counter_for_transfer_region = voxel_counter_for_transfer_region + 1.0;
                }
            }
            v3[0] = v3[0]/voxel_counter_for_transfer_region;
            v3[1] = v3[1]/voxel_counter_for_transfer_region;
            v3[2] = v3[2]/voxel_counter_for_transfer_region;
            
            std::vector<double> Nearest_location = {528,0,0};
            int idx = coarse_well.nearest_voxel_index(Nearest_location);
            
            coarse_well(idx)[0] = v3[0];
            coarse_well(idx)[1] = v3[1];
            coarse_well(idx)[2] = v3[2];
            
            // Write values for "fine" side of transfer region
            for ( int m = 0; m < microenvironment.mesh.voxels.size() ; m++)
            {  
                double mic_cen_x = microenvironment.mesh.voxels[m].center[0];
                if (abs(mic_cen_x - 496) < 0.1 * 32 ) 
                { 
                    int fine_id = transfer_region.nearest_voxel_index(microenvironment.mesh.voxels[m].center);
                    microenvironment(m)[0] = transfer_region(fine_id)[0]; //glucose
                    microenvironment(m)[1] = transfer_region(fine_id)[1]; //glutamine
                    microenvironment(m)[2] = transfer_region(fine_id)[2]; //lactate
                }
            }
			
			if ( intracellular_simulation == true)
			{		
				if( PhysiCell_globals.current_time >= next_intracellular_update )
				{
					//std::cout << "Entering intracellular" << std::endl;
					simulate_DNN(DNN_intracellular_dt);
					next_intracellular_update += DNN_intracellular_dt; 
                    //std::cout << "Exiting intracellular" << std::endl;
				}
			}
            // run PhysiCell 
			((Cell_Container *)microenvironment.agent_container)->update_all_cells( PhysiCell_globals.current_time );
			
			PhysiCell_globals.current_time += diffusion_dt;
		}
		
		if( PhysiCell_settings.enable_legacy_saves == true )
		{			
			log_output(PhysiCell_globals.current_time, PhysiCell_globals.full_output_index, microenvironment, report_file);
			report_file.close();
		}
	}
	catch( const std::exception& e )
	{ // reference to the base of a polymorphic object
		std::cout << e.what(); // information from length_error printed
	}
	
	// save a final simulation snapshot 
	
	sprintf( filename , "%s/final" , PhysiCell_settings.folder.c_str() ); 
	save_PhysiCell_to_MultiCellDS_v2( filename , microenvironment , PhysiCell_globals.current_time ); 
	
	sprintf( filename , "%s/final.svg" , PhysiCell_settings.folder.c_str() ); 
	SVG_plot( filename , microenvironment, 0.0 , PhysiCell_globals.current_time, cell_coloring_function );
	
	// timer 
	
	std::cout << std::endl << "Total simulation runtime: " << std::endl; 
	BioFVM::display_stopwatch_value( std::cout , BioFVM::runtime_stopwatch_value() ); 

	return 0; 
}