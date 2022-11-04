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

#include "./addons/keras/src/model.h"

#include <chrono>
using namespace std::chrono;
	
using namespace BioFVM;
using namespace PhysiCell;

int main( int argc, char* argv[] )
{
    auto start = high_resolution_clock::now();

    
    auto model = keras2cpp::Model::load("WT_DNN.model"); //model input
    
    std::cout << "test" << std::endl;
    keras2cpp::Tensor in{5}; //
    in.data_ = {0.16725,0.003,9.6,0.54,0.78};
    

    #pragma omp parallel for 
    for( int i=0; i < 3125; i++ )
    {
        keras2cpp::Tensor out = model(in);
        out.print();
    } 
    
        
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop - start);
 
    std::cout << "Time taken by function: " << duration.count()/1000000 << " seconds" << std::endl;
 
    return 0;
    
	// load and parse settings file(s)
	bool XML_status = false; 
	char copy_command [1024]; 
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
	
	
	bool whole_well =  parameters.bools( "whole_well_simulation" );
	bool intracellular_simulation = parameters.bools( "intracellular_simulation" );
	Microenvironment coarse_well;
	Microenvironment transfer_region;
	
	if (whole_well == true)
	{
		std::cout << "whole well simulation is active" << std::endl;
		coarse_well.name = "coarse_well";
		coarse_well.spatial_units = "micron";
		coarse_well.mesh.units = "micron";
		coarse_well.time_units = "min";
		
		coarse_well.set_density( 0 , "oxygen", "mM", 108000 , 0.00 );
		coarse_well.add_density( "glucose", "mM", 30000 , 0.0 );
		coarse_well.add_density( "chemokine", "mM", 40000 , 0.0);
		// coarse_well.resize_space( 100, 1 , 1 );
		
		double dx = 32;
		double dy = 2880;
		double dz = 2880;
		
		// coarse_well.resize_space( -dx/2.0+16 , dx/2.0+16, 256.0, 5104.0 , -dz/2.0+16 , dz/2.0+16 , dx, dy, dz );
		coarse_well.resize_space( 256.0, 5120.0, -1440.0, 1440.0, -1440.0, 1440.0, dx, dy, dz );
		std::vector<double> dirichlet_condition = { 0 , 0, 0 };

		coarse_well.set_substrate_dirichlet_activation(0,false);
		coarse_well.set_substrate_dirichlet_activation(1,false);
		coarse_well.set_substrate_dirichlet_activation(2,false);
		
		coarse_well.diffusion_decay_solver = diffusion_decay_solver__constant_coefficients_LOD_1D;   
		

		int coarse_well_voxel_number = coarse_well.mesh.voxels.size();
		
		for ( int m = 0; m < coarse_well_voxel_number ; m++)
		{
			coarse_well(m)[0]=17.5; // oxygen
			coarse_well(m)[1]=5.5; // glucose
			coarse_well(m)[2]=0; //chemokine
		}

		coarse_well.display_information( std::cout );
		coarse_well.write_to_matlab("output/output00000000_microenvironment1.mat");
		
		
		transfer_region.name = "transfer_region";
		transfer_region.spatial_units = "micron";
		transfer_region.mesh.units = "micron";
		transfer_region.time_units = "min";
		
		double tr_dx = 32;
		double tr_dy = 32;
		double tr_dz = 32;

		transfer_region.set_density( 0 , "glucose", "mmHg", 30000 , 0.00 );
		transfer_region.add_density( "glutamine", "mM", 30000 , 0.0 );
		transfer_region.add_density( "lactate", "mM", 30000 , 0.0);
		transfer_region.resize_space( 224.0, 288.0, -1440, 1440, -1440, 1440, tr_dx, tr_dy, tr_dz );
		
		transfer_region.diffusion_decay_solver = diffusion_decay_solver__constant_coefficients_LOD_3D;   
		
		for ( int m = 0; m < transfer_region.mesh.voxels.size() ; m++)
		{
			transfer_region(m)[0]=17.5; // glucose
			transfer_region(m)[1]=5.5; // glucose
			transfer_region(m)[2]=0; //chemokine
		}
		
		transfer_region.display_information( std::cout );
		transfer_region.write_to_matlab("output/output00000000_microenvironment2.mat");    
			
			
	}
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

    double DNN_intracellular_dt = 6.0;
    double last_intracellular_time  = 0.0; 
    double intracellular_dt_tolerance = 0.001 * DNN_intracellular_dt; 
    double next_intracellular_update = DNN_intracellular_dt; 
	/* Users typically stop modifying here. END USERMODS */ 
	
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
					
					if (whole_well == true)
					{
						sprintf( filename , "%s/output%08u_microenvironment1.mat" , PhysiCell_settings.folder.c_str(),  PhysiCell_globals.full_output_index );      
						coarse_well.write_to_matlab(filename);
						
						sprintf( filename , "%s/output%08u_microenvironment2.mat" , PhysiCell_settings.folder.c_str(),  PhysiCell_globals.full_output_index );      
						transfer_region.write_to_matlab(filename);
					}
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
			
			if (whole_well = true)
			{
				
				coarse_well.simulate_diffusion_decay(diffusion_dt);
            
            
				// Obtain coarse well concentrations
				std::vector<double> v1 = {0, 0, 0};
				std::vector<double> v2 = {0, 0, 0};
				for ( int m = 0; m < coarse_well.mesh.voxels.size() ; m++)
				{
					double mic_cen_x = coarse_well.mesh.voxels[m].center[0];
					if (mic_cen_x == 272)
					{ 
						v1[0]+=coarse_well(m)[0]; //oxygen
						v1[1]+=coarse_well(m)[1]; //glucose
						v1[2]+=coarse_well(m)[2]; //chemokine
					}
				}
            
				// Copy coarse well concentrations into "coarse" side of transfer region
				for ( int m = 0; m < transfer_region.mesh.voxels.size() ; m++)
				{
					double mic_cen_x = transfer_region.mesh.voxels[m].center[0];
					if (mic_cen_x == 272)
					{ 
						transfer_region(m)[0]=v1[0]; // oxygen
						transfer_region(m)[1]=v1[1]; // glucose
						transfer_region(m)[2]=v1[2]; // chemokine
					}
				}

				// Obtain fine microenvironment concentrations
				int tr_index = 0;
				int row_length = 90;
				int jump = row_length + 1;
				for ( int m = 0; m < microenvironment.mesh.voxels.size() ; m++)
				{  
					double mic_cen_x = microenvironment.mesh.voxels[m].center[0];
					if (mic_cen_x == 240)
					{
						transfer_region(tr_index)[0]=microenvironment(m)[0]; //oxygen
						transfer_region(tr_index)[1]=microenvironment(m)[1]; //glucose
						transfer_region(tr_index)[2]=microenvironment(m)[2]; //chemokine

						tr_index += 2; 
						// if ((tr_index != 0) && ((tr_index + 1)%row_length == 0))
						// {
						// 	tr_index += jump;
						// }
						// else
						// {
						// 	tr_index += 2;
						// }
					}
				}
				
				// std::vector<double> right_side_before_diffusion = {transfer_region(0)[0], transfer_region(0)[1], transfer_region(0)[2]};
				// std::vector<double> left_side_before_diffusion = {transfer_region(1)[0], transfer_region(1)[1], transfer_region(1)[2]};
				
				transfer_region.simulate_diffusion_decay(diffusion_dt);
				
				
				// std::vector<double> right_side_after_diffusion = {transfer_region(0)[0], transfer_region(0)[1], transfer_region(0)[2]};
				// std::vector<double> left_side_after_diffusion = {transfer_region(1)[0], transfer_region(1)[1], transfer_region(1)[2]};
				
				// left side overwrite
				// coarse_well(0)[0] += left_side_after_diffusion[0] - left_side_before_diffusion[0];
				// coarse_well(0)[1] += left_side_after_diffusion[1] - left_side_before_diffusion[1];
				// coarse_well(0)[2] += left_side_after_diffusion[2] - left_side_before_diffusion[2];
				// Dirichlet Boundary Condition
				// coarse_well(coarse_well_voxel_number-1)[0] = 0.285;
				
				// right side overwrite
				//std::cout << y_240 << std::endl;
				// double oxy_diff = right_side_after_diffusion[0] - right_side_before_diffusion[0];
				// double glu_diff = right_side_after_diffusion[1] - right_side_before_diffusion[1];
				// double chem_diff = right_side_after_diffusion[2] - right_side_before_diffusion[2];
				
				std::vector<double> v3 = {0, 0, 0};

				// Coarsen "coarse" side of transfer region
				for ( int m = 0; m < transfer_region.mesh.voxels.size() ; m++)
				{
					double mic_cen_x = transfer_region.mesh.voxels[m].center[0];
					if (mic_cen_x == 272)
					{ 
						v3[0] += transfer_region(m)[0]*transfer_region.mesh.voxels[m].volume;
						v3[1] += transfer_region(m)[1]*transfer_region.mesh.voxels[m].volume;
						v3[2] += transfer_region(m)[2]*transfer_region.mesh.voxels[m].volume;
						//std::cout << "Glucose difference per voxel  : "  <<glu_diff/y_240 << std::endl; 
						//microenvironment(m)[0] += oxy_diff/y_240; //oxygen
						//microenvironment(m)[1] += glu_diff/y_240; //glucose
						//microenvironment(m)[2] += chem_diff/y_240; //chemokine
					}
				}
				v3[0] /= coarse_well.mesh.voxels[0].volume;
				v3[1] /= coarse_well.mesh.voxels[0].volume;
				v3[2] /= coarse_well.mesh.voxels[0].volume;

				coarse_well(0)[0] = v3[0];
				coarse_well(0)[1] = v3[1];
				coarse_well(0)[2] = v3[2];
				
				// Write values for "fine" side of transfer region
				tr_index = 0;
				for ( int m = 0; m < microenvironment.mesh.voxels.size() ; m++)
				{  
					double mic_cen_x = microenvironment.mesh.voxels[m].center[0];
					if (mic_cen_x == 240)
					{ 
						microenvironment(m)[0] = transfer_region(tr_index)[0]; //oxygen
						microenvironment(m)[1] = transfer_region(tr_index)[1];
						microenvironment(m)[2] = transfer_region(tr_index)[2];

						tr_index += 2;
						// if ((tr_index != 0) && ((tr_index + 1)%row_length == 0))
						// {
						// 	tr_index += jump;
						// }
						// else
						// {
						// 	tr_index += 1;
						// }
					}	
				}
			}
			
			if ( intracellular_simulation == true)
			{		
				if( PhysiCell_globals.current_time >= next_intracellular_update )
				{
					//std::cout << "Entering intracellular" << std::endl;
					simulate_DNN(DNN_intracellular_dt);
					next_intracellular_update += DNN_intracellular_dt; 
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
