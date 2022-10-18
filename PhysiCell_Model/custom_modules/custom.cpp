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
# Copyright (c) 2015-2021, Paul Macklin and the PhysiCell Project             #
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

#include "../addons/keras/src/model.h"
#include "./custom.h"
#include "../modules/PhysiCell_settings.h"

// assume these files; can override in read_DNN()
auto WT_Model = keras2cpp::Model::load("WT_DNN.model");
auto KRAS_Model = keras2cpp::Model::load("KRAS_DNN.model");

void create_cell_types( void )
{
    read_DNN("WT_DNN.model", "KRAS_DNN.model");

	// set the random seed 
	SeedRandom( parameters.ints("random_seed") );  
	
	/* 
	   Put any modifications to default cell definition here if you 
	   want to have "inherited" by other cell types. 
	   
	   This is a good place to set default functions. 
	*/ 
	
	initialize_default_cell_definition(); 
	cell_defaults.phenotype.secretion.sync_to_microenvironment( &microenvironment ); 
	
	cell_defaults.functions.volume_update_function = standard_volume_update_function;
	cell_defaults.functions.update_velocity = standard_update_cell_velocity;

	cell_defaults.functions.update_migration_bias = NULL; 
	cell_defaults.functions.update_phenotype = NULL; // update_cell_and_death_parameters_O2_based; 
	cell_defaults.functions.custom_cell_rule = NULL; 
	cell_defaults.functions.contact_function = NULL; 
	
	cell_defaults.functions.add_cell_basement_membrane_interactions = NULL; 
	cell_defaults.functions.calculate_distance_to_membrane = NULL; 
	
	/*
	   This parses the cell definitions in the XML config file. 
	*/
	
	initialize_cell_definitions_from_pugixml(); 

	/*
	   This builds the map of cell definitions and summarizes the setup. 
	*/
		
	build_cell_definitions_maps(); 

	/*
	   This intializes cell signal and response dictionaries 
	*/

	setup_signal_behavior_dictionaries(); 	

	/* 
	   Put any modifications to individual cell definitions here. 
	   
	   This is a good place to set custom functions. 
	*/ 
	
	cell_defaults.functions.update_phenotype = phenotype_function; 
	cell_defaults.functions.custom_cell_rule = custom_function; 
	cell_defaults.functions.contact_function = contact_function; 
	
	/*
	   This builds the map of cell definitions and summarizes the setup. 
	*/
		
	display_cell_definitions( std::cout ); 
	
	return; 
}

void setup_microenvironment( void )
{
	// set domain parameters 
	
	// put any custom code to set non-homogeneous initial conditions or 
	// extra Dirichlet nodes here. 
	
	// initialize BioFVM 
	
	initialize_microenvironment(); 	
	
	return; 
}

void setup_tissue( void )
{
    static int glucose_substrate_index = microenvironment.find_density_index( "glucose" );
    static int glutamine_substrate_index = microenvironment.find_density_index( "glutamine" ); 
    static int lactate_substrate_index = microenvironment.find_density_index( "lactate");

	// place a cluster of tumor cells at the center 
 
	double tumor_radius = parameters.doubles( "tumor_radius" ); // 250.0; 
	
	int initial_tumor_radius = parameters.doubles( "tumor_radius" ); 
	
    bool Two_Dim_MicEnv =  parameters.bools( "two_dim_seeding" );
    
    std::cout << Two_Dim_MicEnv << std::endl;
    
    
    Cell* pCell; 
	
	

    if ( Two_Dim_MicEnv == true )
    {
        // 2D Cell Seeding
        std::cout << "Seeding Cell as 2D" << std::endl;
        double cell_radius = cell_defaults.phenotype.geometry.radius; 
        double cell_spacing = 0.8 * 2.0 * cell_radius; 
        
        std::vector<std::vector<double>> positions = create_cell_circle_positions(cell_radius,initial_tumor_radius);

        std::cout << "Creating cells" << std::endl;
        
        for( int i=0; i < positions.size(); i++ )
        {
            pCell = create_cell(get_cell_definition("CRC_WT"));
            pCell->functions.volume_update_function = NULL;
            pCell->assign_position( positions[i] );
             
            pCell->phenotype.molecular.internalized_total_substrates[glucose_substrate_index] = 1.56;//(1.56+0.88)/2; // FURKAN : Tomorrow I will make these ones stochastic
            pCell->phenotype.molecular.internalized_total_substrates[glutamine_substrate_index] = 1.08;//(1.08+0.56)/2; // FURKAN : Tomorrow I will make these ones stochastic
            pCell->phenotype.molecular.internalized_total_substrates[lactate_substrate_index] = 19.2;//(19.2+6.4)/2; // FURKAN : Tomorrow I will make these ones stochastic
            
            // Stochastic Volume
            if (parameters.bools("random_initial_volume"))
            {
                double a = uniform_random();
                if ( a > 0.6)
                {
                    pCell->phenotype.volume.multiply_by_ratio(a);
                }
            }
        }
    }
    else
    {
		double xmin=1.e6;
		double ymin=1.e6;
		double zmin=1.e6;
		double xmax= -xmin;
		double ymax= -ymin;
		double zmax= -zmin;
		
		double cell_radius = cell_defaults.phenotype.geometry.radius; 
		// double initial_tumor_radius = 46; // parameters.doubles("initial_tumor_radius");
		double initial_tumor_radius = parameters.doubles("initial_tumor_radius");
		std::cout << "------ initial_tumor_radius = " << initial_tumor_radius << std::endl;

		//rwh
		// double number_of_organoid = 250; //parameters.doubles("number_of_organoid")
		// int number_of_organoid = 250; // parameters.doubles("number_of_organoid")
		int number_of_organoids = parameters.ints("number_of_organoids");
		std::cout << "------ number_of_organoids = " << number_of_organoids << std::endl;		
		
        for (int i = 0; i < number_of_organoids; i++) // seeding number of organoid cells specified in PhysiCell_settings.xml
			    {
                    std::vector<std::vector<double>> positions = create_cell_sphere_positions(cell_radius,initial_tumor_radius); 
                    //std::cout << "creating " << positions.size() << " closely-packed organoid cells ... " << std::endl;
                    // create organoid
                        double xrand = (rand() % 5333) - 2666;
                        double yrand = (rand() % 961) - 480;
                        double zrand = (rand() % 5333) - 2666;
                        if (xrand < xmin) xmin = xrand;
                        if (xrand > xmax) xmax = xrand;
                        if (yrand < ymin) ymin = yrand;
                        if (yrand > ymax) ymax = yrand;
                        if (zrand < zmin) zmin = zrand;
                        if (zrand > zmax) zmax = zrand;
                    //std::cout << positions.size() << std::endl;
                    for( int i=0; i < positions.size(); i++ )
                    {
                        positions[i][0] += xrand;//(rand() % 5333) - 2666;
                        positions[i][1] += yrand;//(rand() % 961) - 480;
                        positions[i][2] += zrand;//(rand() % 5333) - 2666;
                        // pCell = create_cell(KRAS_negative);
                        pCell = create_cell( get_cell_definition("CRC_WT") );
                        pCell->assign_position( positions[i] );
						
						pCell->phenotype.molecular.internalized_total_substrates[glucose_substrate_index] = 1.56;//(1.56+0.88)/2; // FURKAN : Tomorrow I will make these ones stochastic
						pCell->phenotype.molecular.internalized_total_substrates[glutamine_substrate_index] = 1.08;//(1.08+0.56)/2; // FURKAN : Tomorrow I will make these ones stochastic
						pCell->phenotype.molecular.internalized_total_substrates[lactate_substrate_index] = 19.2;//(19.2+6.4)/2; // FURKAN : Tomorrow I will make these ones stochastic
                    }
			    }
    }
    
	
	return; 
}


std::vector<std::vector<double>> create_cell_circle_positions(double cell_radius, double sphere_radius)
{
	std::vector<std::vector<double>> cells;
	int xc=0,yc=0,zc=0;
	double x_spacing= cell_radius*sqrt(3);
	double y_spacing= cell_radius*sqrt(3);

	std::vector<double> tempPoint(3,0.0);
	
	for(double x=-sphere_radius;x<sphere_radius;x+=x_spacing, xc++)
	{
		for(double y=-sphere_radius;y<sphere_radius;y+=y_spacing, yc++)
		{
			tempPoint[1]=y + (xc%2) * cell_radius;
			tempPoint[0]=x;
			tempPoint[2]=0;
			if(sqrt(norm_squared(tempPoint))< sphere_radius)
			{ cells.push_back(tempPoint); }
		}
	}
	return cells;
}


std::vector<std::vector<double>> create_cell_sphere_positions(double cell_radius, double sphere_radius)
{
	std::vector<std::vector<double>> cells;
	int xc=0,yc=0,zc=0;
	double x_spacing= cell_radius*sqrt(3);
	double y_spacing= cell_radius*2;
	double z_spacing= cell_radius*sqrt(3);
	
	std::vector<double> tempPoint(3,0.0);
	
	for(double z=-sphere_radius;z<sphere_radius;z+=z_spacing, zc++)
	{
		for(double x=-sphere_radius;x<sphere_radius;x+=x_spacing, xc++)
		{
			for(double y=-sphere_radius;y<sphere_radius;y+=y_spacing, yc++)
			{
				tempPoint[0]=x + (zc%2) * 0.5 * cell_radius;
				tempPoint[1]=y + (xc%2) * cell_radius;
				tempPoint[2]=z;
				
				if(sqrt(norm_squared(tempPoint))< sphere_radius)
				{ cells.push_back(tempPoint); }
			}
			
		}
	}
	return cells;
}

std::vector<std::string> my_coloring_function( Cell* pCell )
{ return paint_by_number_cell_coloring(pCell); }

void phenotype_function( Cell* pCell, Phenotype& phenotype, double dt )
{ 
	return; 
}

void custom_function( Cell* pCell, Phenotype& phenotype , double dt )
{  
    return;
} 

void contact_function( Cell* pMe, Phenotype& phenoMe , Cell* pOther, Phenotype& phenoOther , double dt )
{ return; } 

void read_DNN(std::string wt_filename, std::string kras_filename)
{
    WT_Model = keras2cpp::Model::load(wt_filename.c_str());
    KRAS_Model = keras2cpp::Model::load(kras_filename.c_str());
}

void simulate_DNN(double intracellular_dt )
{

    static int glc_index = microenvironment.find_density_index( "glucose" );
	static int gln_index = microenvironment.find_density_index( "glutamine" );
	static int lac_index = microenvironment.find_density_index( "lactate" );
	
    static double exp_ave_n_cells = parameters.doubles("experimental_average_number_of_cells" );
    static double exp_vol_well = parameters.doubles("experimental_well_volume" );

    
    #pragma omp parallel for 
    for( int i=0; i < (*all_cells).size(); i++ )
    {
        // Wild type simulation
        if ((*all_cells)[i]->type == 0)
        {
            static int i_Glc_i = (*all_cells)[i]->custom_data.find_variable_index( "int_glc" );
            static int i_Gln_i = (*all_cells)[i]->custom_data.find_variable_index( "int_gln" );
            static int i_Lac_i = (*all_cells)[i]->custom_data.find_variable_index( "int_lac" );
            
            keras2cpp::Tensor in{5};
            keras2cpp::Tensor out;
            double glc_val_int = (*all_cells)[i]->nearest_density_vector()[glc_index];
            double gln_val_int = (*all_cells)[i]->nearest_density_vector()[gln_index];
            
            
            double u_glc = (*all_cells)[i]->custom_data[1] * exp_ave_n_cells / exp_vol_well * glc_val_int;
            double u_gln = (*all_cells)[i]->custom_data[2] * exp_ave_n_cells / exp_vol_well * gln_val_int;
            
            float fl_glc = u_glc;
            float fl_gln = u_gln;
            
            //std::cout << "Glucose = " << fl_glc << std::endl;
            //std::cout << "Glutamine = " << fl_gln << std::endl;    
			
			
			double int_conc_glu = (*all_cells)[i]->phenotype.molecular.internalized_total_substrates[glc_index];
            double int_conc_gln = (*all_cells)[i]->phenotype.molecular.internalized_total_substrates[gln_index];
            double int_conc_lac = (*all_cells)[i]->phenotype.molecular.internalized_total_substrates[lac_index];
			
			float fl_int_conc_glu = (*all_cells)[i]->custom_data[i_Glc_i];
            float fl_int_conc_gln = (*all_cells)[i]->custom_data[i_Gln_i];
			float fl_int_conc_lac = (*all_cells)[i]->custom_data[i_Lac_i];
			
			/* std::cout << "Input Intracellular Glucose = " << fl_int_conc_glu << std::endl;
			std::cout << "Input Intracellular Glutamine = " << fl_int_conc_gln << std::endl;    
			std::cout << "Input Intracellular Lactate = " << fl_int_conc_lac << std::endl;     */
            
			
			
            in.data_ = {fl_glc,fl_gln,fl_int_conc_lac,fl_int_conc_glu,fl_int_conc_gln};
			//in.data_ = {0.1115,0.0025,9.6,0.54,0.39}; test values
            out = WT_Model(in); // model evaluation
            //out.print();
            
			std::vector<double> result;
            result = out.result_vector();
            
            double biomass_creation_flux = result[0]/parameters.doubles("DNN_biomass_normalizer");
            
/*             std::cout << "Biomass = " << biomass_creation_flux << std::endl;
            std::cout << "Intracellular Glucose Change = " << result[1]/10 << std::endl;
			std::cout << "Intracellular Glutamine Change = " << result[2]/10 << std::endl;
			std::cout << "Intracellular Lactate Change = " << result[3]/10 << std::endl; */
            
            //(*all_cells)[i]->custom_data[biomass_vi]  = biomass_creation_flux;
            
            double volume_increase_ratio = 1 + ( biomass_creation_flux / 60 * intracellular_dt);
            (*all_cells)[i]->custom_data[0]  = biomass_creation_flux; // FURKAN to Fix = Manually written indices for custom data - USE dictionaries !!!!!
            (*all_cells)[i]->phenotype.volume.multiply_by_ratio(volume_increase_ratio);
            
			// Note to Furkan : UPTAKE RATES!!!!!!!!
            //(*all_cells)[i]->phenotype.secretion.uptake_rates[gln_index]=0.0;
			
			double glucose_consumption = result[1]/10/60  * intracellular_dt; //DNN multipliers
            double glutamine_consumption = result[2]/10/60 * intracellular_dt;
            double lactate_consumption = result[3]/10/60 * intracellular_dt;
			
			
/* 			std::cout << "Biomass = " << biomass_creation_flux << std::endl;
            std::cout << "Intracellular Glucose Change = " << result[1]/10 << std::endl;
			std::cout << "Intracellular Glutamine Change = " << result[2]/10 << std::endl;
			std::cout << "Intracellular Lactate Change = " << result[3]/10 << std::endl; */
			
			
			//std::cout << "Before : "  << (*all_cells)[i]->custom_data[i_Lac_i] <<std::endl;
            (*all_cells)[i]->custom_data[i_Glc_i] = (*all_cells)[i]->custom_data[i_Glc_i]-glucose_consumption;
            (*all_cells)[i]->custom_data[i_Gln_i] = (*all_cells)[i]->custom_data[i_Gln_i]-glutamine_consumption;
            (*all_cells)[i]->custom_data[i_Lac_i] = (*all_cells)[i]->custom_data[i_Lac_i]-lactate_consumption;
			//std::cout << "After : "  << (*all_cells)[i]->custom_data[i_Lac_i] <<std::endl;

			if ((*all_cells)[i]->custom_data[i_Glc_i] < 0)
			{
				(*all_cells)[i]->custom_data[i_Glc_i] = 0.0;
			}
						
			if ((*all_cells)[i]->custom_data[i_Gln_i] < 0)
			{
				(*all_cells)[i]->custom_data[i_Gln_i] = 0.0;
			}

			if ((*all_cells)[i]->custom_data[i_Lac_i] < 0)
			{
				(*all_cells)[i]->custom_data[i_Lac_i] = 0.0;
			}
/* 
			std::cout << "Biomass Result = " << biomass_creation_flux << std::endl;
            std::cout << "Intracellular Glucose Concentration = " << (*all_cells)[i]->custom_data[i_Glc_i] << std::endl;
			std::cout << "Intracellular Glutamine Concentration = " << (*all_cells)[i]->custom_data[i_Gln_i] << std::endl;
			std::cout << "Intracellular Lactate Concentration = " << (*all_cells)[i]->custom_data[i_Lac_i] << std::endl; */


            
            double cell_pressure = (*all_cells)[i]->state.simple_pressure;
            if ( (*all_cells)[i]->phenotype.volume.total > 2494*2)
                {
                    (*all_cells)[i]->phenotype.cycle.data.transition_rate(0,0) = 9e99;
                }
            else
                {
                    (*all_cells)[i]->phenotype.cycle.data.transition_rate(0,0) = 0.0;
                }
        }
        
        /* // KRAS type simulation
        else if ((*all_cells)[i]->type == 1)
        {  
            keras2cpp::Tensor in{5};
            keras2cpp::Tensor out;
            double glc_val_int = (*all_cells)[i]->nearest_density_vector()[glc_index];
            double gln_val_int = (*all_cells)[i]->nearest_density_vector()[gln_index];
            
            
            double u_glc = (*all_cells)[i]->custom_data[1] * exp_ave_n_cells / exp_vol_well * glc_val_int;
            double u_gln = (*all_cells)[i]->custom_data[2] * exp_ave_n_cells / exp_vol_well * gln_val_int;
            
            float fl_glc = u_glc;
            float fl_gln = u_gln;
            
            //std::cout << "Glucose = " << fl_glc << std::endl;
            //std::cout << "Glutamine = " << fl_gln << std::endl;    
            
            in.data_ = {0.223,0.003,19.2,1.08,1.56};
            out = KRAS_Model(in); // model evaluation
            out.print();
            
            std::vector<double> result;
            result = out.result_vector();
            // std::vector<double> result = out.result_vector();
            
            double biomass_creation_flux = result[0]/parameters.doubles("DNN_biomass_normalizer");
            
            //(*all_cells)[i]->custom_data[biomass_vi]  = biomass_creation_flux;
            
            double volume_increase_ratio = 1 + ( biomass_creation_flux / 60 * intracellular_dt);
            (*all_cells)[i]->custom_data[0]  = biomass_creation_flux; // FURKAN to Fix = Manually written indices for custom data - USE dictionaries !!!!!
            (*all_cells)[i]->custom_data[3]  = fl_glc;
            (*all_cells)[i]->custom_data[4]  = fl_gln;
            (*all_cells)[i]->phenotype.volume.multiply_by_ratio(volume_increase_ratio);
            
            (*all_cells)[i]->phenotype.secretion.uptake_rates[glc_index]=fl_glc;
            (*all_cells)[i]->phenotype.secretion.uptake_rates[gln_index]=fl_gln;
            
            double cell_pressure = (*all_cells)[i]->state.simple_pressure;

            if ( (*all_cells)[i]->phenotype.volume.total > 2494*2)
                {
                    (*all_cells)[i]->phenotype.cycle.data.transition_rate(0,0) = 9e99;
                }
            else
                {
                    (*all_cells)[i]->phenotype.cycle.data.transition_rate(0,0) = 0.0;
                }
        } */
    }
}