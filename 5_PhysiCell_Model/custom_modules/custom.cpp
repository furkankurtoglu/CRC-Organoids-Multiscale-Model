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
auto WT_Model = keras2cpp::Model::load("WT_DNN_including_intracellular.model");
auto KRAS_Model = keras2cpp::Model::load("KRAS_DNN.model");

void create_cell_types( void )
{
    read_DNN("WT_DNN_including_intracellular.model", "KRAS_DNN.model");

	// set the random seed 
	SeedRandom(time(0));  
	
	/* 
	   Put any modifications to default cell definition here if you 
	   want to have "inherited" by other cell types. 
	   
	   This is a good place to set default functions. 
	*/ 
	
	initialize_default_cell_definition(); 
	cell_defaults.phenotype.secretion.sync_to_microenvironment( &microenvironment ); 
	
	cell_defaults.functions.volume_update_function = NULL;
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
	
	//std::cout << "entering setup tissue is done" << std::endl;
    static int glucose_substrate_index = microenvironment.find_density_index( "glucose" );
    static int glutamine_substrate_index = microenvironment.find_density_index( "glutamine" ); 
    static int lactate_substrate_index = microenvironment.find_density_index( "lactate");

	// place a cluster of tumor cells at the center 
 
	double tumor_radius = parameters.doubles( "tumor_radius" ); // 250.0; 
	
	double initial_tumor_radius = parameters.doubles( "tumor_radius" ); 
	
    bool Two_Dim_Seeding =  parameters.bools( "two_dim_seeding" );
    
    std::cout << Two_Dim_Seeding << std::endl;

    double xmin = microenvironment.mesh.bounding_box[0]; 
    double ymin = microenvironment.mesh.bounding_box[1]; 
    double zmin = microenvironment.mesh.bounding_box[2]; 

    double xmax = microenvironment.mesh.bounding_box[3]; 
    double ymax = microenvironment.mesh.bounding_box[4]; 
    double zmax = microenvironment.mesh.bounding_box[5];
    
    Cell* pCell; 
	
	double random_number = uniform_random();

    if ( Two_Dim_Seeding == true )
    {
        // 2D Cell Seeding 
        
        // place CRC_WT 
        //Cell_Definition* pCD = find_cell_definition( "CRC_WT"); 
        srand(time(0));
        double min_number_of_monolayer_cells = parameters.doubles( "minimum_number_of_monolayer_cells" );
        double max_number_of_monolayer_cells = parameters.doubles( "maximum_number_of_monolayer_cells" );
        ////////////////////////////////////////////////////////////////
        double min_g6p_conc = parameters.doubles( "min_g6p_conc" );
        double max_g6p_conc = parameters.doubles( "max_g6p_conc" );
        double min_fbp_conc = parameters.doubles( "min_fbp_conc" );
        double max_fbp_conc = parameters.doubles( "max_fbp_conc" );
        double min_g3p_conc = parameters.doubles( "min_g3p_conc" );
        double max_g3p_conc = parameters.doubles( "max_g3p_conc" );
        double min_pep_conc = parameters.doubles( "min_pep_conc" );
        double max_pep_conc = parameters.doubles( "max_pep_conc" );
        double min_lac_conc = parameters.doubles( "min_lac_conc" );
        double max_lac_conc = parameters.doubles( "max_lac_conc" );
        double min_gln_conc = parameters.doubles( "min_gln_conc" );
        double max_gln_conc = parameters.doubles( "max_gln_conc" );
        double min_glu_conc = parameters.doubles( "min_glu_conc" );
        double max_glu_conc = parameters.doubles( "max_glu_conc" );
        
        
        ////////
        
        double number_of_monolayer_cells = min_number_of_monolayer_cells + (max_number_of_monolayer_cells - min_number_of_monolayer_cells) * UniformRandom();
        for( double n = 0 ; n < number_of_monolayer_cells ; n++ )
        {
            
            //std::cout << "Random Number " << uniform_random() << " ... " << std::endl; 
            
            std::vector<double> positions = {0,0,0}; 
            double xrand = (rand() % int(xmax*2)) - xmax;
            double yrand = (rand() % int(ymax*2)) - ymax;
            double zrand = (rand() % int(zmax*2)) - zmax;
            if (xrand < xmin) xmin = xrand;
            if (xrand > xmax) xmax = xrand;
            if (yrand < ymin) ymin = yrand;
            if (yrand > ymax) ymax = yrand;
            if (zrand < zmin) zmin = zrand;
            if (zrand > zmax) zmax = zrand;
            positions[0] += -500;//(rand() % 5333) - 2666;  ----height
            positions[1] += yrand;//yrand;//(rand() % 961) - 480; --width
            positions[2] += zrand;//(rand() % 5333) - 2666; --widht
            
            pCell = create_cell( get_cell_definition("CRC_WT") ); 
            
            static int i_g6p_i = pCell->custom_data.find_variable_index( "int_g6p" );
            static int i_fbp_i = pCell->custom_data.find_variable_index( "int_fbp" );
            static int i_g3p_i = pCell->custom_data.find_variable_index( "int_g3P" );
            static int i_pep_i = pCell->custom_data.find_variable_index( "int_pep" );
            static int i_lac_i = pCell->custom_data.find_variable_index( "int_lac" );
            static int i_gln_i = pCell->custom_data.find_variable_index( "int_gln" );
            static int i_glu_i = pCell->custom_data.find_variable_index( "int_glu" );
            
            pCell->assign_position( positions );
            pCell->state.orientation = {1,0,0};
            pCell->phenotype.geometry.polarity = 1;
            
            pCell->custom_data[i_g6p_i] = min_g6p_conc + (max_g6p_conc - min_g6p_conc) * UniformRandom();//
            pCell->custom_data[i_fbp_i] = min_fbp_conc + (max_fbp_conc - min_fbp_conc) * UniformRandom();//
            pCell->custom_data[i_g3p_i] = min_g3p_conc + (max_g3p_conc - min_g3p_conc) * UniformRandom();//
            pCell->custom_data[i_pep_i] = min_pep_conc + (max_pep_conc - min_pep_conc) * UniformRandom();//
            pCell->custom_data[i_lac_i] = min_lac_conc + (max_lac_conc-min_lac_conc) * UniformRandom();//
            pCell->custom_data[i_gln_i] = min_gln_conc + (max_gln_conc - min_gln_conc) * UniformRandom();//
            pCell->custom_data[i_glu_i] = min_glu_conc + (max_glu_conc- min_glu_conc) * UniformRandom();// 
            
            
            double a = uniform_random();
                if ( a < 0.2)
                {   
                    pCell->phenotype.volume.multiply_by_ratio(1 + a);
                } 
                if ( a > 0.8)
                {   
                    pCell->phenotype.volume.multiply_by_ratio(1 * a);
                } 
            
        }
    }
    else
    {
		double cell_radius = cell_defaults.phenotype.geometry.radius; 
		// double initial_tumor_radius = 46; // parameters.doubles("initial_tumor_radius");
		double initial_tumor_radius = parameters.doubles("initial_tumor_radius");
		std::cout << "------ initial_tumor_radius = " << initial_tumor_radius << std::endl;

		int number_of_organoids = parameters.ints("number_of_organoids");
		std::cout << "------ number_of_organoids = " << number_of_organoids << std::endl;		
		double tumor_radius = initial_tumor_radius;
        for (int j = 0; j < number_of_organoids; j++) // seeding number of organoid cells specified in PhysiCell_settings.xml
            {
                double spheroid_multiplier = UniformRandom();
                if ( spheroid_multiplier < 0.4 )
                {
                    tumor_radius = initial_tumor_radius * ( 1 + spheroid_multiplier);
                }                    
                if ( spheroid_multiplier > 0.6 )
                {
                    tumor_radius = initial_tumor_radius * spheroid_multiplier;
                }    
                //std::cout << tumor_radius << std::endl;
                std::vector<std::vector<double>> positions = create_cell_sphere_positions(cell_radius,tumor_radius); 
                //std::cout << "creating " << positions.size() << " closely-packed organoid cells ... " << std::endl;
                // create organoid
                    double xrand = (rand() % (int(xmax-100)*2)) - (xmax-100);
                    double yrand = (rand() % (int(ymax-100)*2)) - (ymax-100);
                    double zrand = (rand() % (int(zmax-100)*2)) - (zmax-100);

                    for( int i=0; i < positions.size(); i++ )
                    {
                        positions[i][0] += xrand;
                        positions[i][1] += yrand;
                        positions[i][2] += zrand;
                        pCell = create_cell( get_cell_definition("CRC_WT") );
                        static int i_g6p_i = pCell->custom_data.find_variable_index( "int_g6p" );
                        static int i_fbp_i = pCell->custom_data.find_variable_index( "int_fbp" );
                        static int i_g3p_i = pCell->custom_data.find_variable_index( "int_g3P" );
                        static int i_pep_i = pCell->custom_data.find_variable_index( "int_pep" );
                        static int i_lac_i = pCell->custom_data.find_variable_index( "int_lac" );
                        static int i_gln_i = pCell->custom_data.find_variable_index( "int_gln" );
                        static int i_glu_i = pCell->custom_data.find_variable_index( "int_glu" );
                        //int organ_id = pCell->custom_data.find_variable_index( "int_lac" );
                        //std::cout << pCell->custom_data.find_variable_index( "lactate_result" ) << std::endl;
                        pCell->assign_position( positions[i] );
            
                        pCell->custom_data[i_g6p_i] = 1.68 + (6 - 1.68) * UniformRandom();//
                        pCell->custom_data[i_fbp_i] = 0.16 + (1.08 - 0.16) * UniformRandom();//
                        pCell->custom_data[i_g3p_i] = 0.24 + (0.84 - 0.24) * UniformRandom();//
                        pCell->custom_data[i_pep_i] = 0.24 + (0.84 - 0.24) * UniformRandom();//
                        pCell->custom_data[i_lac_i] = 6.4 + (19.2-6.4) * UniformRandom();//
                        pCell->custom_data[i_gln_i] = 0.56 + (1.08 - 0.56) * UniformRandom();//
                        pCell->custom_data[i_glu_i] = 0.88 + (1.56 - 0.88) * UniformRandom();// 
                        double a = UniformRandom();
                        if ( a < 0.2)
                        {   
                            pCell->phenotype.volume.multiply_by_ratio(1 + a);
                        } 
                        if ( a > 0.8)
                        {   
                            pCell->phenotype.volume.multiply_by_ratio(1 * a);
                        } 
                        double t = j; 
                        int id_number = pCell->custom_data.find_variable_index( "lactate_result" );
                        pCell->custom_data[id_number] = t;
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

    
    #pragma omp parallel for 
    for( int i=0; i < (*all_cells).size(); i++ )
    {
        //std::cout << "Entering the DNN simulations" << std::endl;
        // Wild type simulation
        if ( (*all_cells)[i]->is_out_of_domain == false )
            {  
            if ((*all_cells)[i]->type == 0)
            {
                //housekeeping
                static int i_g6p_i = (*all_cells)[i]->custom_data.find_variable_index( "int_g6p" );
                static int i_fbp_i = (*all_cells)[i]->custom_data.find_variable_index( "int_fbp" );
                static int i_g3p_i = (*all_cells)[i]->custom_data.find_variable_index( "int_g3P" );
                static int i_pep_i = (*all_cells)[i]->custom_data.find_variable_index( "int_pep" );
                static int i_lac_i = (*all_cells)[i]->custom_data.find_variable_index( "int_lac" );
                static int i_gln_i = (*all_cells)[i]->custom_data.find_variable_index( "int_gln" );
                static int i_glu_i = (*all_cells)[i]->custom_data.find_variable_index( "int_glu" );
                static int biomass_result = (*all_cells)[i]->custom_data.find_variable_index( "biomass_flux" );
                
                double k1_glucose = get_single_signal((*all_cells)[i],"custom:k1_glucose");
                double k1_glutamine = get_single_signal((*all_cells)[i],"custom:k1_glutamine" );
                
                
                keras2cpp::Tensor in{9};
                keras2cpp::Tensor out;

                // get pressure
                double cell_pressure = (*all_cells)[i]->state.simple_pressure;
                double cell_pressure_threshold = 0.55;
                    
                if (cell_pressure > cell_pressure_threshold)
                {
                    cell_pressure = cell_pressure_threshold;
                }

                // update exchange rates
                double glc_local_concentration = (*all_cells)[i]->nearest_density_vector()[glc_index];
                double gln_local_concentration = (*all_cells)[i]->nearest_density_vector()[gln_index];
                double glucose_boundary_for_fba = glc_local_concentration * k1_glucose;
                double glutamine_boundary_for_fba = gln_local_concentration * k1_glutamine;

                //****************
                //update exchange boundaries for FBA
                double u_glc = glucose_boundary_for_fba;
                double u_gln = glutamine_boundary_for_fba;
                
                float fl_glc = u_glc;
                float fl_gln = u_gln;
 /*                std::cout << "double : " << u_glc <<std::endl;
                std::cout << "float : " << fl_glc <<std::endl;
                 */
                
                //update intracellular boundaries for FBA

                float fl_int_conc_g6p = (*all_cells)[i]->custom_data[i_g6p_i];
                float fl_int_conc_fbp = (*all_cells)[i]->custom_data[i_fbp_i];
                float fl_int_conc_g3p = (*all_cells)[i]->custom_data[i_g3p_i];
                float fl_int_conc_pep = (*all_cells)[i]->custom_data[i_pep_i];
                float fl_int_conc_lac = (*all_cells)[i]->custom_data[i_lac_i];
                float fl_int_conc_gln = (*all_cells)[i]->custom_data[i_gln_i];
                float fl_int_conc_glu = (*all_cells)[i]->custom_data[i_glu_i];

/*                 std::cout << "FBA MODEL UPDATE" << std::endl; 
                std::cout << "glc exchange boundary : " << u_glc << std::endl;
                std::cout << "gln exchange boundary: " << u_gln << std::endl;           
                std::cout << "Input Intracellular Glucose = " << fl_int_conc_glu << std::endl;
                std::cout << "Input Intracellular Glutamine = " << fl_int_conc_gln << std::endl;    
                std::cout << "Input Intracellular Lactate = " << fl_int_conc_lac << std::endl;    
                std::cout << "-----------------------------" << std::endl;
                std::cout << std::endl; */
                
                in.data_ = {fl_glc,fl_gln,fl_int_conc_g6p,fl_int_conc_fbp,fl_int_conc_g3p,fl_int_conc_pep,fl_int_conc_lac,fl_int_conc_gln,fl_int_conc_glu};
                //in.data_ = {0.223,0.003,6.4,0.56,1.08}; //test values
                
                // model evaluation
                //out = WT_Model(in); 
                //std::cout << "ARRIVED" << std::endl;
                //out.print();
                
                std::vector<double> result;
                result = out.result_vector();
                
                double biomass_creation_flux = result[0]/parameters.doubles("DNN_biomass_normalizer");
                
                (*all_cells)[i]->phenotype.secretion.net_export_rates[glc_index] = glucose_boundary_for_fba * 1000 /60; // 1000 cubic-micron = V_cell, E in mmol/min
                (*all_cells)[i]->phenotype.secretion.net_export_rates[gln_index] = glucose_boundary_for_fba * 1000 /60; // 1000 cubic-micron = V_cell 
                
                
                
                (*all_cells)[i]->custom_data[biomass_result]  = biomass_creation_flux;
                //std::cout << "Biomass Before = " << (*all_cells)[i]->phenotype.volume.total << std::endl;
                
                // update cellular volume
                //double volume_increase_ratio = 1 + ( biomass_creation_flux / 60 * intracellular_dt);
                double volume_increase_ratio = 1 + ( 0.044 / 60 * intracellular_dt);
                (*all_cells)[i]->custom_data[biomass_result]  = biomass_creation_flux;
                //std::cout << "Biomass = " << biomass_creation_flux << std::endl;
                (*all_cells)[i]->phenotype.volume.multiply_by_ratio(volume_increase_ratio);
                //std::cout << "Cell ID : " << (*all_cells)[i]->ID <<  "    Biomass After = " << (*all_cells)[i]->phenotype.volume.total << std::endl;


                
                if ( (*all_cells)[i]->phenotype.volume.total > 1000*2) // 
                    {        
                       //if (cell_pressure < 0.55)
                       //{
                       (*all_cells)[i]->phenotype.cycle.data.transition_rate(0,0) = 9e99;
                       //}
                    }
                else
                    {
                        (*all_cells)[i]->phenotype.cycle.data.transition_rate(0,0) = 0.0;
                    }
            }
        }
    }
}


void create_coarse_microenvironment(Microenvironment* Pcoarse_well)
{
    // Annotate
    Pcoarse_well->name = "coarse_well";
    Pcoarse_well->spatial_units = "micron";
    Pcoarse_well->mesh.units = "micron";
    Pcoarse_well->time_units = "min";
    
    // Add metabolites to coarse well
    Pcoarse_well->set_density( 0 , "glucose", "mM", 30000 , 0.00 ); // Ref: https://www.ncbi.nlm.nih.gov/pubmed/7597991 
    Pcoarse_well->add_density( "glutamine", "mM", 30000 , 0.0 );
    Pcoarse_well->add_density( "lactate", "mM", 30000 , 0.0);
    
    
    // 1D Voxels
    double dx = 32;
    double dy = 5760;
    double dz = 5760;
    
    // 1D Microenvironment creation
    Pcoarse_well->resize_space( 512.0, 5120.0, -2880.0, 2880.0, -2880.0, 2880.0, dx, dy, dz );


    // No Dirichlet Conditions
    Pcoarse_well->set_substrate_dirichlet_activation(0,false);
    Pcoarse_well->set_substrate_dirichlet_activation(1,false);
    Pcoarse_well->set_substrate_dirichlet_activation(2,false);
    
    
    // Set 1D solver
    Pcoarse_well->diffusion_decay_solver = diffusion_decay_solver__constant_coefficients_LOD_1D;   
    
    // Set Initial Conditions
    int coarse_well_voxel_number = Pcoarse_well->mesh.voxels.size();    
    for ( int m = 0; m < coarse_well_voxel_number ; m++)
    {
        (*Pcoarse_well)(m)[0]=17.5; // glucose
        (*Pcoarse_well)(m)[1]=5.5; // glutamine
        (*Pcoarse_well)(m)[2]=0; //lactate
    }

    Pcoarse_well->display_information( std::cout );
    
    // Save the initial data
    Pcoarse_well->write_to_matlab("output/output00000000_microenvironment1.mat"); 

    
    return;
}


void create_transfer_region(Microenvironment* PTransfer)
{
    // Annotate
    PTransfer->name = "transfer_region";
    PTransfer->spatial_units = "micron";
    PTransfer->mesh.units = "micron";
    PTransfer->time_units = "min";

    // Add metabolites to transfer region
    PTransfer->set_density( 0 , "glucose", "mmHg", 30000 , 0.00 );
    PTransfer->add_density( "glutamine", "mM", 30000 , 0.0 );
    PTransfer->add_density( "lactate", "mM", 30000 , 0.0);

      
    // 3D Voxel
    double tr_dx = 32;
    double tr_dy = 32;
    double tr_dz = 32;
    
    // Set 3D meshes
    PTransfer->resize_space( 480.0, 544.0, -2880, 2880, -2880, 2880, tr_dx, tr_dy, tr_dz );
    
    // No Dirichlet Conditions
    PTransfer->set_substrate_dirichlet_activation(0,false);
    PTransfer->set_substrate_dirichlet_activation(1,false);
    PTransfer->set_substrate_dirichlet_activation(2,false);
    
    // Set 3D solver
    PTransfer->diffusion_decay_solver = diffusion_decay_solver__constant_coefficients_LOD_3D;   
    
    // Set Initial Conditions
    int transfer_region_voxel_number = PTransfer->mesh.voxels.size(); 
    for ( int m = 0; m < transfer_region_voxel_number ; m++)
    {
        (*PTransfer)(m)[0]=17.5; // glucose
        (*PTransfer)(m)[1]=5.5; // glutamine
        (*PTransfer)(m)[2]=0; //lactate
    }
    
    PTransfer->display_information( std::cout );
    // Save the initial data
    PTransfer->write_to_matlab("output/output00000000_microenvironment2.mat");    
    
    return;
}