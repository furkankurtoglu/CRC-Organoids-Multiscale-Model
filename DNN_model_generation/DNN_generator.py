# first neural network with keras tutorial
from numpy import loadtxt
from keras.models import Sequential
from keras.layers import Dense
import numpy as np
import matplotlib.pyplot as plt
from keras2cpp import export_model


cell_type = 'KRAS'
train_model = 'Y'
save_model = 'Y'
draw_convergence = 'Y'
plot_verification_results = 'Y'






# load the dataset
dataname = cell_type + '_in_silico_data.csv'
dataset = loadtxt(dataname, delimiter=',')
sorted_data = np.array(dataset)
data = np.array(dataset)
np.random.shuffle(data)

percent_training_data = 0.8
splitter = round(data.shape[0] * percent_training_data)

training_data = data[:splitter]
test_data = data[splitter:]


# glu_vals = dataset[:,0]
# glt_vals = dataset[:,1]
# biomass_vals = dataset[:,3]


# X, Y = np.meshgrid(glu_vals, glt_vals)
# biomass_vals3d = np.array([biomass_vals,biomass_vals])

# fig = plt.figure(figsize=(12, 12))
# ax = fig.add_subplot(projection='3d')
# ax.scatter(glu_vals, glt_vals, biomass_vals)
# ax.set_xlabel('glucose')
# ax.set_ylabel('glutamine')
# ax.set_zlabel('biomass')
# ax.set_zlim3d(0,0.060)
# plt.show()



biomass_multiplier = 100


# split into input (X) and output (y) variables
X = training_data[:,0:5]
y = training_data[:,5:9]
y[:,0] *= biomass_multiplier
y[:,2] *= -1



if (train_model == 'Y'):
    # define the keras model
    model = Sequential()
    model.add(Dense(10, input_shape=(5,), activation='relu'))
    model.add(Dense(20, activation='relu'))
    model.add(Dense(4, activation='relu'))
    # compile the keras model
    model.compile(loss='mse', optimizer='adam', metrics=['mae'])
    # fit the keras model on the dataset
    history = model.fit(X, y, epochs=20, batch_size=5)

if (draw_convergence == 'Y'):
    plt.plot(history.history['mae'],'o',color='black')
    plt.title('mean absolutue error')
    plt.ylabel('mea')
    plt.xlabel('number of epochs')
    plt.show()



#%%
x_test_data = test_data[:,0:5]
y_test_data = test_data[:,5:9]
y_test_data[:,0] *= biomass_multiplier
y_test_data[:,2] *= -1
# evaluate the keras model
print('EVALUATION')



glucose_value = 0.1115
glutamine_value_001 = 0.00225
lactate_conc = 9.6
glucose_conc = 0.54
glutamine_conc = 0.39


testing = model.predict([[glucose_value,glutamine_value_001,lactate_conc,glucose_conc,glutamine_conc]])
print(testing[0]/biomass_multiplier)



#%%

if (plot_verification_results == 'Y'):
    if (cell_type == 'WT'):
        print('cell type is ' + cell_type)
        glucose_value = 0.16725
        glutamine_value_001 = 0.003
        lactate_conc = 9.6
        glucose_conc = 0.54
        glutamine_conc = 0.78


    if (cell_type == 'KRAS'):
        print('cell type is ' + cell_type)
        glucose_value = 0.105
        glutamine_value_001 = 0.0015
        lactate_conc = 9.6
        glucose_conc = 0.54
        glutamine_conc = 0.78

    
    matching_indices = np.where((sorted_data[:,1] == glutamine_value_001) & (sorted_data[:,2] == lactate_conc) & (sorted_data[:,3] == glucose_conc) & (sorted_data[:,4] == glutamine_conc))
    only_glutamine_001 = sorted_data[matching_indices,:]
    unique_rows_001 = np.unique(only_glutamine_001, axis=1)
    
    glucose_concentrations = unique_rows_001[:,0]
    
    glutamine_values_001 = []
    biomass_reals_001 = []
    biomass_predictions_gln_001 = []
    
    for g in range(0,np.size(unique_rows_001[0,:,0])):
        glucose_value = unique_rows_001[0,g,0]
        biomass_value = unique_rows_001[0,g,5]
        biomass_predicted = model.predict([[glucose_value,glutamine_value_001,lactate_conc,glucose_conc,glutamine_conc]])
        glutamine_values_001.append(glucose_value)
        biomass_reals_001.append(biomass_value)
        biomass_predictions_gln_001.append(biomass_predicted[0][0]/biomass_multiplier)
    
    fig = plt.figure()
    plt.plot(glutamine_values_001,biomass_reals_001,'ko')
    plt.plot(glutamine_values_001,biomass_predictions_gln_001,'r')
    plt.xlabel('glucose uptake rate (mM/hr)')
    plt.ylabel('biomass growth (1/hr)')
    
    #%%
    
    
    
    matching_indices = np.where((sorted_data[:,0] == glucose_value) & (sorted_data[:,2] == lactate_conc) & (sorted_data[:,3] == glucose_conc) & (sorted_data[:,4] == glutamine_conc))
    
    
    only_glutamine_001 = sorted_data[matching_indices,:]
    unique_rows_001 = np.unique(only_glutamine_001, axis=1)
    
    glucose_concentrations = unique_rows_001[:,1]
    
    glutamine_values_001 = []
    biomass_reals_001 = []
    biomass_predictions_gln_001 = []
    
    for g in range(0,np.size(unique_rows_001[0,:,0])):
        glutamine_value = unique_rows_001[0,g,1]
        biomass_value = unique_rows_001[0,g,5]
        biomass_predicted = model.predict([[glucose_value,glutamine_value,lactate_conc,glucose_conc,glutamine_conc]])
        glutamine_values_001.append(glutamine_value)
        biomass_reals_001.append(biomass_value)
        biomass_predictions_gln_001.append(biomass_predicted[0][0]/biomass_multiplier)
    
    fig = plt.figure()
    plt.plot(glutamine_values_001,biomass_reals_001,'ko')
    plt.plot(glutamine_values_001,biomass_predictions_gln_001,'r')
    plt.xlabel('glutamine uptake rate (mM/hr)')
    plt.ylabel('biomass growth (1/hr)')
    
    
    
    
    
    #%%
    
    
    matching_indices = np.where((sorted_data[:,1] == glutamine_value_001) & (sorted_data[:,0] == glucose_value) & (sorted_data[:,3] == glucose_conc) & (sorted_data[:,4] == glutamine_conc))
    
    
    only_glutamine_001 = sorted_data[matching_indices,:]
    unique_rows_001 = np.unique(only_glutamine_001, axis=1)
    
    lac_concentrations = unique_rows_001[:,2]
    
    lac_values_001 = []
    biomass_reals_001 = []
    biomass_predictions_lac_001 = []
    
    for g in range(0,np.size(unique_rows_001[0,:,0])):
        lac_value = unique_rows_001[0,g,2]
        biomass_value = unique_rows_001[0,g,5]
        biomass_predicted = model.predict([[glucose_value,glutamine_value_001,lac_value,glucose_conc,glutamine_conc]])
        lac_values_001.append(lac_value)
        biomass_reals_001.append(biomass_value)
        biomass_predictions_lac_001.append(biomass_predicted[0][0]/biomass_multiplier)
    
    fig = plt.figure()
    plt.plot(lac_values_001,biomass_reals_001,'ko')
    plt.plot(lac_values_001,biomass_predictions_lac_001,'r')
    plt.xlabel('initial intracellular lactate concentration (mM)')
    plt.ylabel('biomass growth (1/hr)')
    
    
    
    #%%
    
    
    
    matching_indices = np.where((sorted_data[:,1] == glutamine_value_001) & (sorted_data[:,0] == glucose_value) & (sorted_data[:,2] == lactate_conc) & (sorted_data[:,4] == glutamine_conc))
    
    
    only_glutamine_001 = sorted_data[matching_indices,:]
    unique_rows_001 = np.unique(only_glutamine_001, axis=1)
    
    glu_concentrations = unique_rows_001[:,3]
    
    glu_values_001 = []
    biomass_reals_001 = []
    biomass_predictions_glc_c_001 = []
    
    for g in range(0,np.size(unique_rows_001[0,:,0])):
        glu_value = unique_rows_001[0,g,3]
        biomass_value = unique_rows_001[0,g,5]
        biomass_predicted = model.predict([[glucose_value,glutamine_value_001,lactate_conc,glu_value,glutamine_conc]])
        glu_values_001.append(glu_value)
        biomass_reals_001.append(biomass_value)
        biomass_predictions_glc_c_001.append(biomass_predicted[0][0]/biomass_multiplier)
    
    fig = plt.figure()
    plt.plot(glu_values_001,biomass_reals_001,'ko')
    plt.plot(glu_values_001,biomass_predictions_glc_c_001,'r')
    plt.xlabel('initial intracellular glucose concentration (mM)')
    plt.ylabel('biomass growth (1/hr)')
    
    
    
    
    #%%
    
    
    
    
    matching_indices = np.where((sorted_data[:,1] == glutamine_value_001) & (sorted_data[:,0] == glucose_value) & (sorted_data[:,2] == lactate_conc) & (sorted_data[:,3] == glucose_conc))
    
    
    only_glutamine_001 = sorted_data[matching_indices,:]
    unique_rows_001 = np.unique(only_glutamine_001, axis=1)
    
    gln_concentrations = unique_rows_001[:,3]
    
    gln_values_001 = []
    lac_ex_reals_001 = []
    lac_ex_gln_c_001 = []
    
    for g in range(0,np.size(unique_rows_001[0,:,0])):
        gln_value = unique_rows_001[0,g,4]
        lac_value = unique_rows_001[0,g,8]
        lac_ex_predicted = model.predict([[glucose_value,glutamine_value_001,lactate_conc,glucose_conc,gln_value]])
        gln_values_001.append(gln_value)
        lac_ex_reals_001.append(lac_value)
        lac_ex_gln_c_001.append(lac_ex_predicted[0][3])
    
    fig = plt.figure()
    plt.plot(gln_values_001,lac_ex_reals_001,'ko')
    plt.plot(gln_values_001,lac_ex_gln_c_001,'r')
    plt.xlabel('initial intracellular glutamine concentration (mM)')
    plt.ylabel(' lac creation rate (1/hr)')
    
    #%%
    
    
    
    matching_indices = np.where((sorted_data[:,1] == glutamine_value_001) & (sorted_data[:,2] == lactate_conc) & (sorted_data[:,3] == glucose_conc) & (sorted_data[:,4] == glutamine_conc))
    
    
    only_glutamine_001 = sorted_data[matching_indices,:]
    unique_rows_001 = np.unique(only_glutamine_001, axis=1)
    
    glucose_concentrations = unique_rows_001[:,0]
    
    glutamine_values_001 = []
    biomass_reals_001 = []
    biomass_predictions_gln_001 = []
    
    for g in range(0,np.size(unique_rows_001[0,:,0])):
        glucose_value = unique_rows_001[0,g,0]
        biomass_value = unique_rows_001[0,g,6]
        biomass_predicted = model.predict([[glucose_value,glutamine_value_001,lactate_conc,glucose_conc,glutamine_conc]])
        glutamine_values_001.append(glucose_value)
        biomass_reals_001.append(biomass_value)
        biomass_predictions_gln_001.append(biomass_predicted[0][1]/biomass_multiplier)
    
    fig = plt.figure()
    plt.plot(glutamine_values_001,biomass_reals_001,'ko')
    plt.plot(glutamine_values_001,biomass_predictions_gln_001,'r')
    plt.xlabel('glucose uptake rate (mM/hr)')
    plt.ylabel('glucose consumption rate (1/hr)')
    
    
    
    
    
    #%%
    
    
    
    
    matching_indices = np.where((sorted_data[:,1] == glutamine_value_001) & (sorted_data[:,0] == glucose_value) & (sorted_data[:,2] == lactate_conc) & (sorted_data[:,4] == glutamine_conc))
    
    
    only_glutamine_001 = sorted_data[matching_indices,:]
    unique_rows_001 = np.unique(only_glutamine_001, axis=1)
    
    glu_concentrations = unique_rows_001[:,3]
    
    glu_values_001 = []
    biomass_reals_001 = []
    biomass_predictions_glc_c_001 = []
    
    for g in range(0,np.size(unique_rows_001[0,:,0])):
        glu_value = unique_rows_001[0,g,3]
        biomass_value = unique_rows_001[0,g,7]*-1
        biomass_predicted = model.predict([[glucose_value,glutamine_value_001,lactate_conc,glu_value,glutamine_conc]])
        glu_values_001.append(glu_value)
        biomass_reals_001.append(biomass_value)
        biomass_predictions_glc_c_001.append(biomass_predicted[0][2])
    
    fig = plt.figure()
    plt.plot(glu_values_001,biomass_reals_001,'ko')
    plt.plot(glu_values_001,biomass_predictions_glc_c_001,'r')
    plt.xlabel('initial intracellular glucose concentration (mM)')
    plt.ylabel('intracellular glutamine consumption (1/hr)')
    

    #%%
    
    
    
    matching_indices = np.where((sorted_data[:,1] == glutamine_value_001) & (sorted_data[:,0] == glucose_value) & (sorted_data[:,3] == glucose_conc) & (sorted_data[:,4] == glutamine_conc))
    
    
    only_glutamine_001 = sorted_data[matching_indices,:]
    unique_rows_001 = np.unique(only_glutamine_001, axis=1)
    
    lac_concentrations = unique_rows_001[:,2]
    
    lac_values_001 = []
    biomass_reals_001 = []
    biomass_predictions_lac_001 = []
    
    for g in range(0,np.size(unique_rows_001[0,:,0])):
        lac_value = unique_rows_001[0,g,2]
        biomass_value = unique_rows_001[0,g,8]
        biomass_predicted = model.predict([[glucose_value,glutamine_value_001,lac_value,glucose_conc,glutamine_conc]])
        lac_values_001.append(lac_value)
        biomass_reals_001.append(biomass_value)
        biomass_predictions_lac_001.append(biomass_predicted[0][3])
    
    fig = plt.figure()
    plt.plot(lac_values_001,biomass_reals_001,'ko')
    plt.plot(lac_values_001,biomass_predictions_lac_001,'r')
    plt.xlabel('initial intracellular lactate concentration (mM)')
    plt.ylabel('intracellular lactate creation (1/hr)')


#%%




#%%
#save model
if (save_model == 'Y'):
    model_name = cell_type + '_DNN' + '.model'
    export_model(model, model_name)