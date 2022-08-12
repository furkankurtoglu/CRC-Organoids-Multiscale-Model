# first neural network with keras tutorial
from numpy import loadtxt
from keras.models import Sequential
from keras.layers import Dense
import numpy as np
import matplotlib.pyplot as plt

# load the dataset
dataset = loadtxt('WT_in_silico_data.csv', delimiter=',')
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



# split into input (X) and output (y) variables
X = training_data[:,0:2]
y = training_data[:,3]


# define the keras model
model = Sequential()
model.add(Dense(4, input_shape=(2,), activation='relu'))
model.add(Dense(12, activation='relu'))
model.add(Dense(1, activation='relu'))
# compile the keras model
model.compile(loss='mse', optimizer='adam', metrics=['mae'])
# fit the keras model on the dataset
history = model.fit(X, y, epochs=150, batch_size=100)

plt.plot(history.history['mae'],'o',color='black')
plt.title('mean absolutue error')
plt.ylabel('mea')
plt.xlabel('number of epochs')
plt.show()



#%%
x_test_data = test_data[:,0:2]
y_test_data = test_data[:,3]


# evaluate the keras model
print('EVALUATION')
_, accuracy = model.evaluate(x_test_data, y_test_data)


#%% glucose = 0.223 mM/hr
glutamine_values_0_223 = []
biomass_reals_0_223 = []
biomass_predictions_0_223 = []
for i in range(0,100):
    glutamine_value = sorted_data[i,1]
    biomass_real = sorted_data[i,3]
    biomass_predicted =  model.predict([[0.223, glutamine_value]])
  #  print(biomass_predicted)
    glutamine_values_0_223.append(glutamine_value)
    biomass_reals_0_223.append(biomass_real)
    biomass_predictions_0_223.append(biomass_predicted[0])

fig = plt.figure()
plt.plot(glutamine_values_0_223,biomass_reals_0_223,'ko')
plt.plot(glutamine_values_0_223,biomass_reals_0_223,'b')
plt.xlabel('glutamine lower boundaries')
plt.ylabel('biomass growth')
plt.title('Glucose lb = 0.223 mM')


# glucose = 0.162181818181818 mM/hr
glutamine_values_0_162 = []
biomass_reals_0_162 = []
biomass_predictions_0_162 = []
for i in range(2700,2800):
    glutamine_value = sorted_data[i,1]
    biomass_real = sorted_data[i,3]
    biomass_predicted =  model.predict([[0.162181818181818, glutamine_value]])
  #  print(biomass_predicted)
    glutamine_values_0_162.append(glutamine_value)
    biomass_reals_0_162.append(biomass_real)
    biomass_predictions_0_162.append(biomass_predicted[0])

fig = plt.figure()
plt.plot(glutamine_values_0_162,biomass_reals_0_162,'ko')
plt.plot(glutamine_values_0_162,biomass_reals_0_162,'b')
plt.xlabel('glutamine lower boundaries')
plt.ylabel('biomass growth')
plt.title('Glucose lb = 0.162 mM')


# glucose = 0.0743333333333333 mM/hr
glutamine_values_0_074 = []
biomass_reals_0_074 = []
biomass_predictions_0_074 = []
for i in range(6600,6700):
    glutamine_value = sorted_data[i,1]
    biomass_real = sorted_data[i,3]
    biomass_predicted =  model.predict([[0.0743333333333333, glutamine_value]])
  #  print(biomass_predicted)
    glutamine_values_0_074.append(glutamine_value)
    biomass_reals_0_074.append(biomass_real)
    biomass_predictions_0_074.append(biomass_predicted[0])

fig = plt.figure()
plt.plot(glutamine_values_0_074,biomass_reals_0_074,'ko')
plt.plot(glutamine_values_0_074,biomass_reals_0_074,'b')
plt.xlabel('glutamine lower boundaries')
plt.ylabel('biomass growth')
plt.title('Glucose lb = 0.074 mM')


# glucose = 0.0 mM/hr
glutamine_values_0_162 = []
biomass_reals_0_162 = []
biomass_predictions_0_162 = []
for i in range(2700,2800):
    glutamine_value = sorted_data[i,1]
    biomass_real = sorted_data[i,3]
    biomass_predicted =  model.predict([[0.162181818181818, glutamine_value]])
  #  print(biomass_predicted)
    glutamine_values_0_162.append(glutamine_value)
    biomass_reals_0_162.append(biomass_real)
    biomass_predictions_0_162.append(biomass_predicted[0])

fig = plt.figure()
plt.plot(glutamine_values_0_162,biomass_reals_0_162,'ko')
plt.plot(glutamine_values_0_162,biomass_reals_0_162,'b')
plt.xlabel('glutamine lower boundaries')
plt.ylabel('biomass growth')
plt.title('Glucose lb = 0.162 mM')





#%%

# oxygen_value_100 = 100;
# matching_indices = np.where(sorted_data[:,2] == oxygen_value_100)

# only_oxygen_100 = sorted_data[matching_indices,:]
# unique_rows_100 = np.unique(only_oxygen_100, axis=1)

# glucose_concentrations = unique_rows_100[0,:,0]

# glucose_values_100 = []
# biomass_reals_100 = []
# biomass_predictions_ox_100 = []

# for g in range(0,np.size(unique_rows_100[0,:,0])):
#     glucose_value = unique_rows_100[0,g,0]
#     biomass_value = unique_rows_100[0,g,3]
#     biomass_predicted = model.predict([[glucose_value, 0 ,oxygen_value_100]])
#     glucose_values_100.append(glucose_value)
#     biomass_reals_100.append(biomass_value)
#     biomass_predictions_ox_100.append(biomass_predicted[0])

# fig = plt.figure()
# plt.plot(glucose_concentrations,biomass_reals_100,'ko')
# plt.plot(glucose_concentrations,biomass_predictions_ox_100,'b')
# plt.xlabel('glucose lower boundaries')
# plt.ylabel('biomass growth')
# plt.title('oxygen lb = 100')




# oxygen_value_63 = 63.1578947368421;
# matching_indices = np.where(sorted_data[:,2] == oxygen_value_63)

# only_oxygen_63 = sorted_data[matching_indices,:]
# unique_rows_63 = np.unique(only_oxygen_63, axis=1)

# glucose_concentrations = unique_rows_63[0,:,0]

# glucose_values_63 = []
# biomass_reals_63 = []
# biomass_predictions_ox_63 = []

# for g in range(0,np.size(unique_rows_63[0,:,0])):
#     glucose_value = unique_rows_63[0,g,0]
#     biomass_value = unique_rows_63[0,g,3]
#     biomass_predicted = model.predict([[glucose_value, 0 ,oxygen_value_63]])
#     glucose_values_63.append(glucose_value)
#     biomass_reals_63.append(biomass_value)
#     biomass_predictions_ox_63.append(biomass_predicted[0])

# fig = plt.figure()
# plt.plot(glucose_concentrations,biomass_reals_63,'ko')
# plt.plot(glucose_concentrations,biomass_predictions_ox_63,'b')
# plt.xlabel('glucose lower boundaries')
# plt.ylabel('biomass growth')
# plt.title('oxygen lb = 63')




# oxygen_value_26 = 26.3157894736842;
# matching_indices = np.where(sorted_data[:,2] == oxygen_value_26)

# only_oxygen_26 = sorted_data[matching_indices,:]
# unique_rows_26 = np.unique(only_oxygen_26, axis=1)

# glucose_concentrations = unique_rows_26[0,:,0]

# glucose_values_26 = []
# biomass_reals_26 = []
# biomass_predictions_ox_26 = []

# for g in range(0,np.size(unique_rows_63[0,:,0])):
#     glucose_value = unique_rows_26[0,g,0]
#     biomass_value = unique_rows_26[0,g,3]
#     biomass_predicted = model.predict([[glucose_value, 0 ,oxygen_value_26]])
#     glucose_values_26.append(glucose_value)
#     biomass_reals_26.append(biomass_value)
#     biomass_predictions_ox_26.append(biomass_predicted[0])

# fig = plt.figure()
# plt.plot(glucose_concentrations,biomass_reals_26,'ko')
# plt.plot(glucose_concentrations,biomass_predictions_ox_26,'b')
# plt.xlabel('glucose lower boundaries')
# plt.ylabel('biomass growth')
# plt.title('oxygen lb = 26')



# oxygen_value_0 = -0.;
# matching_indices = np.where(sorted_data[:,2] == oxygen_value_0)

# only_oxygen_0 = sorted_data[matching_indices,:]
# unique_rows_0 = np.unique(only_oxygen_0, axis=1)

# glucose_concentrations = unique_rows_0[0,:,0]

# glucose_values_0 = []
# biomass_reals_0 = []
# biomass_predictions_ox_0 = []

# for g in range(0,np.size(unique_rows_0[0,:,0])):
#     glucose_value = unique_rows_0[0,g,0]
#     biomass_value = unique_rows_0[0,g,3]
#     biomass_predicted = model.predict([[glucose_value, 0 ,oxygen_value_0]])
#     glucose_values_0.append(glucose_value)
#     biomass_reals_0.append(biomass_value)
#     biomass_predictions_ox_0.append(biomass_predicted[0])

# fig = plt.figure()
# plt.plot(glucose_concentrations,biomass_reals_0,'ko')
# plt.plot(glucose_concentrations,biomass_predictions_ox_0,'b')
# plt.xlabel('glucose lower boundaries')
# plt.ylabel('biomass growth')
# plt.title('oxygen lb = 0')





#%%
#save model
from keras2cpp import export_model
export_model(model, 'Wild_Type.model')