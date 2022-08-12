# first neural network with keras tutorial
from numpy import loadtxt
from keras.models import Sequential
from keras.layers import Dense
import numpy as np
import matplotlib.pyplot as plt

# load the dataset
dataset = loadtxt('KRAS_in_silico_data.csv', delimiter=',')
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
    biomass_predicted =  model.predict([[0.21, glutamine_value]])
  #  print(biomass_predicted)
    glutamine_values_0_223.append(glutamine_value)
    biomass_reals_0_223.append(biomass_real)
    biomass_predictions_0_223.append(biomass_predicted[0])

fig = plt.figure()
plt.plot(glutamine_values_0_223,biomass_reals_0_223,'ko')
plt.plot(glutamine_values_0_223,biomass_reals_0_223,'b')
plt.xlabel('glutamine lower boundaries')
plt.ylabel('biomass growth')
plt.title('Glucose lb = 0.21 mM')


# glucose = 0.162181818181818 mM/hr
glutamine_values_0_162 = []
biomass_reals_0_162 = []
biomass_predictions_0_162 = []
for i in range(2700,2800):
    glutamine_value = sorted_data[i,1]
    biomass_real = sorted_data[i,3]
    biomass_predicted =  model.predict([[0.154848484848485, glutamine_value]])
  #  print(biomass_predicted)
    glutamine_values_0_162.append(glutamine_value)
    biomass_reals_0_162.append(biomass_real)
    biomass_predictions_0_162.append(biomass_predicted[0])

fig = plt.figure()
plt.plot(glutamine_values_0_162,biomass_reals_0_162,'ko')
plt.plot(glutamine_values_0_162,biomass_reals_0_162,'b')
plt.xlabel('glutamine lower boundaries')
plt.ylabel('biomass growth')
plt.title('Glucose lb = 0.155 mM')


# glucose = 0.0743333333333333 mM/hr
glutamine_values_0_074 = []
biomass_reals_0_074 = []
biomass_predictions_0_074 = []
for i in range(6600,6700):
    glutamine_value = sorted_data[i,1]
    biomass_real = sorted_data[i,3]
    biomass_predicted =  model.predict([[0.08484848, glutamine_value]])
  #  print(biomass_predicted)
    glutamine_values_0_074.append(glutamine_value)
    biomass_reals_0_074.append(biomass_real)
    biomass_predictions_0_074.append(biomass_predicted[0])

fig = plt.figure()
plt.plot(glutamine_values_0_074,biomass_reals_0_074,'ko')
plt.plot(glutamine_values_0_074,biomass_reals_0_074,'b')
plt.xlabel('glutamine lower boundaries')
plt.ylabel('biomass growth')
plt.title('Glucose lb = 0.085 mM')


# glucose = 0.0 mM/hr
glutamine_values_0 = []
biomass_reals_0 = []
biomass_predictions_0 = []
for i in range(9900,10000):
    glutamine_value = sorted_data[i,1]
    biomass_real = sorted_data[i,3]
    biomass_predicted =  model.predict([[0.0, glutamine_value]])
  #  print(biomass_predicted)
    glutamine_values_0.append(glutamine_value)
    biomass_reals_0.append(biomass_real)
    biomass_predictions_0.append(biomass_predicted[0])

fig = plt.figure()
plt.plot(glutamine_values_0,biomass_reals_0,'ko')
plt.plot(glutamine_values_0,biomass_reals_0,'b')
plt.xlabel('glutamine lower boundaries')
plt.ylabel('biomass growth')
plt.title('Glucose lb = 0.0 mM')





#%%

glutamine_value_003 = 0.003;
matching_indices = np.where(sorted_data[:,1] == glutamine_value_003)

only_glutamine_003 = sorted_data[matching_indices,:]
unique_rows_003 = np.unique(only_glutamine_003, axis=1)

glucose_concentrations = unique_rows_003[:,0]

glucose_values_003 = []
biomass_reals_003 = []
biomass_predictions_gln_003 = []

for g in range(0,np.size(unique_rows_003[0,:,0])):
    glucose_value = unique_rows_003[0,g,0]
    biomass_value = unique_rows_003[0,g,3]
    biomass_predicted = model.predict([[glucose_value,glutamine_value_003]])
    glucose_values_003.append(glucose_value)
    biomass_reals_003.append(biomass_value)
    biomass_predictions_gln_003.append(biomass_predicted[0])

fig = plt.figure()
plt.plot(glucose_values_003,biomass_reals_003,'ko')
plt.plot(glucose_values_003,biomass_predictions_gln_003,'b')
plt.xlabel('glucose lower boundaries')
plt.ylabel('biomass growth')
plt.title('glutamine = 0.003 mM/hr')





glutamine_value_002 = 0.002;
matching_indices = np.where(sorted_data[:,1] == glutamine_value_002)

only_glutamine_002 = sorted_data[matching_indices,:]
unique_rows_002 = np.unique(only_glutamine_002, axis=1)

glucose_concentrations = unique_rows_003[:,0]

glucose_values_002 = []
biomass_reals_002 = []
biomass_predictions_gln_002 = []

for g in range(0,np.size(unique_rows_002[0,:,0])):
    glucose_value = unique_rows_002[0,g,0]
    biomass_value = unique_rows_002[0,g,3]
    biomass_predicted = model.predict([[glucose_value,glutamine_value_002]])
    glucose_values_002.append(glucose_value)
    biomass_reals_002.append(biomass_value)
    biomass_predictions_gln_002.append(biomass_predicted[0])

fig = plt.figure()
plt.plot(glucose_values_002,biomass_reals_002,'ko')
plt.plot(glucose_values_002,biomass_predictions_gln_002,'b')
plt.xlabel('glucose lower boundaries')
plt.ylabel('biomass growth')
plt.title('glutamine = 0.002 mM/hr')




glutamine_value_001 = 0.001;
matching_indices = np.where(sorted_data[:,1] == glutamine_value_001)

only_glutamine_001 = sorted_data[matching_indices,:]
unique_rows_001 = np.unique(only_glutamine_001, axis=1)

glucose_concentrations = unique_rows_001[:,0]

glucose_values_001 = []
biomass_reals_001 = []
biomass_predictions_gln_001 = []

for g in range(0,np.size(unique_rows_001[0,:,0])):
    glucose_value = unique_rows_001[0,g,0]
    biomass_value = unique_rows_001[0,g,3]
    biomass_predicted = model.predict([[glucose_value,glutamine_value_001]])
    glucose_values_001.append(glucose_value)
    biomass_reals_001.append(biomass_value)
    biomass_predictions_gln_001.append(biomass_predicted[0])

fig = plt.figure()
plt.plot(glucose_values_001,biomass_reals_001,'ko')
plt.plot(glucose_values_001,biomass_predictions_gln_001,'b')
plt.xlabel('glucose lower boundaries')
plt.ylabel('biomass growth')
plt.title('glutamine = 0.001 mM/hr')



glutamine_value_000 = 0.000;
matching_indices = np.where(sorted_data[:,1] == glutamine_value_000)

only_glutamine_000 = sorted_data[matching_indices,:]
unique_rows_000 = np.unique(only_glutamine_000, axis=1)

glucose_concentrations = unique_rows_000[:,0]

glucose_values_000 = []
biomass_reals_000 = []
biomass_predictions_gln_000 = []

for g in range(0,np.size(unique_rows_000[0,:,0])):
    glucose_value = unique_rows_000[0,g,0]
    biomass_value = unique_rows_000[0,g,3]
    biomass_predicted = model.predict([[glucose_value,glutamine_value_000]])
    glucose_values_000.append(glucose_value)
    biomass_reals_000.append(biomass_value)
    biomass_predictions_gln_000.append(biomass_predicted[0])

fig = plt.figure()
plt.plot(glucose_values_000,biomass_reals_000,'ko')
plt.plot(glucose_values_000,biomass_predictions_gln_000,'b')
plt.xlabel('glucose lower boundaries')
plt.ylabel('biomass growth')
plt.title('glutamine = 0.0 mM/hr')



#%%
#save model
from keras2cpp import export_model
export_model(model, 'KRAS.model')