# first neural network with keras tutorial
from numpy import loadtxt
from keras.models import Sequential
from keras.layers import Dense
import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import r2_score
from keras2cpp import export_model
import csv
from sklearn.metrics import mean_squared_error
from scipy import stats

cell_type = 'WT'
train_model = 'Y'
save_model = 'Y'
draw_convergence = 'N'
Regression_Testing = 'N'
save_test_CSV = 'N'






# load the dataset
dataname = cell_type + '_in_silico_data_GLC_GLN_and_Seven_Metabolites_5_even.csv'

dataset = []
with open(dataname,encoding="utf-8-sig") as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=',')
    line_count = 0
    for row in csv_reader:
        # print(row[:6])
        if (row[0] != 'Glucose_exchange_input'):
            dataset.append(row[:10])
        
dataset = np.asarray(dataset,dtype=np.float32)
# dataset_prune_ind = np.where(dataset[:,2] > 6.3)
# dataset = dataset[dataset_prune_ind]
# dataset_prune_ind = np.where(dataset[:,3] > 0.55)
# dataset = dataset[dataset_prune_ind]
# dataset_prune_ind = np.where(dataset[:,4] > 0.87)
# dataset = dataset[dataset_prune_ind]

#plt.hist(dataset[:,5],bins=80)
# plt.xticks(range(10))


#%%
sorted_data = np.array(dataset)
data = np.array(dataset)
np.random.shuffle(data)

percent_training_data = 0.8
splitter = round(data.shape[0] * percent_training_data)

training_data = data[:splitter]
test_data = data[splitter:]



biomass_multiplier = 100


# split into input (X) and output (y) variables
X = training_data[:,0:9]
y = training_data[:,9]
multiplier = [100]

x_test = test_data[:,0:9]
y_test = test_data[:,9]

y = y*multiplier



if (train_model == 'Y'):
    # define the keras model
    model = Sequential()
    model.add(Dense(20, input_shape=(9,), activation='relu'))
    model.add(Dense(200, activation='relu'))
    model.add(Dense(20, activation='relu'))
    model.add(Dense(1, activation='relu'))
    # compile the keras model
    model.compile(loss='mse', optimizer='adam', metrics=['mae'])
    # fit the keras model on the dataset
    history = model.fit(X, y, epochs=50, batch_size=20000)
    history2 = model.predict(x_test)   
    mse_test = mean_squared_error(y_test, history2)

if (draw_convergence == 'Y'):
    plt.plot(history.history['loss'],'o',color='black')
    plt.plot(mse_test,'o',color='red')
    plt.title('mean absolulte error')
    plt.ylabel('mae')
    plt.xlabel('number of epochs')
    plt.show()


#%% 
# evaluate the keras model
if Regression_Testing == 'Y':
    print('Testing')
    eval_list = []
    for row in test_data:
        data_list = row.tolist()
        eval_result = model.predict([data_list[0:9]])
        eval_list.append(eval_result)

    True_values = test_data[0:len(eval_list),9]
    Predicted_values = np.asarray(eval_list)
    Predicted_values.resize(len(eval_list))
    
    slope, intercept, r_value, p_value, std_err = stats.linregress(True_values, Predicted_values)



#%% 
if (save_test_CSV == 'Y'):
    import csv
    np.savetxt("testdata_for_WT_DNN_including_intracellular_seven_metabolites.csv", test_data, delimiter=",")


#%%
#save model
if (save_model == 'Y'):
    import csv 
    model_name = cell_type + '_DNN_including_intracellular_seven_metabolites_5_even' + '.model'
    export_model(model, model_name)
    
    with open('multipliers.csv', 'w') as file:
        writer = csv.writer(file)
        writer.writerow(multiplier)