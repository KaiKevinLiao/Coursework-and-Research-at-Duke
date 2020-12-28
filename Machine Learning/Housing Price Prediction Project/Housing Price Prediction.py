#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from statistics import mean 

pd.set_option('display.max_columns', None)


# In[2]:


# Import Data
data_df = pd.read_csv('sf_data.csv', header = None)


# # 1 Data Cleaning
# ##  1.1 Feature Space
# 
# Each data set contains the following Features (in order):
# 
# - House ID \#
# - Price (deflated to year 2000 dollars)
# - county ID # (see below)
# - Year Built
# - Square Footage
# - \# Bathrooms
# - \# Bedrooms
# - \# Total Rooms
# - \# Stories
# - Violent Crime Rate (Cases per 100,000)
# - Property Crime Rate (Cases Per 100,000 )
# - Year of Sale (1993-2008)

# In[3]:


# Rename the columns
data_df.rename(columns={0: 'House_ID', 1: 'Price', 2: 'County_ID', 
                       3: 'Year_Built', 4: 'Square_Footage', 5: 'Bathrooms', 
                       6: 'Bedrooms', 7: 'Rooms', 8: 'Stories', 
                       9: 'VC', 10: 'PC', 11: 'Year_of_Sale'}, inplace=True)


# ##  1.2 Construct Dummies for Counties
# - $1 \quad$  Alameda
# - $13 \quad$ Contra costa
# - $75 \quad$ San Francisco
# - $81 \quad$ San Mateo
# - $85 \quad$ Santa Clara 
# 
# For San Francisco, I include dummies for counties 13, 75, 81, and 85. Alameda is the baseline.

# In[4]:


for label, county in enumerate(set(data_df.iloc[:]['County_ID'])):
    for i in range(len(data_df.index)):
        if data_df.iloc[i]['County_ID'] == county:
            data_df.at[i,str(county)] = 1


# In[5]:


data_df = data_df.fillna(0)


# In[6]:


data_df.to_csv('data_modified.csv')


# ## 1.3 Discriptive Statistic

# In[7]:


from scipy.stats import norm
from sklearn.neighbors import KernelDensity
from matplotlib.pyplot import MultipleLocator
from pylab import rcParams
from matplotlib.ticker import FuncFormatter
rcParams['figure.figsize'] = 16, 9
plt.rcParams.update({'font.size': 22})


# In[12]:


plt.hist(y, bins=150, range = (1,2000000), density = True,
         label = 'Price')
plt.ylabel('Density')
plt.xlabel('House Prices')
plt.title('Distribution of House Prices')

x_major_locator=MultipleLocator(500000)
ax=plt.gca()
ax.xaxis.set_major_locator(x_major_locator)
def formatnum(x, pos):
    return '$%.1f$x$10^{4}$' % (x*1000000)
def formatnum_x(x, pos):
    return '$%.2f \pi$' % (x)
formatter1 = FuncFormatter(formatnum)
formatter2 = FuncFormatter(formatnum_x)
ax.yaxis.set_major_formatter(formatter1)
plt.savefig('Distribution of House Prices.png')
plt.show()


# #  2 Models
# The model assumption and economic theoretical base for the model is given in the write-up.
# ## 2.1 Linear Regression
# ### 2.1.1 Construct Variables for Regression
# 
# The linear regression model 
# 
# - Constant
# - \# Bathrooms
# - \# Bedrooms
# - \# Stories
# - Property Crime Rate
# - (Property Crime Rate) $^{2}$
# - Year Built
# - (Year Built) $^{2}$
# - Square Footage
# - (Square Footage) $^{2}$
# - \# Total Rooms
# - (# Total Rooms) $^{2}$
# - Violent Crime Rate 
# - (Violent Crime Rate) $^{2}$ 
# - Vector of Year Dummies (omit 1999) 
# - Vector of Dummies for Certain Counties $^{2}$
# 
# ## 2.1.2 Conduct Linear Regression on Training Set
# 
# The first two parts are written in my Fortran code. They following python code aim to read the regression result and make prediction on testing set.
# 
# ## 2.1.3 Make Prediction on Testing Set

# In[9]:


# Read dataset
X_df = pd.read_csv('X_SF.csv')
y_df = pd.read_csv('y_SF.csv', header = None)
y_df.rename(columns={0: 'Price'}, inplace=True)


# In[10]:


# Read regression results
reg_result = pd.read_csv('Hedonic_Price_Function_Regression_Results_SF.csv')


# In[11]:


X = X_df.to_numpy()
y = y_df.to_numpy()
beta = reg_result.to_numpy()[0][:]

# Construct Testing Set (10% of the all samples)
X_test = X[340426:][:]
y_test = y[340426:][:]


# In[14]:


# Make predictions
y_pred = np.matmul(X_test, beta.T)


# ## 2.2 Keras Regressions
# ### 2.2.1 Seperate Dataset

# In[64]:


from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense, Activation
from tensorflow.keras.optimizers import Adam
from tensorflow.keras.wrappers.scikit_learn import KerasClassifier
from sklearn.model_selection import GridSearchCV
from sklearn import preprocessing
from sklearn.model_selection import KFold


# In[59]:


# Nomalized the features

X = preprocessing.scale(X)


# In[55]:


X_test = X[340426:][:]
y_test = y[340426:]
X_train = X[:340426][:]
y_train = y[:340426]


# ### 2.2.2 Train the Neural Network Model with Cross Validation
# I cross validate the number of epochs times. 

# In[ ]:


def creat_model(layer = 5):
    clt_nn = Sequential()
    for i in range(layer):
        clt_nn.add(Dense(33,activation='relu'))
    clt_nn.add(Dense(1))
    clt_nn.compile(optimizer='Adam',loss='mse')
    return clt_nn


# In[104]:


X_train_df = pd.DataFrame(data=X_train)
y_train_df = pd.DataFrame(data=y_train)


# In[106]:


EPOCH_MAX = 100
epochs = list(range(1,EPOCH_MAX+1))
val_loss = []
val_loss_mean = []

for train_index, test_index in kf.split(X_train):
    X_train_cross, X_test_cross = X_train_df.loc[train_index], X_train_df.loc[test_index]
    y_train_cross, y_test_cross = y_train_df.loc[train_index], y_train_df.loc[test_index]
    clt_nn = creat_model(5)
    history = clt_nn.fit(x=X_train_cross,y=y_train_cross, 
          validation_data=(X_test_cross,y_test_cross),
          batch_size=256,epochs=EPOCH_MAX)
    val_loss.append(history.history['val_loss'])
val_loss_mean.append(np.mean(val_loss, axis=0))


# In[107]:


val_loss_mean


# In[115]:


# find the best epoch
best_epoch = np.argmin(val_loss_mean[0])+1
best_epoch


# In[117]:


# Train the model again on the entire testing set with the best parameters.
clt_nn = creat_model(5)
clt_nn.fit(x=X_train,y=y_train, 
          validation_data=(X_test,y_test),
          batch_size=256,epochs=best_epoch)


# ### 2.2.3 Make Prediction on Testing Set

# In[118]:


y_pred_nn = clt_nn.predict(X_test)


# ## 2.3 Decision Trees
# ### 2.3.1 Seperate Dataset

# In[119]:


from sklearn.model_selection import GridSearchCV
from sklearn.tree import DecisionTreeRegressor


# ### 2.3.2 Train Model

# In[120]:


params_grid = [{'max_depth':[5,10,15,20,25]}]

clf_dt_grid = GridSearchCV(DecisionTreeRegressor(), params_grid, cv=5)
clf_dt_grid.fit(X_train, y_train)


# In[121]:


clf_dt = clf_dt_grid.best_estimator_
print('Best max_depth:',clf_dt_grid.best_estimator_.max_depth,"\n") 


# ### 2.3.3 Make Prediction on Testing Set

# In[122]:


y_pred_dt = clf_dt.predict(X_test)


# ## 2.4 Random Forests
# ### 2.4.1 Seperate Dataset

# In[123]:


from sklearn.model_selection import GridSearchCV
from sklearn.ensemble import RandomForestRegressor


# ### 2.3.2 Train Model

# In[124]:


params_grid = [{'max_depth':[5,10,15,20,25]}]

clf_rf_grid = GridSearchCV(RandomForestRegressor(), params_grid, cv=5)
clf_rf_grid.fit(X_train, y_train)


# In[125]:


clf_rf = clf_rf_grid.best_estimator_
print('Best max_depth:',clf_rf_grid.best_estimator_.max_depth,"\n") 


# ### 2.3.3 Make Prediction on Testing Set

# In[126]:


y_pred_rf = clf_rf.predict(X_test)


# # 3 Evaluation

# In[127]:


# Number of Algorithms (Prediction using the mean price is included)
ALGORITM_NUMBER = 5
# Name of Algorithms
ALGORITHM = ['Mean Price', 'Linear Regression', 'Neural Network', 'Decision Trees', 'Random Forest']
# Color of Plot
COLOR_OF_PLOT = ['y','b','r', 'g','m']

# Calculate the differences between true values of prices and our predictions.
y_result = np.zeros((len(y_pred), 4 + ALGORITM_NUMBER))

# Calculate the mean price
sum = 0
for i in range(len(y_pred)):
    sum = sum + y_test[i]
y_mean = sum/len(y_pred)

# Calculate the results of predictions from different methods
for i in range(len(y_pred)):
    y_result[i][0] = y_test[i]
    y_result[i][1] = y_pred[i]
    y_result[i][2] = y_pred_nn[i]
    y_result[i][3] = (y_mean - y_test[i]) / y_test[i]
    y_result[i][4] = (y_pred[i] - y_test[i]) / y_test[i]
    y_result[i][5] = (y_pred_nn[i] - y_test[i]) / y_test[i]
    y_result[i][6] = (y_pred_dt[i] - y_test[i]) / y_test[i]
    y_result[i][7] = (y_pred_rf[i] - y_test[i]) / y_test[i]

y_result_df = pd.DataFrame(y_result)


# In[131]:


# Plot the distribution of error

for alg in range(ALGORITM_NUMBER):
    label_plt = ALGORITHM[alg]
    color_plt = COLOR_OF_PLOT[alg]
    plt.hist(y_result_df.iloc[:][alg + 3], bins=150, range = (-1,2), density = True,
             label = label_plt, color = color_plt)
    plt.ylabel('Density')
    plt.xlabel('The Error Rate')
    plt.title(str(label_plt) + ' Prediction Error')
    plt.show()

for alg in range(ALGORITM_NUMBER):
    label_plt = ALGORITHM[alg]
    color_plt = COLOR_OF_PLOT[alg]
    if alg == 0:
        plt.hist(y_result_df.iloc[:][alg + 3], bins=150, range = (-1,2), density = True,
                 label = label_plt, color = color_plt)
    else:
        plt.hist(y_result_df.iloc[:][alg + 3], bins=150, range = (-1,2), density = True,
                 alpha=0.7, label = label_plt, color = color_plt)
plt.ylabel('Density')
plt.legend()
plt.title('Camparison of The Error Rate')
plt.show()


# In[129]:


# Calcualte and plot the accuracy rate

critiria = []
accuracy = []*(ALGORITM_NUMBER)
for i in range(1,21):
    critiria.append(i*0.05)
    

for alg in range(ALGORITM_NUMBER):
    accuracy_temp = []
    for cri in critiria:
        correct_count = 0
        for i in range(len(y_result_df.index)):
            if y_result[i][alg+3] < cri and y_result[i][3] > -cri:
                correct_count = correct_count + 1
        accuracy_temp.append(correct_count/len(y_result_df.index))
    accuracy.append(accuracy_temp)
    


# In[132]:


for alg in range(ALGORITM_NUMBER):
    color_plt = COLOR_OF_PLOT[alg]
    label_plt = ALGORITHM[alg]
    plt.plot(critiria,accuracy[alg][:],color_plt,label=label_plt, marker='.')

plt.ylim(top=1.2) 
plt.ylim(bottom=0) 


for alg in range(ALGORITM_NUMBER):
    counter = 0
    for x,y in zip(critiria,accuracy[alg][:]):
        label = "{:.2f}".format(y)
        if counter % 4 == 3:
            plt.annotate(label, (x,y), textcoords="offset points", xytext=(0,10), ha='center', fontsize=7) 
        counter = counter + 1

plt.legend()
plt.xlabel('alpha')
plt.ylabel('Accuracy Rate')
plt.show()


# In[45]:


ALGORITM_NUMBER = 2

critiria = []
accuracy = []*(ALGORITM_NUMBER+1)


# In[49]:


accuracy[0][:]


# In[27]:


loss_df = pd.DataFrame(model.history.history)
loss_df.plot(figsize=(12,8))


# In[36]:


y_result


# In[ ]:




