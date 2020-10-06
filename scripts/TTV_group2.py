import numpy as np
from sklearn.linear_model import LinearRegression
import matplotlib.pyplot as plt

transit_time = []
transit_number = []
TTV = []

#for i in range(1512): # don't need this
transit_time=np.loadtxt('lhs3844b_history_transit_times_sync_005.txt',usecols=1)
transit_number=np.loadtxt('lhs3844b_history_transit_times_sync_005.txt',usecols=0)
model=LinearRegression().fit(transit_number.reshape((-1,1)),transit_time)
T_initial=model.intercept_
mean_period=model.coef_
#print(T_initial,mean_period*365.25)
for j in range(len(transit_time)):
    tmp_calc=T_initial+transit_number[j]*mean_period
    tmp_obs=transit_time[j]
    TTV.append((tmp_obs-tmp_calc)*365.25*86400.0)
    TTV_each=TTV[j]
    #print(j,TTV_each)

fig = plt.figure(figsize=(10, 10))
ax = fig.add_subplot(2,2,1)
line, = ax.plot(transit_number, TTV, '-', label = "No effects")
ax.set_ylabel("TTV (s)", fontsize = 16)

transit_time = []
transit_number = []
TTV = []

#for i in range(1512): # don't need this
transit_time=np.loadtxt('lhs3844b_history_transit_times_32_005.txt',usecols=1)
transit_number=np.loadtxt('lhs3844b_history_transit_times_32_005.txt',usecols=0)
model=LinearRegression().fit(transit_number.reshape((-1,1)),transit_time)
T_initial=model.intercept_
mean_period=model.coef_
#print(T_initial,mean_period*365.25)
for j in range(len(transit_time)):
    tmp_calc=T_initial+transit_number[j]*mean_period
    tmp_obs=transit_time[j]
    TTV.append((tmp_obs-tmp_calc)*365.25*86400.0)
    TTV_each=TTV[j]
    #print(j,TTV_each)

ax = fig.add_subplot(2,2,2)
line, = ax.plot(transit_number, TTV, '-', label = "No effects")

transit_time = []
transit_number = []
TTV = []

#for i in range(1512): # don't need this
transit_time=np.loadtxt('lhs3844b_history_transit_times_21_005.txt',usecols=1)
transit_number=np.loadtxt('lhs3844b_history_transit_times_21_005.txt',usecols=0)
model=LinearRegression().fit(transit_number.reshape((-1,1)),transit_time)
T_initial=model.intercept_
mean_period=model.coef_
#print(T_initial,mean_period*365.25)
for j in range(len(transit_time)):
    tmp_calc=T_initial+transit_number[j]*mean_period
    tmp_obs=transit_time[j]
    TTV.append((tmp_obs-tmp_calc)*365.25*86400.0)
    TTV_each=TTV[j]
    #print(j,TTV_each)

ax = fig.add_subplot(2,2,3)
line, = ax.plot(transit_number, TTV, '-', label = "No effects")
ax.set_xlabel("Transit number", fontsize = 16)
ax.set_ylabel("TTV (s)", fontsize = 16)

plt.savefig('TTV_l1689b.png')
plt.show()
