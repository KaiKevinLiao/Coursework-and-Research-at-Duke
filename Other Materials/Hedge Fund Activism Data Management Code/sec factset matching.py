#!/usr/bin/env python
# coding: utf-8

# In[68]:


import pandas as pd
import numpy  as np
import os
import re
from difflib import SequenceMatcher
from datetime import datetime

pd.set_option('display.max_columns', 500)
os.chdir("C:\\Users\\kaike\\OneDrive\\Desktop\\Serena\\2020.8.23")
df = pd.read_csv('unmatched_to_check.csv',encoding='latin1')


# In[69]:


for i in range(0,len(df.index)):
    df.at[i,'index'] = i
    df.at[i,'ticker_year_combine'] = str(df.iloc[i]['ticker_combined'])+str(df.iloc[i]['year_combined'])


# In[70]:


i = 0
word_list = ['Environmental','Environment','Sustainability','Sustainable','Renewable','Climate','Energy','Emissions']
while i < len(df.index):
    ticker_year_temp = str(df.iloc[i]['ticker_year_combine'])
    df_temp = df[df['ticker_year_combine'].str.startswith(ticker_year_temp)]
    for ii in range(0, len(df_temp.index)):
        string_temp_1 = re.sub('[^A-Za-z0-9]+', ' ', str(df_temp.iloc[ii]['proxyproposal']) )
        string_temp_2 = re.sub('[^A-Za-z0-9]+', ' ', str(df_temp.iloc[ii]['resolution']) )
        if df_temp.iloc[ii]['sec'] == 1:
            if str(df_temp.iloc[ii]['proposalcategory']).find('Environmental') != -1:
                df.at[df_temp.iloc[ii]['index'],'environment_sec'] = 1
        elif df_temp.iloc[ii]['factset_iss'] == 1:
            # print(string_temp)
            if any(ele in string_temp_1 for ele in word_list) or any(ele in string_temp_2 for ele in word_list) :
                # print('ctmd', string_temp_1, string_temp_2)
                df.at[df_temp.iloc[ii]['index'],'environment_factset'] = 1
    i = i + len(df_temp.index)


# In[71]:


i = 0
while i < len(df.index):
    ticker_year_temp = ''
    ticker_year_temp = str(df.iloc[i]['ticker_year_combine'])
    df_temp = df[df['ticker_year_combine'].str.startswith(ticker_year_temp)]
    kk = 0
    oo = 0
    for ii in range(0, len(df_temp.index)):
        if df_temp.iloc[ii]['environment_sec'] == 1:
            kk = kk +1
        if df_temp.iloc[ii]['environment_factset'] == 1:
            oo = oo +1
    # print(kk, oo)
    if oo == 1:
        for ii in range(0, len(df_temp.index)):
            if df_temp.iloc[ii]['environment_sec'] == 1:
                index_temp = df_temp.iloc[ii]['index']
        for ii in range(0, len(df_temp.index)):
            if df_temp.iloc[ii]['environment_factset'] == 1:
                df.at[df_temp.iloc[ii]['index'],'match_with'] = index_temp
    elif oo > 1:
        print(kk,oo,ticker_year_temp)
    i = i + len(df_temp.index)


# In[59]:


df2 = pd.read_csv('unmatched_to_check - Copy (2).csv',encoding='latin1')


# In[60]:


df2 = df2.astype(str) 
for i in range(0, len(df2.index)):
    #print(df2.at[i,'match_with'])
    index = df2.at[i,'match_with'].split('.', 1)[0]
    if index != 'nan':
        match_from = int(index)
        print(match_from)
        for column in df2:
            #print(df2.iloc[i][column])
            if (df2.iloc[i][column] == 'nan'):
                print(column)
                df2.at[i, column] = df2.iloc[match_from][column]
                
for i in range(0, len(df2.index)):
    for column in df2:
        if df2.at[i, column] == 'nan':
            df2.at[i, column] = ''

delete = []
for i in range(0, len(df2.index)):
    index = df2.at[i,'match_with'].split('.', 1)[0]
    if index != '':
        delete.append(int(index))

df2 = df2.drop(delete)


# In[62]:


df2.to_csv('unmatched_to_check_results3.csv')


# In[61]:


df2


# In[ ]:




