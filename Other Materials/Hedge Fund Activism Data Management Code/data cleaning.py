#!/usr/bin/env python
# coding: utf-8

# In[177]:


import pandas as pd
import numpy  as np
import os
import re
from difflib import SequenceMatcher
from datetime import datetime

pd.set_option('display.max_columns', 500)

df_rus = pd.read_csv('intelligize_russell3000.csv',encoding='latin1')
df_act = pd.read_csv('no_action_v2.csv',encoding='latin1')


# In[178]:


for i in range(len(df_act.index)):
    if str(df_act.at[i,'v58']) != 'nan':
        date_temp = df_act.at[i,'v58']
        year = int(int(date_temp)/10000)
        month = int((int(date_temp)-int(year)*10000)/100)
        if month < 10:
            str_month = '0' + str(month)
        else:
            str_month = str(month)
        df_act.at[i,'data_formerge'] = str(df_act.at[i,'ticker_combined'])+ str(year)
    else:
        date_temp = str(df_act.at[i,'date'])
        year = '20' + date_temp[-2:]
        str_month_temp = ''.join([i for i in date_temp if not i.isdigit()])
        str_month_temp = ''.join([e for e in str_month_temp if e.isalnum()])
        if(str_month_temp == 'jan'):
            str_month = '01'
        elif(str_month_temp == 'feb'):
            str_month = '02'
        elif(str_month_temp == 'mar'):
            str_month = '03'
        elif(str_month_temp == 'apr'):
            str_month = '04'
        elif(str_month_temp == 'may'):
            str_month = '05'
        elif(str_month_temp == 'jun'):
            str_month = '06'
        elif(str_month_temp == 'jul'):
            str_month = '07'
        elif(str_month_temp == 'aug'):
            str_month = '08'
        elif(str_month_temp == 'sep'):
            str_month = '09'
        elif(str_month_temp == 'oct'):
            str_month = '10'
        elif(str_month_temp == 'nov'):
            str_month = '11'
        elif(str_month_temp == 'dec'):
            str_month = '12'
        df_act.at[i,'data_formerge'] = str(df_act.at[i,'ticker_combined']) + year + str_month
        
for i in range(len(df_act.index)):
    if str(df_act.at[i,'proxyproposal']) != 'nan':
        df_act.at[i,'proposal_formerge'] = df_act.at[i,'proxyproposal']
    else:
        df_act.at[i,'proposal_formerge'] = df_act.at[i,'resolution']
    
        


# In[179]:


for i in range(len(df_rus.index)):
    date_temp = str(df_rus.at[i, 'date'])
    date_temp = ''.join([e for e in date_temp if e.isalnum()])
    month = date_temp[0:2]
    year = date_temp[4:]
    ticker = df_rus.at[i,'ticker']
    ticker_split = ticker.split("\t")
    ticker = ticker_split[0]
    df_rus.at[i,'data_formerge'] = str(ticker) + str(year)


# In[180]:


for i in range(len(df_act.index)):
    str_temp = str(df_act.at[i,'proposal_formerge']).lower()
    if ( ( str_temp.find('ratify') != -1 or str_temp.find('ratify') != -1           or str_temp.find('ratification') != -1 or str_temp.find('ratification') != -1 ) 
        and ( str_temp.find('accounting') != -1 or str_temp.find('audit') != -1 or str_temp.find('audit') != -1) \
       or ( str_temp.find('auditor') != -1) ):
        df_act.at[i,'special_proposal_label'] = 'ratify auditors'
        
    elif ( (str_temp.find('audit') != -1 ) ):
        df_rus.at[i,'special_proposal_label'] = 'request for audit'
        
    elif str_temp.find('frequency') != -1 or str_temp.find('frequency')!= -1:
        df_act.at[i,'special_proposal_label'] = 'pay frequency'
        
    elif str_temp.find('pay') != -1 :
        df_act.at[i,'special_proposal_label'] = 'other pay related'
        
    elif ( (str_temp.find('executive') != -1 or str_temp.find('executive') != -1)         and str_temp.find('compensation') != -1 or str_temp.find('compensation') != -1 or str_temp.find('bonus') != -1):
        df_act.at[i,'special_proposal_label'] = 'executive compensation'
        
    elif ( (str_temp.find('employee') != -1 or str_temp.find('employee') != -1)           and (str_temp.find('purchase') != -1 or str_temp.find('purchase') != -1) ):
        df_act.at[i,'special_proposal_label'] = 'employee stock purchase'
        
    elif ( (str_temp.find('omnibus') != -1) ):
        df_act.at[i,'special_proposal_label'] = 'omnibus stock plan'
        
    elif ( (str_temp.find('merge') != -1) or (str_temp.find('poison pill') != -1) or (str_temp.find('proxy') != -1)):
        df_act.at[i,'special_proposal_label'] = 'merge'
    
    elif ( (str_temp.find('vote') != -1) ):
        df_act.at[i,'special_proposal_label'] = 'vote policy change'
        
    elif ( (str_temp.find('adjourn') != -1) ):
        df_act.at[i,'special_proposal_label'] = 'adjourn meeting'
        
    elif ( (str_temp.find('director') != -1) ):
        df_act.at[i,'special_proposal_label'] = 'elect director'
    
    elif ( (str_temp.find('other business') != -1) ):
        df_act.at[i, 'special_proposal_label'] = 'other business' 
    
    elif ( (str_temp.find('social') != -1) or (str_temp.find('human') != -1)):
        df_act.at[i, 'special_proposal_label'] = 'social issues' 
        
    elif ( (str_temp.find('enviornment') != -1) or (str_temp.find('energy') != -1) or (str_temp.find('recycling') != -1)):
        df_act.at[i, 'special_proposal_label'] = 'enviornmental issues' 
        
    elif ( (str_temp.find('animal') != -1) ):
        df_act.at[i, 'special_proposal_label'] = 'animal issues' 
        
    elif ( (str_temp.find('separate') != -1) or (str_temp.find('independent') != -1)):
        df_act.at[i, 'special_proposal_label'] = 'independent chairman' 
        
    elif ( (str_temp.find('requirements') != -1) or (str_temp.find('requirement') != -1)):
        df_act.at[i, 'special_proposal_label'] = 'requirements change' 
        
    elif ( (str_temp.find('health') != -1) ):
        df_act.at[i, 'special_proposal_label'] = 'health issues' 
    
    elif ( (str_temp.find('board') != -1) ):
        df_act.at[i, 'special_proposal_label'] = 'other board related'  
    
    elif ( (str_temp.find('governance') != -1) ):
        df_rus.at[i, 'special_proposal_label'] = 'governance related' 
    
    elif( (str_temp.find('cumulative') != -1) or (str_temp.find('major') != -1) or (str_temp.find('dual class structure') != -1)):
        df_act.at[i, 'special_proposal_label'] = 'cumulative voting'
        
    elif ( (str_temp.find('governance') != -1) or (str_temp.find('meeting') != -1)):
        df_rus.at[i, 'special_proposal_label'] = 'governance related' 
        
    else:
        df_act.at[i,'special_proposal_label'] = str_temp


# In[181]:


for i in range(len(df_act.index)):
    print(df_act.at[i,'proposal_formerge'])
    print(df_act.at[i,'special_proposal_label'])
    print("")


# In[182]:


for i in range(len(df_rus.index)):
    str_temp = str(df_rus.at[i,'proposalcategory']).lower()
    if ( ( str_temp.find('ratify') != -1 or str_temp.find('ratify') != -1           or str_temp.find('ratification') != -1 or str_temp.find('ratification') != -1 ) 
        and ( str_temp.find('accounting') != -1 or str_temp.find('audit') != -1 or str_temp.find('audit') != -1) \
       or ( str_temp.find('auditor') != -1) or ( str_temp.find('auditor') != -1) or ( str_temp.find('rotation') != -1)):
        df_rus.at[i,'special_proposal_label'] = 'ratify auditors'
    
    elif ( (str_temp.find('audit') != -1 ) ):
        df_rus.at[i,'special_proposal_label'] = 'request for audit'
    
    elif str_temp.find('pay') != -1 :
        df_act.at[i,'special_proposal_label'] = 'other pay related'      
    
    elif ( str_temp.find('frequency') != -1 or str_temp.find('frequency')!= -1 or str_temp.find('pay') != -1) :
        df_rus.at[i,'special_proposal_label'] = 'pay frequency'
        
    elif ( (str_temp.find('executive') != -1 or str_temp.find('executive') != -1) or (str_temp.find('compensation') != -1)):
        df_rus.at[i,'special_proposal_label'] = 'executive compensation'
        
    elif ( (str_temp.find('employee') != -1 or str_temp.find('employee') != -1)           and (str_temp.find('purchase') != -1 or str_temp.find('purchase') != -1) ):
        df_rus.at[i,'special_proposal_label'] = 'employee stock purchase'
        
    elif ( (str_temp.find('omnibus') != -1) ):
        df_rus.at[i,'special_proposal_label'] = 'omnibus stock plan'
        
    elif ( (str_temp.find('merge') != -1) or (str_temp.find('poison pill') != -1) or (str_temp.find('proxy') != -1) or (str_temp.find('vesting') != -1)):
        df_rus.at[i,'special_proposal_label'] = 'merge'
        
    elif ( (str_temp.find('adjourn') != -1) ):
        df_rus.at[i,'special_proposal_label'] = 'adjourn meeting'
    
    elif ( (str_temp.find('vote') != -1) ):
        df_act.at[i,'special_proposal_label'] = 'vote policy change'
        
    elif ( (str_temp.find('director') != -1) and (str_temp.find('candidate') != -1)):
        df_rus.at[i,'special_proposal_label'] = 'elect director'
    
    elif ( (str_temp.find('other business') != -1) ):
        df_rus.at[i, 'special_proposal_label'] = 'other business' 
        
    elif ( (str_temp.find('enviornment') != -1) or (str_temp.find('energy') != -1) or (str_temp.find('recycling') != -1)):
        df_rus.at[i, 'special_proposal_label'] = 'enviornmental issues' 
        
    elif ( (str_temp.find('animal') != -1) ):
        df_rus.at[i, 'special_proposal_label'] = 'animal issues' 
        
    elif ( (str_temp.find('separate') != -1) or (str_temp.find('independent') != -1)):
        df_rus.at[i, 'special_proposal_label'] = 'independent chairman' 
        
    elif ( (str_temp.find('requirements') != -1) or (str_temp.find('requirement') != -1)):
        df_rus.at[i, 'special_proposal_label'] = 'requirements change' 
        
    elif ( (str_temp.find('health') != -1) ):
        df_rus.at[i, 'special_proposal_label'] = 'health issues' 
    
    elif ( (str_temp.find('board') != -1) ):
        df_rus.at[i, 'special_proposal_label'] = 'other board related'
        
    elif ( (str_temp.find('social') != -1) or (str_temp.find('human') != -1)):
        df_rus.at[i, 'special_proposal_label'] = 'social issues' 
        
    elif( (str_temp.find('vote') != -1) ):
        df_rus.at[i, 'special_proposal_label'] = 'cumulative voting'
        
    elif ( (str_temp.find('governance') != -1) ):
        df_rus.at[i, 'special_proposal_label'] = 'governance related' 
        
    elif( (str_temp.find('political') != -1) or (str_temp.find('lobby') != -1) ):
        df_rus.at[i, 'special_proposal_label'] = 'political'
        
    else:
        df_rus.at[i,'special_proposal_label'] = str_temp


# In[183]:


for i in range(len(df_rus.index)):
    print(df_rus.at[i,'proposalcategory'])
    print(df_rus.at[i,'special_proposal_label'])
    print("")


# In[184]:


for i in range(0, len(df_act.index)):
    df_act.at[i,'index'] = i


# In[185]:


i = 0
while i < len(df_rus.index):
    Ticker_temp = str(df_rus.iloc[i]['data_formerge'])
    df_rus_temp = df_rus.loc[df_rus['data_formerge']==Ticker_temp]
    df_act_temp = df_act.loc[df_act['data_formerge']==Ticker_temp]
    for ii in range(0, len(df_rus_temp.index)):
        rus_label = str( df_rus_temp.iloc[ii]['special_proposal_label'] )
        for jj in range(0, len(df_act_temp.index)):
            act_label = str( df_act_temp.iloc[jj]['special_proposal_label'] )
            if str(rus_label) == str(act_label):
                df_act.at[int(df_act_temp.iloc[jj]['index']),'mergefrom'] = df_rus_temp.iloc[ii]['counter']-1
                
    i = i+len(df_rus_temp.index)


# In[186]:


Ticker_temp = str(df_rus.iloc[100]['data_formerge'])
df_rus_temp = df_rus.loc[df_rus['data_formerge']==Ticker_temp]
df_act_temp = df_act.loc[df_act['data_formerge']==Ticker_temp]


# In[187]:


for i in range(len(df_rus.index)):
    print(df_rus.at[i,'proposal'])


# In[188]:


df_rus = df_rus.applymap(str)
df_act = df_act.applymap(str)
for i in range(0, len(df_act.index)):
    if str(df_act.at[i,'mergefrom']) != 'nan':
        print(i)
        matching_position = df_act.at[i,'mergefrom']
        df_rus.at[int(float(matching_position)),'merged'] = 1
        # print(matching_position)
        for column in df_rus:
            #if column not in df_factset:
            # print(column, df_factset.iloc[int(matching_position)][column])
            df_act.at[i,column+'_rus'] = df_rus.iloc[int(float(matching_position))][column]


# In[189]:


k = len(df_act.index)           
for i in range(0, len(df_rus.index)):
    if df_rus.at[i,'merged'] == 1:
        print(str(df_rus.at[i,'merged']))
    if df_rus.at[i,'merged'] != 1:
        # print(str(df_rus.at[i,'merged']))
        for column in df_rus:
            #if column not in df_factset:
            # print(column, df_factset.iloc[int(matching_position)][column])
            df_act.at[i+k,column+'_rus'] = df_rus.iloc[i][column]


# In[190]:


df_act = df_act.replace('nan', np.nan)


# In[86]:


df_rus.to_csv('rus_adjusted2.csv')


# In[191]:


df_act.to_csv('act_adjusted8.csv')


# In[207]:


df_act['data_formerge']


# In[66]:


df_rus_temp


# In[68]:


df_act_temp


# In[ ]:




