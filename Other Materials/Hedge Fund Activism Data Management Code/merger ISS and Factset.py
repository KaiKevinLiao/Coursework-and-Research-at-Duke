#!/usr/bin/env python
# coding: utf-8

# In[267]:


import pandas as pd
import numpy  as np
import os

os.chdir("C:\\Users\\kaike\\OneDrive\\Desktop\\Serena") 
df_factset = pd.read_csv('Factset_potential_match_multiple_votedfor.csv')


# In[268]:


df_ISS = pd.read_csv('ISS_potential_match_multiple_votedfor.csv', encoding='latin1')


# In[269]:


# Generate variable "DirectorName" for the cases with proposal to elect a director.
df_factset.assign(DirectorName="")
for i in range(0, len(df_factset.index)):
    standard_str = "Elect Management's Director Nominee - "
    position = str(df_factset.iloc[i]['ProposalText']).find(standard_str)
    if position != -1:
        df_factset.at[i,'Special_Proposal_Label'] = 'Elect Director'
        str_temp = str(df_factset.iloc[i]['ProposalText'])        [position + len(standard_str):len(str(df_factset.iloc[i]['ProposalText']))-1]
        if str_temp.find('Class') != -1:
            str_temp = str_temp.split(',')[1]  # Some special situations
        else:
            str_temp = str_temp.split(',')[0] # delete suffix, e.g., PhD.
        str_temp = str_temp.replace('Dr. ', '') # delete prefix Dr.
        df_factset.at[i,'DirectorName'] = str_temp


# In[270]:


df_ISS.assign(DirectorName="")
for i in range(0, len(df_ISS.index)):
    standard_str = "Elect Director "
    position = str(df_ISS.iloc[i]['ItemDesc']).find(standard_str)
    if position != -1:
        df_ISS.at[i,'Special_Proposal_Label'] = 'Elect Director'
        str_temp = str(df_ISS.iloc[i]['ItemDesc'])        [position + len(standard_str):len(str(df_ISS.iloc[i]['ItemDesc']))]
        if str_temp.find('Class') != -1:
            str_temp = str_temp.split(',')[1]
        else:
            str_temp = str_temp.split(',')[0]
        df_ISS.at[i,'DirectorName'] = str_temp


# In[271]:


for i in range(0, len(df_factset.index)):
    str_temp = str(df_factset.iloc[i]['ProposalText'])
    if ( ( str_temp.find('ratify') != -1 or str_temp.find('Ratify') != -1           or str_temp.find('ratification') != -1 or str_temp.find('Ratification') != -1 ) 
        and ( str_temp.find('accounting') != -1 or str_temp.find('audit') != -1) ):
        df_factset.at[i,'Special_Proposal_Label'] = 'Ratify Auditors'
    elif ( (str_temp.find('executive') != -1 or str_temp.find('Executive') != -1 )         and (str_temp.find('compensation') != -1 or str_temp.find('Compensation') != -1) ):
        df_factset.at[i,'Special_Proposal_Label'] = 'Executive Compensation'
    elif ( (str_temp.find('employee') != -1 or str_temp.find('Employee') != -1)           and (str_temp.find('purchase') != -1 or str_temp.find('Purchase') != -1) ):
        df_factset.at[i,'Special_Proposal_Label'] = 'Employee Stock Purchase'
    elif ( (str_temp.find('Omnibus') != -1) ):
        df_factset.at[i,'Special_Proposal_Label'] = 'Omnibus Stock Plan'


# In[272]:


for i in range(0, len(df_ISS.index)):
    str_temp = str(df_ISS.iloc[i]['ItemDesc'])
    if ( ( str_temp.find('ratify') != -1 or str_temp.find('Ratify') != -1           or str_temp.find('ratification') != -1 or str_temp.find('Ratification') != -1 ) 
        and ( str_temp.find('accounting') != -1 or str_temp.find('audit') != -1 or str_temp.find('Audit') != -1) ):
        df_ISS.at[i,'Special_Proposal_Label'] = 'Ratify Auditors'
    elif ( (str_temp.find('executive') != -1 or str_temp.find('Executive') != -1)         and str_temp.find('compensation') != -1 or str_temp.find('Compensation') != -1):
        df_ISS.at[i,'Special_Proposal_Label'] = 'Executive Compensation'
    elif ( (str_temp.find('employee') != -1 or str_temp.find('Employee') != -1)           and (str_temp.find('purchase') != -1 or str_temp.find('Purchase') != -1) ):
        df_ISS.at[i,'Special_Proposal_Label'] = 'Employee Stock Purchase'
    elif ( (str_temp.find('Omnibus') != -1) ):
        df_ISS.at[i,'Special_Proposal_Label'] = 'Omnibus Stock Plan'


# In[273]:


df_factset = df_factset.astype(str) 
for i in range(0, len(df_factset.index)):
    factset_cusip = str( df_factset.iloc[i]['CUSIP'] )
    factset_meetingdate = str( df_factset.iloc[i]['MeetingDate'] )
    factset_label = str( df_factset.iloc[i]['Special_Proposal_Label'] )
    factset_directorname = str( df_factset.iloc[i]['DirectorName'] )
    for j in range(0, len(df_ISS.index)):
        ISS_cusip = str( df_ISS.iloc[j]['CUSIP'] )
        ISS_meetingdate = str( df_ISS.iloc[j]['MeetingDate'] )
        ISS_label = str( df_ISS.iloc[j]['Special_Proposal_Label'] )
        ISS_directorname = str( df_ISS.iloc[j]['DirectorName'] )
        if ( (factset_cusip == ISS_cusip) and (factset_meetingdate == ISS_meetingdate)            and (factset_label == ISS_label) and (factset_directorname == ISS_directorname) ):
            df_factset.at[i, '_merge_from'] = j
            for column in df_factset:
                if (df_factset.iloc[i][column] == 'nan'):
                    df_factset.at[i, column] = df_ISS.iloc[j][column]               
            break
            


# In[ ]:


df_factset_adjusted = pd.read_csv('factset_adjusted3.csv')
for i in range(0, len(df_factset.index)):
    if( df_factset.at[i, '_merge_from'] != np.nan ):
        df_ISS.at[df_factset.at[i, '_merge_from'], '_merge_to'] = i


# In[266]:


df_ISS.to_csv('ISS_adjusted.csv')


# In[ ]:




