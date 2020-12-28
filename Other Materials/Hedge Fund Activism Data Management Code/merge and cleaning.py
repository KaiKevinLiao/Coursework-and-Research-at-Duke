#!/usr/bin/env python
# coding: utf-8

# In[189]:


import pandas as pd
import numpy  as np
import os
import re
from difflib import SequenceMatcher
from datetime import datetime

pd.set_option('display.max_columns', 500)

def similar(a, b):
    similarity = 0
    if ( a == 'Ratify Auditors' or a == 'Pay Frequency' or a == 'Executive Compensation'        or a == 'Employee Stock Purchase' or a == 'Omnibus Stock Plan' or a == 'Merge'        or a == 'Adjourn Meeting' or a == 'Elect Director' or a == 'Other Business')        and ( b == 'Ratify Auditors' or b == 'Pay Frequency' or b == 'Executive Compensation'        or b == 'Employee Stock Purchase' or b == 'Omnibus Stock Plan' or b == 'Merge'        or b == 'Adjourn Meeting' or b == 'Elect Director' or b == 'Other Business'):
        similarity = 0
    else:
        similarity = SequenceMatcher(None, a, b).ratio()
    return similarity

def checkIfDuplicates_1(listOfElems):
    if len(listOfElems) == len(set(listOfElems)):
        return False
    else:
        return True

os.chdir("C:\\Users\\kaike\\OneDrive\\Desktop\\Serena\\2020.8.16") 
df_factset = pd.read_csv('unmatched_to_check.csv',encoding='latin1')
df_russell = pd.read_csv('ISS.csv',encoding='latin1')


# In[190]:


i = 0

##########
# Initial variable to label the position of matching by sponsors
##########
df_factset.at[i,'sponsor_match_with_russell'] = -1
df_russell.at[i,'sponsor_match_with_factset'] = -1
df_factset['sponsor_match_with_russell'] = df_factset['sponsor_match_with_russell'].astype('object')
df_russell['sponsor_match_with_factset'] = df_russell['sponsor_match_with_factset'].astype('object')

##########
# Initial variable to label the similarity of matching by proposals
##########
df_factset.at[i,'proposal_similarity_russell'] = -1
df_russell.at[i,'proposal_similarity_factset'] = -1
df_factset['proposal_similarity_russell'] = df_factset['proposal_similarity_russell'].astype('object')
df_russell['proposal_similarity_factset'] = df_russell['proposal_similarity_factset'].astype('object')

##########
# Initial variable to label the distance of meetingdate
##########
df_factset.at[i,'meetingdate_distance_russell'] = -1
df_russell.at[i,'meetingdate_distance_factset'] = -1
df_factset['meetingdate_distance_russell'] = df_factset['meetingdate_distance_russell'].astype('object')
df_russell['meetingdate_distance_factset'] = df_russell['meetingdate_distance_factset'].astype('object')

##########
# Initial variable to label the matching results
##########
df_factset.at[i,'matching_russell'] = -1
df_russell.at[i,'matching_factset'] = -1
df_factset['matching_russell'] = df_factset['matching_russell'].astype('object')
df_russell['matching_factset'] = df_russell['matching_factset'].astype('object')

for i in range(0,len(df_factset.index)):
    df_factset.at[i,'index'] = i
    df_factset.at[i,'sponsor_match_with_russell'] = [-1]
    df_factset.at[i,'proposal_similarity_russell'] = [-1]
    df_factset.at[i,'meetingdate_distance_russell'] = [-1]
    df_factset.at[i,'matching_russell'] = [-1]
for i in range(0,len(df_russell.index)):
    df_russell.at[i,'index'] = i
    df_russell.at[i,'sponsor_match_with_factset'] = [-1]
    df_russell.at[i,'proposal_similarity_factset'] = [-1]
    df_russell.at[i,'meetingdate_distance_factset'] = [-1]
    df_russell.at[i,'matching_factset'] = [-1]


# In[191]:


for i in range(0, len(df_russell.index)):
    str_temp = str(df_russell.iloc[i]['MeetingDate'])
    str_date = str_temp[0:2]
    str_month_temp = str_temp[2:5]
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
    str_year = str_temp[5:9]
    df_russell.at[i,'string_MeetingDate'] = str_year + str_month +str_date


# In[192]:


for i in range(0, len(df_russell.index)):
    if (str(df_russell.iloc[i]['ItemDesc']).find('Elect Trustee ') != -1):
        df_russell.at[i,'ItemDesc'] = str(df_russell.iloc[i]['ItemDesc']).replace('Trustee','Director')
    standard_str = "Elect Director "
    position = str(df_russell.iloc[i]['ItemDesc']).find(standard_str)
    if position != -1:
        df_russell.at[i,'Special_Proposal_Label'] = 'Elect Director'
        str_temp = str(df_russell.iloc[i]['ItemDesc'])        [position + len(standard_str):len(str(df_russell.iloc[i]['ItemDesc']))]
        if str_temp.find('Class') != -1:
            str_temp = str_temp.split(',')[1]
        else:
            str_temp = str_temp.split(',')[0]
        df_russell.at[i,'DirectorName'] = str_temp
        str_temp = re.sub('[^a-zA-Z0-9 \n]', '', str_temp)
        # print(str_temp)

        
for i in range(0, len(df_factset.index)):
    str_temp = str(df_factset.iloc[i]['ProxyProposal'])
    str_temp2 = str(df_factset.iloc[i]['ProposalText'])
    
    if ( str_temp.find('Elect Management\'s Director Nominee') != -1         and str_temp2.find('Elect Management\'s Director Nominee - ')!= -1):
        director_name = str_temp2.split(' - ')[1]
        director_name = director_name.split(',')[0]
        director_name = re.sub('[^a-zA-Z0-9 \n]', '', director_name)
        df_factset.at[i,'DirectorName'] = director_name
        df_factset.at[i, 'Special_Proposal_Label'] = 'Elect Director'
        # print(df_factset.at[i,'DirectorName'])
    
    elif ( ( str_temp.find('ratify') != -1 or str_temp.find('Ratify') != -1           or str_temp.find('ratification') != -1 or str_temp.find('Ratification') != -1 ) 
        and ( str_temp.find('accounting') != -1 or str_temp.find('audit') != -1) ):
        df_factset.at[i,'Special_Proposal_Label'] = 'Ratify Auditors'
        
    elif str_temp.find('frequency') != -1 or str_temp.find('Frequency')!= -1:
        df_factset.at[i,'Special_Proposal_Label'] = 'Pay Frequency'
        
    elif ( (str_temp.find('executive') != -1 or str_temp.find('Executive') != -1 or str_temp.find('Management') != -1)         and (str_temp.find('compensation') != -1 or str_temp.find('Compensation') != -1) ):
        df_factset.at[i,'Special_Proposal_Label'] = 'Executive Compensation'
        
    elif ( (str_temp.find('employee') != -1 or str_temp.find('Employee') != -1)           and (str_temp.find('purchase') != -1 or str_temp.find('Purchase') != -1) ):
        df_factset.at[i,'Special_Proposal_Label'] = 'Employee Stock Purchase'
        
    elif ( (str_temp.find('Omnibus') != -1) ):
        df_factset.at[i,'Special_Proposal_Label'] = 'Omnibus Stock Plan'
        
    elif ( (str_temp.find('Merge') != -1) ):
        df_factset.at[i,'Special_Proposal_Label'] = 'Merge'
        
    elif ( (str_temp.find('Adjourn') != -1) ):
        df_factset.at[i,'Special_Proposal_Label'] = 'Adjourn Meeting'
        
    elif ( (str_temp.find('Director') != -1) ):
        df_factset.at[i,'Special_Proposal_Label'] = 'Elect Director'
        
    elif ( (str_temp.find('Miscellaneous') != -1) ):
        df_factset.at[i, 'Special_Proposal_Label'] = 'Other Business'     
        
    else:
        df_factset.at[i,'Special_Proposal_Label'] = str_temp
    

    
    
for i in range(0, len(df_russell.index)):
    str_temp = str(df_russell.iloc[i]['ItemDesc'])
    
    if ( ( str_temp.find('ratify') != -1 or str_temp.find('Ratify') != -1           or str_temp.find('ratification') != -1 or str_temp.find('Ratification') != -1 ) 
        and ( str_temp.find('accounting') != -1 or str_temp.find('audit') != -1 or str_temp.find('Audit') != -1) ):
        df_russell.at[i,'Special_Proposal_Label'] = 'Ratify Auditors'
        
    elif str_temp.find('frequency') != -1 or str_temp.find('Frequency')!= -1:
        df_russell.at[i,'Special_Proposal_Label'] = 'Pay Frequency'
        
    elif ( (str_temp.find('executive') != -1 or str_temp.find('Executive') != -1)         and str_temp.find('compensation') != -1 or str_temp.find('Compensation') != -1 or str_temp.find('Bonus') != -1):
        df_russell.at[i,'Special_Proposal_Label'] = 'Executive Compensation'
        
    elif ( (str_temp.find('employee') != -1 or str_temp.find('Employee') != -1)           and (str_temp.find('purchase') != -1 or str_temp.find('Purchase') != -1) ):
        df_russell.at[i,'Special_Proposal_Label'] = 'Employee Stock Purchase'
        
    elif ( (str_temp.find('Omnibus') != -1) ):
        df_russell.at[i,'Special_Proposal_Label'] = 'Omnibus Stock Plan'
        
    elif ( (str_temp.find('Merge') != -1) ):
        df_russell.at[i,'Special_Proposal_Label'] = 'Merge'
        
    elif ( (str_temp.find('Adjourn') != -1) ):
        df_russell.at[i,'Special_Proposal_Label'] = 'Adjourn Meeting'
        
    elif ( (str_temp.find('Director') != -1) ):
        df_russell.at[i,'Special_Proposal_Label'] = 'Elect Director'
    
    elif ( (str_temp.find('Other Business') != -1) ):
        df_russell.at[i, 'Special_Proposal_Label'] = 'Other Business' 
        
    else:
        df_russell.at[i,'Special_Proposal_Label'] = str_temp


# In[193]:


threshold_similarity_proposal = 0.4
threshold_similarity_directorname = 0.5

i = 0
while i < len(df_russell.index):
    Ticker_temp = str(df_russell.iloc[i]['CUSIP'])
    df_russell_temp = df_russell[df_russell['CUSIP'].str.startswith(Ticker_temp)]
    df_factset_temp = df_factset[df_factset['CUSIP'].str.startswith(Ticker_temp)]
    for ii in range(0, len(df_factset_temp.index)):
        factset_cusip = str( df_factset_temp.iloc[ii]['CUSIP'] )
        factset_meetingdate = str( df_factset_temp.iloc[ii]['MeetingDate'] )
        factset_label = str( df_factset_temp.iloc[ii]['Special_Proposal_Label'] )
        factset_director_name = str( df_factset_temp.iloc[ii]['DirectorName'] )
        
        for jj in range(0, len(df_russell_temp.index)):
            russell_cusip = str( df_russell_temp.iloc[jj]['CUSIP'] )
            russell_meetingdate = str( df_russell_temp.iloc[jj]['MeetingDate'] )
            russell_label = str( df_russell_temp.iloc[jj]['Special_Proposal_Label'] )
            russell_director_name = str( df_russell_temp.iloc[jj]['DirectorName'] )
            
            if len(factset_director_name) < 5 or len(russell_director_name) < 5:
                if ( (factset_cusip == russell_cusip) and (factset_meetingdate == russell_meetingdate)                    and (factset_label == russell_label) ):
                    position = int(df_factset_temp.iloc[ii]['index'])
                    #print(position)
                    matching_temp = df_factset.at[position,'matching_russell']
                    matching_temp.append(int(df_russell_temp.iloc[jj]['index']))
                    df_factset.at[position,'matching_russell'] = matching_temp
                elif  (factset_cusip == russell_cusip) and (factset_meetingdate == russell_meetingdate)                    and (factset_label == 'Omnibus Stock Plan' or factset_label == 'Executive Compensation')                    and (russell_label == 'Omnibus Stock Plan' or russell_label == 'Executive Compensation'):
                    position = int(df_factset_temp.iloc[ii]['index'])
                    #print(position)
                    matching_temp = df_factset.at[position,'matching_russell']
                    matching_temp.append(int(df_russell_temp.iloc[jj]['index']))
                    df_factset.at[position,'matching_russell'] = matching_temp
                elif (factset_cusip == russell_cusip) and (factset_meetingdate == russell_meetingdate)                    and (factset_label == 'Other Business' or factset_label == 'Merge')                    and (russell_label == 'Other Business' or russell_label == 'Merge'):
                    position = int(df_factset_temp.iloc[ii]['index'])
                    #print(position)
                    matching_temp = df_factset.at[position,'matching_russell']
                    matching_temp.append(int(df_russell_temp.iloc[jj]['index']))
                    df_factset.at[position,'matching_russell'] = matching_temp
                elif (factset_cusip == russell_cusip) and (factset_meetingdate == russell_meetingdate)                    and similar(factset_label, russell_label) > threshold_similarity_proposal:
                    position = int(df_factset_temp.iloc[ii]['index'])
                    #print(position)
                    matching_temp = df_factset.at[position,'matching_russell']
                    matching_temp.append(int(df_russell_temp.iloc[jj]['index']))
                    df_factset.at[position,'matching_russell'] = matching_temp
            else:
                if ( (factset_cusip == russell_cusip) and (factset_meetingdate == russell_meetingdate)                    and (factset_label == russell_label) )                    and similar(factset_director_name, russell_director_name) > threshold_similarity_directorname:
                    position = int(df_factset_temp.iloc[ii]['index'])
                    #print(position)
                    matching_temp = df_factset.at[position,'matching_russell']
                    matching_temp.append(int(df_russell_temp.iloc[jj]['index']))
                    df_factset.at[position,'matching_russell'] = matching_temp
                    
    i = i + len(df_russell_temp.index)

list_all = []

for i in range(0, len(df_factset.index)):
    matching_temp = df_factset.at[i,'matching_russell']
    # print(i, df_factset.at[i, 'matching_russell'])
    matching_temp.pop(0)
    df_factset.at[i,'matching_russell'] = matching_temp
    print(i, df_factset.at[i, 'matching_russell'])
    list_all = list_all + df_factset.at[i, 'matching_russell']


# In[194]:


k = 0
for i in range(0, len(df_factset.index)):
    if len(df_factset.at[i, 'matching_russell']) == 0:
        print(i,df_factset.at[i, 'matching_russell'],df_factset.at[i, 'Special_Proposal_Label'])
        k = k+1
        
import collections
print(k)
print([item for item, count in collections.Counter(list_all).items() if count > 1])

# print(df_factset.at[13, 'Special_Proposal_Label'], df_factset.at[14, 'Special_Proposal_Label'], df_russell.at[13, 'Special_Proposal_Label'])


# In[195]:


for i in range(0, len(df_factset.index)):
    matching_temp = df_factset.at[i,'matching_russell']
    matching_temp_2 = []
    for element in matching_temp:
        matching_temp_2 = df_russell.at[element, 'matching_factset']
        matching_temp_2.append(i)
        df_russell.at[element, 'matching_factset'] = matching_temp_2

for i in range(0, len(df_russell.index)):
    matching_temp = df_russell.at[i,'matching_factset']
    matching_temp.pop(0)
    df_russell.at[i,'matching_factset'] = matching_temp
    print(i, df_russell.at[i,'matching_factset'])


# In[196]:


for i in range(0, len(df_russell.index)):
    if df_russell.at[i,'matching_factset']:
        matching_position = df_russell.at[i,'matching_factset'][-1]
        print(matching_position)
        for column in df_factset:
            #if column not in df_factset:
            # print(column, df_factset.iloc[int(matching_position)][column])
            if column == 'length_cusip':
                break
            df_russell.at[i,column] = df_factset.iloc[int(matching_position)][column]


# In[197]:


k = len(df_russell.index)           
for i in range(0, len(df_factset.index)):
    if len(df_factset.at[i,'matching_russell']) == 0:
        print(i)
        for column in df_factset:
            #if column not in df_factset:
            # print(column, df_factset.iloc[int(matching_position)][column])
            if column == 'length_cusip':
                break
            df_russell.at[i+k,column] = df_factset.iloc[i][column]


# In[195]:


i = 1
while i < len(df_russell.index):
    Ticker_temp = str(df_russell.iloc[i]['CUSIP'])
    df_russell_temp = df_russell[df_russell['CUSIP'].str.startswith(Ticker_temp)]
    df_factset_temp = df_factset[df_factset['CUSIP'].str.startswith(Ticker_temp)]
    for ii in range(0, len(df_factset_temp.index)):
        for jj in range(0, len(df_russell_temp.index)):
            ##########
            # write the potential position of matching
            ##########

            # print (ii, jj, proponent_label_factset, str(df_russell_temp.iloc[jj]['Proponent']))
            sponsor_match_with_russell_temp = df_factset.at[int(df_factset_temp.iloc[ii]['index']),'sponsor_match_with_russell']
            sponsor_match_with_russell_temp.append(int(df_russell_temp.iloc[jj]['index']))
            # print(sponsor_match_with_russell_temp)
            df_factset.at[int(df_factset_temp.iloc[ii]['index']),'sponsor_match_with_russell'] = sponsor_match_with_russell_temp

            sponsor_match_with_factset_temp = df_russell.at[int(df_russell_temp.iloc[jj]['index']),'sponsor_match_with_factset']
            sponsor_match_with_factset_temp.append(int(df_factset_temp.iloc[ii]['index']))
            # print(sponsor_match_with_russell_temp)
            df_russell.at[int(df_russell_temp.iloc[jj]['index']),'sponsor_match_with_factset'] = sponsor_match_with_factset_temp

    for ii in range(0, len(df_factset_temp.index)):
        proposal_factset = str(df_factset_temp.iloc[ii]['Special_Proposal_Label'])
        meetingdate_factset = str(df_factset_temp.iloc[ii]['string_MeetingDate'])
        if proposal_factset != '':
            sponsor_match_with_russell = str(df_factset_temp.iloc[ii]['sponsor_match_with_russell'])
            if len(sponsor_match_with_russell) > 1:
                for kk in range(1, len(sponsor_match_with_russell)):

                    ##########
                    # calculate the similarity of proposal
                    ##########

                    # print(sponsor_match_with_russell[kk])
                    # print(df_russell.loc[int(sponsor_match_with_russell[kk])]['ProposalCategory'])
                    proposal_russell = str(df_russell.loc[int(sponsor_match_with_russell[kk])]['Special_Proposal_Label'])
                    print(proposal_russell)
                    # proposal_russell = re.sub('[^A-Za-z0-9]+', ' ', proposal_russell)
                    # print(proposal_russell)
                    # print(similar(proposal_factset,proposal_russell), len(proposal_russell))
                    if proposal_factset == ''
                    similarity = similar(proposal_factset,proposal_russell) * (int((len(proposal_russell)+len(proposal_russell))/80)+1)
                    # print(similarity)

                    ##########
                    # write the similarity of proposal
                    ##########

                    proposal_similarity_russell_temp = df_factset.at[int(df_factset_temp.iloc[ii]['index']),'proposal_similarity_russell']
                    proposal_similarity_russell_temp.append(similarity)
                    # print(proposal_similarity_russell_temp)
                    df_factset.at[int(df_factset_temp.iloc[ii]['index']),'proposal_similarity_russell'] = proposal_similarity_russell_temp

                    proposal_similarity_factset_temp = df_russell.at[int(sponsor_match_with_russell[kk]),'proposal_similarity_factset']
                    proposal_similarity_factset_temp.append(similarity)
                    # print(sponsor_match_with_russell_temp)
                    df_russell.at[int(sponsor_match_with_russell[kk]),'proposal_similarity_factset'] = proposal_similarity_factset_temp

                    ########
                    # calculate the distance of meetingdate
                    ########

                    meetingdate_russell = str(df_russell.loc[int(sponsor_match_with_russell[kk])]['string_MeetingDate'])

                    # print('russell', meetingdate_russell)
                    if len(meetingdate_russell) > 5 and len(meetingdate_factset) > 5:
                        date_format = "%Y%m%d"
                        # print('factset', meetingdate_factset)
                        # print('russell', meetingdate_russell)
                        datetime_factset = datetime.strptime(meetingdate_factset, date_format)
                        datetime_russell = datetime.strptime(meetingdate_russell, date_format)
                        diff = datetime_factset - datetime_russell
                        meetingdate_diff = int(diff.days)
                        # print(meetingdate_diff)
                    else:
                        meetingdate_diff = -999

                    ########
                    # write the distance of meetingdate
                    ########

                    meetingdate_russell_temp = df_factset.at[int(df_factset_temp.iloc[ii]['index']),'meetingdate_distance_russell']
                    meetingdate_russell_temp.append(meetingdate_diff)
                    # print(meetingdate_russell_temp)
                    df_factset.at[int(df_factset_temp.iloc[ii]['index']),'meetingdate_distance_russell'] = meetingdate_russell_temp

                    meetingdate_factset_temp = df_russell.at[int(sponsor_match_with_russell[kk]),'meetingdate_distance_factset']
                    meetingdate_factset_temp.append(meetingdate_diff)
                    # print(sponsor_match_with_russell_temp)
                    df_russell.at[int(sponsor_match_with_russell[kk]),'meetingdate_distance_factset'] = meetingdate_factset_temp
    print(i)
    i = i + len(df_russell_temp.index)


# In[150]:


##########
# Select appopriate matches
##########

threshold_proposal = 0.5
threshold_meetingdate = 1
for i in range(0,len(df_russell.index)):
    sponsor_matching = df_russell.iloc[i]['sponsor_match_with_factset']
    proposal_similarity = df_russell.iloc[i]['proposal_similarity_factset']
    meetingdate_difference = df_russell.iloc[i]['meetingdate_distance_factset']
    # print(sponsor_match_with_russell)
    if len(sponsor_matching) = 2:
        matching_temp = []
        matching_temp.append(sponsor_matching[1])
        df_russell.at[i,'matching_factset'] = matching_temp
    elif len(sponsor_matching) > 1:
        matching_temp = []
        for kk in range(1, len(sponsor_matching)):
            if proposal_similarity[kk] > threshold_proposal             and meetingdate_difference[kk] < threshold_meetingdate             and meetingdate_difference[kk] > -threshold_meetingdate:
                matching_temp.append(sponsor_matching[kk])
        # print(matching_temp)
        df_russell.at[i,'matching_factset'] = matching_temp


# In[200]:


df_russell.to_csv('ISS_adjusted10.csv')


# In[198]:


df_factset.to_csv('factset_adjusted10.csv')


# In[ ]:




