import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import random
from mpl_toolkits.basemap import Basemap

"""
This is a program that takes in a list of properly formatted CSV files(toAdd), 
to append to the Global Seawater d18O database (infile) using Pandas.

This program also has robust plotting capacities using Basemaps, and 
dynamic querrying based on any of the 11 fields, including seasonality,
reference, and even location.

At this time, this program is intended to update these files...

Ref                 File Name               # of Recs       Compiler    Dates
Strauss 2015        'Strauss2015.csv'       271 records     L           7.7.15
Cases 2012          'CASES-Database.csv'    485 records     C   
COX 2010,2001       'Cox-Database.csv'      88  records     C
Conroy 2012         'Conroy2012.csv'        191 records     S

Written by Leon Yin for NASA GISS NYCRI 2015
Please consult NASA GISS or the author for use.

v - 1.0
"""


filein      ='dbO18.txt' #this file is downloaded from the d18o website.
toAdd       =('Strauss2015.csv','CASES-Database.csv','Cox-Database.csv',
              'Conroy2012.csv')
fileout     ='dbO18_2015.txt'

# Logic switches -- use these guys to turn functions on-and-off.
convert         =   2       #set to 1: convert the database (DB) into a CSV, 
                            #set to 2: convert CSV input to DB
append2master   =   1       #set to 1: to add new CSV input to the DB
printOut        =  -1       #set to 1: prints out DB
plot            =  -1       #set to 1: plots a scatter for the DB
search          =  -1       #set to 1: search specfic querry
seasons         =   1       #set to 1: seasonal distribution of observtions
# Global Seawaster d18O database format
col=['lon','lat','month','year','depth','temp', 'sal', 'd18o','dD','notes','ref'] 
fWid = ( 7,    6,      2,     4,      5,     6,     6,      6,   6,    15,   60)

# Read in Database into a Pandas Dataf rame data structure
df=pd.read_fwf(filein,widths=fWid,header=None,index_col=False,names=col,skip_blank_lines=True)

# some variables...
dfRows    = len(df.index-1)
files2Add = len(toAdd)

# Convert DF to to CVS
if (convert == 1):
    toCSV=pd.DataFrame.to_csv(df,sep=',',line_terminator='\n')

# Automated reading of each CSV, and appending to existing DB.
# To troubleshoot adjust the range and double-check 'csv' in the shell.
if (convert == 2):
    for i in range(0,files2Add):
        csv=pd.read_csv(toAdd[i],sep=',',skiprows=1,names=col,na_filter=True)
        #              File Name,Delimiter, skip header
        if (append2master == 1):
            df = df.append(csv)     #add the contents of new csv into DB
            dfRows= len(df.index-1) #update the number of records in the DB

# Print the output dataframe
if(printOut==1):
    df=df.fillna('')
    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")        
    print "Opening the file..."
    dbOut = open(fileout, 'w')    
    v = df.values.T.tolist()
    for i in range(0,dfRows):
#        print ( '{:7.2f}{:>6.2f}{:>2.0f}{:>4.0f}{:>5.0f}{:6.2f}{:6.2f}{:6.2f}{:6.2f}{:>15}{:>60}'\
#        .format(v[0][i],v[1][i],v[2][i],v[3][i],v[4][i],v[5][i],v[6][i],v[7][i],v[8][i],\
#        v[9][i],v[10][i]) )
        dbOut.write(( \
        '{:7.2f}{:>6.2f}{:>2.0f}{:>4.0f}{:>5.0f}{:6.2f}{:6.2f}{:6.2f}{:6.2f}{:>15}{:>60}'\
        .format(v[0][i],v[1][i],v[2][i],v[3][i],v[4][i],v[5][i],v[6][i],v[7][i],v[8][i],\
        v[9][i],v[10][i]) ))
        dbOut.write("\n")
    print "Closing the file..."
    dbOut.close
"""
Use this switch to segregate specific values per field and plot them.
"""
if(search==1):
    # search for specific values
    fromCol     = 'month'       # What column/field are we interested in?
    values      = (12,1,2)  # Put what you're looking for here!
    uniqueVal   = len(values)   # This is for automation
    uniqueQuar  = 1
        
    querry = df.loc[df[fromCol] == values[0]]
    for i in range(1,uniqueVal):
        querry=querry.append(df.loc[df[fromCol] == values[i]])
    print"found",len(querry),"instances of",fromCol,"with values of:",values
    
    # plot the results!
    refList = querry.ref.unique()
    refSize = refList.size
    
    plt.figure(figsize=(12,12))
    map = Basemap(projection='mill', lat_0=0, lon_0=0)
    map.drawcoastlines()
    map.fillcontinents()
    map.drawmapboundary()
    
    x1 = querry.lon.values.T.tolist()
    y1 = querry.lat.values.T.tolist() 
    x, y = map(x1, y1)
    
    map.scatter(x, y, marker='o',color='k',label=querry.ref.unique())
    
    plt.title('Winter Distribution of Observed d18O')
    plt.show
"""
Use this function to section off the observed datapoints seasonally.
"""   
if(seasons==1):
    # Seasonal split!
    fromCol     = 'month'       # What column/field are we interested in?
    winter      = (12,1,2)  # Put what you're looking for here!
    spring      = (3,4,5)  # Put what you're looking for here!
    summer      = (6,7,8)  # Put what you're looking for here!
    fall        = (9,10,11)  # Put what you're looking for here!
    uniqueVal   = len(winter)   # This is for automation
    
    #   Section off winter months  
    querryW = df.loc[df[fromCol] == winter[0]]
    for i in range(1,uniqueVal):
        querryW=querryW.append(df.loc[df[fromCol] == winter[i]])
    refListW = querryW.ref.unique()
    refSize = refListW.size
    print"found",len(querryW),"instances of",fromCol,"with values of:",winter,"from",refSize,"unique references"
    # Section off spring months  
    querryS = df.loc[df[fromCol] == spring[0]]
    for i in range(1,uniqueVal):
        querryS=querryS.append(df.loc[df[fromCol] == spring[i]])
    refListS = querryS.ref.unique()
    refSize = refListS.size
    print"found",len(querryS),"instances of",fromCol,"with values of:",spring,"from",refSize,"unique references"
    # section off summer months
    querrySM = df.loc[df[fromCol] == summer[0]]
    for i in range(1,uniqueVal):
        querrySM=querrySM.append(df.loc[df[fromCol] == summer[i]])
    refListSM = querrySM.ref.unique()
    refSize = refListSM.size    
    print"found",len(querrySM),"instances of",fromCol,"with values of:",summer,"from",refSize,"unique references"
    # section off fall months
    querryF = df.loc[df[fromCol] == fall[0]]
    for i in range(1,uniqueVal):
        querryF=querryF.append(df.loc[df[fromCol] == fall[i]])
    refListF = querryF.ref.unique()
    refSize = refListF.size
    print"found",len(querryF),"instances of",fromCol,"with values of:",fall,"from",refSize,"unique references"   
    
    # plot the dictributions...
    plt.figure(figsize=(12,12))
    map = Basemap(projection='mill', lat_0=0, lon_0=0)
    map.drawcoastlines()
    map.fillcontinents()
    map.drawmapboundary()
    
    #establish lon and Lat per season...    
    x1 = querryS.lon.values.T.tolist()
    y1 = querryS.lat.values.T.tolist() 
    xS, yS = map(x1, y1)
    map.scatter(xS, yS, marker='o',color='g',label='Mar Apr May')
    plt.title('Spring Distribution of Observed d18O')
    plt.show
    
    x1 = querrySM.lon.values.T.tolist()
    y1 = querrySM.lat.values.T.tolist() 
    xSM, ySM = map(x1, y1)
    map.scatter(xSM, ySM, marker='o',color='r',label='Jun Jul Aug')
#    plt.title('Summer Distribution of Observed d18O')
#    plt.show
    
    x1 = querryF.lon.values.T.tolist()
    y1 = querryF.lat.values.T.tolist() 
    xF, yF = map(x1, y1)
    map.scatter(xF, yF, marker='o',color='orange',label='Sept Oct Nov')
#    plt.title('Fall Distribution of Observed d18O')
#    plt.show
    
    x1 = querryW.lon.values.T.tolist()
    y1 = querryW.lat.values.T.tolist() 
    xW, yW = map(x1, y1)
    map.scatter(xW, yW, marker='o',color='k',label='Dec Jan Feb')
#    plt.title('Winter Distribution of Observed d18O')
#    plt.show    
    
    plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.05),ncol=4)
    plt.title('Seaonal Distribution of Observed d18O')
    plt.show
    plt.savefig('seasonalDistribution', dpi = 300)
    
# Scatter plot for all points   WIP!!! 
if (plot==1):
    refList = df.ref.unique()
    refSize = refList.size
    plt.figure(figsize=(12,12))
    map = Basemap(projection='mill', lat_0=0, lon_0=0)
    map.drawcoastlines()
    map.fillcontinents()
    map.drawmapboundary()
    # Make a plot according to reference
    for i in range(0,refSize):
        flavor="#%06x" % random.randint(0,0xFFFFFF)
        #breaks dataframe into subdataframes of unique references    
        uniqueRef=refList[i]
        dfSub=df.loc[df['ref']==uniqueRef]
        x1 = dfSub.lon.values.T.tolist()
        y1 = dfSub.lat.values.T.tolist()
#        for j in range (0,len(dfSub.index)):
        x, y = map(x1, y1)
        map.scatter(x, y, marker='o',color=flavor,label=dfSub.ref.unique())
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.0)
    plt.title('Observed d18O data')
    plt.show
    plt.savefig('observedD18O', dpi = 300)