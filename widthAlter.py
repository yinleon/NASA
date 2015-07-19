import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import random
from mpl_toolkits.basemap import Basemap

"""
This is a program that takes in a string list of properly formatted CSV files
(toAdd), to add to the Global Seawater d18O database (infile).

Added...
Ref                 File Name               # of Recs       Compiler    Dates
Strauss 2015        'Strauss2015.csv'       271 records     L           7.7.15
Cases 2012          'CASES-Database.csv'    485 records     C   
COX 2010,2001       'Cox-Database.csv'      88  records     C
Conroy 2012         'Conroy2012.csv'        191 rows        S
"""


filein      ='dbO18.txt'
toAdd       =('Strauss2015.csv','CASES-Database.csv','Cox-Database.csv',
              'Conroy2012.csv')
fileout     ='dbO18_2015.txt'

# Logic switches -- use these guys to turn functions on-and-off.
convert         =   2       #set to 1: convert the database (DB) into a CSV, 
                            #set to 2: convert CSV input to DB
append2master   =   1       #set to 1: to add new CSV input to the DB
printOut        =   1       #set to 1: prints out DB
plot            =   1       #set to 1: plots a scatter for the DB
search          =  -1       #set to 1: search spefici querry

# Global Seawaster d18) database format
col=['lon','lat','month','year','depth','temp', 'sal', 'd18o','dD','notes','ref'] 
fWid = ( 7,    6,      2,     4,      5,     6,     6,      6,   6,    15,   60)

# Read in Database into a Pandas Dataf rame data structure
df=pd.read_fwf(filein,widths=fWid,header=None,index_col=False,names=col,skip_blank_lines=True)

dfRows    = len(df.index-1)
files2Add = len(toAdd)

# Convert DF to to CVS
if (convert == 1):
    toCSV=pd.DataFrame.to_csv(df,sep=',',line_terminator='\n')

# Automated reading of each CSV.
# To troubleshoot adjust the range and double-check 'csv' in the shell.
if (convert == 2):
    for i in range(0,files2Add):
        csv=pd.read_csv(toAdd[i],sep=',',skiprows=1,names=col,na_filter=True)
        #              File Name,Delimiter, skip header
        if (append2master == 1):
            df = df.append(csv)
            dfRows= len(df.index-1)

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
        dbOut.write(( '{:7.2f}{:>6.2f}{:>2.0f}{:>4.0f}{:>5.0f}{:6.2f}{:6.2f}{:6.2f}{:6.2f}{:>15}{:>60}'\
        .format(v[0][i],v[1][i],v[2][i],v[3][i],v[4][i],v[5][i],v[6][i],v[7][i],v[8][i],\
        v[9][i],v[10][i]) ))
        dbOut.write("\n")
    print "Closing the file..."
    dbOut.close

if(search==1):
    # search for specific values
    fromCol     = 'month'       # What column/field are we interested in?
    values      = (9,10,11,12)  # Put what you're looking for here!
    uniqueVal   = len(values)   # This is for automation
    uniqueQuar  = 1
    # 
    querry = df.loc[df[fromCol]== values[0]]
    for i in range(1,uniqueVal):
        querry=querry.append(df.loc[df[fromCol]== values[i]])
    print"found",len(querry),"instances of",fromCol,"with values of:",values
    
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
    plt.savefig('sampleScatter', dpi = 300)