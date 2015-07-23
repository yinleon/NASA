import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import random
from mpl_toolkits.basemap import Basemap
#from pylab import *
import colorsys
from random import shuffle

"""
This is a program that takes in a list of properly formatted CSV files(toAdd), 
to append to the Global Seawater d18O database (infile) using Pandas.

This thanks to Pandas, there are dynamic database (DB) querrying for any of 
the 11 fields-- including seasonality, reference, and location.

The database is converted to Pandas dataframe(DF) data structure, 
The DF can be accessed by typing 'df' into the command line after running.

Pandas has thorough documentation, and great examples on their website below:
--> http://pandas.pydata.org/pandas-docs/stable/tutorials.html

The querry results can be plotted globally using Basemaps. 
There is an example available by turning the 'seasons' switch to 1.

the integer seasons is a logic switch(much like a boolean), used to access
specific functions in this program for writing out files, searching the DF

At this time, this program is intended to update these files...

Ref                 File Name               # of Recs       Compiler    Dates
Strauss 2015        'Strauss2015.csv'       271 records     L           7.17.15
Cases 2012          'Cases2012.csv'         485 records     C   
COX 2010,2001       'COX2010_2001.csv'      88  records     C
Conroy 2012         'Conroy2012.csv'        191 records     S
Tiwari 2013         'Tiwari2013.csv'        76  records     S           7.20.15
Brown 2014          'BROWN2014.csv'                                     7.25.15
Geotraces 2014...
Rijkenberg 2011     'GeoRijkenberg2011.csv' 432 records     C
Jenkins 2010        'GeoJenkins2010.csv'    192 records     C  

Missing...
Geotraces 2010-onwards          lacks longitude
Brown 2014
srep                            lacks coordinates
Munksgaard 2012                 from paper

Written by Leon Yin for NASA GISS NYCRI 2015
Please consult NASA GISS or the author for use.

v - 1.0
"""

# Files to be read/written in the parent directory
filein      = 'dbO18.txt'       # this file is downloaded from the d18o website.
fileout     = 'dbO18_2015.txt'  # output of the entire DB as a fixed-width txt
csvout      = 'db18O_2015.csv'  # replace 2015 with specifics!
toAdd       = ('Strauss2015.csv','Cases2012.csv','COX2010_2001.csv',
              'Conroy2012.csv','Tiwari2013.csv','BROWN2014.csv',
              'GeoRijkenberg2011.csv','GeoJenkins2010.csv')

# Logic switches -- use these guys to turn functions on-and-off.
convert2DF      =   1               #set to 1: convert CSV input to DB
convert2CSV     =  -1               #set to 1: convert the DB into a CSV, 
append2master   =   1               #set to 1: to add new CSV input to the DB
writeOut        =  -1               #set to 1: writes out DB to fileout
plotALL         =  -1               #set to 1: plots a scatter for the DB
search          =  -1               #set to 1: search specfic querry
seasons         =  -1               #set to 1: seasonal distribution of observtions
uniqueArc       =   1

# Global Seawater d18O database format.
col=\
  ['lon','lat','month','year','depth','temp', 'sal', 'd18o','dD','notes','ref'] 
fWid=( 7,    6,      2,     4,      5,     6,     6,      6,   6,    15,   60)

# Read in Database into a Pandas Dataframe data structure
df=pd.read_fwf(filein,widths=fWid,header=None,names=col,skip_blank_lines=True)

# some variables...
dfRows    = len(df.index-1)
ogLen     = len(df.index-1)

files2Add = len(toAdd)
"""****************************************************************************
Automated reading of each CSV, and appending to existing DB.
To troubleshoot adjust the range and double-check 'csv' in the shell.
****************************************************************************"""
if (convert2DF == 1):
    for i in range(0,files2Add):
        print "opening",toAdd[i],"..."
        csv=pd.read_csv(toAdd[i],sep=',',skiprows=1,names=col,na_filter=True)
        #              File Name,Delimiter, skip header
        dfRows= len(df.index-1)
        if (append2master == 1):
            df = df.append(csv)     # add the contents of new csv into DB.
            dfRows= len(df.index-1) # update the number of records in the DB.
    print dfRows-ogLen,"records added to the database."
"""****************************************************************************
Use this function to convert DF to to CVS, and write it out as a csv.
****************************************************************************"""
if (convert2CSV == 1):
    toCSV=pd.DataFrame.to_csv(df,sep=',',line_terminator='\n')
    csvout = open(csvout, 'w')      # opens the file to write into
    csvout.write(toCSV)             # writes df to csv... 
    print "Database transferred to csv..."
"""****************************************************************************
Use this function to write the database into a fixed-width text file.
****************************************************************************"""
if (writeOut == 1):
    df=df.fillna('')
    print("******************************************************************")        
    dbOut = open(fileout, 'w')
    v = df.values.T.tolist()        # The DF must be converted to a list first.
    for i in range(0,dfRows):       # Right each row into the file.
#        print ( '{:7.2f}{:>6.2f}{:>2.0f}{:>4.0f}{:>5.0f}{:6.2f}{:6.2f}{:6.2f}{:6.2f}{:>15}{:>60}'\
#        .format(v[0][i],v[1][i],v[2][i],v[3][i],v[4][i],v[5][i],v[6][i],v[7][i],v[8][i],\
#        v[9][i],v[10][i]) )
        dbOut.write(( \
        '{:7.2f}{:>6.2f}{:>2.0f}{:>4.0f}{:>5.0f}{:6.2f}{:6.2f}{:6.2f}{:6.2f}{:>15}{:>60}'\
        .format(v[0][i],v[1][i],v[2][i],v[3][i],v[4][i],v[5][i],v[6][i],v[7][i],v[8][i],\
        v[9][i],v[10][i]) ))
        dbOut.write("\n")
    print "Database transferred to txt..."
    dbOut.close
"""****************************************************************************
Use this switch to segregate specific values per field and plot them.
I've used Winter distribution as an example.
****************************************************************************"""
if (search == 1):
    # search for specific values
    fromCol     = 'month'           # What column/field are you interested in?
    values      = (12,1,2)          # Put what you're looking for here!
    uniqueVal   = len(values)       # How many automations are required.
        
    querry = df.loc[df[fromCol] == values[0]]
    for i in range(1,uniqueVal):
        querry=querry.append(df.loc[df[fromCol] == values[i]])
    print"found",len(querry),"instances of",fromCol,"with values of:",values
    
    # Book keeping...
    refList = querry.ref.unique()   # List of search-relevant references.
    refSize = refList.size
    
    # plot the results!
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
"""****************************************************************************
Use this function to section off the observed datapoints seasonally.
****************************************************************************"""
if (seasons == 1):
    # Seasonal split!
    fromCol     = 'month'           # What column/field are we interested in?
    winter      = (12,1,2)          # Winter Months
    spring      = (3,4,5)           # Spring Months
    summer      = (6,7,8)           # Summer Months
    fall        = (9,10,11)         # Fall Months (my favorite...)
    uniqueVal   = len(winter)       # This is for automation'ssake
    
    #   Section off winter months  
    querryW = df.loc[df[fromCol] == winter[0]]
    for i in range(1,uniqueVal):
        querryW=querryW.append(df.loc[df[fromCol] == winter[i]])
    refListW = querryW.ref.unique()
    refSize = refListW.size
    print"found",len(querryW),"instances of",fromCol,"with values of:",winter,\
    "from",refSize,"unique references"
    
    # Section off spring months  
    querryS = df.loc[df[fromCol] == spring[0]]
    for i in range(1,uniqueVal):
        querryS=querryS.append(df.loc[df[fromCol] == spring[i]])
    refListS = querryS.ref.unique()
    refSize = refListS.size
    print"found",len(querryS),"instances of",fromCol,"with values of:",spring,\
    "from",refSize,"unique references"
    
    # Section off summer months
    querrySM = df.loc[df[fromCol] == summer[0]]
    for i in range(1,uniqueVal):
        querrySM=querrySM.append(df.loc[df[fromCol] == summer[i]])
    refListSM = querrySM.ref.unique()
    refSize = refListSM.size    
    print"found",len(querrySM),"instances of",fromCol,"with values of:",summer,\
    "from",refSize,"unique references"
    
    # section off fall months
    querryF = df.loc[df[fromCol] == fall[0]]
    for i in range(1,uniqueVal):
        querryF=querryF.append(df.loc[df[fromCol] == fall[i]])
    refListF = querryF.ref.unique()
    refSize = refListF.size
    print"found",len(querryF),"instances of",fromCol,"with values of:",fall,\
    "from",refSize,"unique references"   
    
    # plot the dictributions...
    plt.figure(figsize=(12,12))
    map = Basemap(projection='npstere',boundinglat=60,lon_0=330,resolution='l')
    parallels = np.arange(-90.,90.,4.)
    map.drawparallels(parallels,labels=[False,False,False,True])
    meridians = np.arange(-180.,181.,5.)
    map.drawmeridians(meridians,labels=[True,False,False,False])
    map.drawcoastlines()
    map.fillcontinents()
    map.drawmapboundary()
    
    #establish lon and Lat per season...    
    x1 = querryS.lon.values.T.tolist()
    y1 = querryS.lat.values.T.tolist() 
    xS, yS = map(x1, y1)
    map.scatter(xS, yS, marker='o',color='g',label='Mar Apr May',s=4)
    plt.title('Spring Distribution of Observed d18O')
    plt.show
    
    x1 = querrySM.lon.values.T.tolist()
    y1 = querrySM.lat.values.T.tolist() 
    xSM, ySM = map(x1, y1)
    map.scatter(xSM, ySM, marker='o',color='r',label='Jun Jul Aug',s=4)
#    plt.title('Summer Distribution of Observed d18O')
#    plt.show
    
    x1 = querryF.lon.values.T.tolist()
    y1 = querryF.lat.values.T.tolist() 
    xF, yF = map(x1, y1)
    map.scatter(xF, yF, marker='o',color='orange',label='Sept Oct Nov',s=4)
#    plt.title('Fall Distribution of Observed d18O')
#    plt.show
    
    x1 = querryW.lon.values.T.tolist()
    y1 = querryW.lat.values.T.tolist() 
    xW, yW = map(x1, y1)
    map.scatter(xW, yW, marker='o',color='k',label='Dec Jan Feb',s=4)
#    plt.title('Winter Distribution of Observed d18O')
#    plt.show    
    
    plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.05),ncol=4)
    plt.title('Seaonal Distribution of Observed d18O')
    plt.show
    plt.savefig('seasonalDistribution', dpi = 300)
"""****************************************************************************
Use this to generate a scatter plot for all points by unique reference.
****************************************************************************"""
if (plotALL == 1):
    refList = df.ref.unique()
    refSize = refList.size
    plt.figure(figsize=(12,12))
    map = Basemap(projection='robin', lat_0=0, lon_0=0)
    map.drawcoastlines()
    map.fillcontinents()
    map.drawmapboundary()
    # Make a plot according to reference
    for i in range(0,refSize):
        flavor="#%06x" % random.randint(0,0xFFFFFF)
        #breaks dataframe into subdataframes of unique references    
        uniqueRef=refList[i]
        dfSub=df.loc[df['ref']==uniqueRef]
        x1 = dfSub.lon.values.T.tolist()    # convert longitude values to list
        y1 = dfSub.lat.values.T.tolist()    # convert latitude valus to list
        x, y = map(x1, y1)
        map.scatter(x, y, marker='o',color=flavor,label=dfSub.ref.unique(),s=3)
#    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.0)
    plt.title('Observed d18O data')
    plt.show
    plt.savefig('observedD18O', dpi = 300)
    plt.show
    
if (uniqueArc == 1):
    def _get_colors(num_colors):
        colors=[]
        for i in np.arange(0., 360., 360. / num_colors):
            hue = i/360.
            lightness = (50 + np.random.rand() * 10)/100.
            saturation = (90 + np.random.rand() * 10)/100.
            colors.append(colorsys.hls_to_rgb(hue, lightness, saturation))
        return colors
    # Define the Arctic...
    newArc = df.loc[df.lat>=66]
    
    refList = newArc.ref.unique()
#    shuffle(refList)    
    
    refSize = refList.size

#    Various Schemes of getting colors
#    z=0
#    N = refSize
#    HSV_tuples = [(z*1.0/N, 0.8, 0.8) for z in range(N)]
#    colorSwatch = map(lambda z: colorsys.hsv_to_rgb(*z), HSV_tuples)
# 
    colorSwatch=_get_colors(refSize)
    
    mark = ['<','>','v','^']
    
    plt.figure(figsize=(12,12))
    map = Basemap(projection='npstere',boundinglat=60,lon_0=330,resolution='l')
    parallels = np.arange(-90.,90.,4.)
    map.drawparallels(parallels,labels=[False,False,False,True])
    meridians = np.arange(-180.,181.,5.)
    map.drawmeridians(meridians,labels=[True,False,False,False])
    map.drawcoastlines()
    map.fillcontinents()
    map.drawmapboundary()
    # Make a plot according to reference
    for i in range(0,refSize):
#        flavor="#%06x" % random.randint(0,0xFFFFFF)
        #breaks dataframe into subdataframes of unique references    
        uniqueRef=refList[i]
        dfSub=newArc.loc[newArc['ref']==uniqueRef]
        x1 = dfSub.lon.values.T.tolist()    # convert longitude values to list
        y1 = dfSub.lat.values.T.tolist()    # convert latitude valus to list
        x, y = map(x1, y1)
        map.scatter(x, y, marker=random.choice(mark),color=colorSwatch[i],label=str(dfSub.ref.unique()),s=7)
    art = []
    lgd = plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.0,prop={'size':10})
    art.append(lgd)    
    plt.title('Spatial distribution of Arctic d18O Observations')
    plt.savefig('arcDist.png', bbox_inches="tight" ,additional_artists=art,format='png', dpi=1000)
#    plt.savefig('arcDist.eps', format='eps', dpi=1000)
    plt.show
    

#####
