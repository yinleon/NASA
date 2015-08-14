import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import os
from mpl_toolkits.basemap import Basemap
import colorsys

"""
This is a program that takes in a list of properly formatted CSV files(toAdd), 
to append to the Global Seawater d18O database (infile) using Pandas.

This thanks to Pandas, there are dynamic database (DB) querrying for any of 
the 11 fields-- including seasonality, reference, and location.

The database is converted to Pandas dataframe(DF) data structure, 
The DF can be accessed by typing 'df' into the command line after running.

Pandas has thorough documentation on their website below:
--> http://pandas.pydata.org/pandas-docs/stable/tutorials.html

The querry results can be plotted globally using Basemaps. 
There is an example available by turning the 'seasons' switch to 1.

the integer seasons is a logic switch(much like a boolean), used to access
specific functions in this program for writing out files, searching the DF

At this time, this program is intended to update these files...

Ref                 File Name               # of Recs       Compiler    Dates
Strauss 2015        'Strauss2015.csv'       271  records    L           7.17.15
Cases 2012          'Cases2012.csv'         485  records    C   
COX 2010,2001       'COX2010_2001.csv'      88   records    C
Conroy 2012         'Conroy2012.csv'        191  records    S
Tiwari 2013         'Tiwari2013.csv'        76   records    S           7.20.15
Brown 2014          'BROWN2014.csv'         1207 records    C           7.25.15
Geotraces 2014...
Rijkenberg 2011     'GeoRijkenberg2011.csv' 432 records     C
Jenkins 2010        'GeoJenkins2010.csv'    192 records     C  

total:                                      2942

Written by Leon Yin for NASA GISS NYCRI 2015
Please consult NASA GISS or the author(leon.yin@nyu.edu)for use.

v - 1.0
"""
# Files to be read/written in the parent directory
filein      = 'dbO18.txt'       # this file is downloaded from the d18o website.
fileout     = 'dbO18_2015.txt'  # output of the entire DB as a fixed-width txt
csvout      = 'eurasianBasin.csv'  # replace 2015 with specifics!
toAdd       = ('Strauss2015.csv','Cases2012.csv','COX2010_2001.csv',
              'Conroy2012.csv','Tiwari2013.csv','BROWN2014.csv',
              'GeoRijkenberg2011.csv','GeoJenkins2010.csv')

# Logic switches -- use these guys to turn functions on-and-off.
convert2DF      =   1               #set to 1: convert CSV input to DB
convert2CSV     =  -1               #set to 1: convert the DB into a CSV, 
append2master   =   1               #set to 1: to add new CSV input to the DB
writeOut        =   1               #set to 1: writes out DB to fileout
plotALL         =   1               #set to 1: plots a scatter for the DB
search          =   1               #set to 1: search specfic querry
seasons         =   1               #set to 1: seasonal distribution of observtions
uniqueArc       =   1
outlierELIM     =   1
annual          =   1
melt            =   1
writeOutIce     =  -1
save            =  -1

"""****************************************************************************
Some necessary functions
****************************************************************************"""
def _get_colors(num_colors):
    colors=[]
    for i in np.arange(0., 360., 360. / num_colors):
        hue = i/360.
        lightness   = (50 + np.random.rand() * 10) / 100.
        saturation  = (90 + np.random.rand() * 10) / 100.
        colors.append(colorsys.hls_to_rgb(hue, lightness, saturation))
    return colors

"""~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The following section should be used to update and export the database
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"""
# Global Seawater d18O database format.
col=\
  ['lon','lat','month','year','depth','temp', 'sal', 'd18o','dD','notes','ref'] 
fWid=( 7,    6,      2,     4,      5,     6,     6,      6,   6,    15,   60)

# Read in Database into a Pandas Dataframe data structure
df=pd.read_fwf(filein,widths=fWid,header=None,names=col,skip_blank_lines=True)

# some variables...
dfRows    = len(df.index-1) # number of records in current DF
ogLen     = len(df.index-1) # number of records in DF from filein
files2Add = len(toAdd)      # number of sources to add


"""****************************************************************************
Automates reading each CSV, and appending to existing DB.
To troubleshoot adjust the range and double-check 'csv' in the shell.
****************************************************************************"""
if (convert2DF == 1):
    print("******************************************************************")        
    print"Integrating CSV into Global d18o Seawater Database..."    
    for i in range(0,files2Add):
        print "opening",toAdd[i],""
        csv2DF=pd.read_csv(toAdd[i],sep=',',skiprows=1,names=col,na_filter=True)
        #              File Name,Delimiter, skip header
        dfRows= len(df.index-1)
        if (append2master == 1):
            df = df.append(csv2DF)  # add contents of new csv into DB.
            dfRows= len(df.index-1) # update the number of records in DB.
    print dfRows-ogLen,"records added to the database."
    print("******************************************************************")        

"""****************************************************************************
Use this function to convert DF to to CVS, and write it out as a csv.
****************************************************************************"""
if (convert2CSV == 1):
    df2csv = 'dbo18_2015.csv'
    toCSV=pd.DataFrame.to_csv(df2csv,sep=',',line_terminator='\n')
    csvout = open(csvout, 'w')      # opens the file to write into
    csvout.write(toCSV)             # writes df to csv... 
    print "Database transferred to csv..."
    csvout.close()
"""****************************************************************************
Seek vengence on outliers from variable high latitutes.
ID's determined using Jack knife statstical analysis in osSlopes.py...
Known troublsome data...

id     lon  lat  month  year  depth   temp    sal  d18o    dD                ref 
18067    3   80      6  1984      0  -1.20  20.00  0.14 -99.9  Ostlund and Grall (1993)  

replace hard coded idOdeath with killList from osSlopes.py for easy(but mindless)
removal of outliers.
****************************************************************************"""
if (outlierELIM == 1):
    print("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX")
    print "Destroying outliers..."
    print("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX")

    idOdeath    = [18067,519,825,613,18383]
    lenOdeath   = len(idOdeath)
    for toKill in range(lenOdeath):     # filters DF by natural index.
        df = df[df.index != idOdeath[toKill]] 
    print lenOdeath,"outliers ELIMINATED!"
    dfRows    = len(df.index-1)
    print("\m/ -___- \m/")
"""****************************************************************************
Use this function to write the database into a fixed-width text file.
****************************************************************************"""
if (writeOut == 1):
    df=df.fillna('')                #replace NaN values with spaces
    temp = 'temp.txt'
    print("******************************************************************")        
    dbOut = open(temp, 'w')
    v = df.values.T.tolist()        # The DF must be converted to a list first.
    for i in range(0,dfRows):       # Right each row into the file.
        dbOut.write(( \
        '{:7.2f}{:>6.2f}{:>2.0f}{:>4.0f}{:>5.0f}{:6.2f}{:6.2f}{:6.2f}{:6.1f}{:>15}{:>60}'\
        .format(v[0][i],v[1][i],v[2][i],v[3][i],v[4][i],v[5][i],v[6][i],v[7][i],v[8][i],\
        v[9][i],v[10][i]) ))
        dbOut.write("\n")
    dbOut.close
    # swap -99.90 with -99.9
    with open(fileout, "wt") as fout:
        with open(temp, "rt") as fin:
            for line in fin:
                fout.write(line.replace('-99.90', ' -99.9'))
    print "Database transferred to txt..."
    fout.close
    fin.close
    os.remove(temp)                 # Remove the temporary file.
"""~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The following section should be used to search and the database.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"""

"""****************************************************************************
Use this switch to segregate specific values per field and plot them.
I've used Winter distribution as an example.
****************************************************************************"""
if (search == 1):
    # search for specific values
    fromCol     = 'month'           # What column are you interested in?(ln:72)
    values      = (12,1,2)          # Put what you're looking for here!
    uniqueVal   = len(values)       # How many automations are required?
    
    # Collect values from specified col into new DF.
    querry = df.loc[df[fromCol] == values[0]]
    for i in range(1,uniqueVal):
        querry=querry.append(df.loc[df[fromCol] == values[i]])
    print"found",len(querry),"instances of",fromCol,"with values of:",values
    
    # Book keeping...
    refList = querry.ref.unique()   # List of search-relevant references.
    refSize = refList.size
    
    # plot the results! Below is optional
    plt.figure(figsize=(12,12))
    map = Basemap(projection='mill', lat_0=0, lon_0=0)
    map.drawcoastlines()
    map.fillcontinents()
    map.drawmapboundary()
    
    x1 = querry.lon.values.T.tolist()
    y1 = querry.lat.values.T.tolist() 
    x, y = map(x1, y1)
    
    map.scatter(x, y, marker='o',color='k',label=querry.ref.unique())
    
    plt.title(r'Winter Distribution of Observed $\delta^{18}$O')
    plt.show
"""****************************************************************************
Use this function to section off the observed datapoints seasonally.
This is an extension of the above function.
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
    plt.figure(figsize=(7,7))
    # choose one of the projections--
    # global projections
#    map = Basemap(projection='mill',llcrnrlat=-90,urcrnrlat=90,
#    llcrnrlon=-180,urcrnrlon=180,resolution='c')
    # Arctic projection    
    map = Basemap(projection='npstere',boundinglat=68,lon_0=330,resolution='l')
    parallels = np.arange(-90.,90.,4.)
    map.drawparallels(parallels,labels=[False,False,False,True])
    meridians = np.arange(-180.,181.,10.)
    map.drawmeridians(meridians,labels=[True,False,False,False])
    
    map.drawcoastlines()
    map.fillcontinents()
    map.drawmapboundary()
    
    #establish lon and Lat per season...    
    x1 = querryS.lon.values.T.tolist()
    y1 = querryS.lat.values.T.tolist() 
    xS, yS = map(x1, y1)
    map.scatter(xS, yS, marker='o',color='g',label='Mar Apr May',s=16)

    x1 = querrySM.lon.values.T.tolist()
    y1 = querrySM.lat.values.T.tolist() 
    xSM, ySM = map(x1, y1)
    map.scatter(xSM, ySM, marker='o',color='r',label='Jun Jul Aug',s=16)

    x1 = querryF.lon.values.T.tolist()
    y1 = querryF.lat.values.T.tolist() 
    xF, yF = map(x1, y1)
    map.scatter(xF, yF, marker='o',color='orange',label='Sept Oct Nov',s=16)

    x1 = querryW.lon.values.T.tolist()
    y1 = querryW.lat.values.T.tolist() 
    xW, yW = map(x1, y1)
    map.scatter(xW, yW, marker='o',color='k',label='Dec Jan Feb',s=16)
    # stuff so the for legend is saved w/ the fig!
    art = []
    lgd = plt.legend(loc='upper center',scatterpoints=1,ncol=4,
         bbox_to_anchor=(0.5, -0.05),prop={'size':12})
    art.append(lgd)      
    
    plt.title(r'Tempoal Distribution of Observed $\delta^{18}$O',size=16)
    plt.show
    if (save ==1):
        plt.savefig('seasonality.png', bbox_inches="tight",
                    additional_artists=art,format='png',dpi = 300)
    # Bookkeeping...    
    seasonDF = (querryF,querryW,querryS,querrySM)
"""****************************************************************************
Inter-annual distribution...
Something to be explored -LY
****************************************************************************"""
if (annual == 1):
    plt.figure(figsize=(9,7))  
    map = Basemap(projection='npstere',boundinglat=68,lon_0=330,resolution='l')
    parallels = np.arange(-90.,90.,4.)
    meridians = np.arange(-180.,181.,10.)
    map.drawparallels(parallels,labels=[False,False,False,True])
    map.drawmeridians(meridians,labels=[True,False,False,False])
#    map = Basemap(projection='mill', lat_0=0, lon_0=0)
    map.drawcoastlines()
    map.fillcontinents()
    map.drawmapboundary()
    
    df=df[df.year!=-999]
    annie = df.year.unique()
    annie.sort()
    annieHall = len(annie)
    colorSwatch=_get_colors(annieHall)
    for i in range(annieHall):
#        print annie[i]
        annualDF = df.loc[df.year==annie[i]]
        annualDF = annualDF[annualDF.depth<=10]
        

        x1 = annualDF.lon.values.T.tolist()
        y1 = annualDF.lat.values.T.tolist()
        x, y = map(x1, y1)
        
        map.scatter(x, y, marker=mark[i%4],color=colorSwatch[i],
                    label=annie[i],s=8)
    art = []
    lgd = plt.legend(loc='upper center',scatterpoints=1,ncol=6,
         bbox_to_anchor=(0.5, -0.05),prop={'size':10})
    plt.title('Interannual Distribution of $\delta^{18}$O',fontsize=14)
    if(save == 1):
        plt.savefig('interannual', dpi = 300)
    plt.show()    
"""****************************************************************************
Use this to generate a scatter plot for all points by unique reference.
****************************************************************************"""
if (plotALL == 1):
    refList = df.ref.unique()       # creates list of unique references from DB
    refSize = refList.size          # how large is the list?
    
    colorSwatch=_get_colors(refSize)# get enough colors for each point. 
    mark = ['o','v','s','d']        # pool of shapes for scatter
    
    plt.figure(figsize=(10,10))     # initiate plot w/ dimensions.
    map = Basemap(projection='robin', lat_0=0, lon_0=0)
    map.drawcoastlines()
    map.fillcontinents()
    map.drawmapboundary()
    
    for i in range(0,refSize):
        #breaks dataframe into subdataframes of unique references    
        uniqueRef=refList[i]
        dfSub=df.loc[df['ref']==uniqueRef]
        x1 = dfSub.lon.values.T.tolist()    # convert longitude values to list
        y1 = dfSub.lat.values.T.tolist()    # convert latitude valus to list
        x, y = map(x1, y1)                  # convert the lon and lat into map
        map.scatter(x, y, marker=mark[i%4],color=colorSwatch[i],
                    label=str(refList[i]),s=3)
    art = []
    lgd = plt.legend(scatterpoints=1,bbox_to_anchor=(1.05, 1), 
                     loc=2, borderaxespad=0.0,prop={'size':8})
#    art.append(lgd)     
   
    plt.title(r'Global Observed $\delta18O$ data')
    plt.show
    if (save == 1):
        plt.savefig('observedD18O.png', bbox_inches="tight",
                    additional_artists=art,format='png',dpi = 300)

"""****************************************************************************
Use this to generate a scatter plot for all arctic points by unique reference.
****************************************************************************""" 
if (uniqueArc == 1):
    # Define the Arctic...
    newArc = df.loc[df.lat>=65.8]
    refList = newArc.ref.unique()
    refList.sort()
    refSize = refList.size

    # Various Schemes of getting colors
    colorSwatch=_get_colors(refSize)
    
    mark = ['o','v','s','d']
    
    plt.figure(figsize=(10,10))
    map = Basemap(projection='npstere',boundinglat=66,lon_0=330,resolution='l')
    parallels = np.arange(-90.,90.,2.)
    map.drawparallels(parallels,labels=[False,False,False,True])
    meridians = np.arange(-180.,181.,5.)
    map.drawmeridians(meridians,labels=[True,False,False,False])
    map.drawcoastlines()
    map.fillcontinents()
    map.drawmapboundary()
    
    # Make a plot according to reference
    for i in range(0,refSize):
        #breaks dataframe into subdataframes of unique references    
        uniqueRef=refList[i]
        dfSub=newArc.loc[newArc['ref']==uniqueRef]
        x1 = dfSub.lon.values.T.tolist()    # convert longitude values to list
        y1 = dfSub.lat.values.T.tolist()    # convert latitude valus to list
        x, y = map(x1, y1)
        map.scatter(x, y, marker=mark[i%4],color=colorSwatch[i],label=str(uniqueRef),s=7)
    art = []
    lgd = plt.legend(scatterpoints=1,bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.0,prop={'size':9})
    art.append(lgd)    
    plt.title(r'Spatial Distribution of Observed Arctic $\delta^{18}$O')
    if (save == 1):
        plt.savefig('arcDist.png', bbox_inches="tight" ,additional_artists=art,format='png', dpi=200)
#       plt.savefig('arcDist.eps', format='eps', dpi=1000) #for vectors
    plt.show
"""****************************************************************************
Separates DB into two DF's one for melt and one for frozen
****************************************************************************"""
if (melt == 1):
    print("******************************************************************")        
    fromCol     = 'month'           # What column/field are we interested in?
    ice     = (1,2,3,4,10,11,12)    # Winter Months
    water   = (5,6,7,8,9)           # Spring Months may-sept
    iLen    = len(ice)
    wLen    = len(water)            # This is for automation'ssake
    
    #   Section off frozen months  
    iceDF = df.loc[df[fromCol] == ice[0]]
    for i in range(1,iLen):
        iceDF=iceDF.append(df.loc[df[fromCol] == ice[i]])
    refListIce = iceDF.ref.unique()
    refSize = refListIce.size
    print"found",len(iceDF),"instances of",fromCol,"with values of:",ice,\
    "from",refSize,"unique references"
    
    # Section off melt months  
    waterDF = df.loc[df[fromCol] == water[0]]
    for i in range(1,wLen):
        waterDF=waterDF.append(df.loc[df[fromCol] == water[i]])
    refListWater = waterDF.ref.unique()
    refSize = refListWater.size
    print"found",len(waterDF),"instances of",fromCol,"with values of:",water,\
    "from",refSize,"unique references"
   
   # write the DB for ice and melt into .dat files for global gridded dataset.
    if (writeOutIce == 1):
        # Ice season
        iceDf=iceDF.fillna('')
        temp = 'temp.txt'
        fileout = 'ice.dat'
        print("******************************************************************")        
        dbOut = open(temp, 'w')
        v = iceDf.values.T.tolist()        # The DF must be converted to a list first.
        for i in range(0,len(iceDF)):       # Right each row into the file.
            dbOut.write(( \
            '{:7.2f}{:>6.2f}{:>2.0f}{:>4.0f}{:>5.0f}{:6.2f}{:6.2f}{:6.2f}{:6.1f}{:>15}{:>60}'\
            .format(v[0][i],v[1][i],v[2][i],v[3][i],v[4][i],v[5][i],v[6][i],v[7][i],v[8][i],\
            v[9][i],v[10][i]) ))
            dbOut.write("\n")
        dbOut.close
        # swap -99.90 with -99.9
        with open(fileout, "wt") as fout:
            with open(temp, "rt") as fin:
                for line in fin:
                    fout.write(line.replace('-99.90', ' -99.9'))
        print "Database transferred to txt..."
        fout.close
        fin.close
        os.remove(temp) 
        
        # Time for meltseason
        waterDf=waterDF.fillna('')
    
        fileout = 'melt.dat'
        dbOut = open(temp, 'w')
        v = waterDf.values.T.tolist()        # The DF must be converted to a list first.
        for i in range(0,len(waterDF)):       # Right each row into the file.
            dbOut.write(( \
            '{:7.2f}{:>6.2f}{:>2.0f}{:>4.0f}{:>5.0f}{:6.2f}{:6.2f}{:6.2f}{:6.1f}{:>15}{:>60}'\
            .format(v[0][i],v[1][i],v[2][i],v[3][i],v[4][i],v[5][i],v[6][i],v[7][i],v[8][i],\
            v[9][i],v[10][i]) ))
            dbOut.write("\n")
        dbOut.close
        # swap -99.90 with -99.9
        with open(fileout, "wt") as fout:
            with open(temp, "rt") as fin:
                for line in fin:
                    fout.write(line.replace('-99.90', ' -99.9'))
        print "Database transferred to txt..."
        fout.close
        fin.close
        os.remove(temp) 

    #surface only
    surfDep     =   10 
    iceDF       =  iceDF.loc[iceDF.depth<=surfDep]
    waterDF     =  waterDF.loc[waterDF.depth<=surfDep]
    
    # Plot 1
    plt.figure(figsize=(7,7))  
    map = Basemap(projection='npstere',boundinglat=68,lon_0=330,resolution='l')
    parallels = np.arange(-90.,90.,4.)
    meridians = np.arange(-180.,181.,10.)
    map.drawparallels(parallels,labels=[False,False,False,True])
    map.drawmeridians(meridians,labels=[True,False,False,False])
    map.drawcoastlines()
    map.fillcontinents()
    map.drawmapboundary()
    
    x1 = iceDF.lon.values.T.tolist()
    y1 = iceDF.lat.values.T.tolist()
    z1 = iceDF.d18o.values.T.tolist()
    x2, y2 = map(x1, y1)
    map.scatter(x2,y2,c=z1,marker='o',s=16,linewidth=.2)
    plt.title('Arctic Frozen Season Observed $\delta^{18}$O',fontsize=14)
    plt.clim(-5.6,2.8)
    cbar = plt.colorbar(orientation='vertical')
    cbar.ax.set_xlabel('$\delta^{18}$O')
    if(save == 1):
        plt.savefig('frozenObs', dpi = 300)
    plt.show()
    
    # Plot 2
    plt.figure(figsize=(7,7))  
    map = Basemap(projection='npstere',boundinglat=68,lon_0=330,resolution='l')
    map.drawparallels(parallels,labels=[False,False,False,True])
    map.drawmeridians(meridians,labels=[True,False,False,False])
    map.drawcoastlines()
    map.fillcontinents()
    map.drawmapboundary()
    
    x1 = waterDF.lon.values.T.tolist()
    y1 = waterDF.lat.values.T.tolist()
    z1 = waterDF.d18o.values.T.tolist()
    x2, y2 = map(x1, y1)
    map.scatter(x2,y2,c=z1,marker='o',s=16,linewidth=.2)
    plt.title('Arctic Melt Season Observed $\delta^{18}$O',fontsize=14)
    plt.clim(-5.6,2.8)
    cbar = plt.colorbar(orientation='vertical')
    cbar.ax.set_xlabel('$\delta^{18}$O')
    if(save == 1):
        plt.savefig('meltObs', dpi = 300)
    plt.show()

