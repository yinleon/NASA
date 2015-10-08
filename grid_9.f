      PROGRAM fd18os3
      !to run in terminal : gfortran -ffixed-form grid_d18o_7_9_15.f followed by ./a.out
      IMPLICIT NONE
      INTEGER,PARAMETER :: skip1=-99,nfil=25514,FILTNT=1,test=1
     *     ,ii=10000,NTTEST=102,REGCT=29,levdep=0,imWOA=360
     *     ,jmWOA=180,lmWOA=33, arcCount=5135
      REAL,PARAMETER :: skip2=-99.99,saltyNA=36.5
      CHARACTER*100 FILEIN,d18OFILE,NOTEST,REFT,phosfile,salfile,ESTFILE
     *     ,tempfile,o2file,mldfile,salout,testout,FILEOUT
      CHARACTER*80 titfort
      character*4 clev(lmWOA)
      character*8 basincode(regct)
      REAL*4 tempWOA(imWOA,jmWOA),WOAs(imWOA,jmWOA,lmWOA),WOAt(imWOA
     *     ,jmWOA,lmWOA),WOAo(imWOA,jmWOA,lmWOA),WOAp(imWOA,jmWOA,lmWOA)
     *     ,WOApstar(imWOA,jmWOA,lmWOA),Levmld(imWOA,jmWOA,13),WOAskip
     *     ,fortest(1:imWOA,1:jmWOA),Levskip,pstar
     *     ,fortest1(1:imWOA,1:jmWOA)
      INTEGER I,J,K,ll,L,P,COUNT(REGCT),as(REGCT,6),grid,gridlp
     *     ,DRAD,MONT,YEART,DEPTHT,klon,klat,lmin,lmax
     *     ,TCT,WorkM(imWOA,jmWOA),WOAlev(lmWOA)
     *     ,levcount,kfil,GCOUNT,KCOUNT,PCOUNT
     *     ,K1,KC,M,N,GS,GI,GJ,LAGO,wdep,geotct(regct),n_northatl
     *     ,arcLon,arcLat
      REAL*4 LAT(REGCT,II),LON(REGCT,II),SAL(REGCT,II),d18O(REGCT,II)
     *     ,DEP(REGCT,II),SLOPE(REGCT),ICEPT(REGCT),DIST
     *     ,SIGSL(REGCT),sigic(regct),rsq(regct),NLON,MLAT
     *     ,GLON,GLAT,TLON,TLAT,LONT,LATT,PTEMPT,SALT,d18OT,dDT,ESTd18O
     *     ,WT,GDISTxy,GDISTz,SUME,SUMO,SUMWTE,SUMWTO,KEEPEST,CALCD18O
     *     ,DISTxy,DISTz,geormsS(regct,ii),geormsO(regct,ii)
     *     ,geormsERR(regct)
     *     ,NAsalmax(2),NAo18max(2),NAsalSLOPE,NAsalICEPT
     *     ,WANT,TTVEC(NTTEST),bdWOAlev(lmWOA+1),phosNAtl,phosmix
     *     ,phosAABW,phosNADW,omit,omit2,omit3,omit4
c     *     ,SYX(REGCT),QFIT(REGCT),TVAL,g1,g2,gln1,x,a,gammq,junk
      integer n_basin(8),n_mostNADW,n_mostAABW,n_NADWAABW,n_deeparco
     *     ,n_dpmedite,n_depacind,n_intnopac,n_AACW_CDW
c     Functions
      REAL*4 smoother,WTSC,wtscz,CALCDIST
c     *     CAL2DIST
      INTEGER waterms

      DATA WOAlev /0, 10, 20, 30, 50, 75, 100, 125, 150, 200, 
     *     250, 300, 400, 500, 600, 700, 800, 900,1000,1100,
     *     1200, 1300, 1400, 1500, 1750, 2000, 2500, 3000, 3500,
     *     4000, 4500, 5000, 5500/    

      DATA clev /"0","10","20","30","50","75","100","125","150","200"
     *     ,"250","300","400","500","600","700","800","900","1000"
     *     ,"1100","1200","1300","1400","1500","1750","2000","2500"
     *     ,"3000","3500","4000","4500","5000","5500"/    

      DATA bdWOAlev /0.,5.,15.,25.,40.,62.5,87.5,112.5,137.5,175.,225.,
     *     275.,350.,450., 550., 650.,750.,850.,950.,1050.,
     *     1150.,1250.,1350.,1450.,1625.,1875.,2250.,2750.,3250.,3750.,
     *     4250.,4750.,5250.,6000./
      phosNADW=0.80 ! minus .14
      phosAABW=1.88 ! plus .14

c      a=194. ;  x=4.82040787
c      call gcf(g1,a,x,gln)
c      call gser(g2,a,x,gln1)
c      print*,a,x,gammq(a,x),g1,g2,gln,gln1
c      stop

**********************GIVE INPUT INSTRUCTIONS***************************
      ! Read in startup files:
      FILEIN="dbO18.txt"
      FILEOUT= "dbO18_CALC_OG"
      ESTFILE = "dbO18_EST_OG"
      salfile="salinityWOA01.s4"
      tempfile="temperatureWOA01.s4"
      o2file="oxygenWOA01.s4"
      phosfile="phosphateWOA01.s4" ! 1x1 0.5E 89.5S
!      phosfile="po4starLev94.s4"
      mldfile="mixedlayerLevitus94.s4" ! 5X2  30E 60N

      if(test.EQ.1) then
         testout='testout'
         open(900,FILE=testout,FORM="UNFORMATTED")
      endif
      
      print*,"Opening salinity file : ",trim(salfile)
      OPEN(200,FILE=trim(salfile),FORM="UNFORMATTED",access=
     *     "SEQUENTIAL",convert="big_endian")
      print*,"Opening temperature file file : ",trim(tempfile)
      OPEN(201,FILE=trim(tempfile),FORM="UNFORMATTED",access=
     *     "SEQUENTIAL",convert="big_endian")
      print*,"Opening o2 file : ",trim(o2file)
      OPEN(202,FILE=trim(o2file),FORM="UNFORMATTED",access=
     *     "SEQUENTIAL",convert="big_endian")
      print*,"Opening phosphate file : ",trim(phosfile)
      OPEN(203,FILE=trim(phosfile),FORM="UNFORMATTED",access=
     *     "SEQUENTIAL",convert="big_endian")
      print*,"Opening MLD file : ",trim(mldfile)
      OPEN(204,FILE=trim(mldfile),FORM="UNFORMATTED",access=
     *     "SEQUENTIAL",convert="big_endian")
      print*,"Opening salinity file : WOAtemperature"
!      OPEN(205,FILE="WOAtemperature",FORM="UNFORMATTED",access=
!     *     "SEQUENTIAL",convert="big_endian")

      do k=1,lmWOA
        read(200,end=9999) tempWOA(1:imWOA,1:jmWOA) !WOAs(1:imWOA,1:jmWOA,k)
        !print*,"A"
        !goto 9999

        WOAs(1:180,1:jmWOA,k)=tempWOA(181:imWOA,1:jmWOA)
        WOAs(181:imWOA,1:jmWOA,k)=tempWOA(1:180,1:jmWOA)
        read(201,end=9999) tempWOA(1:imWOA,1:jmWOA) !WOAt(1:imWOA,1:jmWOA,k)
        WOAt(1:180,1:jmWOA,k)=tempWOA(181:imWOA,1:jmWOA)
        WOAt(181:imWOA,1:jmWOA,k)=tempWOA(1:180,1:jmWOA)
        !print*,"B"
        read(202,end=9999) tempWOA(1:imWOA,1:jmWOA) !WOAo(1:imWOA,1:jmWOA,k)
        WOAo(1:180,1:jmWOA,k)=tempWOA(181:imWOA,1:jmWOA)
        WOAo(181:imWOA,1:jmWOA,k)=tempWOA(1:180,1:jmWOA)
       ! print*,"C"
        read(203,end=9999) tempWOA(1:imWOA,1:jmWOA) !WOAp(1:imWOA,1:jmWOA,k)
        WOAp(1:180,1:jmWOA,k)=tempWOA(181:imWOA,1:jmWOA)
        WOAp(181:imWOA,1:jmWOA,k)=tempWOA(1:180,1:jmWOA)
        !print*,"D"
!        salout="WOA 2001 Temperature Level "//trim(clev(k))
!        write(205) salout,WOAt(1:imWOA,1:jmWOA,k)
      enddo
!      goto 9999
      print*,"Going"
      !goto 9999



      print*,"Captured WOA 2001 salinity,temperature,oxygen,phosphate"
      do k=1,12
        read(204,end=9999) tempWOA(1:imWOA,1:jmWOA)
        Levmld(1:180,1:jmWOA,k)=tempWOA(181:imWOA,1:jmWOA)
        Levmld(181:imWOA,1:jmWOA,k)=tempWOA(1:180,1:jmWOA)
      enddo
      print*,"Captured Levitus 1994 mixed layer depth"
      close(200)
      close(201)
      close(202)
      close(203)
      close(204)
!      OPEN(205,FILE="LevAnnMLD",FORM="UNFORMATTED")
ccccc Convert PO4 to PO4* (Broecker,1998) & UNESCO ml/l -> umol
ccccc Also calculate ann average mld for general sal based est
      WOAskip=WOAs(1,1,1)
      levskip=levmld(1,1,1)
c      print*,"WOA skip set to : ",WOAskip      

      do i=1,imWOA
        do j=1,jmWOA
          do k=1,lmWOA
            if((WOAp(i,j,k).eq.WOAskip).or.(WOAo(i,j,k).eq.WOAskip))
     *           then
              WOApstar(i,j,k)=WOAskip
            else
              WOApstar(i,j,k)=WOAp(i,j,k) + (44.61*WOAo(i,j,k))/175. - 
     *             1.95
            endif
          enddo
          levcount=0
          levmld(i,j,13)=0
          do k=1,12
            if((WOAs(i,j,1).eq.WOAskip) .or. .not.((levmld(i,j,k).lt.0)
     *           .or.(levmld(i,j,k).ge.0)) ) then 
c              Levmld(i,j,13)=WOAskip    
              levmld(i,j,k)=WOAskip
            endif
            if(levmld(i,j,k).ne.WOAskip) then
              levcount=levcount+1
              Levmld(i,j,13)=Levmld(i,j,13)+Levmld(i,j,k)
            else
              levcount=levcount
              Levmld(i,j,13)=Levmld(i,j,13)
c     print*,i,j,k,levmld(i,j,k)
            endif
            
          enddo 
          if(levcount.gt.0) then
            levmld(i,j,13)=levmld(i,j,13)/levcount
          elseif(WOAs(i,j,1).ne.WOAskip) then
            levmld(i,j,13)=0.
          else
            levmld(i,j,13)=WOAskip
          endif
          if(levmld(i,j,13).lt.0) levmld(i,j,13)=WOAskip
          if((WOAs(i,j,1).ne.WOAskip).and.(levcount.eq.0)) then
            levmld(i,j,13)=0
          endif
c          if((levmld(i,j,13).ne.WOAskip).and.(levmld(i,j,13).lt.0)) 
c     *         then
c            print*,i,j,levmld(i,j,13),levcount
c            do k=1,12
c              print*,k,levmld(i,j,k)
c            enddo
c          endif
          if((levmld(i,j,13).gt.100).and.(j.lt.55)) levmld(i,j,13)=100
        enddo
      enddo

      do i=1,imWOA
        do j=1,jmWOA
          if(levmld(i,j,13).ne.WOAskip) then
            levmld(i,j,13)=smoother(i,j,imWOA,jmWOA,levmld(1:imWOA
     *           ,1:jmWOA,13),5,WOAskip)
          else
            levmld(i,j,13)=WOAskip
          endif
        enddo
      enddo

!      print*,"MLD ",levmld(30,30,13)
      salout="Levitus94 ANN MLD"
!      write(205) salout,Levmld(1:imWOA,1:jmWOA,13)
!      close(205)
!      goto 9999

      call basindef(workm,as,basincode,n_basin,n_northatl)
      n_mostNADW=n_basin(1)
      n_mostAABW=n_basin(2)
      n_NADWAABW=n_basin(3)
      n_deeparco=n_basin(4)
      n_dpmedite=n_basin(5)
      n_depacind=n_basin(6)
      n_intnopac=n_basin(7)
      n_aacw_cdw=n_basin(8)

!      print*,workm(1,1),workm(imWOA,jmWOA),workm(1,90),workm(160,90)
!      goto 9999
************************************************************************
****  ONLY KEEP FILES WITH NO 'NOTES',SORT DOWN TO A DEPTH
c     RANGE, SORT DOWN TO A PARTICULAR OCEAN BASIN  ****************
      PRINT*,"Opening and filtering database file ...", FILEIN
      phosNAtl=1.1 ! mostly NADW end member
      phosmix=1.7 !1.4 = equal NADW&AABW      

      do i=1,regct
        count(i)=0
        geotct(i)=0
      enddo

 10   FORMAT (F7.2, F6.2, I2, I4, I5, 3F6.2, F6.1, A25, A50)
 20   FORMAT (F7.2, F6.2, I2, I4, I5, 3F6.2, F6.1, A50)
 17   FORMAT (F8.2, F6.2, I3, I5, I5, 3F6.2, F6.1, 2F9.3, A50)
      OPEN(1,FILE=FILEIN,FORM="FORMATTED")
      OPEN(206,FILE="dbsort.dat",FORM="FORMATTED")
      do k=1,regct
        kfil=210+k
        open(kfil,file="dbsort.dat."//basincode(K),FORM="FORMATTED")
      enddo

cccc Do special salty north atlantic
      NAsalmax(1)=0
      NAo18max(1)=0
      NAsalmax(2)=0
      NAo18max(2)=0

      DO I=1,nfil
        READ (UNIT=1,FMT=10) LONT,LATT,MONT,YEART,DEPTHT,PTEMPT,SALT,
     *       d18OT,dDT, NOTEST, REFT
        if(abs(lont).gt.180) goto 152
! skip records with notes or missing d18O or missing salinity
        IF( (INDEX(NOTEST,'L').GE.1).or.(INDEX(NOTEST,'X').GE.1).or.
     *       (INDEX(NOTEST,'K').GE.1).or.(INDEX(NOTEST,'A').GE.1).or.
     *       (INDEX(NOTEST,'B').GE.1).or.(INDEX(NOTEST,'S').GE.1).or.
     *       (d18OT.LE.-90.0).or.(SALT.LE.0.0) ) GOTO 152 
        
        KLON=INT(180+AINT(LONT))
        KLAT=INT(90+AINT(LATT))         
        if(klat.eq.181) klat=180 ! just in case its on the north pole

c        if((lont.eq.-60.18).and.(latt.eq.46.91)) print*,klon,klat
c     *       ,workm(klon,klat)

        K=workm(klon,klat)

        IF(K.EQ.0) then ! try out surrounding boxes
          K=workm(klon+1,klat)
          if(K.ne.0) goto 159
          K=workm(klon-1,klat)
          if(K.ne.0) goto 159
          K=workm(klon,klat-1)
          if(K.ne.0) goto 159
          K=workm(klon,klat+1)
          if(K.ne.0) goto 159
          GOTO 152
        endif

 159    continue

        if((k.eq.n_northatl).and.(SALT.gt.NAsalmax(2))) then
          NAsalmax(2)=SALT
          NAo18max(2)=d18OT
        endif

        if((deptht.le.50).or.(salt.lt.30).or.
     *         (index(basincode(K),'BaffinBa').gt.0).or.
     *         (index(basincode(K),'HudsonBa').gt.0).or.
     *         (index(basincode(K),'BalticSe').gt.0).or.
     *         (index(basincode(K),'LabdorSe').gt.0)) then
          K=K  
        elseif ((mont.lt.1).or.(mont.gt.12)) then ! need month to det mld
          if(deptht.ge.1000) then ! for really deep, month doesn't matter 
            mont=6
          else
            goto 152
          endif
        else
          wdep=1
          do m=2,lmWOA
            if(deptht.gt.bdWOAlev(m)) wdep=m
          enddo
          pstar=WOApstar(klon,klat,wdep)
          if(pstar.le.-10) then
            pstar=WOApstar(klon+1,klat,wdep)
            if(pstar.gt.-10) goto 160
            pstar=WOApstar(klon-1,klat,wdep)
            if(pstar.gt.-10) goto 160
            pstar=WOApstar(klon,klat+1,wdep)
            if(pstar.gt.-10) goto 160
            pstar=WOApstar(klon,klat-1,wdep)
            if(pstar.gt.-10) goto 160
            pstar=WOApstar(klon,klat,wdep-1)
            if(pstar.gt.-10) goto 160
            GOTO 152
          endif
 160      continue
          K=waterms(K,deptht,Levmld(klon,klat,mont)
     *         ,pstar,phosnatl,phosmix,basincode,n_basin)
        endif

c        elseif((deptht*1.).le.Levmld(klon,klat,mont)) then
c          K=K                   ! surface regions
c        elseif((index(basincode(K),'ArcticOc').gt.0).or
c     *         .(index(basincode(K),'GINSeAtl').gt.0)) then
c          K=n_deeparco          ! Deep Arctic
c        elseif(index(basincode(K),'Mediterr').gt.0) then
c          K=n_dpmedite          ! Mediterrannean (deep)
c        elseif(index(basincode(K),'NorthAtl').gt.0) then
c          K=n_mostnadw
c        elseif( (WOApstar(klon,klat,wdep).ge.phosmix).and
c     *         .(index(basincode(K),'NorthPac').le.0)) then
c          K=n_mostAABW
c        else ! decide which endmember it is closest to if its below mld
c          if(index(basincode(K),'Atl') .gt.0) then
c            if(WOApstar(klon,klat,wdep).le.phosNatl) then
c              K=n_mostnadw      ! NADW endmember
c            else
c              K=n_nadwaabw
c            endif
c          elseif(index(basincode(K),'SouthOcn').gt.0) then
c            if((deptht.gt.3000).or.(WOApstar(klon,klat,wdep).ge.phosmix)
c     *           ) then
c              K=n_mostAABW
c            else
c              K=n_aacw_cdw
c              K=n_nadwaabw
c            endif
c          elseif (deptht.gt.2500) then
c            k=n_dePacInd
c          elseif (index(basincode(K),'NorthPac').gt.0) then
c            K=n_intnopac        ! North Pacific Intermediate Water
c          else                  ! Intermediate Pacific Ocean, Indian Ocean
c            K=n_nadwaabw        ! Mixed NADW and AABW 
c          endif                
c        endif                  

 151    continue
        kfil=210+k
        WRITE(UNIT=kfil,FMT=17)
     *       LONT,LATT,MONT,YEART,DEPTHT,PTEMPT,SALT,d18OT,dDT,
     *       levmld(klon,klat,mont),pstar,REFT
        
        COUNT(K)=COUNT(K)+1
        TCT=COUNT(K)
        SAL(K,TCT)=SALT
        d18O(K,TCT)=d18OT
        LON(K,TCT)=LONT
        LAT(K,TCT)=LATT
        DEP(K,TCT)=DEPTHT
        if(index(reft,'GEOSECS').ge.1) then
          geotct(k)=geotct(k)+1
          tct=geotct(k)
          geormsS(k,tct)=SALT
          geormsO(k,tct)=d18OT
        endif
        
 152    CONTINUE

      ENDDO                     !I=1,II

      CLOSE(1)
      CLOSE(206)
!      goto 9999
      print*,"Database file sorted"

!ccc  Load ttest values 

      WANT=0.9750 ! Confidence interval desired (1-alpha/2)
      CALL LDTTEST(WANT,NTTEST,TTVEC)

**********GO THROUGH ALL OCEAN POINTS AND FIND d18O-S RELATIONSHIP******

  45  FORMAT(F6.1,F6.1,F8.4,I6,F8.4,F10.4,F8.4,F10.4)
c            Long  Lat  Sal N Slope Icept estd18O calcd18O 

  49  FORMAT (F6.1, F6.1, F8.4, I4)
  48  FORMAT (F6.1, F6.1, F8.4, F8.4, I5, A50) 
  46  FORMAT (A8,I8,2F8.2,4F8.3,I6)

      PRINT*, "*************Starting Computation****************"

C Zero out counters
      I=0
      J=0
      K=0
      CALCD18O=skip2

************Calculate O18-S Relationship in each region*****************
      TCT=0
      PRINT*," Region   Count   Slope   ICept  SigSlope SigIcept R-sq"
      DO I=1,REGCT
         TCT=COUNT(I)
         if(tct.gt.7) then
           omit=0.
           omit2=-20.
           omit3=50.
           omit4=10.
           if((index(basincode(I),'TropiAtl').gt.0).or
     *          .(index(basincode(I),'IndianOc').gt.0).or
     *          .(index(basincode(I),'GIN').gt.0)) omit=28.
           if(index(basincode(I),'mostNADW').gt.0) then
             omit=34.
             omit2=-0.8
           endif
           if(index(basincode(I),'NorthAtl').gt.0) omit2=-2.5
           if(index(basincode(I),'SouthOcn').gt.0) omit4=0.
           if(index(basincode(I),'AACW_CDW').gt.0) omit3=34.8
           CALL FIT(SAL(I,1:TCT),d18O(I,1:TCT),TCT,ICEPT(I),SLOPE(I)
     *          ,skip2,sigIC(I),sigSL(I),rsq(I),omit,omit2,omit3,omit4)
           CALL calcrms(geormsS(i,1:geotct(i)),geormsO(i,1:geotct(i))
     *          ,geotct(i),icept(i),slope(i),skip2,geormsERR(i))
!         CALL LSQ(SAL(I,1:TCT),d18O(I,1:TCT),TCT,SLOPE(I),ICEPT(I))
!         CALL XQFIT(SAL(I,1:TCT),d18O(I,1:TCT),SLOPE(I),ICEPT(I),TCT,
!     *     QFIT(I))
!         CALL STATS(SAL(I,1:TCT),d18O(I,1:TCT),SLOPE(I),ICEPT(I),TCT,
!     *        TVAL,GOODP,SIGSL(I),SYX(I))
         endif
         PRINT46,basincode(I),COUNT(I),SLOPE(I),ICEPT(I),sigsl(I)
     *        ,sigIC(I),rsq(I),geormsERR(i),geotct(i)
         !endif
      END DO                    ! I=1,REGCT

cccc Do Special salty north atlantic
      NAsalmax(1)=saltyNA
      NAo18max(1)=slope(n_northatl)*saltyNA+icept(n_northatl)
      NAsalSLOPE=(NAo18max(2)-NAo18max(1))/(NAsalmax(2)-NAsalmax(1))
      NAsalICEPT=NAo18max(2)-NAsalSLOPE*NAsalmax(2)
      print*,"Salty NA mins : ",NAsalmax(1),NAo18max(1)
      print*,"Salty NA maxs : ",NAsalmax(2),NAo18max(2)
      print*,"Salty NA sl ic : ",NAsalSLOPE,NAsalICEPT
************************************************************************
!      goto 9999

*************************Make Smooth d18O*******************************
************************************************************************
  47  FORMAT (I4,I4,F7.1,F6.1,I4,I4,F7.2,F8.4,F10.4,F8.4)
 447  FORMAT (I4,I4,F7.1,F6.1,F7.1,F6.1,I4,I4,F7.2,3F8.4,F10.4)

      grid=5
      DRAD=grid*111 ! approximately 5 grid boxes at equator
      gridlp=2*grid+1

      wdep=0

      open(210,file=FILEOUT,form="unformatted")
      open(220,file=ESTFILE,form="unformatted")
!      wdep=5 ! 75m
!      wdep=28 ! 3500m
!      wdep=13 ! 500m
!      wdep=6 ! 100m
!      wdep=7
 110  wdep=wdep+1

      WOAs(imWOA,jmWOA,wdep)=WOAs(imWOA,jmWOA-1,wdep) ! North Pole undefined

c      if (wdep.eq.14) goto 9999
      if (wdep.eq.34) goto 9999

      kfil=300+wdep
      d18OFILE="d18Og3_sv10.txt."//trim(clev(wdep))
      PRINT*, "Output File: ",trim(d18OFILE)
      OPEN(kfil,FILE=d18OFILE,FORM="FORMATTED")
      
      Print*,"Calculating... ... ... ... Level ",clev(wdep)
      Print*,"  Lon   Lat Salinity  PCt   Slope     Icept KeepEst Calcd1
     *8O"
      TLAT=-90.5
      DO J=1,jmWOA
        TLON=-180.5
        TLAT=TLAT+1             ! DEFINE TLAT LIKE IN LEVITUS(94)
c      TLAT=TLAT+174
c      DO J=175,jmWOA
c      DO J=179,179
c      TLAT=TLAT+J

        DO I=1,imWOA
 100      TLON=TLON+1           ! DEFINE TLON LIKE IN LEVITUS(94)
c      DO I=348,348
c      TLON=TLON+I
!*****   Determine what workzone it is...
c          print*,i,j,tlon,tlat,workm(i,j),basincode(workm(i,j))          
          KC=waterms(workm(i,j),WOAlev(wdep),Levmld(i,j,13),
     *         WOApstar(i,j,wdep),phosNatl,phosmix,basincode,n_basin)
          IF((KC.EQ.0).or.(WOAs(i,j,wdep).eq.WOAskip)) goto 945
!          KC=WorkM(I,J)
!          print*,i,j,tlon,tlat,KC,basincode(KC)          

          PCOUNT=0
          GCOUNT=0
          SUME=0.
          SUMO=0.
          SUMWTE=0.
          SUMWTO=0.
          
          if(WOAlev(wdep).lt.(2*levmld(i,j,13))) then
            lmin=1
            lmax=1
c            lmax=wdep
            do m=2,lmWOA
              if((2*levmld(i,j,13)).gt.bdWOAlev(m)) lmax=m
            enddo
          else
            lmin=wdep
            lmax=wdep
          endif
            
          KCOUNT=as(KC,1)+1 ! as(*,1) is the number of compatible reg.
!********
          IF(J.GT.(jmWOA-gridlp)) THEN ! Special Arctic Ocean stuff
c            GS=jmWOA-gridlp
            gs=J-grid
c            GLAT=-90.5+(jmWOA-gridlp)-1
            glat=TLAT-grid-1
            DO GJ=GS,jmWOA
              GLAT=GLAT+1
              GLON=-180.5
              DO GI=1,imWOA
                GLON=GLON+1
                  GDISTxy=CALCDIST(TLON,TLAT,GLON,GLAT)
                  IF(GDISTxy.GT.DRAD) GOTO 155
                DO ll=lmin,lmax
                  if(ll.ne.wdep) then
                    GDISTz=abs(WOAlev(wdep)-WOAlev(ll))
                  else
                    GDISTz=0
                  endif
                  if (WOAs(gi,gj,ll).eq.WOAskip) goto 154
c     PRINT*,tlon,tlat,GLON,GLAT,GDIST
c                  K1=waterms(workm(gi,gj),WOAlev(wdep),Levmld(gi,gj,13)
c     *                 ,WOApstar(gi,gj,wdep),phosNatl,phosmix,basincode
c     *                 ,n_basin)
c                  K1=waterms(workm(gi,gj),WOAlev(ll),Levmld(gi,gj,13)
c     *                 ,WOApstar(gi,gj,ll),phosNatl,phosmix,basincode
c     *                 ,n_basin)
                  if(WOAlev(ll).le.levmld(gi,gj,13)) then
                    K1= WorkM(GI,GJ)
                  else
                    K1=19
                  endif
                  
                  IF (K1.eq.0) GOTO 154
*******Use only compatible regions*****
c                LAGO=0
c                DO L=2,KCOUNT
c                  IF(K1.eq.as(KC,L)) LAGO=1
c                END DO
c                IF(LAGO.eq.0) GOTO 154
***************************************
                  GCOUNT=GCOUNT+1
                  WT=WTSC(GDISTxy)*WTSCz(GDISTz,levmld(gi,gj,13))
                  ESTd18O=WOAs(gi,gj,ll)*SLOPE(K1)+ICEPT(K1)
                  WRITE(900)ESTd18O
                  SUME=SUME+WT*ESTd18O
                  SUMWTE=SUMWTE+WT
c                  PRINT447,GI,GJ,GLON,GLAT,TLON,TLAT,K1,GCOUNT,GDISTxy
c     *                 ,WT,sumwte,WOAs(gi,gj,ll),ESTd18O
                  if((wt.lt.0).or.(wt.gt.1)) PRINT47,GI,GJ,GLON
     *                 ,GLAT,K1,GCOUNT,GDISTxy,WT,ESTd18O,SUMWTE
                  if(sumwte.lt.0) PRINT47,GI,GJ,GLON
     *                 ,GLAT,K1,GCOUNT,GDISTxy,WT,ESTd18O,SUMWTE
 154              CONTINUE
                enddo
 155            continue
              END DO            ! JJ
            END DO              ! II

          else                !(I.GT.180-gridlp) THEN special Arctic
!**************          

!          PRINT*,"  M   N  NLON   MLAT  K1 GCT   DIST      WT   ESTd18 
!     *O   SUMWT"
          DO M=-grid,grid ! Calculate d18O based on nearby sal-d18O slope
            DO N=-gridlp,gridlp
                NLON=TLON+N
                MLAT=TLAT+M
                DISTxy=CALCDIST(TLON,TLAT,NLON,MLAT)
                IF (DISTxy.gt.DRAD) GOTO 157
              do ll=lmin,lmax
                if(ll.ne.wdep) then
                  DISTz=abs(WOAlev(wdep)-WOAlev(ll))
                else
                  DISTz=0
                endif
c              print*,nlon,mlat,dist,drad,workm(i-n,j-m),woas(i-n,j-m
c     *             ,wdep)

                if (WOAs(i-n,j-m,ll).eq.WOAskip) goto 153
                K1=waterms(workm(i-n,j-m),WOAlev(wdep),Levmld(i-n,j-m,13
     *               ),WOApstar(i-n,j-m,wdep),phosNatl,phosmix,basincode
     *               ,n_basin)

                IF (K1.eq.0) GOTO 153

*******Use only compatible regions*****
                LAGO=0
                DO L=2,KCOUNT
                  IF(K1.eq.as(KC,L)) LAGO=1
                END DO
                IF(LAGO.eq.0) GOTO 153
***************************************
                K1=waterms(workm(i-n,j-m),WOAlev(ll),Levmld(i-n,j-m,13)
     *               ,WOApstar(i-n,j-m,ll),phosNatl,phosmix,basincode
     *               ,n_basin)
                GCOUNT=GCOUNT+1
                WT=WTSC(DISTxy)*wtscz(distz,levmld(i-n,j-m,13))
                ESTd18O=WOAs(I-N,J-M,ll)*SLOPE(K1)+ICEPT(K1)
                if((k1.eq.n_northatl).and.(WOAs(I-N,J-M,ll).gt.saltyNA))
     *               ESTd18O=WOAs(i-N,J-M,ll)*NAsalSLOPE+NAsalICEPT
                SUME=SUME+WT*ESTd18O
                SUMWTE=SUMWTE+WT
 153            CONTINUE
              enddo ! ll=lmin,lmax
 157          continue
            END DO              ! N=-gridlp,gridlp
          END DO                ! M=-grid,grid
          endif
c 156      CONTINUE
          KEEPEST=WOAs(I,J,wdep)*SLOPE(KC)+ICEPT(KC)
          if((KC.eq.n_northatl).and.(WOAs(I,J,ll).gt.saltyNA))
     *         KEEPEST=WOAs(I,J,ll)*NAsalSLOPE+NAsalICEPT
ccccccc Blend in actual points
          DO P=2,KCOUNT ! kcount=number of compatible regions 
c            print*,kc,p,as(kc,p)
            K=as(KC,P)          ! Go through compatible regions
c            print*,KC,P,as(KC,P),K
c            goto 9999
            DO L=1,COUNT(K)
c              print*,basincode(l)
              DIST=CALCDIST(TLON,TLAT,LON(K,L),LAT(K,L))
              IF (DIST.GT.DRAD) GOTO 166
              IF ((DEP(K,L).lt.bdWOAlev(wdep)).or.(DEP(K,L).gt.
     *             bdWOAlev(wdep+1))) goto 166
              if((tlon.eq.-67.5).and.(tlat.eq.-79.5)) PRINT*,K,L,LON(K,L
     *             ),LAT(K,L),dist,d18O(K,L)
              PCOUNT=PCOUNT+1
              WT=WTSC(DIST)
              SUMO=SUMO+WT*d18O(K,L)
              SUMWTO=SUMWTO+WT
              if((wt.lt.0).or.(wt.gt.1)) PRINT47,I,J,LON(K,L)
     *             ,LAT(K,L),K,PCOUNT,DIST,WT,d18O(K,L),SUMWTO
              if(sumwto.lt.0) PRINT47,I,J,LON(K,L)
     *             ,LAT(K,L),K,PCOUNT,DIST,WT,d18O(K,L),SUMWTO
 166          CONTINUE
            END DO              !L=1,COUNT
          END DO                ! P=1,KCOUNT
!***** Weight points at 90%, calculated at 10%
          CALCD18O=(0.9*SUMO+0.1*SUME)/(0.9*SUMWTO+0.1*SUMWTE)
c          if(sumo.gt.0) 
c     print*,i,j,tlon,tlat,sume,sumwte,sumo,sumwto          
          WRITE (UNIT=kfil,FMT=45) TLON,TLAT,WOAs(I,J,wdep),PCOUNT
     *         ,SLOPE(KC),ICEPT(KC),KEEPEST,CALCD18O
c          print45,TLON,TLAT,WOAs(I,J,wdep),PCOUNT
c     *         ,SLOPE(KC),ICEPT(KC),KEEPEST,CALCD18O
c     if (calcd18o.lt.-5) print*,tlon,tlat,calcd18o
          fortest(i,j)=calcd18O
          fortest1(i,j)=KEEPEST
 945      CONTINUE
          IF ((KC.EQ. 0).or.(WOAs(i,j,wdep).eq.WOAskip)) THEN
            WRITE (UNIT=kfil,FMT=45) TLON,TLAT,0.,0,0.,0.,0.,0.
            fortest(i,j)=-1.E30
          END IF ! (KC.EQ. 0).or.(WOAs(i,j,wdep).eq.WOAskip)) THEN
        END DO                  !I=1,360
        if ((abs(tlat).eq.70.5).or.(abs(tlat).eq.50.5).or.(abs(tlat).eq
     *       .30.5).or.(abs(tlat).eq.10.5)) PRINT45
     *       ,TLON,TLAT,WOAs(I,J,wdep),PCOUNT,SLOPE(KC),ICEPT(KC)
     *       ,KEEPEST,CALCD18O
      END DO                    !J=1,180

      CLOSE(kfil)
      titfort="O18water based on database S-O18 and WOAs, Level "/
     *     /trim(clev(wdep))
      write(210) titfort,fortest(1:imWOA,1:jmWOA)  !calculated
      write(220) titfort,fortest1(1:imWOA,1:jmWOA) !estimated
      print*,"********* Finished  level ",trim(clev(wdep)),"***********"
      print*," "
      if(wdep.eq.1) then
        open(555,file="calculated_d18O_level0_sv10",form="unformatted")
        write(555) titfort,fortest(1:imWOA,1:jmWOA)
        close(555)
      endif
      goto 110

 9999 continue
      close(210)
      close(220)
      
      END PROGRAM fd18os3

************************************************************************
***********************SUBROUTINES AND FUNCTIONS************************
************************************************************************
      integer function waterms(workm,dep,mld,pstar,phosNatl,phosmix
     *     ,basincode,n_basin)
      implicit none
      integer,parameter :: regct=29
      real*4 mld,pstar,phosNatl,phosmix
      integer dep,workm
      character*8 basincode(regct)
      integer n_basin(8),n_mostNADW,n_mostAABW,n_NADWAABW,n_deeparco
     *     ,n_dpmedite,n_depacind,n_intnopac,n_aacw_cdw

      n_mostNADW=n_basin(1)
      n_mostAABW=n_basin(2)
      n_NADWAABW=n_basin(3)
      n_deeparco=n_basin(4)
      n_dpmedite=n_basin(5)
      n_depacind=n_basin(6)
      n_intnopac=n_basin(7)
      n_aacw_cdw=n_basin(8)

c      print*,"waterms: ",workm,dep,mld,pstar,phosnatl,phosmix
      if((pstar.le.-10).or.(workm.eq.0)) then
        waterms=0
      elseif((dep.le.50).or.(dep.le.mld).or.
     *       (index(basincode(workm),'BaffinBa').gt.0).or. !what??????
     *       (index(basincode(workm),'HudsonBa').gt.0).or.
     *       (index(basincode(workm),'BalticSe').gt.0).or.
     *       (index(basincode(workm),'LabdorSea').gt.0)) then
        waterms=workm           ! surface regions
      elseif((index(basincode(workm),'ArcticOc').gt.0).or
     *       .(index(basincode(workm),'GINSeAtl').gt.0)) then
        waterms=n_deeparco              ! Deep Arctic
      elseif(index(basincode(workm),'Mediterr').gt.0) then
        waterms=n_dpmedite              ! Mediterrannean (deep)
      elseif(index(basincode(workm),'NorthAtl').gt.0) then
        waterms=n_mostnadw
      elseif( (pstar.ge.phosmix).and.(index(basincode(workm),'NorthPac')
     *       .le.0)) then
        waterms=n_mostAABW
      elseif(index(basincode(workm),'Atl').gt.0) then
        if(pstar.le.phosNatl) then
          waterms=n_mostnadw    ! NADW endmember
        else
          waterms=n_nadwaabw
        endif
      elseif(index(basincode(workm),'SouthOcn').gt.0) then
        if( (dep.gt.3000).or.(pstar.ge.phosmix)) then
          waterms=n_mostAABW
        else
          waterms=n_aacw_cdw
        endif
      elseif (dep.gt.2500) then
        waterms=n_dePacInd
      elseif (index(basincode(workm),'NorthPac').gt.0) then
        waterms=n_IntNoPac
      else
        waterms=n_nadwaabw      ! Mixed NADW and AABW
      endif                 
 9989 continue
      
      return
      end function              !waterms
************************************************************************
      real*4 FUNCTION smoother(i,j,im,jm,data,krad,skip)
      implicit none
      integer i,j,iw,jw,im,jm,krad,ii,jj,smct
      real*4 data(360,180),skip,sm
      
      if((im.gt.360).or.(jm.gt.180).or.(krad.gt.180)) print*,"Issue smoo
     *     ther won't work"

      sm=0.
      smct=0
      do ii=i-krad,i+krad
        if (ii.le.0) then
          iw=im+ii
        elseif (ii.gt.im) then
          iw=ii-im
        else
          iw=ii
        endif
        do jj=j-krad,j+krad
          if (jj.le.0) then
            jw=1-jj
          elseif (jj.gt.jm) then
            jw=jm-(jj-jm)
          else
            jw=jj
          endif
          if((iw.gt.im).or.(iw.lt.0).or.(jw.gt.jm).or.(jw.lt.0)) then
            print*,"i ",i,ii,iw,"j ",j,jj,jw
          endif
          if(data(iw,jw).ne.skip) then
            sm=sm+data(iw,jw)
            smct=smct+1
          else
            smct=smct
            sm=sm
          endif
        enddo
      enddo
      if(smct.gt.0) then
        smoother=sm/smct
      else
        smoother=skip
      endif
          
      return
      end function ! smoother
************************DISTANCE CALCULATOR*****************************
      real*4 FUNCTION CALCDIST(LON1,LAT1,LON2,LAT2)
      IMPLICIT NONE
      REAL,PARAMETER :: pi=3.14159,radius=6371.,d2r = pi/180
      real*4 LON1,LON2,LAT1,LAT2,BLAT,BLON,alon,alat,dellat,dellon
      REAL*4 a,b,c
c     *     ,CALCDIST 
      ! convert degree coordinates to radians
      ALAT=LAT1*d2r
      ALON=LON1*d2r
      BLAT=LAT2*d2r
      BLON=LON2*d2r

      dellat=ALAT-BLAT
      dellon=ALON-BLON

      a = (sin((dellat*5d-1)))**2 + cos(alat)
     *     * cos(Blat) * (sin((dellon*5d-1)))**2  
      b=sqrt(a)
      c=sqrt(1d0-a)
      if((lat1.eq.lat2).and.(lon1.eq.lon2)) then
        CALCDIST=1E-5
      else
        CALCDIST=2d0 * atan2(b,c) * radius
      endif

      return
      end FUNCTION
************************************************************************
**********************2 DISTANCE CALCULATOR*****************************
      real*4 FUNCTION CAL2DIST(LON1,LAT1,LON2,LAT2)
      IMPLICIT NONE
      REAL,PARAMETER :: pi=3.14159,radius=6371.,d2r = pi/180
      real*4 LON1,LON2,LAT1,LAT2,BLAT,BLON,alon,alat
      REAL*4 CALCDIST,a,b,c

c      if((lon1.ne.lon2).and.(lat1.ne.lat2)) then
      if((lon1.eq.lon2).and.(lat1.eq.lat2)) then
        CAL2DIST=1E-10
      else
        ALAT=LAT1*d2r
        ALON=LON1*d2r
        BLAT=LAT2*d2r
        BLON=LON2*d2r
        
        a=cos(alat)*cos(alon)*cos(blat)*cos(blon)
        b=cos(alat)*sin(alon)*cos(blat)*sin(blon)
        c=sin(alat)*sin(blat)
        CAL2DIST=ACOS(a+b+c)*radius
      endif

      return
      end FUNCTION
************************************************************************
**************************Scale Vertical********************************
      real*4 FUNCTION WTSCZ(DIST,mld)
      IMPLICIT NONE
      REAL*4 DIST,mld,a,b,c

      if(mld.lt.50) mld=50

      a=DIST/mld
      b=a**4
      c=b/0.01

      WTSCZ=EXP(-c)

      return
      end FUNCTION
************************Scale Horizontal*********************************
      real*4 FUNCTION WTSC(DIST)
      IMPLICIT NONE
      REAL,PARAMETER :: radius=6371.
      REAL*4 DIST,a,b,c

      a=DIST/radius
c      b=a**2 !gav
      b=a**4
      c=b/2E-5 ! grid 5
c      c=0.01 ! gav
c      c=0.0003 ! grid 10

      WTSC=EXP(-c)

      return
      end FUNCTION
************************Scale Horizontal Deep****************************
      real*4 FUNCTION WTSCD(DIST)
      IMPLICIT NONE
      REAL,PARAMETER :: radius=6371.
      REAL*4 DIST,a,b,c

      a=DIST/radius
c      b=a**2 !gav
      b=a**5
      c=b/2E-5 ! grid 5
c      c=0.01 ! gav
c      c=0.0003 ! grid 10

      WTSCD=EXP(-c)

      return
      end FUNCTION

************************************************************************
**************************Subroutine FIT********************************
      SUBROUTINE FIT(X1,X2,N,A,B,skip,siga,sigb,rsq,omit,omit2,omit3
     *     ,omit4)
C     A=ICEPT B=SLOPE
C     CAN ALSO RETURN SIGA,SIGB,CHI2, Q, AND READ IN SIG,MWT
      IMPLICIT NONE
      INTEGER N,PCT,I,MWT
      REAL X(N),Y(N),SIG(N),A,B,sigdat,ss,st2,sx,sx0ss,sy,t,wt,
     *     gammq,X1(N),X2(N),skip,q,siga,sigb,chi2,ax,ay,omit,omit2
     *     ,omit3,omit4,ssy,sse,rsq,sxx,syy,sxy
c     *     ,corrc

c      skip=-1.E+30
      SX=0.
      SY=0.
      ST2=0.
      B=0.
      pct=0
      ay=0.
      ax=0.
      syy=0.
      sxx=0.
      sxy=0.
      ssy=0.
      sse=0.

      do i=1,n
        if( (x1(i).ne.skip).and.(x2(i).ne.skip)) then
          if((x1(i).ge.omit).and.(x2(i).gt.omit2).and.(x1(i).le.omit3)
     *         .and.(x2(i).le.omit4)) then
            pct=pct+1
            X(pct)=X1(i)
            Y(pct)=X2(i)
            sig(pct)=0.2
            ax=ax+x(pct)
            ay=ay+y(pct)
          endif
        endif
      enddo

      ax=ax/pct
      ay=ay/pct

      MWT=0

      IF(MWT.NE.0) THEN
         SS=0.
         DO I=1,PCT
            WT=1./(SIG(I)**2)
            SS=SS+WT
            SX=SX+X(I)*WT
            SY=SY+Y(I)*WT
 11      ENDDO
      ELSE
         DO I=1,PCT
            SX=SX+X(I)
            SY=SY+Y(I)
 12      ENDDO
         SS=FLOAT(PCT)
      ENDIF
      SX0SS=SX/SS
      IF(MWT.NE.0) THEN
         DO I=1,PCT
            T=(X(I)-SX0SS)/SIG(I)
            ST2=ST2+T*T
            B=B+T*Y(I)/SIG(I)
 13      ENDDO
      ELSE
         DO I = 1,PCT
            T=X(I)-SX0SS
            ST2=ST2+T*T
            B=B+T*Y(I)
 14      ENDDO
      ENDIF
      B=B/ST2
      A=(SY-SX*B)/SS
      SIGA=SQRT((1.+SX*SX/(SS*ST2))/SS)
      SIGB=SQRT(1./ST2)
      CHI2=0.
      IF(MWT.EQ.0) THEN
         DO I=1,PCT
            CHI2=CHI2+(Y(I)-A-B*X(I))**2
            ssy=ssy+(Y(I)-AY)**2
            sse=sse+(Y(I)-A-B*X(I))**2
c            syy=syy+(y(i)-ay)**2
c            sxx=sxx+(x(i)-ax)**2
c            sxy=sxy+(y(i)-ay)*(x(i)-ax)
 15      ENDDO
         Q=1.
         SIGDAT=SQRT(CHI2/(PCT-2))
         SIGA=SIGA*SIGDAT
         SIGB=SIGB*SIGDAT
      ELSE
         DO 16 I=1,PCT
            CHI2=CHI2+((Y(I)-A-B*X(I))/SIG(I))**2
            ssy=ssy+(Y(I)-AY)**2
            sse=sse+(Y(I)-A-B*X(I))**2
c            syy=syy+(y(i)-ay)**2
c            sxx=sxx+(x(i)-ax)**2
c            sxy=sxy+(y(i)-ay)*(x(i)-ax)
 16      CONTINUE
       IF(PCT.GT.2) Q=GAMMQ(0.5*(PCT-2),0.5*CHI2)
       ENDIF

       if (ssy.ne.0) rsq=(ssy-sse)/ssy
c       if((sxx.ne.0).and.(syy.ne.0)) rsq=sxy/sqrt(sxx*syy)
c      print*,pct,ay,ssy,sse,ssy-sse

c       if (rsq.ge.0) then
c         corrc=SQRT(rsq)
c       else
c         corrc=0
c       endif
c       if (B .lt. 0 ) corrc=-1*corrc

!       print*,"chi2,pct",chi2,pct,q,GAMMQ(0.5*(PCT-2),0.5*CHI2*25.)
!     *      ,GAMMQ(0.5*(PCT-2),0.5*CHI2)
      RETURN
      END ! subroutine FIT

************************************************************************
**************************RMS Error Calculator**************************    
      SUBROUTINE calcrms(x1,x2,n,a,b,skip,rmserr)
C     A=ICEPT B=SLOPE
      implicit none;
      integer n,pct,i
      real*4 x1(n),x2(n),a,b,x(n),y(n),x3(n),rmserr,skip,sse

      pct=0
      sse=0

      do i=1,n
        if( (x1(i).ne.skip).and.(x2(i).ne.skip)) then
          pct=pct+1
          x(pct)=x1(i)
          y(pct)=x2(i)
          x3(pct)=A+x(pct)*B
          sse=sse+(x3(pct)-y(pct))**2
c          print*,x(pct),y(pct),x3(pct),sse
        endif
      enddo

      if(pct.ge.2) then 
        rmserr=(sse/((pct*1.)-1))**0.5
      else
        rmserr=skip
      endif
 
      end !subroutine calcrms
**************************LINEAR LEAST SQUARES**************************
       SUBROUTINE LSQ(X,Y,PCT,SLOPE,ICEPT)
 
      INTEGER,PARAMETER :: M=1  ! DEGREE OF POLYNOMIAL MAX OF 10
      INTEGER PCT,LENGTH,N,K,L,KP1,MX2,IP1,NUMBER,I,J
      REAL*4 X(PCT),Y(PCT),SLOPE,ICEPT,SUM
      REAL*4 A(M+1,M+1),B(M+1),C(M+1),P(2*M),TEMP,FACTOR

      LENGTH=PCT
      NUMBER=PCT

      MX2 = M*2
      DO 13 I=1,MX2
         P(I) = 0.0
         DO 13 J=1,LENGTH
   13        P(I)=P(I)+X(J)**I
      N=M+1
      DO 30 I = 1,N
         DO 30 J = 1,N
            K=I+J-2
            IF (K) 29,29,28
   28          A(I,J)=P(K)
               GO TO 30
   29          A(1,1)=NUMBER
   30          CONTINUE
      B(1)=0.0
      DO 21 J = 1, NUMBER
   21    B(1)=B(1)+Y(J)
      DO 22 I = 2,N
         B(I)=0.0
         DO 22 J = 1,NUMBER
   22       B(I)=B(I)+Y(J)*X(J)**(I-1)
      DO 300 K = 1,M
         KP1=K+1
         L=K
         DO 400 I = KP1,N
            IF(ABS(A(I,K))-ABS(A(L,K))) 400,400,401
  401          L=I
  400          CONTINUE
         IF(L-K) 500,500,405
  405          DO 410 J=K,N
                  TEMP=A(K,J)
                  A(K,J)=A(L,J)
  410          A(L,J)=TEMP
            TEMP = B(K)
            B(K)=B(L)
            B(L)=TEMP
  500 DO 300 I = KP1,N
               FACTOR = A(I,K)/A(K,K)
               A(I,K)=0.0
               DO 301 J=KP1,N
  301          A(I,J)=A(I,J) - FACTOR*A(K,J)
  300 B(I)=B(I)-FACTOR*B(K)
      C(N)=B(N)/A(N,N)
      I=M
  710 IP1=I+1
      SUM=0.0
      DO 700 J=IP1,N
  700  SUM = SUM+A(I,J)*C(J)
      C(I)=(B(I)-SUM) / A(I,I)
      I=I-1
      IF (I) 800,800,710
  800 ICEPT=C(1)
      SLOPE=C(2)
c      PRINT*,'Slope: Intercept:'
c      PRINT*,SLOPE,ICEPT
      
      END SUBROUTINE  !LSQ
************************************************************************
****************Load Student's T distribution***************************
      SUBROUTINE LDTTEST(WANT,NTTEST,TTVEC)

      INTEGER NTTEST,DOF(NTTEST),TRSH,TI
      CHARACTER*80 TTFILE
      REAL*4 WANT,CONF(1:6),C1(NTTEST),C2(NTTEST),C3(NTTEST),
     *     C4(NTTEST),C5(NTTEST),C6(NTTEST),TTVEC(NTTEST),COUT
      
**** Load up for the t-test
      TTFILE="ttest2.prn"
      PRINT*,"Getting student's t distribution:"
      PRINT*,TTFILE

   15 FORMAT (I8, 6F8.3)
      OPEN(5,FILE=TTFILE,FORM="FORMATTED")
      READ(UNIT=5,FMT=15) TRSH,CONF(1),CONF(2),CONF(3),CONF(4),CONF(5),
     *     CONF(6)
      DO TI=1,NTTEST
         READ(UNIT=5,FMT=15) DOF(TI),C1(TI),C2(TI),C3(TI),C4(TI),
     *        C5(TI),C6(TI)
      END DO ! TI=1,NTTEST
      CLOSE(5)

      IF (WANT .EQ. CONF(1)) THEN
          TTVEC=C1
          PRINT*, "Confidence equals: "
          COUT=1-(1-CONF(1))*2
          PRINT*,COUT
      ELSEIF (WANT .EQ. CONF(2)) THEN
          TTVEC=C2
          PRINT*, "Confidence equals: "
          COUT=1-(1-CONF(2))*2
          PRINT*,COUT
      ELSEIF (WANT .EQ. CONF(3)) THEN
          TTVEC=C3
          PRINT*, "Confidence equals: "
          COUT=1-(1-CONF(3))*2
          PRINT*,COUT
      ELSEIF (WANT .EQ. CONF(4)) THEN
          TTVEC=C4
          PRINT*, "Confidence equals: "
          COUT=1-(1-CONF(4))*2
          PRINT*,COUT
      ELSEIF (WANT .EQ. CONF(5)) THEN
          TTVEC=C5
          PRINT*, "Confidence equals: "
          COUT=1-(1-CONF(5))*2
          PRINT*,COUT
      ELSEIF (WANT .EQ. CONF(6)) THEN
          TTVEC=C6
          PRINT*, "Confidence equals: "
          COUT=1-(1-CONF(6))*2
          PRINT*,COUT
      ELSE
          TTVEC=C2
          PRINT*, "<Default> Confidence equals: "
          COUT=1-(1-CONF(2))*2
          PRINT*,COUT
      END IF ! (WANT .eq. CONF(1))

      END SUBROUTINE ! LDTTEST
************************************************************************
************************************************************************
      SUBROUTINE STATS(X1,X2,SLOPE,ICEPT,PCT,TVAL,GOODP,SIGSL,SYX)

      INTEGER PCT,L,LL,FLAG(PCT),N,GOODP,LI
      REAL*4 X1(PCT),X2(PCT),X3(PCT),TVAL,PREC,Checkh,Checkl
      REAL*4 YHIGH(PCT),YLOW(PCT),XMM,XSQRT,LHSQ
      REAL*4 SLOPE,ICEPT,RMS,MEANX,MEANY
      REAL*4 SSE,S2YX,SXX,SYY,SYX,SX,SY
      REAL*4 SUM2X,SUMX,SUM2Y,SUMY
      REAL*4 X(PCT),Y(PCT),REX(PCT),REY(PCT),CHISQ,SIGSL

C     HERE X1 = INDEPENDENT X
C          X2 = OBSERVED Y
C          X3 = CALCULATED Y
C          YHIGH = (WANT)% CONFIDENCE HIGH BOUND
C          YLOW = (WANT)% CONFIDENCE LOW BOUND
C          SXX = SAMPLE VARIANCE OF OBSERVED X (X1)
C          SYY = SAMPLE VARIANCE OF OBSERVED Y (X2)
C          RMS = ROOT MEANS SQUARED ERROR
C          TVAL = STUDENTS T-TEST VALUE (DOF=PCT-2@WANT=CONFIDENCE)

      N=PCT
      GOODP=PCT
      PREC=0.1   ! Precision
      X=X1
      Y=X2

      MEANX=0
      MEANY=0
      RMS=0
      SSE=0
      S2YX=0
      SXX=0
      SYY=0
      SUM2X=0
      SUMX=0
      SUM2Y=0
      SUMY=0
      SYX=0
      SX=0
      SY=0
      CHISQ=0
      
      DO L=1,N
         X3(L)=X(L)*SLOPE+ICEPT
         MEANX=MEANX+X(L)
         MEANY=MEANY+Y(L)
         SSE=SSE+(X3(L)-Y(L))**2
         CHISQ=CHISQ+((X3(L)-Y(L))**2)/Y(L)
         SUM2X=SUM2X+X(L)**2
         SUM2Y=SUM2Y+Y(L)**2
         SUMX=SUMX+X(L)
         SUMY=SUMY+Y(L)
      END DO ! L = 1,N

      MEANX=MEANX/N
      MEANY=MEANY/N
      RMS=SQRT(SSE/N)
c      S2YX=SSE/(N-2)
      SXX=(SUM2X-(SUMX**2)/N)/(N-1)
      SYY=(SUM2Y-(SUMY**2)/N)/(N-1)
      S2YX=( (N-1)/(N-2) )*(SYY-SLOPE*SLOPE*SXX)

      SYX=SQRT(S2YX)
      SX=SQRT(SXX)
      SY=SQRT(SYY)
      SIGSL=SYX/(SX*(N-1)**0.5)

***** CONSTRUCT THE CONFIDENCE INTERVAL
*     FOR POINT X0, RANGE +-YX0
*     SYX0=SYX*SQRT( (1/N) + ( (X0-MEANX)^2 / ((N-1)*S2X) ))
*     YX0 = MEANY+-TVAL*SYX0
      XMM=0
      LHSQ=0
      XSQRT=0
      Checkh=0
      Checkl=0
      LI=0
      
      DO LL=1,N
         XMM=MEANY+SLOPE*(X(LL)-MEANX)
         SSY=SSY+(Y(LL)-MEANY)**2
         LHSQ=( (X(LL)-MEANX)**2 ) / ( (N-1) * SXX )
         XSQRT=TVAL*SYX*SQRT( (1/N)+ LHSQ )
         YHIGH(LL)=XMM+XSQRT
         YLOW(LL)=XMM-XSQRT
         Checkh=Y(LL)-YHIGH(LL)
         Checkl=YLOW(LL)-Y(LL)
         IF ( Checkh .gt. PREC )THEN
            FLAG(LL)=1
            GOODP=GOODP-1
         ELSEIF ( Checkl .gt. PREC ) THEN
            FLAG(LL)=1
            GOODP=GOODP-1
         ELSE
            LI=LI+1
            FLAG(LL)=0
            REX(LI)=X(LL)
            REY(LI)=Y(LL)
         END IF
         XMM=0
         LHSQ=0
         XSQRT=0
      END DO



      if (SSY.ne.0) RSQ=(SSY-SSE)/SSY

      if (RSQ.ge.0) then
        CORRC=SQRT(RSQ)
      else
        CORRC=0
      endif
      if (SL .lt. 0 ) CORRC=-1*CORRC

c      PRINT*, "Mean X1, Mean X2"
c      PRINT*, MEANX1,MEANX2
c      PRINT*, "Root Mean Squared Error"
c      PRINT*, RMS
c      PRINT*, "Sample Variance X1, Sample Variance X2"
c      PRINT*, SXX, SYY
c      PRINT*, "Sample variance X1|X2"
c      PRINT*, S2YX
c      PRINT*, "Sample X-Value, Y-Value"
c      PRINT*, X1(test),X2(test)
c      PRINT*, "Y-Predicted, + confidence & - confidence bounds"
c      PRINT*, X3(test),YHIGH(test), YLOW(test)
c      PRINT*, "Total number of points, Number of points in bounds"
c      PRINT*, PCT,GOODP

      END SUBROUTINE !STATS
***********************************************************************
      subroutine basindef(workm,as,basincode,n_basin,n_northatl)
      integer, parameter :: im=360, jm=180, regct=29, LeonMask=-1 !set LeonMask to 1 for unique artic Reigons
      integer i,j,k,kt4(im,jm),workm(im,jm),as(regct,6)
      CHARACTER*80 MASKA,MASKB,MASKC,latFile,lonFile
     *     ,beauLonF,beauLatF,canLonF,canLatF,sibLonF,sibLatF
     *     ,eurLonF,eurLatF,chukLonF,chukLatF,lomLonF,lomLatF
     *     ,lapLonF,lapLatF,karLatF,karLonF,norLonF,norLatF
     *     ,barLonF,barLatF
c     *     ,namebin
      character*8 basincode(regct)
      integer n_arcticoc,n_balticse,n_ginseatl,n_northatl,n_mediterr
     *     ,n_tropiatl,n_southatl,n_southocn,n_northpac,n_tropipac
     *     ,n_southpac,n_reds_p_g,n_indianoc,n_labdorse,n_mostnadw
     *     ,n_nadwaabw,n_mostaabw,n_hudsonba,n_deeparco,n_dpmedite
     *     ,n_baffinba,n_depacind,n_intnopac,n_eurasian,n_lomonosov
     *     ,n_beaufort,n_chukchi
     *     ,n_basin(8)
c     *     ,n_intindoc


!      MASKA="SalMaskA.prn.0" ! File to define ocean basins
!      MASKB="SalMaskB.prn.0" ! File to define ocean basins
!      MASKC="SalMaskC.prn.0" ! File to define ocean basins

c      MASKA="SalMaskA.prn"     ! File to define ocean basins
c      MASKB="SalMaskB0.prn"
c      MASKC="SalMaskC.prn"     ! File to define ocean basins

      MASKA="WOAMaskA2.prn"
      MASKB="WOAMaskB2.prn"
      IF (LeonMask.EQ. 1) THEN !unique arctic
         MASKC="LeonArcMask3.prn"

         PRINT*,"Special Arctic Waterbodies"
         chukLonF="chukLon.txt"
         chukLatF="chukLat.txt"
         beauLonF="beauLon.txt"
         beauLatF="beauLat.txt"
         canLonF="canLon.txt"
         canLatF="canLat.txt"
         sibLonF="sibLon.txt"
         sibLatF="sibLat.txt"
         eurLonF="eurLon.txt"
         eurLatF="eurLat.txt"
         lomLonF="lomLon.txt"
         lomLatF="lomLat.txt"
         lapLonF="lapLon.txt"
         lapLatF="lapLat.txt"
         karLonF="karLon.txt"
         karLatF="karLat.txt"
         norLonF="norLon.txt"
         norLatF="norLat.txt"
         barLonF="barLon.txt"
         barLatF="barLat.txt"

         OPEN(800,FILE=chukLonF,action="write",status="replace")
         OPEN(801,FILE=chukLatF,action="write",status="replace")
         OPEN(810,FILE=beauLonF,action="write",status="replace")
         OPEN(811,FILE=beauLatF,action="write",status="replace")
         OPEN(820,FILE=sibLonF,action="write",status="replace")
         OPEN(821,FILE=sibLatF,action="write",status="replace")
         OPEN(830,FILE=canLonF,action="write",status="replace")
         OPEN(831,FILE=canLatF,action="write",status="replace")
         OPEN(840,FILE=eurLonF,action="write",status="replace")
         OPEN(841,FILE=eurLatF,action="write",status="replace")
         OPEN(850,FILE=lomLonF,action="write",status="replace")
         OPEN(851,FILE=lomLatF,action="write",status="replace")
         OPEN(860,FILE=lapLonF,action="write",status="replace")
         OPEN(861,FILE=lapLatF,action="write",status="replace")
         OPEN(870,FILE=karLonF,action="write",status="replace")
         OPEN(871,FILE=karLatF,action="write",status="replace")
         OPEN(890,FILE=norLonF,action="write",status="replace")
         OPEN(891,FILE=norLatF,action="write",status="replace")
         OPEN(880,FILE=barLonF,action="write",status="replace")
         OPEN(881,FILE=barLatF,action="write",status="replace")
      else
         MASKC="WOAMaskC2.prn"
      ENDIF
      
      
   
         
      
         

      PRINT*, "Opening Ocean Basin Definition Files: "
      PRINT*, trim(MASKA)
      PRINT*, trim(MASKB)
      PRINT*, trim(MASKC)

   41 FORMAT (60I3)  
      
      OPEN(61,FILE=MASKA,FORM="FORMATTED")
      OPEN(62,FILE=MASKB,FORM="FORMATTED")
      OPEN(63,FILE=MASKC,FORM="FORMATTED")
      
      READ(UNIT=61,FMT=41) ((KT4(i,j), J=1,60), I=1,im)
      READ(UNIT=62,FMT=41) ((KT4(i,j), J=61,120), I=1,im)
      READ(UNIT=63,FMT=41) ((KT4(i,j), J=121,jm), I=1,im)

      CLOSE(61)
      CLOSE(62)
      CLOSE(63)

!      namebin="Surface Mask"
!      open(208,FILE="binSalMask",FORM="UNFORMATTED")
!      write(208) namebin,KT4(1:im,1:jm)
!      close(208)
!      print*,"Japan 141 54",kt4(181+141,91+54)
      PRINT*, "Ocean basin mask defined."
      
c For SalMask*.prn files:
c     1x1 Grid with Latitudes in COLUMNS and Longitudes in ROWS
c     Each GridBox coded as follows:
c     60 = Arctic Ocean
c     62 = Baltic Sea
c     61 = Labrador Sea
c     66 = Baffin Bay
c     63 = Hudson Bay, Strait
c   
c     16 = North East Atlantic, GIN: Greenland, Iceland, Norwegian Seas 
c     10 = North Atlantic
c     18 = Mediterranean
c     11 = Tropical Atlantic
c     14 = South Atlantic
c     19 = Southern Ocean south of Atlantic
c     
c     48 = Black Sea
c     58 = Caspian Sea
c
c     26 = NONE! (would be Arctic around Bering Straits)
c     20 = North Pacific
c     22 = Tropical Pacific
c     24 = South Pacific
c     29 = Southern Ocean south of Pacific
c     
c     23 = Indian/Pacific waters South of Aus/NZ
c     32 = Indian/Pacific waters around Indonesia
c     
c     38 = Red Sea waters
c     30 = Tropical Indian Ocean
c     34 = Southern Indian Ocean
c     13 = Atlantic/Indian waters
c     39 = Southern Ocean south of Indian Ocean
c
c     77 = LAND
c      
c     For Leon's Mask...
c     80 = Chukchi Sea
c     81 = Beaufort Sea
c     82 = East Siberian Sea
c     83 = Canadian Basin
c     84 = Eurasian Basin
c     85 = Lomonsov Ridge
c     86 = Laptev sea 
c     87 = Kara Sea
c     88 = Barents Sea       
c     89 = Norweigan Sea
      

      n_arcticoc=1
      n_deeparco=2
      n_balticse=3
      n_ginseatl=4
      n_baffinba=5
      n_labdorse=6
      n_hudsonba=7

      n_northatl=8
      n_tropiatl=9
      n_southatl=10

      n_mediterr=11
      n_reds_p_g=12

      n_northpac=13
      n_tropipac=14
      n_southpac=15

      n_indianoc=16

      n_southocn=17

      n_mostnadw=18
      n_nadwaabw=19
      n_mostaabw=20

      n_intnopac=21

      n_dpmedite=22

      n_depacind=23

      n_aacw_cdw=24

      n_basin(1)=n_mostNADW
      n_basin(2)=n_mostAABW
      n_basin(3)=n_NADWAABW
      n_basin(4)=n_deeparco
      n_basin(5)=n_dpmedite
      n_basin(6)=n_depacind
      n_basin(7)=n_intnopac
      n_basin(8)=n_aacw_cdw
      arcCount=0
      
      PRINT*,"Searchig PRN files"
      if(LeonMask.EQ.1) then
         n_eurasian=25
         n_lomonosov=26
         n_beaufort=27
         n_barents=28
         n_chukchi=29
         
         DO I=1,im              ! Lon 360
            DO J=1,jm           ! Lat 180
               if((j.eq.1).and.(kt4(i,j).ne.77)) print*,"problem! ",i,j
     *              ,kt4(i,j)
               IF (KT4(i,j) .EQ. 80) THEN
                  WorkM(I,J)=n_chukchi ! CHUKCHI
                  arcLon=i
                  arcLat=j
                  WRITE(800,*)arcLon
                  WRITE(801,*)arcLat
c     print*,arcLon,arcLat
                  arcCount =  arcCount+1
               ENDIF
               IF (KT4(i,j) .EQ. 81) THEN !BEAUFORT
                  Print*,"Beau"
                  WorkM(I,J)=n_beaufort
                  arcLon=i
                  arcLat=j
                  WRITE(810,*)arcLon !beauLonF
                  WRITE(811,*)arcLat
                  arcCount =  arcCount+1
               ENDIF
               IF (KT4(i,j) .EQ. 82) THEN !E SIBERIAN
                  Print*,"sib"
                  WorkM(I,J)=n_arcticoc
                  arcLon=i
                  arcLat=j
                  WRITE(820,*)arcLon !sibLonF
                  WRITE(821,*)arcLat
                  arcCount =  arcCount+1
               ENDIF
               IF (KT4(i,j) .EQ. 83) THEN !CANADIAN
                  PRINT*,"Can"
                  WorkM(I,J)=n_arcticoc
                  arcLon=i
                  arcLat=j
                  WRITE(830,*)arcLon !canLonF
                  WRITE(831,*)arcLat
                  arcCount =  arcCount+1
               ENDIF
               IF (KT4(i,j) .EQ. 84) THEN !EURASIAN
                  PRINT*,"eur"
                  WorkM(I,J)=n_eurasian
                  arcLon=i
                  arcLat=j
                  WRITE(840,*)arcLon !euroLonF
                  WRITE(841,*)arcLat
                  arcCount =  arcCount+1
               ENDIF
               IF (KT4(i,j) .EQ. 85) THEN !LAM
                  PRINT*,"Lam"
                  WorkM(I,J)= n_lomonosov
                  arcLon=i
                  arcLat=j
                  WRITE(850,*)arcLon
                  WRITE(851,*)arcLat
                  arcCount =  arcCount+1
               ENDIF
               IF (KT4(i,j) .EQ. 86) THEN !Laptev
                  Print*,"Lap"
                  WorkM(I,J)=n_arcticoc
                  arcLon=i
                  arcLat=j
                  WRITE(860,*)arcLon
                  WRITE(861,*)arcLat
                  arcCount =  arcCount+1
               ENDIF
               IF (KT4(i,j) .EQ. 87) THEN !Kara
                  Print*,"kara"
                  WorkM(I,J)=n_arcticoc
                  arcLon=i
                  arcLat=j
                  WRITE(870,*)arcLon
                  WRITE(871,*)arcLat
                  arcCount =  arcCount+1
               ENDIF
               IF (KT4(i,j) .EQ. 88) THEN !Barents Sea
                  Print*,"Bar"
                  WorkM(I,J)= n_barents
                  arcLon=i
                  arcLat=j
                  WRITE(880,*)arcLon
                  WRITE(881,*)arcLat
                  arcCount =  arcCount+1
               ENDIF
               IF (KT4(i,j) .EQ. 89) THEN !Nowway
                  WorkM(I,J)=n_arcticoc
                  arcLon=i
                  arcLat=j
                  WRITE(890,*)arcLon
                  WRITE(891,*)arcLat
                  arcCount =  arcCount+1
               ENDIF
            
               IF (KT4(i,j) .EQ. 60) WorkM(I,J)=n_arcticoc ! Baltic Sea
               IF (KT4(i,j) .EQ. 62) WorkM(I,J)=n_balticse ! Baltic Sea
               IF (KT4(i,j) .EQ. 62) WorkM(I,J)=n_balticse ! Baltic Sea
               IF (KT4(i,j) .EQ. 66) WorkM(I,J)=n_baffinba ! Baffin Bay
               IF (KT4(i,j) .EQ. 63) WorkM(I,J)=n_hudsonba ! Hudson Bay
               IF (KT4(i,j) .EQ. 61) WorkM(I,J)=n_labdorse ! Labrador Sea
               IF (KT4(i,j) .EQ. 16) WorkM(I,J)=n_ginseatl ! N. N. Atl.
               IF (KT4(i,j) .EQ. 10) WorkM(I,J)=n_northatl ! N. Atl
               IF (KT4(i,j) .EQ. 18) WorkM(I,J)=n_mediterr ! Mediterrannean       
               IF (KT4(i,j) .EQ. 11) WorkM(I,J)=n_tropiatl ! Tr. Atl
               IF (KT4(i,j) .EQ. 14) WorkM(I,J)=n_southatl ! S. Atl         
               IF (KT4(i,j) .EQ. 19) WorkM(I,J)=n_southocn ! Atl. Antarctic
               IF (KT4(i,j) .EQ. 20) WorkM(I,J)=n_northpac ! N. Pac
               IF (KT4(i,j) .EQ. 22) WorkM(I,J)=n_tropipac ! Tr. Pac
               IF (KT4(i,j) .EQ. 24) WorkM(I,J)=n_southpac ! S. Pac
               IF (KT4(i,j) .EQ. 29) WorkM(I,J)=n_southocn ! Pac. Antarctic
               IF (KT4(i,j) .EQ. 38) WorkM(I,J)=n_reds_p_g ! Red Sea/Persian Gulf
               IF (KT4(i,j) .EQ. 32) WorkM(I,J)=n_indianoc ! Indonesian Seas
               IF (KT4(i,j) .EQ. 30) WorkM(I,J)=n_indianoc ! N. Indian
               IF (KT4(i,j) .EQ. 34) WorkM(I,J)=n_indianoc ! S. Indian
               IF (KT4(i,j) .EQ. 13) WorkM(I,J)=n_indianoc ! Indian
               IF (KT4(i,j) .EQ. 39) WorkM(I,J)=n_southocn ! Ind. Antarctic
               IF (KT4(i,j) .EQ. 48) WorkM(I,J)=0 ! Black Sea, not used
               IF (KT4(i,j) .EQ. 58) WorkM(I,J)=0 ! Caspian Sea, not used
               IF (KT4(i,j) .EQ. 77) WorkM(I,J)=0 ! Land
            END DO              ! J=1,jm
         END DO                 ! I=1,im
         CLOSE(800)
         CLOSE(801)
         CLOSE(810)
         CLOSE(811)
         CLOSE(820)
         CLOSE(821)
         CLOSE(830)
         CLOSE(831)
         CLOSE(840)
         CLOSE(841)
         CLOSE(850)
         CLOSE(851)
         CLOSE(860)
         CLOSE(861)
         CLOSE(870)
         CLOSE(871)
         CLOSE(880)
         CLOSE(881)
         CLOSE(890)
         CLOSE(891)
      
         PRINT*, arcCount
         PRINT*,"Arctic Points found...",arcCount
         
      ELSE
         DO I=1,im
            DO J=1,jm
               if((j.eq.1).and.(kt4(i,j).ne.77)) print*,"problem! ",i,j
     *              ,kt4(i,j)
               IF (KT4(i,j) .EQ. 60) WorkM(I,J)=n_arcticoc ! Arctic         
               IF (KT4(i,j) .EQ. 62) WorkM(I,J)=n_balticse ! Baltic Sea
               IF (KT4(i,j) .EQ. 66) WorkM(I,J)=n_baffinba ! Baffin Bay
               IF (KT4(i,j) .EQ. 63) WorkM(I,J)=n_hudsonba ! Hudson Bay
               IF (KT4(i,j) .EQ. 61) WorkM(I,J)=n_labdorse ! Labrador Sea
               IF (KT4(i,j) .EQ. 16) WorkM(I,J)=n_ginseatl ! N. N. Atl.
               IF (KT4(i,j) .EQ. 10) WorkM(I,J)=n_northatl ! N. Atl
               IF (KT4(i,j) .EQ. 18) WorkM(I,J)=n_mediterr ! Mediterrannean       
               IF (KT4(i,j) .EQ. 11) WorkM(I,J)=n_tropiatl ! Tr. Atl
               IF (KT4(i,j) .EQ. 14) WorkM(I,J)=n_southatl ! S. Atl         
               IF (KT4(i,j) .EQ. 19) WorkM(I,J)=n_southocn ! Atl. Antarctic
               IF (KT4(i,j) .EQ. 20) WorkM(I,J)=n_northpac ! N. Pac
               IF (KT4(i,j) .EQ. 22) WorkM(I,J)=n_tropipac ! Tr. Pac
               IF (KT4(i,j) .EQ. 24) WorkM(I,J)=n_southpac ! S. Pac
               IF (KT4(i,j) .EQ. 29) WorkM(I,J)=n_southocn ! Pac. Antarctic
               IF (KT4(i,j) .EQ. 38) WorkM(I,J)=n_reds_p_g ! Red Sea/Persian Gulf
               IF (KT4(i,j) .EQ. 32) WorkM(I,J)=n_indianoc ! Indonesian Seas
               IF (KT4(i,j) .EQ. 30) WorkM(I,J)=n_indianoc ! N. Indian
               IF (KT4(i,j) .EQ. 34) WorkM(I,J)=n_indianoc ! S. Indian
               IF (KT4(i,j) .EQ. 13) WorkM(I,J)=n_indianoc ! Indian
               IF (KT4(i,j) .EQ. 39) WorkM(I,J)=n_southocn ! Ind. Antarctic
               IF (KT4(i,j) .EQ. 48) WorkM(I,J)=0 ! Black Sea, not used
               IF (KT4(i,j) .EQ. 58) WorkM(I,J)=0 ! Caspian Sea, not used
               IF (KT4(i,j) .EQ. 77) WorkM(I,J)=0 ! Land
            END DO              ! J=1,jm
         END DO                 ! I=1,im
      END IF   

!     print*,"Japan",workm(181+141,91+54)
      PRINT*, "Ocean Basins Defined"



*******Define Compatible regions 

      do k=1,regct
        as(k,2)=k
        do i=3,6
          as(k,i)=0
        enddo
      enddo

      basincode(n_arcticoc)="ArcticOc" ! Arctic Ocean
      as(n_arcticoc,3)=n_ginseatl ! with GIN

      basincode(n_balticse)="BalticSe" ! Baltic Sea

      basincode(n_ginseatl)="GINSeAtl" ! Greenland/Iceland/Norwegian Seas
      as(n_ginseatl,3)=n_arcticoc ! with Arctic
      as(n_ginseatl,4)=n_northatl ! with N. Atl
      
      basincode(n_northatl)="NorthAtl" ! North Atlantic
      as(n_northatl,3)=n_ginseatl ! with GIN
      as(n_northatl,4)=n_tropiatl ! with Tr. Atl
      as(n_northatl,5)=n_labdorse ! with Labrador Sea
      
      basincode(n_mediterr)="Mediterr" ! Mediterrannean Sea

      basincode(n_tropiatl)="TropiAtl" ! Tropical Atlantic
      as(n_tropiatl,3)=n_northatl ! with N. Atl
      as(n_tropiatl,4)=n_southatl ! with S. Atl

      basincode(n_southatl)="SouthAtl" ! South Atlantic
      as(n_southatl,3)=n_tropiatl ! with Tr. Atl
      as(n_southatl,4)=n_southocn ! with Antarctic
      as(n_southatl,5)=n_indianoc ! with Indian 

      basincode(n_southocn)="SouthOcn" ! Southern Ocean
      as(n_southocn,3)=n_southatl ! with S. Atl
      as(n_southocn,4)=n_southpac ! with S. Pac
      as(n_southocn,5)=n_indianoc ! with  Indian Ocean

      basincode(n_northpac)="NorthPac" ! North Pacific
      as(n_northpac,3)=n_tropipac ! with  Tr. Pac

      basincode(n_tropipac)="TropiPac" ! Tropical Pacific
      as(n_tropipac,3)=n_northpac ! with N. Pac
      as(n_tropipac,4)=n_indianoc ! with Ind
      as(n_tropipac,5)=n_southpac ! with S. Pac
      
      basincode(n_southpac)="SouthPac" ! South Pacific
      as(n_southpac,3)=n_tropipac ! with Tr. Pac
      as(n_southpac,4)=n_southocn ! with Antarctic
      as(n_southpac,5)=n_indianoc ! with Indian Ocean

      basincode(n_reds_p_g)="RedS,P-G" ! Red Sea / Persian Gulf

      basincode(n_indianoc)="IndianOc" ! Indian Ocean
      as(n_indianoc,3)=n_tropipac ! with Tr. Pac
      as(n_indianoc,4)=n_southatl ! with S. Atl.
      as(n_indianoc,5)=n_southocn ! with Antarctic

      basincode(n_labdorse)="LabdorSe" ! Labrador Sea
      as(n_labdorse,3)=n_northatl ! with North Atlantic

      basincode(n_mostnadw)="mostNADW" ! North Atlantic Deep Water
      as(n_mostnadw,3)=n_ginseatl ! with GIN seas
      as(n_mostnadw,4)=n_northatl ! with north atl
c      as(n_mostnadw,5)=n_tropiatl ! with trop atl
c      as(n_mostnadw,6)=n_deeparco ! with deep Arctic 

      basincode(n_nadwaabw)="NADWAABW" ! mix NADW & AABW
c      as(n_nadwaabw,3)=n_mostaabw ! with most AABW

      basincode(n_mostaabw)="mostAABW" ! Antarctic Bottom Water
c      as(n_mostaabw,3)=n_nadwaabw ! with NADWAABW

      basincode(n_hudsonba)="HudsonBa"

      basincode(n_deeparco)="deepArcO" ! Deep Arctic Ocean
      as(n_deeparco,3)=n_mostnadw ! with NADW
      as(n_deeparco,4)=n_arcticoc ! with Shallow Arctic

      basincode(n_dpmedite)="DpMedite" ! Deep Mediterr
      as(n_dpmedite,3)=n_mediterr ! Shallow Mediterr

      basincode(n_baffinba)="BaffinBa"! Baffin Bay
      basincode(n_depacind)="DePacInd" ! Pacific/Indian Ocean deeper than 2500

      basincode(n_intnopac)="IntNoPac" ! Intermediate North Pacific
      as(n_intnopac,3)=n_northpac

      basincode(n_aacw_cdw)="AACW_CDW" ! Antarctic Commonwater/ Common DW

!     Leon added artic
      IF(LeonMask.EQ.1) THEN
         basincode(n_eurasian)= "Eurasian" ! 
         basincode(n_lomonosov)="Lomonoso" ! 
         basincode(n_barents)=  "BarentsS" ! 
         basincode(n_beaufort)= "Beaufort" !  
         basincode(n_chukchi)=  "ChukchiS" !
      endif
      do k=1,regct
        as(k,1)=1
        do i=3,6
          if (as(k,i).gt.0) as(k,1)=as(k,1)+1
        enddo
      enddo

c      do k=1,regct
c        print*,trim(basincode(k)),as(k,1),as(k,2),as(k,3),as(k,4),as(k,5
c     *       ),as(k,6)
c        enddo
      PRINT*,"Compatible regions assigned"
 
      end subroutine            ! basindef

****************************Subroutine XQFIT***************************
! CALCULATE GOODNESS OF FIT
      subroutine XQFIT(X1,X2,SLOPE,ICEPT,PCT,QFIT)
      IMPLICIT NONE
      INTEGER I,PCT,A
      REAL X1(PCT),X2(PCT),SLOPE,ICEPT,Y,CHI2,DIFFSQ,HCHI2,HA,SIGMA2
     *     ,gammq,qfit
      REAL,PARAMETER :: SIGMA=0.2
      INTEGER,PARAMETER :: M=2
	
      A=PCT-M
      DIFFSQ=0
      CHI2=0
      SIGMA2=SIGMA**2
 22   format(f8.2, f8.3, f8.3)
!     ***Calculate predicted Y-values and corresponding chisq
c      PRINT*,"    X1      X2,  YPRED"
      DO I=1,PCT
         Y=X1(I)*SLOPE+ICEPT
c         PRINT22,X1(I),X2(I),Y
	 DIFFSQ=(X2(I)-Y)**2
	 CHI2=CHI2+DIFFSQ/SIGMA2
      END DO  ! I=1,PCT
      HCHI2=0.5*CHI2
      HA=0.5*A
      QFIT = GAMMQ(hA,HCHI2)
c      print*,"ha",hA,hchi2,qfit
c      PRINT*,"Chi-Squared: "
c      PRINT*, CHI2
c      PRINT*,"DOF"
c      PRINT*,A      
c      PRINT*,"QFit "
c      PRINT*,QFIT

      END subroutine ! XQFIT
	
************************************************************************
************************Function GAMMQ*********************************
      FUNCTION gammq(a,x)
C USES gcf,gser
C Returns the incomplete gamma function Q(a,x) = 1 - P(a,x).
      IMPLICIT NONE
      REAL a,gammq,x
      REAL gammcf,gamser,gln
!      print*,a,x
      if(x.lt.0..or.a.le.0.) PRINT*, "bad arguments in gammq" ! pause
      if(x.lt.a+1.) then    ! Use the series representation
         call gser(gamser,a,x,gln)
         gammq=1.-gamser   ! and take its complement.
      else                 ! Use the continued fraction representation.
         call gcf(gammcf,a,x,gln)
         gammq=gammcf
      endif                ! (x.lt.a+1.)then
!      print*,gammq
      return
      END ! Function gammq(a,x)
************************************************************************
***************************Subroutine GSER******************************
      SUBROUTINE gser(gamser,a,x,gln)
! USES gammln
! Returns the incomplete gamma function P(a, x) evaluated by its series 
! representation as gamser. Also returns ln gamma(a) as gln.
      IMPLICIT NONE
      INTEGER,PARAMETER :: ITMAX=100
      REAL, PARAMETER :: EPS=3.e-7
      REAL a,gamser,gln,x
      INTEGER n
      REAL ap,del,sum,gammln
!      print*,"gser",a,x
      gln=gammln(a)
      if(x.le.0.)then
      if(x.lt.0.) print*, "x < 0 in gser" !pause
      gamser=0.
      return
      endif
      ap=a
      sum=1./a
      del=sum
      do n=1,ITMAX
         ap=ap+1.
         del=del*x/ap
         sum=sum+del
         if(abs(del).lt.abs(sum)*EPS)goto 1
      enddo ! n=1,ITMAX
!      print*, "a too large, ITMAX too small in gser" !pause
    1 gamser=sum*exp(-x+a*log(x)-gln)
!      print*,"gser",gln,gamser
      return
      END ! subroutine GSER

************************************************************************
***************************Subroutine GCF*******************************
      SUBROUTINE gcf(gammcf,a,x,gln)
C USES gammln
! Returns the incomplete gamma function Q(a,x) evaluated by its 
! continued fraction representation as gammcf. Also returns ln gamma(a)
! as gln.
! Parameters: ITMAX is the maximum allowed number of iterations;
! EPS is the relative accuracy;
! FPMIN:  number near the smallest representable floating-point number.
      IMPLICIT NONE
      INTEGER ITMAX
      REAL a,gammcf,gln,x,EPS,FPMIN
      PARAMETER (ITMAX=100,EPS=3.e-7,FPMIN=1.e-30)
      INTEGER i
      REAL an,b,c,d,del,h,gammln
!      print*,"gcf",a,x
      gln=gammln(a)
      b=x+1.-a            ! Set up for evaluating continued fraction by 
      c=1./FPMIN          ! modified Lentz's method (5.2) with b0 = 0.
      d=1./b
      h=d
      do i=1,ITMAX     ! Iterate to convergence.
         an=-i*(i-a)
         b=b+2.
         d=an*d+b
         if(abs(d).lt.FPMIN)d=FPMIN
         c=b+an/c
         if(abs(c).lt.FPMIN)c=FPMIN
         d=1./d
         del=d*c
         h=h*del
         if(abs(del-1.).lt.EPS)goto 1
      enddo !i=1,ITMAX 
      print*, "a too large, ITMAX too small in gcf" ! pause
    1 gammcf=exp(-x+a*log(x)-gln)*h ! Put factors in front.
!      print*,"gcf",gln,gammcf
      return
      END ! subroutine GCF

************************************************************************
************************************************************************

      FUNCTION GAMMLN(XX)
! Returns the value ln(gamma(xx)) for xx>0.  Full accuracy is obtained 
! for XX>1.  For 0 < XX < 1, the reflection formula can be used first.
      IMPLICIT NONE
      REAL GAMMLN,XX ! Returns the value
      INTEGER J
      DOUBLE PRECISION SER, STP, TMP, X, Y, COF(6)
! Internal arithmetic will be done in double precision, a nicety that
! you can omit if five figure accuracy is good enough.
      SAVE COF, STP
      DATA COF, STP/76.18009172947146d0,-86.50532032941677D0,
     *  24.01409824083091D0,-1.231739572450155D0,0.1208650973866179D-2,
     *  -0.5395239384953D-5,2.5066282746310005D0/
      
      X=XX
      Y=X
      TMP=X+5.5D0
      TMP=(X+0.5D0)*LOG(TMP)-TMP
      SER=1.000000000190015D0
      DO J=1,6
         Y=Y+1.D0
         SER=SER+COF(J)/Y
      ENDDO

      GAMMLN=TMP+LOG(STP*SER/X)
      RETURN
      END FUNCTION GAMMLN


