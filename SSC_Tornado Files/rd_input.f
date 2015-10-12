c -------------------------------------------------------------------------

      subroutine RdInput ( nInten,  testInten, iPer )

      implicit none
      include 'tornado.h'

      real testInten(MAX_INTEN)
      real minlat,maxlat,minlong,maxlong,maxdist
      integer nInten, ntotal, attentype(MAX_FLT)
      integer nwr, iProb, iPer, dirFlag, j,jj
      real specT, sigTrunc, varADD
      integer jCalc

      character*80 filein, title
      integer nfiles

      integer nProb, nattentype, nGM_Model
      real  checkwt, c1, c2, wtgm

c      pause 'inside read input'


c     Set Data file units
      nwr = 11

      ntotal = 0

c     Read in the number of data files.
c      read (5,*) nfiles
c     Program no longer allowed to read from multiple files.
      nFiles = 1

c     Loop over the number of files.
c      do 111 iii=1,nfiles

c     Open PSHA Run Input File
      read (31,'( a80)') filein
      write (*,'( a80)') filein
      open (20,file=filein,status='old')

c     Open Input PSHA Source/Fault file
      read (20,'( a80)') filein
      open (10,file=filein,status='old')

c     Read in parameters for background grid.
      read (20,*) minlat,maxlat,minlong,maxlong

Cnjg  Added back in read of single maxdist
      read (20,*) maxdist

c     Input Title (not used) 
      read(20,'( a80)') title

c     Number of Spectral Periods and Number of attenuation relations types
      read(20,*) nProb, nattentype

      do iprob=1,nProb

C       Read period, maxeps dir flag and gm intensities
        read (20,*) specT, sigtrunc, dirflag 
        read (20,*) nInten
        call CheckDim ( nInten, MAX_INTEN, 'MAX_INTEN' )
        backspace (20)
        read (20,*) nInten, (testInten(j), j=1,nInten)
        if ( iPer .eq. iProb ) return

C       Read in the suite of attenution models and wts for each attentype
        do j=1,nattentype
          checkwt = 0.0
          read (20,*) nGM_model

          do jj=1,nGM_model
            read (20,*) jcalc, c1, c2, wtgm, Varadd
          enddo
        enddo

      enddo

      close (20)
       

       return
       end
