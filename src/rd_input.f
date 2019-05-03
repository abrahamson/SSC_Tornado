c -------------------------------------------------------------------------

      subroutine RdInput ( nInten,  testInten,  nGM_model, nAttenType, attenType, nProb, iPer,
     5     jcalc, scalc, sigFix, sssCalc, wt_GM1 )

      implicit none
      include 'tornado.h'

      real testInten(MAX_INTEN)
      real minlat,maxlat,minlong,maxlong,maxdist
      integer nInten, ntotal, attentype(MAX_FLT)
      integer nwr, iProb, iPer, dirFlag, j,jj
      real specT, sigTrunc, varADD
      integer jCalc(MAX_ATTENTYPE,MAX_ATTEN),  sssCalc(MAX_ATTENTYPE,MAX_ATTEN), scalc(MAX_ATTENTYPE,MAX_ATTEN)
      real sigFix(MAX_ATTENTYPE,MAX_ATTEN)

      character*80 filein, title
      integer nfiles

      integer nProb, nattentype, nGM_Model(MAX_ATTENTYPE)
      real  checkwt, c1, c2, wtgm(4,MAX_ATTEN), wt_gm1(4,MAX_ATTEN)

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
      flush (6)
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
      write (*,'( 3i5)') nProb, nattentype

      do iprob=1,nProb

C       Read period, maxeps dir flag and gm intensities
        read (20,*) specT, sigtrunc, dirflag 
        read (20,*) nInten
        write (*,'( i5)') nInten
        pause 'nInten'
c        read (20,*) nInten, (testInten(j), j=1,nInten)
        call CheckDim ( nInten, MAX_INTEN, 'MAX_INTEN' )

C       Read in the suite of attenution models and wts for each attentype
        do j=1,nattentype
          write (*,'( 2i5)') iProb, j
          pause 'test j'
          checkwt = 0.0
          read (20,*) nGM_model(j)
c         Check for Max number of attenuation model
          call checkDim ( nGM_model(j), MAX_ATTEN, 'MAX_ATTEN' )

          do jj=1,nGM_model(j)
            read (20,*) jcalc(j,jj), c1, c2, wtgm(j,jj), Varadd
            if ( jCalc(j,jj) .lt. 0 ) then
               backspace (20)
               read (20,*) jcalc(j,jj), c1, c2, wtgm(j,jj), Varadd, sCalc(j,jj), sigfix(j,jj), sssCalc(j,jj)
            endif
          enddo
        enddo

c       keep the weights for the selected problem only
        if ( iper .eq. iProb) then 
          do j=1,nattentype
            do jj=1,nGM_model(j)
              Wt_GM1(j,jj) = Wtgm(j,jj)
            enddo
          enddo
        endif

      enddo

      close (20)
      pause 'test 4'
       

       return
       end
