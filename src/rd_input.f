c -------------------------------------------------------------------------

      subroutine RdInput ( nInten,  
     1           testInten, nGM_model, cumWt_GM, nAttenType, attenType, 
     4           nProb, iPer, SpecT1, version )

      implicit none
      include 'tornado.h'

      integer nInten, ntotal, attentype(MAX_FLT), nfiles, ix(MAX_FILES),
     1        jCalc(MAX_ATTENTYPE,MAX_ATTEN), 
     2        nProb, nattentype, nGM_Model(MAX_ATTENTYPE), 
     3        jj, j
      integer iprob, iPer,  nwr
      real testInten(MAX_INTEN), specT, dirflag,
     2     checkwt, c1, c2, wtgm(4,MAX_ATTEN), 
     3     sigtrunc, Varadd,
     4    cumWt_GM(MAX_ATTENTYPE,MAX_ATTEN)
      real  SpecT1, version, x
      character*80 filein, title , dummy     

c     Set Data file units
      nwr = 11

      ntotal = 0

c     Program no longer allowed to read from multiple files.
      nFiles = 1

c     Open PSHA Run Input File
      read (31,'( a80)') filein
      write (*,'( a80)') filein
      open (20,file=filein,status='old')

c     Open Input PSHA Source/Fault file
      read (20,'( a80)') filein
      write (*,'( a80)') filein
      
      open (10,file=filein,status='old')
      
C     Check for version compatibility with hazard code
      read (20,*,err=1000) version
         if (version .ne. 45.1 .and. version .ne. 45.2 .and. version .ne. 45.3) then
         write (*,*) 'Incompatible version of Haz45, use Haz45.2 or haz45.3'
         stop 99
        endif

c     Read in parameters for background grid.
      read (20,*,err=1001) x 
      read (20,*,err=1002) x

c     Input Title (not used) 
      read(20,'( a80)') title
      write (*,'( a80)') title

c     Number of Spectral Periods and Number of attenuation relations types
      read(20,*,err=1003) nProb, nattentype

      write (*,'( 2i5)') nProb, nattentype
      call checkDim ( nProb, MAX_PROB, 'MAX_PROB' )
      call checkDim ( nattentype, MAX_ATTENTYPE, 'MAX_ATTENTYPE' )
      
      do iprob=1,nProb

C       Read period, maxeps dir flag and gm intensities
        read (20,*,err=1004) specT, sigtrunc, dirflag 
        write (*,'( 2f10.4)') specT, sigtrunc
c        read (20,*,err=1005) nInten
        read (20,*,err=1005) nInten
        call CheckDim ( nInten, MAX_INTEN, 'MAX_INTEN' ) 
        backspace (20)
        read (20,*,err=1005) nInten, (testInten(j), j=1,nInten)

C       Read in the suite of attenution models and wts for each attentype
        do j=1,nattentype
          checkwt = 0.0
          read (20,*,err=1006) nGM_model(j)
c         Check for Max number of attenuation model
          call checkDim ( nGM_model(j), MAX_ATTEN, 'MAX_ATTEN' )

          do jj=1,nGM_model(j)
            read (20,*,err=1007) jcalc(j,jj), c1, c2, wtgm(j,jj), Varadd
            if ( abs(jcalc(j,jj)) .ge. 15000 .and. abs(jcalc(j,jj)) .lt. 15200 ) then
               read(20,'( a80)') dummy
            endif
          enddo
        enddo

c       keep the weights for the selected problem only
        if ( iper .eq. iProb) then 
         specT1 = specT
C         Set up cum wts for attenuation models. 
          
          do j=1,nattentype
            do jj=1,nGM_model(j)
              if (jj .eq. 1) then
                cumWt_GM(j,jj) = Wtgm(j,jj)
              else
                cumWt_GM(j,jj) = Wtgm(j,jj) + cumWt_GM(j,jj-1)
              endif
            enddo
c           check that weights sum to unity
            if ( abs(cumWt_GM(j,nGM_model(j))-1.) .gt. .001 ) then
              write (*,'( 2x,''weights do nut sum to unity for GM models'')')
              write (*,'( 2x,''attenType='',i5)') j
              write (*,'( 10f10.4)')  (Wtgm(j,jj),jj=1, nGM_model(j))
              write (*,'( f12.8)')cumWt_GM(j,nGM_model(j)) 
              stop 99
            else
              cumWt_GM(j,nGM_model(j)) = 1.
            endif           
          enddo
        endif

      enddo

      close (20)

      ix(1) = 0


       return
 1000  write (*,'( 2x,''err reading run file, 1000'')')
       backspace (20)
       read (20,'( a80)') dummy
       write (*,'( a80)') dummy
       stop 99     
 1001  write (*,'( 2x,''err reading run file, 1001'')')
       backspace (20)
       read (20,'( a80)') dummy
       write (*,'( a80)') dummy
       stop 99     
 1002  write (*,'( 2x,''err reading run file 1002'')')
       backspace (20)
       read (20,'( a80)') dummy
       write (*,'( a80)') dummy
       stop 99     
 1003  write (*,'( 2x,''err reading run file 1003'')')
       backspace (20)
       read (20,'( a80)') dummy
       write (*,'( a80)') dummy
       stop 99     
 1004  write (*,'( 2x,''err reading run file 1004'')')
       backspace (20)
       read (20,'( a80)') dummy
       write (*,'( a80)') dummy
       stop 99     
 1005  write (*,'( 2x,''err reading nInten'')')
       backspace (20)
       read (20,'( a80)') dummy
       write (*,'( a80)') dummy
       stop 99     
 1006  write (*,'( 2x,''error in reading nGM_model'')')
       backspace (20)
       read (20,'( a80)') dummy
       write (*,'( a80)') dummy
       stop 99 
 1007  write (*,'( 2x,''err reading run file 1007'',2i5)') j, jj
       backspace (20)
       read (20,'( a80)') dummy
       write (*,'( a80)') dummy       
       stop 99     
       
       end
