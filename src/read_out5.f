c -------------------------------------------------------------------

      subroutine read_out5 (haz, iPer, iInten, nFlt, nProb, iOpen)

      implicit none
      include 'tornado.h'

      real haz(MAX_FLT, MAX_WIDTH, MAXPARAM,MAX_FTYPE)
      real temp(MAX_INTEN), version
      integer  nInten, nAtten, iProb, nWidth(MAX_FLT), nParamVar(MAX_FLT,MAX_WIDTH)
      integer iPer, nwr, iFlt, jProb, jFlt, i, iWidth, iFtype, j
      integer nFtype(MAX_FLT), iInten, iWidth1
      integer iUnit_Out5, nFlt, nPRob, iOpen
      character*80 file1, dummy
      iUnit_Out5 = 32
            
c     Open output file (only for the first intensity)
      if (iOpen .eq. 0) then
        read (31,'( a80)') file1
        write (*,*) 'Opening out5 file from the hazard runs.'
        write (*,*) file1
        open (iUnit_Out5,file=file1,status='old')
        iOpen = 1
      endif
      
c     Rewind for this intensity
      rewind (iUnit_Out5)
      
C     Check for version compatibility with hazard code
      read (iUnit_Out5,*) version
      if (version .ne. 45.2 .and. version .ne. 45.3) then
          write (*,*) 'Incompatible version of Haz45 out1 file, use Haz45.2 or haz45.3'
          backspace (iUnit_Out5)
          read (iUnit_Out5,'( a80)') dummy
          write (*,'( 2x,a80)') dummy
          stop 99
      endif

c     Loop over faults
      do iFLt=1,nFlt
        
c       Initial read to check the dimension
        read (iUnit_Out5,*,err=2001) jFlt, iProb, nAtten, nWidth(iFlt), nFtype(iFlt), 
     1         (nParamVar(iFlt,i),i=1,nWidth(iFlt)), nInten
        backspace (iUnit_Out5)

C       Check for Array Dimension of haz Array.
        if (nWidth(iFlt) .gt. MAX_WIDTH) then
           write (*,*) 'MAX_WIDTH needs to be increased to ', nWidth(iFlt)
           stop 99
        endif
        if (nFtype(iFlt) .gt. MAX_FTYPE) then
           write (*,*) 'MAX_FTYPE needs to be increased to ', nFtype(iFlt)
           stop 99
        endif
        do iWidth1=1,nWidth(iFlt)
          if (nParamVar(iFlt,iWidth1) .gt. MAXPARAM) then
             write (*,*) 'MAXPARAM needs to be increased to ', nParamVar(iFlt,iWidth)
             stop 99
          endif
        enddo

       do iWidth=1,nWidth(iFlt)
        
        do jProb=1,nProb 

c        read header (repeated for each problem)

       
         read (iUnit_Out5,*,err=2001) jFlt, iProb, nAtten, nWidth(iFlt), nFtype(iFlt), 
     1         (nParamVar(iFlt,i),i=1,nWidth(iFlt)), nInten

         do iFtype=1,nFtype(iFlt)
          do i=1,nParamVar(iFlt,iWidth)
 	    read (iUnit_Out5,*) (temp(j),j=1,nInten)

 	       
c           Only keep if it is the desired spectral period
	    if ( iProb .eq. iPer) then
	      haz(iFlt,iWidth,i,iFtype) = temp(iInten)

            endif
          enddo
         enddo
        enddo
       enddo
      enddo

      return
      
 2000 write (*,'( 2x,''Error in version number'')')
      goto 3001
 2001 write (*,'( 2x,''Error in out5, header for iflt, jProb:'',2i5)') iFlt, jProb
      goto 3001
 2002 write (*,'( 2x,''Error in out5, haz values: '')')
      write (*,'( 2x,''iflt, jProb, iWidth, iparam:'',5i5 )') 
     1       iFlt, jProb, iWidth, iFtype, i
      goto 3000

 3000 backspace (iUnit_Out5)
      backspace (iUnit_Out5)
      read (iUnit_Out5,'( a80)') dummy
      write (*,'( a80)') dummy
      read (iUnit_Out5,'( a80)') dummy
      write (*,'( a80)') dummy
      stop 99

 3001 backspace (iUnit_Out5)
      read (iUnit_Out5,'( a80)') dummy
      write (*,'( a80)') dummy
      stop 99

      end

c -------------------------------------------------------------------------

      subroutine read_out3 (haz_mean, haz_total, iPer )

      implicit none
      include 'tornado.h'

      real haz_mean(MAX_FLT,MAX_INTEN), haz_total(MAX_INTEN)
      real version
      real poisson(MAX_INTEN), M_bar(MAX_INTEN), D_bar(MAX_INTEN), eps_bar(MAX_INTEN)
      integer  nInten, nFlt
      integer iPer, iFlt, jProb, i
      integer iUnit_Out3, nInten1, nPer1, iPEr1, k
      real period, pSeg, al_wt, Rmin, x1, x2, x3
      character*80 file1, dummy, fname
      iUnit_Out3 = 13

c     Open output file and read header
      read (31,'( a80)') file1
      write (*,*) 'Opening out3 file from the hazard runs.'
      write (*,*) file1
      open (iUnit_Out3,file=file1,status='old')

C     Check for version compatibility with hazard code
      read (iUnit_Out3,*) version
      if (version .ne. 45.2 .and. version .ne. 45.3) then
        write (*,*) 'Incompatible version of Haz45 out1 file, use Haz45.2 or haz45.3'
        backspace (iUnit_Out3)
        read (iUnit_Out3,'( a80)') dummy
        write (*,'( 2x,a80)') dummy
        stop 99
      endif
      
c     Read out3 header
      read (iUnit_Out3,*) Nflt
      do i=1,5
       read (iUnit_Out3,'( a80)') dummy
      enddo   
      read (iUnit_Out3,*) nPer1

c     Read out3 file for this period 
c     Only read up to this problem         
      do jProb=1,iPer 
       
c      Read header for this period
c       read (iUnit_Out3,'( a80)') dummy
c       read (iUnit_Out3,'( a80)') dummy
c       read (iUnit_Out3,'( a80)') dummy
c       write (*,'( a80)') dummy
       read (iUnit_Out3,'( 14x,i6,f12.4)') iPer1, period
       read (iUnit_Out3,*) nInten
       read (iUnit_Out3,'( a80)') dummy

c      Read hazard by fault
       do iFlt=1,nFlt
         read (iUnit_Out3,'( a40,2f6.2,f8.2,1x,50(e12.4))',err=2002) fname, pSeg, al_wt, Rmin,
     1       (haz_mean(iFlt,k),k=1,nInten)
c         write (*,'( i5,2x,a40,20e12.4)') iFlt, fname,(haz_mean(iFlt,k),k=1,nInten) 
       enddo    

c      Read the total hazard
       Read (iUnit_Out3,'( 40x,2f6.2,f8.2,1x,50(e12.4))',err=2003) x1, x2, x3,
     1       (haz_total(k),k=1,nInten)
       Read (iUnit_Out3,'( 40x,2f6.2,f8.2,1x,50(e12.4))',err=2004) x1, x2, x3,
     1       (Poisson(k),k=1,nInten)
       Read (iUnit_Out3,'( 40x,2f6.2,f8.2,1x,50(e12.4))',err=2005) x1, x2, x3,
     1       (M_Bar(k),k=1,nInten)
       Read (iUnit_Out3,'( 40x,2f6.2,f8.2,1x,50(e12.4))',err=2006) x1, x2, x3,
     1       (D_bar(k),k=1,nInten)
       Read (iUnit_Out3,'( 40x,2f6.2,f8.2,1x,50(e12.4))',err=2007) x1, x2, x3,
     1       (eps_bar(k),k=1,nInten)
     
       read (iUnit_Out3,'( a80)') dummy
       read (iUnit_Out3,'( a80)') dummy
       read (iUnit_Out3,'( a80)') dummy
      enddo

      return
      
 2000 write (*,'( 2x,''Error in version number'')')
      goto 3000
 2002 write (*,'( 2x,''Error in out3, haz values: '')')
      write (*,'( 2x,''jProb:'',5i5 )') jProb
      goto 3000
 2003 write (*,'( 2x,''Error in out3, haz total values: '')')
      write (*,'( 2x,''jProb:'',5i5 )') jProb
      goto 3000
 2004 write (*,'( 2x,''Error in out3, poisson values: '')')
      write (*,'( 2x,''jProb:'',5i5 )') jProb
      goto 3000
 2005 write (*,'( 2x,''Error in out3, M_bar values: '')')
      write (*,'( 2x,''jProb:'',5i5 )') jProb
      goto 3000
 2006 write (*,'( 2x,''Error in out3, D_bar values: '')')
      write (*,'( 2x,''jProb:'',5i5 )') jProb
      goto 3000
 2007 write (*,'( 2x,''Error in out3, eps_bar values: '')')
      write (*,'( 2x,''jProb:'',5i5 )') jProb
      goto 3000

 3000 backspace (iUnit_Out3)
      backspace (iUnit_Out3)
      read (iUnit_Out3,'( a80)') dummy
      write (*,'( a80)') dummy
      read (iUnit_Out3,'( a80)') dummy
      write (*,'( a80)') dummy
      stop 99

      end

c -------------------------------------------------------------------------
