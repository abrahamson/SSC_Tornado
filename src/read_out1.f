c -------------------------------------------------------------------

      subroutine read_logichaz ( nFlt, haz, nAtten, nFtypeTotal, iPer, nProb)

      implicit none
      include 'tornado.h'

      real haz(MAX_INTEN, MAX_ATTEN, MAX_FLT, MAX_WIDTH, MAXPARAM,MAX_FTYPE)
      real temp(MAX_INTEN)
      integer  nInten, nFlt, nAtten(MAX_FLT), iProb, nWidth(MAX_FLT), nParamVar(MAX_FLT,MAX_WIDTH)
      integer iPer, nProb, nwr, iFlt, jProb, jFlt, i, iWidth, iAtten, iFtype, j
      integer nFtypeTotal(MAX_FLT)
      character*80 file1
      nwr = 12

                     
c     Open output 1 file
      read (31,'( a80)') file1

      write (*,*) 'Opening output file from the hazard runs.'
      write (*,'( a80)' ) file1
      open (nwr,file=file1,status='old')

            
      do iFlt=1,nFlt      

        do jProb=1,nProb

          read (nwr,*) jFlt, iProb, nAtten(iFlt), nWidth(iFlt), nFtypeTotal(iFlt), 
     1         (nParamVar(iFlt,i),i=1,nWidth(iFlt)), nInten
          write (*,'( 20i5)') iFlt, nFlt, jProb, nProb,  nAtten(iFlt), nWidth(iFlt), nFtypeTotal(iFlt), 
     1         (nParamVar(iFlt,i),i=1,nWidth(iFlt)), nInten
 

C     Check for Array Dimension of haz Array.
         if (nFlt .gt. MAX_FLT) then
           write (*,*) 'MAX_FLT needs to be increased to ', nFlt
           stop 99
         endif
         if (nWidth(iFlt) .gt. MAX_WIDTH) then
           write (*,*) 'MAX_WIDTH needs to be increased to ', nWidth(iFlt)
           stop 99
         endif
         if (iProb .gt. MAX_PROB) then
           write (*,*) 'MAX_PROB needs to be increased to ', iProb
           stop 99
         endif
         if (nFtypeTotal(iFlt) .gt. MAX_FTYPE) then
           write (*,*) 'MAX_FTYPE needs to be increased to ', nFtypeTotal(iFlt)
           stop 99
         endif
         if (nAtten(iFlt) .gt. MAX_ATTEN) then
           write (*,*) 'MAX_ATTEN needs to be increased to ', nAtten(iFlt)
           stop 99
         endif
         if (nWidth(iFlt) .gt. MAX_WIDTH) then
           write (*,*) 'MAX_WIDTH needs to be increased to ', nWidth(iFlt)
           stop 99
         endif
         do iWidth=1,nWidth(iFlt)
           if (nParamVar(iFlt,iWidth) .gt. MAXPARAM) then
             write (*,*) 'MAXPARAM needs to be increased to ', nParamVar(iFlt,iWidth)
             stop 99
           endif
         enddo
         if (nInten .gt. MAX_INTEN) then
           write (*,*) 'MAX_INTEN needs to be increased to ', nInten
           stop 99
         endif

          do iAtten=1,nAtten(iFlt)
            do iWidth=1,nWidth(iFlt)
               do iFtype=1,nFtypeTotal(iFlt)
                  do i=1,nParamVar(iFlt,iWidth)
 	             read (nwr,*) (temp(j),j=1,nInten)

c                    Only keep if it is the desired spectral period
	             if ( iProb .eq. iPer) then
                       do j=1,nInten
	                 haz(j,iAtten,iFlt,iWidth,i,iFtype) = temp(j)
	               enddo
                     endif
                  enddo
               enddo
	    enddo
          enddo
        enddo
      enddo

      close (nwr)

      return
  

      end

c -------------------------------------------------------------------------
