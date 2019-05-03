      
      subroutine Rd_Fault_Data (version, nFltTotal,nFlt0, probAct, AttenType, 
     2           nSegModel, Wt_segModel, nFtype, wt_ftype, f_start, f_num,
     1           faultFlag, al_Segwt, nBR_all, wt_all, fname, SR_Fact, ftype2,
     3           flt0_flag, jFlt_index, system_name )

      implicit none
      include 'tornado.h'
      
      integer Attentype(MAX_FLT),
     1        nFlt0, nFtype(MAX_FLT, MAXPARAM), faultflag(MAX_FLT1,MAX_SEG,MAX_FLT),
     2        f_start(MAX_FLT), f_num(MAX_FLT), nSegModel(MAX_FLT), nFltTotal
      real al_segWt(MAX_FLT), Wt_segModel(MAX_FLT,MAX_SEG), version, 
     1     ProbAct(MAX_FLT)
      character*80 fname(MAX_FLT) , system_name(MAX_FLT)
      real wt_all(MAX_FLT,12,MAX_THICK,MAX_BR)
      integer nBR_all(MAX_FLT,12,MAX_THICK)
      real wt_ftype(MAX_FLT,5,5), ftype2(MAX_FLT)
      real SR_Fact(MAX_FLT1,MAX_SEG,MAX_FLT)
      integer n1, n2, flt0_flag(MAX_FLT), iFlt0, i, jFLT_index(MAX_FLT)

      if (version .eq. 45.1) then
          call Rd_Fault_Data_45_1 (nFltTotal, nFlt0, Wt_all, nBR_all, 
     1           probAct, Wt_segModel, nSegModel, AttenType, wt_Ftype, 
     2           nFtype, f_start, f_num, faultFlag, al_Segwt, fname, SR_Fact, 
     3           ftype2, system_name )
    
      elseif (version .eq. 45.2 .or. version .eq. 45.3) then
          call Rd_Fault_Data_45_2 (nFltTotal, nFlt0, Wt_all, nBR_all, 
     1           probAct, Wt_segModel, nSegModel, AttenType, wt_Ftype, 
     2           nFtype, f_start, f_num, faultFlag, al_Segwt, fname, SR_Fact,
     3           ftype2, system_name  )
      else
          write (*,*) 'Incompatible fault file, use Haz45.2 or Haz45.1'
          stop 99
      endif
      
      do iFlt0=1,nFlt0
        n1 = f_start(iFlt0)
        n2 = f_start(iFlt0)+f_num(iFlt0)-1
        do i=n1,n2
          flt0_flag(i) = iFlt0
          jFlt_index(i) = i-n1+1
        enddo 
      enddo
              
      return 
      end
      
c ----------------------------------------------------------------------

      subroutine Rd_Fault_Data_45_1 (nFltTotal, nFlt0, Wt_cum_all, nBR_all, 
     1           probAct, 
     1           cumWt_segModel, nSegModel, AttenType, wt_Ftype, 
     2           nFtype, f_start, f_num, faultFlag, al_Segwt, fname, SR_Fact,
     3           ftype2, system_name  )

      implicit none
      include 'tornado.h'
      
      integer Attentype(MAX_FLT),
     1        nFlt0, nFlt2, nFm, nFtype(MAX_FLT,MAXPARAM),
     2        nThick1, directflag, synflag, 
     3        nSR, nActRAte, nRecInt, n_Dip, nRefMag(MAX_WIDTH), 
     4        faultflag(MAX_FLT1,MAX_SEG,MAX_FLT), nsyn, nFltTotal
      integer f_start(MAX_FLT), f_num(MAX_FLT), nSegModel(MAX_FLT),
     1        iflt, iflt0, iflt2, k, i, iCoor, isourceType,
     2        iFM, nRupArea, nRupWidth, iThick, nMoRate,
     3        nMagRecur, n_bValue, ii, ipt, nDownDip, insyn, iRecur,
     4        iDepthModel, nfp 
      real segwt(MAX_FLT,MAX_FLT), al_segWt(MAX_FLT), minmag, magstep, 
     1     magRecurWt(MAXPARAM), dipWt(MAXPARAM), bValueWt(MAXPARAM), minDepth, 
     4     x(MAXPARAM), ProbAct(MAX_FLT) 
      real ftmodelwt(MAXPARAM), ftype1(MAXPARAM,MAXPARAM), 
     1      magsyn, rupsyn, jbsyn, 
     2     seismosyn, hyposyn, wtsyn, wt_srBranch, wt_recIntBranch, wt_SR(MAXPARAM), 
     3     wt_RecInt(MAXPARAM), bValue2(MAXPARAM), actRate(MAXPARAM), 
     4     actRateWt(MAXPARAM), wt_MoRate(MAXPARAM)
      real faultThickWt(MAX_WIDTH), refMagWt(MAX_WIDTH,MAXPARAM), fZ,
     1     hxStep, hyStep,
     2     probAct0, sampleStep, dip1, top, fLong, fLat, wt_ActRateBranch,
     3     wt_MoRateBranch, sigArea, sigWidth, sum, attensyn, hwsyn, ftypesyn
      character*80 fName1, fName(MAX_FLT), system_name(MAX_FLT)
      real SR_Fact(MAX_FLT1,MAX_SEG,MAX_FLT)
      integer iSeg
      
      real wt_cum_all(MAX_FLT,12,MAX_THICK,MAX_BR)
      integer nBR_all(MAX_FLT,12,MAX_THICK), iNode, iBR, iOverRideMag
      real cumWt_segModel(MAX_FLT,MAX_SEG)
      real wt_ftype(MAX_FLT,5,5), ftype2(MAX_FLT)

C     Input Fault Parameters
      read (10,*) iCoor
      read (10,*) NFLT0

      call CheckDim ( NFLT0, MAX_FLT1, 'MAX_FLT1' )
  
      iflt = 0

C.....Loop over each fault in the source file....
      DO 900 iFlt0=1,NFLT0
        read (10,'( a80)') system_name(iFlt0)
        read (10,*) probAct0

c       Read number of segmentation models for this fault system       
        read (10,*) nSegModel(iFlt0)
        read (10,*) (segWt(iFlt0,k),k=1,nSegModel(iFlt0))

C       Set up cum wt for segment models for this fault system
        do k=1,nSegModel(iFlt0)
          cumWt_segModel(iFlt0,k) = segWt(iFlt0,k) 
        enddo

c       Read total number of fault segments for this fault system
        read (10,*) nFlt2
        do i=1, nSegModel(iFlt0)
         read (10,*) (faultFlag(iFlt0,i,k),k=1,nFlt2)
        enddo

c       Set the index for the first fault in this fault system and the number
c       This allows us to find the right fault from the large list
        f_start(iFlt0) = iFlt + 1
        f_num(iFlt0) = nFlt2 

C.......Loop over number of individual fault segments....                        
        do iflt2=1,nflt2

          iFlt = iFlt + 1
          call CheckDim ( iflt, MAX_FLT, 'MAX_FLT   ' )

c         Set Segment slip rate factors (not used in this version)
          do iSeg=1,nSegModel(iFlt0)
            SR_Fact(iFlt0,iSeg,iflt) = 1.0
          enddo

c         Read past name of this segment
          read(10,'( a80)') fname(iFlt)
          read (10,*) isourceType, attenType(iFlt), sampleStep, directflag, synflag

c         Read past the synchronous Rupture parameters
          if (synflag .gt. 0) then
            read (10,*) nsyn, attensyn
            do insyn=1,nsyn
              read (10,*) magsyn, rupsyn, jbsyn, seismosyn, hwsyn,ftypesyn, hyposyn, wtsyn
            enddo
          endif

c         Read aleatory segmentation wts
          read (10,*) al_segWt(iFlt)

c         Check for standard fault source or areal source
          if ( isourceType .eq. 1 .or. isourceType .eq. 2) then
            read (10,*) dip1, top
            read(10,*) nfp     
            do ipt=1,nfp
              read (10,*) fLong, fLat
            enddo
          endif

c         Check for grid source (w/o depth)
          if ( isourceType .eq. 3 .or. isourceType .eq. 7 ) then
            read (10,*)  dip1, top
c           Read past grid filename...
            read (10,*)
          endif

c         Check for grid source (w/ depth)
          if ( isourceType .eq. 4 ) then
            read (10,*)  dip1
c           Read over grid filename...
            read (10,*)
          endif

c         Check for custom fault source
          if ( isourceType .eq. 5) then
            read(10,*) nDownDip, nfp
c........   Only read in the first downdip case - rest is not needed...
            do ipt=1,nfp
              read (10,*) fLong, fLat, fZ
            enddo
          endif

c         Read dip Variation
          if ( isourceType .ne. 5 ) then
            read (10,*) n_Dip
            read (10,*) (x(i),i=1,n_Dip)
            read (10,*) (dipWt(i),i=1,n_Dip)
          else
            n_Dip = 1
            x(1) = 0.
            dipWt(1) = 1.
          endif

c         Read b-values (not for activity rate cases)
          read (10,*) n_bValue
          call CheckDim ( n_bValue, MAX_N1, 'MAX_N1    ' )
          if ( n_bValue .gt. 0 ) then
            read (10,*) (x(i),i=1,n_bValue)
            read (10,*) (bValueWt(i),i=1,n_bValue)
          endif
c          if ( n_bValue .eq. 0 ) n_bValue=1

c         Read activity rate - b-value pairs
          read (10,*) nActRate 
          if ( nActRate .ne. 0 ) then
            do ii=1,nActRate
              read (10,*) bValue2(ii), actRate(ii), actRateWt(ii)
            enddo
          endif

c         Read weights for rate methods
          read (10,*) wt_srBranch, wt_ActRateBranch, wt_recIntBranch, wt_MoRateBranch
          sum = wt_srBranch + wt_ActRateBranch + wt_recIntBranch + wt_MoRateBranch
          if ( sum .lt. 0.999 .or. sum .gt. 1.001 ) then
              write (*,'( 2x,''rate method weights do not sum to unity for fault, '',a30)') fName(iFLt)
              stop 99
          endif

c         Read slip-rates
          read (10,*) nSR
          if ( nSR .gt. 0 ) then
            read (10,*) (x(k),k=1,nSR)
            read (10,*) (wt_sr(k),k=1,nSR)
          endif

c         Read recurrence intervals
          read (10,*) nRecInt
          if ( nRecInt .gt. 0 ) then
            read (10,*) (x(k),k=1,nRecInt)
            read (10,*) (wt_recInt(k),k=1,nRecInt)
          endif

c         Read moment-rates
          read (10,*) nMoRate
          if ( nMoRate .gt. 0 ) then
            read (10,*) (x(k),k=1,nMoRate)
            read (10,*) (x(k),k=1,nMoRate)
            read (10,*) (wt_MoRate(k),k=1,nMoRate)
          endif

c         Read Mag recurrence weights (char and exp)
          read (10,*) nMagRecur
          call CheckDim ( nMagRecur, MAX_N1, 'MAX_N1    ' )
          read (10,*) (x(i),i=1,nMagRecur)
          read (10,*) (magRecurWt(i),i=1,nMagRecur)

c         Read past corresponding magnitude parameters. 
          do iRecur=1,nMagRecur
            read (10,*) x(iRecur), x(iRecur), x(iRecur)
          enddo

c         Read seismogenic thickness
          if ( isourceType .ne. 5) then
            read (10,*) nThick1
            call CheckDim ( nThick1, MAX_WIDTH, 'MAX_WIDTH ' )
            read (10,*) (x(i),i=1,nThick1)
            read (10,*) (faultThickWt(i),i=1,nThick1)
          else
            nThick1 = 1
            faultThickWt(1) = 1.
          endif

c         Read depth pdf
          read (10,*) iDepthModel       

c        Read Mag method (scaling relations or set values)
         read (10,*) iOverRideMag
         if ( iOverRideMag .ne. 1 ) then
           write (*,'( 2x,''iOverRideMag flag option not working'')') 
           stop 99
         endif

c         Read reference mags for each fault thickness
          do iThick=1,nThick1
            read (10,*) nRefMag(iThick)
            read (10,*) (x(i),i=1,nRefMag(iThick))
            read (10,*) (refMagWt(iThick,i),i=1,nRefMag(iThick))
          enddo

c         Read Past remaining input for this fault
          read (10,*) minMag, magStep, hxStep, hyStep, nRupArea, nRupWidth, minDepth
          read (10,*) (x(k),k=1,2), sigArea
          read (10,*) (x(k),k=1,2), sigWidth

c         Read ftype Models
          read (10,*) nFM
          do iFM=1,nFM
            read (10,*) ftmodelwt(iFM)
            read (10,*) nFtype(iFlt,iFM)
            read (10,*) (ftype1(iFM,k),k=1,nFtype(iFlt,iFM))
            read (10,*) (wt_ftype(iFlt,iFM,k), k=1,nFtype(iFlt,iFM))
          enddo

c         set ftype2 to just get the first case
          ftype2(iflt) = ftype1(1,1)

c         Add index on thickness because the refMag is correlated with thickness
c         For other param, just enter them for each of the thickness to make this easier to track
          do iThick=1,nThick1

c          Load seismo thickness wts into global cumulative wt array
           iNode = 1
           nBR_all(iflt,iNode,iThick) = nThick1
           do iBR=1,nBR_all(iflt,iNode,iThick)
            wt_cum_all(iflt,iNode,iThick,iBR) = faultThickWt(iBR)
           enddo

c          Load dip wts into global cumulative wt array (node 1)
           iNode = 2
           nBR_all(iflt,iNode,iThick) = n_Dip
           do iBR=1,nBR_all(iflt,iNode,iThick)
            wt_cum_all(iflt,iNode,iThick,iBR) = dipWt(iBR)
           enddo

c          Load Ftype Model wts into global cumulative wt array
           iNode = 3
           nBR_all(iflt,iNode,iThick) = nFM
           do iBR=1,nBR_all(iflt,iNode,iThick)
             wt_cum_all(iflt,iNode,iThick,iBR) = ftmodelwt(iBR)
           enddo

c          Load Mag Recur  wts into global cumulative wt array
           iNode = 4
           nBR_all(iflt,iNode,iThick) = nMagRecur
           do iBR=1,nBR_all(iflt,iNode,iThick)
             wt_cum_all(iflt,iNode,iThick,iBR) = magRecurWt(iBR)
           enddo

c          Load ref mag  wts into global cumulative wt array
           iNode = 5
           nBR_all(iflt,iNode,iThick) = nRefMag(iThick)
           do iBR=1,nBR_all(iflt,iNode,iThick)
             wt_cum_all(iflt,iNode,iThick,iBR) = refMagWt(ithick,iBR)
           enddo

c          Load rate method wts into global cumulative wt array
           iNode = 6
           nBR_all(iflt,iNode,iThick) = 4
           wt_cum_all(iflt,iNode,iThick,1) = wt_srBranch
           wt_cum_all(iflt,iNode,iThick,2) = wt_ActRateBranch 
           wt_cum_all(iflt,iNode,iThick,3) = wt_recIntBranch 
           wt_cum_all(iflt,iNode,iThick,4) = wt_MoRateBranch 

c          Load slip-rate wts into global cumulative wt array
           iNode = 7
           nBR_all(iflt,iNode,iThick) = nSR
           if ( nSR .ne. 0 ) then
             do iBR=1,nBR_all(iflt,iNode,iThick)
               wt_cum_all(iflt,iNode,iThick,iBR) = wt_sr(iBR)
             enddo
           endif
                
c          Load activity rate - b-value pairs wts into global cumulative wt array
           iNode = 8
           nBR_all(iflt,iNode,iThick) = nActRate
           if ( nActRate .ne. 0 ) then
            do iBR=1,nBR_all(iflt,iNode,iThick)
              wt_cum_all(iflt,iNode,iThick,iBR) = actRateWt(iBR)
            enddo
           endif

c          Load recurrence intervals wts into global cumulative wt array
           iNode = 9
           nBR_all(iflt,iNode,iThick) = nRecInt
           if ( nRecInt .ne. 0 ) then
            do iBR=1,nBR_all(iflt,iNode,iThick)
              wt_cum_all(iflt,iNode,iThick,iBR) = wt_recInt(iBR)
            enddo
           endif

c          Load recurrence intervals wts into global cumulative wt array
           iNode = 10
           nBR_all(iflt,iNode,iThick) = nMoRate
           if ( nMoRate .ne. 0 ) then
            do iBR=1,nBR_all(iflt,iNode,iThick)
              wt_cum_all(iflt,iNode,iThick,iBR) = wt_MoRate(iBR)
            enddo
           endif
                       
c          Load b-values wts into global cumulative wt array (node 11)
           iNode = 11
           nBR_all(iflt,iNode,iThick) = n_bValue
           do iBR=2,nBR_all(iflt,iNode,iThick)
            wt_cum_all(iflt,iNode,iThick,iBR) = bValueWt(iBR)
           enddo

c          End of Loop over iThick1
          enddo
                   
         probAct(iFlt) = probAct0
         
c       End of Loop over iFlt2 - number of segments    
        enddo

c     End of Loop over iFlt
  900 continue
      nFltTotal = iFlt
      
      return
      end

c -------------------

      subroutine Rd_Fault_Data_45_2 (nFltTotal, nFlt0, Wt_cum_all, nBR_all, 
     1           probAct, 
     1           cumWt_segModel, nSegModel, AttenType, wt_Ftype, 
     2           nFtype, f_start, f_num, faultFlag, al_Segwt, fname, SR_Fact,
     3           ftype2, system_name )

      implicit none
      include 'tornado.h'
      
      integer Attentype(MAX_FLT),
     1        nFlt0, nFlt2, nFm, nFtype(MAX_FLT,MAXPARAM),
     2        nThick1, directflag, synflag, 
     3        nSR, nActRAte, nRecInt, n_Dip, nRefMag(MAX_WIDTH), 
     4        faultflag(MAX_FLT1,MAX_SEG,MAX_FLT), nsyn, nFltTotal
      integer f_start(MAX_FLT), f_num(MAX_FLT), nSegModel(MAX_FLT),
     1        iflt, iflt0, iflt2, k, i, iCoor, isourceType,
     2        iFM, nRupArea, nRupWidth, iThick, nMoRate,
     3        nMagRecur, n_bValue, ii, ipt, nDownDip, insyn, iRecur,
     4        iDepthModel, nfp 
      real segwt(MAX_FLT,MAX_FLT), al_segWt(MAX_FLT), minmag, magstep, 
     1     magRecurWt(MAXPARAM), dipWt(MAXPARAM), bValueWt(MAXPARAM), minDepth, 
     4     x(MAXPARAM), ProbAct(MAX_FLT) 
      real ftmodelwt(MAXPARAM), ftype1(MAXPARAM,MAXPARAM), 
     1      magsyn, rupsyn, jbsyn, 
     2     seismosyn, hyposyn, wtsyn, wt_srBranch, wt_recIntBranch, wt_SR(MAXPARAM), 
     3     wt_RecInt(MAXPARAM), bValue2(MAXPARAM), actRate(MAXPARAM), 
     4     actRateWt(MAXPARAM), wt_MoRate(MAXPARAM)
      real faultThickWt(MAX_WIDTH), refMagWt(MAX_WIDTH,MAXPARAM), fZ,
     1     hxStep, hyStep,
     2     probAct0, sampleStep, dip1, top, fLong, fLat, wt_ActRateBranch,
     3     wt_MoRateBranch, sigArea, sigWidth, sum, attensyn, hwsyn, ftypesyn
      character*80 fName1, fName(MAX_FLT), dummy, system_name(MAX_FLT)
      real actRate0
      real SR_Fact(MAX_FLT1,MAX_SEG,MAX_FLT), xx(MAX_SEG,MAX_FLT)
      
      real wt_cum_all(MAX_FLT,12,MAX_THICK,MAX_BR)
      integer nBR_all(MAX_FLT,12,MAX_THICK), iNode, iBR
      real cumWt_segModel(MAX_FLT,MAX_SEG)
      real wt_ftype(MAX_FLT,5,5), ftype2(MAX_FLT)
      integer IST5, iSR_Fact, ix, iST5_flag, iSeg
      
C     Input Fault Parameters
      read (10,*,err=1000) iCoor
      read (10,*,err=1001) NFLT0

      call CheckDim ( NFLT0, MAX_FLT1, 'MAX_FLT1' )
  
      iflt = 0

C.....Loop over each fault in the source file....
      DO 900 iFlt0=1,NFLT0
        read (10,'( a80)',err=1002) system_name(iFlt0)
        read (10,*,err=1003) probAct0

c       Read number of segmentation models for this fault system       
        iSR_Fact = 0
        read (10,*,err=1004) nSegModel(iFlt0)
        if (nSegModel(iFlt0) .eq. -1 ) then
         iSR_Fact = 1
         backspace ( 10)
         read (10,*,err=3005) ix, nSegModel(iFlt0)
        endif

        if (nSegModel(iFlt0) .ge. MAX_seg) then
          write (*,'( 2x,'' dimension error, increase max_seg to '',i5)')nSegModel(iFlt0)
          stop 99
        endif 
        read (10,*,err=1005) (segWt(iFlt0,k),k=1,nSegModel(iFlt0))

C       Set up cum wt for segment models for this fault system
        do k=1,nSegModel(iFlt0)
          cumWt_segModel(iFlt0,k) = segWt(iFlt0,k) 
        enddo

c       Read total number of fault segments for this fault system
        read (10,*,err=1006) nFlt2
        do i=1, nSegModel(iFlt0)
         read (10,*,err=1007) (faultFlag(iFlt0,i,k),k=1,nFlt2)
        enddo

c       If SR factors used, then read these values
        if ( iSR_Fact .eq. 1 ) then
         do i=1,nSegModel(iFlt0)
           read (10,*,err=3059) (xx(i,k),k=1,nFlt2)
         enddo
        else
c         Set to unity it not used
          do i=1,nSegModel(iFlt0)
            do k=1,nFlt2
              xx(i,k) = 1.0
            enddo
          enddo
        endif       

c       Set the index for the first fault in this fault system and the number
c       This allows us to find the right fault from the large list
        f_start(iFlt0) = iFlt + 1
        f_num(iFlt0) = nFlt2 

C.......Loop over number of individual fault segments....                        
        do iflt2=1,nflt2

          iFlt = iFlt + 1
          call CheckDim ( iflt, MAX_FLT, 'MAX_FLT   ' )

c         Set Segment slip rate factors
          do iSeg=1,nSegModel(iFlt0)
            SR_Fact(iFlt0,iSeg,iflt) = xx(iSeg,iflt2)
          enddo

c         Read past name of this segment
          read(10,'( a80)',err=1008) fname(iFlt)
          read (10,*,err=1009) isourceType, attenType(iFlt), sampleStep, directflag, synflag

c       Set flag for source type 5 to allow limit on crustal thickness
        iST5_flag = 0
        if (isourceType .eq. -5 ) then
          isourceType = 5
          iST5_flag = 1
        endif   

c         Read past the synchronous Rupture parameters
          if (synflag .gt. 0) then
            read (10,*) nsyn, attensyn
            do insyn=1,nsyn
              read (10,*,err=1010) magsyn, rupsyn, jbsyn, seismosyn, hwsyn,ftypesyn, hyposyn, wtsyn
            enddo
          endif

c         Read aleatory segmentation wts
          read (10,*,err=1011) al_segWt(iFlt)

c         Check for standard fault source or areal source
          if ( isourceType .eq. 1 .or. isourceType .eq. 2) then
            read (10,*,err=1012) dip1, top
            read(10,*,err=1013) nfp     
            do ipt=1,nfp
              read (10,*,err=1014) fLong, fLat
            enddo
          endif

c         Check for grid source (w/o depth)
          if ( isourceType .eq. 3 .or. isourceType .eq. 7 ) then
            read (10,*,err=1015)  dip1, top
c           Read past grid filename...
            read (10,*,err=1016) dummy
          endif

c         Check for grid source (w/ depth)
          if ( isourceType .eq. 4 ) then
            read (10,*,err=1017)  dip1
c           Read over grid filename...
            read (10,*,err=1018) dummy
          endif

c         Check for custom fault source
          if ( isourceType .eq. 5) then
            read(10,*,err=1019) nDownDip, nfp
c........   Only read in the first downdip case - rest is not needed...
            do ipt=1,nfp
              read (10,*,err=1020) fLong, fLat, fZ
            enddo
          endif

c         Read dip Variation
          if ( isourceType .ne. 5 ) then
            read (10,*,err=1021) n_Dip
            read (10,*,err=1022) (x(i),i=1,n_Dip)
            read (10,*,err=1023) (dipWt(i),i=1,n_Dip)
          else
            n_Dip = 1
            x(1) = 0.
            dipWt(1) = 1.
          endif

c         Read b-values (not for activity rate cases)
          read (10,*,err=1024) n_bValue
          call CheckDim ( n_bValue, MAX_N1, 'MAX_N1    ' )
          if ( n_bValue .gt. 0 ) then
            read (10,*,err=1025) (x(i),i=1,n_bValue)
            read (10,*,err=1026) (bValueWt(i),i=1,n_bValue)
          endif

c         Read activity rate - b-value pairs
          read (10,*) nActRate 
          if ( nActRate .ne. 0 ) then
            do ii=1,nActRate
              read (10,*,err=1027) bValue2(ii), actRate(ii), actRateWt(ii)
            enddo
          endif

c         Read weights for rate methods
          read (10,*,err=1028) wt_srBranch, wt_ActRateBranch, wt_recIntBranch, wt_MoRateBranch
          sum = wt_srBranch + wt_ActRateBranch + wt_recIntBranch + wt_MoRateBranch
          if ( sum .lt. 0.999 .or. sum .gt. 1.001 ) then
              write (*,'( 2x,''rate method weights do not sum to unity for fault, '',a30)') fName(iFLt)
              stop 99
          endif

c         Read slip-rates
          read (10,*) nSR
          if ( nSR .gt. 0 ) then
            read (10,*,err=1029) (x(k),k=1,nSR)

c          Check if this is a vertical slip rate (indicated by SR(1)=-999)
           if ( x(1) .eq. -999. ) then
             read (10,*,err=3058) (x(k),k=1,nSR)
           endif

            read (10,*,err=1030) (wt_sr(k),k=1,nSR)
          endif

c         Read recurrence intervals
          read (10,*) nRecInt
          if ( nRecInt .gt. 0 ) then
            read (10,*,err=1031) (x(k),k=1,nRecInt)
            read (10,*,err=1032) (wt_recInt(k),k=1,nRecInt)
          endif

c         Read moment-rates
          read (10,*) nMoRate
          if ( nMoRate .gt. 0 ) then
            read (10,*,err=1033) (x(k),k=1,nMoRate)
            read (10,*,err=1034) (x(k),k=1,nMoRate)
            read (10,*,err=1035) (wt_MoRate(k),k=1,nMoRate)
          endif

c         Read Mag recurrence weights (char and exp)
          read (10,*,err=1036) nMagRecur
          call CheckDim ( nMagRecur, MAX_N1, 'MAX_N1    ' )
          read (10,*,err=1037) (x(i),i=1,nMagRecur)
          read (10,*,err=1038) (magRecurWt(i),i=1,nMagRecur)

c         Read past corresponding magnitude parameters. 
          do iRecur=1,nMagRecur
            read (10,*,err=1039) x(iRecur), x(iRecur), x(iRecur)
          enddo

c         Read seismogenic thickness
          if ( isourceType .ne. 5 .or. iST5_flag .eq. 1) then
            read (10,*,err=1040) nThick1
            call CheckDim ( nThick1, MAX_THICK, 'MAX_THICK ' )
            read (10,*,err=1041) (x(i),i=1,nThick1)
            read (10,*,err=1042) (faultThickWt(i),i=1,nThick1)
          else
            nThick1 = 1
            faultThickWt(1) = 1.
          endif

c         Read depth pdf
          read (10,*,err=1043) iDepthModel  

c         Read reference mags for each fault thickness
          do iThick=1,nThick1
            read (10,*,err=1044) nRefMag(iThick)
            read (10,*,err=1045) (x(i),i=1,nRefMag(iThick))
            read (10,*,err=1046) (refMagWt(iThick,i),i=1,nRefMag(iThick))
          enddo

c         Read Past remaining input for this fault
          read (10,*,err=1047) minMag, magStep, hxStep, hyStep, nRupArea, nRupWidth, minDepth
          read (10,*,err=1048) (x(k),k=1,2), sigArea
          read (10,*,err=1049) (x(k),k=1,2), sigWidth

c         Read ftype Models
          read (10,*) nFM
          do iFM=1,nFM
            read (10,*,err=1050) ftmodelwt(iFM)
            read (10,*,err=1051) nFtype(iFlt,iFM)
            read (10,*,err=1052) (ftype1(iFM,k),k=1,nFtype(iFlt,iFM))
            read (10,*,err=1053) (wt_ftype(iFlt,iFM,k), k=1,nFtype(iFlt,iFM))
          enddo

c         set ftype2 to just get the first case
          ftype2(iflt) = ftype1(1,1)

c         Add index on thickness because the refMag is correlated with thickness
c         For other param, just enter them for each of the thickness to make this easier to track
          do iThick=1,nThick1
          
c          Load seismo thickness wts into global cumulative wt array
           iNode = 1
           nBR_all(iflt,iNode,iThick) = nThick1
           if ( nBR_all(iflt,iNode,iThick) .gt. MAX_BR ) goto 4000
           do iBR=1,nBR_all(iflt,iNode,iThick)
            wt_cum_all(iflt,iNode,iThick,iBR) = faultThickWt(iBR)
           enddo

c          Load dip wts into global cumulative wt array (node 1)
           iNode = 2
           nBR_all(iflt,iNode,iThick) = n_Dip
           if ( nBR_all(iflt,iNode,iThick) .gt. MAX_BR ) goto 4000
           do iBR=1,nBR_all(iflt,iNode,iThick)
            wt_cum_all(iflt,iNode,iThick,iBR) = dipWt(iBR)
           enddo

c          Load Ftype Model wts into global cumulative wt array
           iNode = 3
           nBR_all(iflt,iNode,iThick) = nFM
           if ( nBR_all(iflt,iNode,iThick) .gt. MAX_BR ) goto 4000
           do iBR=1,nBR_all(iflt,iNode,iThick)
             wt_cum_all(iflt,iNode,iThick,iBR) = ftmodelwt(iBR)
           enddo

c          Load Mag Recur  wts into global cumulative wt array
           iNode = 4
           nBR_all(iflt,iNode,iThick) = nMagRecur
           if ( nBR_all(iflt,iNode,iThick) .gt. MAX_BR ) goto 4000
           do iBR=1,nBR_all(iflt,iNode,iThick)
             wt_cum_all(iflt,iNode,iThick,iBR) = magRecurWt(iBR)
           enddo

c          Load ref mag  wts into global cumulative wt array
           iNode = 5
           nBR_all(iflt,iNode,iThick) = nRefMag(iThick)
           if ( nBR_all(iflt,iNode,iThick) .gt. MAX_BR ) goto 4000
           do iBR=1,nBR_all(iflt,iNode,iThick)
             wt_cum_all(iflt,iNode,iThick,iBR) = refMagWt(ithick,iBR)
           enddo

c          Load rate method wts into global cumulative wt array
           iNode = 6
           nBR_all(iflt,iNode,iThick) = 4
           wt_cum_all(iflt,iNode,iThick,1) = wt_srBranch
           wt_cum_all(iflt,iNode,iThick,2) = wt_ActRateBranch 
           wt_cum_all(iflt,iNode,iThick,3) = wt_recIntBranch 
           wt_cum_all(iflt,iNode,iThick,4) = wt_MoRateBranch 

c          Load slip-rate wts into global cumulative wt array
           iNode = 7
           nBR_all(iflt,iNode,iThick) = nSR
           if ( nBR_all(iflt,iNode,iThick) .gt. MAX_BR ) goto 4000
           if ( nSR .ne. 0 ) then
             do iBR=1,nBR_all(iflt,iNode,iThick)
               wt_cum_all(iflt,iNode,iThick,iBR) = wt_sr(iBR)
             enddo
           endif
                
c          Load activity rate - b-value pairs wts into global cumulative wt array
           iNode = 8
           nBR_all(iflt,iNode,iThick) = nActRate
           if ( nBR_all(iflt,iNode,iThick) .gt. MAX_BR ) goto 4000
           if ( nActRate .ne. 0 ) then
            do iBR=1,nBR_all(iflt,iNode,iThick)
              wt_cum_all(iflt,iNode,iThick,iBR) = actRateWt(iBR)
            enddo
           endif

c          Load recurrence intervals wts into global cumulative wt array
           iNode = 9
           nBR_all(iflt,iNode,iThick) = nRecInt
           if ( nBR_all(iflt,iNode,iThick) .gt. MAX_BR ) goto 4000
           if ( nRecInt .ne. 0 ) then
            do iBR=1,nBR_all(iflt,iNode,iThick)
              wt_cum_all(iflt,iNode,iThick,iBR) = wt_recInt(iBR)
            enddo
           endif

c          Load recurrence intervals wts into global cumulative wt array
           iNode = 10
           nBR_all(iflt,iNode,iThick) = nMoRate
           if ( nBR_all(iflt,iNode,iThick) .gt. MAX_BR ) goto 4000
           if ( nMoRate .ne. 0 ) then
            do iBR=1,nBR_all(iflt,iNode,iThick)
              wt_cum_all(iflt,iNode,iThick,iBR) = wt_MoRate(iBR)
            enddo
           endif
                       
c          Load b-values wts into global cumulative wt array (node 11)
           iNode = 11
           nBR_all(iflt,iNode,iThick) = n_bValue
           if ( nBR_all(iflt,iNode,iThick) .gt. MAX_BR ) goto 4000
           do iBR=1,nBR_all(iflt,iNode,iThick)
            wt_cum_all(iflt,iNode,iThick,iBR) = bValueWt(iBR)
           enddo

c          End of Loop over iThick1
          enddo
                   
         probAct(iFlt) = probAct0
         
c       End of Loop over iFlt2 - number of segments    
        enddo

c     End of Loop over iFlt
  900 continue
      nFltTotal = iFlt
      
      return
 1000 write (*,'( 2x,''Error 1000'')')
      backspace (20)
      read (20,'( a80)') dummy
      write (*,'( a80)') dummy
      stop 99
 1001 write (*,'( 2x,''Error 1001'')')
      backspace (20)
      read (20,'( a80)') dummy
      write (*,'( a80)') dummy
      stop 99
 1002 write (*,'( 2x,''Error 1002'')')
      backspace (20)
      read (20,'( a80)') dummy
      write (*,'( a80)') dummy
      stop 99
 1003 write (*,'( 2x,''Error 1003'')')
      backspace (20)
      read (20,'( a80)') dummy
      write (*,'( a80)') dummy
      stop 99
 1004 write (*,'( 2x,''Error 1004'')')
      backspace (20)
      read (20,'( a80)') dummy
      write (*,'( a80)') dummy
      stop 99
 1005 write (*,'( 2x,''Error 1005'')')
      backspace (20)
      read (20,'( a80)') dummy
      write (*,'( a80)') dummy
      stop 99
 1006 write (*,'( 2x,''Error 1006'')')
      backspace (20)
      read (20,'( a80)') dummy
      write (*,'( a80)') dummy
      stop 99
 1007 write (*,'( 2x,''Error 1007'')')
      backspace (20)
      read (20,'( a80)') dummy
      write (*,'( a80)') dummy
      stop 99
 1008 write (*,'( 2x,''Error 1008'')')
      backspace (20)
      read (20,'( a80)') dummy
      write (*,'( a80)') dummy
      stop 99
 1009 write (*,'( 2x,''Error 1009'')')
      backspace (20)
      read (20,'( a80)') dummy
      write (*,'( a80)') dummy
      stop 99
 1010 write (*,'( 2x,''Error 1010'')')
      backspace (20)
      read (20,'( a80)') dummy
      write (*,'( a80)') dummy
      stop 99
 1011 write (*,'( 2x,''Error 1011'')')
      backspace (20)
      read (20,'( a80)') dummy
      write (*,'( a80)') dummy
      stop 99
 1012 write (*,'( 2x,''Error 1012'')')
      backspace (20)
      read (20,'( a80)') dummy
      write (*,'( a80)') dummy
      stop 99
 1013 write (*,'( 2x,''Error 1013'')')
      backspace (20)
      read (20,'( a80)') dummy
      write (*,'( a80)') dummy
      stop 99
 1014 write (*,'( 2x,''Error 1014'')')
      backspace (20)
      read (20,'( a80)') dummy
      write (*,'( a80)') dummy
      stop 99
 1015 write (*,'( 2x,''Error 1015'')')
      backspace (20)
      read (20,'( a80)') dummy
      write (*,'( a80)') dummy
      stop 99
 1016 write (*,'( 2x,''Error 1016'')')
      backspace (20)
      read (20,'( a80)') dummy
      write (*,'( a80)') dummy
      stop 99
 1017 write (*,'( 2x,''Error 1017'')')
      backspace (20)
      read (20,'( a80)') dummy
      write (*,'( a80)') dummy
      stop 99
 1018 write (*,'( 2x,''Error 1018'')')
      backspace (20)
      read (20,'( a80)') dummy
      write (*,'( a80)') dummy
      stop 99
 1019 write (*,'( 2x,''Error 1019'')')
      backspace (20)
      read (20,'( a80)') dummy
      write (*,'( a80)') dummy
      stop 99
 1020 write (*,'( 2x,''Error 1020'')')
      backspace (20)
      read (20,'( a80)') dummy
      write (*,'( a80)') dummy
      stop 99
 1021 write (*,'( 2x,''Error 1021'')')
      backspace (20)
      read (20,'( a80)') dummy
      write (*,'( a80)') dummy
      stop 99
 1022 write (*,'( 2x,''Error 1022'')')
      backspace (20)
      read (20,'( a80)') dummy
      write (*,'( a80)') dummy
      stop 99
 1023 write (*,'( 2x,''Error 1023'')')
      backspace (20)
      read (20,'( a80)') dummy
      write (*,'( a80)') dummy
      stop 99
 1024 write (*,'( 2x,''Error 1024'')')
      backspace (20)
      read (20,'( a80)') dummy
      write (*,'( a80)') dummy
      stop 99
 1025 write (*,'( 2x,''Error 1025'')')
      backspace (20)
      read (20,'( a80)') dummy
      write (*,'( a80)') dummy
      stop 99
 1026 write (*,'( 2x,''Error 1026'')')
      backspace (20)
      read (20,'( a80)') dummy
      write (*,'( a80)') dummy
      stop 99
 1027 write (*,'( 2x,''Error 1027'')')
      backspace (20)
      read (20,'( a80)') dummy
      write (*,'( a80)') dummy
      stop 99
 1028 write (*,'( 2x,''Error 1028'')')
      backspace (20)
      read (20,'( a80)') dummy
      write (*,'( a80)') dummy
      stop 99
 1029 write (*,'( 2x,''Error 1029'')')
      backspace (20)
      read (20,'( a80)') dummy
      write (*,'( a80)') dummy
      stop 99
 1030 write (*,'( 2x,''Error 1030'')')
      backspace (20)
      read (20,'( a80)') dummy
      write (*,'( a80)') dummy
      stop 99
 1031 write (*,'( 2x,''Error 1031'')')
      backspace (20)
      read (20,'( a80)') dummy
      write (*,'( a80)') dummy
      stop 99
 1032 write (*,'( 2x,''Error 1032'')')
      backspace (20)
      read (20,'( a80)') dummy
      write (*,'( a80)') dummy
      stop 99
 1033 write (*,'( 2x,''Error 1033'')')
      backspace (20)
      read (20,'( a80)') dummy
      write (*,'( a80)') dummy
      stop 99
 1034 write (*,'( 2x,''Error 1034'')')
      backspace (20)
      read (20,'( a80)') dummy
      write (*,'( a80)') dummy
      stop 99
 1035 write (*,'( 2x,''Error 1035'')')
      backspace (20)
      read (20,'( a80)') dummy
      write (*,'( a80)') dummy
      stop 99
 1036 write (*,'( 2x,''Error 1036'')')
      backspace (20)
      read (20,'( a80)') dummy
      write (*,'( a80)') dummy
      stop 99
 1037 write (*,'( 2x,''Error 1037'')')
      backspace (20)
      read (20,'( a80)') dummy
      write (*,'( a80)') dummy
      stop 99
 1038 write (*,'( 2x,''Error 1038'')')
      backspace (20)
      read (20,'( a80)') dummy
      write (*,'( a80)') dummy
      stop 99
 1039 write (*,'( 2x,''Error 1039'')')
      backspace (20)
      read (20,'( a80)') dummy
      write (*,'( a80)') dummy
      stop 99
 1040 write (*,'( 2x,''Error 1040'')')
      backspace (20)
      read (20,'( a80)') dummy
      write (*,'( a80)') dummy
      stop 99
 1041 write (*,'( 2x,''Error 1041'')')
      backspace (20)
      read (20,'( a80)') dummy
      write (*,'( a80)') dummy
      stop 99
 1042 write (*,'( 2x,''Error 1042'')')
      backspace (20)
      read (20,'( a80)') dummy
      write (*,'( a80)') dummy
      stop 99
 1043 write (*,'( 2x,''Error 1043'')')
      backspace (20)
      read (20,'( a80)') dummy
      write (*,'( a80)') dummy
      stop 99
 1044 write (*,'( 2x,''Error 1044'')')
      backspace (20)
      read (20,'( a80)') dummy
      write (*,'( a80)') dummy
      stop 99
 1045 write (*,'( 2x,''Error 1045'')')
      backspace (20)
      read (20,'( a80)') dummy
      write (*,'( a80)') dummy
      stop 99
 1046 write (*,'( 2x,''Error 1046'')')
      backspace (20)
      read (20,'( a80)') dummy
      write (*,'( a80)') dummy
      stop 99
 1047 write (*,'( 2x,''Error 1047'')')
      backspace (20)
      read (20,'( a80)') dummy
      write (*,'( a80)') dummy
      stop 99
 1048 write (*,'( 2x,''Error 1048'')')
      backspace (20)
      read (20,'( a80)') dummy
      write (*,'( a80)') dummy
      stop 99
 1049 write (*,'( 2x,''Error 1049'')')
      backspace (20)
      read (20,'( a80)') dummy
      write (*,'( a80)') dummy
      stop 99
 1050 write (*,'( 2x,''Error 1050'')')
      backspace (20)
      read (20,'( a80)') dummy
      write (*,'( a80)') dummy
      stop 99
 1051 write (*,'( 2x,''Error 1051'')')
      backspace (20)
      read (20,'( a80)') dummy
      write (*,'( a80)') dummy
      stop 99
 1052 write (*,'( 2x,''Error 1052'')')
      backspace (20)
      read (20,'( a80)') dummy
      write (*,'( a80)') dummy
      stop 99
 1053 write (*,'( 2x,''Error 1053'')')
      backspace (20)
      read (20,'( a80)') dummy
      write (*,'( a80)') dummy
      stop 99
 3005 write (*,'( 2x,''Flt file error:  nSegModel'')')
      write (*,'( 2x,''fault sys: '',a80)') fName1
      stop 99
 3058 write (*,'( 2x,''Flt file error: SR rakes for vertical SR'', 2i5)') iFlt0, iFLt2
      write (*,'( 2x,''From fault segment: '',a80)') fName(iFlt)
      stop 99
 3059  write (*,'( 2x,''Flt file error: SR factors for seg model'', i5)') i
       write (*,'( 2x,''From fault: '',a80)') fName1
       stop 99
 4000 write (*,'( 2x,''dimension error, increase MAX_BR'')')
      write (*,'( 4i5)') iFlt, iNode, iThick, nBR_all(iFlt,iNode,iThick)
      write (*,'( 2x,'' iFlt, iNode, iThick, nBR_all'')')
      stop 99      
     
      end

      subroutine Check_sumwt ( wt_cum_all, iFlt, iNode, iTHick, nBR_all )
      implicit none
      include 'tornado.h'
      real wt_cum_all(MAX_FLT,12,5,10)
      integer nBR_all(MAX_FLT,12,5)
      integer iFlt, iNode, iThick, n1, k
      
      n1 = nBR_all(iflt,iNode,iThick)
      if ( n1 .eq. 0 ) return
      
c      write (*,'( 4i5)') iflt,iNode,iThick, nBR_all(iflt,iNode,iThick)
c      pause 'test nBR'
      if ( wt_cum_all(iFlt,iNode,iThick,n1) .ne. 1. ) then
        write (*,'( 2x,''SSC logic tree weights do not sum to unity'')')
        write (*,'( 2x,''iFlt, iNode, iTHick ='',3i5)') iFlt, iNode, iThick
        write (*,'( 2x,''Cumulative weights:'')')
        write (*,'( 10f10.4)') (wt_cum_all(iFlt,iNode,iThick,k),k=1,n1)
        stop 99
      endif
      return
      end
