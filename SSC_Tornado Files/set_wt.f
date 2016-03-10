      subroutine Set_starting_wts ( n_dip, dipWt, dip_wt1, 
     1     nThick, faultThickWt, faultThick_Wt1, nRefMag, refMagWt, refMag_Wt1,
     2     n_bValue, bvalueWt, bvalue_Wt1, nMagRecur, magRecurWt, magRecur_Wt1,
     3     nFtypeModels, nFtype, ftModelwt, ftype_al, ftype_wt1,
     4     nSegModel, segwt, seg_wt1, wt_RateMethod, wt_RateMethod1,
     5     nSR, wt_SR, wt_SR1, nActRate, actRateWt, actRate_Wt1,
     6     nRecInt, wt_recInt, wt_recInt1, nMoRate, wt_MoRate, MoRate_wt1,
     8     segFlag, totalSegWt, jFlt )
      implicit none

      include 'tornado.h'

      integer n_Dip(MAX_FLT),n_bvalue(MAX_FLT), nActRate(MAX_FLT), nSR(MAX_FLT), 
     1        nRecInt(MAX_FLT), nMoRate(MAX_FLT),
     1        nRefMag(MAX_FLT,MAX_WIDTH), nFtypeModels(MAX_FLT)
      integer nFtype(MAX_FLT,MAXPARAM), nThick(MAX_FLT), nMagRecur(MAX_FLT)
      real ftype_al(MAX_FLt,MAXPARAM,5)
      integer nSegModel(MAX_FLT), segFlag(MAX_FLT, MAX_SEG)

      real segwt(MAX_FLT,MAX_SEG), dipWt(MAX_FLT,MAXPARAM), bValueWt(MAX_FLT,MAXPARAM),
     1     actRateWt(MAX_FLT,MAXPARAM), 
     1     wt_SR(MAX_FLT,MAXPARAM), wt_RecInt(MAX_FLT,MAXPARAM), wt_MoRate(MAX_FLT,MAXPARAM),
     1      magRecurWt(MAX_FLT,MAXPARAM),    
     1     faultThickWt(MAX_FLT,MAXPARAM), 
     2     refMagWt(MAX_FLT,MAX_Width, MAXPARAM), ftmodelwt(MAX_FLT,MAXPARAM), 
     1     wt_rateMethod(MAX_FLT,4)

      real dip_Wt1(MAX_FLT,MAXPARAM), bValue_Wt1(MAX_FLT,MAXPARAM), actRate_Wt1(MAX_FLT,MAXPARAM), 
     1     wt_SR1(MAX_FLT,MAXPARAM), wt_RecInt1(MAX_FLT,MAXPARAM), MoRate_wt1(MAX_FLT,MAXPARAM), 
     1     magRecur_Wt1(MAX_FLT,MAXPARAM), faultThick_Wt1(MAX_FLT,MAXPARAM), 
     2     refMag_Wt1(MAX_FLT,MAX_Width, MAXPARAM), ftype_wt1(MAX_FLT,MAXPARAM)
      real seg_wt1(MAX_FLT,MAX_SEG)
      real wt_rateMethod1(MAX_FLT,4)
      
      integer j, jFlt, iThick, k, iFM, i, ib
      real sum, totalSegWt(MAX_FLT)
     
c        Set Dip weights
         do j=1,n_Dip(jFlt)
           dip_Wt1(jFlt,j) = dipWt(jFlt,j)
         enddo

c        Set thick and refMag weights
         do iThick=1,nThick(jFlt)
           faultThick_Wt1(jFlt,iThick) = faultThickWt(jFlt,iThick)
           do j=1,nRefMag(jFlt,iThick) 
             refMag_Wt1(jFlt,iThick,j)  = refMagWt(jFlt,iThick,j) 
           enddo
         enddo

c        Set b-values weights
         do j=1,n_bValue(jFlt)
           if ( n_bValue(jFlt) .eq. 1 ) then
             bvalue_Wt1(jFlt,j) = 1.
             bvalueWt(jFlt,j) = 1.
           else
             bvalue_Wt1(jFlt,j) = bvalueWt(jFlt,j)
           endif
         enddo

c        Set mag recur weights
         do j=1,nMagRecur(jFlt)
           magRecur_Wt1(jFlt,j)  = magRecurWt(jFlt,j) 
         enddo

c        Set ftype weights
         k = 0
         do iFM=1,nFtypeModels(jFlt)
           do j=1,nFtype(jFlt,IFM)
             k = k + 1
             ftype_wt1 (jFlt,k) =  ftModelwt(jFlt, iFM) * ftype_al(jFlt,iFM,k)
           enddo
         enddo

c        set segmentation weights
         do i=1,nSegModel(jFlt)
          seg_wt1(jFlt,i)= segwt(jFlt,i)
         enddo

c        Set Activity Rate Method weights
         do i=1,4
           wt_RateMethod1(jFlt,i) = wt_RateMethod(jFLt,i)
         enddo

c        Set SR weights 
         do i=1,nSR(jFLt)
           if (nSR(jFLt) .eq. 1 ) then
             wt_SR1(jFLT,i) = 1.
             wt_SR(jFlt,i) = 1.
           else
             wt_SR1(jFLT,i) = wt_SR(jFlt,i)
           endif
         enddo

c        Set Activity Rate weights
         do i=1,nActRate(jFLt)
           if ( nActRate(jFlt) .eq. 1 ) then
             actRate_wt1(jFLT,i) = 1.
             actRateWt(jFlt,i) = 1.
           else
             actRate_wt1(jFLT,i) = actRateWt(jFlt,i)
           endif
         enddo

c        Set Rec Int weights
         do i=1,nRecInt(jFLt)
           if ( nRecInt(jFlt) .eq. 1 ) then
             wt_RecInt1(jFLT,i) = 1.
             wt_recInt(jFlt,i) = 1.
           else
             wt_RecInt1(jFLT,i) = wt_recInt(jFlt,i)
           endif
         enddo

c        Set moment rate weights
         do i=1,nMoRate(jFLt)
           if ( nMoRate(jFLt) .eq. 1 ) then
             MoRate_wt1(jFLT,i) = 1.
             wt_MoRate(jFLT,i) = 1.
           else
             MoRate_wt1(jFLT,i) = wt_MoRate(jFLT,i)
           endif 
         enddo

c        Find the total weight for this fault
         sum = 0.
         do i=1,nSegModel(jFlt)
           sum = sum + seg_wt1(jFlt,i) * segFlag(jFlt,i)
         enddo
         totalSegWt(jFlt) = sum
c         write (*,'( i5,f10.4)') jflt, totalSegWt(jFlt) 

      return
      end

c -----------------------------------------------------------

      subroutine Set_new_wts ( n_dip, dipWt, dip_wt1, 
     1     nThick, faultThickWt, faultThick_Wt1, nRefMag, refMagWt,refMag_Wt1,
     2     n_bValue, bvalueWt, bvalue_Wt1, nMagRecur, magRecurWt, magRecur_Wt1,
     3     nFtypeModels, nFtype, ftModelwt, ftype_al, ftype_wt1,
     4     nSegModel, segwt, seg_wt1, wt_RateMethod, wt_RateMethod1,
     5     nSR, wt_SR, wt_SR1, nActRate, actRateWt, actRate_Wt1,
     6     nRecInt, wt_recInt, wt_recInt1, nMoRate, wt_MoRate, MoRate_wt1,
     7     segFlag, totalSegWt, kFlt, iNode, iBR, nBR_SSC, iSkip )
      implicit none

      include 'tornado.h'

      integer n_Dip(MAX_FLT),n_bvalue(MAX_FLT), nActRate(MAX_FLT), nSR(MAX_FLT), 
     1        nRecInt(MAX_FLT), nMoRate(MAX_FLT),
     1        nRefMag(MAX_FLT,MAX_WIDTH), nFtypeModels(MAX_FLT)
      integer nFtype(MAX_FLT,MAXPARAM), nThick(MAX_FLT), nMagRecur(MAX_FLT)
      real ftype_al(MAX_FLt,MAXPARAM,5)
      integer nSegModel(MAX_FLT), segFlag(MAX_FLT, MAX_SEG)

      real segwt(MAX_FLT,MAX_FLT), dipWt(MAX_FLT,MAXPARAM), bValueWt(MAX_FLT,MAXPARAM),
     1     actRateWt(MAX_FLT,MAXPARAM), 
     1     wt_SR(MAX_FLT,MAXPARAM), wt_RecInt(MAX_FLT,MAXPARAM), wt_MoRate(MAX_FLT,MAXPARAM),
     1      magRecurWt(MAX_FLT,MAXPARAM),    
     1     faultThickWt(MAX_FLT,MAXPARAM), 
     2     refMagWt(MAX_FLT,MAX_Width, MAXPARAM), ftmodelwt(MAX_FLT,MAXPARAM), 
     1     wt_rateMethod(MAX_FLT,4)

      real dip_Wt1(MAX_FLT,MAXPARAM), bValue_Wt1(MAX_FLT,MAXPARAM), actRate_Wt1(MAX_FLT,MAXPARAM), 
     1     wt_SR1(MAX_FLT,MAXPARAM), wt_RecInt1(MAX_FLT,MAXPARAM), MoRate_wt1(MAX_FLT,MAXPARAM), 
     1     magRecur_Wt1(MAX_FLT,MAXPARAM), faultThick_Wt1(MAX_FLT,MAXPARAM), 
     2     refMag_Wt1(MAX_FLT,MAX_Width, MAXPARAM), ftype_wt1(MAX_FLT,MAXPARAM)
      real seg_wt1(MAX_FLT,MAX_SEG)
      real wt_rateMethod1(MAX_FLT,4)
      
      integer j, kFlt, iThick, k, iFM, i
      integer iDip, iMagRec, iThick1, iSkip, iSR, iAct, iRecInt, iMo, ib
      integer iNode, iBR, jBR, nBR_SSC(MAX_FLT,MAX_NODE)
      real sum, totalSegWt(MAX_FLT)
      integer j1, iThick2

c     Set skip flag - if return 0, then there is no branch 
c     This is used for the activty rate methods that are not used
      iSkip = 1

c     Reset the dip weight
      if ( iNode .eq. 1 ) then
        do  iDip=1,n_Dip(kFlt)
          dip_wt1(kFlt,iDip) = 0.
        enddo
        dip_wt1(kFlt,iBR) = 1.
      endif

c     Reset the crustal thickness weight
      if ( iNode .eq. 2 ) then
        do  iThick=1,nThick(kFlt)
          faultThick_Wt1(kFlt,iThick) = 0.
        enddo
        faultThick_Wt1(kFlt,iBR) = 1.
      endif

c     Reset the ftype Model weight (This is the product of the ftmodelwt and the aleatory wt)
      if ( iNode .eq. 3 ) then
        k = 0
        do iFM=1,nFtypeModels(kFlt)
          do j=1,nFtype(kFlt,IFM)
            k = k + 1
            if ( iFM .eq. iBR ) then 
              ftype_wt1 (kFlt,k) =  ftype_al(kFlt,iFM,k)
            else
              ftype_wt1 (kFlt,k) =  0.
            endif
          enddo
        enddo
      endif

c     Reset the magpdf Model weight
      if ( iNode .eq. 4 ) then
        do iMagRec=1,nmagRecur(kFlt)
          magRecur_Wt1(kFlt,iMagRec) = 0.
        enddo
        magRecur_Wt1(kFlt,iBR) = 1. 
      endif

c     Reset the maxmag weight
c     combine all of the maxmag (by thickness) into one big branch
      if ( iNode .eq. 5 ) then
        do ithick1=1,nThick(kFlt)
          do j=1,nRefMag(kFlt,ithick1)
            refMag_Wt1(kFlt,iThick1,j) = 0.
          enddo
        enddo
        k = 0
        do ithick1=1,nThick(kFlt)
          do j=1,nRefMag(kFlt,ithick1)
            k = k + 1
            if ( k .eq. iBR) then
              refMag_Wt1(kFlt,iThick1,j) = 1.  
              j1 = j
              iThick2 = iThick1
              jBR = (iBR-1)/3 + 1
              faultThick_Wt1(kFlt,jBR) = 1.
            endif
          enddo
        enddo
      endif
      

c     Reset the rate type weight
      if ( iNode .eq. 6 ) then
        do i=1,nBR_SSC(kFlt,6)
          wt_RateMethod1(kFlt,i) = 0.
        enddo
        wt_RateMethod1(kFlt,iBR) = 1.
           
        if ( iBR .eq. 1 .and. nSR(kFlt) .eq. 0 ) iSkip = 0
        if ( iBR .eq. 2 .and. nActRate(kFlt) .eq. 0 ) iSkip = 0
        if ( iBR .eq. 3 .and. nRecInt(kFlt) .eq. 0 ) iSkip = 0
        if ( iBR .eq. 4 .and. nMoRate(kFlt) .eq. 0 ) iSkip = 0
      endif

c     Reset the SR weight
      if ( iNode .eq. 7 ) then
        do iSR=1,nSR(kFLt)
          wt_sr1(kFlt,iSR) = 0.
        enddo
        wt_SR1(kFlt,iBR) = 1.
      endif

c     Reset the a-value weight
      if ( iNode .eq. 8 ) then
        do iAct=1,nActRate(kFLt)
          actRate_wt1(kFlt,iAct) = 0.
        enddo
        actrate_Wt1(kFlt,iBR) = 1.
      endif

c     Reset the paleo weight
      if ( iNode .eq. 9 ) then
        do iRecInt=1,nRecInt(kFLt)
          wt_recInt1(kFlt,iRecInt) = 0.
        enddo
        wt_recInt1(kFlt,iBR) = 1.
      endif

c     Reset the Moment-rate weight
      if ( iNode .eq. 10 ) then
        do iMo=1,nMoRate(kFLt)
          MoRate_wt1(kFlt,iMo) = 0.
        enddo
        MoRate_wt1(kFlt,iBR) = 1.
      endif

c     Reset the b-value (for fault) weight
      if ( iNode .eq. 11 ) then
        do ib=1,n_bvalue(kFLt)
          bvalue_wt1(kFlt,ib) = 0.
        enddo
        bvalue_wt1(kFlt,iBR) = 1.
      endif

      return
      end
