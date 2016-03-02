      subroutine calcHaz ( haz, haz1, al_segwt, totalSegWt, 
     1     dip_Wt1, bValue_Wt1, actRate_Wt1, wt_SR1, wt_RecInt1, MoRate_wt1, magRecur_Wt1,    
     2     faultThick_Wt1, refMag_Wt1, ftype_wt1,
     3     wt_rateMethod1, 
     4     nInten, nFlt, n_Dip, n_bvalue, nRefMag,
     5     nFtype1, indexrate, nMagRecur, nRate, RateType, nThick, kflt, iPrint )

      implicit none
      include 'tornado.h'
      
c     Passed Variables
      real haz(MAX_INTEN, MAX_FLT, MAX_WIDTH, MAXPARAM,MAX_FTYPE) 
      real haz1(MAX_INTEN)
      real al_segwt(MAX_FLT), totalSegWt(MAX_FLT)
      integer nInten, nFlt, iPrint
      integer attentype(MAX_FLT)
      real dip_Wt1(MAX_FLT,MAXPARAM), bValue_Wt1(MAX_FLT,MAXPARAM), actRate_Wt1(MAX_FLT,MAXPARAM), 
     1     wt_SR1(MAX_FLT,MAXPARAM), wt_RecInt1(MAX_FLT,MAXPARAM), MoRate_wt1(MAX_FLT,MAXPARAM), magRecur_Wt1(MAX_FLT,MAXPARAM),    
     1     faultThick_Wt1(MAX_FLT,MAXPARAM), 
     2     refMag_Wt1(MAX_FLT,MAX_Width, MAXPARAM), ftype_wt1(MAX_FLT,MAXPARAM)
      real wt_rateMethod1(MAX_FLT,4)
      integer n_Dip(MAX_FLT),n_bvalue(MAX_FLT), nRefMag(MAX_FLT,MAX_WIDTH)
      integer nFtype1(MAX_FLT)
      integer indexrate(MAX_FLT,4), nMagRecur(MAX_FLT)
      integer nRate(MAX_FLT), RateType(MAX_Flt,MAXPARAM), nThick(MAX_FLT)


      integer iFlt, i1, iAtten, iInten, i, kflt
      real wt, wt1, sumwt
      integer iFtype, iParam, iThick, iDip
      integer iwidth, irate, iMagpdf, i_bvalue, iRefmag, jtype

c     initialize hazard
      do i=1,nInten
        haz1(i) = 0.
      enddo

      iFlt = kFlt

      
      if ( iPrint .eq. 1 ) then
        write (*,'( 10i5)') kFlt, nThick(iFlt), n_dip(iFlt),nMagRecur(iFlt), 
     1   nRate(iFlt), nFtype1(iFlt)
        pause 'test 1'
      endif

c     Compute hazard for the kflt only
c      do 900 iFlt = 1, nFlt

      sumwt = 0.

c      Loop over fault dips
       iWidth = 0
        do 860 iThick=1,nThick(iFlt)
           
c       Loop over crustal thickness
       do 855 iDip=1,n_Dip(iFlt)
         iWidth = iWidth + 1

c        Loop over parameter variations in same order as in out1 file
         iParam = 0
         do 850 imagpdf=1,nMagRecur(iFlt)
          do 840 iRate=1,nRate(iFlt)
 
c          Set wt for this rate 
c          set shift to get the to weight for the selected method
           i1 = iRate - indexRate(iFlt,rateType(iFlt,iRate))
           if ( rateType(iFlt,iRate) .eq. 1 ) then
             wt1 = wt_sr1(iFlt,i1) * wt_rateMethod1(iFlt,1)
             if (iprint .eq. 1 ) write (*,'( 3f10.4)') wt1 , wt_sr1(iFlt,i1) , wt_rateMethod1(iFlt,1)
            elseif ( rateType(iFlt,iRate) .eq. 2 ) then
             wt1 = actRate_Wt1(iFlt,i1) * wt_rateMethod1(iFlt,2)
           elseif ( rateType(iFlt,iRate) .eq. 3 ) then
             wt1 = wt_recInt1(iFlt,i1) * wt_rateMethod1(iFlt,3)
           elseif ( rateType(iFlt,iRate) .eq. 4 ) then
             wt1 = MoRate_wt1(iFlt,i1) * wt_rateMethod1(iFlt,4)
           endif  

           if ( n_bValue(iFlt) .eq. 0 ) then
             n_bvalue(iFlt) = 1
             bvalue_wt1(iFlt,1) = 1.
           endif

           do 830 i_bValue=1,n_bvalue(iFlt)
            do 820 iRefMag=1,nRefMag(iFlt,iThick)
             iParam = iParam + 1

             do 810 iFtype=1,nFtype1(iFlt)
c        write (*,'( 2i5, 2x,''nGM_Model'')') iFlt, nGM_Model(jtype)
c              do 805 iAtten=1,nGM_Model(jtype)

c              Set the weight for this set of parameters (epistemic)
               wt = wt1 * dip_Wt1(iFlt,iDip) * faultThick_wt1(iFlt,iThick) * bValue_Wt1(iFlt,i_bValue)
     2          * magRecur_Wt1(iFlt,imagpdf) * refMag_Wt1(iFlt,iThick,iRefMag) * ftype_wt1(iFlt,iFtype)
     3          * totalSegWt(iFlt)
               sumwt = sumwt + wt

       if (iPRint .eq. 1 ) then
       
                   write (*,'( 4i5,20f10.4)') iFLt, iWidth, iParam, iFtype, 
     1              wt1, dip_Wt1(iFlt,iDip) , faultThick_wt1(iFlt,iThick) , bValue_Wt1(iFlt,i_bValue)
     2              , magRecur_Wt1(iFlt,imagpdf) , refMag_Wt1(iFlt,iThick,iRefMag) , ftype_wt1(iFlt,iFtype)
     3              ,  totalSegWt(iFlt), wt
                   write (*,'( f10.5, 2x,''wt for branch'')') wt
                   pause
       endif
        
c              Add weight for aleatory rupture segmentation 
               wt = wt * al_segwt(iFlt) 

               if ( iPRint .eq. 1 ) write (65,'( 5i5, f10.4,10e12.4)') iAtten, iFlt, iWidth, iParam, iFtype, wt, 
     1               (haz(iInten, iFlt, iWidth, iParam, iFtype),iInten=1,nInten)

c           write (*, '( i5)') nInten
               do iInten=1,nInten
c              Add to hazard for this set of weights
                 haz1(iInten) =  haz1(iInten) + haz(iInten, iFlt, iWidth, iParam, iFtype) * wt
                 if ( iPrint .eq. 1 ) write (*,'( 3i5,3e12.3)') iFlt, iWidth, 
     1                iParam, haz(iInten, iFlt, iWidth, iParam, iFtype), wt, haz1(iInten)
               enddo
c               pause
c 805          continue

 810         continue
 820        continue
 830       continue
 840      continue
 850     continue
 855    continue
 860   continue
c       write (*,'( i5, f10.4,20e12.4)') iflt, sumwt, (haz1(iInten), iInten=1,nInten)
c       pause 'haz'
 900  continue

      return
      end
