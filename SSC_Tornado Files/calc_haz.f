      subroutine calcHaz ( haz, haz1, al_segwt, totalSegWt, 
     1     dip_Wt1, bValue_Wt1, actRate_Wt1, wt_SR1, wt_RecInt1, MoRate_wt1, magRecur_Wt1,    
     2     faultThick_Wt1, refMag_Wt1, ftype_wt1,
     3     wt_rateMethod1, gm_wt1, 
     4     nInten, nFlt, attentype, nGM_model, n_Dip, n_bvalue, nRefMag,
     5     nFtype1, indexrate, nMagRecur, nRate, RateType, nThick, kflt, iPrint )

      implicit none
      include 'tornado.h'
      
c     Passed Variables
      real haz(MAX_INTEN,MAX_ATTEN, MAX_FLT, MAX_WIDTH, MAXPARAM,MAX_FTYPE) 
      real haz1(MAX_INTEN)
      real al_segwt(MAX_FLT), totalSegWt(MAX_FLT)
      integer nInten, nFlt, iPrint
      integer attentype(MAX_FLT), nGM_model(MAX_ATTENTYPE)
      real dip_Wt1(MAX_FLT,MAXPARAM), bValue_Wt1(MAX_FLT,MAXPARAM), actRate_Wt1(MAX_FLT,MAXPARAM), 
     1     wt_SR1(MAX_FLT,MAXPARAM), wt_RecInt1(MAX_FLT,MAXPARAM), MoRate_wt1(MAX_FLT,MAXPARAM), magRecur_Wt1(MAX_FLT,MAXPARAM),    
     1     faultThick_Wt1(MAX_FLT,MAXPARAM), 
     2     refMag_Wt1(MAX_FLT,MAX_Width, MAXPARAM), ftype_wt1(MAX_FLT,MAXPARAM)
      real wt_rateMethod1(MAX_FLT,4)
      integer n_Dip(MAX_FLT),n_bvalue(MAX_FLT), nRefMag(MAX_FLT,MAX_WIDTH)
      integer nFtype1(MAX_FLT)
      integer indexrate(MAX_FLT,4), nMagRecur(MAX_FLT)
      integer nRate(MAX_FLT), RateType(MAX_Flt,MAXPARAM), nThick(MAX_FLT)
      real gm_wt1(4,MAX_ATTEN)

      integer iFlt, i1, iAtten, iInten, i, kflt
      real wt, wt1, sumwt
      integer iFtype, iParam, iThick, iDip
      integer iwidth, irate, iMagpdf, i_bvalue, iRefmag, jtype

c     initialize hazard
      do i=1,nInten
        haz1(i) = 0.
      enddo

c     Compute hazard for the kflt only
      iFlt = kFlt
c      do 900 iFlt = 1, nFlt

      sumwt = 0.

c      Loop over fault dips
       iWidth = 0
c       write (*,'( 2i5, 2x,''nDip'')') iFlt, n_dip(iFlt)
        do 860 iThick=1,nThick(iFlt)
           
c       Loop over crustal thickness
c       write (*,'( 2i5, 2x,''nThick'')') iFlt, nThick(iFlt)
       do 855 iDip=1,n_Dip(iFlt)
         iWidth = iWidth + 1

c        Loop over parameter variations in same order as in out1 file
         iParam = 0
c       write (*,'( 2i5, 2x,''nMagRecur'')') iFlt, nMagRecur(iFlt)
         do 850 imagpdf=1,nMagRecur(iFlt)
c       write (*,'( 2i5, 2x,''nRate'')') iFlt, nRate(iFlt)
          do 840 iRate=1,nRate(iFlt)
 
c          Set wt for this rate 
c          set shift to get the to weight for the selected method
           i1 = iRate - indexRate(iFlt,rateType(iFlt,iRate))
           if ( rateType(iFlt,iRate) .eq. 1 ) then
             wt1 = wt_sr1(iFlt,i1) * wt_rateMethod1(iFlt,1)
            elseif ( rateType(iFlt,iRate) .eq. 2 ) then
             wt1 = actRate_Wt1(iFlt,i1) * wt_rateMethod1(iFlt,2)
           elseif ( rateType(iFlt,iRate) .eq. 3 ) then
             wt1 = wt_recInt1(iFlt,i1) * wt_rateMethod1(iFlt,3)
           elseif ( rateType(iFlt,iRate) .eq. 4 ) then
             wt1 = MoRate_wt1(iFlt,i1) * wt_rateMethod1(iFlt,4)
           endif  

c       write (*,'( 2i5, 2x,''n_bvalue'')') iFlt, n_bvalue(iFlt)

           do 830 i_bValue=1,n_bvalue(iFlt)
c       write (*,'( 2i5, 2x,''nRefMag'')') iFlt, nRefMag(iFlt,ithick)
            do 820 iRefMag=1,nRefMag(iFlt,iThick)
             iParam = iParam + 1

c       write (*,'( 2i5, 2x,''nFtype1'')') iFlt, nFtype1(iFlt)

             do 810 iFtype=1,nFtype1(iFlt)
              jType = attenType(iFlt)
c       write (*,'( 2i5, 2x,''nGM_Model'')') iFlt, nGM_Model(jtype)
              do 805 iAtten=1,nGM_Model(jtype)

c              Set the weight for this set of parameters (epistemic)
               wt = wt1 * dip_Wt1(iFlt,iDip) * faultThick_wt1(iFlt,iThick) * bValue_Wt1(iFlt,i_bValue)
     2          * magRecur_Wt1(iFlt,imagpdf) * refMag_Wt1(iFlt,iThick,iRefMag) * ftype_wt1(iFlt,iFtype)
     3          * gm_wt1(jType,iAtten) * totalSegWt(iFlt)
               sumwt = sumwt + wt

c                   if ( iflt .eq. 3 ) write (*,'( 5i5,20e12.4)') iFLt, iAtten, iWidth, 
c     1                 iParam, iFtype, haz(1,iAtten, iFlt, iWidth, iParam, iFtype)
c                   write (66,'( 9i5,20f10.4)') iFLt, iAtten, iWidth, iParam, iFtype, 
c     1              wt1, dip_Wt1(iFlt,iDip) , faultThick_wt1(iFlt,iThick) , bValue_Wt1(iFlt,i_bValue)
c     2              , magRecur_Wt1(iFlt,imagpdf) , refMag_Wt1(iFlt,iThick,iRefMag) , ftype_wt1(iFlt,iFtype)
c     3              , gm_wt1(jType,iAtten) ,  totalSegWt(iFlt), wt
c                   write (*,'( f10.5, 2x,''wt for branch'')') wt
c                   pause
        
c              Add weight for aleatory rupture segmentation 
               wt = wt * al_segwt(iFlt) 

               if ( iPRint .eq. 1 ) write (65,'( 5i5, f10.4,10e12.4)') iAtten, iFlt, iWidth, iParam, iFtype, wt, 
     1               (haz(iInten,iAtten, iFlt, iWidth, iParam, iFtype),iInten=1,nInten)

               do iInten=1,nInten
c              Add to hazard for this set of weights
                 haz1(iInten) =  haz1(iInten) + haz(iInten,iAtten, iFlt, iWidth, iParam, iFtype) * wt
c        if ( iflt .eq. 3 .and. iInten .eq. 1 )  then
c            write (*,'( 3e12.3)') haz1(iInten) , haz(iInten,iAtten, iFlt, iWidth, iParam, iFtype) , wt
c         endif
               enddo
 805          continue

 810         continue
 820        continue
 830       continue
 840      continue
 850     continue
 855    continue
 860   continue
c       write (*,'( i5, f10.4,e12.4)') iflt, sumwt, haz1(1)
 900  continue
c      pause

      return
      end
