
      subroutine calc_haz ( iFlt, nBR, haz_out5, wt_all, 
     1     al_segwt, wt_ftype, SR_Fact, faultFlag, haz1, flt0_flag,
     2     nSegMOdel, Wt_segModel, nFtype, jFlt, jNode)

      implicit none
      include 'tornado.h'

      real Haz_out5(MAX_FLT,MAX_WIDTH,MAXPARAM,MAX_FTYPE)
      integer iBR1, iBR2, iBR3, iBR4, IBR5, IBR6, IBR7, IBR8,
     1        iBR9, iBR10, iBR11, iSeg, iFlt
      integer ithick1, k
      integer  nBR(MAX_FLT,12,MAX_THICK), flt0_flag(MAX_FLT)
      integer nSegModel(MAX_FLT)
      real Wt_segModel(MAX_FLT,MAX_SEG)
      real wt_all(MAX_FLT,12,MAX_THICK,MAX_BR)
      real haz1
      real*8 sum, wt_sum
      real SR_Fact(MAX_FLT1,MAX_SEG,MAX_FLT)
      real al_segwt(MAX_FLT)       
      real wt_ftype(MAX_FLT,5,5), wt
      integer faultflag(MAX_FLT1,MAX_SEG,MAX_FLT)
      integer iFlt0, ithick, n7, n8, n9, n10, n11, i2, jFlt
      integer iRate, iParam, nRate, iWidth, iFM, isum, iFtypeIndex, iFtype
      integer nFtype(MAX_FLT,MAXPARAM), jj, jNode, iPrint

c      SSC branches: 1=crustal thick, 2 = Dip, , 3=ftype model, 4=mag pdf, 5=refmag,
c                   6=Rate Method, 7=slip rate for faults, 8=a-value&b-value pair, 
c                   9=paleo recurrence interval, 
c                   10=Moment rate, 11=b_value flt for fault, 12=segModel

c     Initialize hazard sum
      sum = 0.
      wt_sum = 0.
      
      iPrint = 0
      if ( jNode .eq. 7 .and. iflt .eq. 193 ) then
        iPrint = 2
      endif
            
      iFlt0 = flt0_flag(iFlt)
      
      if ( iPrint .eq. 1 ) then
        write (*,'( 3i5)') iFlt0, iFlt, nSegModel(iFlt0)
        write (*,'( 10f10.3)') (Wt_segModel(iflt0,iSeg),iSeg=1,nSegModel(iFlt0))
        pause 'nSegModel'
      endif

      do 100 iSeg=1,nSegModel(iFlt0)
      if ( iPrint .eq. 1 ) write (*,'( 2x,''iseg'',i5,f10.3,i5)') iSeg, Wt_segModel(iflt0,iSeg),
     1  faultflag(iFlt0,iSeg,jFlt)
        

c      only include this fault if it is in the segmentation model
       if ( faultflag(iFlt0,iSeg,jFlt) .eq. 0 ) goto 100     
      
c      Loop over all nodes for this fault
c      write (*,'( i5,2x,''nBR1'')') nBR(iflt,1,1)
       do 101 iBR1=1,nBR(iflt,1,1)
       ithick=iBR1
       do 102 iBR2=1,nBR(iflt,2,ithick)
        do 103 iBR3=1,nBR(iflt,3,ithick)
         do 104 iBR4=1,nBR(iflt,4,ithick)
          do 105 iBR5=1,nBR(iflt,5,ithick)
           do 106 iBR6=1,nBR(iflt,6,ithick)
            if ( wt_all(iFlt,6,ithick,iBR6) .eq. 0. ) goto 106
            if ( nBR(iflt,7,ithick) .eq. 0 ) then
              n7 = 1
              wt_all(iFlt,7,ithick,1) = 1.0
            else
              n7 = nBR(iflt,7,ithick)
            endif

            do 107 iBR7=1,n7
            
             if ( nBR(iflt,8,ithick) .eq. 0 ) then
              n8 = 1
              wt_all(iFlt,8,ithick,1) = 1.0
             else
              n8 = nBR(iflt,8,ithick)
             endif

             do 108 iBR8=1,n8
             
              if ( nBR(iflt,9,ithick) .eq. 0 ) then
               n9 = 1
               wt_all(iFlt,9,ithick,1) = 1.0
              else
               n9 = nBR(iflt,9,ithick)
              endif
             
              do 109 iBR9=1,n9

               if ( nBR(iflt,10,ithick) .eq. 0 ) then
                n10 = 1
                wt_all(iFlt,10,ithick,1) = 1.0
               else
                n10 = nBR(iflt,10,ithick)
               endif

               do 110 iBR10=1,n10

                if ( nBR(iflt,11,ithick) .eq. 0 ) then
                 n11 = 1
                 wt_all(iFlt,11,ithick,1) = 1.0
                else
                 n11 = nBR(iflt,11,ithick)
                endif

                do 111 iBR11=1,n11
                
c                write (*,'( 15i5)') iflt, iSeg, iBR1, iBR2, iBR3, iBR4, 
c     1          iBR5, iBR6, iBR7, iBR8, iBR9, iBR10, iBR11
                      
c                 Now, convert the branch indexes to the iParam used in the hazard code  
c                 First, find the index (iRate) for the method used for the rate     
                  if (iBR6 .eq. 1 ) then
                    iRate = iBR7
                  elseif (iBR6 .eq. 2 ) then
                    iRate = nBR(iFlt,7,ithick) + iBR8
                  elseif (iBR6 .eq. 3 ) then
                    iRate = nBR(iFlt,7,ithick)+ nBR(iFlt,8,ithick) + iBR9
                  elseif (iBR6 .eq. 4 ) then
                    iRate = nBR(iFlt,7,ithick)+ nBR(iFlt,8,ithick)
     1                + nBR(iFlt,9,ithick) + iBR10
                  endif
                  nRate = nBR(iFlt,7,ithick)+ nBR(iFlt,8,ithick)
     1                + nBR(iFlt,9,ithick) + nBR(iFlt,9,ithick) 
          
c                 set iParam = (iMagRecur-1)* Nrate * n_b * nRefMag
c                       + (iRate-1) * N_b * NrefMag 
c                       + (ib-1) * nRefMag 
c                       + iRefMag
c                 This is how the hazard code loads up the iParam cases into the hazard
                  if ( nBR(iflt,11,ithick) .ge. 1 ) then
                    iParam = (iBR4-1) * nRate * nBR(iflt,11,iThick)
     1                   * nBR(iflt,5,iThick)
     1             + (iRate-1) *  nBR(iflt,11,ithick) * nBR(iflt,5,ithick)
     2             + (iBR11-1) * nBR(iflt,5,ithick)
     3             + iBR5 
                  else
                    iParam = (iBR4-1) * nRate * nBR(iflt,5,iThick)
     1             + (iRate-1) * nBR(iflt,5,ithick)
     3             + iBR5 
                  endif

c                 Set width index     
c                 iWidth = (ithick-1) * nDip + iDip
                  iWidth = (ithick-1) *nBR(iflt,2,ithick) + iBR2
              
c                 Rename ftype model selected index
                  iFM = iBR3            

c                 Find the postion of this ftype model in the global array of ftypes
c                 iSum is the number of ftype cases to skip for this selected iFM
                  isum = 0
                  do i2=1,iFM-1
                    isum = isum + nFtype(iflt,i2)
                  enddo
              
c              write (*,'( i5)') nFtype(iflt,iFM)
c              pause 'test 10'
              
c                 sum over aleatory ftypes within this ftype model
                  do iFtype=1, nFtype(iflt,iFM)
                    iFtypeIndex = isum + iFtype
           
                    wt = wt_all(iFlt,1,ithick,iBR1) *
     1                   wt_all(iFlt,2,ithick,iBR2) *
     1                   wt_all(iFlt,3,ithick,iBR3) *
     1                   wt_all(iFlt,4,ithick,iBR4) *
     1                   wt_all(iFlt,5,ithick,iBR5) *
     1                   wt_all(iFlt,6,ithick,iBR6) *
     1                   wt_all(iFlt,7,ithick,iBR7) *
     1                   wt_all(iFlt,8,ithick,iBR8) *
     1                   wt_all(iFlt,9,ithick,iBR9) *
     1                   wt_all(iFlt,10,ithick,iBR10) *
     1                   wt_all(iFlt,11,ithick,iBR11) 

c                  write (*,'( 2x,''iParam, iWidth='',4i5,e12.4)')iFlt,iWidth,iParam,iFtypeIndex,
c     1             Haz_out5(iFlt,iWidth,iParam,iFtypeIndex)
                    
c                   Add the hazard curve to the total hazard for this case
                    sum = sum + Haz_out5(iFlt,iWidth,iParam,iFtypeIndex)
     1                    * al_segwt(iflt) * wt_ftype(iflt,iFM,iFtype)
     1                    * SR_Fact(iFlt0,iSeg,iflt) *wt_segModel(iFlt0,iSeg) * wt
                    wt_sum = wt_sum + wt_ftype(iflt,iFM,iFtype)
     1                     *wt_segModel(iFlt0,iSeg) * wt

           if (iprint .eq. 1 ) then               
                     write (*,'(i5,2e12.4,10f8.3 )') iSeg, Haz_out5(iFlt,iWidth,iParam,iFtypeIndex), sum,
     1 al_segwt(iflt) , wt_ftype(iflt,iFM,iFtype)
     1                    , SR_Fact(iFlt0,iSeg,iflt) , wt_segModel(iFlt0,iSeg) , wt
       write (*,'( 12i5)') iBR1, iBR2, iBR3, iBR4, iBR5, iBR6, iBR7, iBR8, iBR9, iBR10, iBR11, iSeg 
       write (*,'( 5i5)') iFlt, iWidth, iParam, iFtypeIndex 
       
                     write (*,'( 15f10.4)') wt_all(iFlt,1,ithick,iBR1) ,
     1                   wt_all(iFlt,2,ithick,iBR2) ,
     1                   wt_all(iFlt,3,ithick,iBR3) ,
     1                   wt_all(iFlt,4,ithick,iBR4) ,
     1                   wt_all(iFlt,5,ithick,iBR5) ,
     1                   wt_all(iFlt,6,ithick,iBR6) ,
     1                   wt_all(iFlt,7,ithick,iBR7) ,
     1                   wt_all(iFlt,8,ithick,iBR8) ,
     1                   wt_all(iFlt,9,ithick,iBR9) , 
     1                   wt_all(iFlt,10,ithick,iBR10),
     1                   wt_all(iFlt,11,ithick,iBR11) 
             write (*,'( 2i5,10f10.3)') ithick, iBR5, (wt_all(iFlt,1,ithick,k),k=1,nBR(iflt,1,1))
             write (*,'( 2x,''check thick wts'')')
          endif

               
                  enddo
 111            continue
 110           continue
 109          continue
 108         continue
 107        continue
 106       continue
 105      continue
 104     continue
 103    continue
 102   continue
 101  continue
 100  continue
      if (iprint .eq. 1) then
        write (*,'( i5,e12.3,f10.4)') jNode, sum, wt_sum
        pause 'jNode, haz, wt_sum'
      endif
      
      haz1 = sum
      
         
      return
      end


         
