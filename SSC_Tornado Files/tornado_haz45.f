
      program Tornado_Haz45

c     Last modified:  8/15  

      implicit none
      include 'tornado.h'
      
      real haz(MAX_INTEN, MAX_FLT, MAX_WIDTH, MAXPARAM,MAX_FTYPE) 
      real haz1(MAX_INTEN)
      real testInten(MAX_INTEN)
      real al_segwt(MAX_FLT)
      integer nInten, jcalc(MAX_ATTENTYPE,MAX_ATTEN), nFlt, isite, iflt
      integer iInten, iBR
      character*80 filein, file1, fname(MAX_FLT)
      integer attentype(MAX_FLT), nGM_model(MAX_ATTENTYPE), iFlag, iPrint

      integer nFtype1(MAX_FLT)
      real meanhaz(MAX_INTEN) , nonPoisson(3,MAX_INTEN), 
     1     meanHaz1(MAX_FLT,MAX_INTEN), meanHaz2(MAX_INTEN)
      integer nProb      , nSite, jBr
      integer indexrate(MAX_FLT,4)
      real xFact(100), xx
      integer nBR, iSkip
 
      integer sssCalc(MAX_ATTENTYPE,MAX_ATTEN), scalc(MAX_ATTENTYPE,MAX_ATTEN)
      real sigFix(MAX_ATTENTYPE,MAX_ATTEN)

      real*8 haz_SSC(MAX_FLT, MAX_NODE, MAX_BR, MAX_INTEN) 
      real*8 haz_SSC1(MAX_NODE, MAX_BR, MAX_INTEN) 
      real sum, totalSegWt(MAX_FLT), totalSegWt_all(MAX_FLT)
      real*8 sum1

      integer n_Dip(MAX_FLT),n_bvalue(MAX_FLT), nActRate(MAX_FLT), nSR(MAX_FLT), 
     1        nRecInt(MAX_FLT), nMoRate(MAX_FLT),
     1        nRefMag(MAX_FLT,MAX_WIDTH), nFtypeModels(MAX_FLT)
      integer nFtype(MAX_FLT,MAXPARAM)
      real ftype_al(MAX_FLt,MAXPARAM,5)
      real t1, fact1

      real segwt(MAX_FLT,MAX_FLT), dipWt(MAX_FLT,MAXPARAM), bValueWt(MAX_FLT,MAXPARAM), actRateWt(MAX_FLT,MAXPARAM), 
     1     wt_SR(MAX_FLT,MAXPARAM), wt_RecInt(MAX_FLT,MAXPARAM), wt_MoRate(MAX_FLT,MAXPARAM), magRecurWt(MAX_FLT,MAXPARAM),    
     1     faultThickWt(MAX_FLT,MAXPARAM), 
     2     refMagWt(MAX_FLT,MAX_Width, MAXPARAM), ftmodelwt(MAX_FLT,MAXPARAM)
      real wt_rateMethod(MAX_FLT,4)
      integer iPer, nFlt0, f_start(MAX_FLT), f_num(MAX_FLT)
      integer nMagRecur(MAX_FLT)
      real contrib_min
      integer iNode, i,  jFlt, kFlt, nNode_SSC
      integer nBR_SSC(MAX_FLT,MAX_NODE), nRate(MAX_FLT), RateType(MAX_FLt,MAXPARAM)
      integer nBR_SSC1(MAX_FLT,MAX_NODE)
      integer j, k, iThick, nThick(MAX_FLT), iFM, iDip, imagRec, iThick1
      integer iSR, iAct, iRecInt, iMo, ib, j2, iSum
      real seg_wt1(MAX_FLT,MAX_SEG)
      integer nSegModel(MAX_FLT), segFlag(MAX_FLT, MAX_SEG)
      real ratio1, ratio2
      integer iFtype, iWidth, node_12_flag(MAX_FLT), nSegModel1(MAX_FLT)
      real period

      real dip_Wt1(MAX_FLT,MAXPARAM), bValue_Wt1(MAX_FLT,MAXPARAM), actRate_Wt1(MAX_FLT,MAXPARAM), 
     1     wt_SR1(MAX_FLT,MAXPARAM), wt_RecInt1(MAX_FLT,MAXPARAM), MoRate_wt1(MAX_FLT,MAXPARAM), 
     1     magRecur_Wt1(MAX_FLT,MAXPARAM), faultThick_Wt1(MAX_FLT,MAXPARAM), 
     2     refMag_Wt1(MAX_FLT,MAX_Width, MAXPARAM), ftype_wt1(MAX_FLT,MAXPARAM)
      real wt_rateMethod1(MAX_FLT,4), hazLevel, GM_ratio(100), GM0, GM1, GM2(MAX_BR)
      integer iFlt0, fltSys_flag(MAX_FLT), iflt1, iflt2
    
      write (*,*) '*************************'
      write (*,*) '* Tornado Code for use *'
      write (*,*) '*  with Hazard_v45 *'
      write (*,*) '*    Feb, 2016, NAA     *'
      write (*,*) '*************************'

      write (*,*) 'Enter the input filename.'
      
      read (*,'(a80)') filein
      open (31,file=filein,status='old')

      read (31,*) iPer
      read (31,*) contrib_min
      read (31,*) Hazlevel

c     Read Input File
      call RdInput ( nInten,  testInten, iPer, nProb, period )
c      pause 'out of rDinput'

c     Read run file
      call Rd_Fault_Data  (nFlt, nFlt0, f_start, f_num, AttenType, 
     1           n_Dip, n_bValue, nActRate,  nSR,   nRecInt,   nMoRate,   nMagRecur,  nThick,
     1           nRefMag,  nFtypeModels,
     2           dipWt, bValueWt, actRateWt, wt_sr, wt_recInt, wt_MoRate, magRecurWt, 
     3           faultThickWt, refMagWt, ftmodelwt,
     3           nFtype, ftype_al, wt_rateMethod, al_Segwt,
     3           nRate, rateType, nBR_SSC, nSegModel, segwt, segFlag, indexRate, fname )
      write (*,'( i5)') nFlt

c     sources are indexed in two ways
c     When the out5 file is read, the sources are just numbered 1 to NFLT
c     For segmentation, the source are number 1 to nFLT0 for the independent groups of sources
c     THe sources included in each group are inducated f_start and f_num
c        f_start(k) gives the index (1 to NFLT) of the first source in group k
c        f_num(k) is the number of sources in group k

c     Set flag to indicate if this fault is the first in the segmentation  model 
      do iFlt=1,nFlt
        node_12_flag(iFlt) = 0
      enddo  
      do iFlt0=1,nFlt0
        iFlt = f_start(iFlt0)
        node_12_flag(iFlt) = 1 
        nSegModel1(iFlt) = nSegModel(iFlt0)
        fltsys_flag(iFlt) = iFlt0
      enddo   

c     Set the total number of nodes in the SSC
      nNode_SSC = 12

c     Set nBR_SSC1 (for output)
      do iFlt=1,nFlt
       do iNode=1,nNode_SSC
        nBR_SSC1(iFlt,iNode) = nBR_SSC(iFlt,iNode)
        if ( iNode .eq. 6 ) then
          iSum = 0
          do iBR=1,4
           if ( iBR .eq. 1 .and. nSR(iFlt) .gt. 0 ) iSum = iSum + 1
           if ( iBR .eq. 2 .and. nActRate(iFlt) .gt. 0 ) iSum = iSum + 1
           if ( iBR .eq. 3 .and. nRecInt(iFlt) .gt. 0 ) iSum = iSum + 1
           if ( iBR .eq. 4 .and. nMoRate(iFlt) .gt. 0 ) iSum = iSum + 1
          enddo
          nBR_SSC1(iFlt,iNode) = iSum
        endif 
       enddo
      enddo  

      read (31,'( a80)') file1
      write (*,'( a80)') file1
      open (43,file=file1,status='unknown')
      write (43,'( 2x,''period index, period: '',i5, f10.4)') iPer, period
      write (43,'( 2x,''Min contribution to GM uncertainty: '',f10.5)') contrib_min
      write (43,'( 2x,''Hazard level: '',e12.4)') Hazlevel
      write (43,'( 2x,''Site, Flt, Node, ln(max fact), GM ratios...'')')

      read (31,*) nSite

c     Loop Over Number of  sites (iSite is a summary used in haz runs)
      do 1000 iSite = 1,nSite

c      Read the out1 file       
       write (*,'( 2x,''reading logic tree file'')')
       call read_out5 ( nFlt, haz, nFtype1, iPer, nProb )

c      THis is set up using brute force.  It is not efficient, but will be easier to modify to 
c      account for correlations later.

c      First, reset the SSC Weights to starting value for all faults
       do jFlt=1,nFlt           
         call Set_starting_wts ( n_dip, dipWt, dip_wt1, 
     1     nThick, faultThickWt, faultThick_Wt1, nRefMag, refMagWt, refMag_Wt1,
     2     n_bValue, bvalueWt, bvalue_Wt1, nMagRecur, magRecurWt, magRecur_Wt1,
     3     nFtypeModels, nFtype, ftModelwt, ftype_al, ftype_wt1,
     4     nSegModel, segwt, seg_wt1, wt_RateMethod, wt_RateMethod1,
     5     nSR, wt_SR, wt_SR1, nActRate, actRateWt, actRate_Wt1,
     6     nRecInt, wt_recInt, wt_recInt1, nMoRate, wt_MoRate, MoRate_wt1,
     8     segFlag, totalSegWt, totalSegWt_all, jFlt )
       enddo

c      Compute the mean hazard for each fault separately
       do kFlt=1,nFlt
         iPrint = 0
         call calcHaz ( haz, haz1, al_segwt, totalSegWt, 
     1         dip_Wt1, bValue_Wt1, actRate_Wt1, wt_SR1, wt_RecInt1, MoRate_wt1, 
     2         magRecur_Wt1, faultThick_Wt1, refMag_Wt1, ftype_wt1,
     3         wt_rateMethod1,  
     4         nInten, nFlt, n_Dip, n_bvalue, nRefMag,
     5         nFtype1, indexrate, nMagRecur, nRate, RateType, nThick, kflt, iPrint )
        
        write (*,'( i5,10e12.4)') kflt, (haz1(i),i=1,nInten)
c        Save the mean hazard for this flt to new array (meanhaz1)
         do i=1,nInten
           meanhaz1(kFlt,i)= haz1(i)
         enddo
       enddo
        pause

c      Compute the total mean hazard for all faults (meanhaz)
       do i=1,nInten
         meanhaz(i)= 0.
       enddo
       do iFlt=1,nFlt
         do i=1,nInten
           meanhaz(i)= meanhaz1(iFlt,i) + meanhaz(i)
         enddo
       enddo

c      SSC, reset the weights to unity for a given branch, fault and one node at a time
       do 950 kFlt = 1, nFlt
c         write (*,'( 6i5)') iSite, nSite, kFlt, nFlt

        do 940 iNode=1,nNode_SSC 
c        Set special handling for segmentation (Node=12)
         if (iNode .eq. 12 ) goto 935
         
c         write (*,'( 3i5)') kflt, iNode, nBR_SSC(kFlt,iNode)
         do 930, iBR=1,nBR_SSC(kFlt,iNode)

c          Remove the mean hazard from this fault from total mean 
c          The mean haz w/o this fault is called meanhaz2
c          We will add the hazard from this fault back in later after changing the weights
           do iInten=1,nInten
             meanhaz2(iInten) = meanhaz(iInten) - meanhaz1(kFlt,iInten)
           enddo

c         Now, set weights to unity for the Node and Branch of interest for fault kFlt
c         SSC branches: 1 = Dip, 2=crustal thick, 3=ftype, 4=magpdf, 5=maxmag,  6=RateType, 7=SR, 8=a-value, 9=paleo, 
c                       10=Moment, 11=b_value flt, 12=segModel

         call Set_new_wts ( n_dip, dipWt, dip_wt1, 
     1     nThick, faultThickWt, faultThick_Wt1, nRefMag, refMagWt,refMag_Wt1,
     2     n_bValue, bvalueWt, bvalue_Wt1, nMagRecur, magRecurWt, magRecur_Wt1,
     3     nFtypeModels, nFtype, ftModelwt, ftype_al, ftype_wt1,
     4     nSegModel, segwt, seg_wt1, wt_RateMethod, wt_RateMethod1,
     5     nSR, wt_SR, wt_SR1, nActRate, actRateWt, actRate_Wt1,
     6     nRecInt, wt_recInt, wt_recInt1, nMoRate, wt_MoRate, MoRate_wt1,
     7     segFlag, totalSegWt, totalSegWt_all, kFlt, iNode, iBR, nBR_SSC, iSkip )

c        Check if there is no branch for this node (activity method not used)
         if ( iSkip .eq. 0 ) goto 925

c        Compute the hazard for this set of weights
         iPrint = 0
         call calcHaz ( haz, haz1, al_segwt, totalSegWt, 
     1         dip_Wt1, bValue_Wt1, actRate_Wt1, wt_SR1, wt_RecInt1, MoRate_wt1, magRecur_Wt1,    
     2         faultThick_Wt1, refMag_Wt1, ftype_wt1,
     3          wt_rateMethod1, 
     4         nInten, nFlt, n_Dip, n_bvalue, nRefMag,
     5         nFtype1, indexrate, nMagRecur, nRate, RateType, nThick, kflt, iPrint )

c        Add the hazard from this fault (for this branch) to total hazard array
         do iInten=1,nInten
            haz_SSC(kFlt,iNode,iBR,iInten) = meanhaz2(iInten) + haz1(iInten)
         enddo
         if ( iNode .eq. 7 ) then
           if ( kFlt .ge. 31 .and. kFlt .le. 33 ) write (*,'( 3i7,10f10.3)') kFLt, iNode, iBR, 
     1    (haz_SSC(kFlt,iNode,iBR,i)/meanhaz(i),i=1,nInten)
         endif

 925    continue
         
C       Reset the weights for this fault back to the starting weights
        
         call Set_starting_wts ( n_dip, dipWt, dip_wt1, 
     1     nThick, faultThickWt, faultThick_Wt1, nRefMag, refMagWt,refMag_Wt1,
     2     n_bValue, bvalueWt, bvalue_Wt1, nMagRecur, magRecurWt, magRecur_Wt1,
     3     nFtypeModels, nFtype, ftModelwt, ftype_al, ftype_wt1,
     4     nSegModel, segwt, seg_wt1, wt_RateMethod, wt_RateMethod1,
     5     nSR, wt_SR, wt_SR1, nActRate, actRateWt, actRate_Wt1,
     6     nRecInt, wt_recInt, wt_recInt1, nMoRate, wt_MoRate, MoRate_wt1,
     7     segFlag, totalSegWt, totalSegWt_all, kFlt )

 930    continue
c       End loop on iBR 

c       Skip past special coding for node 12
        goto 940

 935    continue
c       Special coding for iNode=12 (segmentation model)

c       Set the number of branches for the segmenation node
c       If 1 branch, then set to mean hazard
        if ( node_12_flag(kFlt) .eq. 0 ) then
          nBR_SSC(kFlt,12) = 1
          nBR_SSC1(kFlt,12) = 1
          do iInten=1,nInten
            haz_SSC(kFlt,iNode,iBR,iInten) = meanhaz(iInten) 
          enddo
          goto 940
        else
          nBR_SSC(kFlt,12) = nSegModel1(kFlt)
          nBR_SSC1(kFlt,12) = nSegModel1(kFlt)
        endif

c       Remove the hazard from faults in this fault group from the total haz
        do i=1,nInten
           meanhaz2(i)= meanhaz(i) 
        enddo
        iflt1 = f_start(fltsys_flag(kFlt))
        iflt2 = f_num(fltsys_flag(kFlt))+iFlt1-1
        do iFlt=iFlt1,iFlt2
          do i=1,nInten
           meanhaz2(i)= meanhaz2(i) - meanhaz1(iFlt,i)
          enddo
        enddo
c        write (*,'( 12e12.4)') (meanhaz2(i),i=1,nInten)
c        pause 'menhaz2'

c       Loop over segmentation branches
        do iBR=1,nBR_SSC(kFlt,12)

c         Reset the weights for this branch (one to unity)
          do i=1,nSegModel(kFlt)
            seg_Wt1(kFlt,i) = 0.
          enddo
          seg_Wt1(kFlt,iBR) = 1.

c         Find the total weight for each segment in this set
          do jFlt=iFlt1,iFlt2
            sum = 0.
            do i=1,nSegModel(kFlt)
              sum = sum + seg_wt1(kFlt,i) * segFlag(jFlt,i)
            enddo
            totalSegWt(jFlt) = sum
c            write (*,'( i5,f10.4)') jFlt, totalSegWt(jFlt)
          enddo

c         Initialize haz_SSC to total w/o this group of faults
          do iInten=1,nInten
            haz_SSC(kFlt,iNode,iBR,iInten) = meanhaz2(iInten) 
          enddo

c         Add the hazard back to the total for this group of faults
          do jFlt=iFlt1,iFlt2
c           Compute the hazard for this set of weights
            iPrint = 0
            call calcHaz ( haz, haz1, al_segwt, totalSegWt, 
     1         dip_Wt1, bValue_Wt1, actRate_Wt1, wt_SR1, wt_RecInt1, MoRate_wt1, magRecur_Wt1,    
     2         faultThick_Wt1, refMag_Wt1, ftype_wt1,
     3         wt_rateMethod1, 
     4         nInten, nFlt, n_Dip, n_bvalue, nRefMag,
     5         nFtype1, indexrate, nMagRecur, nRate, RateType, nThick, jflt, iPrint )

c           Add the hazard from these fault (for this branch) to total hazard array
            do iInten=1,nInten
              haz_SSC(kFlt,iNode,iBR,iInten) = haz_SSC(kFlt,iNode,iBR,iInten) + haz1(iInten)
            enddo
          enddo
c          write (*,'( 5i5)') kFlt,iNode,iBR, iFlt1, iFlt2
c         write (*,'( 3i3, 12e12.4)') kFlt,iNode,iBR,(haz_SSC(kFlt,iNode,iBR,i),i=1,nInten)
         
        enddo
c         pause 'haz_SSC'
 
 940    continue
c       End loop over Nodes 

 950   continue         
c      End loop over faults


c       Write out sensitivity hazard curves for SSC
        read (31,'( a80)') file1
        open (42,file=file1,status='unknown')
        write (42,'( ''SSC Nodes: 1 = Dip, 2=crustal thick, 3=ftype, 4=magpdf, 5=maxmag,  6=RateType, '')')
        write (42,'( ''           7=SR, 8=a-value, 9=paleo, 10=Moment, 11=b_value flt, 12=segModel'')')
        write (42,'( 2x,''period index, period: '',i5, f10.4)') iPer, period
        write (42,'( 2x,'' Z values:'')')
        write(42,'(6x,25f12.4)') (testInten(J2),J2=1,nInten)
        write (42,'( 2x,''Mean hazard'')')
        write(42,'( 25e12.4)') (meanHaz(j),j=1,nInten)

        write (42,'( /,2x,''Sensitivity hazard'')')
        write (42,'( 2x,''iFlt, Number of Nodes, Number of Branches for each Node'')')
        write (42,'( 2x,'' iFlt, iNode, iBranch, hazard(z) '')')
        
        do iFlt=1,nFlt
c        write (42,'( 20i5)') iFlt, nNode_SSC, (nBR_SSC1(iFlt,iNode),iNode=1,nNode_SSC)
         do iNode=1,nNode_SSC
          do iBR=1,nBR_SSC(iFlt,iNode)
            if ( iNode .eq. 6 ) then
              if ( iBR .eq. 1 .and. nSR(iFlt) .eq. 0 ) goto 990
              if ( iBR .eq. 2 .and. nActRate(iFlt) .eq. 0 ) goto 990
              if ( iBR .eq. 3 .and. nRecInt(iFlt) .eq. 0 ) goto 990
              if ( iBR .eq. 4 .and. nMoRate(iFlt) .eq. 0 ) goto 990
            endif
            if ( nBR_SSC1(iFlt,iNode) .eq. 1 ) goto 990
            write (42,'( 6x, 3i5,25e12.4)') iFlt, iNode, iBR, 
     1         (haz_SSC(iFlt,iNode,iBR,iInten),iInten=1,nInten)
 990        continue
          enddo
         enddo
        enddo
        close (42)
c        pause ' test sensitivity'

c      Interpolate the desired hazard level for tornado plot
c      First find the GM for the mean hazard, interpolated to desired haz level
       do iInten=2,nInten
         if ( meanhaz(iInten-1) .ge. hazLevel .and. meanhaz(iInten) .le. hazLevel ) then
          GM0 = exp( alog(hazLevel / meanhaz(iInten-1)) / 
     1                  alog( meanhaz(iInten)/ meanhaz(iInten-1))
     2                  * alog( testInten(iInten)/testInten(iInten-1) ) + alog(testInten(iInten-1)) )
         endif
        enddo

        do iFlt=1,nFlt
         do iNode=1,nNode_SSC

          k = 0
          do iBR=1,nBR_SSC(iFlt,iNode)
            if ( iNode .eq. 6 ) then
              if ( iBR .eq. 1 .and. nSR(iFlt) .eq. 0 ) goto 995
              if ( iBR .eq. 2 .and. nActRate(iFlt) .eq. 0 ) goto 995
              if ( iBR .eq. 3 .and. nRecInt(iFlt) .eq. 0 ) goto 995
              if ( iBR .eq. 4 .and. nMoRate(iFlt) .eq. 0 ) goto 995
            endif
            if ( nBR_SSC1(iFlt,iNode) .eq. 1 ) goto 995
            if (iNode .eq. 3) goto 995
 
c           Interpolate
            k = k + 1
            do iInten=2,nInten
              if ( haz_SSC(iFlt,iNode,iBR,iInten-1) .ge. hazLevel
     1        .and. haz_SSC(iFlt,iNode,iBR,iInten)  .le. hazLevel ) then
                 GM1 = exp( alog(hazLevel / haz_SSC(iFlt,iNode,iBR,iInten-1)) / 
     1                  alog( haz_SSC(iFlt,iNode,iBR,iInten)/ haz_SSC(iFlt,iNode,iBR,iInten-1))
     2                  * alog( testInten(iInten)/testInten(iInten-1) ) + alog(testInten(iInten-1)) )              

                   GM_ratio(k) = GM1 / GM0
                 goto 995 
              endif
            enddo
 995        continue
          enddo
 996      continue

c         Check if the range is large enought to be relevant
          ratio1 = 1 - contrib_min
          ratio2 = 1 + contrib_min
          iFlag = 0
          do iBR=1,k
            if (GM_ratio(iBR) .gt. 0. ) then
              if ( GM_ratio(iBR) .lt. ratio1 .or. GM_ratio(iBR) .gt. ratio2 ) iFlag = 1         
            endif
          enddo
c          write (*,'( 3i5, 10e12.3)')  iFlt, iNode, iFlag, (GM_ratio(iBR),iBR=1,k)

c         Find max factor for this branch
          fact1 = 0.
          do iBR=1,k
            t1 = abs (  alog(GM_ratio(iBR) ) )
            if ( t1 .gt. fact1 ) fact1 = t1
          enddo

          if ( iFlag .eq. 1 ) then
c            write (43,'( 3i5,f10.3, 25f10.4)') iSite, iFlt, iNode,  fact1, 
c     1                      (GM_ratio(iBR), iBR=1,nBR_SSC(iFlt,iNode) )
            xx = 999.
            write (43,'( 3i5,f10.3,$)') iSite, iFlt, iNode,  fact1
            do iBR=1,nBR_SSC(iFlt,iNode) 
              write (43,'( f10.4,$)') GM_ratio(iBR)
            enddo
            write (43,'( 2x,a80)') fname(iFlt)
          endif
         enddo

        enddo

c       Set branch for non-poisson (node 13)
        read (31,*) nBR, (xFact(k),k=1,nBR)

c       Scale the hazard for non-poisson, but only crustal sources
        do iInten=1,nInten
          do iBR=1,nBR
            sum1 = 0.
            do kFlt=1,nFlt    
              if ( attenType(kFlt) .eq. 1 ) then
                sum1 = sum1 + meanhaz1(kFlt,iInten)* xfact(iBR)
              else 
                sum1 = sum1 + meanhaz1(kFlt,iInten)
              endif 
            enddo 
            nonPoisson(iBR,iInten) = sum1
          enddo        
        enddo

        do k=1,9
          GM_ratio(k) = 1.
        enddo

c       interpolate the meanhaz and scaled mean haz 
        do j=1,nBR
          do iInten=2,nInten
            if ( nonPoisson(j,iInten-1) .ge. hazLevel .and. nonPoisson(j,iInten)  .le. hazLevel ) then
                 GM1 = exp( alog(hazLevel / nonPoisson(j,iInten-1)) / 
     1                  alog( nonPoisson(j,iInten)/ nonPoisson(j,iInten-1))
     2                  * alog( testInten(iInten)/testInten(iInten-1) ) + alog(testInten(iInten-1)) )
                 GM_ratio(j) = GM1 / GM0
                 goto 997
            endif
          enddo
 997      continue
        enddo
        iFlt=0
        iNode = 13
        
c       Find max factor for this branch
        fact1 = 0.
        do iBR=1,nBR
         t1 = abs (  alog(GM_ratio(iBR) ) )
         if ( t1 .gt. fact1 ) fact1 = t1
        enddo
          
        write (43,'( 3i5,f10.3, 25f10.4)') iSite, iFlt, iNode, fact1, 
     1      (GM_ratio(iBR), iBR=1,3),
     2       (xx,iBR=nBR+1,9)

        call flush (43)


 1000 continue

      write (*,*) 
      write (*,*) '*** Tornado Code (45) Completed with Normal Termination ***'

      stop
      end

