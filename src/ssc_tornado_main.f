      program SSC_Tornado_Haz45
 
!c     compatible with Haz45.2

      implicit none
      include 'tornado.h'
      
      integer iPer, nFlt0, nattentype, nProb, jFlt, kFlt,
     1        i, j, k, nInten, nFlt,
     3        iInten, nWidth(MAX_FLT),
     4        attentype(MAX_FLT), nGM_model(MAX_ATTENTYPE)   
      integer f_start(MAX_FLT), f_num(MAX_FLT), nSegModel(MAX_FLT),
     2        faultflag(MAX_FLT1,MAX_SEG,MAX_FLT), nPer, jPer   
      integer iFlt, ithick1
      real testInten(MAX_INTEN), probAct(MAX_FLT),
     1     Haz_out5(MAX_FLT,MAX_WIDTH,MAXPARAM,MAX_FTYPE),
     2     al_segwt(MAX_FLT), Wt_segModel(MAX_FLT,MAX_SEG),
     3     x
      real cumWt_GM(MAX_ATTENTYPE,MAX_ATTEN)
      real specT1(MAX_PROB), hazLevel, version      
      character*80 filein, file1, fname(MAX_FLT), dummy, system_name(MAX_FLT)
      real SR_Fact(MAX_FLT1,MAX_SEG,MAX_FLT)
      
      integer nCorr1,
     1        corrNode(100), corrNFlt(100), corrFlt(100,MAX_FLT), iCorr
      integer iFlt0, iSeg, kk
      real wt_seg1(MAX_FLT,MAX_SEG)
      
      real ftype2(MAX_FLT)
      integer nFtype(MAX_FLT,MAXPARAM)
      integer jCorrFlag, NFLT2, jFlt_index(MAX_FLT)
      real hazTotal1, haz1, 
     1     wt_BR(MAX_FLT,12,MAX_THICK,MAX_BR)
      real*8 sum, sum1
      real hazTotal2(MAX_INTEN,MAX_FLT,10,12,MAX_BR1)
      real GM0, GM1, GM_Ratio(100), ratio1, ratio2, contrib_min, max_ratio
      real GM2(10)
      integer nBR1, iFlag, iii, jFlt0

      integer nNode_SSC, iNode, iThick, iSum, iOpen
      integer iBR, nBR_all(MAX_FLT,12,MAX_THICK)
      real wt_ftype(MAX_FLT,5,5)
      real wt_all(MAX_FLT,12,MAX_THICK,MAX_BR)
      real haz_mean(MAX_FLT,MAX_INTEN), haz_total(MAX_INTEN)

      integer index_flt(100), jBR, n1, iAdd
      integer flt0_flag(MAX_FLT)

      iii = 0
      iOpen = 0
            
      write (*,*) '**************************'
      write (*,*) '*   SSC Tornado Code for *'
      write (*,*) '*       HAZ 45.2 or 3    *'
      write (*,*) '*       Feb 2018         *'
      write (*,*) '* Reduced Memory version *'
      write (*,*) '**************************'

c     Open and read the run file
      write (*,*) 'Enter the input filename.'      
      read (*,'(a80)') filein
      open (31,file=filein,status='old')
      read (31,*,err=2002) nPer
      read (31,*,err=2003) HazLevel
      read (31,*,err=2011) contrib_min
      
c     Read correlated nodes from the input file
      read (31,*,err=2004) nCorr1
      write (*,'( i5,'' nCorr'')') nCOrr1
      if ( nCorr1 .gt. 100 ) then
        write (*,'( 2x,''Maximum number of correlations is 100'')')
        stop 99
      endif
      if ( nCorr1 .ne. 0 ) then
        do i=1,nCorr1
          read (31,*,err=2005) corrNode(i), corrNFlt(i)
          if ( corrNFlt(i) .gt. MAX_FLT) then
            write (*,'(2x,''Increase MAX_FLT to'',i5)') corrNFlt(i)
            stop 99
          endif
          read (31,*,err=2006) (corrFlt(i,j),j=1,corrNFlt(i))
        enddo
      endif

c     Read output file name and open file
      read (31,'( a80)',err=2007) file1
      write (*,'( a80)') file1
      open (30,file=file1,status='unknown')
      rewind (30)
      write (30,'(''SSC branches:'')')
      write (30,'(''1=crustal thick'')')       
      write (30,'(''2 = Dip  '')')       
      write (30,'(''3=ftype model   '')')       
      write (30,'(''4=mag pdf   '')')       
      write (30,'(''5=refmag   '')')       
      write (30,'(''6=Rate Method   '')')       
      write (30,'(''7=slip rate for faults   '')')       
      write (30,'(''8=a-value&b-value pair   '')')       
      write (30,'(''9=paleo recurrence interval   '')')       
      write (30,'(''10=Moment rate   '')')       
      write (30,'(''11=b_value flt for fault   '')')       
      write (30,'(''12=segModel   '',/)') 
      write (30,'(44x, '' iFlt  iNode ithick LN_Fac  BR1   BR2  ....   '')')  

c     Read hazard curve output file name and open file
      read (31,'( a80)',err=2007) file1
      write (*,'( a80)') file1
      open (40,file=file1,status='unknown')

c     Loop over each period
      do 1100 jPer=1,nPer
       read (31,*,err=2012) iPer

c      Read Input File
       call RdInput ( nInten,  testInten, nGM_model, cumWt_GM, 
     5     nattentype, attenType, nProb, iPer, specT1(jPer), version)
          
c      SSC branches: 1=crustal thick, 2 = Dip, , 3=ftype model, 4=mag pdf, 5=refmag,
c                   6=Rate Method, 7=slip rate for faults, 8=a-value&b-value pair, 
c                   9=paleo recurrence interval, 
c                   10=Moment rate, 11=b_value flt for fault, 12=segModel
       nNode_SSC = 12

c      Read the fault file 
       call Rd_Fault_Data (version, nFlt, nFlt0, probAct, AttenType, 
     2           nSegModel, Wt_segModel, nFtype, wt_ftype, f_start, f_num,
     1           faultFlag, al_Segwt, nBR_all, wt_all, fname, SR_Fact, ftype2,
     3           flt0_flag, jFlt_index, system_name )
       write (*,'( 2x,''out of Rd_fault_data'')')
                            
c      Read the out3 file (mean hazard by fault) for this period 
       call read_out3 (haz_mean, haz_total, iPer )
       write (*,'( 2x,''Out of read out3 file'')')

c      Compute the total number of widths, combining all dips and thicknesses
       do iFlt=1,nFlt
         isum = 0
         do ithick=1,nBR_all(iFlt,1,1)
c          Add the number of dips for this thickness (currently, independent of thickness)          
           isum = isum + nBR_all(iFlt,2,ithick)
         enddo
         nWidth(iFlt) = isum
       enddo

c      Initialize the weights to the starting weights
       do iNode=1,11
         call Reset_wt ( nFlt, nBR_all, wt_all, wt_BR )  
       enddo
       do iFlt0=1,nFlt0
        do iSeg=1,nSegModel(iFlt0)
         wt_seg1(iFlt0,iSeg) = Wt_segModel(iFlt0,iSeg)
        enddo
       enddo

c      Loop over intensity levels (put this in the outer loop to save memory)
       do iInten=1,nInten

c       Read the hazard from out5 for this  period and Z value      
        call read_out5 (haz_out5, iPer, iInten, nFlt, nProb, iOpen)
        write (*,'( 2x,''out of readout5 for intensity '',i5,2x,'' '',i5)') iInten, nInten

c       Sensitivity for each branch and fault
c       Loop over each fault
        do iFlt=1,nFlt
        
c        Loop over each SSC node for this fault
c        (Nodes: 1 through 11, node 12 for seg model is separate
         do iNode=1,11
         
c         For node 5, (ref mag), there is a refmag range for each thickness
          if ( iNode .eq. 5 ) then
            n1 = nBR_all(iFlt,1,1)
            if ( n1 .gt. 10 ) stop 88
          else
            n1 = 1
          endif         
         
          do iThick=1,n1
         
c          skip Node if there is only one branch for this node
           if ( nBR_all(iFlt,iNode,ithick) .le. 1 .and. iNode .ne. 1 ) goto 100  

c          Check if this is a correlated node
           jCorrFlag = 0
           do iCorr=1,nCorr1
            if (corrNode(iCorr) .eq. iNode ) then
             do jFlt=1,corrNFlt(iCorr)  
              if ( iflt .eq. corrFlt(iCorr,jFlt)) then
               jCorrFlag = iCorr
               goto 40
              endif
             enddo
            endif   
           enddo
 40        continue

c          If not correlated, remove just this fault from the total hazard          
           if (jCorrFlag .eq. 0 ) then
            hazTotal1 = haz_Total(iInten) - haz_mean(iFlt,iInten)
            index_flt(1) = iFlt
           endif
         
c          If correlated, then remove all of the correlated faults from the total hazard
           if (jCorrFlag .gt. 0 ) then
            hazTotal1 = haz_Total(iInten)
            do i=1,corrNFlt(jCorrFlag)
             jFlt = corrFlt(jCorrFlag,i)
             hazTotal1 = hazTotal1 - haz_mean(jFlt,iInten)
            enddo
           endif
           
c          Loop over branches for this node
           do iBR = 1, nBR_all(iFlt,iNode,ithick)
          
c           iF not correlated, reset weights and calc hazard for this fault only
            if (jCorrFlag .eq. 0 ) then

c            Reset the weights (set weight=1 for the selected branch, 0 for all others)
             do k=1,nBR_all(iFlt,1,1)
              do jBR=1,nBR_all(iFlt,iNode,k)
               wt_BR(iFlt,iNode,k,jBR) = 0.
              enddo
              wt_BR(iFlt,iNode,k,iBR) = 1.
             enddo
             
c            For node 5 (ref mag), set the weights for a single thickness
             if (iNode .eq. 5) then
              do k=1,nBR_all(iFlt,1,1)
               do jBR=1,nBR_all(iFlt,1,1)
                wt_BR(iFlt,1,k,jBR) = 0.
               enddo
               wt_BR(iFlt,1,k,ithick) = 1.
              enddo
             endif                          

             call calc_haz ( iFlt, nBR_all, haz_out5, wt_BR,  
     1                  al_segwt, wt_ftype, SR_Fact, faultFlag, haz1, flt0_flag,
     2                  nSegModel, Wt_segModel, nFtype, jFlt_index(iFlt), iNode )
             sum = haz1
            endif

c           iF correlated, reset weights and calc hazard for all correlated faults
            if (jCorrFlag .gt. 0 ) then

c            Set the weights for sensitivity 
c            (set weight=1 for the selected branch, 0 for all others)
             do i=1,corrNFlt(jCorrFlag)
              index_flt(i) = corrFlt(jCorrFlag,i)
              jFlt = corrFlt(jCorrFlag,i)
              do k=1,nBR_all(jFlt,1,1)
               do jBR=1,nBR_all(jFlt,iNode,k)
                wt_BR(jFlt,iNode,k,jBR) = 0.
               enddo
               wt_BR(jFlt,iNode,k,iBR) = 1.
              enddo

c            For node 5 (ref mag), set the weights for a single thickness
             if (iNode .eq. 5) then
              do k=1,nBR_all(jFlt,1,1)
               do jBR=1,nBR_all(jFlt,1,1)
                wt_BR(jFlt,1,k,jBR) = 0.
               enddo
               wt_BR(jFlt,1,k,ithick) = 1.
              enddo
             endif                          


             enddo
             
c            Add the hazard contribution from the correlated faults
             sum = 0.
             do i=1,corrNFlt(jCorrFlag)
              kFlt = index_flt(i) 
              call calc_haz ( kFlt, nBR_all, haz_out5, wt_BR,  
     1                  al_segwt, wt_ftype, SR_Fact, faultFlag, haz1, flt0_flag,
     2                  nSegModel, Wt_segModel, nFtype, jFlt_index(kFlt), iNode )
              sum = sum + haz1
             enddo
            endif

            hazTotal2(iInten,iFlt,iThick,iNode,iBR) = sum + hazTotal1           

c            write (*,'( 5i5,e12.5)') iInten,iFlt,iThick,iNode,iBR,
c     1         hazTotal2(iInten,iFlt,iThick,iNode,iBR)
                   
           enddo
c          end loop over branch          

 100       continue         

c          Reset the weights to the starting weights
           call Reset_wt ( nFlt, nBR_all, wt_all, wt_BR )  

          enddo
c         end loop over thickness

         enddo
c        end loop over node         
         
        enddo
c       end loop over faults    

c --------------
c       Sensitivity to the Segmentation model

c       Reset the weights to the starting weights
        call Reset_wt ( nFlt, nBR_all, wt_all, wt_BR )  

c       evaluate sensitivity for each fault system
        do iFlt0=1,nFlt0
         if ( nSegModel(iFlt0) .eq. 1 ) goto 101

c        set segmodel weights to starting values         
         do jFlt0=1,nFlt0
          do iSeg=1,nSegModel(jFlt0)
            wt_seg1(jFlt0,iSeg) = wt_segModel(jFlt0,iSeg)
          enddo
         enddo

c        Sensitivity for each segmentation model         
         do iSeg=1,nSegModel(iFlt0)

c         Set the seg model weights to unity for this branch
          do jBR=1,nSegModel(iFlt0)
            wt_seg1(iFlt0,jBR) = 0.
          enddo
          wt_seg1(iFlt0,iSeg) = 1.
          
          sum = 0.
          sum1 = 0.
          iNode = 12
          do iFlt=1,nFlt

            call calc_haz ( iFlt, nBR_all, haz_out5, wt_BR,  
     1                  al_segwt, wt_ftype, SR_Fact, faultFlag, haz1, flt0_flag,
     2                  nSegModel, Wt_seg1, nFtype, jFlt_index(iFlt), iNode )
            sum = sum + haz1
            sum1 = sum1 + haz_mean(iflt,iInten)

          enddo
          hazTotal2(iInten,iFlt0,1,iNode,iSeg) = sum           
         enddo
     1      
 101     continue
         
        enddo
c       end loop over source groups

       enddo
c      end loop over intensity       

c      Interpolate the desired hazard level for tornado plot
c      First find the GM for the mean hazard, interpolated to desired haz level
       do iInten=2,nInten
        if ( haz_total(iInten-1) .ge. hazLevel .and. haz_total(iInten) .le. hazLevel ) then
          GM0 = exp( alog(real(hazLevel / haz_total(iInten-1))) / 
     1          alog( haz_total(iInten)/ haz_total(iInten-1))
     2         * alog( testInten(iInten)/testInten(iInten-1) ) + alog(testInten(iInten-1)) )
        endif
       enddo

c      Compute the sensitivity for each node for each fault        
       do iFlt=1,nFlt
        do iNode=1,11

c        For node 5, (ref mag), there is a refmag range for each thickness
         if ( iNode .eq. 5 ) then
            n1 = nBR_all(iFlt,1,1)
         else
            n1 = 1
         endif         
         
         do iThick=1,n1
          
c         Set the number of branches by adding up all cases for ref mag (node5)          
          nBR1 = nBR_all(iFlt,iNode,1)
c          write (*,'( 4i5,2x,'' nBR1'')') iFlt, iNode, iThick, nBR1
          if ( nBR1 .le. 1 .and. iNode .ne. 1 ) goto 995
                   
c         Check if there are at least two different rate methods used
          if ( iNode .eq. 6) then
           if ( nBR_all(iFlt,7,1) .ne. 0 ) then
            if ( nBR_all(iFlt,8,1) .ne. 0) goto 80
            if ( nBR_all(iFlt,9,1) .ne. 0) goto 80
            if ( nBR_all(iFlt,10,1) .ne. 0) goto 80
           endif
           if ( nBR_all(iFlt,8,1) .ne. 0 ) then
            if ( nBR_all(iFlt,9,1) .ne. 0) goto 80
            if ( nBR_all(iFlt,10,1) .ne. 0) goto 80
           endif
           if ( nBR_all(iFlt,9,1) .ne. 0 ) then
            if ( nBR_all(iFlt,10,1) .ne. 0) goto 80
           endif
           goto 995
          endif
  80      continue        

c         Compute Sensitivity for each branch
          do iBR=1,nBR1
  
c          Interpolate to the desired hazard level
           GM_ratio(iBR) = -999.
           do iInten=2,nInten
            if ( hazTotal2(iInten-1,iFlt,ithick,iNode,iBR) .ge. hazLevel
     1         .and. hazTotal2(iInten,iFlt,ithick,iNode,iBR)  .le. hazLevel ) then
              GM1 = exp( alog(hazLevel / hazTotal2(iInten-1,iFlt,iThick,iNode,iBR)) / 
     1                  alog( hazTotal2(iInten,iFlt,ithick,iNode,iBR)/ 
     1                       hazTotal2(iInten-1,iFlt,ithick,iNode,iBR))
     2                  * alog( testInten(iInten)/testInten(iInten-1) ) 
     2                  + alog(testInten(iInten-1)) )
     
c             save the GM for thickness branch for later use with node 5
              if ( iNode .eq. 1 ) then
                GM2(iBR) = GM1
c                write (*,'( i5,f10.3)') iBR, GM2(iBR)
              endif

c             set the GM ratio for this branch.  If iNode=5 (refmag), then
c             use the GM for this thickness.  Otherwise, use the GM for the total haz
              if ( iNode .eq. 5 .and. nBR_all(iFlt,1,1) .gt. 1 ) then
               GM_ratio(iBR) = GM1 / GM2(ithick)             
              else
               GM_ratio(iBR) = GM1 / GM0
              endif
c              write (*,'( 3i5,4f10.3, ''GM0, GM1, GM2'')') iflt, iNode, ithick, GM0, GM1, GM2(ithick), GM_ratio(iBR)
c              if ( iNode .eq. 5 ) pause
              goto 990 
            endif
           enddo
  990      continue

          enddo
c         end loop over branches

c         Check if the range is large enought to be relevant
          ratio1 = 1 - contrib_min
          ratio2 = 1 + contrib_min
          iFlag = 0
          max_ratio = 0.
          do iBR=1,nBR1
           if (GM_ratio(iBR) .gt. 0. ) then
             if (abs(alog(GM_ratio(iBR))) .gt. max_ratio) max_ratio =abs(alog(GM_ratio(iBR))) 
             if ( GM_ratio(iBR) .lt. ratio1 .or. GM_ratio(iBR) .gt. ratio2 ) then
c               write (*,'( 3f10.3)') GM_ratio(iBR), ratio1, ratio2
c               pause
               iFlag = 1       
             endif
           endif
          enddo
         
c         check if the ratio is large enough to keep         
          if ( iFlag .eq. 1 ) then
         
c          Check if this is a correlated node (if so, then list only once)
           jCorrFlag = 0
           do iCorr=1,nCorr1
            if (corrNode(iCorr) .eq. iNode ) then
             do jFlt=1,corrNFlt(iCorr)  
              if ( iflt .eq. corrFlt(iCorr,jFlt) .and. jflt .eq. 1) then
               jCorrFlag=1
              elseif (iflt .eq. corrFlt(iCorr,jFlt) ) then
               goto 41
              endif
             enddo
            endif   
           enddo

           if ( inode .ne. 5 ) then
             ithick1 = 0
           else
             ithick1 = ithick
           endif
                    
           if ( jCorrFlag .eq. 0 ) then
             write (30,'( a30, 10x, 3i7,20f7.2)') fname(iFlt), iFlt, iNode, ithick1, 
     1            max_ratio, (GM_ratio(iBR), iBR=1,nBR1)
           else
             write (30,'( a30,''correlated'', 3i7,20f7.2)') fname(iFlt), iFlt, iNode,
     1         ithick1, max_ratio, (GM_ratio(iBR), iBR=1,nBR1)
           endif
           
 41        continue
          endif

 995      continue

         enddo
c        end loop over thicknesses

        enddo
c       end loop over node
        
       enddo
c      end loop over flt

c --------------
c      Compute the sensitivity for node 12 for each source group
       iNode = 12
       iThick = 1
       write (30,'(/,44x, '' iFlt0 iNode ithick LN_Fac  iSeg1 iSeg2  ....   '')')  
       do iFlt0=1,nFlt0
       
c        do iBR=1,nSegModel(iFlt0)
c         write (*,'(2i5,20e12.4 )') iFlt0, iBR, (hazTotal2(i,iFlt0,iNode,iBR), i=1,nInten)
c        enddo
       
c       Set the number of branches by adding up all cases for ref mag (node5)          
        nBR1 = nSegModel(iFlt0)
        if ( nBR1 .le. 1 ) goto 996

c       Compute sensitivity for each branch of the segmodels         
        do iBR=1,nBR1
  
c        Interpolate to the desired hazard level
         GM_ratio(iBR) = -999.
         do iInten=2,nInten
           if ( hazTotal2(iInten-1,iFlt0,iThick,iNode,iBR) .ge. hazLevel
     1         .and. hazTotal2(iInten,iFlt0,iThick,iNode,iBR)  .le. hazLevel ) then
             GM1 = exp( alog(hazLevel / hazTotal2(iInten-1,iFlt0,iThick,iNode,iBR)) / 
     1                  alog( hazTotal2(iInten,iFlt0,iThick,iNode,iBR)/ 
     1                       hazTotal2(iInten-1,iFlt0,iThick,iNode,iBR))
     2                  * alog( testInten(iInten)/testInten(iInten-1) ) 
     2                  + alog(testInten(iInten-1)) )
             GM_ratio(iBR) = GM1 / GM0
             goto 991 
           endif
         enddo
  991    continue

        enddo
c       end loop over branches

c       Check if the range is large enought to be relevant
        ratio1 = 1 - contrib_min
        ratio2 = 1 + contrib_min
        iFlag = 0
        max_ratio = 0.
        do iBR=1,nBR1
          if (GM_ratio(iBR) .gt. 0. ) then
            if (abs(alog(GM_ratio(iBR))) .gt. max_ratio) max_ratio =abs(alog(GM_ratio(iBR))) 
            if ( GM_ratio(iBR) .lt. ratio1 .or. GM_ratio(iBR) .gt. ratio2 ) then
              iFlag = 1       
            endif
          endif
        enddo
         
c       check if the ratio is large enough to keep         
        if ( iFlag .eq. 1 ) then
          ithick1 = 0
          write (30,'( a40, 3i7,20f7.2)') system_name(iFlt0), iFlt0, iNode, iThick1, max_ratio, 
     1        (GM_ratio(iBR), iBR=1,nBR1)
        endif

 996    continue
       enddo
c      end loop over node
       
 1100 continue
c     end loop over period

c     Write out hazard curves
      do iFlt=1,nFlt
       do iNode=1,11

c       For node 5, (ref mag), there is a refmag range for each thickness
        if ( iNode .eq. 5 ) then
            n1 = nBR_all(iFlt,1,1)
        else
            n1 = 1
        endif   
        

        do iThick=1,n1
         do iBR=1,nBR_all(iFlt,iNode,ithick)
           write (40,'(4i5, 20e12.4)') iFlt, iNode, ithick,iBR, 
     1         (hazTotal2(k,iFlt,iThick,iNode,iBR),k=1,nInten)
c           write (40,'(4i5, 20e12.4)') iFlt, iNode, ithick,iBR
         enddo
        enddo
       enddo
      enddo   
        

      write (*,*) 
      write (*,*) '*** SSC Tornado Code (45) Completed with Normal Termination ***'

      stop
 2000 write (*,'( 2x,''input error, random seed'')')
      goto 3000     
 2001 write (*,'( 2x,''input error, nSample'')')
      goto 3000     
 2002 write (*,'( 2x,''input error, nPer'')')
      goto 3000     
 2003 write (*,'( 2x,''input error, haz level'')')
      goto 3000     
 2004 write (*,'( 2x,''input error, nCOrr'')')
      goto 3000     
 2005 write (*,'( 2x,''input error, corrNode'')')
      goto 3000     
 2006 write (*,'( 2x,''input error, corrFlt'')')
      goto 3000     
 2007 write (*,'( 2x,''input error, output file name'')')
      goto 3000   
 2010 write (*,'( 2x,''input error, select mech for fractile'')')
      goto 3000  
 2011 write (*,'( 2x,''input error, Min GM ratio'')')
      goto 3000  
 2012 write (*,'( 2x,''input error, iPer'')')
      goto 3000  
      
 3000 backspace (31)
      read (31,'( a80)') dummy
      x = 1
      write (*,4000) dummy
      stop 99
 4000 format ( a80)
      end


c -------------------------------------------------------------------------

      subroutine Reset_wt ( nFlt, nBR_all, wt_all, wt_BR )  
      implicit none
      include 'tornado.h'

      real wt_all(MAX_FLT,12,MAX_THICK,MAX_BR)
      real wt_BR(MAX_FLT,12,MAX_THICK,MAX_BR)
      integer iBR, nBR_all(MAX_FLT,12,MAX_THICK)
      integer iFlt, nFlt, nTHick, iTHick, iNode
      
      do iFLt=1,nFlt
       nThick = nBR_all(iFlt,1,1) 
       do iThick=1,nThick
        do iNode=1,11
         do iBR=1,nBR_all(iFlt,iNode,ithick) 
           wt_BR(iFlt,iNode,iThick,iBR) = wt_all(iFlt,iNode,iThick,iBR)
         enddo
        enddo
       enddo
      enddo
      return
      end
