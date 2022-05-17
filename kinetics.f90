
!
! File: kinetics.f90
! Description: Implements Gillespie's SSA algorithm for chemical reaction simulation, as
!              well as direct simulation.
! Author: Scott Habershon
! Start date: 6 July 2015
!
! Change log:
!
! 6 July 2015  - Original implementation
! 13 July 2016 - Added direct kinetics simulation, tidied code.
!

Program SSA
  implicit none
  include 'globals.h'

  integer :: narg
  integer :: log_file
  integer :: nmol
  integer :: status
  integer :: i, idum, ii, jj, id, irxn, j, kk
  integer :: it
  integer :: nrxn
  integer :: rid(3,NRXNMAX)
  integer :: pid(3,NRXNMAX)
  integer :: idtag(NMOLMAX)
  integer :: nt
  integer :: irun
  integer :: ifile
  integer :: ntmax
  integer :: ic, k
  integer :: ir, nrun, ibin,ichange,jdum,nsample

  character (len=50)  :: filename(4)
  character (len=8)   :: c_date
  character (len=10)  :: c_time
  character (len=90)  :: line, comment
  character (len=30)  :: molid(NMOLMAX)
  character (len=1)   :: rxntype(NRXNMAX), prodtype(NRXNMAX)
  character (len=90)  :: rxnid(NRXNMAX)
  character (len=3)   :: c3
  character (len=8)   :: fmt
  character (len=3)   :: GTStype

  real(8) :: stoich(NRXNMAX,NMOLMAX)
  real(8) :: Gmol(NMOLMAX)
  real(8) :: pop(NMOLMAX)
  real(8) :: Gts(NRXNMAX)
  real(8) :: Gts_mean(NRXNMAX)
  real(8) :: Tstop
  real(8) :: timestep
  real(8) :: Temperature
  real(8) :: ktst(NRXNMAX)
  real(8) :: dG, Gi, Gf
  real(8) :: propens(NRXNMAX)
  real(8) :: s1,s2
  real(8) :: tau, sum, sumprop, time
  real(8) :: ran2
  real(8) :: volume, bincount(NMOLMAX, NOUT)
  real(8) :: Kdiff, Keq,del,cond
  real(8) :: pressure, conc(NMOLMAX), pop_init_store(NMOLMAX),psd(NMOLMAX,NOUT),pbin(NMOLMAX,NOUT)
  real(8) :: flux(NRXNMAX), fluxtot, sflux(NMOLMAX), popsd(NMOLMAX,NOUT),Tbin(NOUT)
  real(8) :: Tout(NOUT), popout(NMOLMAX,NOUT), popsave(NMOLMAX,NRUNMAX,NOUT),Tsave(NRUNMAX,NOUT)
  real(8) :: covar(NRXNMAX,NRXNMAX),Fnew,Fold,de,prob,Gold,step,ratio
  real(8), allocatable :: cov_new(:,:), inv_cov_new(:,:)
  logical :: gvconnect(NRXNMAX), dynamic_graph, Readout(NOUT), rfound(NRXNMAX), readrates
  logical :: usecovar,readout_temp(NOUT)

  character*3 :: Method ! 'SSA' or 'DIR'

  do i = 1, NOUT
     Readout(i) = .FALSE.
  enddo


  ! Set the output index of the log file:
  !
  log_file = 11


  ! Get the current date and time.
  !
  call date_and_time(date=c_date, time=c_time)


  ! Read the input filename, the starting structure file and
  ! the target structure file from the command line.
  !
  Call ReadInputs(Method,Tstop,Temperature,irun,volume,ntmax,Kdiff,timestep,filename, &
       c_date,c_time,log_file,GTStype,nrun,readrates,usecovar)


  ! Initialize random numbers.
  !
  Call StartRan( irun )


  ! Open and read the molecules file.
  !
  Call ReadMolFile(log_file,filename,nmol,idtag,molid,Gmol)



  ! Open and read the reactions file - note that we only read the
  ! forward reactions here, then double up to get the reverse reactions.
  !
  Call ReadRxnFile(log_file,filename,nrxn,rxnid,rxntype,prodtype,gvconnect,rid,pid, &
       Gts,Gmol,readrates,ktst)


  ! Open and read the initial state file.
  !
  Call ReadInitialState(log_file,filename,pop,conc,Method,pressure,volume, &
       Temperature,molid,nmol)


  ! Set up the stoichiometry matrix for each reaction..
  !
  Call GetStoichiometryMatrix(stoich,nrxn,rid,pid,rxntype,prodtype)


  ! Set up the rate vector for all reactions.
  !
  if (.not.readrates) then
    Call GetRates(log_file, ktst, rxntype, prodtype, Kdiff, Gmol, Gts, nrxn, &
    pid, rid, Temperature,GTStype)
    call flush(log_file)
  endif


  ! Output a static GraphViz representation of the network.
  !
  ifile = 21
  open( ifile, file = 'network.dat', status = 'unknown' )
  dynamic_graph = .false.
  Call GraphViz( nmol, nrxn, ktst, molid, NRXNMAX, NMOLMAX, dynamic_graph, Gts, ifile, &
       rid, pid,rxntype,prodtype,flux,gvconnect,SMALL)
  close(unit = ifile)



  ! We're now ready to run the dynamics simulations...
  !
  if (Method == 'SSA') then

     ! Zero the averaged populations and output times.
     !
     psd(:,:) = 0.0
     pbin(:,:) = 0.0
     Tbin(:) = 0.0
     popsave(:,:,:) = 0.0
     Tsave(:,:) = 0.0


     ! Setup the output time-bin
     !
     do i = 1, NOUT
        Tbin(i) = dble(i) * Tstop / dble(NOUT)
     enddo
     del = Tstop / dble(Nout-1)   ! in s.


     ! Store the initial concentrations.
     !
     pop_init_store(:) = pop(:)


     ! If using covariance matrix, store the input Gts.
     !
     if (usecovar) then
        Gts_mean(:) = Gts(:)
     endif

     ! Loop over individual runs.
     !
     do ir = 1, nrun
91      continue

        pop(:) = pop_init_store(:)

        ! NEW - if usecovar = .TRUE., then we sample a new Gts in a normal distribution with the
        !       covariance matrix read in from 'covar.dat'. Note that covar.dat MUST contain
        !       the elements of an nrxn x nrxn matrix. After generating a new Gts, we need to
        !       recalculate the rection rates.
        !       The covariance matrix should be given in ATOMIC UNITS (Eh**2).
        !
        if (usecovar) then
          
           allocate(cov_new(nrxn/2,nrxn/2))
           allocate(inv_cov_new(nrxn/2,nrxn/2))

           open(10,file='covar.dat',status='unknown')
           do i = 1, nrxn/2
              do j = 1, nrxn/2
                 read(10,*)idum,jdum,covar(i,j)
                 covar(i,j) = covar(i,j) * (tokjmol * tokjmol)
                 cov_new(i,j) = covar(i,j)
              enddo
           enddo
           close(10)

           ! Calculate inverse covariance matrix:
           !
           Call SVD_inverse(cov_new, inv_cov_new, nrxn/2, cond, 1d-12)
            

           ! Set number of MC samples and step-size (Hartrees)
           !
           nsample = 10000
           step = 1.00 

           ! Reset Gts.
           do i = 1, nrxn
              Gts(i) = Gts_mean(i)
           enddo


           ! Set initial covariance function to be equal to 1.0. In other words,
           ! We start with all Gts = Gts_mean.
           !
           ! Note that this sampling is only performed for the forward reactions!
           ! We deal with the beackward reactions afterwards....

           Fold = 1.0

           do i = 1, nsample

              ! Choose reaction to change.
              !
              ichange = 1 + int( ran2(irun,0.d0,dble(nrxn/2)) )
              if (ichange.gt.nrxn/2)then
                print*,'RANDOM NUMBER ERROR'
                stop
              endif

              Gold = Gts(ichange)

              ! Make a change.
              de = ran2(irun,-step,+step)
              Gts(ichange) = Gts(ichange) + de

              ! Calculate new covariance function
              sum = 0.d0
              do j = 1, nrxn/2
                 do k = 1, nrxn/2
                  sum = sum + 0.5 * (Gts(j) - Gts_mean(j)) *inv_cov_new(j,k)* (Gts(k) - Gts_mean(k))
                 enddo
              enddo
              Fnew = exp(-sum)

              ratio = Fnew / Fold

              ! Calculate acceptance/rejection
              prob = ran2(irun,0.d0,1.d0)
              if (prob.ge.ratio) then      ! Reject !
                 Gts(ichange) = Gts(ichange) - de
              else
                 Fold = Fnew
              endif

           enddo

           ! Now we have to calculate the barriers for the backwards reactions - we do this
           ! by calculating the sampled change in the original forward reaction, and
           ! then adding this onto the original backward barriers:
           !
           do i = 1, nrxn/2
             de = Gts(i) - Gts_mean(i)
             jj = i + (nrxn/2)
             Gts(jj) = Gts_mean(jj) + de
             write(6,*)'BACKWARDS: ',jj,Gts(jj),Gts_mean(jj)
             if (Gts(jj).lt.0.d0)then
               write(6,*) 'ERROR: NEGATIVE BACKWARDS BARRIER AFTER SAMPLING COVARIANCE '
               deallocate(cov_new)
               deallocate(inv_cov_new)
               goto 91
             endif
           enddo

           ! Set up the rate vector for all reactions.
           !
           if (.not.readrates) then
              Call GetRates(log_file, ktst, rxntype, prodtype, Kdiff, Gmol, Gts, nrxn, &
              pid, rid, Temperature,GTStype)
              call flush(log_file)
           else
              write(6,*)'* ERROR: Cannot use covariance with readrates=.false.'
              stop
           endif
           do i = 1, nrxn
             write(66,*)i,Gts(i),Gts_mean(i)
           enddo

           deallocate(cov_new)
           deallocate(inv_cov_new)

        endif

        do i = 1, NOUT
           readout_temp(i) = .FALSE.
        enddo
        write(6,*)'INTO STOCHASTIC';call flush(6)
        Call Stochastic(Tstop,pop,ktst,rid,pid,rxntype,prodtype,volume,nrxn,nmol,irun,stoich,flux, &
             ir,popsave,Tsave,log_file,fluxtot,readout_temp)
        write(6,*)'OUT OF STOCHASTIC';call flush(6)

        do i = 1, NOUT 
           if (readout_temp(i)) then
               readout(i) = .true.
           endif
        enddo

     enddo


     ! Bin averaged populations.
     !
     pbin(:,:) = 0.d0
     bincount(:,:) = 0.d0
     do j = 1, nrun
        do i = 1, NOUT
           do ibin = 1, NOUT
              if (Tsave(j,i) > Tbin(ibin) - 0.5d0 * del .and. Tsave(j,i) < Tbin(ibin)+0.5d0*del) then

                 do k = 1, nmol
                    pbin(k,ibin) = pbin(k,ibin) + popsave(k,j,i)
                    bincount(k,ibin) = bincount(k,ibin) + 1.d0
                 enddo

              endif
           enddo
        enddo
     enddo


     ! Determine the averaged populations in each time bin for each molecule.
     !
     do k = 1, nmol
        do ibin = 1, NOUT
           pbin(k,ibin) = pbin(k,ibin) / bincount(k,ibin)
        enddo
     enddo


     ! Now determine the SD in each bin for each molecule.
     !
     psd(:,:) = 0.d0
     bincount(:,:) = 0.d0
     do j = 1, nrun
        do i = 1, NOUT
           do ibin = 1, NOUT
              if (Tsave(j,i) > Tbin(ibin) - 0.5d0 * del .and. Tsave(j,i) < Tbin(ibin)+0.5d0*del) then
                 do k = 1, nmol
                    psd(k,ibin) = psd(k,ibin) + (popsave(k,j,i) - pbin(k,ibin))**2
                    bincount(k,ibin) = bincount(k,ibin) + 1.d0
                 enddo
              endif
           enddo
        enddo
     enddo

     do k = 1, nmol
        do ibin = 1, NOUT
           if (bincount(k,ibin) > 1.1d0) then
              psd(k,ibin) = psd(k,ibin) / (bincount(k,ibin)-1.d0)
              psd(k,ibin) = sqrt(psd(k,ibin))
              psd(k,ibin) = psd(k,ibin) / sqrt(bincount(k,ibin))
           else
              psd(k,ibin) = 0.d0
           endif
        enddo
     enddo


     ! Output averaged....
     !
     do i = 1, nmol
        fmt = '(I3.3)'
        write(c3,fmt)i
        open(67,file = 'popave'//trim(c3)//'.dat',status='unknown')
        do ibin = 1, NOUT
          if (readout(ibin)) then
          if (bincount(i,ibin) > 1.1d0) then
            write(67,'(3(1x,e14.6))')Tbin(ibin),(pbin(i,ibin)/6.02e23)/(volume*1d3), &
            (psd(i,ibin)/6.02e23)/(volume*1d3)
          endif
          endif
        enddo
        close(67)
     enddo


     ! Output population files - individual runs.
     !
     do i = 1, nmol
        fmt = '(I3.3)'
        write(c3,fmt)i
        open(67,file = 'pop'//trim(c3)//'.dat',status='unknown')

        do ir = 1, nrun
           do j = 1, NOUT
              if (readout(j))then
              write(67,'(2(1x,e14.6))')Tsave(ir,j),(popsave(i,ir,j)/6.02e23)/(volume*1d3)
              endif
           enddo
           write(67,*)
        enddo
        ifile = ifile + 1

     enddo




  else if (Method == 'DIR') then
     stop 'DIRECT simulation method not yet complete'
     ! Call Direct(Tstop,conc,ktst,rid,pid,rxntype,prodtype,volume,nrxn,nmol,irun, &
     !     stoich,flux,timestep)
  endif



  ! Output the reaction flux.
  !
  open(30,file = 'flux_forward.dat',status='unknown')
  open(31,file = 'flux_reverse.dat',status='unknown')
  open(32,file = 'flux_net.dat',status='unknown')
  ic = 0
  do i = 1, nrxn/2
!     write(6,*)'WTF: ',i,flux(i),fluxtot
     flux(i) = flux(i) / fluxtot
     flux(i+nrxn/2) = flux(i+nrxn/2) / fluxtot
     if (flux(i) > 0.d0) then
        ic = ic + 1
        if (rxntype(i) == 'U'.and.prodtype(i) == 'U') then
           write(30,'("#",2x,a,1x,a1,1x,a1,i5,i5)')trim(adjustl(rxnid(i))),rxntype(i),prodtype(i),&
                   rid(1,i),pid(1,i)
           write(31,'("#",2x,a,1x,a1,1x,a1,i5,i5)')trim(adjustl(rxnid(i))),rxntype(i),prodtype(i),&
                   rid(1,i),pid(1,i)
           write(32,'("#",2x,a,1x,a1,1x,a1,i5,i5)')trim(adjustl(rxnid(i))),rxntype(i),prodtype(i),&
                   rid(1,i),pid(1,i)
        else if (rxntype(i) == 'B'.and.prodtype(i) == 'U') then
           write(30,'("#",2x,a,1x,a1,1x,a1,i5,i5,i5)')trim(adjustl(rxnid(i))),rxntype(i),prodtype(i), &
                rid(1,i),rid(2,i),pid(1,i)
           write(31,'("#",2x,a,1x,a1,1x,a1,i5,i5,i5)')trim(adjustl(rxnid(i))),rxntype(i),prodtype(i), &
                rid(1,i),rid(2,i),pid(1,i)
           write(32,'("#",2x,a,1x,a1,1x,a1,i5,i5,i5)')trim(adjustl(rxnid(i))),rxntype(i),prodtype(i), &
                rid(1,i),rid(2,i),pid(1,i)

        else if (rxntype(i) == 'U'.and.prodtype(i) == 'B') then
           write(30,'("#",2x,a,1x,a1,1x,a1,i5,i5,i5)')trim(adjustl(rxnid(i))),rxntype(i),prodtype(i),&
                rid(1,i),pid(1,i),pid(2,i)
           write(31,'("#",2x,a,1x,a1,1x,a1,i5,i5,i5)')trim(adjustl(rxnid(i))),rxntype(i),prodtype(i),&
                rid(1,i),pid(1,i),pid(2,i)
           write(32,'("#",2x,a,1x,a1,1x,a1,i5,i5,i5)')trim(adjustl(rxnid(i))),rxntype(i),prodtype(i),&
                rid(1,i),pid(1,i),pid(2,i)

        else if (rxntype(i) == 'B'.and.prodtype(i) == 'B') then
           write(30,'("#",2x,a,1x,a1,1x,a1,i5,i5,i5,i5)')trim(adjustl(rxnid(i))),rxntype(i),prodtype(i), &
                rid(1,i),rid(2,i),pid(1,i),pid(2,i)
           write(31,'("#",2x,a,1x,a1,1x,a1,i5,i5,i5,i5)')trim(adjustl(rxnid(i))),rxntype(i),prodtype(i), &
                rid(1,i),rid(2,i),pid(1,i),pid(2,i)
           write(32,'("#",2x,a,1x,a1,1x,a1,i5,i5,i5,i5)')trim(adjustl(rxnid(i))),rxntype(i),prodtype(i), &
                rid(1,i),rid(2,i),pid(1,i),pid(2,i)

        endif

        write(30,*)ic,flux(i),ktst(i)
        write(31,*)ic,-flux(i+nrxn/2),ktst(i+nrxn/2)
        write(32,*)ic,flux(i)-flux(i+nrxn/2)

     endif
  enddo
  close(30)
  close(31)
  close(32)


  ! Output the structure flux.
  !
  open(30,file = 'sflux.dat',status='unknown')
  do i = 1, nmol
  !   sflux(i) = sflux(i) / fluxtot
!     if (sflux(i) > 0.d0) then
        write(30,'(i4,2x,f14.8,2x,i4,a)')i,sflux(i),idtag(i),molid(i)
!     endif
  enddo
  close(30)


  ! Output a flux-weighted network graph.
  !
  ifile = 21
  open( ifile, file = 'flux_network.dat', status = 'unknown' )
  dynamic_graph = .true.
  Call GraphViz( nmol, nrxn, ktst, molid, NRXNMAX, NMOLMAX, dynamic_graph, Gts, ifile, &
       rid, pid,rxntype,prodtype,flux,gvconnect,SMALL)
  close(unit = ifile)


end Program SSA



!
!*****************************************************************************************
!
! SUBROUTINE: GraphViz
!
! Outputs a GraphViz file for creating pictures of the network using Graphviz.
!
!*****************************************************************************************
!
Subroutine GraphViz( nmol, nrxn, ktst, molid, NRXNMAX, NMOLMAX, dynamic_graph, Gts, ifile, &
     rid, pid,rxntype,prodtype,flux,gvconnect,SMALL)
  implicit none
  integer :: nmol
  integer :: nrxn
  integer :: NMOLMAX, NRXNMAX
  integer :: ifile
  integer :: rid(3,NRXNMAX)
  integer :: pid(3,NRXNMAX)
  integer :: ii,jj,kk,ll, i
  real(8) :: ktst(NRXNMAX), weight, flux(NRXNMAX), scale
  real(8) :: Gts(NRXNMAX),SMALL
  character (len=1)  :: rxntype(NRXNMAX), prodtype(NRXNMAX)
  character (len=30)  :: molid(NMOLMAX),line,line2,line3,line4,line5
  logical :: dynamic_graph,gvconnect(NRXNMAX)
  logical :: gvc_flag(NMOLMAX), icon(NMOLMAX,NMOLMAX)
  character (len=120) :: node_option,edge_option,TS_option,B_option
  logical, allocatable :: found(:)

  node_option = '[shape=circle,fontname=Helvetica,fontcolor=white,fontsize=48,color=grey,  &
       style=filled,fillcolor="#3366CC"];'
  TS_option = '[height=1.0,width=1.0,label="",shape=square,&
       fontsize=16,color=white,style=filled,fillcolor="#cc0033"];'
  B_option = '[height=1.0,width=1.0,label="",shape=diamond,color=white,style=filled,fillcolor="#669933"];'
  edge_option = '[penwidth=2,weight=2,arrowsize=2,color="#3f4246",dir="both"]'

  ! Write the header.
  !
  write(ifile,'("graph G {")')

  scale = -1.d0

  ! Write temporary options:
  !
  write(ifile,'("overlap=scalexy;")')
  write(ifile,'("splines=true;")')
  write(ifile,'("nodesep=0.8;")')
  !write(ifile,'("node",a)')node_option
  write(ifile,'("edge",a)')edge_option

  flux(:) = log(flux(:))


  ! Write out the node list.
  !
  allocate( found(NMOLMAX) )

  found(1:NMOLMAX) = .false.
  do i = 1, nrxn
     if (gvconnect(i)) then

        if (rxntype(i) == 'B') then
           ii = rid(1,i)
           jj = rid(2,i)

           if (.not.found(ii)) then
              write(line,*)ii
              write(ifile,'("S",a,a)')trim( adjustl(line) ),node_option
              found(ii) = .true.
           endif

           if (.not.found(jj)) then
              write(line,*)jj
              write(ifile,'("S",a,a)')trim( adjustl(line) ),node_option
              found(jj) = .true.
           endif

           if (prodtype(i) == 'U') then
              ii = pid(1,i)

              if (.not.found(ii)) then
                 write(line,*)ii
                 write(ifile,'("S",a,a)')trim( adjustl(line) ),node_option
                 found(ii) = .true.
              endif
           else if (prodtype(i) == 'B') then
              ii = pid(1,i)
              jj = pid(2,i)

              if (.not.found(ii)) then
                 write(line,*)ii
                 write(ifile,'("S",a,a)')trim( adjustl(line) ),node_option
                 found(ii) = .true.
              endif

              if (.not.found(jj)) then
                 write(line,*)jj
                 write(ifile,'("S",a,a)')trim( adjustl(line) ),node_option
                 found(jj) = .true.
              endif
           endif

        else if (rxntype(i) == 'U') then

           ii = rid(1,i)

           if (.not.found(ii)) then
              write(line,*)ii
              write(ifile,'("S",a,a)')trim( adjustl(line) ),node_option
              found(ii) = .true.
           endif

           if (prodtype(i) == 'U') then
              ii = pid(1,i)

              if (.not.found(ii)) then
                 write(line,*)ii
                 write(ifile,'("S",a,a)')trim( adjustl(line) ),node_option
                 found(ii) = .true.
              endif
           else if (prodtype(i) == 'B') then
              ii = pid(1,i)
              jj = pid(2,i)

              if (.not.found(ii)) then
                 write(line,*)ii
                 write(ifile,'("S",a,a)')trim( adjustl(line) ),node_option
                 found(ii) = .true.
              endif

              if (.not.found(jj)) then
                 write(line,*)jj
                 write(ifile,'("S",a,a)')trim( adjustl(line) ),node_option
                 found(jj) = .true.
              endif
           endif
        endif
     endif
  enddo


  ! Loop over reactions.
  !
  do i = 1, nrxn

     if (gvconnect(i)) then

        if (rxntype(i) == 'B') then

           ii = rid(1,i)
           jj = rid(2,i)

           if (prodtype(i) == 'U') then
              kk = pid(1,i)

              if (abs( Gts(i) ) < SMALL) then  ! BARRIERLESS REACTION !

                 write(line,*)ii
                 write(line2,*)jj
                 write(line3,*)i
                 write(line4,*)kk

                 write(ifile,'("B",a,a)')trim( adjustl(line3) ),B_option
                 write(ifile,'("S",a,"--B",a,a,";")')trim(adjustl(line)),trim(adjustl(line3)),edge_option
                 write(ifile,'("S",a,"--B",a,a,";")')trim(adjustl(line2)),trim(adjustl(line3)),edge_option
                 write(ifile,'("B",a,"--S",a,a,";")')trim(adjustl(line3)),trim(adjustl(line4)),edge_option

              else

                 write(line,*)ii
                 write(line2,*)jj
                 write(line3,*)i
                 write(line4,*)kk

                 if (.not.dynamic_graph) then
                    write(ifile,'("TS",a,a)')trim( adjustl(line3) ),TS_option
                    write(ifile,'("S",a,"--TS",a,a,";")')trim(adjustl(line)),trim(adjustl(line3)),edge_option
                    write(ifile,'("S",a,"--TS",a,a,";")')trim(adjustl(line2)),trim(adjustl(line3)),edge_option
                    write(ifile,'("TS",a,"--S",a,a,";")')trim(adjustl(line3)),trim(adjustl(line4)),edge_option
                 else
                    weight = flux(i)
                    write(line3,*)weight
                    write(ifile,'("S",a,"--S",a,2x,"[penwidth=",f15.10,"];")')trim(adjustl(line)), &
                         trim(adjustl(line2)),flux(i)*scale
                 endif
              endif

           else if (prodtype(i) == 'B') then

              kk = pid(1,i)
              ll = pid(2,i)

              if (abs( Gts(i) ) < SMALL) then  ! BARRIERLESS REACTION !

                 write(line,*)ii
                 write(line2,*)jj
                 write(line3,*)i
                 write(line4,*)kk
                 write(line5,*)ll

                 write(ifile,'("B",a,a)')trim( adjustl(line3) ),B_option
                 write(ifile,'("S",a,"--B",a,a,";")')trim(adjustl(line)),trim(adjustl(line3)),edge_option
                 write(ifile,'("S",a,"--B",a,a,";")')trim(adjustl(line2)),trim(adjustl(line3)),edge_option
                 write(ifile,'("B",a,"--S",a,a,";")')trim(adjustl(line3)),trim(adjustl(line4)),edge_option
                 write(ifile,'("B",a,"--S",a,a,";")')trim(adjustl(line3)),trim(adjustl(line5)),edge_option

              else

                 write(line,*)ii
                 write(line2,*)jj
                 write(line3,*)i
                 write(line4,*)kk
                 write(line5,*)ll

                 if (.not.dynamic_graph) then
                    write(ifile,'("TS",a,a)')trim( adjustl(line3) ),TS_option
                    write(ifile,'("S",a,"--TS",a,a,";")')trim(adjustl(line)),trim(adjustl(line3)),edge_option
                    write(ifile,'("S",a,"--TS",a,a,";")')trim(adjustl(line2)),trim(adjustl(line3)),edge_option
                    write(ifile,'("TS",a,"--S",a,a,";")')trim(adjustl(line3)),trim(adjustl(line4)),edge_option
                    write(ifile,'("TS",a,"--S",a,a,";")')trim(adjustl(line3)),trim(adjustl(line5)),edge_option
                 else
                    weight = flux(i)
                    write(line3,*)weight
                    write(ifile,'("S",a,"--S",a,2x,"[penwidth=",f15.10,"];")')trim(adjustl(line)), &
                         trim(adjustl(line2)),flux(i)*scale
                 endif
              endif
           endif

        else if (rxntype(i) == 'U') then
           ii = rid(1,i)

           if (prodtype(i) == 'U') then
              kk = pid(1,i)

              write(line,*)ii
              write(line2,*)kk
              write(line3,*)i

              if (.not.dynamic_graph) then
                 write(ifile,'("TS",a,a)')trim( adjustl(line3) ),TS_option
                 write(ifile,'("S",a,"--TS",a,a,";")')trim(adjustl(line)),trim(adjustl(line3)),edge_option
                 write(ifile,'("TS",a,"--S",a,a,";")')trim(adjustl(line3)),trim(adjustl(line2)),edge_option
              else
                 weight = flux(i)
                 write(line3,*)weight
                 write(ifile,'("S",a,"--S",a,2x,"[penwidth=",f15.10,"];")')trim(adjustl(line)), &
                      trim(adjustl(line2)),flux(i)*scale
              endif

           else if (prodtype(i) == 'B') then

              kk = pid(1,i)
              ll = pid(2,i)

              if (abs( Gts(i) ) < SMALL) then  ! BARRIERLESS REACTION !

                 write(line,*)ii
                 write(line3,*)i
                 write(line4,*)kk
                 write(line5,*)ll

                 write(ifile,'("B",a,a)')trim( adjustl(line3) ),B_option
                 write(ifile,'("S",a,"--B",a,a,";")')trim(adjustl(line)),trim(adjustl(line3)),edge_option
                 write(ifile,'("B",a,"--S",a,a,";")')trim(adjustl(line3)),trim(adjustl(line4)),edge_option
                 write(ifile,'("B",a,"--S",a,a,";")')trim(adjustl(line3)),trim(adjustl(line5)),edge_option

              else

                 write(line,*)ii
                 write(line3,*)i
                 write(line4,*)kk
                 write(line5,*)ll

                 if (.not.dynamic_graph) then
                    write(ifile,'("TS",a,a)')trim( adjustl(line3) ),TS_option
                    write(ifile,'("S",a,"--TS",a,a,";")')trim(adjustl(line)),trim(adjustl(line3)),edge_option
                    write(ifile,'("TS",a,"--S",a,a,";")')trim(adjustl(line3)),trim(adjustl(line4)),edge_option
                    write(ifile,'("TS",a,"--S",a,a,";")')trim(adjustl(line3)),trim(adjustl(line5)),edge_option
                 else
                    weight = flux(i)
                    write(line3,*)weight
                    write(ifile,'("S",a,"--S",a,2x,"[penwidth=",f15.10,"];")')trim(adjustl(line)), &
                         trim(adjustl(line2)),flux(i)*scale
                 endif
              endif

           endif

        endif

     endif

  enddo

  write(ifile,'("}")')


  return
end Subroutine GraphViz


!
!*****************************************************************************************
!
! SUBROUTINE: READINPUTS
!
! Read the input file information and the main input file.
! Also, opens a log file and writes initial information.
!
!*****************************************************************************************
!
Subroutine ReadInputs(Method,Tstop,Temperature,irun,volume,ntmax,Kdiff,timestep, &
     filename,c_date,c_time,log_file,GTStype,nrun,readrates,usecovar)
  implicit none
  integer :: i, irun, ntmax, log_file
  integer :: narg, nrun
  character*3 :: Method
  real(8) :: Tstop, Temperature, volume, Kdiff, timestep
  character (len=50)  :: filename(4)
  character (len=8)   :: c_date
  character (len=10)  :: c_time
  character (len=3)   :: GTStype
  logical :: readrates, usecovar
  namelist /input/ Method, Tstop, Temperature, irun, volume, ntmax, &
       Kdiff, timestep,GTStype,nrun,readrates,usecovar


  ! Read input files from command line.
  !
  narg = iargc()
  if ( narg .ne. 4 ) then
     write(6,'(/"* Program usage: kinetics.x [input file] [molecules file] [reactions file] [initial state]"/)')
     stop
  endif
  do i = 1, 4
     call getarg (i, filename(i) )
  enddo


  ! Open the input file and read the input namelist.
  !
  Open(10, file = filename(1), status = 'unknown')
  read(10,input)
  close(10)


  ! Error check input parameters.
  !
  if (Method.ne.'DIR'.and.Method.ne.'SSA') then
     stop '* ERROR: Method should be DIR or SSA'
  endif



  ! Open the log file.
  !
  Open(log_file, file = trim( filename(1) )//'.log', status = 'unknown')
  write(log_file,'(/"*** Kinetics simulation code ***"/)')
  write(log_file,'("# Start date: ",3x,a2,"/",a2,"/",a4)')c_date(7:8),c_date(5:6),c_date(1:4)
  write(log_file,'("# Start time: ",3x,a2,":",a2,":",a2/)')c_time(1:2),c_time(3:4),c_time(5:6)

  write(log_file,'("# Input file:",3x,a50)')adjustl( filename(1) )
  write(log_file,'("# Molecules file:",3x,a50)')adjustl( filename(2) )
  write(log_file,'("# Reactions file:",3x,a50)')adjustl( filename(3) )
  write(log_file,'("# Initial state file:",3x,a50)')adjustl( filename(4) )

  write(log_file,'(/"*** Input parameters ***"/)')
  write(log_file,'("* Method: ",1x,a3)')
  write(log_file,'("* Stopping time =",1x,f14.8,1x,"s")')Tstop
  write(log_file,'("* Temperature =",1x,f14.8,1x,"K")')Temperature
  write(log_file,'("* Random number seed =",1x,i5)')irun
  write(log_file,'("* Volume =",1x,f14.8,1x,"m^3")')volume
  write(log_file,'("* Maximum number of time-steps =",1x,i5)')ntmax
  write(log_file,'("* Diffusion-limited rate =",1x,f14.8,1x,"s^-1")')Kdiff
  write(log_file,'("* Timestep =",1x,f14.8,1x,"s")')timestep

  if (readrates) then
    write(log_file,'("* Reading RATES from rxn file.....")')
  else
    write(log_file,'("* Reading dG from rxn file...calculating TST rates separately....")')
  endif

  write(log_file,*) ; call flush(log_file)


  ! Convert parameters.
  !
  Kdiff = Kdiff * 1d-12        ! To ps-1
  timestep = timestep * 1d+12  ! To ps.

  return
end Subroutine ReadInputs



!
!*****************************************************************************************
!
! SUBROUTINE: READMOLFILE
!
! Reads the set of molecules comprising the network from filename(2).
!
!*****************************************************************************************
!
Subroutine ReadMolFile(log_file,filename,nmol,idtag,molid,Gmol)
  implicit none
  include 'globals.h'
  integer :: log_file, nmol, status
  integer :: idtag(NMOLMAX)
  real(8) :: Gmol(NMOLMAX)
  character (len=30)  :: molid(NMOLMAX)
  character (len=50)  :: filename(4)
  character (len=90)  :: line, comment


  ! Write report to log file.
  !
  write(log_file,'(/"*** Reading molecules file ***"/)')

  open(10, file = filename(2), status='unknown')
  nmol = 0
  do
     read(10,'(a100)',iostat=status)line

     if (status /=0)exit

     if (line(1:1) == '#') then  ! Comment lines begin with #

        read(line,'(a100)')comment
        write(log_file,'(":: COMMENT ::",1x,a70)')comment

     else                        ! Otherwise assume a molecule found.

        nmol = nmol + 1
        read(line,*)idtag(nmol),molid(nmol), Gmol(nmol)
        write(log_file,'(/"*** MOLECULE:",1x,i4,1x,"***")')nmol
        write(log_file,'(" ID number =",3x,i6)')idtag(nmol)
        write(log_file,'(" Label =",3x,a30)')molid(nmol)
        write(log_file,'(" Gmol =",3x,f14.8,1x,"Hartrees")')Gmol(nmol)

     endif

  enddo
  write(log_file,'(/"* TOTAL NUMBER OF MOLECULES =",1x,i4)')nmol
  write(log_file,'(/" Finished reading molecules file...."/)')
  close(unit = 10)

  return
end Subroutine ReadMolFile



!
!*****************************************************************************************
!
! SUBROUTINE: READRXNFILE
!
! Reads the set of reactions comprising the network from filename(3).
!
!*****************************************************************************************
!

Subroutine ReadRxnFile(log_file,filename,nrxn,rxnid,rxntype,prodtype,gvconnect,rid,pid, &
     Gts,Gmol,readrates,ktst)
  implicit none
  include 'globals.h'
  integer :: log_file, status
  integer :: nrxn, i, j, jj, ii, kk
  integer :: rid(3,NRXNMAX), index_c1,index_c2
  integer :: pid(3,NRXNMAX),index_pp,index_pr,index_rx
  real(8) :: Gts(NRXNMAX), Gmol(NMOLMAX)
  real(8) :: Gf, Gi, ktst(NRXNMAX)
  character (len=1)   :: rxntype(NRXNMAX), prodtype(NRXNMAX)
  character (len=50)  :: filename(4)
  character (len=90)  :: line, comment
  character (len=90)  :: rxnid(NRXNMAX)
  logical :: gvconnect(NRXNMAX), readrates

  ! Write report to log file.
  !
  write(log_file,'(/"*** Reading reactions file ***"/)')

  ! Open reactions file.
  !
  open(10, file = filename(3), status='unknown')
  nrxn = 0
  do
     read(10,'(a100)',iostat=status)line

     if (status /=0)exit

     if (line(1:1) == '#') then  ! Comment lines begin with #

        read(line,'(a100)')comment
        write(log_file,'(":: COMMENT ::",1x,a70)')comment

     else if (line(1:1) == ' ') then

        exit

     else if (line(1:1) == '*') then

         read(line,'(a100)')comment

     else                        ! Otherwise assume a reaction found.

        nrxn = nrxn + 1


        ! We split the line into segments based on :, + and  = delimiters, then
        ! read the relevant parts....
        !
        index_c1 = scan(line,":")
        index_c2 = scan(line,":",.true.)
        index_pr = scan( line, "+" )
        index_pp = scan( line, "+", .true. )
        index_rx = scan( line, "=" )


        ! Up to the first :, read the reaction ID.
        !
        read(line(1:index_c1),*)rxnid(nrxn)
        rxnid(nrxn) = trim(adjustl(rxnid(nrxn)))//'F'


        ! If index_pp == 0 and index_pr == 0, it's unimolecular (there are no + present).
        !
        if (index_pp == 0 .and. index_pr == 0) then
           rxntype(nrxn) = 'U'
           prodtype(nrxn) = 'U'
           read(line(index_c1+1:index_rx-1),*)rid(1,nrxn)
           read(line(index_rx+1:index_c2-1),*)pid(1,nrxn)

        else

           ! if index_pp > index_rx, then the + is on the right and the
           ! products must be bimolecular....
           !
           if (index_pp > index_rx) then
              prodtype(nrxn) = 'B'
              read(line(index_rx+1:index_pp),*)pid(1,nrxn)
              read(line(index_pp+1:index_c2),*)pid(2,nrxn)
           else
              prodtype(nrxn) = 'U'
              read(line(index_rx+1:index_c2),*)pid(1,nrxn)
           endif


           ! If index_pr < index_rx, then the + is on the left and the reactants
           ! must be bimolecular...
           !
           if (index_pr < index_rx) then
              rxntype(nrxn) = 'B'
              read(line(index_c1+1:index_pr-1),*)rid(1,nrxn)
              read(line(index_pr+1:index_rx-1),*)rid(2,nrxn)
           else
              rxntype(nrxn) = 'U'
              read(line(index_c1+1:index_rx-1),*)rid(1,nrxn)
           endif
        endif


        ! Read the rate in s-1 OR the TS energy.
        !
        if (readrates) then
          read(line(index_c2+1:100),*)ktst(nrxn)
          ktst(nrxn) = ktst(nrxn) * 1d-12
        else           ! Read energy in kj/mol
          read(line(index_c2+1:100),*)Gts(nrxn)
        endif


        ! Output to log file.
        !
        write(log_file,'(/"*** REACTION:",1x,i4,1x,"***")')nrxn
        write(log_file,'(" Label =",3x,a30)')rxnid(nrxn)

        if (rxntype(nrxn) == 'U') then
          write(log_file,'(" Reactants type: Unimolecular")')
          write(log_file,'(" Reactants =",3x,i4)')rid(1,nrxn)
          if (readrates) then
            write(log_file,'(" RATE = ",3x,f14.8,1x,"ps-1")')ktst(nrxn)
          else
            write(log_file,'(" Gr =",3x,f14.8,1x,"Hartrees")')Gmol(rid(1,nrxn))
            Gi = (Gmol(rid(1,nrxn))) * tokjmol
          endif
        else if (rxntype(nrxn) == 'B') then
          write(log_file,'(" Reactants type: Bimolecular")')
          write(log_file,'(" Reactants =",3x,i4,3x,i4)')rid(1,nrxn),rid(2,nrxn)
          if (readrates) then
            write(log_file,'(" RATE = ",3x,f14.8,1x,"ps-1")')ktst(nrxn)
          else
            write(log_file,'(" Gr =",3x,f14.8,1x,"Hartrees")')Gmol(rid(1,nrxn))+Gmol(rid(2,nrxn))
            Gi = (Gmol(rid(1,nrxn))+Gmol(rid(2,nrxn))) * tokjmol
          endif
        endif

        ! RELATIVE GTs
        !
        !     if (abs( Gts(nrxn) ) > SMALL) then
        !        Gts(nrxn) = Gts(nrxn)*tokjmol - (Gi)
        !     else
        !        Gts(nrxn) = 0.d0
        !     endif

        if (.not.readrates) then
          if (prodtype(nrxn) == 'U') then
            write(log_file,'(" Product type: Unimolecular")')
            write(log_file,'(" Products =",3x,i4)')pid(1,nrxn)
            write(log_file,'(" Gp =",3x,f14.8,1x,"Hartrees")')Gmol(pid(1,nrxn))
          else if (prodtype(nrxn) == 'B') then
            write(log_file,'(" Product type: Bimolecular")')
            write(log_file,'(" Products =",3x,i4,3x,i4)')pid(1,nrxn),pid(2,nrxn)
            write(log_file,'(" Gp =",3x,f14.8,1x,"Hartrees")')Gmol(pid(1,nrxn))+Gmol(pid(2,nrxn))
          endif
          write(log_file,'(" TS energy in kj/mol=",1x,f14.8)')Gts(nrxn) 
        endif

     endif

  enddo


  ! Now sort out the reverse reactions.
  !
  if (.not.readrates) then
    do i = 1, nrxn

      jj = i + nrxn

      write(6,*)'BACK',jj

      gvconnect(jj) = .false.
      prodtype(jj) = rxntype(i)
      rxnid(jj) = trim(adjustl(rxnid(i)))//'R'
      if (rxntype(i) == 'U') then
        pid(1,jj) = rid(1,i)
        ii = pid(1,jj)
        Gf = Gmol(ii)
      else if (rxntype(i) == 'B') then
        pid(1,jj) = rid(1,i)
        pid(2,jj) = rid(2,i)
        ii = pid(1,jj)
        kk = pid(2,jj)
        Gf = Gmol(ii) + Gmol(kk)
      endif

      rxntype(jj) = prodtype(i)
      if (prodtype(i) == 'U') then
        rid(1,jj) = pid(1,i)
        ii = rid(1,jj)
        Gi = Gmol(ii)
      else if (prodtype(i) == 'B') then
        rid(1,jj) = pid(1,i)
        rid(2,jj) = pid(2,i)
        ii = rid(1,jj)
        kk = rid(2,jj)
        Gi = Gmol(ii) + Gmol(kk)
      endif

      if (Gts(i) > SMALL) then
        Gi = Gi * ToKjmol
        Gf = Gf * toKjmol
        Gts(jj) = Gts(i) - Gi + Gf ! TS in kj/mol 

      else
        Gts(jj) = 0.d0
      endif

      write(log_file,'(/"*** REACTION:",1x,i4,1x,"***")')jj
      write(log_file,'(" Label =",3x,a30)')rxnid(jj)

      if (rxntype(jj) == 'U') then
        write(log_file,'(" Reactants type: Unimolecular")')
        write(log_file,'(" Reactants =",3x,i4)')rid(1,jj)
        write(log_file,'(" Gr =",3x,f14.8,1x,"Hartrees")')Gmol(rid(1,jj))
        Gi = Gmol(rid(1,jj))
      else if (rxntype(jj) == 'B') then
        write(log_file,'(" Reactants type: Bimolecular")')
        write(log_file,'(" Reactants =",3x,i4,3x,i4)')rid(1,jj),rid(2,jj)
        write(log_file,'(" Gr =",3x,f14.8,1x,"Hartrees")')Gmol(rid(1,jj))+Gmol(rid(2,jj))
        Gi = (Gmol(rid(1,jj))+Gmol(rid(2,jj)))
      endif

      if (prodtype(jj) == 'U') then
        write(log_file,'(" Product type: Unimolecular")')
        write(log_file,'(" Products =",3x,i4)')pid(1,jj)
        write(log_file,'(" Gp =",3x,f14.8,1x,"Hartrees")')Gmol(pid(1,jj))
      else if (prodtype(jj) == 'B') then
        write(log_file,'(" Product type: Bimolecular")')
        write(log_file,'(" Products =",3x,i4,3x,i4)')pid(1,jj),pid(2,jj)
        write(log_file,'(" Gp =",3x,f14.8,1x,"Hartrees")')Gmol(pid(1,jj))+Gmol(pid(2,jj))
      endif

      write(log_file,'(" TS energy =",1x,f14.8,1x,"kJ/mol")')Gts(jj)

      ! Sanity check...
      !
      if (Gts(jj) < 0.d0) then
        write(log_file,*)'*** ERROR: NEGATIVE TS ENERGY FOUND ***'
        stop
      endif

    enddo
    nrxn = 2 * nrxn
  endif

  write(log_file,'(/"* TOTAL NUMBER OF REACTIONS =",1x,i4)')nrxn
  write(log_file,'(/" Finished reading reactions file...."/)')
  close(unit = 10)

  return
end Subroutine ReadRxnFile


!
!*****************************************************************************************
!
! SUBROUTINE: READINITIALSTATE
!
! Reads the initial populations and concentrations of the reactants from filename(4).
!
!*****************************************************************************************
!
Subroutine ReadInitialState(log_file,filename,pop,conc,Method,pressure,volume, &
     Temperature,molid,nmol)
  implicit none
  include 'globals.h'
  integer :: log_file, status
  integer :: i, j, id, nmol
  real(8) :: pop(NMOLMAX), conc(NMOLMAX), pressure, volume, Temperature
  character*3 :: Method
  character (len=30)  :: molid(NMOLMAX)
  character (len=50)  :: filename(4)
  character (len=90)  :: line, comment


  ! Zero out the populations and concentrations.
  !
  conc(:) = 0.d0
  pop(:) = 0.d0
  write(log_file,'(/"*** Reading initial state...."/)')
  open(10, file = filename(4), status='unknown')
  do
     read(10,'(a100)',iostat=status)line

     if (status /=0)exit

     if (line(1:1) == '#') then  ! Comment lines begin with #

        read(line,'(a100)')comment
        write(log_file,'(":: COMMENT ::",1x,a70)')comment

     else if (line(1:1) == '*') then

         read(line,'(a100)')comment

     else                        ! Otherwise assume a population found.

        if (Method == 'SSA') then
           read(line,*)id,pop(id)
        else if (Method == 'DIR') then
           read(line,*)id,conc(id)          ! in mol / dm**3
        endif

     endif
  enddo
  close(10)

  write(log_file,'(/"* Initial populations:")')
  do i = 1, nmol
     if (pop(i) > 0) then
        write(log_file,'(/"* Molecule number:",1x,i4)')i
        write(log_file,'("* Molecule ID:",1x,a30)')adjustl( molid(i) )

        if (Method == 'SSA') then
           write(log_file,'("* Initial population =",1x,f14.8)')pop(i)

           pressure = ( (pop(i)/Navogadro) * 8.314d0 * Temperature) / Volume
           pressure = pressure / 1d5
           write(log_file,'("* Initial pressure =",1x,f14.8,1x,"bar")')pressure

           conc(i) = (pop(i)/Navogadro) / (Volume*1d3)
           write(log_file,'("* Concentration = ",1x,f14.8,1x,"mol/dm**3")')conc(i)
        else if (Method == 'DIR') then
           write(log_file,'("* Concentration = ",1x,f14.8,1x,"mol/dm**3")')conc(i)
        endif
     endif
  enddo
  write(log_file,'(/"*** Finished reading initial state...."/)')

  return
end Subroutine ReadInitialState



!
!*****************************************************************************************
!
! SUBROUTINE: GETSTOICHIOMETRYMATRIX
!
! Determines the stoichiometry matrix based on the reaction types and reactant and product
! ID numbers.
!
!*****************************************************************************************
!
Subroutine GetStoichiometryMatrix(stoich,nrxn,rid,pid,rxntype,prodtype)
  implicit none
  include 'globals.h'
  integer :: i, nrxn, ii, jj
  integer :: rid(3,NRXNMAX)
  integer :: pid(3,NRXNMAX)
  real(8) :: stoich(NRXNMAX,NMOLMAX)
  character(len=1) :: rxntype(NRXNMAX), prodtype(NRXNMAX)

  stoich(:,:) = 0
  do i = 1, nrxn

     if (rxntype(i) == 'B') then
        ii = rid(1,i)
        jj = rid(2,i)
        stoich(i,ii) = -1
        stoich(i,jj) = -1
     else if (rxntype(i) == 'U') then
        ii = rid(1,i)
        stoich(i,ii) = -1
     endif

     if (prodtype(i) == 'B') then
        ii = pid(1,i)
        jj = pid(2,i)
        stoich(i,ii) = +1
        stoich(i,jj) = +1
     else if (prodtype(i) == 'U') then
        ii = pid(1,i)
        stoich(i,ii) = 1
     endif

  enddo
  return
end Subroutine GetStoichiometryMatrix





!
!*****************************************************************************************
!
! SUBROUTINE: GETRATES
!
! Calculates the rates ktst(:) for all reactions using TST.
!
! Currently implements equilibrium TST.
!
!*****************************************************************************************
!
Subroutine GetRates(log_file, ktst, rxntype, prodtype, Kdiff, Gmol, Gts, nrxn, &
     pid, rid, Temperature,GTStype)
  implicit none
  include 'globals.h'
  integer :: log_file, nrxn, i, ii, jj
  real(8) :: ktst(NRXNMAX), Kdiff, Gmol(NMOLMAX), Gts(NRXNMAX)
  real(8) :: Temperature, Keq, Gi, Gf, dg
  integer :: rid(3,NRXNMAX)
  integer :: pid(3,NRXNMAX)
  character(len=1) :: rxntype(NRXNMAX), prodtype(NRXNMAX)
  character(len=3) :: Gtstype


  ! Set up the TST rate vector for all reactions. ktst(i) contains the
  ! rate for reaction i. Note that if ktst(i) is less than zero, the
  ! reaction DOES NOT EXIST...
  !
  write(log_file,'("* Calculated TST rates:"/)')

  ktst(:) = -1.d0

  do i = 1, nrxn

     ! BARRIERLESS REACTIONS - signalled by Gts(rxn) = 0.0 (technically, less than SMALL)
     !
     if (abs( Gts(i) ) < SMALL) then

        !stop 'needs correcting for barrierless'

        if (rxntype(i) == 'B') then         ! BIMOLECULAR, BARRIERLESS
           ii = rid(1,i)
           jj = rid(2,i)
           Gi = Gmol(ii) + Gmol(jj)

           if (prodtype(i) == 'U') then
              ii = pid(1,i)
              Gf = Gmol(ii)

              if (Gf < Gi) then

          !       dg = Gf - Gi
          !       dg = dG * ToKjmol * 1000.d0
          !       dG = dG / Navogadro              ! in J
                 ktst(i) = Kdiff

              else if (Gf > Gi) then

                 dg = (Gf - Gi) * tokjmol
                 dg = dG * 1000.d0      ! in J/mol !
                 dG = dG / Navogadro              ! in J
                 Keq = exp(-dG/(kboltz*temperature))

                 ktst(i) = Keq * kdiff

              endif

           else if (prodtype(i) == 'B') then
              ii = pid(1,i)
              jj = pid(2,i)
              Gf = Gmol(ii) + Gmol(jj)
              Stop 'BARRIERLESS RXNS OF TYPE B --> B NOT CODED YET'

           endif

        else if (rxntype(i) == 'U') then    ! UNIMOLECULAR, BARRIERLESS

           ii = rid(1,i)
           Gi = Gmol(ii)

           if (prodtype(i) == 'U') then
              ii = pid(1,i)
              Gf = Gmol(ii)
           else if (prodtype(i) == 'B') then
              ii = pid(1,i)
              jj = pid(2,i)
              Gf = Gmol(ii) + Gmol(jj)
           endif

           if (Gf < Gi) then
              !                 dg = Gf - Gi
              !                 dg = dG * ToKjmol * 1000.d0      ! in J/mol !
              !                 dG = dG / Navogadro              ! in J
              !                 Keq = exp(-dG/(kboltz*temperature))
              !                 ktst(i) = Keq * kdiff

              ktst(i) = kdiff

           else if (Gf > Gi) then

              dg = (Gf - Gi)*tokjmol
              dg = dG * 1000.d0      ! in J/mol !
              dG = dG / Navogadro              ! in J
              Keq = exp(-dG/(kboltz*temperature))

              ktst(i) = Keq * kdiff

           endif

        endif

        write(log_file,'("* Reaction:",1x,i4)')i
        write(log_file,'("* Delta G =",1x,f12.6,1x,"kJ/mol")')dG * Navogadro / 1000.d0
        write(log_file,'("* Rate =",1x,e14.8,1x,"ps^-1"/)')ktst(i)

     ! Normal TST rates...
     !
     else

        if (rxntype(i) == 'B') then
           ii = rid(1,i)
           jj = rid(2,i)
           Gi = Gmol(ii) + Gmol(jj)
        else if (rxntype(i) == 'U') then
           ii = rid(1,i)
           Gi = Gmol(ii)
        endif
        
        ! Convert to kj/mol
        Gi = Gi * tokjmol

        if (GTStype == 'REL') then     ! ASSUME TS ENERGIES GIVEN IN kJ/mol relative to reactants
           dG = Gts(i) * 1000.d0       ! convert to J/mol
           dG = dG / Navogadro         ! in J
        else if (GTStype == 'ABS') then ! ABSOLUTE TS ENERGIES GIVEN
           dG = Gts(i) - Gi                 ! in kj/mol !
           dg = dG * ToKjmol * 1000.d0      ! in J/mol !
           dG = dG / Navogadro              ! in J
        endif

        ktst(i) = (Kboltz * Temperature / h) * exp(-dG / (Kboltz * Temperature) )

        write(log_file,'("* Reaction:",1x,i4)')i
        write(log_file,'("* Delta G =",1x,f12.6,1x,"kJ/mol")')Gts(i) !* Navogadro / 1000.d0

        ktst(i) = ktst(i) * 1d-12 ! Convert to ps-1 !

        ! For a bimolecular reaction, limit rate to the diffusion-limited value Kdiff.
        !
        if (rxntype(i) == 'B') then
            ktst(i) = min(ktst(i),Kdiff)
        endif
!!$           !  ktst(i) = max(ktst(i), Kdiff)
!!$
!!$           if (ktst(i) > kdiff) then
!!$              !              ktst(i) = kdiff
!!$              if (rid(1,i) == 2 .or. rid(2,i) == 2) then
!!$                 !       ktst(i) = ktst(i) * conc(2)
!!$                 !                    ktst(i) = 1.29d9
!!$                 !       ktst(i) = ktst(i) / 34.5
!!$              else if (rid(1,i) == 3 .or. rid(2,i) == 3) then
!!$                 !       ktst(i) = ktst(i) * conc(3)
!!$                 !                    ktst(i) = 1.29d9
!!$                 !      ktst(i) = ktst(i) / 34.5
!!$              endif
!!$
!!$           endif


           !?????
     endif

     write(log_file,'("* Rate =",1x,e14.8,1x,"ps^-1"/)')ktst(i)

  enddo
  write(log_file,'(/"*** Finished calculating rates...."/)')

  return
end Subroutine GetRates



!
!*****************************************************************************************
!
! SUBROUTINE: STOCHASTIC
!
! Implements Gillespie's stochastic simulation algorithm to propagate the
! species populations.
!
!*****************************************************************************************
!
Subroutine Stochastic(Tstop,pop,ktst,rid,pid,rxntype,prodtype,volume,nrxn,nmol,irun, &
     stoich, flux, ir, popsave,Tsave,log_file,fluxtot,readout)

  implicit none
  include 'globals.h'

  integer :: i,j,jj,irxn,irun, ifile, ii, ir
  integer :: rid(3,NRXNMAX), log_file, icounter
  integer :: pid(3,NRXNMAX),nrxn,nmol
  real(8) :: stoich(NRXNMAX,NMOLMAX), Tstop
  real(8) :: Tout(NOUT), flux(NRXNMAX), Sflux(NMOLMAX)
  real(8) :: time, fluxtot, propens(NRXNMAX), ran2
  real(8) :: ktst(NRXNMAX), volume
  real(8) :: pop(NMOLMAX),sum,sumprop,s1,s2,tau
  real(8) :: popsave(NMOLMAX,NRUNMAX,NOUT),Tsave(NRUNMAX,NOUT)
  logical :: Readout(NOUT)
  character(len=1) :: rxntype(NRXNMAX), prodtype(NRXNMAX)


  ! Write to log_file
  !
  write(log_file,'(/"*** SSA simulation number: ",1x,i4/)')ir
  call flush(log_file)


  ! Setup the output times.
  !
  do i = 1, NOUT
     Tout(i) = dble(i) * Tstop / dble(NOUT-1)
     Readout(i) = .false.
  enddo


  ! Start the loop over time-steps.
  !
  time = 0.d0
  Sflux(1:nmol) = 0.d0
  flux(1:nrxn) = 0.d0
  fluxtot = 0.d0
  icounter = 0

  do while ( time*1d-12 < Tstop )

     icounter = icounter + 1


     ! Calculate the reaction propensity for each reaction.
     !
     do i = 1, nrxn
        propens(i) = 0.d0
        if (rxntype(i) == 'B') then
           ii = rid(1,i)
           jj = rid(2,i)

           if (ii /= jj) then
             ! propens(i) = 1d-3*pop(ii) * pop(jj) * ktst(i) / (navogadro*Volume)
              propens(i) = pop(ii) * pop(jj) * ktst(i) / (Volume)
           else if (ii == jj) then
            !  propens(i) = 1d-3*(0.5*pop(ii)*(pop(ii)-1)) * ktst(i) / (navogadro*Volume)   ! Note factor of 2 here
              propens(i) = (0.5*pop(ii)*(pop(ii)-1)) * 2.0 * ktst(i) / (Volume)   ! Note factor of 2 here
           endif

        else if (rxntype(i) == 'U') then
           ii = rid(1,i)
           propens(i) = pop(ii) * ktst(i)
        endif
     enddo


     ! Calculate the sum of reaction propensities.
     !
     sumprop = 0.d0
     do i = 1, nrxn
        sumprop = sumprop + propens(i)
     enddo

     ! Calculate two random numbers.
     !
     s1 = ran2(irun,0.d0,1.d0)
     s2 = ran2(irun,0.d0,1.d0)


     ! Calculate time to next reaction.
     !
     tau = (1.d0/sumprop) * log(1.d0 / s1)


     ! Calculate reaction type.
     !
     sum = 0.d0
     do i = 1, nrxn
        sum = sum + propens(i)
        if (sum > s2 * sumprop)then
           irxn = i
           goto 91
        endif
     enddo
91   continue

     flux(irxn) = flux(irxn) + 1.d0
     fluxtot = fluxtot + 1.d0

!     write(6,*)'HERE: ',time,irxn,flux(irxn), fluxtot



     ! Update the time.
     !
     time = time + tau

!     write(log_file,*)'IRXN: ',irxn

     ! Apply the stoichiometric matrix for reaction irxn.
     !
     do i = 1, nmol
        pop(i) = pop(i) + stoich(irxn,i)
        sflux(i) = sflux(i) + stoich(irxn,i)
     enddo


     ! Store the current time and the populations.
     !
     do j = 1, NOUT
        if (.not.Readout(j))then
           if (time * 1d-12 .gt. Tout(j-1).and.time*1d-12.lt.Tout(j)) then
              Tsave(ir,j) = time * 1d-12
              Readout(j) = .true.
              do i = 1, nmol
                 popsave(i,ir,j) = pop(i)
              enddo
           endif
        endif
     enddo


     ! Stop the calculation if we've got to time Tstop.
     !
     if (time*1d-12 >= Tstop) goto 98


     ! Output to log file every so often...
     !
!     if (mod(icounter,10000)==0) then
!        write(log_file,*)'POP: ',time*1d-12,(pop(26)/6.02e23)/(volume*1d3)
!     endif

  enddo
98 continue

!  write(6,*)'WWW: ',fluxtot

  return
end Subroutine Stochastic



!
!*****************************************************************************************
!
! SUBROUTINE: DIRECT
!
! Implements direct algorithm to propagate the species concentrations.
!
!*****************************************************************************************
!
Subroutine Direct(Tstop,conc,ktst,rid,pid,rxntype,prodtype,volume,nrxn,nmol,irun, &
     stoich,flux,timestep)
  implicit none
  include 'globals.h'

  integer :: i,j,jj,irxn,irun, ifile, ii
  integer :: rid(3,NRXNMAX), ic
  integer :: pid(3,NRXNMAX),nrxn,nmol
  real(8) :: stoich(NRXNMAX,NMOLMAX), Tstop, stoichTr(NMOLMAX,NRXNMAX)
  real(8) :: Tout(NOUT), flux(NRXNMAX), Sflux(NMOLMAX)
  real(8) :: time, fluxtot, ran2, timestep
  real(8) :: ktst(NRXNMAX), volume
  real(8) :: pop(NMOLMAX),sum,sumprop,s1,s2,tau
  real(8) :: dcdt(NMOLMAX), conc_store(NMOLMAX), conc(NMOLMAX)
  real(8) :: k1(NMOLMAX), k2(NMOLMAX), k3(NMOLMAX), k4(NMOLMAX)
  logical :: Readout(NOUT)
  character(len=1) :: rxntype(NRXNMAX), prodtype(NRXNMAX)


  ! Setup the output times.
  !
  do i = 1, NOUT
     Tout(i) = dble(i-1) * Tstop / dble(NOUT-1)
     Readout(i) = .false.
  enddo


  ! Form the transposed stoichiometry matrix.
  !
  do i = 1, nmol
     do j = 1, nrxn
        stoichtr(i,j) = stoich(j,i)
     enddo
  enddo

  ! Change from numbers to mol dm-3.
  !
  conc(:) = (conc(:)/ 6.02e23)  /(volume*1d3)


  ! Start the loop over time-steps.
  !
  time = 0.d0
  Sflux(1:nmol) = 0.d0
  flux(1:nrxn) = 0.d0
  fluxtot = 0.d0

  ic = 0
  do while ( time*1d-12 < Tstop )

     ic = ic + 1

     ! Store current concentrations.
     !
     conc_store(:) = conc(:)

     ! RK4 algorithm.
     !
     Call GetDerivs(k1,stoichtr,ktst,conc,rxntype,prodtype,rid,pid,nmol,nrxn)
     conc(:) = conc_store(:) + 0.5 * timestep * k1(:)
     Call GetDerivs(k2,stoichtr,ktst,conc,rxntype,prodtype,rid,pid,nmol,nrxn)
     conc(:) = conc_store(:) + 0.5 * timestep * k2(:)
     Call GetDerivs(k3,stoichtr,ktst,conc,rxntype,prodtype,rid,pid,nmol,nrxn)
     conc(:) = conc_store(:) + timestep * k3(:)
     Call GetDerivs(k4,stoichtr,ktst,conc,rxntype,prodtype,rid,pid,nmol,nrxn)
     dcdt(:) = (1.d0/6.d0) * (k1(:) + 2.0 * k2(:) + 2.0 * k3(:) + k4(:))
     conc(:) = conc(:) + dcdt(:) * timestep

     time = time + timestep

     do i = 1, nmol
        write(6,*)'WTF: ',i,dcdt(i),conc(i)
     enddo


     ! Output the current time and the concentrations.
     !
!     do j = 1, NOUT
!        if (.not.Readout(j))then
!           if (time * 1d-12 .gt. Tout(j-1).and.time*1d-12.lt.Tout(j)) then
!              Readout(j) = .true.
     if (mod(ic,1000) == 0) then
        ifile = 30
        do i = 1, nmol
           write(ifile,*)time * 1d-12, conc(i)
           !                 write(ifile,*)time * 1d-12, conc(i)
           ifile = ifile + 1
        enddo
     endif
!           endif
!        endif
!     enddo


     ! Stop the calculation if we've got to time Tstop.
     !
     if (time*1d-12 >= Tstop) goto 98

  enddo
98 continue

  return
end Subroutine Direct



!
!*****************************************************************************************
!
! SUBROUTINE: GETDERIVS
!
! Calculates the time-derivatives of concentrations for direct simulations.
!
!*****************************************************************************************
!
Subroutine GetDerivs(derivs,stoichtr,ktst,conc,rxntype,prodtype,rid,pid,nmol,nrxn)
  implicit none
  include 'globals.h'

  integer :: i, nrxn, nmol, ii, jj, j
  integer :: rid(3,NRXNMAX)
  integer :: pid(3,NRXNMAX)
  real(8) :: stoichTr(NMOLMAX,NRXNMAX)
  real(8) :: ktst(NRXNMAX), conc(NMOLMAX)
  real(8) :: derivs(NMOLMAX), yvec(NRXNMAX)
  character(len=1) :: rxntype(NRXNMAX), prodtype(NRXNMAX)


  ! Create vector of rates.
  !
  do i = 1, nrxn
     yvec(i) = 0.0

     if (rxntype(i) == 'U') then
        ii = rid(1,i)
        yvec(i) = ktst(i) * conc(ii)
     else if (rxntype(i) == 'B') then
        ii = rid(1,i)
        jj = rid(2,i)
        yvec(i) = ktst(i) * conc(ii) * conc(jj)
     endif

     write(6,*)'NOW: ',rxntype(i),i,yvec(i),conc(ii)

  enddo


  ! Multiplication to get derivs.
  !
  do i = 1, nmol
     derivs(i) = 0.0
     do j = 1, nrxn
        derivs(i) = derivs(i) + stoichtr(i,j) * yvec(j)
     enddo
  enddo

  return
end Subroutine GetDerivs


Subroutine SVD_inverse(A, Ainv, n, cond, svdeps)
  implicit none
  integer :: n,LWORK,INFO,i,j
  double precision :: A(n,n), Ainv(n,n), cond,svdeps
  double precision, allocatable :: Astore(:,:), sigma(:), s(:,:),work(:)
  double precision, allocatable :: U(:,:), UT(:,:), V(:,:), VT(:,:)

  LWORK = 5 * n

  allocate ( Astore(n,n) )
  allocate ( sigma(n), s(n,n) )
  allocate ( U(n,n), UT(n,n), V(n,n), VT(n,n) )
  allocate ( work(LWORK) )
  Astore = A

  Call DGESVD( 'A', 'A', n, n, A, n, sigma,U,n,VT,n,work,LWORK,INFO)

  ! Check that everything worked correctly with the matrix inversion.
  if (info.ne.0) then
     write(6,*)'WARNING: Singular matrix detected in DGESVD'
     stop
  endif

  A = Astore

  cond = sigma(1) / sigma(n)

  ! Calculate the inverse matrix, A^-1 = V * (1/sigma) * Ut. Note that we
  ! regularize the inversion by removing singular values less than svdeps
  !
  S(:,:) = 0.0
  do i = 1, n
     if (sigma(i).le.svdeps)then
        S(i,i) = 0.d0
     else
        S(i,i) = 1.d0 / Sigma(i)
     endif
  enddo
  do i = 1, n
     do j = 1, n
        UT(i,j) = U(j,i)
        V(i,j) = VT(j,i)
     enddo
  enddo
  U = matmul(S,UT)
  Ainv = matmul(V,U)

  deallocate ( Astore )
  deallocate (s, sigma)
  deallocate ( u,ut,v,vt)
  deallocate ( work )

  return
end Subroutine SVD_inverse

