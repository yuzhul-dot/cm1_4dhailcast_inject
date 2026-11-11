  MODULE restart_read_module

  implicit none

  private
  public :: restart_read

  logical :: dointerp

  integer, parameter :: infofile = 50
  integer, parameter :: rstfile  = 51
  integer, parameter :: sfile    = 52
  integer, parameter :: ufile    = sfile+1
  integer, parameter :: vfile    = sfile+2
  integer, parameter :: wfile    = sfile+3

  integer, parameter :: in_file_no   =   0
  integer, parameter :: in_file_yes  =   1

  integer :: nfiles
  real :: z1,z2

  logical :: dorstfile,doinfo,dobcast

  integer :: niold,njold,nkold
  integer :: myi1old,myi2old,myj1old,myj2old
  integer :: numprocs_old,info1,info2
  logical :: newdecomp
  integer :: ppnode_old,restart_filetype_old,nioldmax,njoldmax,stag
  integer :: nxold1,nyold1

  integer :: dk1i1,dk1i2,dk1j1,dk1j2
  integer :: dk2i1,dk2i2,dk2j1,dk2j2
  integer :: dk3i1,dk3i2,dk3j1,dk3j2

    ! note: only format 1 (binary) can be done, for now !
    ! (netcdf format will be added in near future) !
  integer, parameter ::    restart_read_format   = 1

  CONTAINS

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine restart_read(nstep,srec,sirec,urec,vrec,wrec,nrec,mrec,prec,      &
                              trecs,trecw,arecs,arecw,                             &
                              nwrite,nwritet,nwritea,nwriteh,nrst,nstatout,        &
                              restart_time,spinup_t2p,                             &
                              num_soil_layers,nrad2d,                              &
                              avgsfcu,avgsfcv,avgsfcs,avgsfcsu,avgsfcsv,avgsfct,avgsfcq,avgsfcp, &
                              dt,dtlast,mtime,ndt,adt,adtlast,acfl,dbldt,mass1,    &
                              stattim,taptim,rsttim,radtim,prcltim,hailtrajtim,    &
                              qbudget,asq,bsq,qname,                               &
                              xfref,xhref,yfref,yhref,xh,xf,yh,yf,zh,zf,sigma,sigmaf,zs, &
                              th0,prs0,pi0,rho0,qv0,u0,v0,                         &
                              rain,raincv,hail,sws,svs,sps,srs,sgs,sus,shs,        &
                              tsk,znt,ust,cd,ch,cq,u1,v1,s1,th1,qv1,thflux,qvflux, &
                              prate,ustt,ut,vt,st,                                 &
                              radbcw,radbce,radbcs,radbcn,                         &
                              rho,prs,ua,dumu,va,dumv,wa,ppi,tha,qa,tkea,          &
                              swten,lwten,thraten,radsw,rnflx,radswnet,radlwin,rad2d,      &
                              effc,effi,effs,effr,effg,effis,                      &
                              thcuten,qvcuten,qccuten,qicuten,ucuten,vcuten,       &
                              lu_index,kpbl2d,psfc,u10,v10,s10,hfx,qfx,xland,      &
                              hpbl,wspd,psim,psih,gz1oz0,br,                       &
                              CHS,CHS2,CQS2,CPMM,ZOL,MAVAIL,                       &
                              MOL,RMOL,REGIME,LH,FLHC,FLQC,QGH,                    &
                              CK,CKA,CDA,USTM,QSFC,T2,Q2,TH2,EMISS,THC,ALBD,       &
                              gsw,glw,chklowq,capg,snowc,fm,fh,mznt,swspd,wstar,delta,tslb,    &
                              tmn,tml,t0ml,hml,h0ml,huml,hvml,tmoml,               &
                              qpten,qtten,qvten,qcten,pta,pdata,ploc,              &
                              haildata,hailloc,ppx,phi1,phi2,                      &
                              pbl3d,pbl2d,pbl2dint, &
                              thpten,qvpten,qcpten,qipten,upten,vpten,qnipten,qncpten, &
                              icenter,jcenter,xcenter,ycenter,domainlocx,domainlocy,adaptmovetim,mvrec,nwritemv,ug,vg, &
                              gamk,kmw,ufw,vfw,u1b,v1b,                                &
                              ntavg,rtavg,tavg,timavg,sfctavg,sfctimavg,               &
                              dum1,dum2,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,restarted,restart_prcl,restart_hail,  &
                              sw31,sw32,se31,se32,ss31,ss32,sn31,sn32,reqs_s,          &
                              confignum)
          ! end_restart_read
      use input
      use lsnudge_module
      use goddard_module, only : consat,consat2
      use lfoice_module, only : lfoice_init
      use bc_module, only : bc2d
#ifdef MPI
      use mpi
      use comm_module, only : comm_1s2d_start,comm_1s2d_end
#endif
#ifdef NETCDF
      use netcdf
      use writeout_nc_module, only : disp_err
#endif

      implicit none

      !----------------------------------------------------------
      ! This subroutine organizes the reading of restart files
      !----------------------------------------------------------

      integer, intent(inout) :: nstep,srec,sirec,urec,vrec,wrec,nrec,mrec,prec,trecs,trecw,arecs,arecw
      integer, intent(inout) :: nwrite,nwritet,nwritea,nwriteh,nrst,nstatout
      double precision, intent(inout) :: restart_time
      logical, intent(inout) :: spinup_t2p
      integer, intent(in) :: num_soil_layers,nrad2d
      double precision, intent(inout) :: avgsfcu,avgsfcv,avgsfcs,avgsfcsu,avgsfcsv,avgsfct,avgsfcq,avgsfcp
      real, intent(inout) :: dt,dtlast
      integer, intent(inout) :: ndt
      double precision, intent(inout) :: adt,adtlast,acfl,dbldt
      double precision, intent(inout) :: mass1
      double precision, intent(inout) :: mtime,stattim,taptim,rsttim,radtim,prcltim,hailtrajtim
      double precision, intent(inout), dimension(nbudget) :: qbudget
      double precision, intent(inout), dimension(numq) :: asq,bsq
      character(len=3), intent(in), dimension(maxq) :: qname
      real, dimension(1-ngxy:nx+ngxy+1) :: xfref,xhref
      real, dimension(1-ngxy:ny+ngxy+1) :: yfref,yhref
      real, intent(in), dimension(ib:ie) :: xh
      real, intent(in), dimension(ib:ie+1) :: xf
      real, intent(in), dimension(jb:je) :: yh
      real, intent(in), dimension(jb:je+1) :: yf
      real, dimension(kb:ke) :: sigma
      real, dimension(kb:ke+1) :: sigmaf
      real, intent(in), dimension(ib:ie,jb:je) :: zs
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: zh
      real, intent(in), dimension(ib:ie,jb:je,kb:ke+1) :: zf
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: th0,prs0,pi0,rho0,qv0
      real, intent(inout), dimension(ib:ie+1,jb:je,kb:ke) :: u0
      real, intent(inout), dimension(ib:ie,jb:je+1,kb:ke) :: v0
      real, intent(inout), dimension(ib:ie,jb:je,nrain) :: rain,raincv,hail,sws,svs,sps,srs,sgs,sus,shs
      real, intent(inout), dimension(ib:ie,jb:je) :: tsk,znt,ust,cd,ch,cq,u1,v1,s1,th1,qv1,xland,psfc,thflux,qvflux,prate,ustt,ut,vt,st
      real, intent(inout), dimension(jb:je,kb:ke) :: radbcw,radbce
      real, intent(inout), dimension(ib:ie,kb:ke) :: radbcs,radbcn
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: rho,prs
      real, intent(inout), dimension(ib:ie+1,jb:je,kb:ke) :: ua
      real, intent(inout), dimension(ib:ie+1,jb:je,kb:ke) :: dumu
      real, intent(inout), dimension(ib:ie,jb:je+1,kb:ke) :: va
      real, intent(inout), dimension(ib:ie,jb:je+1,kb:ke) :: dumv
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke+1) :: wa
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: ppi,tha
      real, intent(inout), dimension(ibm:iem,jbm:jem,kbm:kem,numq) :: qa
      real, intent(inout), dimension(ibt:iet,jbt:jet,kbt:ket) :: tkea
      real, intent(inout), dimension(ibr:ier,jbr:jer,kbr:ker) :: swten,lwten,thraten,effc,effi,effs,effr,effg,effis
      real, intent(inout), dimension(ni,nj) :: radsw,rnflx,radswnet,radlwin
      real, intent(inout), dimension(ni,nj,nrad2d) :: rad2d
      real, intent(inout), dimension(ibcu:iecu,jbcu:jecu,kbcu:kecu) :: thcuten,qvcuten,qccuten,qicuten,ucuten,vcuten
      integer, intent(inout), dimension(ibl:iel,jbl:jel) :: lu_index
      integer, intent(inout), dimension(ibl:iel,jbl:jel) :: kpbl2d
      real, intent(inout), dimension(ibl:iel,jbl:jel) :: u10,v10,s10,hfx,qfx, &
                                      hpbl,wspd,psim,psih,gz1oz0,br,          &
                                      CHS,CHS2,CQS2,CPMM,ZOL,MAVAIL,          &
                                      MOL,RMOL,REGIME,LH,FLHC,FLQC,QGH,   &
                                      CK,CKA,CDA,USTM,QSFC,T2,Q2,TH2,EMISS,THC,ALBD,   &
                                      gsw,glw,chklowq,capg,snowc,fm,fh,mznt,swspd,wstar,delta
      real, intent(inout), dimension(ibl:iel,jbl:jel,num_soil_layers) :: tslb
      real, intent(inout), dimension(ibl:iel,jbl:jel) :: tmn,tml,t0ml,hml,h0ml,huml,hvml,tmoml
      real, intent(inout), dimension(ibm:iem,jbm:jem,kbm:kem) :: qpten,qtten,qvten,qcten
      real, intent(inout), dimension(ibp:iep,jbp:jep,kbp:kep,npt) :: pta
      real, intent(inout), dimension(nparcels,npvals) :: pdata
      real, intent(inout), dimension(nparcels,3) :: ploc
      real, intent(inout), dimension(nhailtrajs,nhvals) :: haildata
      real, intent(inout), dimension(nhailtrajs,13) :: hailloc
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: ppx
      real, intent(inout), dimension(ibph:ieph,jbph:jeph,kbph:keph) :: phi1,phi2
      real, intent(inout), dimension(ibpbl:iepbl,jbpbl:jepbl,kbpbl:kepbl,npbl3d) :: pbl3d
      real, intent(inout), dimension(ibpbl:iepbl,jbpbl:jepbl,npbl2d) :: pbl2d
      integer, intent(inout), dimension(ibpbl:iepbl,jbpbl:jepbl,npbl2dint) :: pbl2dint
      real, intent(inout), dimension(ibb:ieb,jbb:jeb,kbb:keb) :: thpten,qvpten,qcpten,qipten,upten,vpten,qnipten,qncpten
      integer, intent(inout) :: icenter,jcenter
      real, intent(inout) :: xcenter,ycenter
      double precision, intent(inout) :: domainlocx,domainlocy
      double precision, intent(inout) :: adaptmovetim
      integer, intent(inout) :: mvrec,nwritemv
      real, intent(inout), dimension(kb:ke) :: ug,vg
      real, intent(inout), dimension(kb:ke) :: gamk,kmw,ufw,vfw,u1b,v1b
      integer, intent(inout), dimension(ntim) :: ntavg
      double precision, intent(inout), dimension(ntim) :: rtavg
      real, intent(inout), dimension(ibta:ieta,jbta:jeta,kbta:keta,ntim,ntavr) :: tavg
      real, intent(inout), dimension(ibta:ieta,jbta:jeta,kbta:keta,ntavr) :: timavg
      real, intent(inout), dimension(ibta:ieta,jbta:jeta,ntim,nsfctavr) :: sfctavg
      real, intent(inout), dimension(ibta:ieta,jbta:jeta,nsfctavr) :: sfctimavg
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: dum1,dum2
      real, intent(inout), dimension(d3i,d3j) :: dat1
      real, intent(inout), dimension(d2i,d2j) :: dat2
      real, intent(inout), dimension(d3i,d3j,d3n) :: dat3
      integer, intent(inout), dimension(d3t) :: reqt
      integer, intent(in), dimension(0:numprocs-1) :: myi1p,myi2p,myj1p,myj2p
      logical, intent(in) :: restarted
      logical, intent(inout) :: restart_prcl
      logical, intent(inout) :: restart_hail
      real, intent(inout), dimension(cmp,jmp,kmp) :: sw31,sw32,se31,se32
      real, intent(inout), dimension(imp,cmp,kmp) :: ss31,ss32,sn31,sn32
      integer, intent(inout), dimension(rmp) :: reqs_s
      integer, intent(inout) :: confignum

      character(len=maxstring) fname
      character(len=8) :: text1
      character(len=6) :: aname
      integer :: i,j,k,kk,n,nn,n2,nfile,np,nt,nvar,nread,orecs,orecu,orecv,orecw,ndum,ifoo,jfoo
      real :: w1,w2,w3,w4,w5,w6,w7,w8
      integer :: ncid,time_index,old_format,old_filetype,old_ptype,old_ipbl,old_sfcmodel,old_sgsmodel
      integer :: nbudgetfoo,numqfoo
      double precision, dimension(nbudget,0:numprocs-1) :: sbudget
      double precision, dimension(numq,0:numprocs-1) :: csq,dsq
      real, dimension(:,:), allocatable :: pfoo
      real, dimension(:,:), allocatable :: hfoo
      real, dimension(:), allocatable :: dumz
      real, dimension(:,:), allocatable :: dumx,dumy
#ifdef MPI
      integer :: proc,index,count,req1,req2,req3,reqp,ierr
#endif
#ifdef NETCDF
      integer :: varid,ncstatus
#endif
      integer :: nxold,nyold,nzold,ngxyold,iold,jold,kold,nxr,nyr,nold
      integer :: wbcold,ebcold,sbcold,nbcold,bbcold,tbcold
      integer :: oorecs,oorecu,oorecv,oorecw
      real :: maxxold,maxyold,maxzold
      integer :: ntim_old,ntavr_old,nsfctavr_old,keta_old,confignum_old
      real, dimension(:), allocatable :: xhref0,xfref0,yhref0,yfref0,sigma0,sigmaf0
      real, dimension(:,:), allocatable :: datk1,datk2,datk3
      integer, dimension(:), allocatable :: ix,iy,iz,ixu,iyv,izw
      real, dimension(:), allocatable :: rx,ry,rz,rxu,ryv,rzw
      double precision, dimension(:), allocatable :: dptmp
      integer :: in_file
      logical :: openrst,opens,openu,openv,openw,openinfo
      integer :: open_nifiles,open_nxfiles,open_sfiles,open_ufiles,open_vfiles,open_wfiles
      real :: cm1rversion_old
      integer, dimension(:,:), allocatable :: info
      integer :: nodes_old
      integer :: nloop,nloop2,i1old,i2old,j1old,j2old
      integer :: iostat,filenametype,restart_read_filetype
      logical :: doit
      integer :: utav_old,vtav_old,wtav_old,ttav_old,qtav_old,etav_old,uutav_old,vvtav_old

!-----------------------------------------------------------------------

  openrst  = .false.
  opens    = .false.
  openu    = .false.
  openv    = .false.
  openw    = .false.
  openinfo = .false.
  open_nifiles = 0
  open_nxfiles = 0
  open_sfiles = 0
  open_ufiles = 0
  open_vfiles = 0
  open_wfiles = 0
  if(myid.eq.0)then
    doinfo = .true.
    dorstfile = .true.
  else
    doinfo = .false.
    dorstfile = .false.
  endif
  dobcast = .true.
  newdecomp = .false.

    if(dowr) write(outfile,*)
    if(dowr) write(outfile,*) 'Inside RESTART_READ'
    if(dowr) write(outfile,*)

  rformat:  &
  IF( restart_read_format.eq.1 )THEN
    ! unformatted direct-access (grads) format:

    do i=1,maxstring
      fname(i:i) = ' '
    enddo

    doit = .true.
    nloop = 0

101 format(i6.6)

  myid0:  &
  IF( myid.eq.0 )THEN

    fileloop:  &
    do while( doit )

      nloop = nloop+1

      if(     nloop.eq.1 )then
        fname = 'cm1rst_i.dat'
      elseif( nloop.eq.2 )then
        fname = 'cm1rst_tXXXXXX_i.dat'
        write(fname( 9:14),101) rstnum
      elseif( nloop.eq.3 )then
        fname = 'cm1rst_nXXXXXX_tYYYYYY_i.dat'
        write(fname( 9:14),101) 0
        write(fname(17:22),101) rstnum
      elseif( nloop.eq.4 )then
        fname = 'cm1rst_pXXXXXX_tYYYYYY_i.dat'
        write(fname( 9:14),101) 0
        write(fname(17:22),101) rstnum
      else
        print *
        print *,'  Could not find restart files '
        print *
        print *,'      stopping cm1 '
        print *
        call stopcm1
      endif
      print *
      print *,'  trying fname = ',trim(fname)
      open(unit=infofile,file=fname,form='unformatted',status='old',iostat=iostat)
      print *,'  iostat = ',iostat

      if( iostat.eq.0 )then
        print *,'  ... success ... '
        print *
        openinfo = .true.
      endif

      if( openinfo )then
        if(     nloop.eq.1 )then
          fname = 'cm1rst_x.dat'
        elseif( nloop.eq.2 )then
          fname = 'cm1rst_tXXXXXX_x.dat'
          write(fname( 9:14),101) rstnum
        elseif( nloop.eq.3 )then
          fname = 'cm1rst_nXXXXXX_tYYYYYY_x.dat'
          write(fname( 9:14),101) 0
          write(fname(17:22),101) rstnum
        elseif( nloop.eq.4 )then
          fname = 'cm1rst_pXXXXXX_tYYYYYY_x.dat'
          write(fname( 9:14),101) 0
          write(fname(17:22),101) rstnum
        else
          print *
          print *,'  Could not find restart files '
          print *
          print *,'      stopping cm1 '
          print *
          call stopcm1
        endif
        print *
        print *,'  trying fname = ',trim(fname)
        open(unit=rstfile,file=fname,form='unformatted',status='old',iostat=iostat)
        print *,'  iostat = ',iostat
        if( iostat.eq.0 )then
          print *,'  ... success ... '
          print *
          openrst = .true.
        endif
      endif

      if( openinfo .and. openrst )then
        print *,'  found info files; proceeding ... '
        print *
        filenametype = nloop
        doit = .false.
      else
        print *,'  did not find info files; try other file names ... '
        print *
        openinfo = .false.
        openrst = .false.
      endif

      print *

    enddo  fileloop

  ENDIF  myid0

    if(myid.eq.0) print *
    if(myid.eq.0) print *,'  filenametype = ',filenametype
    if(myid.eq.0) print *

  ENDIF  rformat

!---------------------------------------------------------------
! metadata:

  interp_on_restart = .false.
  dointerp = .false.
  nxold = 1
  nyold = 1
  nzold = 1
  ngxyold = 1
  orecs = 1
  orecu = 1
  orecv = 1
  orecw = 1
  oorecs = 1
  oorecu = 1
  oorecv = 1
  oorecw = 1

  IF( restart_read_format.eq.1 )THEN

    if(myid.eq.0)then
      read(rstfile,err=7001) cm1rversion_old
      if( cm1rversion_old.lt.21.999 ) go to 7001
      read(rstfile) confignum_old
      read(rstfile) nxold
      read(rstfile) nyold
      read(rstfile) nzold
      read(rstfile) ngxyold
      if(myid.eq.0) print *,'  nxold,nyold,nzold = ',nxold,nyold,nzold
      read(rstfile) restart_filetype_old
      read(rstfile) numprocs_old
      read(rstfile) ppnode_old
      read(rstfile) nodes_old
    endif

#ifdef MPI
    call MPI_BCAST(restart_filetype_old,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(numprocs_old        ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(ppnode_old          ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(nodes_old           ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(nxold               ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(nyold               ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(nzold               ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
#endif

    nxold1 = nxold
    nyold1 = nyold
    nkold = nzold
    info1 = 1
    info2 = 1

    restart_read_filetype = restart_filetype_old

  if( restart_filetype_old.eq.3 .or. restart_filetype_old.eq.4 )then
    info1 = 4
    info2 = numprocs_old-1
    if( myid.eq.0 )then
      print *,'  restart_filetype_old = ',restart_filetype_old
      print *,'  numprocs_old,nodes_old,ppnode_old = ',numprocs_old,nodes_old,ppnode_old
      print *
      if( restart_filetype_old.eq.3 )then
        fname = 'cm1rst_pXXXXXX_tYYYYYY_i.dat'
        nfiles = numprocs_old
        nloop2 = 1
      elseif( restart_filetype_old.eq.4 )then
        fname = 'cm1rst_nXXXXXX_tYYYYYY_i.dat'
        nfiles = nodes_old
        nloop2 = ppnode_old
      endif
      write(fname(17:22),101) rstnum
      print *
      print *,'  nfiles,nloop2 = ',nfiles,nloop2
      print *
      open_nifiles = nfiles-1
      do n=1,nfiles-1
        write(fname( 9:14),101) n
        nfile = 1000 + 0*nfiles + n
        open(unit=nfile,file=trim(fname),form='unformatted',status='old')
        print *,'  fname = ',nfile,trim(fname)
      enddo
      print *
      write(fname(24:24),103) 'x'
103   format(a1)
      open_nxfiles = nfiles-1
      do n=1,nfiles-1
        write(fname( 9:14),101) n
        nfile = 1000 + 1*nfiles + n
        open(unit=nfile,file=trim(fname),form='unformatted',status='old')
        print *,'  fname = ',nfile,trim(fname)
      enddo
      print *
    endif
  endif

    allocate( info(info1,0:info2) )
    info = 0

    if( myid.eq.0 )then
      nioldmax = 0
      njoldmax = 0
      if( restart_filetype_old.eq.3 .or. restart_filetype_old.eq.4 )then
        n = 0
        do nf=0,nfiles-1
          if( n.eq.0 )then
            nfile = rstfile
          else
            nfile = 1000 + 1*nfiles + nf
          endif
          do nloop=1,nloop2
            read(nfile) i1old
            read(nfile) i2old
            read(nfile) j1old
            read(nfile) j2old
            info(1,n) = i1old
            info(2,n) = j1old
            info(3,n) = i2old-i1old+1
            info(4,n) = j2old-j1old+1
            write(outfile,121) n,info(1,n),info(2,n),info(3,n),info(4,n)
121         format('  info = ',5(1x,i6))
!!!            write(outfile,122) n,myi1p(n),myj1p(n),myi2p(n)-myi1p(n)+1,  &
!!!                                                   myj2p(n)-myj1p(n)+1
!!!122         format('         ',5(1x,i6))
            nioldmax = max( nioldmax , info(3,n) )
            njoldmax = max( njoldmax , info(4,n) )
            ! compare subdomain dimensions:
            if( .not. newdecomp )then
              if( info(3,n).ne.(myi2p(n)-myi1p(n)+1) .or.  &
                  info(4,n).ne.(myj2p(n)-myj1p(n)+1) )then
                print *,'    ... different domain decomposition detected ... '
                newdecomp = .true.
              endif
            endif
            n = n+1
          enddo
        enddo
      else
        read(rstfile) i1old
        read(rstfile) i2old
        read(rstfile) j1old
        read(rstfile) j2old
      endif
    endif

    myid0only:  &
    IF(myid.eq.0)THEN
      ! only processor 0 has these variables:

      read(rstfile) noutstat0
      read(rstfile) nout2ds0
      read(rstfile) nout3ds0
      read(rstfile) nout3du0
      read(rstfile) nout3dv0
      read(rstfile) nout3dw0
      read(rstfile) nout3dxv0
      read(rstfile) nout3dyv0
      read(rstfile) nout3dzv0
      read(rstfile) noutdiags0
      read(rstfile) noutdiagw0
      read(rstfile) noutazims0
      read(rstfile) noutazimw0
      read(rstfile) noutparcel0
      read(rstfile) nouthifrq0
      read(rstfile) nconfig
      ! ghb: are comms needed for these variables?

      read(rstfile) wbcold
      read(rstfile) ebcold
      read(rstfile) sbcold
      read(rstfile) nbcold
      read(rstfile) bbcold
      read(rstfile) tbcold
      read(rstfile) old_ptype
      read(rstfile) old_ipbl
      read(rstfile) old_sfcmodel
      read(rstfile) old_sgsmodel
      if( old_ptype .ne. ptype )then
        print *
        print *,'  Cannot change ptype on a restart '
        print *
        print *,'      old_ptype = ',old_ptype
        print *,'      ptype     = ',ptype
        print *
        print *,'   stopping model .... '
        print *
        call stopcm1
      endif
      if( old_ipbl .ne. ipbl )then
        init_pbl = 1
        print *
        print *,'  old_ipbl,ipbl = ',old_ipbl,ipbl
        print *
        print *,'  init_pbl = ',init_pbl
        print *
      endif
      if( old_sfcmodel .ne. sfcmodel )then
        init_sfcmodel = 1
        print *
        print *,'  old_sfcmodel,sfcmodel = ',old_sfcmodel,sfcmodel
        print *
        print *,'  init_sfcmodel = ',init_sfcmodel
        print *
      endif
      if( old_sgsmodel .ne. sgsmodel )then
        print *
        print *,'  old_sgsmodel,sgsmodel = ',old_sgsmodel,sgsmodel
        print *
      endif
      if( nxold.ne.nx .or. nyold.ne.ny .or. nzold.ne.nz )then
        interp_on_restart = .true.
        dointerp = .true.
        print *
        print *,'  Different grids detected: '
        print *
        print *,'  nxold,nx = ',nxold,nx
        print *,'  nyold,ny = ',nyold,ny
        print *,'  nzold,nz = ',nzold,nz
        print *
        print *,'  Interpolating to new grid ... '
        print *
      else
        interp_on_restart = .false.
      endif
      if( numprocs .ne. numprocs_old )then
        newdecomp = .true.
        print *
        print *,'  Different domain decomposition detected: '
        print *
        print *,'  nxold,nx  =  ',nxold,nx
        print *,'  nyold,ny  =  ',nyold,ny
        print *,'  nzold,nz  =  ',nzold,nz
        print *
        print *,'  numprocs      =  ',numprocs
        print *,'  numprocs_old  =  ',numprocs_old
        print *
      endif
      if( restart_filetype.eq.4 )then
        if( (nodes.ne.nodes_old) .or. (ppnode.ne.ppnode_old) )then
          newdecomp = .true.
          print *
          print *,'  Different domain decomposition detected: '
          print *
          print *,'  nxold,nx  =  ',nxold,nx
          print *,'  nyold,ny  =  ',nyold,ny
          print *,'  nzold,nz  =  ',nzold,nz
          print *
          print *,'  numprocs      =  ',numprocs
          print *,'  numprocs_old  =  ',numprocs_old
          print *
          print *,'  nodes         =  ',nodes
          print *,'  nodes_old     =  ',nodes_old
          print *
          print *,'  ppnode        =  ',ppnode
          print *,'  ppnode_old    =  ',ppnode_old
          print *
        endif
      endif
      read(rstfile) maxxold
      read(rstfile) maxyold
      read(rstfile) maxzold
      if( interp_on_restart )then
        if( ((maxx-minx)-maxxold).gt.0.1*min_dx .or.  &
            ((maxy-miny)-maxyold).gt.0.1*min_dy .or.  &
            ( maxz      -maxzold).gt.0.1*min_dz )then
          print *
          print *,'  new domain is too large '
          print *
          print *,'  old,new x size = ',maxxold,maxx-minx
          print *,'  old,new y size = ',maxyold,maxy-miny
          print *,'  old,new z size = ',maxzold,maxz

          print *,'  x:',(maxx-minx),maxxold,(maxx-minx)-maxxold
          print *,'  y:',(maxy-miny),maxyold,(maxy-miny)-maxyold
          print *,'  z:',maxz,maxzold,  maxz     -maxzold

          print *
          print *,'  ...  stopping cm1  ... '
          print *
          call stopcm1
        endif
      endif
      read(rstfile) nstep
      read(rstfile) srec
      read(rstfile) sirec
      read(rstfile) urec
      read(rstfile) vrec
      read(rstfile) wrec
      read(rstfile) nrec
      read(rstfile) mrec
      read(rstfile) prec
      read(rstfile) trecs
      read(rstfile) trecw
      read(rstfile) arecs
      read(rstfile) arecw
      read(rstfile) mvrec
      read(rstfile) nwrite
      read(rstfile) nwritet
      read(rstfile) nwritea
      read(rstfile) nwritemv
      read(rstfile) nwriteh
      read(rstfile) nrst
      read(rstfile) nstatout
      read(rstfile) ndt
      read(rstfile) icenter
      read(rstfile) jcenter
      read(rstfile) old_format
      read(rstfile) old_filetype
      read(rstfile) dt
      read(rstfile) dtlast
      read(rstfile) xcenter
      read(rstfile) ycenter
      read(rstfile) umove
      read(rstfile) vmove
      read(rstfile) umovesfc
      read(rstfile) vmovesfc
      read(rstfile) umovem1
      read(rstfile) vmovem1
      read(rstfile) umovem2
      read(rstfile) vmovem2
      read(rstfile) domainlocx
      read(rstfile) domainlocy
      read(rstfile) adaptmovetim
      read(rstfile) cflmax
      read(rstfile) mtime
      read(rstfile) stattim
      read(rstfile) taptim
      read(rstfile) rsttim
      read(rstfile) radtim
      read(rstfile) prcltim
      read(rstfile) hailtrajtim
      read(rstfile) restart_time
      read(rstfile) spinup_t2p
      read(rstfile) adt
      read(rstfile) adtlast
      read(rstfile) acfl
      read(rstfile) dbldt
      read(rstfile) mass1
      read(rstfile) avgsfcu
      read(rstfile) avgsfcv
      read(rstfile) avgsfcs
      read(rstfile) avgsfcsu
      read(rstfile) avgsfcsv
      read(rstfile) avgsfct
      read(rstfile) avgsfcq
      read(rstfile) avgsfcp
      allocate( xhref0(1-ngxyold:nxold+ngxyold+1) )
      xhref0 = 0.0
      allocate( xfref0(1-ngxyold:nxold+ngxyold+1) )
      xfref0 = 0.0
      allocate( yhref0(1-ngxyold:nyold+ngxyold+1) )
      yhref0 = 0.0
      allocate( yfref0(1-ngxyold:nyold+ngxyold+1) )
      yfref0 = 0.0
      allocate( sigma0(0:nzold+1) )
      sigma0 = 0.0
      allocate( sigmaf0(0:nzold+2) )
      sigmaf0 = 0.0
      read(rstfile) xhref0
      read(rstfile) xfref0
      read(rstfile) yhref0
      read(rstfile) yfref0
      read(rstfile) sigma0
      read(rstfile) sigmaf0

      ! ghb: some day, may need to use zh instead of sigma
      z1 = sigma0(1)
      z2 = sigma0(2)

      z1old = sigma0(1)

      print *
      print *,'  Depth of old domain: '
      print *,'  sigma0 = ',sigma0(1),sigma0(nzold)
      print *,'  sigmaf0 = ',sigmaf0(1),sigmaf0(nzold+1)
      print *
      print *,'  z1old = ',z1old
      print *
      print *,'  dointerp = ',dointerp
      print *
    ENDIF  myid0only

#ifdef MPI
      call MPI_BCAST(interp_on_restart,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(dointerp         ,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(newdecomp        ,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(init_pbl         ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(init_sfcmodel    ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(ngxyold          ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(nxold            ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(nyold            ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(nzold            ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(z1old            ,1,MPI_REAL   ,0,MPI_COMM_WORLD,ierr)
#endif
      if( myid.ne.0 )then
        allocate( xhref0(1-ngxyold:nxold+ngxyold+1) )
        xhref0 = 0.0
        allocate( xfref0(1-ngxyold:nxold+ngxyold+1) )
        xfref0 = 0.0
        allocate( yhref0(1-ngxyold:nyold+ngxyold+1) )
        yhref0 = 0.0
        allocate( yfref0(1-ngxyold:nyold+ngxyold+1) )
        yfref0 = 0.0
        allocate( sigma0(0:nzold+1) )
        sigma0 = 0.0
        allocate( sigmaf0(0:nzold+2) )
        sigmaf0 = 0.0
      endif
#ifdef MPI
      call MPI_BCAST(xhref0,nxold+2*ngxyold+1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(xfref0,nxold+2*ngxyold+1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(yhref0,nyold+2*ngxyold+1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(yfref0,nyold+2*ngxyold+1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(sigma0 ,nzold+2,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(sigmaf0,nzold+3,MPI_REAL,0,MPI_COMM_WORLD,ierr)
#endif

      dk1i1 = 0
      dk1j1 = 0
      dk1i2 = 0
      dk1j2 = 0

      dk2i1 = 0
      dk2j1 = 0
      dk2i2 = 0
      dk2j2 = 0

      if( dointerp )then
        dk1i2 = nxold+2
        dk1j2 = nyold+2
        dk2i2 = nxold+2
        dk2j2 = nyold+2
      endif

      dk3i1 = 0
      dk3j1 = 0
      dk3i2 = 0
      dk3j2 = 0

#ifdef MPI
      call MPI_BCAST(nioldmax,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(njoldmax,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
#endif

      if( restart_read_filetype.eq.3 .and. myid.eq.0 )then
        dk3i2 = nioldmax
        dk3j2 = njoldmax
      endif
      if( restart_read_filetype.eq.4 .and. myid.eq.0 )then
        dk3i2 = (nioldmax+1) * ppnode_old
        dk3j2 = (njoldmax+1) * ppnode_old
      endif

      if(myid.eq.0) print *,'  dk1i1,dk1i2,dk1j1,dk1j2 = ',dk1i1,dk1i2,dk1j1,dk1j2
      if(myid.eq.0) print *,'  dk2i1,dk2i2,dk2j1,dk2j2 = ',dk2i1,dk2i2,dk2j1,dk2j2
      if(myid.eq.0) print *,'  dk3i1,dk3i2,dk3j1,dk3j2 = ',dk3i1,dk3i2,dk3j1,dk3j2

      if(myid.eq.0) print *,'  allocating datk1 ... '
      allocate( datk1(dk1i1:dk1i2,dk1j1:dk1j2) )
      datk1 = 0.0
      if(myid.eq.0) print *,'  allocating datk2 ... '
      allocate( datk2(dk2i1:dk2i2,dk2j1:dk2j2) )
      datk2 = 0.0
      if(myid.eq.0) print *,'  allocating datk3 ... '
      allocate( datk3(dk3i1:dk3i2,dk3j1:dk3j2) )
      datk3 = 0.0


      allocate( ix(ni) )
      ix = 0
      allocate( iy(nj) )
      iy = 0
      allocate( iz(nk) )
      iz = 0
      allocate( ixu(ni+1) )
      ixu = 0
      allocate( iyv(nj+1) )
      iyv = 0
      allocate( izw(nk+1) )
      izw = 0

      allocate( rx(ni) )
      rx = 0.0
      allocate( ry(nj) )
      ry = 0.0
      allocate( rz(nk) )
      rz = 0.0
      allocate( rxu(ni+1) )
      rxu = 0.0
      allocate( ryv(nj+1) )
      ryv = 0.0
      allocate( rzw(nk+1) )
      rzw = 0.0

      if(myid.eq.0) print *,'  ... done allocating data '
      if(myid.eq.0) print *
      if(myid.eq.0) print *,'    interp_on_restart = ',interp_on_restart
      if(myid.eq.0) print *,'    newdecomp         = ',newdecomp
      if(myid.eq.0) print *

      if( interp_on_restart )then
        !...................................................................
        ! Get interpolation info for scalar points:
        if(myid.eq.0) print *
        if(myid.eq.0) print *,'  xold,xnew,xold,rx: '
        iold = 0
        do i=1,ni
          do while( xhref0(iold+1).lt.xh(i) )
            iold = iold + 1
          enddo
          ix(i) = iold
          rx(i) = sngl( (dble(xh(i))-dble(xhref0(iold)))/(dble(xhref0(iold+1))-dble(xhref0(iold))) )
          if(myid.eq.0) print *,iold,xhref0(iold),xh(i),xhref0(iold+1),rx(i)
        enddo
        if(myid.eq.0) print *
        if(myid.eq.0) print *,'  yold,ynew,yold,ry: '
        jold = 0
        do j=1,nj
          do while( yhref0(jold+1).lt.yh(j) )
            jold = jold + 1
          enddo
          iy(j) = jold
          ry(j) = sngl( (dble(yh(j))-dble(yhref0(jold)))/(dble(yhref0(jold+1))-dble(yhref0(jold))) )
          if(myid.eq.0) print *,jold,yhref0(jold),yh(j),yhref0(jold+1),ry(j)
        enddo
        if(myid.eq.0) print *
        if(myid.eq.0) print *,'  zold,znew,zold,rz: '
        kold = 0
        do k=1,nk
          do while( sigma0(kold+1).lt.sigma(k) )
            kold = kold + 1
          enddo
          iz(k) = kold
          rz(k) = (sigma(k)-sigma0(kold))/(sigma0(kold+1)-sigma0(kold))
          if(myid.eq.0) print *,kold,sigma0(kold),sigma(k),sigma0(kold+1),rz(k)
        enddo
        if(myid.eq.0) print *
        !...................................................................
        ! Get interpolation info for velocity points:
        if(myid.eq.0) print *
        if(myid.eq.0) print *,'  xold,xnew,xold,rx: '
        iold = 0
        do i=1,ni+1
          do while( xfref0(iold+1).lt.xf(i) )
            iold = iold + 1
          enddo
          ixu(i) = iold
          rxu(i) = sngl( (dble(xf(i))-dble(xfref0(iold)))/(dble(xfref0(iold+1))-dble(xfref0(iold))) )
          if(myid.eq.0) print *,iold,xfref0(iold),xf(i),xfref0(iold+1),rxu(i)
        enddo
        if(myid.eq.0) print *
        if(myid.eq.0) print *,'  yold,ynew,yold,ry: '
        jold = 0
        do j=1,nj+1
          do while( yfref0(jold+1).lt.yf(j) )
            jold = jold + 1
          enddo
          iyv(j) = jold
          ryv(j) = sngl( (dble(yf(j))-dble(yfref0(jold)))/(dble(yfref0(jold+1))-dble(yfref0(jold))) )
          if(myid.eq.0) print *,jold,yfref0(jold),yf(j),yfref0(jold+1),ryv(j)
        enddo
        if(myid.eq.0) print *
        if(myid.eq.0) print *,'  zold,znew,zold,rz: '
        kold = 0
        do k=1,nk+1
          do while( sigmaf0(kold+1).lt.sigmaf(k) )
            kold = kold + 1
          enddo
          izw(k) = kold
          rzw(k) = (sigmaf(k)-sigmaf0(kold))/(sigmaf0(kold+1)-sigmaf0(kold))
          if(myid.eq.0) print *,kold,sigmaf0(kold),sigmaf(k),sigmaf0(kold+1),rzw(k)
        enddo
        if(myid.eq.0) print *
        !...................................................................
      endif

#ifdef NETCDF
  ELSEIF( restart_read_format.eq.2 )THEN

    call disp_err( nf90_inq_varid(ncid,"cm1rversion",varid) , .true. )
    call disp_err( nf90_get_var(ncid,varid,cm1rversion,(/time_index/)) , .true. )

    IF(myid.eq.0)THEN

      call disp_err( nf90_inq_varid(ncid,"nstep",varid) , .true. )
      call disp_err( nf90_get_var(ncid,varid,nstep,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"srec",varid) , .true. )
      call disp_err( nf90_get_var(ncid,varid,srec,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"sirec",varid) , .true. )
      call disp_err( nf90_get_var(ncid,varid,sirec,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"urec",varid) , .true. )
      call disp_err( nf90_get_var(ncid,varid,urec,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"vrec",varid) , .true. )
      call disp_err( nf90_get_var(ncid,varid,vrec,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"wrec",varid) , .true. )
      call disp_err( nf90_get_var(ncid,varid,wrec,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"nrec",varid) , .true. )
      call disp_err( nf90_get_var(ncid,varid,nrec,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"mrec",varid) , .true. )
      call disp_err( nf90_get_var(ncid,varid,mrec,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"prec",varid) , .true. )
      call disp_err( nf90_get_var(ncid,varid,prec,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"trecs",varid) , .true. )
      call disp_err( nf90_get_var(ncid,varid,trecs,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"trecw",varid) , .true. )
      call disp_err( nf90_get_var(ncid,varid,trecw,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"arecs",varid) , .true. )
      call disp_err( nf90_get_var(ncid,varid,arecs,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"arecw",varid) , .true. )
      call disp_err( nf90_get_var(ncid,varid,arecw,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"mvrec",varid) , .true. )
      call disp_err( nf90_get_var(ncid,varid,mvrec,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"nwrite",varid) , .true. )
      call disp_err( nf90_get_var(ncid,varid,nwrite,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"nwritet",varid) , .true. )
      call disp_err( nf90_get_var(ncid,varid,nwritet,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"nwritea",varid) , .true. )
      call disp_err( nf90_get_var(ncid,varid,nwritea,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"nwritemv",varid) , .true. )
      call disp_err( nf90_get_var(ncid,varid,nwritemv,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"nwriteh",varid) , .true. )
      call disp_err( nf90_get_var(ncid,varid,nwriteh,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"nrst",varid) , .true. )
      call disp_err( nf90_get_var(ncid,varid,nrst,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"nstatout",varid) , .true. )
      call disp_err( nf90_get_var(ncid,varid,nstatout,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"ndt",varid) , .true. )
      call disp_err( nf90_get_var(ncid,varid,ndt,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"icenter",varid) , .true. )
      call disp_err( nf90_get_var(ncid,varid,icenter,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"jcenter",varid) , .true. )
      call disp_err( nf90_get_var(ncid,varid,jcenter,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"old_format",varid) , .true. )
      call disp_err( nf90_get_var(ncid,varid,old_format,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"dt",varid) , .true. )
      call disp_err( nf90_get_var(ncid,varid,dt,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"dtlast",varid) , .true. )
      call disp_err( nf90_get_var(ncid,varid,dtlast,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"xcenter",varid) , .true. )
      call disp_err( nf90_get_var(ncid,varid,xcenter,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"ycenter",varid) , .true. )
      call disp_err( nf90_get_var(ncid,varid,ycenter,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"umove",varid) , .true. )
      call disp_err( nf90_get_var(ncid,varid,umove,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"vmove",varid) , .true. )
      call disp_err( nf90_get_var(ncid,varid,vmove,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"domainlocx",varid) , .true. )
      call disp_err( nf90_get_var(ncid,varid,domainlocx,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"domainlocy",varid) , .true. )
      call disp_err( nf90_get_var(ncid,varid,domainlocy,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"adaptmovetim",varid) , .true. )
      call disp_err( nf90_get_var(ncid,varid,adaptmovetim,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"cflmax",varid) , .true. )
      call disp_err( nf90_get_var(ncid,varid,cflmax,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"mtime",varid) , .true. )
      call disp_err( nf90_get_var(ncid,varid,mtime,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"stattim",varid) , .true. )
      call disp_err( nf90_get_var(ncid,varid,stattim,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"taptim",varid) , .true. )
      call disp_err( nf90_get_var(ncid,varid,taptim,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"rsttim",varid) , .true. )
      call disp_err( nf90_get_var(ncid,varid,rsttim,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"radtim",varid) , .true. )
      call disp_err( nf90_get_var(ncid,varid,radtim,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"prcltim",varid) , .true. )
      call disp_err( nf90_get_var(ncid,varid,prcltim,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"hailtrajtim",varid) , .true. )
      call disp_err( nf90_get_var(ncid,varid,hailtrajtim,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"adt",varid) , .true. )
      call disp_err( nf90_get_var(ncid,varid,adt,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"adtlast",varid) , .true. )
      call disp_err( nf90_get_var(ncid,varid,adtlast,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"acfl",varid) , .true. )
      call disp_err( nf90_get_var(ncid,varid,acfl,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"dbldt",varid) , .true. )
      call disp_err( nf90_get_var(ncid,varid,dbldt,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"mass1",varid) , .true. )
      call disp_err( nf90_get_var(ncid,varid,mass1,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"avgsfcu",varid) , .true. )
      call disp_err( nf90_get_var(ncid,varid,avgsfcu,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"avgsfcv",varid) , .true. )
      call disp_err( nf90_get_var(ncid,varid,avgsfcv,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"avgsfcs",varid) , .true. )
      call disp_err( nf90_get_var(ncid,varid,avgsfcs,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"avgsfcsu",varid) , .true. )
      call disp_err( nf90_get_var(ncid,varid,avgsfcsu,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"avgsfcsv",varid) , .true. )
      call disp_err( nf90_get_var(ncid,varid,avgsfcsv,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"avgsfct",varid) , .true. )
      call disp_err( nf90_get_var(ncid,varid,avgsfct,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"avgsfcq",varid) , .true. )
      call disp_err( nf90_get_var(ncid,varid,avgsfcq,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"avgsfcp",varid) , .true. )
      call disp_err( nf90_get_var(ncid,varid,avgsfcp,(/time_index/)) , .true. )

    ENDIF
#endif
  ENDIF

#ifdef MPI
      ! communicate to all other processors:
      call MPI_BCAST(nstep  ,1,MPI_INTEGER         ,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(srec   ,1,MPI_INTEGER         ,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(sirec  ,1,MPI_INTEGER         ,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(urec   ,1,MPI_INTEGER         ,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(vrec   ,1,MPI_INTEGER         ,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(wrec   ,1,MPI_INTEGER         ,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(nrec   ,1,MPI_INTEGER         ,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(mrec   ,1,MPI_INTEGER         ,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(prec   ,1,MPI_INTEGER         ,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(trecs  ,1,MPI_INTEGER         ,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(trecw  ,1,MPI_INTEGER         ,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(arecs  ,1,MPI_INTEGER         ,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(arecw  ,1,MPI_INTEGER         ,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(mvrec  ,1,MPI_INTEGER         ,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(nwrite ,1,MPI_INTEGER         ,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(nwritet,1,MPI_INTEGER         ,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(nwritea,1,MPI_INTEGER         ,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(nwritemv,1,MPI_INTEGER         ,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(nwriteh ,1,MPI_INTEGER         ,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(nrst   ,1,MPI_INTEGER         ,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(nstatout,1,MPI_INTEGER         ,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(ndt    ,1,MPI_INTEGER         ,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(icenter,1,MPI_INTEGER         ,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(jcenter,1,MPI_INTEGER         ,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(old_format,1,MPI_INTEGER         ,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(dt     ,1,MPI_REAL            ,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(dtlast ,1,MPI_REAL            ,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(xcenter,1,MPI_REAL            ,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(ycenter,1,MPI_REAL            ,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(umove  ,1,MPI_REAL            ,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(vmove  ,1,MPI_REAL            ,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(umovesfc,1,MPI_REAL            ,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(vmovesfc,1,MPI_REAL            ,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(umovem1,1,MPI_REAL            ,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(vmovem1,1,MPI_REAL            ,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(umovem2,1,MPI_REAL            ,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(vmovem2,1,MPI_REAL            ,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(domainlocx  ,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(domainlocy  ,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(adaptmovetim,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(cflmax ,1,MPI_REAL            ,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(mtime  ,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(stattim,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(taptim ,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(rsttim ,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(radtim ,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(prcltim,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(hailtrajtim,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(restart_time,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(spinup_t2p,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(adt    ,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(adtlast,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(acfl   ,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(dbldt  ,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(mass1  ,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(avgsfcu,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(avgsfcv,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(avgsfcs,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(avgsfcsu,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(avgsfcsv,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(avgsfct,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(avgsfcq,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(avgsfcp,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
#endif

!---------------------------------------------------------------
! budget variables:

    IF( myid.eq.0 )THEN

      IF( restart_read_format.eq.1 )THEN
        read(rstfile) nbudgetfoo
        read(rstfile) numqfoo
        if(myid.eq.0) print *
        if(myid.eq.0) print *,'  nbudget (old,new) = ',nbudgetfoo,nbudget
        if(myid.eq.0) print *,'  numq    (old,new) = ',numqfoo,numq
        if(myid.eq.0) print *
        if( nbudgetfoo.eq.nbudget )then
          read(rstfile) qbudget
        else
          allocate( dptmp(nbudgetfoo) )
          dptmp = 0.0
          read(rstfile) dptmp
          deallocate( dptmp )
          qbudget = 0.0
        endif
        if( numqfoo.eq.numq )then
          read(rstfile) asq
          read(rstfile) bsq
        else
          allocate( dptmp(numqfoo) )
          dptmp = 0.0
          read(rstfile) dptmp
          read(rstfile) dptmp
          deallocate( dptmp )
          asq = 0.0
          bsq = 0.0
        endif
#ifdef NETCDF
      ELSEIF( restart_read_format.eq.2 )THEN
        call disp_err( nf90_inq_varid(ncid,"qbudget",varid) , .true. )
        call disp_err( nf90_get_var(ncid,varid,qbudget,(/1,time_index/),(/nbudget,1/)) , .true. )
        call disp_err( nf90_inq_varid(ncid,"asq",varid) , .true. )
        call disp_err( nf90_get_var(ncid,varid,asq,(/1,time_index/),(/numq,1/)) , .true. )
        call disp_err( nf90_inq_varid(ncid,"bsq",varid) , .true. )
        call disp_err( nf90_get_var(ncid,varid,bsq,(/1,time_index/),(/numq,1/)) , .true. )
#endif
      ENDIF

    ELSE

      qbudget = 0.0
      asq = 0.0
      bsq = 0.0

    ENDIF

!---------------------------------------------------------------

  call readinfo(in_file,myid,oorecs,oorecu,oorecv,oorecw,'  ug,vg:')

  if( in_file.eq.in_file_yes )then

    if( myid.eq.0 )then
    if( interp_on_restart )then
      allocate( dumz(0:nzold+1) )
      dumz = 0.0
      read(rstfile) dumz
      do k=1,nk
        kk = iz(k)
        ug(k) = (1.0-rz(k))*dumz(kk) + rz(k)*dumz(kk+1)
      enddo
      read(rstfile) dumz
      do k=1,nk
        kk = iz(k)
        vg(k) = (1.0-rz(k))*dumz(kk) + rz(k)*dumz(kk+1)
      enddo
      deallocate( dumz )
    else
      read(rstfile) ug
      read(rstfile) vg
    endif
    endif

#ifdef MPI
    call MPI_BCAST(ug(kb),ke-kb+1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(vg(kb),ke-kb+1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
#endif
  endif

  if(myid.eq.0) print *,'  check1a:', orecs, orecu, orecv, orecw
  if(myid.eq.0) print *,'  check2a:',oorecs,oorecu,oorecv,oorecw

!-----------------------------------------------------------------------
!  open more files:

  rf2:  &
  IF( restart_read_format.eq.1 )THEN
    ! unformatted direct-access (grads) format:

    if( interp_on_restart )then
      nxr = nxold
      nyr = nyold
    else
      nxr = nx
      nyr = ny
    endif

    if( myid.eq.0 ) print *
    if( myid.eq.0 ) print *,'  nxr,nyr = ',nxr,nyr
    if( myid.eq.0 ) print *
    if( myid.eq.0 ) print *,'  restart_read_filetype = ',restart_read_filetype
    if( myid.eq.0 ) print *

  rft2:  &
  IF( restart_read_filetype.eq.1 )THEN

    !------------------
    ! one restart file (per stagger type):
    IF(myid.eq.nodeleader)THEN
      if(dowr) write(outfile,*)
      if(dowr) write(outfile,*) '  restart_read_filetype = 1 is no longer supported with restart_format = 1 '
      if(dowr) write(outfile,*)
      if(dowr) write(outfile,*) '      ... stopping cm1 ... '
      if(dowr) write(outfile,*)
      call stopcm1
    ENDIF

  ELSEIF( restart_read_filetype.eq.2 )THEN  rft2

    !------------------
    ! one restart file (per restart time):
    IF(myid.eq.0)THEN

      do i=1,maxstring
        fname(i:i) = ' '
      enddo

      fname = 'cm1rst_tXXXXXX_s.dat'
      write(fname( 9:14),101) rstnum
      if(dowr) write(outfile,*) '  fname=',fname
      open(unit=sfile,file=fname,form='unformatted',access='direct',recl=4*nxr*nyr,status='old',err=778)
      opens = .true.
      orecs = 1

      fname = 'cm1rst_tXXXXXX_u.dat'
      write(fname( 9:14),101) rstnum
      if(dowr) write(outfile,*) '  fname=',fname
      open(unit=ufile,file=fname,form='unformatted',access='direct',recl=4*(nxr+1)*nyr,status='old',err=778)
      openu = .true.
      orecu = 1

      fname = 'cm1rst_tXXXXXX_v.dat'
      write(fname( 9:14),101) rstnum
      if(dowr) write(outfile,*) '  fname=',fname
      open(unit=vfile,file=fname,form='unformatted',access='direct',recl=4*nxr*(nyr+1),status='old',err=778)
      openv = .true.
      orecv = 1

      fname = 'cm1rst_tXXXXXX_w.dat'
      write(fname( 9:14),101) rstnum
      if(dowr) write(outfile,*) '  fname=',fname
      open(unit=wfile,file=fname,form='unformatted',access='direct',recl=4*nxr*nyr,status='old',err=778)
      openw = .true.
      orecw = 1

      if(dowr) write(outfile,*)

    ENDIF

  ELSEIF( restart_read_filetype.eq.3 )THEN

    !------------------
    ! one restart file per processor:

    if(myid.eq.0) print *,'  one restart file per processor ... '

  if( newdecomp .or. interp_on_restart )then

      ! only processor 0 opens files:

    if( myid.eq.0 )then

      print *,'  processor 0 will read all files ... '
      print *,'      (this may take a while) '

      do i=1,maxstring
        fname(i:i) = ' '
      enddo

      fname = 'cm1rst_pXXXXXX_tYYYYYY_s.dat'
      write(fname(17:22),101) rstnum
      do n=0,numprocs_old-1
        write(fname( 9:14),101) n
        nfile = 1000 + 2*nfiles + n
        open(unit=nfile,file=trim(fname),form='unformatted',access='direct',recl=4*info(3,n)*info(4,n))
        open_sfiles = n
        print *,'  fname = ',nfile,trim(fname)
      enddo
      print *

      fname = 'cm1rst_pXXXXXX_tYYYYYY_u.dat'
      write(fname(17:22),101) rstnum
      do n=0,numprocs_old-1
        write(fname( 9:14),101) n
        nfile = 1000 + 3*nfiles + n
        open(unit=nfile,file=trim(fname),form='unformatted',access='direct',recl=4*(info(3,n)+1)*info(4,n))
        open_ufiles = n
        print *,'  fname = ',nfile,trim(fname)
      enddo
      print *

      fname = 'cm1rst_pXXXXXX_tYYYYYY_v.dat'
      write(fname(17:22),101) rstnum
      do n=0,numprocs_old-1
        write(fname( 9:14),101) n
        nfile = 1000 + 4*nfiles + n
        open(unit=nfile,file=trim(fname),form='unformatted',access='direct',recl=4*info(3,n)*(info(4,n)+1))
        open_vfiles = n
        print *,'  fname = ',nfile,trim(fname)
      enddo
      print *

      fname = 'cm1rst_pXXXXXX_tYYYYYY_w.dat'
      write(fname(17:22),101) rstnum
      do n=0,numprocs_old-1
        write(fname( 9:14),101) n
        nfile = 1000 + 5*nfiles + n
        open(unit=nfile,file=trim(fname),form='unformatted',access='direct',recl=4*info(3,n)*info(4,n))
        open_wfiles = n
        print *,'  fname = ',nfile,trim(fname)
      enddo
      print *

    endif

  else

      ! everybody opens their own files:

      if(myid.eq.0) print *,'  every processor will read data ... '

      do i=1,maxstring
        fname(i:i) = ' '
      enddo

      fname = 'cm1rst_pXXXXXX_tYYYYYY_s.dat'
      write(fname( 9:14),101) myid
      write(fname(17:22),101) rstnum
      open(unit=sfile,file=trim(fname),form='unformatted',access='direct',recl=4*ni*nj)
      opens = .true.
      if(myid.eq.0) print *,myid,sfile,trim(fname)

      fname = 'cm1rst_pXXXXXX_tYYYYYY_u.dat'
      write(fname( 9:14),101) myid
      write(fname(17:22),101) rstnum
      open(unit=ufile,file=trim(fname),form='unformatted',access='direct',recl=4*(ni+1)*nj)
      openu = .true.
      if(myid.eq.0) print *,myid,ufile,trim(fname)

      fname = 'cm1rst_pXXXXXX_tYYYYYY_v.dat'
      write(fname( 9:14),101) myid
      write(fname(17:22),101) rstnum
      open(unit=vfile,file=trim(fname),form='unformatted',access='direct',recl=4*ni*(nj+1))
      openv = .true.
      if(myid.eq.0) print *,myid,vfile,trim(fname)

      fname = 'cm1rst_pXXXXXX_tYYYYYY_w.dat'
      write(fname( 9:14),101) myid
      write(fname(17:22),101) rstnum
      open(unit=wfile,file=trim(fname),form='unformatted',access='direct',recl=4*ni*nj)
      openw = .true.
      if(myid.eq.0) print *,myid,wfile,trim(fname)

  endif

  ELSEIF( restart_read_filetype.eq.4 )THEN

    !------------------
    ! one restart file per node:

  if( newdecomp .or. interp_on_restart )then

      ! only processor 0 opens files:

    if( myid.eq.0 )then

      if(myid.eq.0) print *,'  processor 0 will read all data '
      if(myid.eq.0) print *,'  nioldmax,njoldmax    = ',nioldmax,njoldmax
      if(myid.eq.0) print *,'  nodes_old,ppnode_old = ',nodes_old,ppnode_old

      do i=1,maxstring
        fname(i:i) = ' '
      enddo

      fname = 'cm1rst_nXXXXXX_tYYYYYY_s.dat'
      write(fname(17:22),101) rstnum
      do n=0,nodes_old-1
        write(fname( 9:14),101) n
        nfile = 1000 + 2*nfiles + n
        open(unit=nfile,file=trim(fname),form='unformatted',access='direct',recl=4*(nioldmax*ppnode_old)*njoldmax)
        open_sfiles = n
        print *,'  fname = ',nfile,trim(fname)
      enddo
      print *

      fname = 'cm1rst_nXXXXXX_tYYYYYY_u.dat'
      write(fname(17:22),101) rstnum
      do n=0,nodes_old-1
        write(fname( 9:14),101) n
        nfile = 1000 + 3*nfiles + n
        open(unit=nfile,file=trim(fname),form='unformatted',access='direct',recl=4*((nioldmax+1)*ppnode_old)*njoldmax)
        open_ufiles = n
        print *,'  fname = ',nfile,trim(fname)
      enddo
      print *

      fname = 'cm1rst_nXXXXXX_tYYYYYY_v.dat'
      write(fname(17:22),101) rstnum
      do n=0,nodes_old-1
        write(fname( 9:14),101) n
        nfile = 1000 + 4*nfiles + n
        open(unit=nfile,file=trim(fname),form='unformatted',access='direct',recl=4*(nioldmax*ppnode_old)*(njoldmax+1))
        open_vfiles = n
        print *,'  fname = ',nfile,trim(fname)
      enddo
      print *

      fname = 'cm1rst_nXXXXXX_tYYYYYY_w.dat'
      write(fname(17:22),101) rstnum
      do n=0,nodes_old-1
        write(fname( 9:14),101) n
        nfile = 1000 + 5*nfiles + n
        open(unit=nfile,file=trim(fname),form='unformatted',access='direct',recl=4*(nioldmax*ppnode_old)*njoldmax)
        open_wfiles = n
        print *,'  fname = ',nfile,trim(fname)
      enddo
      print *

    endif

  else

      ! every nodeleader opens files:

    if( myid.eq.nodeleader )then

      if(myid.eq.0) print *,'  every nodeleader will read data ... '
      if(myid.eq.0) print *,'  nioldmax,njoldmax = ',nioldmax,njoldmax

      do i=1,maxstring
        fname(i:i) = ' '
      enddo

      fname = 'cm1rst_nXXXXXX_tYYYYYY_s.dat'
      write(fname( 9:14),101) mynode
      write(fname(17:22),101) rstnum
      open(unit=sfile,file=trim(fname),form='unformatted',access='direct',recl=4*(nioldmax*ppnode_old)*njoldmax)
      opens = .true.
      if(myid.eq.0) print *,myid,sfile,trim(fname)

      fname = 'cm1rst_nXXXXXX_tYYYYYY_u.dat'
      write(fname( 9:14),101) mynode
      write(fname(17:22),101) rstnum
      open(unit=ufile,file=trim(fname),form='unformatted',access='direct',recl=4*((nioldmax+1)*ppnode_old)*njoldmax)
      openu = .true.
      if(myid.eq.0) print *,myid,ufile,trim(fname)

      fname = 'cm1rst_nXXXXXX_tYYYYYY_v.dat'
      write(fname( 9:14),101) mynode
      write(fname(17:22),101) rstnum
      open(unit=vfile,file=trim(fname),form='unformatted',access='direct',recl=4*(nioldmax*ppnode_old)*(njoldmax+1))
      openv = .true.
      if(myid.eq.0) print *,myid,vfile,trim(fname)

      fname = 'cm1rst_nXXXXXX_tYYYYYY_w.dat'
      write(fname( 9:14),101) mynode
      write(fname(17:22),101) rstnum
      open(unit=wfile,file=trim(fname),form='unformatted',access='direct',recl=4*(nioldmax*ppnode_old)*njoldmax)
      openw = .true.
      if(myid.eq.0) print *,myid,wfile,trim(fname)
      if(myid.eq.0) print *

    endif

  endif

  ENDIF  rft2

  ENDIF  rf2


!---------------------------------------------------------------
! standard 2D:

      call  readr(ni,nj,1,1,nx,ny,tsk(ib,jb),'tsk     ',           &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                  ncid,time_index,restart_read_format,restart_read_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),datk3(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),info,reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)

      n = 1
      call  readr(ni,nj,1,1,nx,ny,rain(ib,jb,n),'rain    ',        &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                  ncid,time_index,restart_read_format,restart_read_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),datk3(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),info,reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call  readr(ni,nj,1,1,nx,ny,raincv(ib,jb,n),'raincv  ',        &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                  ncid,time_index,restart_read_format,restart_read_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),datk3(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),info,reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call  readr(ni,nj,1,1,nx,ny,hail(ib,jb,n),'hail    ',        &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                  ncid,time_index,restart_read_format,restart_read_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),datk3(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),info,reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call  readr(ni,nj,1,1,nx,ny,sws(ib,jb,n),'sws     ',         &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                  ncid,time_index,restart_read_format,restart_read_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),datk3(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),info,reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call  readr(ni,nj,1,1,nx,ny,svs(ib,jb,n),'svs     ',         &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                  ncid,time_index,restart_read_format,restart_read_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),datk3(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),info,reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call  readr(ni,nj,1,1,nx,ny,sps(ib,jb,n),'sps     ',         &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                  ncid,time_index,restart_read_format,restart_read_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),datk3(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),info,reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call  readr(ni,nj,1,1,nx,ny,srs(ib,jb,n),'srs     ',         &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                  ncid,time_index,restart_read_format,restart_read_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),datk3(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),info,reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call  readr(ni,nj,1,1,nx,ny,sgs(ib,jb,n),'sgs     ',         &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                  ncid,time_index,restart_read_format,restart_read_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),datk3(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),info,reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call  readr(ni,nj,1,1,nx,ny,sus(ib,jb,n),'sus     ',         &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                  ncid,time_index,restart_read_format,restart_read_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),datk3(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),info,reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call  readr(ni,nj,1,1,nx,ny,shs(ib,jb,n),'shs     ',         &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                  ncid,time_index,restart_read_format,restart_read_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),datk3(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),info,reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)

!---------------------------------------------------------------

      call readinfo(in_file,myid,oorecs,oorecu,oorecv,oorecw,'  nrain2:')
    if( in_file.eq.in_file_yes )then
    if( nrain.eq.2 )then
      n = 2
      call  readr(ni,nj,1,1,nx,ny,rain(ib,jb,n),'rain2   ',        &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                  ncid,time_index,restart_read_format,restart_read_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),datk3(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),info,reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call  readr(ni,nj,1,1,nx,ny,raincv(ib,jb,n),'raincv2 ',        &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                  ncid,time_index,restart_read_format,restart_read_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),datk3(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),info,reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call  readr(ni,nj,1,1,nx,ny,hail(ib,jb,n),'hail2   ',        &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                  ncid,time_index,restart_read_format,restart_read_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),datk3(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),info,reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call  readr(ni,nj,1,1,nx,ny,sws(ib,jb,n),'sws2    ',         &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                  ncid,time_index,restart_read_format,restart_read_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),datk3(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),info,reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call  readr(ni,nj,1,1,nx,ny,svs(ib,jb,n),'svs2    ',         &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                  ncid,time_index,restart_read_format,restart_read_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),datk3(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),info,reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call  readr(ni,nj,1,1,nx,ny,sps(ib,jb,n),'sps2    ',         &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                  ncid,time_index,restart_read_format,restart_read_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),datk3(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),info,reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call  readr(ni,nj,1,1,nx,ny,srs(ib,jb,n),'srs2    ',         &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                  ncid,time_index,restart_read_format,restart_read_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),datk3(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),info,reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call  readr(ni,nj,1,1,nx,ny,sgs(ib,jb,n),'sgs2    ',         &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                  ncid,time_index,restart_read_format,restart_read_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),datk3(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),info,reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call  readr(ni,nj,1,1,nx,ny,sus(ib,jb,n),'sus2    ',         &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                  ncid,time_index,restart_read_format,restart_read_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),datk3(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),info,reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call  readr(ni,nj,1,1,nx,ny,shs(ib,jb,n),'shs2    ',         &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                  ncid,time_index,restart_read_format,restart_read_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),datk3(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),info,reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
    else
      orecs = oorecs
      orecu = oorecu
      orecv = oorecv
      orecw = oorecw
    endif
    endif

  if(myid.eq.0) print *,'  check1b:', orecs, orecu, orecv, orecw
  if(myid.eq.0) print *,'  check2b:',oorecs,oorecu,oorecv,oorecw

!---------------------------------------------------------------
! standard 3D:

      call  readr(ni,nj,1,nk,nx,ny,rho(ib,jb,1),'rho     ',        &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                  ncid,time_index,restart_read_format,restart_read_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),datk3(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),info,reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call  readr(ni,nj,1,nk,nx,ny,prs(ib,jb,1),'prs     ',        &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                  ncid,time_index,restart_read_format,restart_read_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),datk3(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),info,reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call  readr(ni+1,nj,1,nk,nx+1,ny,ua(ib,jb,1),'ua      ',     &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecu,ufile,   &
                  ncid,time_index,restart_read_format,restart_read_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),datk3(0,0),nx,ny,nz,ixu,iy,iz,rxu,ry,rz,nxold+1,nyold,nzold,sigma(0),sigma0(0),info,reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2iu,d2ju,d3iu,d3ju)
      call  readr(ni,nj+1,1,nk,nx,ny+1,va(ib,jb,1),'va      ',     &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecv,vfile,   &
                  ncid,time_index,restart_read_format,restart_read_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),datk3(0,0),nx,ny,nz,ix,iyv,iz,rx,ryv,rz,nxold,nyold+1,nzold,sigma(0),sigma0(0),info,reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2iv,d2jv,d3iv,d3jv)
      call  readr(ni,nj,1,nk+1,nx,ny,wa(ib,jb,1),'wa      ',       &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecw,wfile,   &
                  ncid,time_index,restart_read_format,restart_read_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),datk3(0,0),nx,ny,nz+1,ix,iy,izw,rx,ry,rzw,nxold,nyold,nzold+1,sigma(0),sigma0(0),info,reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call  readr(ni,nj,1,nk,nx,ny,ppi(ib,jb,1),'ppi     ',        &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                  ncid,time_index,restart_read_format,restart_read_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),datk3(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),info,reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      if( interp_on_restart )then
        call  readr(ni,nj,1,nk,nx,ny,ppi(ib,jb,1),'pptot   ',        &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_read_format,restart_read_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),datk3(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),info,reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        do k=1,nk
        do j=1,nj
        do i=1,ni
          ppi(i,j,k) = ppi(i,j,k)-pi0(i,j,k)
        enddo
        enddo
        enddo
      else
        orecs = orecs + nkold
      endif
      call  readr(ni,nj,1,nk,nx,ny,tha(ib,jb,1),'tha     ',        &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                  ncid,time_index,restart_read_format,restart_read_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),datk3(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),info,reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      if( interp_on_restart )then
        call  readr(ni,nj,1,nk,nx,ny,tha(ib,jb,1),'thtot   ',        &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_read_format,restart_read_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),datk3(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),info,reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        do k=1,nk
        do j=1,nj
        do i=1,ni
          tha(i,j,k) = tha(i,j,k)-th0(i,j,k)
        enddo
        enddo
        enddo
      else
        orecs = orecs + nkold
      endif
      call  readr(ni,nj,1,nk,nx,ny,ppx(ib,jb,1),'ppx     ',        &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                  ncid,time_index,restart_read_format,restart_read_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),datk3(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),info,reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)

!---------------------------------------------------------------

      call readinfo(in_file,myid,oorecs,oorecu,oorecv,oorecw,'  phi:')
    if( in_file.eq.in_file_yes )then
    if( psolver.eq.6 )then
      call  readr(ni,nj,1,nk,nx,ny,phi1(ib,jb,1),'phi1    ',       &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                  ncid,time_index,restart_read_format,restart_read_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),datk3(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),info,reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call  readr(ni,nj,1,nk,nx,ny,phi2(ib,jb,1),'phi2    ',       &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                  ncid,time_index,restart_read_format,restart_read_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),datk3(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),info,reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
    else
      orecs = oorecs
      orecu = oorecu
      orecv = oorecv
      orecw = oorecw
    endif
    endif

  if(myid.eq.0) print *,'  check1c:', orecs, orecu, orecv, orecw
  if(myid.eq.0) print *,'  check2c:',oorecs,oorecu,oorecv,oorecw

!---------------------------------------------------------------
!  moisture

      call readinfo(in_file,myid,oorecs,oorecu,oorecv,oorecw,'  qa:',nold)
#ifdef MPI
      call MPI_BCAST(nold,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
#endif
      if( numq.ne.nold )then
        if(myid.eq.0)then
        print *
        print *,'  numq = ',numq
        print *,'  nold = ',nold
        print *
        print *,'  cannot change number of moisture variables on a restart '
        print *
        print *,'  ... stopping cm1 ... '
        print *
        endif
#ifdef MPI
        call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif
        call stopcm1
      endif
    if( in_file.eq.in_file_yes )then
    IF(imoist.eq.1)THEN
    do n=1,numq
      text1 = '        '
      write(text1(1:3),156) qname(n)
156   format(a3)
      call  readr(ni,nj,1,nk,nx,ny,qa(ib,jb,1,n),text1     ,       &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                  ncid,time_index,restart_read_format,restart_read_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),datk3(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),info,reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call pdef(ni,nj,nk,qa(ib,jb,kb,n))
    enddo
    ELSE
      orecs = oorecs
      orecu = oorecu
      orecv = oorecv
      orecw = oorecw
    ENDIF
    endif

  if(myid.eq.0) print *,'  check1d:', orecs, orecu, orecv, orecw
  if(myid.eq.0) print *,'  check2d:',oorecs,oorecu,oorecv,oorecw

!---------------------------------------------------------------
!  vars for eqtset=2

      call readinfo(in_file,myid,oorecs,oorecu,oorecv,oorecw,'  qxten:')
    if( in_file.eq.in_file_yes )then
    if(imoist.eq.1.and.eqtset.eq.2)then
      call  readr(ni,nj,1,nk,nx,ny,qpten(ib,jb,1),'qpten   ',      &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                  ncid,time_index,restart_read_format,restart_read_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),datk3(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),info,reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call  readr(ni,nj,1,nk,nx,ny,qtten(ib,jb,1),'qtten   ',      &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                  ncid,time_index,restart_read_format,restart_read_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),datk3(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),info,reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call  readr(ni,nj,1,nk,nx,ny,qvten(ib,jb,1),'qvten   ',      &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                  ncid,time_index,restart_read_format,restart_read_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),datk3(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),info,reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call  readr(ni,nj,1,nk,nx,ny,qcten(ib,jb,1),'qcten   ',      &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                  ncid,time_index,restart_read_format,restart_read_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),datk3(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),info,reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
    else
      orecs = oorecs
      orecu = oorecu
      orecv = oorecv
      orecw = oorecw
    endif
    endif

  if(myid.eq.0) print *,'  check1e:', orecs, orecu, orecv, orecw
  if(myid.eq.0) print *,'  check2e:',oorecs,oorecu,oorecv,oorecw

!---------------------------------------------------------------
!  subgrid tke for LES:

      call readinfo(in_file,myid,oorecs,oorecu,oorecv,oorecw,'  tkea:')
    if( in_file.eq.in_file_yes )then
    if( idoles .and. iusetke )then
      call  readr(ni,nj,1,nk+1,nx,ny,tkea(ib,jb,1),'tkea    ',       &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecw,wfile,   &
                  ncid,time_index,restart_read_format,restart_read_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),datk3(0,0),nx,ny,nz+1,ix,iy,izw,rx,ry,rzw,nxold,nyold,nzold+1,sigma(0),sigma0(0),info,reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call pdef(ni,nj,nk+1,tkea)
    else
      orecs = oorecs
      orecu = oorecu
      orecv = oorecv
      orecw = oorecw
    endif
    endif

  if(myid.eq.0) print *,'  check1f:', orecs, orecu, orecv, orecw
  if(myid.eq.0) print *,'  check2f:',oorecs,oorecu,oorecv,oorecw

!---------------------------------------------------------------
!  radiation:

      call readinfo(in_file,myid,oorecs,oorecu,oorecv,oorecw,'  radiation:')

    if( in_file.eq.in_file_yes )then

      radoptge1:  &
      if(radopt.ge.1)then
        call  readr(ni,nj,1,nk,nx,ny,lwten(ib,jb,1),'lwten   ',      &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_read_format,restart_read_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),datk3(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),info,reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,nk,nx,ny,swten(ib,jb,1),'swten   ',      &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_read_format,restart_read_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),datk3(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),info,reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,nk,nx,ny,thraten(ib,jb,1),'thraten ',    &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_read_format,restart_read_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),datk3(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),info,reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,dum1(ib,jb,1),'radsw   ',        &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_read_format,restart_read_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),datk3(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),info,reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
!$omp parallel do default(shared)  &
!$omp private(i,j)
        do j=1,nj
        do i=1,ni
          radsw(i,j) = dum1(i,j,1)
        enddo
        enddo
        call  readr(ni,nj,1,1,nx,ny,dum1(ib,jb,1),'rnflx   ',        &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_read_format,restart_read_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),datk3(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),info,reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
!$omp parallel do default(shared)  &
!$omp private(i,j)
        do j=1,nj
        do i=1,ni
          rnflx(i,j) = dum1(i,j,1)
        enddo
        enddo
        call  readr(ni,nj,1,1,nx,ny,dum1(ib,jb,1),'radswnet',        &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_read_format,restart_read_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),datk3(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),info,reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
!$omp parallel do default(shared)  &
!$omp private(i,j)
        do j=1,nj
        do i=1,ni
          radswnet(i,j) = dum1(i,j,1)
        enddo
        enddo
        call  readr(ni,nj,1,1,nx,ny,dum1(ib,jb,1),'radlwin ',        &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_read_format,restart_read_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),datk3(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),info,reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
!$omp parallel do default(shared)  &
!$omp private(i,j)
        do j=1,nj
        do i=1,ni
          radlwin(i,j) = dum1(i,j,1)
        enddo
        enddo
        do n=1,nrad2d
        if( n.lt.10 )then
          text1 = 'radX    '
          write(text1(4:4),181) n
181       format(i1.1)
        elseif( n.lt.100 )then
          text1 = 'radXX   '
          write(text1(4:5),182) n
182       format(i2.2)
        elseif( n.lt.1000 )then
          text1 = 'radXXX  '
          write(text1(4:6),183) n
183       format(i3.3)
        else
          stop 11611
        endif
        call  readr(ni,nj,1,1,nx,ny,dum1(ib,jb,1),text1,             &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_read_format,restart_read_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),datk3(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),info,reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
!$omp parallel do default(shared)  &
!$omp private(i,j)
        do j=1,nj
        do i=1,ni
          rad2d(i,j,n) = dum1(i,j,1)
        enddo
        enddo
        enddo

      else

        orecs = oorecs
        orecu = oorecu
        orecv = oorecv
        orecw = oorecw

      endif  radoptge1
    endif

  if(myid.eq.0) print *,'  check1g:', orecs, orecu, orecv, orecw
  if(myid.eq.0) print *,'  check2g:',oorecs,oorecu,oorecv,oorecw

!---------------------------------------------------------------

      call readinfo(in_file,myid,oorecs,oorecu,oorecv,oorecw,'  effx:')

    if( in_file.eq.in_file_yes )then

      if( radopt.ge.1 .and. ptype.eq.5 )then
        call  readr(ni,nj,1,nk,nx,ny,effc(ib,jb,1),'effc    ',       &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_read_format,restart_read_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),datk3(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),info,reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,nk,nx,ny,effi(ib,jb,1),'effi    ',       &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_read_format,restart_read_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),datk3(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),info,reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,nk,nx,ny,effs(ib,jb,1),'effs    ',       &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_read_format,restart_read_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),datk3(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),info,reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,nk,nx,ny,effr(ib,jb,1),'effr    ',       &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_read_format,restart_read_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),datk3(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),info,reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,nk,nx,ny,effg(ib,jb,1),'effg    ',       &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_read_format,restart_read_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),datk3(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),info,reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,nk,nx,ny,effis(ib,jb,1),'effis   ',      &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_read_format,restart_read_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),datk3(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),info,reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      else
        orecs = oorecs
        orecu = oorecu
        orecv = oorecv
        orecw = oorecw
      endif

    endif

  if(myid.eq.0) print *,'  check1h:', orecs, orecu, orecv, orecw
  if(myid.eq.0) print *,'  check2h:',oorecs,oorecu,oorecv,oorecw

!---------------------------------------------------------------
!  cu param:

      call readinfo(in_file,myid,oorecs,oorecu,oorecv,oorecw,'  cuparam:')

    if( in_file.eq.in_file_yes )then

      if(cuparam.ge.1)then
        call  readr(ni,nj,1,nk,nx,ny,thcuten(ib,jb,1),'thcuten ',      &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_read_format,restart_read_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),datk3(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),info,reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,nk,nx,ny,qvcuten(ib,jb,1),'qvcuten ',      &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_read_format,restart_read_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),datk3(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),info,reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,nk,nx,ny,qccuten(ib,jb,1),'qccuten ',      &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_read_format,restart_read_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),datk3(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),info,reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,nk,nx,ny,qicuten(ib,jb,1),'qicuten ',      &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_read_format,restart_read_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),datk3(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),info,reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,nk,nx,ny,ucuten(ib,jb,1),'ucuten  ',      &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_read_format,restart_read_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),datk3(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),info,reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,nk,nx,ny,vcuten(ib,jb,1),'vcuten  ',      &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_read_format,restart_read_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),datk3(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),info,reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      else
        orecs = oorecs
        orecu = oorecu
        orecv = oorecv
        orecw = oorecw
      endif

    endif

  if(myid.eq.0) print *,'  check1i:', orecs, orecu, orecv, orecw
  if(myid.eq.0) print *,'  check2i:',oorecs,oorecu,oorecv,oorecw

!---------------------------------------------------------------
!  surface:
!     I don't know how many of these are really needed in restart
!     files, but let's include them all for now ... just to be safe

      call readinfo(in_file,myid,oorecs,oorecu,oorecv,oorecw,'  sfcvars1:')

    if( in_file.eq.in_file_yes )then

      sfcvars1:  &
      if( ((oceanmodel.eq.2).or.(ipbl.ge.1).or.(sfcmodel.ge.1)) .and. sfcmodel.ge.1 )then

        call  readr(ni,nj,1,1,nx,ny,ust(ib,jb),'ust     ',           &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_read_format,restart_read_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),datk3(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),info,reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,znt(ib,jb),'znt     ',           &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_read_format,restart_read_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),datk3(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),info,reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,cd(ib,jb),'cd      ',            &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_read_format,restart_read_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),datk3(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),info,reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,ch(ib,jb),'ch      ',            &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_read_format,restart_read_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),datk3(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),info,reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,cq(ib,jb),'cq      ',            &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_read_format,restart_read_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),datk3(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),info,reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,u1(ib,jb),'u1      ',            &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_read_format,restart_read_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),datk3(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),info,reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,v1(ib,jb),'v1      ',            &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_read_format,restart_read_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),datk3(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),info,reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,s1(ib,jb),'s1      ',            &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_read_format,restart_read_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),datk3(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),info,reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,th1(ib,jb),'th1     ',           &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_read_format,restart_read_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),datk3(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),info,reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,qv1(ib,jb),'qv1     ',           &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_read_format,restart_read_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),datk3(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),info,reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,u10(ib,jb),'u10     ',           &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_read_format,restart_read_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),datk3(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),info,reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,v10(ib,jb),'v10     ',           &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_read_format,restart_read_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),datk3(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),info,reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,s10(ib,jb),'s10     ',           &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_read_format,restart_read_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),datk3(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),info,reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,xland(ib,jb),'xland   ',         &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_read_format,restart_read_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),datk3(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),info,reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,thflux(ib,jb),'thflux  ',        &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_read_format,restart_read_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),datk3(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),info,reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,qvflux(ib,jb),'qvflux  ',        &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_read_format,restart_read_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),datk3(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),info,reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,psfc(ib,jb),'psfc    ',          &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_read_format,restart_read_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),datk3(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),info,reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)

      else

        orecs = oorecs
        orecu = oorecu
        orecv = oorecv
        orecw = oorecw

      endif  sfcvars1

    endif

  if(myid.eq.0) print *,'  check1j:', orecs, orecu, orecv, orecw
  if(myid.eq.0) print *,'  check2j:',oorecs,oorecu,oorecv,oorecw



      call readinfo(in_file,myid,oorecs,oorecu,oorecv,oorecw,'  sfcvars2:')

    if( in_file.eq.in_file_yes )then

      sfcvars2:  &
      if( ((oceanmodel.eq.2).or.(ipbl.ge.1).or.(sfcmodel.ge.1)) .and. tbc.eq.3 )then

        call  readr(ni,nj,1,1,nx,ny,ustt(ib,jb),'ustt    ',          &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_read_format,restart_read_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),datk3(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),info,reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,  ut(ib,jb),'ut      ',          &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_read_format,restart_read_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),datk3(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),info,reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,  vt(ib,jb),'vt      ',          &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_read_format,restart_read_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),datk3(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),info,reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,  st(ib,jb),'st      ',          &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_read_format,restart_read_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),datk3(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),info,reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)

      else

        orecs = oorecs
        orecu = oorecu
        orecv = oorecv
        orecw = oorecw

      endif  sfcvars2

    endif


  if(myid.eq.0) print *,'  check1k:', orecs, orecu, orecv, orecw
  if(myid.eq.0) print *,'  check2k:',oorecs,oorecu,oorecv,oorecw


      call readinfo(in_file,myid,oorecs,oorecu,oorecv,oorecw,'  sfcvars3:')

    if( in_file.eq.in_file_yes )then

      sfcvars3:  &
      if( ((oceanmodel.eq.2).or.(ipbl.ge.1).or.(sfcmodel.ge.1)) .and. sfcmodel.ge.1 )then

        call  readr(ni,nj,1,1,nx,ny,dum1(ib,jb,1),'lu_index',        &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_read_format,restart_read_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),datk3(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),info,reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
!$omp parallel do default(shared)  &
!$omp private(i,j)
        do j=1,nj
        do i=1,ni
          lu_index(i,j) = nint(dum1(i,j,1))
        enddo
        enddo
        call  readr(ni,nj,1,1,nx,ny,dum1(ib,jb,1),'kpbl2d  ',        &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_read_format,restart_read_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),datk3(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),info,reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
!$omp parallel do default(shared)  &
!$omp private(i,j)
        do j=1,nj
        do i=1,ni
          kpbl2d(i,j) = nint(dum1(i,j,1))
        enddo
        enddo
        call  readr(ni,nj,1,1,nx,ny,hfx(ib,jb),'hfx     ',           &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_read_format,restart_read_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),datk3(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),info,reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,qfx(ib,jb),'qfx     ',           &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_read_format,restart_read_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),datk3(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),info,reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,hpbl(ib,jb),'hpbl    ',          &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_read_format,restart_read_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),datk3(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),info,reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,wspd(ib,jb),'wspd    ',          &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_read_format,restart_read_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),datk3(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),info,reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,psim(ib,jb),'psim    ',          &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_read_format,restart_read_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),datk3(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),info,reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,psih(ib,jb),'psih    ',          &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_read_format,restart_read_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),datk3(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),info,reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,gz1oz0(ib,jb),'gz1oz0  ',        &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_read_format,restart_read_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),datk3(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),info,reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,br(ib,jb),'br      ',            &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_read_format,restart_read_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),datk3(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),info,reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,CHS(ib,jb),'chs     ',           &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_read_format,restart_read_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),datk3(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),info,reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,CHS2(ib,jb),'chs2    ',          &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_read_format,restart_read_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),datk3(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),info,reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,CQS2(ib,jb),'cqs2    ',          &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_read_format,restart_read_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),datk3(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),info,reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,CPMM(ib,jb),'cpmm    ',          &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_read_format,restart_read_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),datk3(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),info,reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,ZOL(ib,jb),'zol     ',           &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_read_format,restart_read_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),datk3(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),info,reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,MAVAIL(ib,jb),'mavail  ',        &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_read_format,restart_read_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),datk3(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),info,reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,MOL(ib,jb),'mol     ',           &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_read_format,restart_read_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),datk3(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),info,reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,RMOL(ib,jb),'rmol    ',          &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_read_format,restart_read_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),datk3(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),info,reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,REGIME(ib,jb),'regime  ',        &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_read_format,restart_read_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),datk3(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),info,reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,LH(ib,jb),'lh      ',            &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_read_format,restart_read_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),datk3(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),info,reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,tmn(ib,jb),'tmn     ',           &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_read_format,restart_read_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),datk3(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),info,reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,FLHC(ib,jb),'flhc    ',          &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_read_format,restart_read_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),datk3(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),info,reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,FLQC(ib,jb),'flqc    ',          &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_read_format,restart_read_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),datk3(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),info,reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,QGH(ib,jb),'qgh     ',           &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_read_format,restart_read_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),datk3(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),info,reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,CK(ib,jb),'ck      ',            &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_read_format,restart_read_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),datk3(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),info,reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,CKA(ib,jb),'cka     ',           &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_read_format,restart_read_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),datk3(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),info,reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,CDA(ib,jb),'cda     ',           &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_read_format,restart_read_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),datk3(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),info,reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,USTM(ib,jb),'ustm    ',          &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_read_format,restart_read_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),datk3(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),info,reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,QSFC(ib,jb),'qsfc    ',          &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_read_format,restart_read_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),datk3(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),info,reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,T2(ib,jb),'t2      ',            &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_read_format,restart_read_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),datk3(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),info,reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,Q2(ib,jb),'q2      ',            &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_read_format,restart_read_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),datk3(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),info,reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,TH2(ib,jb),'th2     ',           &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_read_format,restart_read_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),datk3(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),info,reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,EMISS(ib,jb),'emiss   ',         &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_read_format,restart_read_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),datk3(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),info,reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,THC(ib,jb),'thc     ',           &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_read_format,restart_read_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),datk3(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),info,reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,ALBD(ib,jb),'albd    ',          &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_read_format,restart_read_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),datk3(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),info,reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,gsw(ib,jb),'gsw     ',           &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_read_format,restart_read_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),datk3(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),info,reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,glw(ib,jb),'glw     ',           &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_read_format,restart_read_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),datk3(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),info,reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,chklowq(ib,jb),'chklowq ',       &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_read_format,restart_read_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),datk3(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),info,reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,capg(ib,jb),'capg    ',          &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_read_format,restart_read_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),datk3(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),info,reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,snowc(ib,jb),'snowc   ',         &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_read_format,restart_read_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),datk3(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),info,reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,fm(ib,jb),'fm      ',            &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_read_format,restart_read_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),datk3(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),info,reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,fh(ib,jb),'fh      ',            &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_read_format,restart_read_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),datk3(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),info,reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,mznt(ib,jb),'mznt    ',          &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_read_format,restart_read_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),datk3(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),info,reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,swspd(ib,jb),'swspd   ',         &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_read_format,restart_read_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),datk3(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),info,reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,wstar(ib,jb),'wstar   ',         &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_read_format,restart_read_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),datk3(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),info,reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,delta(ib,jb),'delta   ',         &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_read_format,restart_read_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),datk3(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),info,reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        do n=1,num_soil_layers
          if( n.lt.10 )then
            text1 = 'tslbX   '
            write(text1(5:5),171) n
171         format(i1.1)
          elseif( n.lt.100 )then
            text1 = 'tslbXX  '
            write(text1(5:6),172) n
172         format(i2.2)
          else
            stop 22122
          endif
        call  readr(ni,nj,1,1,nx,ny,tslb(ib,jb,n),text1,             &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_read_format,restart_read_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),datk3(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),info,reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        enddo

      else

        orecs = oorecs
        orecu = oorecu
        orecv = oorecv
        orecw = oorecw

      endif  sfcvars3

    endif

  if(myid.eq.0) print *,'  check1l:', orecs, orecu, orecv, orecw
  if(myid.eq.0) print *,'  check2l:',oorecs,oorecu,oorecv,oorecw


      call readinfo(in_file,myid,oorecs,oorecu,oorecv,oorecw,'  omlvars:')

    if( in_file.eq.in_file_yes )then

      omlvars:  &
      if(oceanmodel.eq.2)then
        call  readr(ni,nj,1,1,nx,ny,tml(ib,jb),'tml     ',           &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_read_format,restart_read_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),datk3(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),info,reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,t0ml(ib,jb),'t0ml    ',          &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_read_format,restart_read_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),datk3(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),info,reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,hml(ib,jb),'hml     ',           &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_read_format,restart_read_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),datk3(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),info,reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,h0ml(ib,jb),'h0ml    ',          &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_read_format,restart_read_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),datk3(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),info,reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,huml(ib,jb),'huml    ',          &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_read_format,restart_read_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),datk3(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),info,reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,hvml(ib,jb),'hvml    ',          &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_read_format,restart_read_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),datk3(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),info,reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,tmoml(ib,jb),'tmoml   ',         &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_read_format,restart_read_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),datk3(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),info,reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)

      else

        orecs = oorecs
        orecu = oorecu
        orecv = oorecv
        orecw = oorecw

      endif  omlvars

    endif

  if(myid.eq.0) print *,'  check1m:', orecs, orecu, orecv, orecw
  if(myid.eq.0) print *,'  check2m:',oorecs,oorecu,oorecv,oorecw

!---------------------------------------------------------------
!  MYNN and MYJ:

      call readinfo(in_file,myid,oorecs,oorecu,oorecv,oorecw,'  myx:',n)

      if(myid.eq.0) print *,'  n,npbl3d = ',n,npbl3d

    if( in_file.eq.in_file_yes )then

    IF( ( ipbl.eq.4 .or. ipbl.eq.5 .or. ipbl.eq.6 ) .and. ( init_pbl.eq.0 ) )THEN
      if(myid.eq.0) print *,'  reading myx vars ... '
      do n=1,npbl3d
        call  readr(ni,nj,1,nk,nx,ny,pbl3d(ib,jb,1,n),'myx-3d  ',        &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_read_format,restart_read_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),datk3(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),info,reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      enddo

    ELSE

        orecs = oorecs
        orecu = oorecu
        orecv = oorecv
        orecw = oorecw

    ENDIF

    endif

    !----------

      call readinfo(in_file,myid,oorecs,oorecu,oorecv,oorecw,'  myx-2d:',n)

      if(myid.eq.0) print *,'  n,npbl2d = ',n,npbl2d

    if( in_file.eq.in_file_yes )then

    IF( ( sfcmodel.eq.6 .or. sfcmodel.eq.7 ) .and. ( init_sfcmodel.eq.0 ) )THEN
    do n=1,npbl2d
      call  readr(ni,nj,1,1,nx,ny,pbl2d(ib,jb,n),'myx-2d  ',        &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                  ncid,time_index,restart_read_format,restart_read_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),datk3(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),info,reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
    enddo

    ELSE

        orecs = oorecs
        orecu = oorecu
        orecv = oorecv
        orecw = oorecw

    ENDIF

    endif

    !----------

      call readinfo(in_file,myid,oorecs,oorecu,oorecv,oorecw,'  myx-2d-int:',n)

      if(myid.eq.0) print *,'  n,npbl2dint = ',n,npbl2dint

    if( in_file.eq.in_file_yes )then
    IF( ( sfcmodel.eq.6 .or. sfcmodel.eq.7 ) .and. ( init_sfcmodel.eq.0 ) )THEN
    do n=1,npbl2dint
      call  readr(ni,nj,1,1,nx,ny,dum1(ib,jb,1),'myx-2di ',        &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                  ncid,time_index,restart_read_format,restart_read_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),datk3(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),info,reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      do j=1,nj
      do i=1,ni
        pbl2dint(i,j,n) = nint(dum1(i,j,1))
      enddo
      enddo
    enddo
      if(myid.eq.0) print *,'  ... done reading myx vars '

    ELSE

        orecs = oorecs
        orecu = oorecu
        orecv = oorecv
        orecw = oorecw

    ENDIF
    endif

  if(myid.eq.0) print *,'  check1n:', orecs, orecu, orecv, orecw
  if(myid.eq.0) print *,'  check2n:',oorecs,oorecu,oorecv,oorecw

!---------------------------------------------------------------
!  PBL tendencies:

      call readinfo(in_file,myid,oorecs,oorecu,oorecv,oorecw,'  pbl tends:')

    if( in_file.eq.in_file_yes )then

    pbltends:  &
    IF( ( ipbl.eq.4 .or. ipbl.eq.5 .or. ipbl.eq.6 ) .and. ( init_pbl.eq.0 ) )THEN
      if(myid.eq.0) print *,'  reading pbl tends ... '
      call  readr(ni,nj,1,nk,nx,ny,thpten(ib,jb,1),'thpten  ',     &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                  ncid,time_index,restart_read_format,restart_read_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),datk3(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),info,reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call  readr(ni,nj,1,nk,nx,ny,qvpten(ib,jb,1),'qvpten  ',     &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                  ncid,time_index,restart_read_format,restart_read_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),datk3(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),info,reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call  readr(ni,nj,1,nk,nx,ny,qcpten(ib,jb,1),'qcpten  ',     &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                  ncid,time_index,restart_read_format,restart_read_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),datk3(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),info,reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call  readr(ni,nj,1,nk,nx,ny,qipten(ib,jb,1),'qipten  ',     &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                  ncid,time_index,restart_read_format,restart_read_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),datk3(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),info,reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call  readr(ni,nj,1,nk,nx,ny,upten(ib,jb,1),'upten   ',      &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                  ncid,time_index,restart_read_format,restart_read_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),datk3(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),info,reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call  readr(ni,nj,1,nk,nx,ny,vpten(ib,jb,1),'vpten   ',      &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                  ncid,time_index,restart_read_format,restart_read_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),datk3(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),info,reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call  readr(ni,nj,1,nk,nx,ny,qnipten(ib,jb,1),'qnipten ',    &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                  ncid,time_index,restart_read_format,restart_read_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),datk3(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),info,reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call  readr(ni,nj,1,nk,nx,ny,qncpten(ib,jb,1),'qncpten ',    &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                  ncid,time_index,restart_read_format,restart_read_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),datk3(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),info,reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)

    ELSE

        orecs = oorecs
        orecu = oorecu
        orecv = oorecv
        orecw = oorecw

    ENDIF  pbltends

    endif

!---------------------------------------------------------------
!  2-part turbulence param:

      call readinfo(in_file,myid,oorecs,oorecu,oorecv,oorecw,'  t2p_avg1:')

    if( in_file.eq.in_file_yes )then

      if( myid.eq.0 .and. t2p_avg.eq.1 )then
        if( interp_on_restart )then
          allocate( dumz(0:nzold+1) )
          dumz = 0.0
          !---
          read(rstfile) dumz
          do k=1,nk
            kk = iz(k)
            gamk(k) = (1.0-rz(k))*dumz(kk) + rz(k)*dumz(kk+1)
          enddo
          !---
          read(rstfile) dumz
          do k=1,nk
            kk = iz(k)
            kmw(k) = (1.0-rz(k))*dumz(kk) + rz(k)*dumz(kk+1)
          enddo
          !---
          read(rstfile) dumz
          do k=1,nk
            kk = iz(k)
            ufw(k) = (1.0-rz(k))*dumz(kk) + rz(k)*dumz(kk+1)
          enddo
          !---
          read(rstfile) dumz
          do k=1,nk
            kk = iz(k)
            vfw(k) = (1.0-rz(k))*dumz(kk) + rz(k)*dumz(kk+1)
          enddo
          !---
          read(rstfile) dumz
          do k=1,nk
            kk = iz(k)
            u1b(k) = (1.0-rz(k))*dumz(kk) + rz(k)*dumz(kk+1)
          enddo
          !---
          read(rstfile) dumz
          do k=1,nk
            kk = iz(k)
            v1b(k) = (1.0-rz(k))*dumz(kk) + rz(k)*dumz(kk+1)
          enddo
          !---
          deallocate( dumz )
        else
         read(rstfile) gamk
         read(rstfile) kmw
         read(rstfile) ufw
         read(rstfile) vfw
         read(rstfile) u1b
         read(rstfile) v1b
        endif
      endif
#ifdef MPI
      call MPI_BCAST(gamk,ke-kb+1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(kmw ,ke-kb+1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(ufw ,ke-kb+1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(vfw ,ke-kb+1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(u1b ,ke-kb+1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(v1b ,ke-kb+1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
#endif
    endif

  if(myid.eq.0) print *,'  check1o:', orecs, orecu, orecv, orecw
  if(myid.eq.0) print *,'  check2o:',oorecs,oorecu,oorecv,oorecw

!---------------------------------------------------------------
!  passive tracers:

      call readinfo(in_file,myid,oorecs,oorecu,oorecv,oorecw,'  pt:',nvar)

    if( in_file.eq.in_file_yes )then

      if( myid.eq.0 )then
#ifdef NETCDF
        if( restart_read_format.eq.2 )then
          call disp_err( nf90_inq_varid(ncid,"npt",varid) , .true. )
          call disp_err( nf90_get_var(ncid,varid,nvar,(/time_index/)) , .true. )
          if( iptra.eq.0 ) nvar = 0
        endif
#endif
        print *,'  nvar_npt = ',nvar
      endif

      readpt:  &
      if( iptra.eq.1 .and. nvar.gt.0 )then
            do n=1,min(nvar,npt)
              if( n.lt.10 )then
                text1 = 'ptX     '
                write(text1(3:3),161) n
161             format(i1.1)
              elseif( n.lt.100 )then
                text1 = 'ptXX    '
                write(text1(3:4),162) n
162             format(i2.2)
              else
                stop 11512
              endif
              call  readr(ni,nj,1,nk,nx,ny,pta(ib,jb,1,n),text1,           &
                          ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                          ncid,time_index,restart_read_format,restart_read_filetype,myi1p,myi2p,myj1p,myj2p, &
                          dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),datk3(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),info,reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
            enddo
            if( nvar.ne.npt )then
              orecs = oorecs
              orecu = oorecu
              orecv = oorecv
              orecw = oorecw
            endif
      else
        orecs = oorecs
        orecu = oorecu
        orecv = oorecv
        orecw = oorecw
      endif  readpt

    endif

  if(myid.eq.0) print *,'  check1p:', orecs, orecu, orecv, orecw
  if(myid.eq.0) print *,'  check2p:',oorecs,oorecu,oorecv,oorecw

!---------------------------------------------------------------
!  time-average vars:

      call readinfo(in_file,myid,oorecs,oorecu,oorecv,oorecw,'  timavg:')

    if( in_file.eq.in_file_yes )then

        if( myid.eq.0 )then
           read(rstfile) ntim_old
           read(rstfile) ntavr_old
           read(rstfile) nsfctavr_old
           read(rstfile) keta_old
           read(rstfile) utav_old
           read(rstfile) vtav_old
           read(rstfile) wtav_old
           read(rstfile) ttav_old
           read(rstfile) qtav_old
           read(rstfile) etav_old
           read(rstfile) uutav_old
           read(rstfile) vvtav_old
           read(rstfile) ntavg
           read(rstfile) rtavg
          print *
          print *,'  ntim,ntim_old         = ',ntim,ntim_old
          print *,'  ntavr,ntavr_old       = ',ntavr,ntavr_old
          print *,'  nsfctavr,nsfctavr_old = ',nsfctavr,nsfctavr_old
          print *,'  keta,keta_old         = ',keta,keta_old
          print *,'  utav,utav_old         = ',utav,utav_old
          print *,'  vtav,vtav_old         = ',vtav,vtav_old
          print *,'  wtav,wtav_old         = ',wtav,wtav_old
          print *,'  ttav,ttav_old         = ',ttav,ttav_old
          print *,'  qtav,qtav_old         = ',qtav,qtav_old
          print *,'  etav,etav_old         = ',etav,etav_old
          print *,'  uutav,uutav_old       = ',uutav,uutav_old
          print *,'  vvtav,vvtav_old       = ',vvtav,vvtav_old
          print *
        endif

      readtimavg:  &
      IF( dotimeavg )THEN

#ifdef MPI
        call MPI_BCAST(ntavr_old      ,1,MPI_INTEGER         ,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(nsfctavr_old   ,1,MPI_INTEGER         ,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(keta_old       ,1,MPI_INTEGER         ,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(utav_old       ,1,MPI_INTEGER         ,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(vtav_old       ,1,MPI_INTEGER         ,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(wtav_old       ,1,MPI_INTEGER         ,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(ttav_old       ,1,MPI_INTEGER         ,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(qtav_old       ,1,MPI_INTEGER         ,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(etav_old       ,1,MPI_INTEGER         ,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(uutav_old      ,1,MPI_INTEGER         ,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(vvtav_old      ,1,MPI_INTEGER         ,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(ntavg       ,ntim,MPI_INTEGER         ,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(rtavg       ,ntim,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
#endif

        if(myid.eq.0 )then
          if( ntim.ne.ntim_old )then
            print *
            print *,'  with time-avg code, cannot change ntim on a restart (for now) '
            print *
            print *,'  ... stopping cm1 ... '
            print *
            call stopcm1
          endif
        endif
!!!        if(myid.eq.0 )then
!!!          if( keta.ne.keta_old )then
!!!            print *
!!!            print *,'  with time-avg code, cannot change number of levels on a restart (for now) '
!!!            print *
!!!            print *,'  ... stopping cm1 ... '
!!!            print *
!!!            call stopcm1
!!!          endif
!!!        endif
        if( (utav_old.ne.utav) .or.  &
            (vtav_old.ne.vtav) .or.  &
            (wtav_old.ne.wtav) .or.  &
            (ttav_old.ne.ttav) .or.  &
            (qtav_old.ne.qtav) .or.  &
            (etav_old.ne.etav) .or.  &
            (uutav_old.ne.uutav) .or.  &
            (vvtav_old.ne.vvtav) )then
          restart_time = mtime
          spinup_t2p = .true.
          if(myid.eq.0)then
            print *
            print *,'  using t2p-timeavg code, but some variables are not in restart file '
            print *,'    so, changing settings to allow for spinup of t2p vars: '
            print *
            print *,'      restart_time = ',restart_time
            print *,'      spinup_t2p   = ',spinup_t2p
            print *
          endif
        endif


        !----------
        !  3d timeavg vars:
        nloopa:  &
        do n=1,ntavr_old
          if(myid.eq.0) print *,'  n = ',n
          if(     n.eq.utav_old .or. n.eq.uutav_old )then
            call  readr(ni+1,nj,kbta,keta,nx+1,ny,dumu(ib,jb,1),'timavgu ',     &
                        ni,nj,ngxy,myid,numprocs,nodex,nodey,orecu,ufile,   &
                        ncid,time_index,restart_read_format,restart_read_filetype,myi1p,myi2p,myj1p,myj2p, &
                        dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),datk3(0,0),nx,ny,keta,ixu,iy,iz(1),rxu,ry,rz(1),nxold+1,nyold,keta_old,sigma(0),sigma0(0),info,reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2iu,d2ju,d3iu,d3ju)
            nn = 0
            if( n.eq.utav_old  .and. utav.ge.1  ) nn = utav
            if( n.eq.uutav_old .and. uutav.ge.1 ) nn = uutav
            if( nn.ge.1 )then
              do k=kbta,keta
              do j=1,nj
              do i=1,ni+1
                timavg(i,j,k,nn) = dumu(i,j,k)
              enddo
              enddo
              enddo
            endif
          elseif( n.eq.vtav_old .or. n.eq.vvtav_old )then
            call  readr(ni,nj+1,kbta,keta,nx,ny+1,dumv(ib,jb,1),'timavgv ',     &
                        ni,nj,ngxy,myid,numprocs,nodex,nodey,orecv,vfile,   &
                        ncid,time_index,restart_read_format,restart_read_filetype,myi1p,myi2p,myj1p,myj2p, &
                        dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),datk3(0,0),nx,ny,keta,ix,iyv,iz(1),rx,ryv,rz(1),nxold,nyold+1,keta_old,sigma(0),sigma0(0),info,reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2iv,d2jv,d3iv,d3jv)
            nn = 0
            if( n.eq.vtav_old  .and. vtav.ge.1  ) nn = vtav
            if( n.eq.vvtav_old .and. vvtav.ge.1 ) nn = vvtav
            if( nn.ge.1 )then
              do k=kbta,keta
              do j=1,nj+1
              do i=1,ni
                timavg(i,j,k,nn) = dumv(i,j,k)
              enddo
              enddo
              enddo
            endif
          else
            call  readr(ni,nj,kbta,keta,nx,ny,dum1(ib,jb,kbta),'timavg  ',               &
                        ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,                           &
                        ncid,time_index,restart_read_format,restart_read_filetype,myi1p,myi2p,myj1p,myj2p, &
                        dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),datk3(0,0),nx,ny,keta,ix,iy,iz(1),rx,ry,rz(1),nxold,nyold,keta_old,sigma(0),sigma0(0),info,reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
            nn = 0
            if( n.eq.wtav_old  .and. wtav.ge.1  ) nn = wtav
            if( n.eq.ttav_old  .and. ttav.ge.1  ) nn = ttav
            if( n.eq.qtav_old  .and. qtav.ge.1  ) nn = qtav
            if( n.eq.etav_old  .and. etav.ge.1  ) nn = etav
            if( nn.ge.1 )then
              do k=kbta,keta
              do j=1,nj
              do i=1,ni
                timavg(i,j,k,nn) = dum1(i,j,k)
              enddo
              enddo
              enddo
            endif
          endif
        enddo  nloopa
        !----------
        !  sfc timeavg vars:
        nloopb:  &
        do n=1,nsfctavr_old
          if(myid.eq.0) print *,'  n = ',n
          call  readr(ni,nj,1,1,nx,ny,dum1(ib,jb,1),'sfctimav',                       &
                      ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,                           &
                      ncid,time_index,restart_read_format,restart_read_filetype,myi1p,myi2p,myj1p,myj2p, &
                      dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),datk3(0,0),nx,ny,keta,ix,iy,iz(1),rx,ry,rz(1),nxold,nyold,keta_old,sigma(0),sigma0(0),info,reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
          if( n.le.nsfctavr )then
            do j=1,nj
            do i=1,ni
              sfctimavg(i,j,n) = dum1(i,j,1)
            enddo
            enddo
          endif
        enddo  nloopb
        !----------
        !  3d tavg arrays:
        nloopc:  &
        do n=1,ntavr_old
          if(myid.eq.0) print *,'  n = ',n
        do nt=1,ntim
          if(     n.eq.utav_old .or. n.eq.uutav_old )then
            call  readr(ni+1,nj,kbta,keta,nx+1,ny,dumu(ib,jb,1),'tavgu   ',     &
                        ni,nj,ngxy,myid,numprocs,nodex,nodey,orecu,ufile,   &
                        ncid,time_index,restart_read_format,restart_read_filetype,myi1p,myi2p,myj1p,myj2p, &
                        dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),datk3(0,0),nx,ny,keta,ixu,iy,iz(1),rxu,ry,rz(1),nxold+1,nyold,keta_old,sigma(0),sigma0(0),info,reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2iu,d2ju,d3iu,d3ju)
            nn = 0
            if( n.eq.utav_old  .and. utav.ge.1  ) nn = utav
            if( n.eq.uutav_old .and. uutav.ge.1 ) nn = uutav
            if( nn.ge.1 )then
              do k=kbta,keta
              do j=1,nj
              do i=1,ni+1
                tavg(i,j,k,nt,nn) = dumu(i,j,k)
              enddo
              enddo
              enddo
            endif
          elseif( n.eq.vtav_old .or. n.eq.vvtav_old )then
            call  readr(ni,nj+1,kbta,keta,nx,ny+1,dumv(ib,jb,1),'tavgv   ',     &
                        ni,nj,ngxy,myid,numprocs,nodex,nodey,orecv,vfile,   &
                        ncid,time_index,restart_read_format,restart_read_filetype,myi1p,myi2p,myj1p,myj2p, &
                        dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),datk3(0,0),nx,ny,keta,ix,iyv,iz(1),rx,ryv,rz(1),nxold,nyold+1,keta_old,sigma(0),sigma0(0),info,reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2iv,d2jv,d3iv,d3jv)
            nn = 0
            if( n.eq.vtav_old  .and. vtav.ge.1  ) nn = vtav
            if( n.eq.vvtav_old .and. vvtav.ge.1 ) nn = vvtav
            if( nn.ge.1 )then
              do k=kbta,keta
              do j=1,nj+1
              do i=1,ni
                tavg(i,j,k,nt,nn) = dumv(i,j,k)
              enddo
              enddo
              enddo
            endif
          else
            call  readr(ni,nj,kbta,keta,nx,ny,dum1(ib,jb,kbta),'tavg    ',                       &
                        ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,                           &
                        ncid,time_index,restart_read_format,restart_read_filetype,myi1p,myi2p,myj1p,myj2p, &
                        dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),datk3(0,0),nx,ny,keta,ix,iy,iz(1),rx,ry,rz(1),nxold,nyold,keta_old,sigma(0),sigma0(0),info,reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
            nn = 0
            if( n.eq.wtav_old  .and. wtav.ge.1  ) nn = wtav
            if( n.eq.ttav_old  .and. ttav.ge.1  ) nn = ttav
            if( n.eq.qtav_old  .and. qtav.ge.1  ) nn = qtav
            if( n.eq.etav_old  .and. etav.ge.1  ) nn = etav
            if( nn.ge.1 )then
              do k=kbta,keta
              do j=1,nj
              do i=1,ni
                tavg(i,j,k,nt,nn) = dum1(i,j,k)
              enddo
              enddo
              enddo
            endif
          endif
        enddo
        enddo  nloopc
        !----------
        !  sfc tavg arrays:
        nloopd:  &
        do n=1,nsfctavr_old
          if(myid.eq.0) print *,'  n = ',n
        do nt=1,ntim
          call  readr(ni,nj,1,1,nx,ny,dum1(ib,jb,1),'sfctavg ',                                &
                      ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,                           &
                      ncid,time_index,restart_read_format,restart_read_filetype,myi1p,myi2p,myj1p,myj2p, &
                      dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),datk3(0,0),nx,ny,keta,ix,iy,iz(1),rx,ry,rz(1),nxold,nyold,keta_old,sigma(0),sigma0(0),info,reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
          if( n.le.nsfctavr )then
            do j=1,nj
            do i=1,ni
              sfctavg(i,j,nt,n) = dum1(i,j,1)
            enddo
            enddo
          endif
        enddo
        enddo  nloopd
        !----------

      ELSE

        orecs = oorecs
        orecu = oorecu
        orecv = oorecv
        orecw = oorecw

      ENDIF  readtimavg

    endif

  if(myid.eq.0) print *,'  check1q:', orecs, orecu, orecv, orecw
  if(myid.eq.0) print *,'  check2q:',oorecs,oorecu,oorecv,oorecw

!---------------------------------------------------------------
!  parcels:

      call readinfo(in_file,myid,oorecs,oorecu,oorecv,oorecw,'  parcels:')

    if( in_file.eq.in_file_yes )then

      if( myid.eq.0 )then
        if( restart_read_format.eq.1 )then
          read(rstfile) nvar
#ifdef NETCDF
        elseif( restart_read_format.eq.2 )then
          call disp_err( nf90_inq_varid(ncid,"numparcels",varid) , .true. )
          call disp_err( nf90_get_var(ncid,varid,nvar,(/time_index/)) , .true. )
          if( iprcl.eq.0 ) nvar = 0
#endif
        endif
        print *,'  nvar_parcels = ',nvar
      endif

#ifdef MPI
      call MPI_BCAST(nvar,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
#endif

      readparcels:  &
      if( nvar.gt.0 )then
          ! only read position info:
          if( myid.eq.0 )then
            IF( nvar.eq.nparcels )THEN
              ! easy:  restart file matches current config
              if( restart_read_format.eq.1 )then
                read(rstfile) ploc
#ifdef NETCDF
              elseif( restart_read_format.eq.2 )then
                call disp_err( nf90_inq_varid(ncid,"ploc",varid) , .true. )
                n = 3
                call disp_err( nf90_get_var(ncid,varid,ploc,(/1,1,time_index/),(/nparcels,n,1/)) , .true. )
#endif
              endif
            ELSE
              ! annoying:  restart file has different nparcels than current config
              IF( iprcl.eq.1 )THEN
                if( .not. terrain_flag )then
                  do np=1,nparcels
                    ploc(np,1) = pdata(np,prx)
                    ploc(np,2) = pdata(np,pry)
                    ploc(np,3) = pdata(np,prz)
                  enddo
                else
                  do np=1,nparcels
                    ploc(np,1) = pdata(np,prx)
                    ploc(np,2) = pdata(np,pry)
                    ploc(np,3) = pdata(np,prsig)
                  enddo
                endif
              ENDIF
              if( myid.eq.0 ) print *,'  start pfoo ' 
              allocate( pfoo(nvar,3) )
              if( restart_read_format.eq.1 )then
                read(rstfile) pfoo
#ifdef NETCDF
              elseif( restart_read_format.eq.2 )then
                call disp_err( nf90_inq_varid(ncid,"ploc",varid) , .true. )
                n = 3
                call disp_err( nf90_get_var(ncid,varid,pfoo,(/1,1,time_index/),(/nvar,n,1/)) , .true. )
#endif
              endif
              IF( iprcl.eq.1 )THEN
                do n=1,3
                do np=1,min(nvar,nparcels)
                  ploc(np,n) = pfoo(np,n)
                enddo
                enddo
              ENDIF
              deallocate( pfoo )
              if( myid.eq.0 ) print *,'  end pfoo ' 
            ENDIF
          endif
          IF( iprcl.eq.1 )THEN
#ifdef MPI
            call MPI_BCAST(ploc,3*nparcels,MPI_REAL,0,MPI_COMM_WORLD,ierr)
#endif
            if( .not. terrain_flag )then
              DO np=1,nparcels
                pdata(np,prx)=ploc(np,1)
                pdata(np,pry)=ploc(np,2)
                pdata(np,prz)=ploc(np,3)
              ENDDO
            else
              DO np=1,nparcels
                pdata(np,prx)=ploc(np,1)
                pdata(np,pry)=ploc(np,2)
                pdata(np,prsig)=ploc(np,3)
              ENDDO
            endif
            restart_prcl = .true.
          ENDIF
      else
        if( myid.eq.0 ) print *
        if( myid.eq.0 ) print *,'  Note:  no parcel data in the restart file '
        if( myid.eq.0 ) print *
        restart_prcl = .false.
      endif readparcels

    endif

!---------------------------------------------------------------
!  hailstones:

      call readinfo(in_file,myid,oorecs,oorecu,oorecv,oorecw,'  hailstones:')

    if( in_file.eq.in_file_yes )then

      if( myid.eq.0 )then
        if( restart_read_format.eq.1 )then
          read(rstfile) nvar
#ifdef NETCDF
        elseif( restart_read_format.eq.2 )then
          call disp_err( nf90_inq_varid(ncid,"numhailtrajs",varid) , .true. )
          call disp_err( nf90_get_var(ncid,varid,nvar,(/time_index/)) , .true. )
          if( ihailtraj.eq.0 ) nvar = 0
#endif
        endif
        print *,'  nvar_hailtraj = ',nvar
      endif

#ifdef MPI
      call MPI_BCAST(nvar,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
#endif

      readhail:  &
      if( nvar.gt.0 )then
          ! unlike parcels, need position info and internal hailstone characteristics:
          !  x,y,z(or sigma),u,v,w,time,d,dense,vt,ts,fw,itype.  13 in total
          if( myid.eq.0 )then
            IF( nvar.eq.nhailtrajs )THEN
              ! easy:  restart file matches current config
              if( restart_read_format.eq.1 )then
                read(rstfile) hailloc
#ifdef NETCDF
              elseif( restart_read_format.eq.2 )then
                call disp_err( nf90_inq_varid(ncid,"hailloc",varid) , .true. )
                n = 13
                call disp_err( nf90_get_var(ncid,varid,hailloc,(/1,1,time_index/),(/nhailtrajs,n,1/)) , .true. )
#endif
              endif
              print *,'  restart max d: ', MAXVAL(hailloc(:,8))
            ELSE
              ! annoying:  restart file has different nhailtrajs than current config
              ! (hailstones beyond those in the restart file keep their init3d values)
              IF( ihailtraj.eq.1 )THEN
                do np=1,nhailtrajs
                  hailloc(np,1) = haildata(np,hailx)
                  hailloc(np,2) = haildata(np,haily)
                  if( .not. terrain_flag )then
                    hailloc(np,3) = haildata(np,hailz)
                  else
                    hailloc(np,3) = haildata(np,hailsig)
                  endif
                  hailloc(np,4) = haildata(np,hailu)
                  hailloc(np,5) = haildata(np,hailv)
                  hailloc(np,6) = haildata(np,hailw)
                  hailloc(np,7) = haildata(np,hailtime)
                  hailloc(np,8) = haildata(np,haild)
                  hailloc(np,9) = haildata(np,haildense)
                  hailloc(np,10) = haildata(np,hailvt)
                  hailloc(np,11) = haildata(np,hailts)
                  hailloc(np,12) = haildata(np,hailfw)
                  hailloc(np,13) = haildata(np,hailitype)
                enddo
              ENDIF
              if( myid.eq.0 ) print *,'  start hfoo '
              allocate( hfoo(nvar,13) )
              if( restart_read_format.eq.1 )then
                read(rstfile) hfoo
#ifdef NETCDF
              elseif( restart_read_format.eq.2 )then
                call disp_err( nf90_inq_varid(ncid,"hailloc",varid) , .true. )
                n = 13
                call disp_err( nf90_get_var(ncid,varid,hfoo,(/1,1,time_index/),(/nvar,n,1/)) , .true. )
#endif
              endif
              IF( ihailtraj.eq.1 )THEN
                do n=1,13
                do np=1,min(nvar,nhailtrajs)
                  hailloc(np,n) = hfoo(np,n)
                enddo
                enddo
              ENDIF
              deallocate( hfoo )
              if( myid.eq.0 ) print *,'  end hfoo '
            ENDIF
          endif
          IF( ihailtraj.eq.1 )THEN
#ifdef MPI
            call MPI_BCAST(hailloc,13*nhailtrajs,MPI_REAL,0,MPI_COMM_WORLD,ierr)
#endif
            DO np=1,nhailtrajs
              haildata(np,hailx)=hailloc(np,1)
              haildata(np,haily)=hailloc(np,2)
              if( .not. terrain_flag )then
                haildata(np,hailz)=hailloc(np,3)
              else
                haildata(np,hailsig)=hailloc(np,3)
              endif
              haildata(np,hailu)=hailloc(np,4)
              haildata(np,hailv)=hailloc(np,5)
              haildata(np,hailw)=hailloc(np,6)
              haildata(np,hailtime)=hailloc(np,7)
              haildata(np,haild)=hailloc(np,8)
              haildata(np,haildense)=hailloc(np,9)
              haildata(np,hailvt)=hailloc(np,10)
              haildata(np,hailts)=hailloc(np,11)
              haildata(np,hailfw)=hailloc(np,12)
              haildata(np,hailitype)=hailloc(np,13)
            ENDDO
            restart_hail = .true.
          ENDIF
      else
        if( myid.eq.0 ) print *
        if( myid.eq.0 ) print *,'  Note:  no hailstone data in the restart file '
        if( myid.eq.0 ) print *
        restart_hail = .false.
      endif readhail

    else
      if( myid.eq.0 ) print *
      if( myid.eq.0 ) print *,'  Note:  no hailstone data in the restart file '
      if( myid.eq.0 ) print *
      restart_hail = .false.
    endif

  if(myid.eq.0) print *,'  check1r:', orecs, orecu, orecv, orecw
  if(myid.eq.0) print *,'  check2r:',oorecs,oorecu,oorecv,oorecw

!---------------------------------------------------------------
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!---------------------------------------------------------------
!  open bc:

      call readinfo(in_file,myid,oorecs,oorecu,oorecv,oorecw,'  openbc:')

    if( in_file.eq.in_file_yes )then

        !----------------------
        !cccccccccccccccccccccc
        !----------------------
!        if(myid.eq.0)then
          ndum = nyold
!        else
!          ndum = 1
!        endif
        allocate( dumy(0:ndum+1,0:nzold+1) )
        !----------------------
        if(myid.eq.0)then
          print *,'  radbcw, ndum,nzold = ',ndum,nzold
          aname = 'radbcw'
          call readbc(ndum=ndum,                      &
                      nzold=nzold,                    &
                      dum=dumy,                       &
                      restart_read_format=restart_read_format,  &
                      ncid=ncid,                      &
                      aname=aname,                    &
                      interp_on_restart=interp_on_restart)
        endif
        jfoo = 0
#ifdef MPI
        jfoo = myj1p(myid)-1
        call MPI_BCAST(dumy,(ndum+2)*(nzold+2),MPI_REAL,0,MPI_COMM_WORLD,ierr)
#endif
        if( wbc.eq.2 .and. ibw.eq.1 )then
          if( interp_on_restart )then
            do k=1,nk
            do j=1,nj
              w1=(1.0-ry(j))*(1.0-rz(k))
              w3=ry(j)*(1.0-rz(k))
              w4=(1.0-ry(j))*rz(k)
              w6=ry(j)*rz(k)
              radbcw(j,k) = ( dumy(iy(j)  ,iz(k)  )*w1    &
                             +dumy(iy(j)+1,iz(k)  )*w3    &
                             +dumy(iy(j)  ,iz(k)+1)*w4    &
                             +dumy(iy(j)+1,iz(k)+1)*w6 )
            enddo
            enddo
          else
            do k=1,nk
            do j=1,nj
              radbcw(j,k) = dumy(jfoo+j,k)
            enddo
            enddo
          endif
        endif
        !----------------------
        !cccccccccccccccccccccc
        !----------------------
        if(myid.eq.0)then
          print *,'  radbce, ndum = ',ndum
          aname = 'radbce'
          call readbc(ndum=ndum,                      &
                      nzold=nzold,                    &
                      dum=dumy,                       &
                      restart_read_format=restart_read_format,  &
                      ncid=ncid,                      &
                      aname=aname,                    &
                      interp_on_restart=interp_on_restart)
        endif
        jfoo = 0
#ifdef MPI
        jfoo = myj1p(myid)-1
        call MPI_BCAST(dumy,(ndum+2)*(nzold+2),MPI_REAL,0,MPI_COMM_WORLD,ierr)
#endif
        if( ebc.eq.2 .and. ibe.eq.1 )then
          if( interp_on_restart )then
            do k=1,nk
            do j=1,nj
              w1=(1.0-ry(j))*(1.0-rz(k))
              w3=ry(j)*(1.0-rz(k))
              w4=(1.0-ry(j))*rz(k)
              w6=ry(j)*rz(k)
              radbce(j,k) = ( dumy(iy(j)  ,iz(k)  )*w1    &
                             +dumy(iy(j)+1,iz(k)  )*w3    &
                             +dumy(iy(j)  ,iz(k)+1)*w4    &
                             +dumy(iy(j)+1,iz(k)+1)*w6 )
            enddo
            enddo
          else
            do k=1,nk
            do j=1,nj
              radbce(j,k) = dumy(jfoo+j,k)
            enddo
            enddo
          endif
        endif
        deallocate( dumy )
        !----------------------
        !cccccccccccccccccccccc
        !----------------------
!        if(myid.eq.0)then
          ndum = nxold
!        else
!          ndum = 1
!        endif
        allocate( dumx(0:ndum+1,0:nzold+1) )
        !----------------------
        if(myid.eq.0)then
          print *,'  radbcs, ndum = ',ndum
          aname = 'radbcs'
          call readbc(ndum=ndum,                      &
                      nzold=nzold,                    &
                      dum=dumx,                       &
                      restart_read_format=restart_read_format,  &
                      ncid=ncid,                      &
                      aname=aname,                    &
                      interp_on_restart=interp_on_restart)
        endif
        ifoo = 0
#ifdef MPI
        ifoo = myi1p(myid)-1
        call MPI_BCAST(dumx,(ndum+2)*(nzold+2),MPI_REAL,0,MPI_COMM_WORLD,ierr)
#endif
        if( sbc.eq.2 .and. ibs.eq.1 )then
          if( interp_on_restart )then
            do k=1,nk
            do i=1,ni
              w1=(1.0-rx(i))*(1.0-rz(k))
              w2=rx(i)*(1.0-rz(k))
              w4=(1.0-rx(i))*rz(k)
              w5=rx(i)*rz(k)
              radbcs(i,k) = ( dumx(ix(i)  ,iz(k)  )*w1    &
                             +dumx(ix(i)+1,iz(k)  )*w2    &
                             +dumx(ix(i)  ,iz(k)+1)*w4    &
                             +dumx(ix(i)+1,iz(k)+1)*w5 )
            enddo
            enddo
          else
            do k=1,nk
            do i=1,ni
              radbcs(i,k) = dumx(ifoo+i,k)
            enddo
            enddo
          endif
        endif
        !----------------------
        !cccccccccccccccccccccc
        !----------------------
        if(myid.eq.0)then
          print *,'  radbcn, ndum = ',ndum
          aname = 'radbcn'
          call readbc(ndum=ndum,                      &
                      nzold=nzold,                    &
                      dum=dumx,                       &
                      restart_read_format=restart_read_format,  &
                      ncid=ncid,                      &
                      aname=aname,                    &
                      interp_on_restart=interp_on_restart)
        endif
        ifoo = 0
#ifdef MPI
        ifoo = myi1p(myid)-1
        call MPI_BCAST(dumx,(ndum+2)*(nzold+2),MPI_REAL,0,MPI_COMM_WORLD,ierr)
#endif
        if( nbc.eq.2 .and. ibn.eq.1 )then
          if( interp_on_restart )then
            do k=1,nk
            do i=1,ni
              w1=(1.0-rx(i))*(1.0-rz(k))
              w2=rx(i)*(1.0-rz(k))
              w4=(1.0-rx(i))*rz(k)
              w5=rx(i)*rz(k)
              radbcn(i,k) = ( dumx(ix(i)  ,iz(k)  )*w1    &
                             +dumx(ix(i)+1,iz(k)  )*w2    &
                             +dumx(ix(i)  ,iz(k)+1)*w4    &
                             +dumx(ix(i)+1,iz(k)+1)*w5 )
            enddo
            enddo
          else
            do k=1,nk
            do i=1,ni
              radbcn(i,k) = dumx(ifoo+i,k)
            enddo
            enddo
          endif
        endif
        deallocate( dumx )
        !----------------------
        !cccccccccccccccccccccc
        !----------------------

    endif

  if(myid.eq.0) print *,'  check1s:', orecs, orecu, orecv, orecw
  if(myid.eq.0) print *,'  check2s:',oorecs,oorecu,oorecv,oorecw

!---------------------------------------------------------------
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!---------------------------------------------------------------
!  151001:  use theta  (over-rides perturbation value that was read-in above)


      call readinfo(in_file,myid,oorecs,oorecu,oorecv,oorecw,'  theta:')

    if( in_file.eq.in_file_yes )then

    IF( restart_use_theta )THEN
      call  readr(ni,nj,1,nk,nx,ny,dum1(ib,jb,1),'theta   ',       &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                  ncid,time_index,restart_read_format,restart_read_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),datk3(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),info,reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      do k=1,nk
      do j=1,nj
      do i=1,ni
        tha(i,j,k) = dum1(i,j,k)-th0(i,j,k)
      enddo
      enddo
      enddo
    ELSE
      orecs = oorecs
      orecu = oorecu
      orecv = oorecv
      orecw = oorecw
    ENDIF

    endif

  if(myid.eq.0) print *,'  check1t:', orecs, orecu, orecv, orecw
  if(myid.eq.0) print *,'  check2t:',oorecs,oorecu,oorecv,oorecw

!---------------------------------------------------------------
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!---------------------------------------------------------------


    call readinfo(in_file,myid,oorecs,oorecu,oorecv,oorecw,'  u0:')
    if( in_file.eq.in_file_yes )then
    IF( do_lsnudge .or. do_adapt_move )THEN
      call  readr(ni+1,nj,1,nk,nx+1,ny,u0(ib,jb,1),'u0      ',     &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecu,ufile,   &
                  ncid,time_index,restart_read_format,restart_read_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),datk3(0,0),nx,ny,nz,ixu,iy,iz,rxu,ry,rz,nxold+1,nyold,nzold,sigma(0),sigma0(0),info,reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2iu,d2ju,d3iu,d3ju)

    ELSE
      orecs = oorecs
      orecu = oorecu
      orecv = oorecv
      orecw = oorecw
    ENDIF
    endif

  if(myid.eq.0) print *,'  check1u:', orecs, orecu, orecv, orecw
  if(myid.eq.0) print *,'  check2u:',oorecs,oorecu,oorecv,oorecw

    call readinfo(in_file,myid,oorecs,oorecu,oorecv,oorecw,'  v0:')
    if( in_file.eq.in_file_yes )then
    IF( do_lsnudge .or. do_adapt_move )THEN
      call  readr(ni,nj+1,1,nk,nx,ny+1,v0(ib,jb,1),'v0      ',     &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecv,vfile,   &
                  ncid,time_index,restart_read_format,restart_read_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),datk3(0,0),nx,ny,nz,ix,iyv,iz,rx,ryv,rz,nxold,nyold+1,nzold,sigma(0),sigma0(0),info,reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2iv,d2jv,d3iv,d3jv)
    ELSE
      orecs = oorecs
      orecu = oorecu
      orecv = oorecv
      orecw = oorecw
    ENDIF
    endif

  if(myid.eq.0) print *,'  check1v:', orecs, orecu, orecv, orecw
  if(myid.eq.0) print *,'  check2v:',oorecs,oorecu,oorecv,oorecw

!---------------------------------------------------------------
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!---------------------------------------------------------------

    if(myid.eq.0) print *
    if(myid.eq.0) print *,'  Done reading restart file '
    if(myid.eq.0) print *

    IF( restart_read_format.eq.1 )THEN
      if(openrst)  close(unit=rstfile)
      if(openinfo) close(unit=infofile)
      if(opens)    close(unit=sfile)
      if(openu)    close(unit=ufile)
      if(openv)    close(unit=vfile)
      if(openw)    close(unit=wfile)
      if(open_nifiles.ge.1)then
        do n=1,open_nifiles
          nfile = 1000 + 0*nfiles + n
          close(unit=nfile)
        enddo
      endif
      if(open_nxfiles.ge.1)then
        do n=1,open_nxfiles
          nfile = 1000 + 1*nfiles + n
          close(unit=nfile)
        enddo
      endif
      if(open_sfiles.ge.1)then
        do n=0,open_sfiles
          nfile = 1000 + 2*nfiles + n
          close(unit=nfile)
        enddo
      endif
      if(open_ufiles.ge.1)then
        do n=0,open_ufiles
          nfile = 1000 + 3*nfiles + n
          close(unit=nfile)
        enddo
      endif
      if(open_vfiles.ge.1)then
        do n=0,open_vfiles
          nfile = 1000 + 4*nfiles + n
          close(unit=nfile)
        enddo
      endif
      if(open_wfiles.ge.1)then
        do n=0,open_wfiles
          nfile = 1000 + 5*nfiles + n
          close(unit=nfile)
        enddo
      endif
#ifdef NETCDF
    ELSEIF( restart_read_format.eq.2 )THEN
      if( myid.eq.0 )then
        call disp_err( nf90_close(ncid) , .true. )
      endif
#endif
    ENDIF

    if( restarted ) nrst = nrst+1

!---------

      if(dowr) write(outfile,*)
      if(dowr) write(outfile,*) '  From restart file: '
      if(dowr) write(outfile,*)
      if(dowr) write(outfile,*) '   mtime   = ',mtime
      if(dowr) write(outfile,*) '   stattim = ',stattim
      if(dowr) write(outfile,*) '   taptim  = ',taptim
      if(dowr) write(outfile,*) '   rsttim  = ',rsttim
      if(dowr) write(outfile,*) '   radtim  = ',radtim
      if(dowr) write(outfile,*) '   prcltim = ',prcltim
      if(dowr) write(outfile,*)
      if(dowr) write(outfile,*) '   nstep   = ',nstep
      if(dowr) write(outfile,*) '   srec    = ',srec
      if(dowr) write(outfile,*) '   sirec   = ',sirec
      if(dowr) write(outfile,*) '   urec    = ',urec
      if(dowr) write(outfile,*) '   vrec    = ',vrec
      if(dowr) write(outfile,*) '   wrec    = ',wrec
      if(dowr) write(outfile,*) '   nrec    = ',nrec
      if(dowr) write(outfile,*) '   mrec    = ',mrec
      if(dowr) write(outfile,*) '   prec    = ',prec
      if(dowr) write(outfile,*) '   nwrite  = ',nwrite
      if(dowr) write(outfile,*) '   nrst    = ',nrst
      if(dowr) write(outfile,*) '   nstatout= ',nstatout
      if( dodomaindiag )then
      if(dowr) write(outfile,*) '   trecs   = ',trecs
      if(dowr) write(outfile,*) '   trecw   = ',trecw
      if(dowr) write(outfile,*) '   nwritet = ',nwritet
      endif
      if( doazimavg )then
      if(dowr) write(outfile,*) '   arecs   = ',arecs
      if(dowr) write(outfile,*) '   arecw   = ',arecw
      if(dowr) write(outfile,*) '   nwritea = ',nwritea
      endif
      if( dohifrq )then
      if(dowr) write(outfile,*) '   nwriteh = ',nwriteh
      endif
      if( do_adapt_move )then
      if(dowr) write(outfile,*) '   nwritemv      = ',nwritemv
      if(dowr) write(outfile,*) '   mvrec         = ',mvrec
      if(dowr) write(outfile,*) '   umove         = ',umove
      if(dowr) write(outfile,*) '   vmove         = ',vmove
      if(dowr) write(outfile,*) '   adaptmovetime = ',adaptmovetim
      endif
      if(dowr) write(outfile,*)
      if(dowr) write(outfile,*) '   mass1 = ',mass1
      if(dowr) write(outfile,*)
      if(dowr) write(outfile,*) '   dt      = ',dt
      if(dowr) write(outfile,*) '   dtlast  = ',dtlast
      if(dowr) write(outfile,*) '   dbldt   = ',dbldt
      if(dowr) write(outfile,*) '   adt     = ',adt
      if(dowr) write(outfile,*) '   adtlast = ',adtlast

!---------------------------------------------------------------
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!---------------------------------------------------------------
!  For interp_on_restart, use log-profile winds below lowest model level
!    251231: setting this off by default; need a way to adjust the pressure field

    doit = .false.

    IF( doit .and. bbc.eq.3 .and. interp_on_restart .and. (.not.terrain_flag) )THEN

      nf=0

      if(myid.eq.0) print *
      if(myid.eq.0) print *,'  log-profile winds below previous lowest model level '
      if(myid.eq.0) print *

      do k=1,nk
        if( zh(1,1,k).lt.sigma0(1) )then
          if(myid.eq.0) print *,k,zh(1,1,k),sigma0(1)
          do j=1,nj
          do i=1,ni
            dum1(i,j,k) = u1(i,j)*alog(zh(1,1,k)/znt(i,j))/alog(sigma0(1)/znt(i,j))
            dum2(i,j,k) = v1(i,j)*alog(zh(1,1,k)/znt(i,j))/alog(sigma0(1)/znt(i,j))
          enddo
          enddo
          call bc2d(dum1(ib,jb,k))
          call bc2d(dum2(ib,jb,k))
#ifdef MPI
          call comm_1s2d_start(dum1(ib,jb,k),sw31(1,1,1),sw32(1,1,1),se31(1,1,1),se32(1,1,1),  &
                                             ss31(1,1,1),ss32(1,1,1),sn31(1,1,1),sn32(1,1,1),reqs_s)
          call comm_1s2d_end(  dum1(ib,jb,k),sw31(1,1,1),sw32(1,1,1),se31(1,1,1),se32(1,1,1),  &
                                             ss31(1,1,1),ss32(1,1,1),sn31(1,1,1),sn32(1,1,1),reqs_s)
          call comm_1s2d_start(dum2(ib,jb,k),sw31(1,1,1),sw32(1,1,1),se31(1,1,1),se32(1,1,1),  &
                                             ss31(1,1,1),ss32(1,1,1),sn31(1,1,1),sn32(1,1,1),reqs_s)
          call comm_1s2d_end(  dum2(ib,jb,k),sw31(1,1,1),sw32(1,1,1),se31(1,1,1),se32(1,1,1),  &
                                             ss31(1,1,1),ss32(1,1,1),sn31(1,1,1),sn32(1,1,1),reqs_s)
#endif
          do j=1,nj+1
          do i=1,ni+1
            ua(i,j,k) = 0.5*(dum1(i-1,j,k)+dum1(i,j,k))
            va(i,j,k) = 0.5*(dum2(i,j-1,k)+dum2(i,j,k))
          enddo
          enddo
        endif
      enddo

      if(myid.eq.0) print *

    ENDIF

!---------------------------------------------------------------
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!---------------------------------------------------------------

!---------

      if( adapt_dt.eq.0 )then
        dt = dtl
        dbldt = dtl
      endif

      ! this is needed for stats files:
      nrec=nrec-1
      nstatout=nstatout-1

!      ! need to think about this
!      IF( ( output_format .ne. old_format )  .or.  &
!          ( output_filetype .ne. old_filetype )   )THEN
!        srec = 1
!        sirec = 1
!        urec = 1
!        vrec = 1
!        wrec = 1
!        nrec = 1
!        mrec = 1
!        nwrite = 1
!        prec = 1
!      ENDIF

!---------

        IF( (imoist.eq.1).and.(ptype.eq.2) )then
          if(timestats.ge.1) time_misc=time_misc+mytime()
          call consat2(dt)
          if(timestats.ge.1) time_microphy=time_microphy+mytime()
        ENDIF
        IF( (imoist.eq.1).and.(ptype.eq.4) )then
          if(timestats.ge.1) time_misc=time_misc+mytime()
          call lfoice_init(dt)
          if(timestats.ge.1) time_microphy=time_microphy+mytime()
        ENDIF

!---------

      if( allocated(xhref0) )then
        deallocate( xhref0 )
      endif
      if( allocated(xfref0) )then
        deallocate( xfref0 )
      endif
      if( allocated(yhref0) )then
        deallocate( yhref0 )
      endif
      if( allocated(yfref0) )then
        deallocate( yfref0 )
      endif
      if( allocated(sigma0) )then
        deallocate( sigma0 )
      endif
      if( allocated(sigmaf0) )then
        deallocate( sigmaf0 )
      endif
      if( allocated(datk1) )then
        deallocate( datk1 )
      endif
      if( allocated(datk2) )then
        deallocate( datk2 )
      endif
      if( allocated(datk3) )then
        deallocate( datk3 )
      endif
      if( allocated(ix) )then
        deallocate( ix )
      endif
      if( allocated(iy) )then
        deallocate( iy )
      endif
      if( allocated(iz) )then
        deallocate( iz )
      endif
      if( allocated(ixu) )then
        deallocate( ixu )
      endif
      if( allocated(iyv) )then
        deallocate( iyv )
      endif
      if( allocated(izw) )then
        deallocate( izw )
      endif
      if( allocated(rx) )then
        deallocate( rx )
      endif
      if( allocated(ry) )then
        deallocate( ry )
      endif
      if( allocated(rz) )then
        deallocate( rz )
      endif
      if( allocated(rxu) )then
        deallocate( rxu )
      endif
      if( allocated(ryv) )then
        deallocate( ryv )
      endif
      if( allocated(rzw) )then
        deallocate( rzw )
      endif
      if( allocated(info) )then
        deallocate( info )
      endif

      if( interp_on_restart )then
        if(myid.eq.0) print *
        if(myid.eq.0) print *,'  setting qpten to zero '
        qpten = 0.0
        qtten = 0.0
        qvten = 0.0
        qcten = 0.0
      endif

      if(dowr) write(outfile,*)
      if(dowr) write(outfile,*) 'Leaving RESTART_READ'
      if(dowr) write(outfile,*)

#ifdef MPI
      ! seems to help with memory on some hardware (somehow):
      call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif

      return

778   print *,'  error opening restart file '
      print *,'    ... stopping cm1 ... '
      call stopcm1

7001  print *
      print *,'  7001: error opening restart files '
      print *
      print *,'    Error attempting to read cm1rversion '
      print *
      print *,'  Note that cm1r22.0 and later versions require restart files '
      print *,'  written by cm1r22.0 or later. '
      print *
      print *,'    ... stopping cm1 ... '
      print *
      call stopcm1

      end subroutine restart_read


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


    subroutine readbc(ndum,nzold,dum,restart_read_format,ncid,aname,interp_on_restart)
#ifdef NETCDF
    use netcdf
#endif
    implicit none

    integer, intent(in) :: ndum,nzold
    real, intent(inout), dimension(0:ndum+1,0:nzold+1) :: dum
    integer, intent(in) :: restart_read_format,ncid
    character(len=6), intent(in) :: aname
    logical, intent(in) :: interp_on_restart

    integer :: k,kk,n
#ifdef NETCDF
    integer :: varid,ncstatus
#endif

    if( restart_read_format.eq.1 )then
      do k=1,nzold
        read(rstfile) dum(1:ndum,k)
      enddo
#ifdef NETCDF
    elseif( restart_read_format.eq.2 )then
      ncstatus = nf90_inq_varid(ncid,aname,varid)
      if(ncstatus.ne.nf90_noerr)then
        print *,'  Error1 in readbc, aname = ',aname
        print *,nf90_strerror(ncstatus)
        call stopcm1
      endif
      stop 33328
#endif
    endif

    if( interp_on_restart )then

      ! interp to new grid:
      !$omp parallel do default(shared)   &
      !$omp private(n)
      do n=1,ndum
        dum(n,0) = dum(n,1)
        dum(n,nzold+1) = dum(n,nzold)
      enddo

      !$omp parallel do default(shared)   &
      !$omp private(k)
      do k=0,nzold+1
        dum(0,k) = dum(1,k)
        dum(ndum+1,k) = dum(ndum,k)
      enddo

    endif

    end subroutine readbc


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


    subroutine  readr(numi,numj,numk1,numk2,nxr,nyr,var,aname,           &
                      ni,nj,ngxy,myid,numprocs,nodex,nodey,orec,nfile,   &
                      ncid,time_index,restart_read_format,restart_read_filetype,myi1p,myi2p,myj1p,myj2p,   &
                      dat1,dat2,dat3,datk1,datk2,datk3,nx,ny,nzr,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma,sigma0,info,reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2i,d2j,d3i,d3j)
#ifdef MPI
    use mpi
#endif
#ifdef NETCDF
    use netcdf
#endif
    implicit none

    !-------------------------------------------------------------------
    ! This subroutine reads restart files and then passes data 
    ! to other processors if this is a MPI run. 
    !-------------------------------------------------------------------

    integer, intent(in) :: numi,numj,numk1,numk2,nxr,nyr,nxold,nyold,nzold,nx,ny,nzr
    integer, intent(in) :: ppnode,d3n,d3t,d2i,d2j,d3i,d3j
    real, intent(inout), dimension(1-ngxy:numi+ngxy,1-ngxy:numj+ngxy,numk1:numk2) :: var
    character(len=8), intent(in) :: aname
    integer, intent(in) :: ni,nj,ngxy,myid,numprocs,nodex,nodey
    integer, intent(inout) :: orec,ncid
    integer, intent(in) :: time_index,restart_read_format,restart_read_filetype
    real, intent(inout), dimension(d3i,d3j) :: dat1
    real, intent(inout), dimension(d2i,d2j) :: dat2
    real, intent(inout), dimension(d3i,d3j,0:d3n-1) :: dat3
    real, intent(inout), dimension(dk1i1:dk1i2,dk1j1:dk1j2) :: datk1
    real, intent(inout), dimension(dk2i1:dk2i2,dk2j1:dk2j2) :: datk2
    real, intent(inout), dimension(dk3i1:dk3i2,dk3j1:dk3j2) :: datk3
    integer, intent(inout), dimension(d3t) :: reqt
    integer, intent(in), dimension(info1,0:info2) :: info
    integer, intent(in) :: mynode,nodeleader,nodes,nfile
    integer, intent(in), dimension(0:numprocs-1) :: myi1p,myi2p,myj1p,myj2p
    real, intent(in), dimension(0:numk2) :: sigma
    real, intent(in), dimension(0:nzold) :: sigma0
    integer, intent(in), dimension(numi) :: ix
    integer, intent(in), dimension(numj) :: iy
    integer, intent(in), dimension(nzr)  :: iz
    real,    intent(in), dimension(numi) :: rx
    real,    intent(in), dimension(numj) :: ry
    real,    intent(in), dimension(nzr)  :: rz

    integer :: i,j,k,msk
#ifdef MPI
    integer :: ierr
#endif

!-------------------------------------------------------------------------------

    msk = 0

    stag = 1
    if( nfile.eq.ufile ) stag = 2
    if( nfile.eq.vfile ) stag = 3

    if(myid.eq.0) print *,aname,numk1,numk2
!!!    if(myid.eq.0) print *,'    stag,nxold,nyold = ',stag,nxold,nyold

    rf1:  &
    IF( restart_read_filetype.eq.1 .or. restart_read_filetype.eq.2 .or. restart_read_filetype.eq.3 .or. restart_read_filetype.eq.4 )THEN

#ifndef MPI

    !  single processor 

    IF( dointerp )THEN
      ! with interpolation:
      call readr_interp(ppnode,d2i,d2j,d3i,d3j,d3n,d3t,ngxy,numi,numj,numk1,numk2,reqt,dat1,dat2,dat3,var,  &
                        nxold,nyold,nzold,nxr,nyr,nzr,ix,iy,iz,rx,ry,rz,restart_read_format,nfile,orec,  &
                        msk,myid,numprocs,mynode,nodes,ni,nj,myi1p,myi2p,myj1p,myj2p,datk1,datk2,datk3,info,restart_read_filetype)
    ELSE
      ! no interpolation:
        call   readr_sp(numk1,numk2,nxold,nyold,nxr,nyr,d2i,d2j,numi,numj,ngxy,dat2,var,orec,restart_read_format,nfile,ncid,aname,time_index,d3i,d3j,dat1,info,restart_read_filetype,datk3)
    ENDIF

#else

    !  multiple processors

    IF( dointerp )THEN
      ! with interpolation:
      call readr_interp(ppnode,d2i,d2j,d3i,d3j,d3n,d3t,ngxy,numi,numj,numk1,numk2,reqt,dat1,dat2,dat3,var,  &
                        nxold,nyold,nzold,nxr,nyr,nzr,ix,iy,iz,rx,ry,rz,restart_read_format,nfile,orec,  &
                        msk,myid,numprocs,mynode,nodes,ni,nj,myi1p,myi2p,myj1p,myj2p,datk1,datk2,datk3,info,restart_read_filetype)
    ELSE
      ! no interpolation:
      if( .not. newdecomp .and. (restart_filetype_old.eq.3.or.restart_filetype_old.eq.4) )then
        ! same grid, same domain decomposition ... read data directly (easy peasy):
        if( restart_filetype_old.eq.3 )then
          do k=numk1,numk2
            read(nfile,rec=orec) ((var(i,j,k),i=1,numi),j=1,numj)
            orec = orec+1
          enddo
        elseif( restart_filetype_old.eq.4 )then
          if(myid.ne.nodeleader)then
            call readr_comm1(numk1,numk2,d3i,d3j,ngxy,numi,numj,nodeleader,dat1,var)
          else
            call readr_comm2(ppnode,d2i,d2j,d3i,d3j,d3n,d3t,ngxy,numi,numj,numk1,numk2,msk,myid,reqt,dat2(1,1),dat3,var,nfile,orec,.false.)
          endif
        endif
      else
        ! same grid, new decomposition ... complicated (aka fun):
        ! ... and, restart_filetype = 1,2 (one file for entire domain)
        IF(myid.ne.nodeleader)THEN
          call   readr_comm1(numk1,numk2,d3i,d3j,ngxy,numi,numj,nodeleader,dat1,var)
        ELSE
          IF(myid.ne.msk)THEN
            call readr_comm2(ppnode,d2i,d2j,d3i,d3j,d3n,d3t,ngxy,numi,numj,numk1,numk2,msk,myid,reqt,dat2(1,1),dat3,var,nfile,orec,.true.)
          ELSE
            call readr_comm3(ppnode,d2i,d2j,d3i,d3j,d3n,d3t,ngxy,numi,numj,numk1,numk2,reqt,dat1,dat2,dat3,var,  &
                             nxold,nyold,nxr,nyr,restart_read_format,nfile,orec,                                 &
                             msk,myid,numprocs,mynode,nodes,ni,nj,myi1p,myi2p,myj1p,myj2p,datk3,                 &
                             ncid,time_index,aname,info,restart_read_filetype)
          ENDIF
        ENDIF
      endif
    ENDIF

#endif


    ENDIF  rf1

!-------------------------------------------------------------------------------
!ccccc  done  cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!-------------------------------------------------------------------------------

    end subroutine  readr


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
#ifndef MPI


    ! code for single processor runs:

    subroutine readr_sp(numk1,numk2,nxold,nyold,nxr,nyr,d2i,d2j,numi,numj,ngxy,dat2,var,orec,restart_read_format,nfile,ncid,aname,time_index,d3i,d3j,dat1,info,restart_read_filetype,datk3)
#ifdef NETCDF
    use netcdf
#endif
    implicit none

    integer, intent(in) :: numk1,numk2,nxold,nyold,nxr,nyr,d2i,d2j,numi,numj,ngxy
    real, intent(inout), dimension(d2i,d2j) :: dat2
    real, intent(inout), dimension(1-ngxy:numi+ngxy,1-ngxy:numj+ngxy,numk1:numk2) :: var
    integer, intent(inout) :: orec,ncid
    integer, intent(in) :: restart_read_format,nfile
    character(len=8), intent(in) :: aname
    integer, intent(in) :: time_index
    integer, intent(in) :: d3i,d3j
    real, intent(inout), dimension(d3i,d3j) :: dat1
    integer, intent(in), dimension(info1,0:info2) :: info
    integer, intent(in) :: restart_read_filetype
    real, intent(inout), dimension(dk3i1:dk3i2,dk3j1:dk3j2) :: datk3

    integer :: i,j,k

#ifdef NETCDF
    integer :: varid,status
#endif

#ifdef NETCDF
    if( restart_read_format.eq.2 )then
      status = nf90_inq_varid(ncid,aname,varid)
      if(status.ne.nf90_noerr)then
        print *,'  Error1 in  readr, aname = ',aname
        print *,nf90_strerror(status)
        call stopcm1
      endif
    endif
#endif


    kloop:  &
    DO k=numk1,numk2

      IF( restart_read_format.eq.1 )THEN
        ! ..... binary format .....
        call read_binary(numi,numj,nfile,nxr,nyr,1,d2i,1,d2j,dat2(1,1),orec,nxold,nyold,info,datk3,restart_read_filetype)
#ifdef NETCDF
      ELSEIF( restart_read_format.eq.2 )THEN
        ! ..... netcdf format .....
        if(numk1.eq.numk2)then
          status = nf90_get_var(ncid,varid,dat2,(/1,1,time_index/),(/d2i,d2j,1/))
        else
          status = nf90_get_var(ncid,varid,dat2,(/1,1,k,time_index/),(/d2i,d2j,1,1/))
        endif
        if(status.ne.nf90_noerr)then
          print *,'  Error2 in  readr, aname = ',aname
          print *,nf90_strerror(status)
          call stopcm1
        endif
#endif
      ENDIF


      ! all done ... put contents of dat2 into var array and move on !

      !$omp parallel do default(shared)   &
      !$omp private(i,j)
      do j=1,numj
      do i=1,numi
        var(i,j,k)=dat2(i,j)
      enddo
      enddo

    ENDDO  kloop



    end subroutine readr_sp


#endif
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
#ifdef MPI


! yo comm1a
    ! ordinary processor ... send data to nodeleader:
    subroutine readr_comm1(numk1,numk2,d3i,d3j,ngxy,numi,numj,nodeleader,dat1,var)
    use mpi
    implicit none

    integer, intent(in) :: numk1,numk2,d3i,d3j,ngxy,numi,numj,nodeleader
    real, intent(inout), dimension(d3i,d3j) :: dat1
    real, intent(inout), dimension(1-ngxy:numi+ngxy,1-ngxy:numj+ngxy,numk1:numk2) :: var

    integer :: i,j,k,tag,reqs,ierr

    tag = 1

    kloop:  &
    DO k=numk1,numk2

      call MPI_IRECV(dat1(1,1),d3i*d3j,MPI_REAL,nodeleader,tag,MPI_COMM_WORLD,reqs,ierr)
      call MPI_WAIT(reqs,MPI_STATUS_IGNORE,ierr)

      !$omp parallel do default(shared)   &
      !$omp private(i,j)
      do j=1,numj
      do i=1,numi
        var(i,j,k)=dat1(i,j)
      enddo
      enddo

      tag = tag+2

    ENDDO  kloop

    end subroutine readr_comm1

#endif
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
#ifdef MPI


    ! nodeleader but not proc 0:
    subroutine readr_comm2(ppnode,d2i,d2j,d3i,d3j,d3n,d3t,ngxy,numi,numj,numk1,numk2,msk,myid,reqt,dat2,dat3,var,nfile,orec,get_from_msk)
    use mpi
    implicit none

    integer, intent(in) :: ppnode,d2i,d2j,d3i,d3j,d3n,d3t,ngxy,numi,numj,numk1,numk2,msk,myid
    integer, intent(inout), dimension(d3t) :: reqt
    real, intent(inout), dimension(d3i*ppnode_old,d3j) :: dat2
    real, intent(inout), dimension(d3i,d3j,0:d3n-1) :: dat3
    real, intent(inout), dimension(1-ngxy:numi+ngxy,1-ngxy:numj+ngxy,numk1:numk2) :: var
    integer, intent(in)    :: nfile
    integer, intent(inout) :: orec
    logical, intent(in)    :: get_from_msk

    integer :: i,j,k,tag,reqs,ierr,proc,fooi
    integer, dimension(mpi_status_size,ppnode-1) :: status1

    tag = 1

    kloop:  &
    DO k=numk1,numk2

      if( get_from_msk )then
        ! get data from msk:
        call MPI_IRECV(dat3(1,1,myid),d3i*d3j*ppnode,MPI_REAL,msk,tag+1,MPI_COMM_WORLD,reqs,ierr)
        ! wait for data to arrive:
        call MPI_WAIT(reqs,MPI_STATUS_IGNORE,ierr)
        do proc=myid+1,myid+(ppnode-1)
          call MPI_ISEND(dat3(1,1,proc),d3i*d3j,MPI_REAL,proc,tag,MPI_COMM_WORLD,reqt(proc-myid),ierr)
        enddo
        !$omp parallel do default(shared)  &
        !$omp private(i,j)
        do j=1,numj
        do i=1,numi
          var(i,j,k)=dat3(i,j,myid)
        enddo
        enddo
      else
        read(nfile,rec=orec) ((dat2(i,j),i=1,d3i*ppnode),j=1,d3j)
        orec = orec+1
        do proc=myid+1,myid+(ppnode-1)
          fooi = d3i*(proc-myid)
          !$omp parallel do default(shared)   &
          !$omp private(i,j)
          do j=1,d3j
          do i=1,d3i
            dat3(i,j,proc)=dat2(fooi+i,j)
          enddo
          enddo
          call MPI_ISEND(dat3(1,1,proc),d3i*d3j,MPI_REAL,proc,tag,MPI_COMM_WORLD,reqt(proc-myid),ierr)
        enddo
        !$omp parallel do default(shared)  &
        !$omp private(i,j)
        do j=1,numj
        do i=1,numi
          var(i,j,k)=dat2(i,j)
        enddo
        enddo
      endif
      ! wait for sends to finish:
      call mpi_waitall(ppnode-1,reqt(1:ppnode-1),status1,ierr)
      ! DONE

      tag = tag+2

    ENDDO  kloop

    end subroutine readr_comm2
    ! done, nodeleader but not proc 0:

! yo comm2b

#endif
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
#ifdef MPI


    ! nodeleader and proc 0:
! yo comm3a
    subroutine   readr_comm3(ppnode,d2i,d2j,d3i,d3j,d3n,d3t,ngxy,numi,numj,numk1,numk2,reqt,dat1,dat2,dat3,var,  &
                             nxold,nyold,nxr,nyr,restart_read_format,nfile,orec,                                 &
                             msk,myid,numprocs,mynode,nodes,ni,nj,myi1p,myi2p,myj1p,myj2p,datk3,                 &
                             ncid,time_index,aname,info,restart_read_filetype)
    use mpi
#ifdef NETCDF
    use netcdf
#endif
    implicit none

    integer, intent(in) :: ppnode,d2i,d2j,d3i,d3j,d3n,d3t,ngxy,numi,numj,numk1,numk2
    integer, intent(inout), dimension(d3t) :: reqt
    real, intent(inout), dimension(d3i,d3j) :: dat1
    real, intent(inout), dimension(d2i,d2j) :: dat2
    real, intent(inout), dimension(d3i,d3j,0:d3n-1) :: dat3
    real, intent(inout), dimension(1-ngxy:numi+ngxy,1-ngxy:numj+ngxy,numk1:numk2) :: var
    integer, intent(in) :: nxold,nyold,nxr,nyr
    integer, intent(in) :: restart_read_format,nfile
    integer, intent(inout) :: orec
    integer, intent(in) :: msk,myid,numprocs,mynode,nodes,ni,nj
    integer, intent(in), dimension(0:numprocs-1) :: myi1p,myi2p,myj1p,myj2p
    real, intent(inout), dimension(dk3i1:dk3i2,dk3j1:dk3j2) :: datk3
    integer, intent(in) :: ncid,time_index
    character(len=8), intent(in) :: aname
    integer, intent(in), dimension(info1,0:info2) :: info
    integer, intent(in) :: restart_read_filetype

    integer :: i,j,k,tag,reqs,ierr
    integer :: fooi,fooj,n,nn,nnn,ntot,n1,n2,nitmp,njtmp,proc,index,index2
    integer, dimension(mpi_status_size,ppnode-1 + nodes-1) :: status2
#ifdef NETCDF
    integer :: varid,status
#endif

#ifdef NETCDF
    if( restart_read_format.eq.2 )then
      status = nf90_inq_varid(ncid,aname,varid)
      if(status.ne.nf90_noerr)then
        print *,'  Error1 in  readr, aname = ',aname
        print *,nf90_strerror(status)
        call stopcm1
      endif
    endif
#endif

    tag = 1

    kloop:  &
    DO k=numk1,numk2

          ! read data:
      IF( restart_read_format.eq.1 )THEN
        ! ..... binary format .....
        call read_binary(numi,numj,nfile,nxr,nyr,1,d2i,1,d2j,dat2(1,1),orec,nxold,nyold,info,datk3,restart_read_filetype)
#ifdef NETCDF
      ELSEIF( restart_read_format.eq.2 )THEN
        ! ..... netcdf format .....
        if(numk1.eq.numk2)then
          status = nf90_get_var(ncid,varid,dat2,(/1,1,time_index/),(/d2i,d2j,1/))
        else
          status = nf90_get_var(ncid,varid,dat2,(/1,1,k,time_index/),(/d2i,d2j,1,1/))
        endif
        if(status.ne.nf90_noerr)then
          print *,'  Error2 in  readr, aname = ',aname
          print *,nf90_strerror(status)
          call stopcm1
        endif
#endif
      ENDIF

          ! send data:
          do nn=1,( nodes-1 )
              ! send data to other nodeleaders:
              index2 = nn
              if( index2.le.mynode )then
                index2 = index2-1
              endif
              n1 = index2*ppnode
              n2 = (index2+1)*ppnode-1
              do nnn=n1,n2
                proc = nnn
                fooi = myi1p(proc)-1
                fooj = myj1p(proc)-1
                nitmp = myi2p(proc)-myi1p(proc)+1
                njtmp = myj2p(proc)-myj1p(proc)+1
                if( numi.gt.ni ) nitmp = nitmp+1
                if( numj.gt.nj ) njtmp = njtmp+1
                !$omp parallel do default(shared)  &
                !$omp private(i,j)
                do j=1,njtmp
                do i=1,nitmp
                  dat3(i,j,proc) = dat2(fooi+i,fooj+j)
                enddo
                enddo
              enddo
              proc = index2*ppnode
              call MPI_ISEND(dat3(1,1,proc),d3i*d3j*ppnode,MPI_REAL,proc,tag+1,MPI_COMM_WORLD,reqt(ppnode-1+nn),ierr)
          enddo
          do nn=1,( ppnode-1 )
              ! send data to ordinary procs on this node:
              proc = myid+nn
              fooi = myi1p(proc)-1
              fooj = myj1p(proc)-1
              nitmp = myi2p(proc)-myi1p(proc)+1
              njtmp = myj2p(proc)-myj1p(proc)+1
              if( numi.gt.ni ) nitmp = nitmp+1
              if( numj.gt.nj ) njtmp = njtmp+1
              !$omp parallel do default(shared)  &
              !$omp private(i,j)
              do j=1,njtmp
              do i=1,nitmp
                dat3(i,j,proc) = dat2(fooi+i,fooj+j)
              enddo
              enddo
              call MPI_ISEND(dat3(1,1,proc),d3i*d3j,MPI_REAL,proc,tag,MPI_COMM_WORLD,reqt(nn),ierr)
          enddo

          if( myid.eq.0 )then
            !$omp parallel do default(shared)  &
            !$omp private(i,j)
            do j=1,numj
            do i=1,numi
              var(i,j,k) = dat2(i,j)
            enddo
            enddo
          else
            proc = myid
            fooi = myi1p(proc)-1
            fooj = myj1p(proc)-1
            nitmp = myi2p(proc)-myi1p(proc)+1
            njtmp = myj2p(proc)-myj1p(proc)+1
            if( numi.gt.ni ) nitmp = nitmp+1
            if( numj.gt.nj ) njtmp = njtmp+1
            !$omp parallel do default(shared)  &
            !$omp private(i,j)
            do j=1,njtmp
            do i=1,nitmp
              var(i,j,k) = dat2(fooi+i,fooj+j)
            enddo
            enddo
          endif
          ntot = ppnode-1 + nodes-1
          call mpi_waitall(ntot,reqt(1:ntot),status2,ierr)


      tag = tag+2

    ENDDO  kloop

    end subroutine readr_comm3
    ! done nodeleader and proc 0
! yo comm3b


#endif
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


    subroutine readr_interp(ppnode,d2i,d2j,d3i,d3j,d3n,d3t,ngxy,numi,numj,numk1,numk2,reqt,dat1,dat2,dat3,var,  &
                           nxold,nyold,nzold,nxr,nyr,nzr,ix,iy,iz,rx,ry,rz,restart_read_format,nfile,orec,  &
                           msk,myid,numprocs,mynode,nodes,ni,nj,myi1p,myi2p,myj1p,myj2p,datk1,datk2,datk3,info,restart_read_filetype)
#ifdef MPI
    use mpi
#endif
#ifdef NETCDF
    use netcdf
#endif
    implicit none

    integer, intent(in) :: ppnode,d2i,d2j,d3i,d3j,d3n,d3t,ngxy,numi,numj,numk1,numk2
    integer, intent(inout), dimension(d3t) :: reqt
    real, intent(inout), dimension(d3i,d3j) :: dat1
    real, intent(inout), dimension(d2i,d2j) :: dat2
    real, intent(inout), dimension(d3i,d3j,0:d3n-1) :: dat3
    real, intent(inout), dimension(1-ngxy:numi+ngxy,1-ngxy:numj+ngxy,numk1:numk2) :: var
    integer, intent(in) :: nxold,nyold,nzold,nxr,nyr,nzr
    integer, intent(in), dimension(numi) :: ix
    integer, intent(in), dimension(numj) :: iy
    integer, intent(in), dimension(nzr)  :: iz
    real,    intent(in), dimension(numi) :: rx
    real,    intent(in), dimension(numj) :: ry
    real,    intent(in), dimension(nzr)  :: rz
    integer, intent(in) :: restart_read_format,nfile
    integer, intent(inout) :: orec
    integer, intent(in) :: msk,myid,numprocs,mynode,nodes,ni,nj
    integer, intent(in), dimension(0:numprocs-1) :: myi1p,myi2p,myj1p,myj2p
    real, intent(inout), dimension(dk1i1:dk1i2,dk1j1:dk1j2) :: datk1
    real, intent(inout), dimension(dk2i1:dk2i2,dk2j1:dk2j2) :: datk2
    real, intent(inout), dimension(dk3i1:dk3i2,dk3j1:dk3j2) :: datk3
    integer, intent(in), dimension(info1,0:info2) :: info
    integer, intent(in) :: restart_read_filetype

    integer :: i,j,k,tag,reqs,ierr,kold1,kold2
    integer :: fooi,fooj,n,nn,nnn,ntot,n1,n2,nitmp,njtmp,proc,index,index2
    logical :: alldone
#ifdef MPI
    integer, dimension(mpi_status_size,ppnode-1) :: status1
    integer, dimension(mpi_status_size,ppnode-1 + nodes-1) :: status2
#endif
#ifdef NETCDF
    integer :: varid,status
#endif

#ifdef NETCDF
!    if( myid.eq.0 )then
!    if( restart_read_format.eq.2 )then
!      status = nf90_inq_varid(ncid,aname,varid)
!      if(status.ne.nf90_noerr)then
!        print *,'  Error1 in  readr, aname = ',aname
!        print *,nf90_strerror(status)
!        call stopcm1
!      endif
!    endif
!    endif
#endif

    tag = 1
    kold1 = -1
    kold2 =  0
    alldone = .false.

    kloop:  &
    DO k=numk1,numk2

    if(myid.eq.msk)then
!!!      print *,'  k = ',k
          ! read data:
      IF( restart_read_format.eq.1 )THEN
        ! ..... binary format .....
          call   prep_for_interp(k,numi,numj,nzr,nxold,nyold,nzold,numk1,numk2,nfile,myid,iz,kold1,kold2,orec,alldone,datk1,datk2,datk3,d3i,d3j,dat1,info,restart_read_filetype)
#ifdef NETCDF
      ELSEIF( restart_read_format.eq.2 )THEN
        ! ..... netcdf format .....
        stop 12871
#endif
      ENDIF
    endif

#ifdef MPI
        call MPI_BCAST(datk1,(dk1i2-dk1i1+1)*(dk1j2-dk1j1+1),MPI_REAL,msk,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(datk2,(dk2i2-dk2i1+1)*(dk2j2-dk2j1+1),MPI_REAL,msk,MPI_COMM_WORLD,ierr)
#endif

      if( dointerp )then
        call     do_interp(k,myid,nxr,nyr,nzr,nxold,nyold,d2i,d2j,ix,iy,iz,rx,ry,rz,datk1,datk2,numi,numj,numk1,numk2,ngxy,var)
      endif

    ENDDO  kloop

!!!    if(myid.eq.0)then
!!!      print *,'  numk1,numk2 = ',numk1,numk2
!!!      print *,'  nkold,kold2 = ',nkold,kold2
!!!    endif

    ! if new domain is shallower than old domain:
    if( (numk2.gt.numk1) .and. (kold2.lt.nzold) )then
      if(myid.eq.0) print *,'  incrementing orec:',nzold-kold2
      orec = orec + (nzold-kold2)
    endif

    end subroutine readr_interp


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      ! only proc 0 calls this 
    subroutine readr_newdecomp(nxold,nyold,numi,numj,i1,i2,j1,j2,dat2,datk3,orec,info,skip,restart_read_filetype)
    implicit none

    integer, intent(in) :: nxold,nyold,numi,numj,i1,i2,j1,j2,skip
    real, intent(inout), dimension(i1:i2,j1:j2) :: dat2
    real, intent(inout), dimension(dk3i1:dk3i2,dk3j1:dk3j2) :: datk3
    integer, intent(in) :: orec
    integer, intent(in), dimension(info1,0:info2) :: info
    integer, intent(in) :: restart_read_filetype

    integer :: ierr,i,j,n,nfile,ifoo,jfoo,nitmp,njtmp,np

    if( restart_read_filetype.eq.3 )then

      do n=0,nfiles-1
        nfile = 1000 + skip*nfiles + n
        ifoo = info(1,n)-1
        jfoo = info(2,n)-1
        nitmp = info(3,n)
        njtmp = info(4,n)
        if(skip.eq.3) nitmp = nitmp+1
        if(skip.eq.4) njtmp = njtmp+1
        call readfile_ft3(nfile,orec,ifoo,jfoo,nxold,nyold,nitmp,njtmp,i1,i2,j1,j2,dat2(i1,j1),datk3(dk3i1,dk3j1))
      enddo

    elseif( restart_read_filetype.eq.4 )then

      np = 0
      nitmp = nioldmax
      njtmp = njoldmax
      if(skip.eq.3) nitmp = nitmp+1
      if(skip.eq.4) njtmp = njtmp+1

      do n=0,nfiles-1
        nfile = 1000 + skip*nfiles + n
        call readfile_ft4(nfile,orec,skip,nxold,nyold,nitmp,njtmp,i1,i2,j1,j2,dat2(i1,j1),datk3(dk3i1,dk3j1),np,info)
      enddo

    endif

    end subroutine readr_newdecomp

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    subroutine readfile_ft3(nfile,orec,ifoo,jfoo,nxold,nyold,nitmp,njtmp,i1,i2,j1,j2,dat2,dat3)
    implicit none

    integer, intent(in) :: nfile,orec,ifoo,jfoo,nxold,nyold,nitmp,njtmp,i1,i2,j1,j2
    real, intent(inout), dimension(i1:i2,j1:j2) :: dat2
    real, intent(inout), dimension(nitmp,njtmp) :: dat3

    integer :: i,j

    read(nfile,rec=orec) ((dat3(i,j),i=1,nitmp),j=1,njtmp)

    !$omp parallel do default(shared)   &
    !$omp private(i,j)
    do j=1,njtmp
    do i=1,nitmp
      dat2(ifoo+i,jfoo+j) = dat3(i,j)
    enddo
    enddo

    end subroutine readfile_ft3

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    subroutine readfile_ft4(nfile,orec,skip,nxold,nyold,nitmp,njtmp,i1,i2,j1,j2,dat2,datk3,np,info)
    implicit none

    integer, intent(in) :: nfile,orec,skip,nxold,nyold,nitmp,njtmp,i1,i2,j1,j2
    real, intent(inout), dimension(i1:i2,j1:j2) :: dat2
    real, intent(inout), dimension(1:nitmp*ppnode_old,1:njtmp) :: datk3
    integer, intent(inout) :: np
    integer, intent(in), dimension(info1,0:info2) :: info

    integer :: i,j,n,itmp,jtmp,ni2,nj2

    read(nfile,rec=orec) ((datk3(i,j),i=1,nitmp*ppnode_old),j=1,njtmp)

    do n=0,ppnode_old-1
      itmp = info(1,np)-1
      jtmp = info(2,np)-1
      ni2 = info(3,np)
      nj2 = info(4,np)
      if(skip.eq.3) ni2 = ni2+1
      if(skip.eq.4) nj2 = nj2+1
      !$omp parallel do default(shared)   &
      !$omp private(i,j)
      do j=1,nj2
      do i=1,ni2
        dat2(itmp+i,jtmp+j) = datk3(n*nitmp+i,j)
      enddo
      enddo
      np = np+1
    enddo

    end subroutine readfile_ft4

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine read_binary(numi,numj,nfile,nxr,nyr,i1,i2,j1,j2,dat2,orec,nxold,nyold,info,datk3,restart_read_filetype)
      implicit none

      integer, intent(in) :: numi,numj,nfile,nxr,nyr,i1,i2,j1,j2
      real, intent(inout), dimension(i1:i2,j1:j2) :: dat2
      integer, intent(inout) :: orec
      integer, intent(in) :: nxold,nyold
      integer, intent(in), dimension(info1,0:info2) :: info
      real, intent(inout), dimension(dk3i1:dk3i2,dk3j1:dk3j2) :: datk3
      integer, intent(in) :: restart_read_filetype

      integer :: i,j,n,skip

      if( restart_read_filetype.eq.1 .or. restart_read_filetype.eq.2 )then

        read(nfile,rec=orec) ((dat2(i,j),i=1,nxr),j=1,nyr)
        orec = orec+1

      elseif( restart_read_filetype.eq.3 .or. restart_read_filetype.eq.4 )then

        skip = 2 + nfile-sfile
        call   readr_newdecomp(nxold,nyold,numi,numj,i1,i2,j1,j2,dat2,datk3,orec,info,skip,restart_read_filetype)
        orec = orec+1

      else

        stop 13097

      endif

      end subroutine read_binary


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


    subroutine setbc(dki1,dki2,dkj1,dkj2,nxold,nyold,datk)
    use input, only : wbc,ebc,sbc,nbc
    implicit none

    integer, intent(in) :: dki1,dki2,dkj1,dkj2,nxold,nyold
    real, intent(inout), dimension(dki1:dki2,dkj1:dkj2) :: datk

    integer :: i,j

    !----------

    if( wbc.eq.1 )then
      ! periodic:
      if( stag.eq.2 )then
        !$omp parallel do default(shared)   &
        !$omp private(j)
        do j=0,nyold1+2
          datk(0,j) = datk(nxold1,j)
        enddo
      else
        !$omp parallel do default(shared)   &
        !$omp private(j)
        do j=0,nyold1+2
          datk(0,j) = datk(nxold1,j)
        enddo
      endif
    else
      ! not periodic:
      if( stag.eq.2 )then
        !$omp parallel do default(shared)   &
        !$omp private(j)
        do j=0,nyold1+2
          datk(0,j) = datk(1,j)
        enddo
      else
        !$omp parallel do default(shared)   &
        !$omp private(j)
        do j=0,nyold1+2
          datk(0,j) = datk(1,j)
        enddo
      endif
    endif

    !----------

    if( ebc.eq.1 )then
      ! periodic:
      if( stag.eq.2 )then
        !$omp parallel do default(shared)   &
        !$omp private(j)
        do j=0,nyold1+2
          datk(nxold1+2,j) = datk(2,j)
        enddo
      else
        !$omp parallel do default(shared)   &
        !$omp private(j)
        do j=0,nyold1+2
          datk(nxold1+1,j) = datk(1,j)
        enddo
      endif
    else
      ! not periodic:
      if( stag.eq.2 )then
        !$omp parallel do default(shared)   &
        !$omp private(j)
        do j=0,nyold1+2
          datk(nxold1+2,j) = datk(nxold1+1,j)
        enddo
      else
        !$omp parallel do default(shared)   &
        !$omp private(j)
        do j=0,nyold1+2
          datk(nxold1+1,j) = datk(nxold1,j)
        enddo
      endif
    endif

    !----------

    if( sbc.eq.1 )then
      ! periodic:
      if( stag.eq.3 )then
        !$omp parallel do default(shared)   &
        !$omp private(i)
        do i=0,nxold1+2
          datk(i,0) = datk(i,nyold1)
        enddo
      else
        !$omp parallel do default(shared)   &
        !$omp private(i)
        do i=0,nxold1+2
          datk(i,0) = datk(i,nyold1)
        enddo
      endif
    else
      ! not periodic:
      if( stag.eq.3 )then
        !$omp parallel do default(shared)   &
        !$omp private(i)
        do i=0,nxold1+2
          datk(i,0) = datk(i,1)
        enddo
      else
        !$omp parallel do default(shared)   &
        !$omp private(i)
        do i=0,nxold1+2
          datk(i,0) = datk(i,1)
        enddo
      endif
    endif

    !----------

    if( nbc.eq.1 )then
      ! periodic:
      if( stag.eq.3 )then
        !$omp parallel do default(shared)   &
        !$omp private(i)
        do i=0,nxold1+2
          datk(i,nyold1+2) = datk(i,2)
        enddo
      else
        !$omp parallel do default(shared)   &
        !$omp private(i)
        do i=0,nxold1+2
          datk(i,nyold1+1) = datk(i,1)
        enddo
      endif
    else
      ! not periodic:
      if( stag.eq.3 )then
        !$omp parallel do default(shared)   &
        !$omp private(i)
        do i=0,nxold1+2
          datk(i,nyold1+2) = datk(i,nyold1+1)
        enddo
      else
        !$omp parallel do default(shared)   &
        !$omp private(i)
        do i=0,nxold1+2
          datk(i,nyold1+1) = datk(i,nyold1)
        enddo
      endif
    endif

    !----------

    end subroutine setbc


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine prep_for_interp(k,numi,numj,nzr,nxold,nyold,nzold,numk1,numk2,nfile,myid,iz,kold1,kold2,orec,alldone,datk1,datk2,datk3,d3i,d3j,dat1,info,restart_read_filetype)
      implicit none

      integer, intent(in) :: k,numi,numj,nzr,nxold,nyold,nzold,numk1,numk2,nfile,myid
      integer, intent(in), dimension(nzr) :: iz
      integer, intent(inout) :: kold1,kold2,orec
      logical, intent(inout) :: alldone
      real, intent(inout), dimension(dk1i1:dk1i2,dk1j1:dk1j2) :: datk1
      real, intent(inout), dimension(dk2i1:dk2i2,dk2j1:dk2j2) :: datk2
      real, intent(inout), dimension(dk3i1:dk3i2,dk3j1:dk3j2) :: datk3
      integer, intent(in) :: d3i,d3j
      real, intent(inout), dimension(d3i,d3j) :: dat1
      integer, intent(in), dimension(info1,0:info2) :: info
      integer, intent(in) :: restart_read_filetype

      integer :: i,j,kk,kk1,kk2,orectmp

    !   Layout:
    !
    !    old grid     new grid     old grid
    !      kold1                   kold2
    !        iz          k          iz+1
    !

!!!            print *,'    iz(k)+1,kold2,nzold = ',iz(k)+1,kold2,nzold

            ! check to see if we need to increment kold:
            neednew:  &
            IF( iz(k)+1.gt.kold2 )THEN

              ! yes, we need a new level

              ! is more data available?
              inside_old_domain:  &
              IF( (iz(k)+1).le.nzold )THEN

                  ! yes, more levels are available ... so, get to work:

                  ! determine how many levels to read from the old grid:
                  if( numk1.eq.numk2 )then
                    ! 2D variable:
                    kk1 = 1
                    kk2 = 1
                  else
                    ! 3D variable:
                    kk1 = kold2+1
                    kk2 = iz(k)+1
                  endif

!!!                  print *,'    kk1,kk2 = ',kk1,kk2

                  ! read through old levels:
                  do kk=kk1,kk2

                    kold1 = kold1+1
                    kold2 = kold2+1

                    ! save old level:
                    if( numk2.gt.numk1 )then
                      !$omp parallel do default(shared)   &
                      !$omp private(i,j)
                      do j=0,nyold1+2
                      do i=0,nxold1+2
                        datk1(i,j) = datk2(i,j)
                      enddo
                      enddo
                    endif

!!!                    print *,'      ...... reading kold ...... ',kk
                    call read_binary(numi,numj,nfile,nxold,nyold,dk2i1,dk2i2,dk2j1,dk2j2,datk2,orec,nxold,nyold,info,datk3,restart_read_filetype)
                    call setbc(dk2i1,dk2i2,dk2j1,dk2j2,nxold,nyold,datk2)

                  enddo

                  if( numk2.eq.numk1 )then
                    !$omp parallel do default(shared)   &
                    !$omp private(i,j)
                    do j=0,nyold1+2
                    do i=0,nxold1+2
                      datk1(i,j) = datk2(i,j)
                    enddo
                    enddo
                  endif

                  if( (kold1.eq.0) .and. (numk2.gt.numk1) )then
                    ! we just read the first level ... extrapolate to lower levels

!!!                    print *,'    below first level: '
                    ! (NOTE: need to replace with extrapolation)
                    ! temporarily:  k=1 in datk2,  k=2 in datk1
                    orectmp = orec+1
                    call read_binary(numi,numj,nfile,nxold,nyold,dk1i1,dk1i2,dk1j1,dk1j2,datk1,orectmp,nxold,nyold,info,datk3,restart_read_filetype)

                    ! someday: need actual model levels here (zh arrays) !

                    ! now extrapolate to z=1
!!!                    print *,'  nzr = ',nzr
!!!                    print *,'  z1,z2 = ',z1,z2
                    !$omp parallel do default(shared)   &
                    !$omp private(i,j)
                    do j=1,nyold
                    do i=1,nxold
                      datk1(i,j) = datk2(i,j) - z1*(datk1(i,j)-datk2(i,j))/(z2-z1)
                    enddo
                    enddo

                    call setbc(dk1i1,dk1i2,dk1j1,dk1j2,nxold,nyold,datk1)

                  endif


              ELSE  inside_old_domain

!!!      print *,'  here '

                  ! we have already read all levels from old grid:

                  if( .not. alldone )then
                    ! First time here ... equate the two old levels (for now)

                    ! (NOTE: need to replace with extrapolation)

!!!                    print *,'    above last level: '
                    !$omp parallel do default(shared)   &
                    !$omp private(i,j)
                    do j=0,nyold1+2
                    do i=0,nxold1+2
                      datk1(i,j) = datk2(i,j)
                    enddo
                    enddo

                    call setbc(dk1i1,dk1i2,dk1j1,dk1j2,nxold,nyold,datk1)

                  endif

                  alldone = .true.

              ENDIF  inside_old_domain

            ENDIF  neednew

      end subroutine prep_for_interp


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine do_interp(k,myid,nxr,nyr,nzr,nxold,nyold,d2i,d2j,ix,iy,iz,rx,ry,rz,datk1,datk2,numi,numj,numk1,numk2,ngxy,var)
      implicit none

      integer, intent(in) :: k,myid,nxr,nyr,nzr,nxold,nyold,d2i,d2j
      integer, intent(in), dimension(numi) :: ix
      integer, intent(in), dimension(numj) :: iy
      integer, intent(in), dimension(nzr)  :: iz
      real, intent(in), dimension(numi) :: rx
      real, intent(in), dimension(numj) :: ry
      real, intent(in), dimension(nzr)  :: rz
      real, intent(inout), dimension(dk1i1:dk1i2,dk1j1:dk1j2) :: datk1
      real, intent(inout), dimension(dk2i1:dk2i2,dk2j1:dk2j2) :: datk2
      integer, intent(in) :: numi,numj,numk1,numk2,ngxy
      real, intent(inout), dimension(1-ngxy:numi+ngxy,1-ngxy:numj+ngxy,numk1:numk2) :: var

      double precision :: w1,w2,w3,w4,w5,w6,w7,w8,wsum

      integer :: i,j

        ! this is where the interpolation happens:

        !$omp parallel do default(shared)   &
        !$omp private(i,j,w1,w2,w3,w4,w5,w6,w7,w8,wsum)
        do j=1,numj
        do i=1,numi

          w1=(1.0d0-dble(rx(i)))*(1.0d0-dble(ry(j)))*(1.0d0-dble(rz(k)))
          w2=dble(rx(i))*(1.0d0-dble(ry(j)))*(1.0d0-dble(rz(k)))
          w3=(1.0d0-dble(rx(i)))*dble(ry(j))*(1.0d0-dble(rz(k)))
          w4=(1.0d0-dble(rx(i)))*(1.0d0-dble(ry(j)))*dble(rz(k))
          w5=dble(rx(i))*(1.0d0-dble(ry(j)))*dble(rz(k))
          w6=(1.0d0-dble(rx(i)))*dble(ry(j))*dble(rz(k))
          w7=dble(rx(i))*dble(ry(j))*(1.0d0-dble(rz(k)))
          w8=dble(rx(i))*dble(ry(j))*dble(rz(k))
          wsum = w1+w2+w3+w4+w5+w6+w7+w8

!          ! debug !
!          if( rx(i).lt.-0.0001 .or. rx(i).gt.1.0001 .or.  &
!              ry(j).lt.-0.0001 .or. ry(j).gt.1.0001 .or.  &
!              rz(k).lt.-0.0001 .or. rz(k).gt.1.0001 .or.  &
!              wsum.le.0.99999 .or.                  &
!              wsum.ge.1.00001 )then
!            print *,'  ack:',rx(i),ry(j),rz(k),wsum
!            print *,'  a2:',k,rz(k)
!            print *,'  245987 '
!            call stopcm1
!          endif
!          ! debug !

          var(i,j,k)=( ( ( dble(datk1(ix(i)  ,iy(j)  ))*w1       &
                          +dble(datk1(ix(i)+1,iy(j)+1))*w7 )     &
                        +( dble(datk1(ix(i)+1,iy(j)  ))*w2       &
                          +dble(datk1(ix(i)  ,iy(j)+1))*w3 ) )   &
                      +( ( dble(datk2(ix(i)  ,iy(j)  ))*w4       &
                          +dble(datk2(ix(i)+1,iy(j)+1))*w8 )     &
                        +( dble(datk2(ix(i)+1,iy(j)  ))*w5       &
                          +dble(datk2(ix(i)  ,iy(j)+1))*w6 ) ) )

        enddo
        enddo


      end subroutine do_interp


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    subroutine readinfo(in_file,myid,oorecs,oorecu,oorecv,oorecw,text,num)
#ifdef MPI
    use mpi
#endif
    implicit none

    integer, intent(out) :: in_file
    integer, intent(in) :: myid
    integer, intent(inout) :: oorecs,oorecu,oorecv,oorecw
    character(len=*), intent(in) :: text
    integer, intent(out), optional :: num

    integer :: ierr

    if(doinfo)then
      read(infofile) in_file
      read(infofile) oorecs
      read(infofile) oorecu
      read(infofile) oorecv
      read(infofile) oorecw
      if(myid.eq.0)then
      print *,'  inside readinfo:'
      print *,'    text    = ',trim(text)
      print *,'    in_file = ',in_file
      print *,'    oorecs  = ',oorecs
      print *,'    oorecu  = ',oorecu
      print *,'    oorecv  = ',oorecv
      print *,'    oorecw  = ',oorecw
      endif
      if(present(num))then
        read(infofile) num
        if(myid.eq.0) print *,'    num     = ',num
      endif
    endif
#ifdef MPI
    call MPI_BCAST(in_file  ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    if(present(num))  &
    call MPI_BCAST(num      ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    if( dobcast )then
      call MPI_BCAST(oorecs   ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(oorecu   ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(oorecv   ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(oorecw   ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    endif
#endif

    end subroutine readinfo

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    subroutine pdef(i2,j2,k2,var)
    use input, only : ngxy,ngz
    implicit none

    integer, intent(in) :: i2,j2,k2
    real, intent(inout), dimension(1-ngxy:i2+ngxy,1-ngxy:j2+ngxy,1-ngz:k2+ngz) :: var

    integer :: i,j,k

    !$omp parallel do default(shared)   &
    !$omp private(i,j,k)
    do k=1,k2
    do j=1,j2
    do i=1,i2
      var(i,j,k) = max( 0.0 , var(i,j,k) )
    enddo
    enddo
    enddo

    end subroutine pdef

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  END MODULE restart_read_module
