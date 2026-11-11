  MODULE restart_write_module

  implicit none

  private
  public :: restart_write

  integer, parameter :: rstfile  = 50
  integer, parameter :: sfile    = 51
  integer, parameter :: ufile    = 52
  integer, parameter :: vfile    = 53
  integer, parameter :: wfile    = 54
  integer, parameter :: infofile = 55

  integer, parameter :: in_file_no   =   0
  integer, parameter :: in_file_yes  =   1

  logical :: dorstfile

  CONTAINS

      subroutine restart_write(nstep,srec,sirec,urec,vrec,wrec,nrec,mrec,prec,      &
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
                               haildata,hailloc,ppx,                                &
                               tdiag,qdiag,phi1,phi2,                               &
                               pbl3d,pbl2d,pbl2dint, &
                               thpten,qvpten,qcpten,qipten,upten,vpten,qnipten,qncpten, &
                               icenter,jcenter,xcenter,ycenter,domainlocx,domainlocy,adaptmovetim,mvrec,nwritemv,ug,vg, &
                               gamk,kmw,ufw,vfw,u1b,v1b,                                &
                               ntavg,rtavg,tavg,timavg,sfctavg,sfctimavg,               &
                               dum1,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,confignum)
          ! end_restart_write
      use input
      use lsnudge_module
#ifdef MPI
      use mpi
#endif
#ifdef NETCDF
      use netcdf
      use writeout_nc_module, only : restart_prelim,disp_err
#endif
      implicit none

      !----------------------------------------------------------
      ! This subroutine organizes the writing of restart files
      !----------------------------------------------------------

      integer, intent(in) :: nstep,srec,sirec,urec,vrec,wrec,nrec,mrec,prec,trecs,trecw,arecs,arecw
      integer, intent(in) :: nwrite,nwritet,nwritea,nwriteh,nrst,nstatout
      double precision, intent(in) :: restart_time
      logical, intent(in) :: spinup_t2p
      integer, intent(in) :: num_soil_layers,nrad2d
      double precision, intent(in)    :: avgsfcu,avgsfcv,avgsfcs,avgsfcsu,avgsfcsv,avgsfct,avgsfcq,avgsfcp
      real, intent(in) :: dt,dtlast
      integer, intent(in) :: ndt
      double precision, intent(in) :: adt,adtlast,acfl,dbldt
      double precision, intent(in) :: mass1
      double precision, intent(in) :: mtime,stattim,taptim,rsttim,radtim,prcltim,hailtrajtim
      double precision, intent(inout), dimension(nbudget) :: qbudget
      double precision, intent(inout), dimension(numq) :: asq,bsq
      character(len=3), intent(in), dimension(maxq) :: qname
      real, dimension(1-ngxy:nx+ngxy+1) :: xfref,xhref
      real, dimension(1-ngxy:ny+ngxy+1) :: yfref,yhref
      real, intent(in), dimension(ib:ie) :: xh
      real, intent(in), dimension(ib:ie+1) :: xf
      real, intent(in), dimension(jb:je) :: yh
      real, intent(in), dimension(jb:je+1) :: yf
      real, intent(in), dimension(kb:ke) :: sigma
      real, intent(in), dimension(kb:ke+1) :: sigmaf
      real, intent(in), dimension(ib:ie,jb:je) :: zs
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: zh
      real, intent(in), dimension(ib:ie,jb:je,kb:ke+1) :: zf
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: th0,prs0,pi0,rho0,qv0
      real, intent(in), dimension(ib:ie+1,jb:je,kb:ke) :: u0
      real, intent(in), dimension(ib:ie,jb:je+1,kb:ke) :: v0
      real, intent(in), dimension(ib:ie,jb:je,nrain) :: rain,raincv,hail,sws,svs,sps,srs,sgs,sus,shs
      real, intent(in), dimension(ib:ie,jb:je) :: tsk,znt,ust,cd,ch,cq,u1,v1,s1,th1,qv1,xland,psfc,thflux,qvflux,prate,ustt,ut,vt,st
      real, intent(in), dimension(jb:je,kb:ke) :: radbcw,radbce
      real, intent(in), dimension(ib:ie,kb:ke) :: radbcs,radbcn
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: rho,prs
      real, intent(in),    dimension(ib:ie+1,jb:je,kb:ke) :: ua
      real, intent(inout), dimension(ib:ie+1,jb:je,kb:ke) :: dumu
      real, intent(in),    dimension(ib:ie,jb:je+1,kb:ke) :: va
      real, intent(inout), dimension(ib:ie,jb:je+1,kb:ke) :: dumv
      real, intent(in), dimension(ib:ie,jb:je,kb:ke+1) :: wa
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: ppi,tha
      real, intent(in), dimension(ibm:iem,jbm:jem,kbm:kem,numq) :: qa
      real, intent(in), dimension(ibt:iet,jbt:jet,kbt:ket) :: tkea
      real, intent(in), dimension(ibr:ier,jbr:jer,kbr:ker) :: swten,lwten,thraten,effc,effi,effs,effr,effg,effis
      real, intent(in), dimension(ni,nj) :: radsw,rnflx,radswnet,radlwin
      real, intent(in), dimension(ni,nj,nrad2d) :: rad2d
      real, intent(in), dimension(ibcu:iecu,jbcu:jecu,kbcu:kecu) :: thcuten,qvcuten,qccuten,qicuten,ucuten,vcuten
      integer, intent(in), dimension(ibl:iel,jbl:jel) :: lu_index
      integer, intent(in), dimension(ibl:iel,jbl:jel) :: kpbl2d
      real, intent(in), dimension(ibl:iel,jbl:jel) :: u10,v10,s10,hfx,qfx,    &
                                      hpbl,wspd,psim,psih,gz1oz0,br,          &
                                      CHS,CHS2,CQS2,CPMM,ZOL,MAVAIL,          &
                                      MOL,RMOL,REGIME,LH,FLHC,FLQC,QGH,   &
                                      CK,CKA,CDA,USTM,QSFC,T2,Q2,TH2,EMISS,THC,ALBD,   &
                                      gsw,glw,chklowq,capg,snowc,fm,fh,mznt,swspd,wstar,delta
      real, intent(in), dimension(ibl:iel,jbl:jel,num_soil_layers) :: tslb
      real, intent(in), dimension(ibl:iel,jbl:jel) :: tmn,tml,t0ml,hml,h0ml,huml,hvml,tmoml
      real, intent(in), dimension(ibm:iem,jbm:jem,kbm:kem) :: qpten,qtten,qvten,qcten
      real, intent(in), dimension(ibp:iep,jbp:jep,kbp:kep,npt) :: pta
      real, intent(in), dimension(nparcels,npvals) :: pdata
      real, intent(inout), dimension(nparcels,3) :: ploc
      real, intent(in), dimension(nhailtrajs,nhvals) :: haildata
      real, intent(inout), dimension(nhailtrajs,13) :: hailloc
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: ppx
      real, intent(in),    dimension(ibph:ieph,jbph:jeph,kbph:keph) :: phi1,phi2
      real, intent(in),    dimension(ibpbl:iepbl,jbpbl:jepbl,kbpbl:kepbl,npbl3d) :: pbl3d
      real, intent(inout), dimension(ibpbl:iepbl,jbpbl:jepbl,npbl2d) :: pbl2d
      integer, intent(inout), dimension(ibpbl:iepbl,jbpbl:jepbl,npbl2dint) :: pbl2dint
      real, intent(inout), dimension(ibb:ieb,jbb:jeb,kbb:keb) :: thpten,qvpten,qcpten,qipten,upten,vpten,qnipten,qncpten
      real, intent(in   ) , dimension(ibdt:iedt,jbdt:jedt,kbdt:kedt,ntdiag) :: tdiag
      real, intent(in   ) , dimension(ibdq:iedq,jbdq:jedq,kbdq:kedq,nqdiag) :: qdiag
      integer, intent(in   ) :: icenter,jcenter
      real, intent(in   ) :: xcenter,ycenter
      double precision, intent(in) :: domainlocx,domainlocy
      double precision, intent(in   ) :: adaptmovetim
      integer, intent(in   ) :: mvrec,nwritemv
      real, intent(in   ), dimension(kb:ke) :: ug,vg
      real, intent(in   ), dimension(kb:ke) :: gamk,kmw,ufw,vfw,u1b,v1b
      integer, intent(in   ), dimension(ntim) :: ntavg
      double precision, intent(in   ), dimension(ntim) :: rtavg
      real, intent(in   ), dimension(ibta:ieta,jbta:jeta,kbta:keta,ntim,ntavr) :: tavg
      real, intent(in   ), dimension(ibta:ieta,jbta:jeta,kbta:keta,ntavr) :: timavg
      real, intent(in   ), dimension(ibta:ieta,jbta:jeta,ntim,nsfctavr) :: sfctavg
      real, intent(in   ), dimension(ibta:ieta,jbta:jeta,nsfctavr) :: sfctimavg
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: dum1
      real, intent(inout), dimension(d3i,d3j) :: dat1
      real, intent(inout), dimension(d2i,d2j) :: dat2
      real, intent(inout), dimension(d3i,d3j,d3n) :: dat3
      integer, intent(inout), dimension(d3t) :: reqt
      integer, intent(in), dimension(0:numprocs-1) :: myi1p,myi2p,myj1p,myj2p
      integer, intent(in) :: confignum

      character(len=maxstring) :: fname
      character(len=8) :: text1
      character(len=6) :: aname
      integer :: i,j,k,n,np,nt,nvar,reqs,orecs,orecu,orecv,orecw,ndum
      integer :: ncid,time_index
      real, dimension(:,:), allocatable :: dumx,dumy
#ifdef MPI
      integer :: proc,index,count,req1,req2,req3,reqp,ierr
      double precision, dimension(nbudget) :: cfoo
      double precision, dimension(numq) :: afoo,bfoo
#endif
#ifdef NETCDF
      integer :: varid,ncstatus
#endif
      integer :: in_file,num
      logical :: openrst,opens,openu,openv,openw,openinfo
      logical :: doinfo

!-----------------------------------------------------------------------

  openrst  = .false.
  opens    = .false.
  openu    = .false.
  openv    = .false.
  openw    = .false.
  openinfo = .false.
  if(myid.eq.0)then
    doinfo = .true.
    dorstfile = .true.
  else
    doinfo = .false.
    dorstfile = .false.
  endif

101 format(i6.6)

  IF( restart_format.eq.1 )THEN
    ! unformatted direct-access (grads) format:

  IF( restart_filetype.eq.1 )THEN

    !------------------
    ! one restart file (per stagger type):
    IF(myid.eq.nodeleader)THEN
      if(dowr) write(outfile,*)
      if(dowr) write(outfile,*) '  restart_filetype = 1 is not available with restart_format = 1 '
      if(dowr) write(outfile,*)
      if(dowr) write(outfile,*) '      ... stopping cm1 ... '
      if(dowr) write(outfile,*)
      call stopcm1
    ENDIF

  ELSEIF( restart_filetype.eq.2 )THEN

    !------------------
    ! one restart file (per restart time):
    IF(myid.eq.0)THEN
      do i=1,maxstring
        fname(i:i) = ' '
      enddo

      fname = 'cm1rst_tXXXXXX_i.dat'
      write(fname( 9:14),101) nrst
      if( myid.eq.0 )  &
      open(unit=infofile,file=fname,form='unformatted',status='unknown')
      openinfo = .true.

      fname = 'cm1rst_tXXXXXX_x.dat'
      write(fname( 9:14),101) nrst

      if(dowr) write(outfile,*)
      if(dowr) write(outfile,*) '  Writing to restart file!'
      if(dowr) write(outfile,*) '  fname=',fname

      if( myid.eq.0 )  &
      open(unit=rstfile,file=fname,form='unformatted',status='unknown')
      openrst  = .true.

      fname = 'cm1rst_tXXXXXX_s.dat'
      write(fname( 9:14),101) nrst
      if(dowr) write(outfile,*) '  fname=',fname
      open(unit=sfile,file=fname,form='unformatted',access='direct',recl=4*nx*ny)
      opens = .true.
      orecs = 1

      fname = 'cm1rst_tXXXXXX_u.dat'
      write(fname( 9:14),101) nrst
      if(dowr) write(outfile,*) '  fname=',fname
      open(unit=ufile,file=fname,form='unformatted',access='direct',recl=4*(nx+1)*ny)
      openu = .true.
      orecu = 1

      fname = 'cm1rst_tXXXXXX_v.dat'
      write(fname( 9:14),101) nrst
      if(dowr) write(outfile,*) '  fname=',fname
      open(unit=vfile,file=fname,form='unformatted',access='direct',recl=4*nx*(ny+1))
      openv = .true.
      orecv = 1

      fname = 'cm1rst_tXXXXXX_w.dat'
      write(fname( 9:14),101) nrst
      if(dowr) write(outfile,*) '  fname=',fname
      open(unit=wfile,file=fname,form='unformatted',access='direct',recl=4*nx*ny)
      openw = .true.
      orecw = 1

      if(dowr) write(outfile,*)
    ENDIF

  ELSEIF( restart_filetype.eq.3 )THEN

    !------------------
    ! one restart file per processor:

      do i=1,maxstring
        fname(i:i) = ' '
      enddo

      fname = 'cm1rst_pXXXXXX_tYYYYYY_i.dat'
      write(fname( 9:14),101) myid
      write(fname(17:22),101) nrst

      ! all processors:
      open(unit=infofile,file=fname,form='unformatted',status='unknown')
      openinfo = .true.
      doinfo = .true.

      fname = 'cm1rst_pXXXXXX_tYYYYYY_x.dat'
      write(fname( 9:14),101) myid
      write(fname(17:22),101) nrst

      if(dowr) write(outfile,*)
      if(dowr) write(outfile,*) '  Writing to restart file!'
      if(dowr) write(outfile,*) '  fname=',fname

      open(unit=rstfile,file=fname,form='unformatted',status='unknown')
      openrst  = .true.
      dorstfile = .true.

      fname = 'cm1rst_pXXXXXX_tYYYYYY_s.dat'
      write(fname( 9:14),101) myid
      write(fname(17:22),101) nrst
      if(dowr) write(outfile,*) '  fname=',fname
      open(unit=sfile,file=fname,form='unformatted',access='direct',recl=4*ni*nj)
      opens = .true.
      orecs = 1

      fname = 'cm1rst_pXXXXXX_tYYYYYY_u.dat'
      write(fname( 9:14),101) myid
      write(fname(17:22),101) nrst
      if(dowr) write(outfile,*) '  fname=',fname
      open(unit=ufile,file=fname,form='unformatted',access='direct',recl=4*(ni+1)*nj)
      openu = .true.
      orecu = 1

      fname = 'cm1rst_pXXXXXX_tYYYYYY_v.dat'
      write(fname( 9:14),101) myid
      write(fname(17:22),101) nrst
      if(dowr) write(outfile,*) '  fname=',fname
      open(unit=vfile,file=fname,form='unformatted',access='direct',recl=4*ni*(nj+1))
      openv = .true.
      orecv = 1

      fname = 'cm1rst_pXXXXXX_tYYYYYY_w.dat'
      write(fname( 9:14),101) myid
      write(fname(17:22),101) nrst
      if(dowr) write(outfile,*) '  fname=',fname
      open(unit=wfile,file=fname,form='unformatted',access='direct',recl=4*ni*nj)
      openw = .true.
      orecw = 1

      if(dowr) write(outfile,*)


  ELSEIF( restart_filetype.eq.4 )THEN

    !------------------
    ! one restart file per node:

    IF(myid.eq.nodeleader)THEN

      do i=1,maxstring
        fname(i:i) = ' '
      enddo

      fname = 'cm1rst_nXXXXXX_tYYYYYY_i.dat'

      write(fname( 9:14),101) mynode
      write(fname(17:22),101) nrst

      open(unit=infofile,file=fname,form='unformatted',status='unknown')
      openinfo = .true.
      doinfo = .true.

      fname = 'cm1rst_nXXXXXX_tYYYYYY_x.dat'
      write(fname( 9:14),101) mynode
      write(fname(17:22),101) nrst

      if(dowr) write(outfile,*)
      if(dowr) write(outfile,*) '  Writing to restart file!'
      if(dowr) write(outfile,*) '  fname=',fname

      open(unit=rstfile,file=fname,form='unformatted',status='unknown')
      openrst  = .true.
      dorstfile = .true.

      fname = 'cm1rst_nXXXXXX_tYYYYYY_s.dat'
      write(fname( 9:14),101) mynode
      write(fname(17:22),101) nrst
      if(dowr) write(outfile,*) '  fname=',fname
      open(unit=sfile,file=fname,form='unformatted',access='direct',recl=4*(nimax*ppnode)*njmax)
      opens = .true.
      orecs = 1

      fname = 'cm1rst_nXXXXXX_tYYYYYY_u.dat'
      write(fname( 9:14),101) mynode
      write(fname(17:22),101) nrst
      if(dowr) write(outfile,*) '  fname=',fname
      open(unit=ufile,file=fname,form='unformatted',access='direct',recl=4*((nimax+1)*ppnode)*njmax)
      openu = .true.
      orecu = 1

      fname = 'cm1rst_nXXXXXX_tYYYYYY_v.dat'
      write(fname( 9:14),101) mynode
      write(fname(17:22),101) nrst
      if(dowr) write(outfile,*) '  fname=',fname
      open(unit=vfile,file=fname,form='unformatted',access='direct',recl=4*(nimax*ppnode)*(njmax+1))
      openv = .true.
      orecv = 1

      fname = 'cm1rst_nXXXXXX_tYYYYYY_w.dat'
      write(fname( 9:14),101) mynode
      write(fname(17:22),101) nrst
      if(dowr) write(outfile,*) '  fname=',fname
      open(unit=wfile,file=fname,form='unformatted',access='direct',recl=4*(nimax*ppnode)*njmax)
      openw = .true.
      orecw = 1

      if(dowr) write(outfile,*)

    ENDIF

  ELSE
    stop 12388
  ENDIF

#ifdef NETCDF
  ELSEIF( restart_format.eq.2 )THEN
    ! netcdf format:

    ! broken (at the moment)

    if( myid.eq.0 )then
      print *
      print *,'      restart_format = ',restart_format
      print *
      print *,'  netcdf-format restart files are not available right now '
      print *,'    (but will be fixed someday soon)   '
      print *
    endif
#ifdef MPI
    call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif
    call stopcm1

!    if( myid.eq.0 )then
!      call     restart_prelim(nrst,ncid,mtime,xfref,yfref,zh,zf,sigma,sigmaf,  &
!                              qname,num_soil_layers,nrad2d,dat2(1,1),dat2(1,2),dum1(ib,jb,kb),time_index)
!    endif

#endif
  ELSE

    if( myid.eq.0 )then
      print *
      print *,'  unrecognized value for restart_format '
      print *
      print *,'      restart_format = ',restart_format
      print *
    endif
#ifdef MPI
    call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif
    call stopcm1

  ENDIF

!---------------------------------------------------------------
! metadata:

  IF( restart_format.eq.1 )THEN

    if(myid.eq.0)then
      ! only processor 0 does this:
      write(rstfile) cm1rversion
      write(rstfile) confignum
      write(rstfile) nx
      write(rstfile) ny
      write(rstfile) nz
      write(rstfile) ngxy
      print *,'  nx,ny,nz,ngxy = ',nx,ny,nz,ngxy
      write(rstfile) restart_filetype
      write(rstfile) numprocs
      write(rstfile) ppnode
      write(rstfile) nodes
    endif

    if( dorstfile )then
      write(rstfile) myi1
      write(rstfile) myi2
      write(rstfile) myj1
      write(rstfile) myj2
      if(myid.eq.0) print *,'  myi1,myi2 = ',myi1,myi2
      if(myid.eq.0) print *,'  myj1,myj2 = ',myj1,myj2
      if( restart_filetype.eq.4 )then
        do n=1,ppnode-1
          if(myid.eq.0) print *,'  n = ',n
          write(rstfile) myi1p(myid+n)
          write(rstfile) myi2p(myid+n)
          write(rstfile) myj1p(myid+n)
          write(rstfile) myj2p(myid+n)
          if(myid.eq.0) print *,'  myi1,myi2 = ',myi1p(myid+n),myi2p(myid+n)
          if(myid.eq.0) print *,'  myj1,myj2 = ',myj1p(myid+n),myj2p(myid+n)
        enddo
      endif
    endif

    IF(myid.eq.0)THEN
      ! only processor 0 does this:
      write(rstfile) noutstat
      write(rstfile) nout2ds
      write(rstfile) nout3ds
      write(rstfile) nout3du
      write(rstfile) nout3dv
      write(rstfile) nout3dw
      write(rstfile) nout3dxv
      write(rstfile) nout3dyv
      write(rstfile) nout3dzv
      write(rstfile) noutdiags
      write(rstfile) noutdiagw
      write(rstfile) noutazims
      write(rstfile) noutazimw
      write(rstfile) noutparcel
      write(rstfile) nouthifrq
      write(rstfile) nconfig

      write(rstfile) wbc
      write(rstfile) ebc
      write(rstfile) sbc
      write(rstfile) nbc
      write(rstfile) bbc
      write(rstfile) tbc
      write(rstfile) ptype
      write(rstfile) ipbl
      write(rstfile) sfcmodel
      write(rstfile) sgsmodel










      write(rstfile) (maxx-minx)
      write(rstfile) (maxy-miny)
      write(rstfile) maxz
















      write(rstfile) nstep
      write(rstfile) srec
      write(rstfile) sirec
      write(rstfile) urec
      write(rstfile) vrec
      write(rstfile) wrec
      write(rstfile) nrec
      write(rstfile) mrec
      write(rstfile) prec
      write(rstfile) trecs
      write(rstfile) trecw
      write(rstfile) arecs
      write(rstfile) arecw
      write(rstfile) mvrec
      write(rstfile) nwrite
      write(rstfile) nwritet
      write(rstfile) nwritea
      write(rstfile) nwritemv
      write(rstfile) nwriteh
      write(rstfile) nrst
      write(rstfile) nstatout
      write(rstfile) ndt
      write(rstfile) icenter
      write(rstfile) jcenter
      write(rstfile) output_format
      write(rstfile) output_filetype
      write(rstfile) dt
      write(rstfile) dtlast
      write(rstfile) xcenter
      write(rstfile) ycenter
      write(rstfile) umove
      write(rstfile) vmove
      write(rstfile) umovesfc
      write(rstfile) vmovesfc
      write(rstfile) umovem1
      write(rstfile) vmovem1
      write(rstfile) umovem2
      write(rstfile) vmovem2
      write(rstfile) domainlocx
      write(rstfile) domainlocy
      write(rstfile) adaptmovetim
      write(rstfile) cflmax
      write(rstfile) mtime
      write(rstfile) stattim
      write(rstfile) taptim
      write(rstfile) rsttim
      write(rstfile) radtim
      write(rstfile) prcltim
      write(rstfile) hailtrajtim
      write(rstfile) restart_time
      write(rstfile) spinup_t2p
      write(rstfile) adt
      write(rstfile) adtlast
      write(rstfile) acfl
      write(rstfile) dbldt
      write(rstfile) mass1
      write(rstfile) avgsfcu
      write(rstfile) avgsfcv
      write(rstfile) avgsfcs
      write(rstfile) avgsfcsu
      write(rstfile) avgsfcsv
      write(rstfile) avgsfct
      write(rstfile) avgsfcq
      write(rstfile) avgsfcp
      write(rstfile) xhref
      write(rstfile) xfref
      write(rstfile) yhref
      write(rstfile) yfref
      write(rstfile) sigma
      write(rstfile) sigmaf
      print *
      print *,'  Depth of domain: '
      print *,'  sigma = ',sigma(1),sigma(nz)
      print *,'  sigmaf = ',sigmaf(1),sigmaf(nz+1)
      print *
    ENDIF
#ifdef NETCDF
  ELSEIF( restart_format.eq.2 )THEN

    call disp_err( nf90_inq_varid(ncid,"cm1rversion",varid) , .true. )
    call disp_err( nf90_put_var(ncid,varid,cm1rversion,(/time_index/)) , .true. )

    IF(myid.eq.0)THEN

      call disp_err( nf90_inq_varid(ncid,"nstep",varid) , .true. )
      call disp_err( nf90_put_var(ncid,varid,nstep,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"srec",varid) , .true. )
      call disp_err( nf90_put_var(ncid,varid,srec,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"sirec",varid) , .true. )
      call disp_err( nf90_put_var(ncid,varid,sirec,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"urec",varid) , .true. )
      call disp_err( nf90_put_var(ncid,varid,urec,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"vrec",varid) , .true. )
      call disp_err( nf90_put_var(ncid,varid,vrec,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"wrec",varid) , .true. )
      call disp_err( nf90_put_var(ncid,varid,wrec,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"nrec",varid) , .true. )
      call disp_err( nf90_put_var(ncid,varid,nrec,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"mrec",varid) , .true. )
      call disp_err( nf90_put_var(ncid,varid,mrec,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"prec",varid) , .true. )
      call disp_err( nf90_put_var(ncid,varid,prec,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"trecs",varid) , .true. )
      call disp_err( nf90_put_var(ncid,varid,trecs,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"trecw",varid) , .true. )
      call disp_err( nf90_put_var(ncid,varid,trecw,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"arecs",varid) , .true. )
      call disp_err( nf90_put_var(ncid,varid,arecs,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"arecw",varid) , .true. )
      call disp_err( nf90_put_var(ncid,varid,arecw,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"mvrec",varid) , .true. )
      call disp_err( nf90_put_var(ncid,varid,mvrec,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"nwrite",varid) , .true. )
      call disp_err( nf90_put_var(ncid,varid,nwrite,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"nwritet",varid) , .true. )
      call disp_err( nf90_put_var(ncid,varid,nwritet,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"nwritea",varid) , .true. )
      call disp_err( nf90_put_var(ncid,varid,nwritea,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"nwritemv",varid) , .true. )
      call disp_err( nf90_put_var(ncid,varid,nwritemv,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"nwriteh",varid) , .true. )
      call disp_err( nf90_put_var(ncid,varid,nwriteh,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"nrst",varid) , .true. )
      call disp_err( nf90_put_var(ncid,varid,nrst,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"nstatout",varid) , .true. )
      call disp_err( nf90_put_var(ncid,varid,nstatout,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"ndt",varid) , .true. )
      call disp_err( nf90_put_var(ncid,varid,ndt,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"icenter",varid) , .true. )
      call disp_err( nf90_put_var(ncid,varid,icenter,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"jcenter",varid) , .true. )
      call disp_err( nf90_put_var(ncid,varid,jcenter,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"old_format",varid) , .true. )
      call disp_err( nf90_put_var(ncid,varid,output_format,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"dt",varid) , .true. )
      call disp_err( nf90_put_var(ncid,varid,dt,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"dtlast",varid) , .true. )
      call disp_err( nf90_put_var(ncid,varid,dtlast,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"xcenter",varid) , .true. )
      call disp_err( nf90_put_var(ncid,varid,xcenter,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"ycenter",varid) , .true. )
      call disp_err( nf90_put_var(ncid,varid,ycenter,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"umove",varid) , .true. )
      call disp_err( nf90_put_var(ncid,varid,umove,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"vmove",varid) , .true. )
      call disp_err( nf90_put_var(ncid,varid,vmove,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"domainlocx",varid) , .true. )
      call disp_err( nf90_put_var(ncid,varid,domainlocx,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"domainlocy",varid) , .true. )
      call disp_err( nf90_put_var(ncid,varid,domainlocy,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"adaptmovetim",varid) , .true. )
      call disp_err( nf90_put_var(ncid,varid,adaptmovetim,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"cflmax",varid) , .true. )
      call disp_err( nf90_put_var(ncid,varid,cflmax,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"mtime",varid) , .true. )
      call disp_err( nf90_put_var(ncid,varid,mtime,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"stattim",varid) , .true. )
      call disp_err( nf90_put_var(ncid,varid,stattim,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"taptim",varid) , .true. )
      call disp_err( nf90_put_var(ncid,varid,taptim,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"rsttim",varid) , .true. )
      call disp_err( nf90_put_var(ncid,varid,rsttim,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"radtim",varid) , .true. )
      call disp_err( nf90_put_var(ncid,varid,radtim,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"prcltim",varid) , .true. )
      call disp_err( nf90_put_var(ncid,varid,prcltim,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"hailtrajtim",varid) , .true. )
      call disp_err( nf90_put_var(ncid,varid,hailtrajtim,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"adt",varid) , .true. )
      call disp_err( nf90_put_var(ncid,varid,adt,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"adtlast",varid) , .true. )
      call disp_err( nf90_put_var(ncid,varid,adtlast,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"acfl",varid) , .true. )
      call disp_err( nf90_put_var(ncid,varid,acfl,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"dbldt",varid) , .true. )
      call disp_err( nf90_put_var(ncid,varid,dbldt,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"mass1",varid) , .true. )
      call disp_err( nf90_put_var(ncid,varid,mass1,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"avgsfcu",varid) , .true. )
      call disp_err( nf90_put_var(ncid,varid,avgsfcu,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"avgsfcv",varid) , .true. )
      call disp_err( nf90_put_var(ncid,varid,avgsfcv,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"avgsfcs",varid) , .true. )
      call disp_err( nf90_put_var(ncid,varid,avgsfcs,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"avgsfcsu",varid) , .true. )
      call disp_err( nf90_put_var(ncid,varid,avgsfcsu,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"avgsfcsv",varid) , .true. )
      call disp_err( nf90_put_var(ncid,varid,avgsfcsv,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"avgsfct",varid) , .true. )
      call disp_err( nf90_put_var(ncid,varid,avgsfct,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"avgsfcq",varid) , .true. )
      call disp_err( nf90_put_var(ncid,varid,avgsfcq,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"avgsfcp",varid) , .true. )
      call disp_err( nf90_put_var(ncid,varid,avgsfcp,(/time_index/)) , .true. )

    ENDIF
#endif
  ENDIF

!-----------------------------------------------------------------------
#ifdef MPI
      cfoo = 0.0
      call MPI_REDUCE(qbudget(1),cfoo(1),nbudget,MPI_DOUBLE_PRECISION,MPI_SUM,0,  &
                      MPI_COMM_WORLD,ierr)
      if( myid.eq.0 )then
        do n=1,nbudget
          qbudget(n)=cfoo(n)
        enddo
      else
        qbudget = 0.0
      endif
      if( imoist.eq.1 )then
        afoo = 0.0
        call MPI_REDUCE(asq(1),afoo(1),numq,MPI_DOUBLE_PRECISION,MPI_SUM,0,  &
                        MPI_COMM_WORLD,ierr)
        if( myid.eq.0 )then
          do n=1,numq
            asq(n)=afoo(n)
          enddo
        else
          asq = 0.0
        endif
        bfoo = 0.0
        call MPI_REDUCE(bsq(1),bfoo(1),numq,MPI_DOUBLE_PRECISION,MPI_SUM,0,  &
                        MPI_COMM_WORLD,ierr)
        if( myid.eq.0 )then
          do n=1,numq
            bsq(n)=bfoo(n)
          enddo
        else
          bsq = 0.0
        endif
      endif
#endif
!-----------------------------------------------------------------------
! budget variables:

    IF( myid.eq.0 )THEN

      IF( restart_format.eq.1 )THEN
        write(rstfile) nbudget
        write(rstfile) numq
        print *,'  nbudget,numq = ',nbudget,numq
        write(rstfile) qbudget
        write(rstfile) asq
        write(rstfile) bsq
#ifdef NETCDF
      ELSEIF( restart_format.eq.2 )THEN
        call disp_err( nf90_inq_varid(ncid,"qbudget",varid) , .true. )
        call disp_err( nf90_put_var(ncid,varid,qbudget,(/1,time_index/),(/nbudget,1/)) , .true. )
        call disp_err( nf90_inq_varid(ncid,"asq",varid) , .true. )
        call disp_err( nf90_put_var(ncid,varid,asq,(/1,time_index/),(/numq,1/)) , .true. )
        call disp_err( nf90_inq_varid(ncid,"bsq",varid) , .true. )
        call disp_err( nf90_put_var(ncid,varid,bsq,(/1,time_index/),(/numq,1/)) , .true. )
#endif
      ENDIF

    ENDIF

!---------------------------------------------------------------

  num = 0

  IF( do_lsnudge .or. do_adapt_move )THEN

    in_file = in_file_yes

    if( myid.eq.0 )then
      write(rstfile) ug
      write(rstfile) vg
    endif

  ELSE

    in_file = in_file_no

  ENDIF

  if(doinfo) call writeinfo(in_file,orecs,orecu,orecv,orecw,'  ug,vg:')

!---------------------------------------------------------------
! standard 2D:

      call writer(ni,nj,1,1,nx,ny,tsk(ib,jb),'tsk     ',           &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)

      n = 1
      call writer(ni,nj,1,1,nx,ny,rain(ib,jb,n),'rain    ',        &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call writer(ni,nj,1,1,nx,ny,raincv(ib,jb,n),'raincv  ',        &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call writer(ni,nj,1,1,nx,ny,hail(ib,jb,n),'hail    ',        &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call writer(ni,nj,1,1,nx,ny,sws(ib,jb,n),'sws     ',         &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call writer(ni,nj,1,1,nx,ny,svs(ib,jb,n),'svs     ',         &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call writer(ni,nj,1,1,nx,ny,sps(ib,jb,n),'sps     ',         &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call writer(ni,nj,1,1,nx,ny,srs(ib,jb,n),'srs     ',         &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call writer(ni,nj,1,1,nx,ny,sgs(ib,jb,n),'sgs     ',         &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call writer(ni,nj,1,1,nx,ny,sus(ib,jb,n),'sus     ',         &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call writer(ni,nj,1,1,nx,ny,shs(ib,jb,n),'shs     ',         &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)

!---------------------------------------------------------------
!  for "moving" swaths:

    if( nrain.eq.2 )then
      in_file = in_file_yes
      n = 2
      call writer(ni,nj,1,1,nx,ny,rain(ib,jb,n),'rain2   ',        &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call writer(ni,nj,1,1,nx,ny,raincv(ib,jb,n),'raincv2 ',        &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call writer(ni,nj,1,1,nx,ny,hail(ib,jb,n),'hail2   ',        &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call writer(ni,nj,1,1,nx,ny,sws(ib,jb,n),'sws2    ',         &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call writer(ni,nj,1,1,nx,ny,svs(ib,jb,n),'svs2    ',         &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call writer(ni,nj,1,1,nx,ny,sps(ib,jb,n),'sps2    ',         &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call writer(ni,nj,1,1,nx,ny,srs(ib,jb,n),'srs2    ',         &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call writer(ni,nj,1,1,nx,ny,sgs(ib,jb,n),'sgs2    ',         &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call writer(ni,nj,1,1,nx,ny,sus(ib,jb,n),'sus2    ',         &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call writer(ni,nj,1,1,nx,ny,shs(ib,jb,n),'shs2    ',         &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
    else
      in_file = in_file_no
    endif

  if(doinfo) call writeinfo(in_file,orecs,orecu,orecv,orecw,'  nrain2:')

!---------------------------------------------------------------
! standard 3D:

      call writer(ni,nj,1,nk,nx,ny,rho(ib,jb,1),'rho     ',        &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call writer(ni,nj,1,nk,nx,ny,prs(ib,jb,1),'prs     ',        &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call writer(ni+1,nj,1,nk,nx+1,ny,ua(ib,jb,1),'ua      ',     &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecu,ufile,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2iu,d2ju,d3iu,d3ju)
      call writer(ni,nj+1,1,nk,nx,ny+1,va(ib,jb,1),'va      ',     &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecv,vfile,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2iv,d2jv,d3iv,d3jv)
      call writer(ni,nj,1,nk+1,nx,ny,wa(ib,jb,1),'wa      ',       &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecw,wfile,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call writer(ni,nj,1,nk,nx,ny,ppi(ib,jb,1),'ppi     ',        &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      do k=1,nk
      do j=1,nj
      do i=1,ni
        dum1(i,j,k) = pi0(i,j,k)+ppi(i,j,k)
      enddo
      enddo
      enddo
      call writer(ni,nj,1,nk,nx,ny,dum1(ib,jb,1),'pptot   ',        &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call writer(ni,nj,1,nk,nx,ny,tha(ib,jb,1),'tha     ',        &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      do k=1,nk
      do j=1,nj
      do i=1,ni
        dum1(i,j,k) = th0(i,j,k)+tha(i,j,k)
      enddo
      enddo
      enddo
      call writer(ni,nj,1,nk,nx,ny,dum1(ib,jb,1),'thtot   ',        &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call writer(ni,nj,1,nk,nx,ny,ppx(ib,jb,1),'ppx     ',        &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)

!---------------------------------------------------------------

    if( psolver.eq.6 )then
      in_file = in_file_yes
      call writer(ni,nj,1,nk,nx,ny,phi1(ib,jb,1),'phi1    ',       &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call writer(ni,nj,1,nk,nx,ny,phi2(ib,jb,1),'phi2    ',       &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
    else
      in_file = in_file_no
    endif

    if(doinfo) call writeinfo(in_file,orecs,orecu,orecv,orecw,'  phi:')

!---------------------------------------------------------------
!  moisture:

    IF(imoist.eq.1)THEN
      in_file = in_file_yes
    do n=1,numq
      text1 = '        '
      write(text1(1:3),156) qname(n)
156   format(a3)
      call writer(ni,nj,1,nk,nx,ny,qa(ib,jb,1,n),text1     ,       &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
    enddo
    ELSE
      in_file = in_file_no
    ENDIF

    if(doinfo) call writeinfo(in_file,orecs,orecu,orecv,orecw,'  qa:',numq)

!---------------------------------------------------------------
!  vars for eqtset=2

    if(imoist.eq.1.and.eqtset.eq.2)then
      in_file = in_file_yes
      call writer(ni,nj,1,nk,nx,ny,qpten(ib,jb,1),'qpten   ',      &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call writer(ni,nj,1,nk,nx,ny,qtten(ib,jb,1),'qtten   ',      &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call writer(ni,nj,1,nk,nx,ny,qvten(ib,jb,1),'qvten   ',      &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call writer(ni,nj,1,nk,nx,ny,qcten(ib,jb,1),'qcten   ',      &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
    else
      in_file = in_file_no
    endif

    if(doinfo) call writeinfo(in_file,orecs,orecu,orecv,orecw,'  qxten:')

!---------------------------------------------------------------
!  subgrid tke for LES:

    if( idoles .and. iusetke )then
      in_file = in_file_yes
      call writer(ni,nj,1,nk+1,nx,ny,tkea(ib,jb,1),'tkea    ',     &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecw,wfile,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
    else
      in_file = in_file_no
    endif

    if(doinfo) call writeinfo(in_file,orecs,orecu,orecv,orecw,'  tkea:')

!---------------------------------------------------------------
!  radiation:

      radoptge1:  &
      if(radopt.ge.1)then
        in_file = in_file_yes
        call writer(ni,nj,1,nk,nx,ny,lwten(ib,jb,1),'lwten   ',      &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,nk,nx,ny,swten(ib,jb,1),'swten   ',      &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,nk,nx,ny,thraten(ib,jb,1),'thraten ',    &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
!$omp parallel do default(shared)  &
!$omp private(i,j)
        do j=1,nj
        do i=1,ni
          dum1(i,j,1) = radsw(i,j)
        enddo
        enddo
        call writer(ni,nj,1,1,nx,ny,dum1(ib,jb,1),'radsw   ',        &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
!$omp parallel do default(shared)  &
!$omp private(i,j)
        do j=1,nj
        do i=1,ni
          dum1(i,j,1) = rnflx(i,j)
        enddo
        enddo
        call writer(ni,nj,1,1,nx,ny,dum1(ib,jb,1),'rnflx   ',        &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
!$omp parallel do default(shared)  &
!$omp private(i,j)
        do j=1,nj
        do i=1,ni
          dum1(i,j,1) = radswnet(i,j)
        enddo
        enddo
        call writer(ni,nj,1,1,nx,ny,dum1(ib,jb,1),'radswnet',        &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
!$omp parallel do default(shared)  &
!$omp private(i,j)
        do j=1,nj
        do i=1,ni
          dum1(i,j,1) = radlwin(i,j)
        enddo
        enddo
        call writer(ni,nj,1,1,nx,ny,dum1(ib,jb,1),'radlwin ',        &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        do n=1,nrad2d
!$omp parallel do default(shared)  &
!$omp private(i,j)
        do j=1,nj
        do i=1,ni
          dum1(i,j,1) = rad2d(i,j,n)
        enddo
        enddo
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
        call writer(ni,nj,1,1,nx,ny,dum1(ib,jb,1),text1,             &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        enddo
      else
        in_file = in_file_no
      endif  radoptge1

      if(doinfo) call writeinfo(in_file,orecs,orecu,orecv,orecw,'  radiation:')

!---------------------------------------------------------------

      if( radopt.ge.1 .and. ptype.eq.5 )then
        in_file = in_file_yes
        call writer(ni,nj,1,nk,nx,ny,effc(ib,jb,1),'effc    ',       &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,nk,nx,ny,effi(ib,jb,1),'effi    ',       &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,nk,nx,ny,effs(ib,jb,1),'effs    ',       &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,nk,nx,ny,effr(ib,jb,1),'effr    ',       &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,nk,nx,ny,effg(ib,jb,1),'effg    ',       &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,nk,nx,ny,effis(ib,jb,1),'effis   ',      &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      else
        in_file = in_file_no
      endif

    if(doinfo) call writeinfo(in_file,orecs,orecu,orecv,orecw,'  effx:')

!---------------------------------------------------------------
!  cu param:

      if(cuparam.ge.1)then
        in_file = in_file_yes
        call writer(ni,nj,1,nk,nx,ny,thcuten(ib,jb,1),'thcuten ',      &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,nk,nx,ny,qvcuten(ib,jb,1),'qvcuten ',      &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,nk,nx,ny,qccuten(ib,jb,1),'qccuten ',      &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,nk,nx,ny,qicuten(ib,jb,1),'qicuten ',      &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,nk,nx,ny,ucuten(ib,jb,1),'ucuten  ',      &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,nk,nx,ny,vcuten(ib,jb,1),'vcuten  ',      &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      else
        in_file = in_file_no
      endif

    if(doinfo) call writeinfo(in_file,orecs,orecu,orecv,orecw,'  cuparam:')

!---------------------------------------------------------------
!  surface:
!     I don't know how many of these are really needed in restart
!     files, but let's include them all for now ... just to be safe

      sfcvars1:  &
      if( ((oceanmodel.eq.2).or.(ipbl.ge.1).or.(sfcmodel.ge.1)) .and. sfcmodel.ge.1 )then
        in_file = in_file_yes
        call writer(ni,nj,1,1,nx,ny,ust(ib,jb),'ust     ',           &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,znt(ib,jb),'znt     ',           &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,cd(ib,jb),'cd      ',            &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,ch(ib,jb),'ch      ',            &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,cq(ib,jb),'cq      ',            &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,u1(ib,jb),'u1      ',            &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,v1(ib,jb),'v1      ',            &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,s1(ib,jb),'s1      ',            &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,th1(ib,jb),'th1     ',           &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,qv1(ib,jb),'qv1     ',           &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,u10(ib,jb),'u10     ',           &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,v10(ib,jb),'v10     ',           &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,s10(ib,jb),'s10     ',           &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,xland(ib,jb),'xland   ',         &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,thflux(ib,jb),'thflux  ',        &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,qvflux(ib,jb),'qvflux  ',        &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,psfc(ib,jb),'psfc    ',          &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      else
        in_file = in_file_no
      endif  sfcvars1

    if(doinfo) call writeinfo(in_file,orecs,orecu,orecv,orecw,'  sfcvars1:')

  !.............................................................................

      sfcvars2:  &
      if( ((oceanmodel.eq.2).or.(ipbl.ge.1).or.(sfcmodel.ge.1)) .and. tbc.eq.3 )then
        in_file = in_file_yes
        call writer(ni,nj,1,1,nx,ny,ustt(ib,jb),'ustt    ',          &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,  ut(ib,jb),'ut      ',          &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,  vt(ib,jb),'vt      ',          &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,  st(ib,jb),'st      ',          &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      else
        in_file = in_file_no
      endif  sfcvars2

      if(doinfo) call writeinfo(in_file,orecs,orecu,orecv,orecw,'  sfcvars2:')

  !.............................................................................

      sfcvars3:  &
      if( ((oceanmodel.eq.2).or.(ipbl.ge.1).or.(sfcmodel.ge.1)) .and. sfcmodel.ge.1 )then
        in_file = in_file_yes

!$omp parallel do default(shared)  &
!$omp private(i,j)
        do j=1,nj
        do i=1,ni
          dum1(i,j,1) = lu_index(i,j)
        enddo
        enddo
        call writer(ni,nj,1,1,nx,ny,dum1(ib,jb,1),'lu_index',        &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
!$omp parallel do default(shared)  &
!$omp private(i,j)
        do j=1,nj
        do i=1,ni
          dum1(i,j,1) = kpbl2d(i,j)
        enddo
        enddo
        call writer(ni,nj,1,1,nx,ny,dum1(ib,jb,1),'kpbl2d  ',        &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,hfx(ib,jb),'hfx     ',           &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,qfx(ib,jb),'qfx     ',           &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,hpbl(ib,jb),'hpbl    ',          &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,wspd(ib,jb),'wspd    ',          &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,psim(ib,jb),'psim    ',          &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,psih(ib,jb),'psih    ',          &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,gz1oz0(ib,jb),'gz1oz0  ',        &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,br(ib,jb),'br      ',            &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,CHS(ib,jb),'chs     ',           &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,CHS2(ib,jb),'chs2    ',          &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,CQS2(ib,jb),'cqs2    ',          &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,CPMM(ib,jb),'cpmm    ',          &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,ZOL(ib,jb),'zol     ',           &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,MAVAIL(ib,jb),'mavail  ',        &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,MOL(ib,jb),'mol     ',           &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,RMOL(ib,jb),'rmol    ',          &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,REGIME(ib,jb),'regime  ',        &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,LH(ib,jb),'lh      ',            &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,tmn(ib,jb),'tmn     ',           &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,FLHC(ib,jb),'flhc    ',          &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,FLQC(ib,jb),'flqc    ',          &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,QGH(ib,jb),'qgh     ',           &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,CK(ib,jb),'ck      ',            &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,CKA(ib,jb),'cka     ',           &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,CDA(ib,jb),'cda     ',           &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,USTM(ib,jb),'ustm    ',          &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,QSFC(ib,jb),'qsfc    ',          &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,T2(ib,jb),'t2      ',            &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,Q2(ib,jb),'q2      ',            &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,TH2(ib,jb),'th2     ',           &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,EMISS(ib,jb),'emiss   ',         &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,THC(ib,jb),'thc     ',           &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,ALBD(ib,jb),'albd    ',          &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,gsw(ib,jb),'gsw     ',           &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,glw(ib,jb),'glw     ',           &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,chklowq(ib,jb),'chklowq ',       &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,capg(ib,jb),'capg    ',          &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,snowc(ib,jb),'snowc   ',         &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,fm(ib,jb),'fm      ',            &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,fh(ib,jb),'fh      ',            &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,mznt(ib,jb),'mznt    ',          &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,swspd(ib,jb),'swspd   ',         &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,wstar(ib,jb),'wstar   ',         &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,delta(ib,jb),'delta   ',         &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
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
          call writer(ni,nj,1,1,nx,ny,tslb(ib,jb,n),text1,             &
                      ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                      ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                      dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        enddo
      else
        in_file = in_file_no
      endif  sfcvars3

      if(doinfo) call writeinfo(in_file,orecs,orecu,orecv,orecw,'  sfcvars3:')

  !.............................................................................

      omlvars:  &
      if(oceanmodel.eq.2)then
        in_file = in_file_yes
        call writer(ni,nj,1,1,nx,ny,tml(ib,jb),'tml     ',           &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,t0ml(ib,jb),'t0ml    ',          &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,hml(ib,jb),'hml     ',           &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,h0ml(ib,jb),'h0ml    ',          &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,huml(ib,jb),'huml    ',          &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,hvml(ib,jb),'hvml    ',          &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,tmoml(ib,jb),'tmoml   ',         &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      else
        in_file = in_file_no
      endif  omlvars

      if(doinfo) call writeinfo(in_file,orecs,orecu,orecv,orecw,'  omlvars:')

!---------------------------------------------------------------
!  MYNN and MYJ:

    myxvars:  &
    IF( ipbl.eq.4 .or. ipbl.eq.5 .or. ipbl.eq.6 )THEN
      if(myid.eq.0) print *,'  writing myx-3d vars ... '
      in_file = in_file_yes
      do n=1,npbl3d
        call writer(ni,nj,1,nk,nx,ny,pbl3d(ib,jb,1,n),'myx-3d  ',        &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      enddo
    ELSE
      in_file = in_file_no
    ENDIF  myxvars
    if(doinfo) call writeinfo(in_file,orecs,orecu,orecv,orecw,'  myx-3d:',npbl3d)


    IF( sfcmodel.eq.6 .or. sfcmodel.eq.7 )THEN
      if(myid.eq.0) print *,'  writing myx-2d vars ... '
      in_file = in_file_yes
      do n=1,npbl2d
        call writer(ni,nj,1,1,nx,ny,pbl2d(ib,jb,n),'myx-2d  ',        &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      enddo
    ELSE
      in_file = in_file_no
    ENDIF
    if(doinfo) call writeinfo(in_file,orecs,orecu,orecv,orecw,'  myx-2d:',npbl2d)


    IF( sfcmodel.eq.6 .or. sfcmodel.eq.7 )THEN
      if(myid.eq.0) print *,'  writing myx-2d vars ... '
      in_file = in_file_yes
      do n=1,npbl2dint
        do j=1,nj
        do i=1,ni
          dum1(i,j,1) = pbl2dint(i,j,n)
        enddo
        enddo
        call writer(ni,nj,1,1,nx,ny,dum1(ib,jb,1),'myx-2di ',        &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      enddo

      if(myid.eq.0) print *,'  ... done writing myx vars '
    ELSE
      in_file = in_file_no
    ENDIF
    if(doinfo) call writeinfo(in_file,orecs,orecu,orecv,orecw,'  myx:',npbl2dint)

!---------------------------------------------------------------
!  PBL tendencies:

    pbltends:  &
    IF( ipbl.eq.4 .or. ipbl.eq.5 .or. ipbl.eq.6 )THEN
      if(myid.eq.0) print *,'  writing pbl tends ... '
      in_file = in_file_yes
      call writer(ni,nj,1,nk,nx,ny,thpten(ib,jb,1),'thpten  ',     &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call writer(ni,nj,1,nk,nx,ny,qvpten(ib,jb,1),'qvpten  ',     &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call writer(ni,nj,1,nk,nx,ny,qcpten(ib,jb,1),'qcpten  ',     &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call writer(ni,nj,1,nk,nx,ny,qipten(ib,jb,1),'qipten  ',     &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call writer(ni,nj,1,nk,nx,ny,upten(ib,jb,1),'upten   ',      &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call writer(ni,nj,1,nk,nx,ny,vpten(ib,jb,1),'vpten   ',      &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call writer(ni,nj,1,nk,nx,ny,qnipten(ib,jb,1),'qnipten ',    &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call writer(ni,nj,1,nk,nx,ny,qncpten(ib,jb,1),'qncpten ',    &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
    ELSE
      in_file = in_file_no
    ENDIF  pbltends

    if(doinfo) call writeinfo(in_file,orecs,orecu,orecv,orecw,'  pbl tends:')

!---------------------------------------------------------------
!  2-part turbulence param:

      if( (sgsmodel.eq.3.or.sgsmodel.eq.4) .and. t2p_avg.eq.1 )then
        in_file = in_file_yes
        if(myid.eq.0)then
          write(rstfile) gamk
          write(rstfile) kmw
          write(rstfile) ufw
          write(rstfile) vfw
          write(rstfile) u1b
          write(rstfile) v1b
        endif
      else
        in_file = in_file_no
      endif

      if(doinfo) call writeinfo(in_file,orecs,orecu,orecv,orecw,'  t2p_avg1:')

!---------------------------------------------------------------
!  passive tracers:

      ptvars:  &
      if(iptra.eq.1)then
        in_file = in_file_yes

        if(myid.eq.0)then
          if( restart_format.eq.1 )then
#ifdef NETCDF
          elseif( restart_format.eq.2 )then
            call disp_err( nf90_inq_varid(ncid,"npt",varid) , .true. )
            call disp_err( nf90_put_var(ncid,varid,npt,(/time_index/)) , .true. )
#endif
          endif
        endif
        do n=1,npt
          if( n.lt.10 )then
            text1 = 'ptX     '
            write(text1(3:3),161) n
161         format(i1.1)
          elseif( n.lt.100 )then
            text1 = 'ptXX    '
            write(text1(3:4),162) n
162         format(i2.2)
          else
            stop 11512
          endif
          call writer(ni,nj,1,nk,nx,ny,pta(ib,jb,1,n),text1,           &
                      ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                      ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                      dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        enddo

      else
        in_file = in_file_no
      endif  ptvars

      if(doinfo) call writeinfo(in_file,orecs,orecu,orecv,orecw,'  pt:',npt)

!---------------------------------------------------------------
!  time-average vars:

      IF( dotimeavg )THEN
        in_file = in_file_yes

        if( myid.eq.0 )then
          write(rstfile) ntim
          write(rstfile) ntavr
          write(rstfile) nsfctavr
          write(rstfile) keta
          write(rstfile) utav
          write(rstfile) vtav
          write(rstfile) wtav
          write(rstfile) ttav
          write(rstfile) qtav
          write(rstfile) etav
          write(rstfile) uutav
          write(rstfile) vvtav
          write(rstfile) ntavg
          write(rstfile) rtavg
        endif

        !----------
        !  3d timeavg vars:
        nloopa:  &
        do n=1,ntavr
          if(     n.eq.utav .or. n.eq.uutav )then
            do k=kbta,keta
            do j=1,nj
            do i=1,ni+1
              dumu(i,j,k) = timavg(i,j,k,n)
            enddo
            enddo
            enddo
            call writer(ni+1,nj,kbta,keta,nx+1,ny,dumu(ib,jb,1),'timavgu ',     &
                        ni,nj,ngxy,myid,numprocs,nodex,nodey,orecu,ufile,   &
                        ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                        dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2iu,d2ju,d3iu,d3ju)
          elseif( n.eq.vtav .or. n.eq.vvtav )then
            do k=kbta,keta
            do j=1,nj+1
            do i=1,ni
              dumv(i,j,k) = timavg(i,j,k,n)
            enddo
            enddo
            enddo
            call writer(ni,nj+1,kbta,keta,nx,ny+1,dumv(ib,jb,1),'timavgv ',     &
                        ni,nj,ngxy,myid,numprocs,nodex,nodey,orecv,vfile,   &
                        ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                        dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2iv,d2jv,d3iv,d3jv)
          else
            call writer(ni,nj,kbta,keta,nx,ny,timavg(ibta,jbta,kbta,n),'timavg  ',               &
                        ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,                           &
                        ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                        dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
          endif
        enddo  nloopa
        !----------
        !  sfc timeavg vars:
        nloopb:  &
        do n=1,nsfctavr
          call writer(ni,nj,1,1,nx,ny,sfctimavg(ibta,jbta,n),'sfctimav',                       &
                      ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,                           &
                      ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                      dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        enddo  nloopb
        !----------
        !  3d tavg arrays:
        nloopc:  &
        do n=1,ntavr
        do nt=1,ntim
          if(     n.eq.utav .or. n.eq.uutav )then
            do k=kbta,keta
            do j=1,nj
            do i=1,ni+1
              dumu(i,j,k) = tavg(i,j,k,nt,n)
            enddo
            enddo
            enddo
            call writer(ni+1,nj,kbta,keta,nx+1,ny,dumu(ib,jb,1),'tavgu   ',     &
                        ni,nj,ngxy,myid,numprocs,nodex,nodey,orecu,ufile,   &
                        ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                        dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2iu,d2ju,d3iu,d3ju)
          elseif( n.eq.vtav .or. n.eq.vvtav )then
            do k=kbta,keta
            do j=1,nj+1
            do i=1,ni
              dumv(i,j,k) = tavg(i,j,k,nt,n)
            enddo
            enddo
            enddo
            call writer(ni,nj+1,kbta,keta,nx,ny+1,dumv(ib,jb,1),'tavgv   ',     &
                        ni,nj,ngxy,myid,numprocs,nodex,nodey,orecv,vfile,   &
                        ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                        dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2iv,d2jv,d3iv,d3jv)
          else
            do k=kbta,keta
            do j=1,nj
            do i=1,ni
              dum1(i,j,k) = tavg(i,j,k,nt,n)
            enddo
            enddo
            enddo
            call writer(ni,nj,kbta,keta,nx,ny,dum1(ib,jb,kbta),'tavg    ',                       &
                        ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,                           &
                        ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                        dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
          endif
        enddo
        enddo  nloopc
        !----------
        !  sfc tavg arrays:
        nloopd:  &
        do n=1,nsfctavr
        do nt=1,ntim
          do j=1,nj
          do i=1,ni
            dum1(i,j,1) = sfctavg(i,j,nt,n)
          enddo
          enddo
          call writer(ni,nj,1,1,nx,ny,dum1(ib,jb,1),'sfctavg ',                                &
                      ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,                           &
                      ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                      dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        enddo
        enddo  nloopd
        !----------

      ELSE
        in_file = in_file_no
      ENDIF

      if(doinfo) call writeinfo(in_file,orecs,orecu,orecv,orecw,'  timavg:')

!---------------------------------------------------------------
!  parcels:

      if(iprcl.eq.1)then
        in_file = in_file_yes
        !-----------------
        ! with parcels:
        if(myid.eq.0)then
          if( restart_format.eq.1 )then
            write(rstfile) nparcels
#ifdef NETCDF
          elseif( restart_format.eq.2 )then
            call disp_err( nf90_inq_varid(ncid,"numparcels",varid) , .true. )
            call disp_err( nf90_put_var(ncid,varid,nparcels,(/time_index/)) , .true. )
#endif
          endif
        endif
        ! only write position info:
        if(myid.eq.0)then
          if( .not. terrain_flag )then
            DO np=1,nparcels
              ploc(np,1)=pdata(np,prx)
              ploc(np,2)=pdata(np,pry)
              ploc(np,3)=pdata(np,prz)
            ENDDO
          else
            DO np=1,nparcels
              ploc(np,1)=pdata(np,prx)
              ploc(np,2)=pdata(np,pry)
              ploc(np,3)=pdata(np,prsig)
            ENDDO
          endif
          if( restart_format.eq.1 )then
            write(rstfile) ploc
#ifdef NETCDF
          elseif( restart_format.eq.2 )then
            call disp_err( nf90_inq_varid(ncid,"ploc",varid) , .true. )
            n = 3
            call disp_err( nf90_put_var(ncid,varid,ploc,(/1,1,time_index/),(/nparcels,n,1/)) , .true. )
#endif
          endif
        endif
      else
        in_file = in_file_no
      endif

      if(doinfo) call writeinfo(in_file,orecs,orecu,orecv,orecw,'  parcels:')

!---------------------------------------------------------------
!  hailstones:

      if(ihailtraj.eq.1)then
        in_file = in_file_yes
        !-----------------
        ! with hailstones:
        if(myid.eq.0)then
          if( restart_format.eq.1 )then
            write(rstfile) nhailtrajs
#ifdef NETCDF
          elseif( restart_format.eq.2 )then
            call disp_err( nf90_inq_varid(ncid,"numhailtrajs",varid) , .true. )
            call disp_err( nf90_put_var(ncid,varid,nhailtrajs,(/time_index/)) , .true. )
#endif
          endif
        endif
        ! need position info as well as 10 additional internal characteristics:
        if(myid.eq.0)then
          DO np=1,nhailtrajs
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
          ENDDO
          if( restart_format.eq.1 )then
            write(rstfile) hailloc
#ifdef NETCDF
          elseif( restart_format.eq.2 )then
            call disp_err( nf90_inq_varid(ncid,"hailloc",varid) , .true. )
            n = 13
            call disp_err( nf90_put_var(ncid,varid,hailloc,(/1,1,time_index/),(/nhailtrajs,n,1/)) , .true. )
#endif
          endif
        endif
      else
        in_file = in_file_no
      endif

      if(doinfo) call writeinfo(in_file,orecs,orecu,orecv,orecw,'  hailstones:')

!---------------------------------------------------------------
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!---------------------------------------------------------------
!  open bc:

      irbc4:  &
      if(wbc.eq.2.or.ebc.eq.2.or.sbc.eq.2.or.nbc.eq.2)then
        in_file = in_file_yes

        !----------------------
        !cccccccccccccccccccccc
        !----------------------
        if(myid.eq.0)then
          ndum = ny
        else
          ndum = 1
        endif
        allocate( dumy(ndum,nk) )
        !----------------------
        if(myid.eq.0) print *,'  radbcw, ndum = ',ndum
        aname = 'radbcw'
        dumy = 0.0
        if(wbc.eq.2)  &
        call writerbcwe(radbcw,aname,ndum,dumy,ibw,jb,je,kb,ke,ny,ni,nj,nk,  &
                        nodex,nodey,restart_format,myid,numprocs,myj1p)
        if( myid.eq.0 )then
          call writebc(ndum=ndum,                      &
                       nk=nk,                          &
                       dum=dumy,                       &
                       restart_format=restart_format,  &
                       ncid=ncid,                      &
                       time_index=time_index,          &
                       aname=aname)
        endif
        !----------------------
        !cccccccccccccccccccccc
        !----------------------
        if(myid.eq.0) print *,'  radbce, ndum = ',ndum
        aname = 'radbce'
        dumy = 0.0
        if(ebc.eq.2)  &
        call writerbcwe(radbce,aname,ndum,dumy,ibe,jb,je,kb,ke,ny,ni,nj,nk,  &
                        nodex,nodey,restart_format,myid,numprocs,myj1p)
        if( myid.eq.0 )then
          call writebc(ndum=ndum,                      &
                       nk=nk,                          &
                       dum=dumy,                       &
                       restart_format=restart_format,  &
                       ncid=ncid,                      &
                       time_index=time_index,          &
                       aname=aname)
        endif
        deallocate( dumy )
        !----------------------
        !cccccccccccccccccccccc
        !----------------------
        if(myid.eq.0)then
          ndum = nx
        else
          ndum = 1
        endif
        allocate( dumx(ndum,nk) )
        !----------------------
        if(myid.eq.0) print *,'  radbcs, ndum = ',ndum
        aname = 'radbcs'
        dumx = 0.0
        if(sbc.eq.2)  &
        call writerbcsn(radbcs,aname,ndum,dumx,ibs,ib,ie,kb,ke,nx,ni,nj,nk,  &
                        nodex,nodey,restart_format,myid,numprocs,myi1p)
        if( myid.eq.0 )then
          call writebc(ndum=ndum,                      &
                       nk=nk,                          &
                       dum=dumx,                       &
                       restart_format=restart_format,  &
                       ncid=ncid,                      &
                       time_index=time_index,          &
                       aname=aname)
        endif
        !----------------------
        !cccccccccccccccccccccc
        !----------------------
        if(myid.eq.0) print *,'  radbcn, ndum = ',ndum
        aname = 'radbcn'
        dumx = 0.0
        if(nbc.eq.2)  &
        call writerbcsn(radbcn,aname,ndum,dumx,ibn,ib,ie,kb,ke,nx,ni,nj,nk,  &
                        nodex,nodey,restart_format,myid,numprocs,myi1p)
        if( myid.eq.0 )then
          call writebc(ndum=ndum,                      &
                       nk=nk,                          &
                       dum=dumx,                       &
                       restart_format=restart_format,  &
                       ncid=ncid,                      &
                       time_index=time_index,          &
                       aname=aname)
        endif
        deallocate( dumx )
        !----------------------
        !cccccccccccccccccccccc
        !----------------------

      else
        in_file = in_file_no
      endif  irbc4

      if(doinfo) call writeinfo(in_file,orecs,orecu,orecv,orecw,'  openbc:')

!---------------------------------------------------------------
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!---------------------------------------------------------------
!  150820:  optional variables


    IF( restart_file_theta )THEN
      in_file = in_file_yes
      do k=1,nk
      do j=1,nj
      do i=1,ni
        dum1(i,j,k) = th0(i,j,k)+tha(i,j,k)
      enddo
      enddo
      enddo
      call writer(ni,nj,1,nk,nx,ny,dum1(ib,jb,1),'theta   ',       &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
    ELSE
      in_file = in_file_no
    ENDIF

    if(doinfo) call writeinfo(in_file,orecs,orecu,orecv,orecw,'  theta:')

    !...................................................................

    IF( restart_file_u0 .or. do_lsnudge .or. do_adapt_move )THEN
      in_file = in_file_yes
      call writer(ni+1,nj,1,nk,nx+1,ny,u0(ib,jb,1),'u0      ',     &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecu,ufile,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2iu,d2ju,d3iu,d3ju)
    ELSE
      in_file = in_file_no
    ENDIF

    if(doinfo) call writeinfo(in_file,orecs,orecu,orecv,orecw,'  u0:')

    !...................................................................

    IF( restart_file_v0 .or. do_lsnudge .or. do_adapt_move )THEN
      in_file = in_file_yes
      call writer(ni,nj+1,1,nk,nx,ny+1,v0(ib,jb,1),'v0      ',     &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecv,vfile,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2iv,d2jv,d3iv,d3jv)
    ELSE
      in_file = in_file_no
    ENDIF

    if(doinfo) call writeinfo(in_file,orecs,orecu,orecv,orecw,'  v0:')

    !...................................................................

    IF( restart_file_dbz .and. qd_dbz.ge.1 )THEN
      ! cm1r19:
      call writer(ni,nj,1,nk,nx,ny,qdiag(ibdq,jbdq,1,qd_dbz),'dbz     ',        &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
    ENDIF
    !-----
    IF( restart_file_th0 )THEN
      call writer(ni,nj,1,nk,nx,ny,th0(ib,jb,1),'th0     ',        &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
    ENDIF
    IF( restart_file_prs0 )THEN
      call writer(ni,nj,1,nk,nx,ny,prs0(ib,jb,1),'prs0    ',       &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
    ENDIF
    IF( restart_file_pi0 )THEN
      call writer(ni,nj,1,nk,nx,ny,pi0(ib,jb,1),'pi0     ',        &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
    ENDIF
    IF( restart_file_rho0 )THEN
      call writer(ni,nj,1,nk,nx,ny,rho0(ib,jb,1),'rho0    ',       &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
    ENDIF
    IF( restart_file_qv0 )THEN
      call writer(ni,nj,1,nk,nx,ny,qv0(ib,jb,1),'qv0     ',        &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
    ENDIF
    !-----
    IF( restart_file_zs )THEN
      call writer(ni,nj,1,1,nx,ny,zs(ib,jb),'zs      ',            &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
    ENDIF
    IF( restart_file_zh )THEN
      call writer(ni,nj,1,nk,nx,ny,zh(ib,jb,1),'zhalf   ',         &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
    ENDIF
    IF( restart_file_zf )THEN
      call writer(ni,nj,1,nk+1,nx,ny,zf(ib,jb,1),'zfull   ',       &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecw,wfile,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
    ENDIF

    IF( restart_file_diags .and. td_diss.gt.0 )THEN
      call writer(ni,nj,1,nk,nx,ny,tdiag(ib,jb,1,td_diss),'dissheat',    &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,         &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p,       &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
    ENDIF
    IF( restart_file_diags .and.  td_mp.gt.0 )THEN
      call writer(ni,nj,1,nk,nx,ny,tdiag(ib,jb,1,td_mp),'mptend  ',      &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,         &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p,       &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
    ENDIF
    IF( restart_file_diags .and. qd_vtc.gt.0 )THEN
      call writer(ni,nj,1,nk,nx,ny,qdiag(ib,jb,1,qd_vtc),'vtc     ',     &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,         &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p,       &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
    ENDIF
    IF( restart_file_diags .and. qd_vtr.gt.0 )THEN
      call writer(ni,nj,1,nk,nx,ny,qdiag(ib,jb,1,qd_vtr),'vtr     ',     &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,         &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p,       &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
    ENDIF
    IF( restart_file_diags .and. qd_vts.gt.0 )THEN
      call writer(ni,nj,1,nk,nx,ny,qdiag(ib,jb,1,qd_vts),'vts     ',     &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,         &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p,       &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
    ENDIF
    IF( restart_file_diags .and. qd_vtg.gt.0 )THEN
      call writer(ni,nj,1,nk,nx,ny,qdiag(ib,jb,1,qd_vtg),'vtg     ',     &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,         &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p,       &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
    ENDIF
    IF( restart_file_diags .and. qd_vti.gt.0 )THEN
      call writer(ni,nj,1,nk,nx,ny,qdiag(ib,jb,1,qd_vti),'vti     ',     &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,sfile,         &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p,       &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
    ENDIF

!---------------------------------------------------------------
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!---------------------------------------------------------------


    IF( restart_format.eq.1 )THEN
      if(openrst)  close(unit=rstfile)
      if(openinfo) close(unit=infofile)
      if(opens)    close(unit=sfile)
      if(openu)    close(unit=ufile)
      if(openv)    close(unit=vfile)
      if(openw)    close(unit=wfile)
#ifdef NETCDF
    ELSEIF( restart_format.eq.2 )THEN
      if( myid.eq.0 )then
        call disp_err( nf90_close(ncid) , .true. )
      endif
#endif
    ENDIF

    if(myid.eq.0)then
      print *
      print *,' noutstat   = ',noutstat
      print *,' nout2ds    = ',nout2ds
      print *,' nout3ds    = ',nout3ds
      print *,' nout3du    = ',nout3du
      print *,' nout3dv    = ',nout3dv
      print *,' nout3dw    = ',nout3dw
      print *,' nout3dxv   = ',nout3dxv
      print *,' nout3dyv   = ',nout3dyv
      print *,' nout3dzv   = ',nout3dzv
      print *,' noutdiags  = ',noutdiags
      print *,' noutdiagw  = ',noutdiagw
      print *,' noutazims  = ',noutazims
      print *,' noutazimw  = ',noutazimw
      print *,' noutparcel = ',noutparcel
      print *,' nouthifrq  = ',nouthifrq
      print *,' nconfig    = ',nconfig
      print *
    endif

#ifdef MPI
      ! seems to help with memory on some hardware (somehow):
      call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif

      end subroutine restart_write


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


    subroutine writebc(ndum,nk,dum,restart_format,ncid,time_index,aname)
#ifdef NETCDF
    use netcdf
#endif
    implicit none

    integer, intent(in) :: ndum,nk
    real, intent(in) :: dum(ndum,nk)
    integer, intent(in) :: restart_format,ncid,time_index
    character(len=6), intent(in) :: aname

    integer :: k
#ifdef NETCDF
    integer :: varid,ncstatus
#endif

    if( restart_format.eq.1 )then
      do k=1,nk
        write(rstfile) dum(:,k)
      enddo
#ifdef NETCDF
    elseif( restart_format.eq.2 )then
      ncstatus = nf90_put_var(ncid,varid,dum,(/1,k,time_index/),(/ndum,1,1/))
      if(ncstatus.ne.nf90_noerr)then
        print *,'  Error2 in writebc, aname = ',aname
        print *,nf90_strerror(ncstatus)
        call stopcm1
      endif
      stop 33329
#endif
    endif

    end subroutine writebc


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine writerbcwe(radbc,aname,ndum,dumy,ibndy,jb,je,kb,ke,ny,ni,nj,nk,  &
                            nodex,nodey,restart_format,myid,numprocs,myj1p)
#ifdef MPI
      use mpi
#endif
      implicit none

      integer, intent(in) :: ndum,ibndy,jb,je,kb,ke,ny,ni,nj,nk,nodex,nodey
      real, intent(in), dimension(jb:je,kb:ke) :: radbc
      character(len=6), intent(in) :: aname
      real, intent(inout), dimension(ndum,nk) :: dumy
      integer, intent(in) :: restart_format,myid
      integer, intent(in) :: numprocs
      integer, intent(in), dimension(0:numprocs-1) :: myj1p

      integer :: j,j1,j2,k
#ifdef MPI
      integer :: fooj,proc,reqs,ierr,req1,njrecv
#endif

    kloop:  &
    do k=1,nk

#ifndef MPI
      do j=1,ny
        dumy(j,k) = radbc(j,k)
      enddo
#else
      IF(myid.ne.0)THEN
        if( ibndy.eq.1 )then
          call MPI_ISEND(nj,1,MPI_INTEGER,0,30,MPI_COMM_WORLD,req1,ierr)
          call MPI_ISEND(radbc(1,k),nj,MPI_REAL,0,31,MPI_COMM_WORLD,reqs,ierr)
          call MPI_WAIT(req1,mpi_status_ignore,ierr)
          call MPI_WAIT(reqs,mpi_status_ignore,ierr)
        endif
      ELSE
        if( (aname.eq.'radbcw') .or. (aname.eq.'radbce' .and. nodex.eq.1) )then
          do j=1,nj
            dumy(j,k) = radbc(j,k)
          enddo
          j1 = 2
          j2 = nodey
        else
          j1 = 1
          j2 = nodey
        endif
        do j=j1,j2
          if( aname.eq.'radbcw' )then
            proc = (j-1)*nodex
          else
            proc = (j-1)*nodex + (nodex-1)
          endif
          call MPI_IRECV(njrecv,1,MPI_INTEGER,proc,30,MPI_COMM_WORLD,reqs,ierr)
          call MPI_WAIT(reqs,mpi_status_ignore,ierr)
          fooj = myj1p(proc)
          call MPI_IRECV(dumy(fooj,k),njrecv,MPI_REAL,proc,31,MPI_COMM_WORLD,reqs,ierr)
          call MPI_WAIT(reqs,mpi_status_ignore,ierr)
        enddo
      ENDIF
#endif

    enddo  kloop

      end subroutine writerbcwe


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine writerbcsn(radbc,aname,ndum,dumx,ibndy,ib,ie,kb,ke,nx,ni,nj,nk,  &
                            nodex,nodey,restart_format,myid,numprocs,myi1p)
#ifdef MPI
      use mpi
#endif
      implicit none

      integer, intent(in) :: ndum,ibndy,ib,ie,kb,ke,nx,ni,nj,nk,nodex,nodey
      real, intent(in), dimension(ib:ie,kb:ke) :: radbc
      character(len=6), intent(in) :: aname
      real, intent(inout), dimension(ndum,nk) :: dumx
      integer, intent(in) :: restart_format,myid
      integer, intent(in) :: numprocs
      integer, intent(in), dimension(0:numprocs-1) :: myi1p

      integer :: i,i1,i2,k
#ifdef MPI
      integer :: fooi,proc,reqs,ierr,req1,nirecv
#endif

    kloop:  &
    do k=1,nk

#ifndef MPI
      do i=1,nx
        dumx(i,k) = radbc(i,k)
      enddo
#else
      IF(myid.ne.0)THEN
        if( ibndy.eq.1 )then
          call MPI_ISEND(ni,1,MPI_INTEGER,0,33,MPI_COMM_WORLD,req1,ierr)
          call MPI_ISEND(radbc(1,k),ni,MPI_REAL,0,32,MPI_COMM_WORLD,reqs,ierr)
          call MPI_WAIT(req1,mpi_status_ignore,ierr)
          call MPI_WAIT(reqs,mpi_status_ignore,ierr)
        endif
      ELSE
        if( (aname.eq.'radbcs') .or. (aname.eq.'radbcn' .and. nodey.eq.1) )then
          do i=1,ni
            dumx(i,k) = radbc(i,k)
          enddo
          i1 = 2
          i2 = nodex
        else
          i1 = 1
          i2 = nodex
        endif
        do i=i1,i2
          if( aname.eq.'radbcs' )then
            proc = (i-1)
          else
            proc = (i-1) + nodex*(nodey-1)
          endif
          call MPI_IRECV(nirecv,1,MPI_INTEGER,proc,33,MPI_COMM_WORLD,reqs,ierr)
          call MPI_WAIT(reqs,mpi_status_ignore,ierr)
          fooi = myi1p(proc)
          call MPI_IRECV(dumx(fooi,k),nirecv,MPI_REAL,proc,32,MPI_COMM_WORLD,reqs,ierr)
          call MPI_WAIT(reqs,mpi_status_ignore,ierr)
        enddo
      ENDIF
#endif

    enddo  kloop

      end subroutine writerbcsn


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


    subroutine writer(numi,numj,numk1,numk2,nxr,nyr,var,aname,           &
                      ni,nj,ngxy,myid,numprocs,nodex,nodey,orec,nfile,   &
                      ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p,   &
                      dat1,dat2,dat3,reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2i,d2j,d3i,d3j)
#ifdef MPI
    use mpi
#endif
#ifdef NETCDF
    use netcdf
#endif
    implicit none

    !-------------------------------------------------------------------
    ! This subroutine collects data (from other processors if this is a
    ! MPI run) and does the actual writing of restart files.
    !-------------------------------------------------------------------

    integer, intent(in) :: numi,numj,numk1,numk2,nxr,nyr
    integer, intent(in) :: ppnode,d3n,d3t,d2i,d2j,d3i,d3j
    real, intent(in   ), dimension(1-ngxy:numi+ngxy,1-ngxy:numj+ngxy,numk1:numk2) :: var
    character(len=8), intent(in) :: aname
    integer, intent(in) :: ni,nj,ngxy,myid,numprocs,nodex,nodey
    integer, intent(inout) :: orec,ncid
    integer, intent(in) :: time_index,restart_format,restart_filetype
    real, intent(inout), dimension(d3i,d3j) :: dat1
    real, intent(inout), dimension(d2i,d2j) :: dat2
    real, intent(inout), dimension(d3i,d3j,0:d3n-1) :: dat3
    integer, intent(inout), dimension(d3t) :: reqt
    integer, intent(in) :: mynode,nodeleader,nodes,nfile
    integer, intent(in), dimension(0:numprocs-1) :: myi1p,myi2p,myj1p,myj2p

    integer :: i,j,k,msk,varid
#ifdef MPI
    integer :: ierr
#endif
#ifdef NETCDF
    integer :: status
#endif

!-------------------------------------------------------------------------------

#ifdef NETCDF
    if( restart_format.eq.2 )then
      status = nf90_inq_varid(ncid,aname,varid)
      if(status.ne.nf90_noerr)then
        print *,'  Error1 in writer, aname = ',aname
        print *,nf90_strerror(status)
        call stopcm1
      endif
    endif
#endif

    msk = 0

    if(myid.eq.0) print *,aname,numk1,numk2

    rf1:  &
    IF( restart_filetype.eq.1 .or. restart_filetype.eq.2 )THEN
      ! one file:

#ifndef MPI

        call   writer_sp(numk1,numk2,d2i,d2j,numi,numj,ngxy,dat2,var,orec,nfile,ncid,varid,aname,time_index,restart_format)

#else

      if(myid.ne.nodeleader)then
        call   writer_comm1(numk1,numk2,d3i,d3j,ngxy,numi,numj,nodeleader,dat1,var)
      else
        if(myid.ne.msk)then
          call writer_comm2(ppnode,d2i,d2j,d3i,d3j,d3n,d3t,ngxy,numi,numj,numk1,numk2,msk,myid,reqt,dat2(1,1),dat3,var,nfile,orec)
        else
          call  writer_comm3(ppnode,d2i,d2j,d3i,d3j,d3n,d3t,ngxy,numi,numj,numk1,numk2,nfile,orec,reqt,dat1,dat2,dat3,var,  &
                             restart_format,myid,ni,nj,numprocs,mynode,nodes,myi1p,myi2p,myj1p,myj2p,ncid,varid,time_index,aname)
        endif
      endif

#endif

    ENDIF  rf1

#ifdef MPI
!-------------------------------------------------------------------------------

    rf3:  &
    IF( restart_filetype.eq.3 )THEN
      ! one file per processor:

      IF( restart_format.eq.1 )THEN

        kloop4:  &
        DO k=numk1,numk2
          write(nfile,rec=orec) ((var(i,j,k),i=1,numi),j=1,numj)
          orec = orec+1
        ENDDO  kloop4

      ELSE

        stop 95642

      ENDIF

    ENDIF  rf3

!-------------------------------------------------------------------------------

    rf4:  &
    IF( restart_filetype.eq.4 )THEN
      ! one file per node:

      call     writer2(numi,numj,numk1,numk2,nxr,nyr,var,aname,          &
                      ni,nj,ngxy,myid,numprocs,nodex,nodey,orec,nfile,   &
                      ncid,time_index,restart_format,myi1p,myi2p,myj1p,myj2p,   &
                      dat1(1,1),dat2(1,1),dat3(1,1,0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2i,d2j,d3i,d3j)

    ENDIF  rf4
#endif

!-------------------------------------------------------------------------------
!ccccc  done  cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!-------------------------------------------------------------------------------

    end subroutine writer


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    subroutine writeit_restart(d2i,d2j,dat2,k,numk1,numk2,restart_format,nfile,orec,ncid,varid,time_index)
#ifdef NETCDF
    use netcdf
#endif
    implicit none

    integer, intent(in) :: d2i,d2j
    real, intent(in), dimension(d2i,d2j) :: dat2
    integer, intent(in) :: k,numk1,numk2,restart_format,nfile
    integer, intent(inout) :: orec
    integer, intent(in) :: ncid,varid,time_index

    integer :: i,j
#ifdef NETCDF
    integer :: status
#endif

    IF( restart_format.eq.1 )THEN

      write(nfile,rec=orec) ((dat2(i,j),i=1,d2i),j=1,d2j)
      orec = orec+1

#ifdef NETCDF
    ELSEIF( restart_format.eq.2 )THEN

      if(numk1.eq.numk2)then
        status = nf90_put_var(ncid,varid,dat2,(/1,1,time_index/),(/d2i,d2j,1/))
      else
        status = nf90_put_var(ncid,varid,dat2,(/1,1,k,time_index/),(/d2i,d2j,1,1/))
      endif
      if(status.ne.nf90_noerr)then
        print *,'  Error2 in writeit_restart '
        print *,nf90_strerror(status)
        call stopcm1
      endif
#endif

    ENDIF

    end subroutine writeit_restart

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
#ifndef MPI

    subroutine writer_sp(numk1,numk2,d2i,d2j,numi,numj,ngxy,dat2,var,orec,nfile,ncid,varid,aname,time_index,restart_format)
#ifdef NETCDF
    use netcdf
#endif
    implicit none

    integer, intent(in) :: numk1,numk2,d2i,d2j,numi,numj,ngxy
    real, intent(inout), dimension(d2i,d2j) :: dat2
    real, intent(in   ), dimension(1-ngxy:numi+ngxy,1-ngxy:numj+ngxy,numk1:numk2) :: var
    integer, intent(inout) :: orec,ncid,varid
    integer, intent(in) :: nfile
    character(len=8), intent(in) :: aname
    integer, intent(in) :: time_index,restart_format

    integer :: i,j,k

    kloop:  &
    DO k=numk1,numk2

      !$omp parallel do default(shared)   &
      !$omp private(i,j)
      do j=1,numj
      do i=1,numi
        dat2(i,j)=var(i,j,k)
      enddo
      enddo

      call     writeit_restart(d2i,d2j,dat2,k,numk1,numk2,restart_format,nfile,orec,ncid,varid,time_index)

    ENDDO  kloop

    end subroutine writer_sp

#endif
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
#ifdef MPI


    ! ordinary processor ... send data to nodeleader:
    subroutine writer_comm1(numk1,numk2,d3i,d3j,ngxy,numi,numj,nodeleader,dat1,var)
    use mpi
    implicit none

    integer, intent(in) :: numk1,numk2,d3i,d3j,ngxy,numi,numj,nodeleader
    real, intent(inout), dimension(d3i,d3j) :: dat1
    real, intent(in   ), dimension(1-ngxy:numi+ngxy,1-ngxy:numj+ngxy,numk1:numk2) :: var

    integer :: i,j,k,tag,reqs,ierr

    tag = 1

    kloop:  &
    DO k=numk1,numk2

      !$omp parallel do default(shared)   &
      !$omp private(i,j)
      do j=1,numj
      do i=1,numi
        dat1(i,j)=var(i,j,k)
      enddo
      enddo

      call MPI_ISEND(dat1(1,1),d3i*d3j,MPI_REAL,nodeleader,tag,MPI_COMM_WORLD,reqs,ierr)
      call MPI_WAIT(reqs,MPI_STATUS_IGNORE,ierr)

      tag = tag+2

    ENDDO  kloop

    end subroutine writer_comm1
    ! done, ordinary processor

#endif
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
#ifdef MPI

    subroutine writer_comm2(ppnode,d2i,d2j,d3i,d3j,d3n,d3t,ngxy,numi,numj,numk1,numk2,msk,myid,reqt,dat2     ,dat3,var,nfile,orec)
    use mpi
    implicit none

    integer, intent(in) :: ppnode,d2i,d2j,d3i,d3j,d3n,d3t,ngxy,numi,numj,numk1,numk2,msk,myid
    integer, intent(inout), dimension(d3t) :: reqt
    real, intent(inout), dimension(d3i*ppnode,d3j) :: dat2
    real, intent(inout), dimension(d3i,d3j,0:d3n-1) :: dat3
    real, intent(in   ), dimension(1-ngxy:numi+ngxy,1-ngxy:numj+ngxy,numk1:numk2) :: var
    integer, intent(in)    :: nfile
    integer, intent(inout) :: orec

    integer :: i,j,k,tag,reqs,ierr,proc,fooi
    integer, dimension(mpi_status_size,ppnode-1) :: status1

    tag = 1

    kloop:  &
    DO k=numk1,numk2

      do proc=myid+1,myid+(ppnode-1)
        call MPI_IRECV(dat3(1,1,proc),d3i*d3j,MPI_REAL,proc,tag,MPI_COMM_WORLD,reqt(proc-myid),ierr)
      enddo

      !$omp parallel do default(shared)  &
      !$omp private(i,j)
      do j=1,numj
      do i=1,numi
        dat3(i,j,myid)=var(i,j,k)
      enddo
      enddo

      call mpi_waitall(ppnode-1,reqt(1:ppnode-1),status1,ierr)
      call MPI_ISEND(dat3(1,1,myid),d3i*d3j*ppnode,MPI_REAL,msk,tag+1,MPI_COMM_WORLD,reqs,ierr)
      call MPI_WAIT(reqs,MPI_STATUS_IGNORE,ierr)

      tag = tag+2

    ENDDO  kloop

    end subroutine writer_comm2

#endif
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
#ifdef MPI


    ! nodeleader and proc 0:
    subroutine  writer_comm3(ppnode,d2i,d2j,d3i,d3j,d3n,d3t,ngxy,numi,numj,numk1,numk2,nfile,orec,reqt,dat1,dat2,dat3,var,  &
                             restart_format,myid,ni,nj,numprocs,mynode,nodes,myi1p,myi2p,myj1p,myj2p,ncid,varid,time_index,aname)
    use mpi
#ifdef NETCDF
    use netcdf
#endif
    implicit none

    integer, intent(in) :: ppnode,d2i,d2j,d3i,d3j,d3n,d3t,ngxy,numi,numj,numk1,numk2,nfile
    integer, intent(inout) :: orec
    integer, intent(inout), dimension(d3t) :: reqt
    real, intent(inout), dimension(d3i,d3j) :: dat1
    real, intent(inout), dimension(d2i,d2j) :: dat2
    real, intent(inout), dimension(d3i,d3j,0:d3n-1) :: dat3
    real, intent(in),    dimension(1-ngxy:numi+ngxy,1-ngxy:numj+ngxy,numk1:numk2) :: var
    integer, intent(in) :: restart_format,myid,ni,nj,numprocs,mynode,nodes
    integer, intent(in), dimension(0:numprocs-1) :: myi1p,myi2p,myj1p,myj2p
    integer, intent(in) :: ncid,varid,time_index
    character(len=8), intent(in) :: aname

    integer :: i,j,k,tag,reqs,ierr
    integer :: fooi,fooj,n,nn,nnn,ntot,n1,n2,nitmp,njtmp,proc,index,index2
    integer, dimension(mpi_status_size) :: status0
    integer, dimension(mpi_status_size,ppnode-1) :: status1
    integer, dimension(mpi_status_size,ppnode-1 + nodes-1) :: status2
!!!#ifdef NETCDF
!!!    integer :: status
!!!#endif

    tag = 1

    kloop:  &
    DO k=numk1,numk2

      do proc=myid+1,myid+(ppnode-1)
        call MPI_IRECV(dat3(1,1,proc),d3i*d3j,MPI_REAL,proc,tag,MPI_COMM_WORLD,reqt(proc-myid),ierr)
      enddo
      do n = 1,(nodes-1)
        proc = n*ppnode
        call MPI_IRECV(dat3(1,1,proc),d3i*d3j*ppnode,MPI_REAL,proc,tag+1,MPI_COMM_WORLD,reqt(ppnode-1+n),ierr)
      enddo


        !$omp parallel do default(shared)  &
        !$omp private(i,j)
        do j=1,numj
        do i=1,numi
          dat2(i,j)=var(i,j,k)
        enddo
        enddo

          ntot = ppnode-1 + nodes-1
          do nn=1,ntot
            index = nn
            call mpi_wait(reqt(nn),status0,ierr)
!!!            call mpi_waitany(ntot,reqt(1:ntot),index,MPI_STATUS_IGNORE,ierr)
            if( index.le.(ppnode-1) )then
              ! data from ordinary procs on node:
              proc = myid+index
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
                dat2(fooi+i,fooj+j) = dat3(i,j,proc)
              enddo
              enddo
            else
              ! data from other nodeleaders:
              index2 = index-(ppnode-1)
              if( index2.le.mynode )then
                index2 = index2-1
              endif
              n1 = index2*ppnode
              n2 = (index2+1)*ppnode-1
              do nnn = n1,n2
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
                  dat2(fooi+i,fooj+j) = dat3(i,j,proc)
                enddo
                enddo
              enddo
            endif
          enddo

      call     writeit_restart(d2i,d2j,dat2,k,numk1,numk2,restart_format,nfile,orec,ncid,varid,time_index)

      tag = tag+2

    ENDDO  kloop

    end subroutine writer_comm3

#endif
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
#ifdef MPI

    ! cm1r17-format restart files !
    subroutine writer2(numi,numj,numk1,numk2,nxr,nyr,var,aname,          &
                      ni,nj,ngxy,myid,numprocs,nodex,nodey,orec,nfile,   &
                      ncid,time_index,restart_format,myi1p,myi2p,myj1p,myj2p,   &
                      dat1     ,dat2     ,dat3       ,reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2i,d2j,d3i,d3j)
    use mpi
    implicit none

    !-------------------------------------------------------------------
    ! This subroutine collects data (from other processors if this is a
    ! MPI run) and does the actual writing of restart files.
    !-------------------------------------------------------------------

    integer, intent(in) :: numi,numj,numk1,numk2,nxr,nyr
    integer, intent(in) :: ppnode,d3n,d3t,d2i,d2j,d3i,d3j
    real, intent(in   ), dimension(1-ngxy:numi+ngxy,1-ngxy:numj+ngxy,numk1:numk2) :: var
    character(len=8), intent(in) :: aname
    integer, intent(in) :: ni,nj,ngxy,myid,numprocs,nodex,nodey
    integer, intent(inout) :: orec,ncid
    integer, intent(in) :: time_index,restart_format
    real, intent(inout), dimension(d3i,d3j) :: dat1
    real, intent(inout), dimension(d3i*ppnode,d3j) :: dat2
    real, intent(inout), dimension(d3i,d3j,0:d3n-1) :: dat3
    integer, intent(inout), dimension(d3t) :: reqt
    integer, intent(in) :: mynode,nodeleader,nodes,nfile
    integer, intent(in), dimension(0:numprocs-1) :: myi1p,myi2p,myj1p,myj2p

    integer :: i,j,k,msk
    integer :: reqs,index,index2,n,nn,nnn,fooi,fooj,proc,ierr,ntot,n1,n2,tag
    logical :: recv1,recv2

    DO k=numk1,numk2
      IF(myid.ne.nodeleader)THEN
        !$omp parallel do default(shared)   &
        !$omp private(i,j)
        do j=1,numj
        do i=1,numi
          dat1(i,j) = var(i,j,k)
        enddo
        enddo
        call MPI_ISEND(dat1(1,1),d3i*d3j,MPI_REAL,nodeleader,k,MPI_COMM_WORLD,reqs,ierr)
        call MPI_WAIT(reqs,mpi_status_ignore,ierr)
      ELSE
        do proc=myid+1,myid+(ppnode-1)
          call MPI_IRECV(dat3(1,1,proc),d3i*d3j,MPI_REAL,proc,k,MPI_COMM_WORLD,reqt(proc-myid),ierr)
        enddo
        !$omp parallel do default(shared)   &
        !$omp private(i,j)
        do j=1,numj
        do i=1,numi
          dat2(i,j)=var(i,j,k)
        enddo
        enddo
        nn = 1
        do while( nn.le.(ppnode-1) )
          nn = nn + 1
          call mpi_waitany(ppnode-1,reqt(1:ppnode-1),index,MPI_STATUS_IGNORE,ierr)
          fooi = d3i*index
          !$omp parallel do default(shared)   &
          !$omp private(i,j)
          do j=1,numj
          do i=1,numi
            dat2(fooi+i,j)=dat3(i,j,nodeleader+index)
          enddo
          enddo
        enddo
        write(nfile,rec=orec) ((dat2(i,j),i=1,d3i*ppnode),j=1,d3j)
        orec = orec+1
      ENDIF
    ENDDO

    end subroutine writer2


#endif
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    subroutine writeinfo(in_file,orecs,orecu,orecv,orecw,text,num)
    use input, only : myid
    implicit none

    integer, intent(in) :: in_file,orecs,orecu,orecv,orecw
    character(len=*), intent(in) :: text
    integer, intent(in), optional :: num

    write(infofile) in_file
    write(infofile) orecs
    write(infofile) orecu
    write(infofile) orecv
    write(infofile) orecw
    if(myid.eq.0) print *,'  inside writeinfo:'
    if(myid.eq.0) print *,'    text    = ',trim(text)
    if(myid.eq.0) print *,'    in_file = ',in_file
    if(myid.eq.0) print *,'    orecs   = ',orecs
    if(myid.eq.0) print *,'    orecu   = ',orecu
    if(myid.eq.0) print *,'    orecv   = ',orecv
    if(myid.eq.0) print *,'    orecw   = ',orecw
    if(present(num))then
      write(infofile) num
      if(myid.eq.0) print *,'    num     = ',num
    endif

    end subroutine writeinfo

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  END MODULE restart_write_module
