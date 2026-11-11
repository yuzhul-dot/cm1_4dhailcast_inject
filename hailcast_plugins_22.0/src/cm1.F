
      program cm1

!-----------------------------------------------------------------------------
!
!  CM1 Numerical Model, Release 22.0  (cm1r22.0)
!  21 September 2026
!  https://github.com/NCAR/CM1
!
!  (c)2026 - University Corporation for Atmospheric Research 
!
!-----------------------------------------------------------------------------
!
!  Please see documentation at the top of the "solve.F" file.
!
!  See also documentation at the cm1 website, such as:
!
!    "The governing equations for CM1"
!        https://github.com/NCAR/CM1/blob/main/docs/cm1_equations.pdf
!
!-----------------------------------------------------------------------------

      use input
      use constants, only : pi
      use param_module, only : param,wenocheck
      use base_module, only : base
      use init3d_module, only : init3d,getset
      use misclibs, only : calccflquick,calcksquick,getnewdt,dtsmall,pdefq,  &
                           set_time_to_zero,tavg_plus1,tavg_shift
      use solve1_module, only : solve1
      use solve2_module, only : solve2
      use solve3_module, only : solve3
      use pdcomp_module, only : pidcomp
      use mp_driver_module, only : mp_driver
      use diff2_module, only : diff2def
      use turb_module, only : sfc_and_turb
      use statpack_module, only : statpack,setup_stat_vars
      use writeout_module, only : writeout
      use restart_write_module, only : restart_write
      use restart_read_module, only : restart_read
      use radiation_module, only : radiation_driver
      use radtrns3d_module, only : nrad2d,n2d_radiat,n3d_radiat
      use domaindiag_module, only : domaindiag
      use azimavg_module, only : azimavg,getrangle,getpsmth,getpcent,update_adapt_move
      use hifrq_module, only : writeout_hifrq
      use hailcast_module
      use parcel_module, only : setup_parcel_vars,parcel_interp,parcel_write
      use init_physics_module, only : init_physics
      use init_surface_module, only : init_surface
      use ib_module, only : ibib,ieib,jbib,jeib,kbib,keib,kmaxib, &
                            ib_flx_init,init_immersed_boundaries
      use eddy_recycle, only : irecywe,irecysn,jrecywe,jrecysn,krecy,nrecy,do_recycle
      use lsnudge_module, only : lsnudge_unit,lsnudge_count,lsnudge_maxlevels,        &
                                 lsnudge_time1,lsnudge_time2,lsnudge_nmax,lsnudgefac, &
                                 read_lsnudge
#ifdef MPI
      use mpi
#endif
#ifdef OPENMP
      use omp_lib
#endif
      implicit none

      integer :: nstep,nstep0
      integer :: srec,sirec,urec,vrec,wrec,xrec,yrec,zrec,nrec,mrec,prec,trecs,trecw,arecs,arecw,hrec
      integer :: nrst,nstatout,nwrite,nwritet,nwritea,nwriteh
      integer :: rbufsz,num_soil_layers,ndt
      real :: dt,dtlast
      double precision :: mtime,stattim,taptim,rsttim,radtim,prcltim,diagtim,azimavgtim,hifrqtim,restart_time,hailtrajtim
      double precision :: adt,adtlast,acfl,dbldt
      double precision :: mass1,mass2
      double precision :: avgsfcu,avgsfcv,avgsfcs,avgsfcsu,avgsfcsv,avgsfct,avgsfcq,avgsfcp
      logical :: dosfcflx
      logical, dimension(maxq) :: cloudvar,rhovar
      real, dimension(maxq) :: qmag
      character(len=3), dimension(maxq) :: qname
      character(len=20), dimension(maxq) :: qunit
      character(len=6), dimension(maxq) :: budname
      character(len=60), dimension(maxvars) :: desc_output
      character(len=40), dimension(maxvars) :: name_output,unit_output
      character(len=1),  dimension(maxvars) :: grid_output
      integer,           dimension(maxvars) :: levs_output
      logical, dimension(maxvars) :: cmpr_output
      character(len=40), dimension(maxvars) :: name_stat,desc_stat,unit_stat
      character(len=40), dimension(maxvars) :: name_prcl,desc_prcl,unit_prcl
      character(len=40), dimension(maxvars) :: name_hail,desc_hail,unit_hail
      double precision, dimension(:), allocatable :: bud,bud2
      double precision, dimension(:), allocatable :: qbudget
      double precision, dimension(:), allocatable :: asq,bsq
      real, dimension(:), allocatable :: xh,rxh,arh1,arh2,uh,ruh
      real, dimension(:), allocatable :: xf,rxf,arf1,arf2,uf,ruf
      real, dimension(:), allocatable :: yh,vh,rvh
      real, dimension(:), allocatable :: yf,vf,rvf
      real, dimension(:), allocatable :: xfref,xhref,yfref,yhref
      double precision, dimension(:), allocatable :: dumk1,dumk2,dumk3,dumk4
      real, dimension(:), allocatable :: dumk5
      real, dimension(:), allocatable :: rds,sigma,rdsf,sigmaf
      real, dimension(:), allocatable :: wprof,ufrc,vfrc,thfrc,qvfrc,ug,vg,dvdr,  &
                                         uavg,vavg,savg,thavg,pavg,ulspg,vlspg
      real, dimension(:,:), allocatable :: qavg
      double precision, dimension(:,:), allocatable :: cavg
      real, dimension(:,:,:), allocatable :: tauh,taus,zh,mh,rmh,c1,c2
      real, dimension(:,:,:), allocatable :: tauf,zf,mf,rmf
      real, dimension(:), allocatable :: rstat
      real, dimension(:,:), allocatable :: rho0s,pi0s,prs0s,rth0s
      real, dimension(:,:,:), allocatable :: pi0,rho0,prs0,thv0,th0,rth0,qv0
      real, dimension(:,:,:), allocatable :: qc0,qi0,rr0,rf0,rrf0,u0,v0,thrd
      real, dimension(:,:,:), allocatable :: dum1,dum2,dum3,dum4,dum5,dum6,dum7,dum8,dum9,dum10
      real, dimension(:,:), allocatable :: zs,gz,rgz,gzu,rgzu,gzv,rgzv,dzdx,dzdy
      real, dimension(:,:,:), allocatable :: gx,gxu,gy,gyv
      real, dimension(:,:,:), allocatable :: rain,raincv,hail,sws,svs,sps,srs,sgs,sus,shs
      real, dimension(:,:), allocatable :: tsk,znt,rznt,zntmp,ust,stau,tst,qst,z0t,z0q,thflux,qvflux,  &
                                           cd,ch,cq,z1,u1,v1,s1,th1,qv1,xland,psfc,tlh,f2d,psmth,prate,pratecv,ustt,ut,vt,st,cm0,eimaskw,eimaske,eimasks,eimaskn
      real, dimension(:,:), allocatable :: radbcw,radbce
      real, dimension(:,:), allocatable :: radbcs,radbcn
      real, dimension(:,:), allocatable :: dtu,dtu0,dtv,dtv0
      real, dimension(:,:,:), allocatable :: divx,rho,rr,rf,prs
      real, dimension(:,:,:), allocatable :: t11,t12,t13,t22,t23,t33
      real, dimension(:,:,:), allocatable :: m11,m12,m13,m22,m23,m33
      real, dimension(:,:,:), allocatable :: rru,ua,u3d,uten,uten1
      real, dimension(:,:,:), allocatable :: rrv,va,v3d,vten,vten1
      real, dimension(:,:,:), allocatable :: rrw,wa,w3d,wten,wten1
      real, dimension(:,:,:), allocatable :: ppi,pp3d,ppten,sten,sadv,ppx,phi1,phi2
      real, dimension(:,:,:), allocatable :: tha,th3d,thten,thten1,thterm
      real, dimension(:,:,:), allocatable :: qpten,qtten,qvten,qcten
      real, dimension(:,:,:,:), allocatable :: qa,q3d,qten
      real, dimension(:,:,:), allocatable :: p3a
      real, dimension(:,:,:,:), allocatable :: p3o
      real, dimension(:,:,:), allocatable :: kmh,kmv,khh,khv,cme,csm,ce1,ce2,kmv2,khv2
      real, dimension(:,:,:), allocatable :: tkea,tke3d,tketen
      real, dimension(:,:,:), allocatable :: nm,defv,defh,lenscl,dissten,epst
      real, dimension(:,:,:), allocatable :: thpten,qvpten,qcpten,qipten,upten,vpten,qnipten,qncpten,qrpten,qspten,qgpten
      real, dimension(:,:,:), allocatable :: xkzh,xkzq,xkzm
      real, dimension(:,:,:), allocatable :: swten,lwten,swtenc,lwtenc,thraten,cldfra,o30
      real, dimension(:,:), allocatable :: zir,radsw,rnflx,radswnet,radlwin,dsr,olr
      real, dimension(:,:,:), allocatable :: rad2d,effc,effi,effs,effr,effg,effis
      real, dimension(:,:), allocatable :: lwupt,lwuptc,lwdnt,lwdntc,lwupb,lwupbc,lwdnb,lwdnbc
      real, dimension(:,:), allocatable :: swupt,swuptc,swdnt,swdntc,swupb,swupbc,swdnb,swdnbc
      real, dimension(:,:), allocatable :: lwcf,swcf,coszr
      real, dimension(:,:), allocatable :: xice,xsnow,xlat,xlong,coszen,swddir,swddni,swddif,hrang
      integer, dimension(:,:,:), allocatable :: cldfra1_flag
      integer, dimension(:,:), allocatable :: lu_index,kpbl2d
      real, dimension(:,:), allocatable :: u10,v10,s10,hfx,qfx,               &
                                      hpbl,wspd,phim,phih,psim,psih,psiq,gz1oz0,br,brcr, &
                                      CHS,CHS2,CQS2,CPMM,ZOL,MAVAIL,          &
                                      MOL,RMOL,REGIME,LH,FLHC,FLQC,QGH,       &
                                      CK,CKA,CDA,USTM,QSFC,T2,Q2,TH2,EMISS,THC,ALBD,   &
                                      gsw,glw,chklowq,capg,snowc,snowh,qcg,dsxy,wstar,delta,prkpp,fm,fh
      real, dimension(:,:), allocatable :: charn,msang,scurx,scury,zkmax,cd_out,ch_out,wscale,wscaleu
      real, dimension(:,:), allocatable :: mznt,swspd,smois,taux,tauy,hpbl2d,evap2d,heat2d

      ! cm1r21.1: for cu param
      real, dimension(:,:,:), allocatable :: thcuten,qvcuten,qccuten,qicuten,ucuten,vcuten,thften,qvften

      ! for pbl:
      real, dimension(:,:,:,:), allocatable :: pbl3d
      real, dimension(:,:,:), allocatable :: pbl2d
      integer, dimension(:,:,:), allocatable :: pbl2dint

      ! for physics:
      real, dimension(:,:,:), allocatable :: tmp2d
      real, dimension(:,:), allocatable :: tmp1d

      real, dimension(:), allocatable :: slab_zs,slab_dzs
      real, dimension(:,:,:), allocatable :: tslb
      real, dimension(:,:), allocatable :: tmn,tml,t0ml,hml,h0ml,huml,hvml,tmoml
      real, dimension(:,:,:,:),  allocatable :: pta,pt3d,ptten
      real, dimension(:,:), allocatable :: dat1,dat2
      real, dimension(:,:,:), allocatable :: dat3
      integer, dimension(:), allocatable :: reqt
      integer :: reqc,reqk
      real, dimension(:,:), allocatable :: pdata,ploc
      real, dimension(:,:), allocatable :: haildata,hailloc
      logical, dimension(:,:,:), allocatable :: flag

!--- arrays for MPI ---
      integer, dimension(:), allocatable :: myi1p,myi2p,myj1p,myj2p
      integer, dimension(:), allocatable :: myi1pp,myi2pp,myj1pp,myj2pp
      integer, dimension(:), allocatable :: reqs_u,reqs_v,reqs_w,reqs_s,reqs_p,reqs_p2,reqs_p3,reqs_x,reqs_y,reqs_z,reqs_tk,reqs_tk2
      integer, dimension(:,:),  allocatable :: reqs_q,reqs_t
      real, dimension(:), allocatable :: nw1,nw2,ne1,ne2,sw1,sw2,se1,se2
      real, dimension(:,:,:), allocatable :: n3w1,n3w2,n3e1,n3e2,s3w1,s3w2,s3e1,s3e2
      real, dimension(:,:), allocatable :: ww1,ww2,we1,we2
      real, dimension(:,:), allocatable :: ws1,ws2,wn1,wn2
      real, dimension(:,:), allocatable :: pw1,pw2,pe1,pe2
      real, dimension(:,:), allocatable :: ps1,ps2,pn1,pn2
      real, dimension(:,:), allocatable :: p2w1,p2w2,p2e1,p2e2
      real, dimension(:,:), allocatable :: p2s1,p2s2,p2n1,p2n2
      real, dimension(:,:), allocatable :: p3w1,p3w2,p3e1,p3e2
      real, dimension(:,:), allocatable :: p3s1,p3s2,p3n1,p3n2
      real, dimension(:,:), allocatable :: vw1,vw2,ve1,ve2
      real, dimension(:,:), allocatable :: vs1,vs2,vn1,vn2
      real, dimension(:,:), allocatable :: zw1,zw2,ze1,ze2
      real, dimension(:,:), allocatable :: zs1,zs2,zn1,zn2
      real, dimension(:,:,:), allocatable :: uw31,uw32,ue31,ue32
      real, dimension(:,:,:), allocatable :: us31,us32,un31,un32
      real, dimension(:,:,:), allocatable :: vw31,vw32,ve31,ve32
      real, dimension(:,:,:), allocatable :: vs31,vs32,vn31,vn32
      real, dimension(:,:,:), allocatable :: ww31,ww32,we31,we32
      real, dimension(:,:,:), allocatable :: ws31,ws32,wn31,wn32
      real, dimension(:,:,:), allocatable :: sw31,sw32,se31,se32
      real, dimension(:,:,:), allocatable :: ss31,ss32,sn31,sn32
      real, dimension(:,:,:), allocatable :: rw31,rw32,re31,re32
      real, dimension(:,:,:), allocatable :: rs31,rs32,rn31,rn32
      real, dimension(:,:,:,:), allocatable :: qw31,qw32,qe31,qe32
      real, dimension(:,:,:,:), allocatable :: qs31,qs32,qn31,qn32
      real, dimension(:,:,:), allocatable :: tkw1,tkw2,tke1,tke2
      real, dimension(:,:,:), allocatable :: tks1,tks2,tkn1,tkn2
      real, dimension(:,:,:), allocatable :: kw1,kw2,ke1,ke2
      real, dimension(:,:,:), allocatable :: ks1,ks2,kn1,kn2
      real, dimension(:,:,:,:), allocatable :: tw1,tw2,te1,te2
      real, dimension(:,:,:,:), allocatable :: ts1,ts2,tn1,tn2

      ! arrays for elliptic solver:
      real, dimension(:,:,:),    allocatable :: cfb
      real, dimension(:),        allocatable :: cfa,cfc,d1,d2
      complex, dimension(:,:,:), allocatable :: pdt,lgbth,lgbph

      ! arrays for parallelized elliptic solver:
      integer, dimension(:), allocatable :: npxa,npya,npza
      integer, dimension(:), allocatable :: npxb,npyb,npzb
      integer, dimension(:), allocatable :: pmyi1a,pmyi2a,pmyj1a,pmyj2a,pmyk1a,pmyk2a
      integer, dimension(:), allocatable :: pmyi1b,pmyi2b,pmyj1b,pmyj2b,pmyk1b,pmyk2b
      real, dimension(:,:,:), allocatable :: pgrid1
      complex, dimension(:,:,:), allocatable :: cgrid1,cgrid2
      real, dimension(:,:,:,:), allocatable :: psenda,precva
      complex, dimension(:,:,:,:), allocatable :: psendb,precvb
      complex, dimension(:,:,:,:), allocatable :: psendc,precvc
      real, dimension(:,:,:,:), allocatable :: psendd,precvd
      double complex, dimension(:), allocatable :: rhsx,transx
      double complex, dimension(:), allocatable :: rhsy,transy
      integer, dimension(:), allocatable :: npcomm
      integer, dimension(:,:), allocatable :: nplist,npi1,npi2,npj1,npj2,npk1,npk2

      ! diagnostic arrays:
      real, dimension(:,:,:,:), allocatable :: tdiag,qdiag,udiag,vdiag,wdiag,kdiag,pdiag

      ! miscellaneous output:
      real, dimension(:,:,:),   allocatable :: out2d
      real, dimension(:,:,:,:), allocatable :: out3d

      ! for output using cylindrical coords:
      integer, dimension(:,:), allocatable :: cir
      real, dimension(:,:), allocatable :: crr,cangle
      real, dimension(:), allocatable :: crref

      ! t2p arrays:
      real, dimension(:), allocatable :: gamk,gamwall,kmw,ufw,vfw,u1b,v1b,l2p,s2p,s2b,t2pm1,t2pm2,t2pm3,t2pm4
      real, dimension(:,:,:), allocatable :: u2pt,v2pt,kmwk,ufwk,vfwk

      ! cm1ib arrays:
      logical, dimension(:,:,:), allocatable :: bndy
      integer, dimension(:,:), allocatable :: kbdy
      integer, dimension(:,:,:), allocatable :: hflxw,hflxe,hflxs,hflxn

      logical :: dorestart,dowriteout,dostat,doprclout,dotdwrite,doazimwrite,dohifrqwrite,dohailout

      ! eddy recycling:
      integer, dimension(:,:), allocatable :: recy_cap_we,recy_inj_we
      integer, dimension(:,:), allocatable :: recy_cap_sn,recy_inj_sn
      real, dimension(:,:), allocatable :: recy_wsp_w,recy_wsp_e,recy_tau_w,recy_tau_e
      real, dimension(:,:), allocatable :: recy_wsp_s,recy_wsp_n,recy_tau_s,recy_tau_n
      real, dimension(:,:,:,:), allocatable :: recywe
      real, dimension(:,:,:,:), allocatable :: recysn

      ! for large-scale nudging:
      real, dimension(:,:), allocatable :: lsnudge_u,lsnudge_v,lsnudge_th,lsnudge_qv

      ! arrays for time-averaging:
      integer, dimension(:), allocatable :: ntavg
      double precision, dimension(:), allocatable :: rtavg
      real, dimension(:,:,:,:,:), allocatable :: tavg
      real, dimension(:,:,:,:), allocatable :: timavg
      real, dimension(:,:,:,:), allocatable :: sfctavg
      real, dimension(:,:,:), allocatable :: sfctimavg

!-----

      integer*8 :: count,rate,maxr
      real*8 :: xtime,time_solve,time_solve0,sum
      real :: rtime,cfl_target,ks_target,kminit,khinit,dttmp
      real :: steptime1,steptime2
      integer :: i,j,k,k1,k2,n,nn,fnum,frec,confignum
      real :: tem0,tem
      logical :: getsfc,getpbl,update_sfc,startup,restarted,restart_prcl,reset,getrad,restart_hail
      logical :: dosolve,dorad,getdbz,getvt,getsed,getqdiags,doit,dotbud,doqbud,doubud,dovbud,dowbud,donudge,increment

      integer :: icrs,icenter,jcenter
      real :: xcenter,ycenter
      double precision :: domainlocx,domainlocy

      double precision :: adaptmovetim
      integer :: mvrec,nwritemv

      integer :: ntot
      double precision :: rtot,rdenom
      real(kind=qp) :: qtem
      logical :: printit
      logical :: spinup_t2p

#ifdef MPI
      integer :: reqs,rc,ii,jj,id,itmp1,itmp2,jtmp1,jtmp2,ierr
      integer, dimension(MPI_STATUS_SIZE) :: status
      real :: mp_total,minval,temx,temy,temni,temnj
      double precision :: tstart,tend
      integer, dimension(:), allocatable :: isum,jsum
#endif

!----------------------------------------------------------------------
!  Initialize MPI

      myid=0
      numprocs=1

#ifdef MPI
      call MPI_INIT( ierr )
      call MPI_COMM_RANK( MPI_COMM_WORLD, myid, ierr )
      call MPI_COMM_SIZE( MPI_COMM_WORLD, numprocs, ierr )
#endif

      nomp = 1
#ifdef OPENMP
      nomp = omp_get_max_threads()
#endif

!----------------------------------------------------------------------

    ! set a bunch of parameters to some default value 
    ! (note: many of these can/will be over-written by namelist values)

      nf=0
      nu=0
      nv=0
      nw=0
      nstep = 0
      nstep0 = 0
      mtime = 0.0d0
      srec=1
      sirec=1
      urec=1
      vrec=1
      wrec=1
      xrec=1
      yrec=1
      zrec=1
      nrec=1
      mrec=1
      prec=1
      hrec=1
      trecs=1
      trecw=1
      arecs=1
      arecw=1
      nrst=1
      nstatout=1
      nwrite = 1
      nwritet = 1
      nwritea = 1
      nwriteh = 1
      outfile=6
      outunits=1
      startup = .true.
      stopit = .false.
      stopcondition = 0
      smeps = 1.0e-30
      tsmall = 0.0001
      ! (should be same as qsmall in morrison scheme)
      qsmall = 1.0e-14
#ifdef DP
      smeps = 1.0e-60
      tsmall = 0.000000000001
      qsmall = 1.0e-24
#endif
      big_num = 0.0
      big_num = huge(big_num)
      ndcnst = 0.0
      nt_c = 0.0
      qmag = 0.01
      np3a = 1
      np3o = 1
      cflmax = 0.0
      ksmax = 0.0
      csound = 0.0
      cstar = 0.0
      ndt = 0
      adt = 0.0
      adtlast = 0.0
      acfl = 0.0
      reqc = 0
      reqk = 0
      mass1 = 0.0
      mass2 = 0.0
      avgsfcu = 0.0
      avgsfcv = 0.0
      avgsfcs = 0.0
      avgsfcsu = 0.0
      avgsfcsv = 0.0
      avgsfct = 0.0
      avgsfcq = 0.0
      avgsfcp = 0.0
      nudgeobc = 0
      alphobc = 60.0
      diagtim = 0.0
      domainlocx = 0.0
      domainlocy = 0.0
      top_cd = 0.0
      side_cd = 0.0
      brkup_db = 0.0
      solconrce = 0.0
      zenanglrce = 0.0
      trn_swirl = 0.0
      trn_wfrc = 0.0
      trn_lr = 0.0
      trn_lz = 0.0
      trn_xc = 0.0
      trn_yc = 0.0
      trn_zc = 0.0
      z1old = 0.0
      timeavg_num = 1
      timeavg_len = 0.0
      gamma_cp = 100.0
      gamma_dr = 18.0
      cfl_limit = 1.0
      dtfrac = 0.1
      adapt_dt = 0
      ramp_down = 1
      kmaxib = 1

      getsfc = .true.
      getpbl = .true.
      ramp_sgs = 0
      ramp_time = 0.0
      t2p_avg = 1
      dotimeavg = .false.
      restarted = .false.
      restart_prcl = .false.
      restart_hail = .false.
      restart_format = 1
      restart_filetype = 2
      restart_reset_frqtim = .true.
      interp_on_restart = .false.
      restart_time = 0.0
      run_time = -999.0
      output_z_max = -999.9
      hurr_rad = 1.0e30
      dolsw = .false.
      cmpr_output = .false.
      dohifrq = .false.
      iusetke = .false.
      do_adapt_move = .false.
      do_recycle = .false.
      cuparam = 0
      nssl_3moment = .false.
      nssl_density_on = .true.
      irandp = 0
      irandprst = 0
      init_pbl = 1
      init_sfcmodel = 1

      dorestart = .false.
      dowriteout = .false.
      dostat = .false.
      doprclout = .false.
      dohailout = .false.
      dotdwrite = .false.
      doazimwrite = .false.
      dohifrqwrite = .false.
      dotbud = .true.
      doqbud = .false.
      doubud = .false.
      dovbud = .false.
      dowbud = .false.
      getrad = .false.
      writeit = .true.

      confignum = 1

      noutcheck = .false.
      nconfig = 1

      adaptmovetim = 0.0
      nwritemv = 1
      mvrec = 1

      idoles = .false.
      idopbl = .false.
      use_avg_sfc = .false.
      hpbl_grdm = .false.
      spinup_t2p = .false.

      getdbz = .false.
      getvt = .false.
      getsed = .false.
      getqdiags = .false.

      cm1version = 'cm1r22.0'
      cm1rversion = 22.0

!----------------------------------------------------------------------
!  For arbitrary 3d output array:  
!
!    nout3d is the total number of output variables.
!
!    User must then "populate" the out3d array.  (That is, you must
!    fill in the out3d(i,j,k,1) and out3d(i,j,k,2) etc arrays with the
!    information you wish to write out.)
!
!    No other changes are necessary:  output file will contain the 
!    variables "out1", "out2", "out3", ... etc.

      nout3d =   0

!----------------------------------------------------------------------
!  For arbitrary 2d output array:  
!
!    nout2d is the total number of output variables.
!
!    User must then "populate" the out2d array.  (That is, you must
!    fill in the out2d(i,j,1) and out3d(i,j,2) etc arrays with the
!    information you wish to write out.)
!
!    No other changes are necessary:  output file will contain the 
!    variables "out2d1", "out2d2", "out2d3", ... etc.

      nout2d =   0

!----------------------------------------------------------------------

      if( myid.eq.0 )then
        print *
        print *,'|---------------------------------------------------------------|'
        print *,'|                                                               |'
        print *,'|          You are using CM1, Release 22.0 (cm1r22.0)           |'
        print *,'|                                                               |'
        print *,'|---------------------------------------------------------------|'

#ifdef MPI
        print *,'|                                                               |'
        print *,'|    compiled using distributed-memory parallelization (MPI)    |'
        print *,'|                                                               |'
        print *,'|---------------------------------------------------------------|'
#endif

#ifdef OPENMP
        print *,'|                                                               |'
        print *,'|     compiled using shared-memory parallelization (OpenMP)     |'
        print *,'|                                                               |'
        print *,'|---------------------------------------------------------------|'
#endif

        print *
      endif

!----------------------------------------------------------------------

      ! This initializes timer
      time_last = 0.0d0
      call system_clock(count=count,count_rate=rate,count_max=maxr)
      clock_rate = 1.0d0/float(rate)
      xtime = mytime()

      call set_time_to_zero()

!----------------------------------------------------------------------
!  Get domain dimensions, allocate some arrays, then call PARAM

      open(unit=20,file='namelist.input',form='formatted',status='old',    &
           access='sequential',err=8778)
      read(20,nml=param0)
      read(20,nml=param2)
      read(20,nml=param8)
      close(unit=20)

      ! note:  read remainder of namelist sections in param.F !

!----------------------------------------------------------------------

      IF( procfiles )THEN
#ifdef MPI
        outfile=10
#endif
        dowr = .true.
      ELSE
        dowr = .false.
      ENDIF

      IF( myid.eq.0 ) dowr = .true.

!----------------------------------------------------------------------
#ifdef MPI
  !......
  ! cm1r20.1:
  ! determine values for nodex and nodey:
    if(myid.eq.0)then
      minval = 1.0e30
      open(unit=10,file='proc.info',status='unknown')
      write(10,*)
    IF( nx.eq.1 .and. ny.eq.1 )THEN
      print *
      print *,'  nx = ',nx
      print *,'  ny = ',ny
      print *
      print *,'  Cannot use MPI when nx = ny = 1 '
      print *
      call stopcm1
    ELSEIF( nx.le.3 )THEN
      nodex = 1
      nodey = numprocs
    ELSEIF( ny.le.3 )THEN
      nodex = numprocs
      nodey = 1
    ELSE
      do n = numprocs , 1 , -1
        temy = float(n)
        temx = float(numprocs)/temy
        temni = float(nx)/temx
        temnj = float(ny)/temy
!!!        write(10,*) temx,temy,temx*temy,numprocs,temx/temy

        ! target: ni = 2 x nj
        tem0 = abs(temni/temnj-2.0)

        write(10,*) temx,temy,temni,temnj,temni/temnj
        ! only consider cases in which temx,temy are both integers:
        if( abs(float(int(temx))-temx).lt.smeps )then
          write(10,*) '      temx,temy   = ',temx,temy
          write(10,*) '      temni,temnj = ',temni,temnj
          if( tem0.lt.minval )then
            minval = tem0
            nodex = temx
            nodey = temy
            write(10,*) '        minval = ',minval
          endif
        endif
      enddo
      write(10,*)
      write(10,*) '  nodex,nodey = ',nodex,nodey,nodex*nodey,numprocs
      write(10,*)
    ENDIF
    endif

    ! just to be sure:
    nodex  = max(1,nodex)
    nodey  = max(1,nodey)
    ppnode = max(1,ppnode)

    call MPI_BCAST(nodex,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(nodey,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

  !......
#else
      ! serial (i.e. single-processor) run:
      nodex = 1
      nodey = 1
      ppnode = 1
#endif

      ni = nx / nodex
      nj = ny / nodey
      nk = nz
      nkp1 = nk+1

      nimax = ni
      njmax = nj

      myi=1
      myj=1

      myi1 =  1
      myi2 = nx
      myj1 =  1
      myj2 = ny

      if( .not. allocated( myi1p ) ) allocate( myi1p(0:numprocs-1) )
      myi1p = 1
      if( .not. allocated( myi2p ) ) allocate( myi2p(0:numprocs-1) )
      myi2p = nx
      if( .not. allocated( myj1p ) ) allocate( myj1p(0:numprocs-1) )
      myj1p = 1
      if( .not. allocated( myj2p ) ) allocate( myj2p(0:numprocs-1) )
      myj2p = ny

      if( .not. allocated( myi1pp ) ) allocate( myi1pp(0:numprocs-1) )
      myi1pp = 0
      if( .not. allocated( myi2pp ) ) allocate( myi2pp(0:numprocs-1) )
      myi2pp = 0
      if( .not. allocated( myj1pp ) ) allocate( myj1pp(0:numprocs-1) )
      myj1pp = 0
      if( .not. allocated( myj2pp ) ) allocate( myj2pp(0:numprocs-1) )
      myj2pp = 0

#ifdef MPI
!-----------------------------------------------
!  Determine number of grid points for each processor.
!
!  cm1r20:  each processor is no longer required
!           to have the same number of grid points

      myj = myid / nodex + 1
      myi = myid - (myj-1)*nodex  + 1

      if( .not. allocated( isum ) ) allocate( isum(nodex) )
      isum = 0
      if( .not. allocated( jsum ) ) allocate( jsum(nodey) )
      jsum = 0

      ! determine number of "i" points for each processor:

      ni = 0
      id = 0

      do i = 1,nx
        id = id+1
        isum(id) = isum(id)+1
        if( id.eq.myi ) ni = ni+1
        if( id.eq.nodex ) id = 0
      enddo

      ! determine number of "j" points for each processor:

      nj = 0
      id = 0

      do j = 1,ny
        id = id+1
        jsum(id) = jsum(id)+1
        if( id.eq.myj ) nj = nj+1
        if( id.eq.nodey ) id = 0
      enddo

!......

      ! first gridpoint for each processor (from the total grid):

      myi1 = 1
      do i=1,(myi-1)
        myi1 = myi1+isum(i)
      enddo
      myi2 = myi1+ni-1

      myj1 = 1
      do j=1,(myj-1)
        myj1 = myj1+jsum(j)
      enddo
      myj2 = myj1+nj-1

!!!      write(10,*) '  myid,myi1,myi2,myj1,myj2 = ',myid,myi1,myi2,myj1,myj

!......

      ! write processor info:

      n = 0

      if( myid.eq.0 )then
        write(10,*)
        write(10,*) '  n,i,j,ni,nj:'
        do j = 1,nodey
        do i = 1,nodex
          if( i.eq.1 .and. j.eq.1 )then
            write(10,*) n,i,j,ni,nj
          else
            call mpi_recv(ii   ,1,mpi_integer,n,1,MPI_COMM_WORLD,status,ierr)
            call mpi_recv(jj   ,1,mpi_integer,n,2,MPI_COMM_WORLD,status,ierr)
            write(10,*) n,i,j,ii,jj
          endif
          n = n+1
        enddo
        enddo
      else
        call mpi_send(ni  ,1,mpi_integer,0,1,MPI_COMM_WORLD,ierr)
        call mpi_send(nj  ,1,mpi_integer,0,2,MPI_COMM_WORLD,ierr)
      endif

!......

      ! write more processor info:

      n = 0

      if( myid.eq.0 )then
        write(10,*)
        write(10,*) '  n,myi1,myi2,myj1,myj2:'
        do j = 1,nodey
        do i = 1,nodex
          if( i.eq.1 .and. j.eq.1 )then
            write(10,*) n,myi1,myi2,myj1,myj2
            myi1p(0) = myi1
            myi2p(0) = myi2
            myj1p(0) = myj1
            myj2p(0) = myj2
          else
            call mpi_recv(itmp1,1,mpi_integer,n,3,MPI_COMM_WORLD,status,ierr)
            call mpi_recv(itmp2,1,mpi_integer,n,4,MPI_COMM_WORLD,status,ierr)
            call mpi_recv(jtmp1,1,mpi_integer,n,5,MPI_COMM_WORLD,status,ierr)
            call mpi_recv(jtmp2,1,mpi_integer,n,6,MPI_COMM_WORLD,status,ierr)
            write(10,*) n,itmp1,itmp2,jtmp1,jtmp2
            myi1p(n) = itmp1
            myi2p(n) = itmp2
            myj1p(n) = jtmp1
            myj2p(n) = jtmp2
          endif
          n = n+1
        enddo
        enddo
        write(10,*)
      else
        call mpi_send(myi1,1,mpi_integer,0,3,MPI_COMM_WORLD,ierr)
        call mpi_send(myi2,1,mpi_integer,0,4,MPI_COMM_WORLD,ierr)
        call mpi_send(myj1,1,mpi_integer,0,5,MPI_COMM_WORLD,ierr)
        call mpi_send(myj2,1,mpi_integer,0,6,MPI_COMM_WORLD,ierr)
      endif

      call MPI_BCAST(myi1p(0),numprocs,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(myi2p(0),numprocs,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(myj1p(0),numprocs,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(myj2p(0),numprocs,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

      nimax = 0
      do n=0,numprocs-1
        nimax = max( nimax , myi2p(n)-myi1p(n)+1 )
      enddo

      njmax = 0
      do n=0,numprocs-1
        njmax = max( njmax , myj2p(n)-myj1p(n)+1 )
      enddo

!......

      deallocate( isum )
      deallocate( jsum )

      close(unit=10)

!!!      call MPI_BARRIER (MPI_COMM_WORLD,ierr)

!-----------------------------------------------
#endif


      call wenocheck

      ! number of 'ghost' points in the horizontal directions:
      ngxy  = 2
      if( hadvordrs.eq.3 .or. hadvordrv.eq.3 .or. hadvordrs.eq.4 .or. hadvordrv.eq.4 )   ngxy = 2
      if( hadvordrs.eq.5 .or. hadvordrv.eq.5 .or. hadvordrs.eq.6 .or. hadvordrv.eq.6 )   ngxy = 3
      if( hadvordrs.eq.7 .or. hadvordrv.eq.7 .or. hadvordrs.eq.8 .or. hadvordrv.eq.8 )   ngxy = 4
      if( hadvordrs.eq.9 .or. hadvordrv.eq.9 .or. hadvordrs.eq.10 .or. hadvordrv.eq.10 ) ngxy = 5

      if( advwenos.ge.1 .or. advwenov.ge.1 )then
        if( weno_order.eq.3 ) ngxy = max(ngxy,2)
        if( weno_order.eq.5 ) ngxy = max(ngxy,3)
        if( weno_order.eq.7 ) ngxy = max(ngxy,4)
        if( weno_order.eq.9 ) ngxy = max(ngxy,5)
      endif

      ! number of 'ghost' points in the vertical direction:
      ngz   = 1

!---------------------------------------------------------------------
!      For ZVD:
!      ngz   = 3
!      IF( ngz.eq.3 )THEN
!        kb =  1 - ngz
!        ke = nk + ngz
!      ENDIF
!---------------------------------------------------------------------

      ib =  1 - ngxy
      ie = ni + ngxy
      jb =  1 - ngxy
      je = nj + ngxy
      kb =  1 - ngz
      ke = nk + ngz

      if( .not. allocated( xh ) ) allocate(    xh(ib:ie) )
      xh = 0.0
      if( .not. allocated( rxh ) ) allocate(   rxh(ib:ie) )
      rxh = 0.0
      if( .not. allocated( arh1 ) ) allocate(  arh1(ib:ie) )
      arh1 = 0.0
      if( .not. allocated( arh2 ) ) allocate(  arh2(ib:ie) )
      arh2 = 0.0
      if( .not. allocated( uh ) ) allocate(    uh(ib:ie) )
      uh = 0.0
      if( .not. allocated( ruh ) ) allocate(   ruh(ib:ie) )
      ruh = 0.0
      if( .not. allocated( xf ) ) allocate(    xf(ib:ie+1) )
      xf = 0.0
      if( .not. allocated( rxf ) ) allocate(   rxf(ib:ie+1) )
      rxf = 0.0
      if( .not. allocated( arf1 ) ) allocate(  arf1(ib:ie+1) )
      arf1 = 0.0
      if( .not. allocated( arf2 ) ) allocate(  arf2(ib:ie+1) )
      arf2 = 0.0
      if( .not. allocated( uf ) ) allocate(    uf(ib:ie+1) )
      uf = 0.0
      if( .not. allocated( ruf ) ) allocate(   ruf(ib:ie+1) )
      ruf = 0.0
      if( .not. allocated( yh ) ) allocate(    yh(jb:je) )
      yh = 0.0
      if( .not. allocated( vh ) ) allocate(    vh(jb:je) )
      vh = 0.0
      if( .not. allocated( rvh ) ) allocate(   rvh(jb:je) )
      rvh = 0.0
      if( .not. allocated( yf ) ) allocate(    yf(jb:je+1) )
      yf = 0.0
      if( .not. allocated( vf ) ) allocate(    vf(jb:je+1) )
      vf = 0.0
      if( .not. allocated( rvf ) ) allocate(   rvf(jb:je+1) )
      rvf = 0.0
      if( .not. allocated( xfref ) ) allocate( xfref(1-ngxy:nx+ngxy+1) )
      xfref = 0.0
      if( .not. allocated( xhref ) ) allocate( xhref(1-ngxy:nx+ngxy+1) )
      xhref = 0.0
      if( .not. allocated( yfref ) ) allocate( yfref(1-ngxy:ny+ngxy+1) )
      yfref = 0.0
      if( .not. allocated( yhref ) ) allocate( yhref(1-ngxy:ny+ngxy+1) )
      yhref = 0.0
      if( .not. allocated( dumk1 ) ) allocate( dumk1(kb:ke) )
      dumk1 = 0.0
      if( .not. allocated( dumk2 ) ) allocate( dumk2(kb:ke) )
      dumk2 = 0.0
      if( .not. allocated( dumk3 ) ) allocate( dumk3(nk) )
      dumk3 = 0.0
      if( .not. allocated( dumk4 ) ) allocate( dumk4(nk) )
      dumk4 = 0.0
      if( .not. allocated( dumk5 ) ) allocate( dumk5(nk) )
      dumk5 = 0.0
      if( .not. allocated( rds ) ) allocate(   rds(kb:ke) )
      rds = 0.0
      if( .not. allocated( sigma ) ) allocate( sigma(kb:ke) )
      sigma = 0.0
      if( .not. allocated( rdsf ) ) allocate(   rdsf(kb:ke+1) )
      rdsf = 0.0
      if( .not. allocated( sigmaf ) ) allocate( sigmaf(kb:ke+1) )
      sigmaf = 0.0
      if( .not. allocated( tauh ) ) allocate(  tauh(ib:ie,jb:je,kb:ke) )
      tauh = 0.0
      if( .not. allocated( taus ) ) allocate(  taus(ib:ie,jb:je,kb:ke) )
      taus = 0.0
      if( .not. allocated( zh ) ) allocate(    zh(ib:ie,jb:je,kb:ke) )
      zh = 0.0
      if( .not. allocated( mh ) ) allocate(    mh(ib:ie,jb:je,kb:ke) )
      mh = 0.0
      if( .not. allocated( rmh ) ) allocate(   rmh(ib:ie,jb:je,kb:ke) )
      rmh = 0.0
      if( .not. allocated( c1 ) ) allocate(    c1(ib:ie,jb:je,kb:ke) )
      c1 = 0.0
      if( .not. allocated( c2 ) ) allocate(    c2(ib:ie,jb:je,kb:ke) )
      c2 = 0.0
      if( .not. allocated( tauf ) ) allocate(  tauf(ib:ie,jb:je,kb:ke+1) )
      tauf = 0.0
      if( .not. allocated( mf ) ) allocate(    mf(ib:ie,jb:je,kb:ke+1) )
      mf = 0.0
      if( .not. allocated( rmf ) ) allocate(   rmf(ib:ie,jb:je,kb:ke+1) )
      rmf = 0.0
      if( .not. allocated( f2d ) ) allocate(    f2d(ib:ie,jb:je) )
      f2d = 0.0
      if( .not. allocated( psmth ) ) allocate(  psmth(ib:ie,jb:je) )
      psmth = 0.0
      if( .not. allocated( prate ) ) allocate(  prate(ib:ie,jb:je) )
      prate = 0.0
      if( .not. allocated( pratecv ) ) allocate(  pratecv(ib:ie,jb:je) )
      pratecv = 0.0

      if( .not. allocated( gamk ) ) allocate( gamk(kb:ke) )
      gamk = 1.0
      if( .not. allocated( gamwall ) ) allocate( gamwall(kb:ke) )
      gamwall = 0.0

      if(terrain_flag)then
        itb=ib
        ite=ie
        jtb=jb
        jte=je
        ktb=kb
        kte=ke
      else
        itb=1
        ite=1
        jtb=1
        jte=1
        ktb=1
        kte=1
      endif

      if( .not. allocated( zs ) ) allocate(   zs(ib:ie,jb:je) )
      zs = 0.0
      if( .not. allocated( gz ) ) allocate(   gz(itb:ite,jtb:jte) )
      gz = 0.0
      if( .not. allocated( rgz ) ) allocate(  rgz(itb:ite,jtb:jte) )
      rgz = 0.0
      if( .not. allocated( gzu ) ) allocate(  gzu(itb:ite,jtb:jte) )
      gzu = 0.0
      if( .not. allocated( rgzu ) ) allocate( rgzu(itb:ite,jtb:jte) )
      rgzu = 0.0
      if( .not. allocated( gzv ) ) allocate(  gzv(itb:ite,jtb:jte) )
      gzv = 0.0
      if( .not. allocated( rgzv ) ) allocate( rgzv(itb:ite,jtb:jte) )
      rgzv = 0.0
      if( .not. allocated( dzdx ) ) allocate( dzdx(itb:ite,jtb:jte) )
      dzdx = 0.0
      if( .not. allocated( dzdy ) ) allocate( dzdy(itb:ite,jtb:jte) )
      dzdy = 0.0
      if( .not. allocated( gx ) ) allocate(   gx(itb:ite,jtb:jte,ktb:kte) )
      gx = 0.0
      if( .not. allocated( gxu ) ) allocate(  gxu(itb:ite,jtb:jte,ktb:kte) )
      gxu = 0.0
      if( .not. allocated( gy ) ) allocate(   gy(itb:ite,jtb:jte,ktb:kte) )
      gy = 0.0
      if( .not. allocated( gyv ) ) allocate(  gyv(itb:ite,jtb:jte,ktb:kte) )
      gyv = 0.0
      if( .not. allocated( zf ) ) allocate(   zf(ib:ie,jb:je,kb:ke+1) )
      zf = 0.0

      if( .not. allocated( dtu ) ) allocate( dtu(ib:ie,jb:je) )
      dtu = 1.0
      if( .not. allocated( dtu0 ) ) allocate( dtu0(ib:ie,jb:je) )
      dtu0 = 1.0
      if( .not. allocated( dtv ) ) allocate( dtv(ib:ie,jb:je) )
      dtv = 1.0
      if( .not. allocated( dtv0 ) ) allocate( dtv0(ib:ie,jb:je) )
      dtv0 = 1.0

!----------------------------------------------------------------------
!  New for cm1r19.6  (not thoroughly tested)

      dohifrq = .false.     ! do high-frequency output?

      hifrqfrq =    60.0    ! frequency (seconds)

!----------------------------------------------------------------------

!------
! allocate the MPI arrays

#ifdef MPI
      imp = max(1,ni)
      jmp = max(1,nj)
      kmp = max(2,nk)
      kmt = max(2,nk+1)
      rmp = 8
      cmp = ngxy
#else
      imp = 1
      jmp = 1
      kmp = 2
      kmt = 2
      rmp = 1
      cmp = 1
#endif

      if( .not. allocated( reqs_u ) ) allocate( reqs_u(rmp) )
      reqs_u = 0
      if( .not. allocated( reqs_v ) ) allocate( reqs_v(rmp) )
      reqs_v = 0
      if( .not. allocated( reqs_w ) ) allocate( reqs_w(rmp) )
      reqs_w = 0
      if( .not. allocated( reqs_s ) ) allocate( reqs_s(rmp) )
      reqs_s = 0
      if( .not. allocated( reqs_p ) ) allocate( reqs_p(rmp) )
      reqs_p = 0
      if( .not. allocated( reqs_p2 ) ) allocate( reqs_p2(rmp) )
      reqs_p2 = 0
      if( .not. allocated( reqs_p3 ) ) allocate( reqs_p3(rmp) )
      reqs_p3 = 0
      if( .not. allocated( reqs_x ) ) allocate( reqs_x(rmp) )
      reqs_x = 0
      if( .not. allocated( reqs_y ) ) allocate( reqs_y(rmp) )
      reqs_y = 0
      if( .not. allocated( reqs_z ) ) allocate( reqs_z(rmp) )
      reqs_z = 0
      if( .not. allocated( reqs_tk ) ) allocate( reqs_tk(rmp) )
      reqs_tk = 0
      if( .not. allocated( reqs_tk2 ) ) allocate( reqs_tk2(rmp) )
      reqs_tk2 = 0

      if( .not. allocated( nw1 ) ) allocate( nw1(kmt) )
      nw1 = 0.0
      if( .not. allocated( nw2 ) ) allocate( nw2(kmt) )
      nw2 = 0.0
      if( .not. allocated( ne1 ) ) allocate( ne1(kmt) )
      ne1 = 0.0
      if( .not. allocated( ne2 ) ) allocate( ne2(kmt) )
      ne2 = 0.0
      if( .not. allocated( sw1 ) ) allocate( sw1(kmt) )
      sw1 = 0.0
      if( .not. allocated( sw2 ) ) allocate( sw2(kmt) )
      sw2 = 0.0
      if( .not. allocated( se1 ) ) allocate( se1(kmt) )
      se1 = 0.0
      if( .not. allocated( se2 ) ) allocate( se2(kmt) )
      se2 = 0.0

      if( .not. allocated( n3w1 ) ) allocate( n3w1(cmp,cmp,kmt+1) )
      n3w1 = 0.0
      if( .not. allocated( n3w2 ) ) allocate( n3w2(cmp,cmp,kmt+1) )
      n3w2 = 0.0
      if( .not. allocated( n3e1 ) ) allocate( n3e1(cmp,cmp,kmt+1) )
      n3e1 = 0.0
      if( .not. allocated( n3e2 ) ) allocate( n3e2(cmp,cmp,kmt+1) )
      n3e2 = 0.0
      if( .not. allocated( s3w1 ) ) allocate( s3w1(cmp,cmp,kmt+1) )
      s3w1 = 0.0
      if( .not. allocated( s3w2 ) ) allocate( s3w2(cmp,cmp,kmt+1) )
      s3w2 = 0.0
      if( .not. allocated( s3e1 ) ) allocate( s3e1(cmp,cmp,kmt+1) )
      s3e1 = 0.0
      if( .not. allocated( s3e2 ) ) allocate( s3e2(cmp,cmp,kmt+1) )
      s3e2 = 0.0

      if( .not. allocated( ww1 ) ) allocate( ww1(jmp,kmp-1) )
      ww1 = 0.0
      if( .not. allocated( ww2 ) ) allocate( ww2(jmp,kmp-1) )
      ww2 = 0.0
      if( .not. allocated( we1 ) ) allocate( we1(jmp,kmp-1) )
      we1 = 0.0
      if( .not. allocated( we2 ) ) allocate( we2(jmp,kmp-1) )
      we2 = 0.0
      if( .not. allocated( ws1 ) ) allocate( ws1(imp,kmp-1) )
      ws1 = 0.0
      if( .not. allocated( ws2 ) ) allocate( ws2(imp,kmp-1) )
      ws2 = 0.0
      if( .not. allocated( wn1 ) ) allocate( wn1(imp,kmp-1) )
      wn1 = 0.0
      if( .not. allocated( wn2 ) ) allocate( wn2(imp,kmp-1) )
      wn2 = 0.0

      if( .not. allocated( pw1 ) ) allocate( pw1(jmp,kmp) )
      pw1 = 0.0
      if( .not. allocated( pw2 ) ) allocate( pw2(jmp,kmp) )
      pw2 = 0.0
      if( .not. allocated( pe1 ) ) allocate( pe1(jmp,kmp) )
      pe1 = 0.0
      if( .not. allocated( pe2 ) ) allocate( pe2(jmp,kmp) )
      pe2 = 0.0
      if( .not. allocated( ps1 ) ) allocate( ps1(imp,kmp) )
      ps1 = 0.0
      if( .not. allocated( ps2 ) ) allocate( ps2(imp,kmp) )
      ps2 = 0.0
      if( .not. allocated( pn1 ) ) allocate( pn1(imp,kmp) )
      pn1 = 0.0
      if( .not. allocated( pn2 ) ) allocate( pn2(imp,kmp) )
      pn2 = 0.0

      if( .not. allocated( p2w1 ) ) allocate( p2w1(jmp,kmp) )
      p2w1 = 0.0
      if( .not. allocated( p2w2 ) ) allocate( p2w2(jmp,kmp) )
      p2w2 = 0.0
      if( .not. allocated( p2e1 ) ) allocate( p2e1(jmp,kmp) )
      p2e1 = 0.0
      if( .not. allocated( p2e2 ) ) allocate( p2e2(jmp,kmp) )
      p2e2 = 0.0
      if( .not. allocated( p2s1 ) ) allocate( p2s1(imp,kmp) )
      p2s1 = 0.0
      if( .not. allocated( p2s2 ) ) allocate( p2s2(imp,kmp) )
      p2s2 = 0.0
      if( .not. allocated( p2n1 ) ) allocate( p2n1(imp,kmp) )
      p2n1 = 0.0
      if( .not. allocated( p2n2 ) ) allocate( p2n2(imp,kmp) )
      p2n2 = 0.0

      if( .not. allocated( p3w1 ) ) allocate( p3w1(jmp,kmp) )
      p3w1 = 0.0
      if( .not. allocated( p3w2 ) ) allocate( p3w2(jmp,kmp) )
      p3w2 = 0.0
      if( .not. allocated( p3e1 ) ) allocate( p3e1(jmp,kmp) )
      p3e1 = 0.0
      if( .not. allocated( p3e2 ) ) allocate( p3e2(jmp,kmp) )
      p3e2 = 0.0
      if( .not. allocated( p3s1 ) ) allocate( p3s1(imp,kmp) )
      p3s1 = 0.0
      if( .not. allocated( p3s2 ) ) allocate( p3s2(imp,kmp) )
      p3s2 = 0.0
      if( .not. allocated( p3n1 ) ) allocate( p3n1(imp,kmp) )
      p3n1 = 0.0
      if( .not. allocated( p3n2 ) ) allocate( p3n2(imp,kmp) )
      p3n2 = 0.0

      if( .not. allocated( vw1 ) ) allocate( vw1(jmp,kmp) )
      vw1 = 0.0
      if( .not. allocated( vw2 ) ) allocate( vw2(jmp,kmp) )
      vw2 = 0.0
      if( .not. allocated( ve1 ) ) allocate( ve1(jmp,kmp) )
      ve1 = 0.0
      if( .not. allocated( ve2 ) ) allocate( ve2(jmp,kmp) )
      ve2 = 0.0
      if( .not. allocated( vs1 ) ) allocate( vs1(imp,kmp) )
      vs1 = 0.0
      if( .not. allocated( vs2 ) ) allocate( vs2(imp,kmp) )
      vs2 = 0.0
      if( .not. allocated( vn1 ) ) allocate( vn1(imp,kmp) )
      vn1 = 0.0
      if( .not. allocated( vn2 ) ) allocate( vn2(imp,kmp) )
      vn2 = 0.0

      if( .not. allocated( zw1 ) ) allocate( zw1(jmp,kmp) )
      zw1 = 0.0
      if( .not. allocated( zw2 ) ) allocate( zw2(jmp,kmp) )
      zw2 = 0.0
      if( .not. allocated( ze1 ) ) allocate( ze1(jmp,kmp) )
      ze1 = 0.0
      if( .not. allocated( ze2 ) ) allocate( ze2(jmp,kmp) )
      ze2 = 0.0
      if( .not. allocated( zs1 ) ) allocate( zs1(imp,kmp) )
      zs1 = 0.0
      if( .not. allocated( zs2 ) ) allocate( zs2(imp,kmp) )
      zs2 = 0.0
      if( .not. allocated( zn1 ) ) allocate( zn1(imp,kmp) )
      zn1 = 0.0
      if( .not. allocated( zn2 ) ) allocate( zn2(imp,kmp) )
      zn2 = 0.0

      if( .not. allocated( uw31 ) ) allocate( uw31(cmp,jmp,kmp) )
      uw31 = 0.0
      if( .not. allocated( uw32 ) ) allocate( uw32(cmp,jmp,kmp) )
      uw32 = 0.0
      if( .not. allocated( ue31 ) ) allocate( ue31(cmp,jmp,kmp) )
      ue31 = 0.0
      if( .not. allocated( ue32 ) ) allocate( ue32(cmp,jmp,kmp) )
      ue32 = 0.0
      if( .not. allocated( us31 ) ) allocate( us31(imp+1,cmp,kmp) )
      us31 = 0.0
      if( .not. allocated( us32 ) ) allocate( us32(imp+1,cmp,kmp) )
      us32 = 0.0
      if( .not. allocated( un31 ) ) allocate( un31(imp+1,cmp,kmp) )
      un31 = 0.0
      if( .not. allocated( un32 ) ) allocate( un32(imp+1,cmp,kmp) )
      un32 = 0.0

      if( .not. allocated( vw31 ) ) allocate( vw31(cmp,jmp+1,kmp) )
      vw31 = 0.0
      if( .not. allocated( vw32 ) ) allocate( vw32(cmp,jmp+1,kmp) )
      vw32 = 0.0
      if( .not. allocated( ve31 ) ) allocate( ve31(cmp,jmp+1,kmp) )
      ve31 = 0.0
      if( .not. allocated( ve32 ) ) allocate( ve32(cmp,jmp+1,kmp) )
      ve32 = 0.0
      if( .not. allocated( vs31 ) ) allocate( vs31(imp,cmp,kmp) )
      vs31 = 0.0
      if( .not. allocated( vs32 ) ) allocate( vs32(imp,cmp,kmp) )
      vs32 = 0.0
      if( .not. allocated( vn31 ) ) allocate( vn31(imp,cmp,kmp) )
      vn31 = 0.0
      if( .not. allocated( vn32 ) ) allocate( vn32(imp,cmp,kmp) )
      vn32 = 0.0

      if( .not. allocated( ww31 ) ) allocate( ww31(cmp,jmp,kmp-1) )
      ww31 = 0.0
      if( .not. allocated( ww32 ) ) allocate( ww32(cmp,jmp,kmp-1) )
      ww32 = 0.0
      if( .not. allocated( we31 ) ) allocate( we31(cmp,jmp,kmp-1) )
      we31 = 0.0
      if( .not. allocated( we32 ) ) allocate( we32(cmp,jmp,kmp-1) )
      we32 = 0.0
      if( .not. allocated( ws31 ) ) allocate( ws31(imp,cmp,kmp-1) )
      ws31 = 0.0
      if( .not. allocated( ws32 ) ) allocate( ws32(imp,cmp,kmp-1) )
      ws32 = 0.0
      if( .not. allocated( wn31 ) ) allocate( wn31(imp,cmp,kmp-1) )
      wn31 = 0.0
      if( .not. allocated( wn32 ) ) allocate( wn32(imp,cmp,kmp-1) )
      wn32 = 0.0

      if( .not. allocated( sw31 ) ) allocate( sw31(cmp,jmp,kmp) )
      sw31 = 0.0
      if( .not. allocated( sw32 ) ) allocate( sw32(cmp,jmp,kmp) )
      sw32 = 0.0
      if( .not. allocated( se31 ) ) allocate( se31(cmp,jmp,kmp) )
      se31 = 0.0
      if( .not. allocated( se32 ) ) allocate( se32(cmp,jmp,kmp) )
      se32 = 0.0
      if( .not. allocated( ss31 ) ) allocate( ss31(imp,cmp,kmp) )
      ss31 = 0.0
      if( .not. allocated( ss32 ) ) allocate( ss32(imp,cmp,kmp) )
      ss32 = 0.0
      if( .not. allocated( sn31 ) ) allocate( sn31(imp,cmp,kmp) )
      sn31 = 0.0
      if( .not. allocated( sn32 ) ) allocate( sn32(imp,cmp,kmp) )
      sn32 = 0.0

      if( .not. allocated( rw31 ) ) allocate( rw31(cmp,jmp,kmp) )
      rw31 = 0.0
      if( .not. allocated( rw32 ) ) allocate( rw32(cmp,jmp,kmp) )
      rw32 = 0.0
      if( .not. allocated( re31 ) ) allocate( re31(cmp,jmp,kmp) )
      re31 = 0.0
      if( .not. allocated( re32 ) ) allocate( re32(cmp,jmp,kmp) )
      re32 = 0.0
      if( .not. allocated( rs31 ) ) allocate( rs31(imp,cmp,kmp) )
      rs31 = 0.0
      if( .not. allocated( rs32 ) ) allocate( rs32(imp,cmp,kmp) )
      rs32 = 0.0
      if( .not. allocated( rn31 ) ) allocate( rn31(imp,cmp,kmp) )
      rn31 = 0.0
      if( .not. allocated( rn32 ) ) allocate( rn32(imp,cmp,kmp) )
      rn32 = 0.0

      if( .not. allocated( tkw1 ) ) allocate( tkw1(cmp,jmp,kmt) )
      tkw1 = 0.0
      if( .not. allocated( tkw2 ) ) allocate( tkw2(cmp,jmp,kmt) )
      tkw2 = 0.0
      if( .not. allocated( tke1 ) ) allocate( tke1(cmp,jmp,kmt) )
      tke1 = 0.0
      if( .not. allocated( tke2 ) ) allocate( tke2(cmp,jmp,kmt) )
      tke2 = 0.0
      if( .not. allocated( tks1 ) ) allocate( tks1(imp,cmp,kmt) )
      tks1 = 0.0
      if( .not. allocated( tks2 ) ) allocate( tks2(imp,cmp,kmt) )
      tks2 = 0.0
      if( .not. allocated( tkn1 ) ) allocate( tkn1(imp,cmp,kmt) )
      tkn1 = 0.0
      if( .not. allocated( tkn2 ) ) allocate( tkn2(imp,cmp,kmt) )
      tkn2 = 0.0

      if( .not. allocated( kw1 ) ) allocate( kw1(jmp,kmt,4) )
      kw1 = 0.0
      if( .not. allocated( kw2 ) ) allocate( kw2(jmp,kmt,4) )
      kw2 = 0.0
      if( .not. allocated( ke1 ) ) allocate( ke1(jmp,kmt,4) )
      ke1 = 0.0
      if( .not. allocated( ke2 ) ) allocate( ke2(jmp,kmt,4) )
      ke2 = 0.0
      if( .not. allocated( ks1 ) ) allocate( ks1(imp,kmt,4) )
      ks1 = 0.0
      if( .not. allocated( ks2 ) ) allocate( ks2(imp,kmt,4) )
      ks2 = 0.0
      if( .not. allocated( kn1 ) ) allocate( kn1(imp,kmt,4) )
      kn1 = 0.0
      if( .not. allocated( kn2 ) ) allocate( kn2(imp,kmt,4) )
      kn2 = 0.0

      if( .not. allocated( npxa ) ) allocate( npxa(0:numprocs-1) )
      npxa = 0
      if( .not. allocated( npya ) ) allocate( npya(0:numprocs-1) )
      npya = 0
      if( .not. allocated( npza ) ) allocate( npza(0:numprocs-1) )
      npza = 0
      if( .not. allocated( npxb ) ) allocate( npxb(0:numprocs-1) )
      npxb = 0
      if( .not. allocated( npyb ) ) allocate( npyb(0:numprocs-1) )
      npyb = 0
      if( .not. allocated( npzb ) ) allocate( npzb(0:numprocs-1) )
      npzb = 0

      if( .not. allocated( pmyi1a ) ) allocate( pmyi1a(0:numprocs-1) )
      pmyi1a = 0
      if( .not. allocated( pmyi2a ) ) allocate( pmyi2a(0:numprocs-1) )
      pmyi2a = 0
      if( .not. allocated( pmyj1a ) ) allocate( pmyj1a(0:numprocs-1) )
      pmyj1a = 0
      if( .not. allocated( pmyj2a ) ) allocate( pmyj2a(0:numprocs-1) )
      pmyj2a = 0
      if( .not. allocated( pmyk1a ) ) allocate( pmyk1a(0:numprocs-1) )
      pmyk1a = 0
      if( .not. allocated( pmyk2a ) ) allocate( pmyk2a(0:numprocs-1) )
      pmyk2a = 0

      if( .not. allocated( pmyi1b ) ) allocate( pmyi1b(0:numprocs-1) )
      pmyi1b = 0
      if( .not. allocated( pmyi2b ) ) allocate( pmyi2b(0:numprocs-1) )
      pmyi2b = 0
      if( .not. allocated( pmyj1b ) ) allocate( pmyj1b(0:numprocs-1) )
      pmyj1b = 0
      if( .not. allocated( pmyj2b ) ) allocate( pmyj2b(0:numprocs-1) )
      pmyj2b = 0
      if( .not. allocated( pmyk1b ) ) allocate( pmyk1b(0:numprocs-1) )
      pmyk1b = 0
      if( .not. allocated( pmyk2b ) ) allocate( pmyk2b(0:numprocs-1) )
      pmyk2b = 0

      if( .not. allocated( npcomm ) ) allocate( npcomm(8) )
      npcomm = 0
      if( .not. allocated( nplist ) ) allocate( nplist(8,numprocs) )
      nplist = -1
      if( .not. allocated( npi1 ) ) allocate( npi1(8,numprocs) )
      npi1 =  1e9
      if( .not. allocated( npi2 ) ) allocate( npi2(8,numprocs) )
      npi2 = -1e9
      if( .not. allocated( npj1 ) ) allocate( npj1(8,numprocs) )
      npj1 =  1e9
      if( .not. allocated( npj2 ) ) allocate( npj2(8,numprocs) )
      npj2 = -1e9
      if( .not. allocated( npk1 ) ) allocate( npk1(8,numprocs) )
      npk1 =  1e9
      if( .not. allocated( npk2 ) ) allocate( npk2(8,numprocs) )
      npk2 = -1e9

#ifdef MPI
          call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif
      time_misc = mytime()
      time_param = 0.0

      call       param(dt,dtlast,stattim,taptim,rsttim,radtim,prcltim,  &
                       hailtrajtim,                                     &
                       cloudvar,rhovar,qmag,qname,qunit,budname,        &
                       xh,rxh,arh1,arh2,uh,ruh,xf,rxf,arf1,arf2,uf,ruf, &
                       yh,vh,rvh,yf,vf,rvf,xfref,xhref,yfref,yhref,     &
                       rds,sigma,rdsf,sigmaf,tauh,taus,                 &
                       zh,mh,rmh,c1 ,c2 ,tauf,zf,mf,rmf,f2d,dtu0,dtv0,  &
                       gamk,gamwall,                                    &
                       zs,gz,rgz,gzu,rgzu,gzv,rgzv,dzdx,dzdy,gx,gxu,gy,gyv,  &
                       reqs_u,reqs_v,reqs_s,reqs_p,                     &
                       nw1,nw2,ne1,ne2,sw1,sw2,se1,se2,                 &
                       n3w1,n3w2,n3e1,n3e2,s3w1,s3w2,s3e1,s3e2,         &
                       sw31,sw32,se31,se32,ss31,ss32,sn31,sn32,         &
                       uw31,uw32,ue31,ue32,us31,us32,un31,un32,         &
                       vw31,vw32,ve31,ve32,vs31,vs32,vn31,vn32,         &
                       ww31,ww32,we31,we32,ws31,ws32,wn31,wn32,         &
                       npxa,npya,npza,npxb,npyb,npzb,                   &
                       myi1p,myi2p,myj1p,myj2p,                         &
                       myi1pp,myi2pp,myj1pp,myj2pp,                     &
                       pmyi1a,pmyi2a,pmyj1a,pmyj2a,pmyk1a,pmyk2a,       &
                       pmyi1b,pmyi2b,pmyj1b,pmyj2b,pmyk1b,pmyk2b,       &
                       npcomm,nplist,npi1,npi2,npj1,npj2,npk1,npk2)

#ifdef MPI
          call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif
        time_param = time_param + mytime()
        if(dowr) write(outfile,*)
        if(dowr) write(outfile,*) '  time spent in subroutine param = ',time_param,' seconds '
        if(dowr) write(outfile,*)

!----------------------------------------------------------------------

      dbldt = dble(dt)

#ifdef MPI
          call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif
          if(myid.eq.0) print *
          if(myid.eq.0) print *,'  allocating MPI comm arrays ... '
#ifdef MPI
          call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif

      if( .not. allocated( reqs_q ) ) allocate( reqs_q(rmp,numq) )
      reqs_q = 0
      if( .not. allocated( reqs_t ) ) allocate( reqs_t(rmp,npt) )
      reqs_t = 0

      if( .not. allocated( qw31 ) ) allocate( qw31(cmp,jmp,kmp,numq) )
      qw31 = 0.0
      if( .not. allocated( qw32 ) ) allocate( qw32(cmp,jmp,kmp,numq) )
      qw32 = 0.0
      if( .not. allocated( qe31 ) ) allocate( qe31(cmp,jmp,kmp,numq) )
      qe31 = 0.0
      if( .not. allocated( qe32 ) ) allocate( qe32(cmp,jmp,kmp,numq) )
      qe32 = 0.0
      if( .not. allocated( qs31 ) ) allocate( qs31(imp,cmp,kmp,numq) )
      qs31 = 0.0
      if( .not. allocated( qs32 ) ) allocate( qs32(imp,cmp,kmp,numq) )
      qs32 = 0.0
      if( .not. allocated( qn31 ) ) allocate( qn31(imp,cmp,kmp,numq) )
      qn31 = 0.0
      if( .not. allocated( qn32 ) ) allocate( qn32(imp,cmp,kmp,numq) )
      qn32 = 0.0

      if( .not. allocated( tw1 ) ) allocate( tw1(cmp,jmp,kmp,npt) )
      tw1 = 0.0
      if( .not. allocated( tw2 ) ) allocate( tw2(cmp,jmp,kmp,npt) )
      tw2 = 0.0
      if( .not. allocated( te1 ) ) allocate( te1(cmp,jmp,kmp,npt) )
      te1 = 0.0
      if( .not. allocated( te2 ) ) allocate( te2(cmp,jmp,kmp,npt) )
      te2 = 0.0
      if( .not. allocated( ts1 ) ) allocate( ts1(imp,cmp,kmp,npt) )
      ts1 = 0.0
      if( .not. allocated( ts2 ) ) allocate( ts2(imp,cmp,kmp,npt) )
      ts2 = 0.0
      if( .not. allocated( tn1 ) ) allocate( tn1(imp,cmp,kmp,npt) )
      tn1 = 0.0
      if( .not. allocated( tn2 ) ) allocate( tn2(imp,cmp,kmp,npt) )
      tn2 = 0.0


      ! 1d arrays for test cases:
      if( .not. allocated( wprof ) ) allocate( wprof(kb:ke) )
      wprof = 0.0
      if( .not. allocated( ufrc ) ) allocate( ufrc(kb:ke) )
      ufrc = 0.0
      if( .not. allocated( vfrc ) ) allocate( vfrc(kb:ke) )
      vfrc = 0.0
      if( .not. allocated( thfrc ) ) allocate( thfrc(kb:ke) )
      thfrc = 0.0
      if( .not. allocated( qvfrc ) ) allocate( qvfrc(kb:ke) )
      qvfrc = 0.0
      if( .not. allocated( ug ) ) allocate( ug(kb:ke) )
      ug = 0.0
      if( .not. allocated( vg ) ) allocate( vg(kb:ke) )
      vg = 0.0
      if( .not. allocated( dvdr ) ) allocate( dvdr(kb:ke) )
      dvdr = 0.0

      if( .not. allocated( uavg ) ) allocate( uavg(kb:ke) )
      uavg = 0.0
      if( .not. allocated( vavg ) ) allocate( vavg(kb:ke) )
      vavg = 0.0
      if( .not. allocated( savg ) ) allocate( savg(kb:ke) )
      savg = 0.0
      if( .not. allocated( thavg ) ) allocate( thavg(kb:ke) )
      thavg = 0.0
      if( .not. allocated( pavg ) ) allocate( pavg(kb:ke) )
      pavg = 0.0

      if( .not. allocated( ulspg ) ) allocate( ulspg(kb:ke) )
      ulspg = 0.0
      if( .not. allocated( vlspg ) ) allocate( vlspg(kb:ke) )
      vlspg = 0.0

      if( .not. allocated( qavg ) ) allocate( qavg(kb:ke,numq) )
      qavg = 0.0

      if( .not. allocated( cavg ) ) allocate( cavg(kb:ke,3+numq) )
      cavg = 0.0

      if( .not. allocated( kmw ) ) allocate( kmw(kb:ke) )
      kmw = 0.0
      if( .not. allocated( ufw ) ) allocate( ufw(kb:ke) )
      ufw = 0.0
      if( .not. allocated( vfw ) ) allocate( vfw(kb:ke) )
      vfw = 0.0
      if( .not. allocated( u1b ) ) allocate( u1b(kb:ke) )
      u1b = 0.0
      if( .not. allocated( v1b ) ) allocate( v1b(kb:ke) )
      v1b = 0.0
      if( .not. allocated( l2p ) ) allocate( l2p(kb:ke) )
      l2p = 0.0
      if( .not. allocated( s2p ) ) allocate( s2p(kb:ke) )
      s2p = 0.0
      if( .not. allocated( s2b ) ) allocate( s2b(kb:ke) )
      s2b = 0.0
      if( .not. allocated( t2pm1 ) ) allocate( t2pm1(kb:ke) )
      t2pm1 = 0.0
      if( .not. allocated( t2pm2 ) ) allocate( t2pm2(kb:ke) )
      t2pm2 = 0.0
      if( .not. allocated( t2pm3 ) ) allocate( t2pm3(kb:ke) )
      t2pm3 = 0.0
      if( .not. allocated( t2pm4 ) ) allocate( t2pm4(kb:ke) )
      t2pm4 = 0.0

!----------------------------------------------------------------------
!  allocate the base state arrays, then call BASE

#ifdef MPI
          call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif
          if(myid.eq.0) print *,'  allocating base state arrays ... '
#ifdef MPI
          call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif

      if( .not. allocated( rho0s ) ) allocate( rho0s(ib:ie,jb:je) )
      rho0s = 0.0
      if( .not. allocated( pi0s ) ) allocate(  pi0s(ib:ie,jb:je) )
      pi0s = 0.0
      if( .not. allocated( prs0s ) ) allocate( prs0s(ib:ie,jb:je) )
      prs0s = 0.0
      if( .not. allocated( rth0s ) ) allocate( rth0s(ib:ie,jb:je) )
      rth0s = 0.0
      if( .not. allocated( pi0 ) ) allocate(  pi0(ib:ie,jb:je,kb:ke) )
      pi0 = 0.0
      if( .not. allocated( rho0 ) ) allocate( rho0(ib:ie,jb:je,kb:ke) )
      rho0 = 0.0
      if( .not. allocated( prs0 ) ) allocate( prs0(ib:ie,jb:je,kb:ke) )
      prs0 = 0.0
      if( .not. allocated( thv0 ) ) allocate( thv0(ib:ie,jb:je,kb:ke) )
      thv0 = 0.0
      if( .not. allocated( th0 ) ) allocate(  th0(ib:ie,jb:je,kb:ke) )
      th0 = 0.0
      if( .not. allocated( rth0 ) ) allocate( rth0(ib:ie,jb:je,kb:ke) )
      rth0 = 0.0
      if( .not. allocated( qv0 ) ) allocate(  qv0(ib:ie,jb:je,kb:ke) )
      qv0 = 0.0
      if( .not. allocated( qc0 ) ) allocate(  qc0(ib:ie,jb:je,kb:ke) )
      qc0 = 0.0
      if( .not. allocated( qi0 ) ) allocate(  qi0(ib:ie,jb:je,kb:ke) )
      qi0 = 0.0
      if( .not. allocated( rr0 ) ) allocate(  rr0(ib:ie,jb:je,kb:ke) )
      rr0 = 0.0
      if( .not. allocated( rf0 ) ) allocate(  rf0(ib:ie,jb:je,kb:ke) )
      rf0 = 0.0
      if( .not. allocated( rrf0 ) ) allocate( rrf0(ib:ie,jb:je,kb:ke) )
      rrf0 = 0.0
      if( .not. allocated( u0 ) ) allocate(   u0(ib:ie+1,jb:je,kb:ke) )
      u0 = 0.0
      if( .not. allocated( v0 ) ) allocate(   v0(ib:ie,jb:je+1,kb:ke) )
      v0 = 0.0

#ifdef MPI
          call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif
          if(myid.eq.0) print *,'  allocating dum, out2d, out3d arrays ... '
#ifdef MPI
          call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif

      if( .not. allocated( thrd ) ) allocate( thrd(ibb2:ibe2,jbb2:jbe2,kbb2:kbe2) )
      thrd = 0.0

      if( .not. allocated( dum1 ) ) allocate( dum1(ib:ie,jb:je,kb:ke) )
      dum1 = 0.0
      if( .not. allocated( dum2 ) ) allocate( dum2(ib:ie,jb:je,kb:ke) )
      dum2 = 0.0
      if( .not. allocated( dum3 ) ) allocate( dum3(ib:ie,jb:je,kb:ke) )
      dum3 = 0.0
      if( .not. allocated( dum4 ) ) allocate( dum4(ib:ie,jb:je,kb:ke) )
      dum4 = 0.0
      if( .not. allocated( dum5 ) ) allocate( dum5(ib:ie,jb:je,kb:ke) )
      dum5 = 0.0
      if( .not. allocated( dum6 ) ) allocate( dum6(ib:ie,jb:je,kb:ke) )
      dum6 = 0.0
      if( .not. allocated( dum7 ) ) allocate( dum7(ib:ie,jb:je,kb:ke) )
      dum7 = 0.0
      if( .not. allocated( dum8 ) ) allocate( dum8(ib:ie,jb:je,kb:ke) )
      dum8 = 0.0
      if( .not. allocated( dum9 ) ) allocate( dum9(ib:ie,jb:je,kb:ke) )
      dum9 = 0.0
      if( .not. allocated( dum10 ) ) allocate( dum10(ib:ie+1,jb:je+1,kb:ke) )
      dum10 = 0.0

      if( .not. allocated( out2d ) ) allocate( out2d(ib2d:ie2d,jb2d:je2d,nout2d) )
      out2d = 0.0

      if( .not. allocated( out3d ) ) allocate( out3d(ib3d:ie3d,jb3d:je3d,kb3d:ke3d,nout3d) )
      out3d = 0.0

      if(myid.eq.0) print *
#ifdef MPI
      call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif
      time_base=mytime()

      call       base(zh,mh,rmh,c1,c2,zf,mf,rho0s,pi0s,prs0s,rth0s,         &
                      wprof,ufrc,vfrc,thfrc,qvfrc,ug,vg,dvdr,               &
                      uavg,vavg,thavg,pavg,qavg,ulspg,vlspg,                &
                      pi0,prs0,rho0,thv0,th0,rth0,qv0,u0,v0,thrd,           &
                      qc0,qi0,rr0,rf0,rrf0,dum1,dum2,                       &
                      reqs_u,reqs_v,reqs_s,nw1,nw2,ne1,ne2,sw1,sw2,se1,se2, &
                      n3w1,n3w2,n3e1,n3e2,s3w1,s3w2,s3e1,s3e2,              &
                      uw31,uw32,ue31,ue32,us31,us32,un31,un32,              &
                      vw31,vw32,ve31,ve32,vs31,vs32,vn31,vn32,              &
                      sw31,sw32,se31,se32,ss31,ss32,sn31,sn32)

#ifdef MPI
          call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif
        time_base = mytime()
        if(dowr) write(outfile,*)
        if(dowr) write(outfile,*) '  time spent in subroutine base = ',time_base,' seconds '
        if(dowr) write(outfile,*)

!----------------------------------------------------------------------
!  Now, allocate the mother lode, then call INIT3D

#ifdef MPI
          call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif
          if(myid.eq.0) print *
          if(myid.eq.0) print *,'  allocating 2d vars ... '
#ifdef MPI
          call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif
      if( .not. allocated( rain ) ) allocate(   rain(ib:ie,jb:je,nrain) )
      rain = 0.0
      if( .not. allocated( raincv ) ) allocate( raincv(ib:ie,jb:je,nrain) )
      raincv = 0.0
      if( .not. allocated( hail ) ) allocate( hail(ib:ie,jb:je,nrain) )
      hail = 0.0
      if( .not. allocated( sws ) ) allocate(    sws(ib:ie,jb:je,nrain) )
      sws = 0.0
      if( .not. allocated( svs ) ) allocate(    svs(ib:ie,jb:je,nrain) )
      svs = 0.0
      if( .not. allocated( sps ) ) allocate(    sps(ib:ie,jb:je,nrain) )
      sps = 0.0
      if( .not. allocated( srs ) ) allocate(    srs(ib:ie,jb:je,nrain) )
      srs = 0.0
      if( .not. allocated( sgs ) ) allocate(    sgs(ib:ie,jb:je,nrain) )
      sgs = 0.0
      if( .not. allocated( sus ) ) allocate(    sus(ib:ie,jb:je,nrain) )
      sus = 0.0
      if( .not. allocated( shs ) ) allocate(    shs(ib:ie,jb:je,nrain) )
      shs = 0.0

      if( .not. allocated( tsk ) ) allocate(    tsk(ib:ie,jb:je) )
      tsk = 0.0
      if( .not. allocated( znt ) ) allocate(    znt(ib:ie,jb:je) )
      znt = 0.0
      if( .not. allocated( rznt ) ) allocate(   rznt(ib:ie,jb:je) )
      rznt = 0.0
      if( .not. allocated( zntmp ) ) allocate(  zntmp(ib:ie,jb:je) )
      zntmp = 0.0
      if( .not. allocated( ust ) ) allocate(    ust(ib:ie,jb:je) )
      ust = 0.0
      if( .not. allocated( stau ) ) allocate(   stau(ib:ie,jb:je) )
      stau = 0.0
      if( .not. allocated( tst ) ) allocate(    tst(ib:ie,jb:je) )
      tst = 0.0
      if( .not. allocated( qst ) ) allocate(    qst(ib:ie,jb:je) )
      qst = 0.0
      if( .not. allocated( z0t ) ) allocate(    z0t(ib:ie,jb:je) )
      z0t = 0.0
      if( .not. allocated( z0q ) ) allocate(    z0q(ib:ie,jb:je) )
      z0q = 0.0
      if( .not. allocated( thflux ) ) allocate( thflux(ib:ie,jb:je) )
      thflux = 0.0
      if( .not. allocated( qvflux ) ) allocate( qvflux(ib:ie,jb:je) )
      qvflux = 0.0
      if( .not. allocated( cd ) ) allocate(     cd(ib:ie,jb:je) )
      cd = 0.0
      if( .not. allocated( ch ) ) allocate(     ch(ib:ie,jb:je) )
      ch = 0.0
      if( .not. allocated( cq ) ) allocate(     cq(ib:ie,jb:je) )
      cq = 0.0
      if( .not. allocated( z1 ) ) allocate(     z1(ib:ie,jb:je) )
      z1 = 0.0
      if( .not. allocated( u1 ) ) allocate(     u1(ib:ie,jb:je) )
      u1 = 0.0
      if( .not. allocated( v1 ) ) allocate(     v1(ib:ie,jb:je) )
      v1 = 0.0
      if( .not. allocated( s1 ) ) allocate(     s1(ib:ie,jb:je) )
      s1 = 0.0
      if( .not. allocated( th1 ) ) allocate(     th1(ib:ie,jb:je) )
      th1 = 0.0
      if( .not. allocated( qv1 ) ) allocate(     qv1(ib:ie,jb:je) )
      qv1 = 0.0
      if( .not. allocated( xland ) ) allocate(  xland(ib:ie,jb:je) )
      xland = 0.0
      if( .not. allocated( psfc ) ) allocate(   psfc(ib:ie,jb:je) )
      psfc = 0.0
      if( .not. allocated( tlh ) ) allocate(    tlh(ib:ie,jb:je) )
      tlh = l_h
      if( .not. allocated( ustt ) ) allocate(   ustt(ib:ie,jb:je) )
      ustt = 0.0
      if( .not. allocated( ut ) ) allocate(     ut(ib:ie,jb:je) )
      ut = 0.0
      if( .not. allocated( vt ) ) allocate(     vt(ib:ie,jb:je) )
      vt = 0.0
      if( .not. allocated( st ) ) allocate(     st(ib:ie,jb:je) )
      st = 0.0
      if( .not. allocated( cm0 ) ) allocate(    cm0(ib:ie,jb:je) )
      cm0 = 0.0
      if( .not. allocated( eimaskw ) ) allocate( eimaskw(ib:ie,jb:je) )
      eimaskw = 0.0
      if( .not. allocated( eimaske ) ) allocate( eimaske(ib:ie,jb:je) )
      eimaske = 0.0
      if( .not. allocated( eimasks ) ) allocate( eimasks(ib:ie,jb:je) )
      eimasks = 0.0
      if( .not. allocated( eimaskn ) ) allocate( eimaskn(ib:ie,jb:je) )
      eimaskn = 0.0

      if( .not. allocated( radbcw ) ) allocate( radbcw(jb:je,kb:ke) )
      radbcw = 0.0
      if( .not. allocated( radbce ) ) allocate( radbce(jb:je,kb:ke) )
      radbce = 0.0
      if( .not. allocated( radbcs ) ) allocate( radbcs(ib:ie,kb:ke) )
      radbcs = 0.0
      if( .not. allocated( radbcn ) ) allocate( radbcn(ib:ie,kb:ke) )
      radbcn = 0.0

#ifdef MPI
          call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif
          if(myid.eq.0) print *,'  allocating tau + etc ... '
#ifdef MPI
          call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif
      if( .not. allocated( divx ) ) allocate( divx(ib:ie,jb:je,kb:ke) )
      divx = 0.0
      if( .not. allocated( rho ) ) allocate(  rho(ib:ie,jb:je,kb:ke) )
      rho = 0.0
      if( .not. allocated( rr ) ) allocate(   rr(ib:ie,jb:je,kb:ke) )
      rr = 0.0
      if( .not. allocated( rf ) ) allocate(   rf(ib:ie,jb:je,kb:ke) )
      rf = 0.0
      if( .not. allocated( prs ) ) allocate(  prs(ib:ie,jb:je,kb:ke) )
      prs = 0.0

#ifdef MPI
          call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif
          if(myid.eq.0) print *,'  allocating u,v,w ... '
#ifdef MPI
          call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif
      if( .not. allocated( rru ) ) allocate(   rru(ib:ie+1,jb:je,kb:ke) )
      rru = 0.0
      if( .not. allocated( ua ) ) allocate(    ua(ib:ie+1,jb:je,kb:ke) )
      ua = 0.0
      if( .not. allocated( u3d ) ) allocate(   u3d(ib:ie+1,jb:je,kb:ke) )
      u3d = 0.0
      if( .not. allocated( uten ) ) allocate(  uten(ib:ie+1,jb:je,kb:ke) )
      uten = 0.0
      if( .not. allocated( uten1 ) ) allocate( uten1(ib:ie+1,jb:je,kb:ke) )
      uten1 = 0.0

      if( .not. allocated( rrv ) ) allocate(   rrv(ib:ie,jb:je+1,kb:ke) )
      rrv = 0.0
      if( .not. allocated( va ) ) allocate(    va(ib:ie,jb:je+1,kb:ke) )
      va = 0.0
      if( .not. allocated( v3d ) ) allocate(   v3d(ib:ie,jb:je+1,kb:ke) )
      v3d = 0.0
      if( .not. allocated( vten ) ) allocate(  vten(ib:ie,jb:je+1,kb:ke) )
      vten = 0.0
      if( .not. allocated( vten1 ) ) allocate( vten1(ib:ie,jb:je+1,kb:ke) )
      vten1 = 0.0

      if( .not. allocated( rrw ) ) allocate(   rrw(ib:ie,jb:je,kb:ke+1) )
      rrw = 0.0
      if( .not. allocated( wa ) ) allocate(    wa(ib:ie,jb:je,kb:ke+1) )
      wa = 0.0
      if( .not. allocated( w3d ) ) allocate(   w3d(ib:ie,jb:je,kb:ke+1) )
      w3d = 0.0
      if( .not. allocated( wten ) ) allocate(  wten(ib:ie,jb:je,kb:ke+1) )
      wten = 0.0
      if( .not. allocated( wten1 ) ) allocate( wten1(ib:ie,jb:je,kb:ke+1) )
      wten1 = 0.0

#ifdef MPI
          call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif
          if(myid.eq.0) print *,'  allocating p,th ... '
#ifdef MPI
          call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif
      if( .not. allocated( ppi ) ) allocate(   ppi(ib:ie,jb:je,kb:ke) )
      ppi = 0.0
      if( .not. allocated( pp3d ) ) allocate(  pp3d(ib:ie,jb:je,kb:ke) )
      pp3d = 0.0
      if( .not. allocated( ppten ) ) allocate( ppten(ib:ie,jb:je,kb:ke) )
      ppten = 0.0
      if( .not. allocated( sten ) ) allocate(  sten(ib:ie,jb:je,kb:ke) )
      sten = 0.0
      if( .not. allocated( sadv ) ) allocate(  sadv(ib:ie,jb:je,kb:ke) )
      sadv = 0.0
      if( .not. allocated( ppx ) ) allocate(   ppx(ib:ie,jb:je,kb:ke) )
      ppx = 0.0

      if( .not. allocated( phi1 ) ) allocate(  phi1(ibph:ieph,jbph:jeph,kbph:keph) )
      phi1 = 0.0
      if( .not. allocated( phi2 ) ) allocate(  phi2(ibph:ieph,jbph:jeph,kbph:keph) )
      phi2 = 0.0

      if( .not. allocated( tha ) ) allocate(   tha(ib:ie,jb:je,kb:ke) )
      tha = 0.0
      if( .not. allocated( th3d ) ) allocate(  th3d(ib:ie,jb:je,kb:ke) )
      th3d = 0.0
      if( .not. allocated( thten ) ) allocate( thten(ib:ie,jb:je,kb:ke) )
      thten = 0.0
      if( .not. allocated( thten1 ) ) allocate(thten1(ib:ie,jb:je,kb:ke) )
      thten1 = 0.0
      if( .not. allocated( thterm ) ) allocate(thterm(ib:ie,jb:je,kb:ke) )
      thterm = 0.0

      if( .not. allocated( qpten ) ) allocate( qpten(ibm:iem,jbm:jem,kbm:kem) )
      qpten = 0.0
      if( .not. allocated( qtten ) ) allocate( qtten(ibm:iem,jbm:jem,kbm:kem) )
      qtten = 0.0
      if( .not. allocated( qvten ) ) allocate( qvten(ibm:iem,jbm:jem,kbm:kem) )
      qvten = 0.0
      if( .not. allocated( qcten ) ) allocate( qcten(ibm:iem,jbm:jem,kbm:kem) )
      qcten = 0.0

      if( .not. allocated( bud ) ) allocate(   bud(nk) )
      bud = 0.0
      if( .not. allocated( bud2 ) ) allocate(  bud2(nj) )
      bud2 = 0.0
      if( .not. allocated( qbudget ) ) allocate( qbudget(nbudget) )
      qbudget = 0.0
      if( .not. allocated( asq ) ) allocate(    asq(numq) )
      asq = 0.0
      if( .not. allocated( bsq ) ) allocate(    bsq(numq) )
      bsq = 0.0

#ifdef MPI
          call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif
          if(myid.eq.0) print *,'  allocating qa,q3d,qten ... '
#ifdef MPI
          call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif
      if( .not. allocated( qa ) ) allocate(     qa(ibm:iem,jbm:jem,kbm:kem,numq) )
      qa = 0.0
      if( .not. allocated( q3d ) ) allocate(    q3d(ibm:iem,jbm:jem,kbm:kem,numq) )
      q3d = 0.0
      if( .not. allocated( qten ) ) allocate(   qten(ibm:iem,jbm:jem,kbm:kem,numq) )
      qten = 0.0

      if( .not. allocated( p3a ) ) allocate( p3a(ibp3:iep3,kbp3:kep3,np3a) )
      p3a = 0.0
      if( .not. allocated( p3o ) ) allocate( p3o(ibp3:iep3,jbp3:jep3,kbp3:kep3,np3o) )
      p3o = 0.0

#ifdef MPI
          call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif
          if(myid.eq.0) print *,'  allocating k,tke ... '
#ifdef MPI
          call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif

      if( cm1setup.eq.3 )then
        kminit = viscosity
        khinit = viscosity/pr_num
      else
        kminit = 0.0
        khinit = 0.0
      endif

      if( .not. allocated( kmh ) ) allocate(    kmh(ibc:iec,jbc:jec,kbc:kec) )
      kmh = kminit
      if( .not. allocated( kmv ) ) allocate(    kmv(ibc:iec,jbc:jec,kbc:kec) )
      kmv = kminit
      if( .not. allocated( khh ) ) allocate(    khh(ibc:iec,jbc:jec,kbc:kec) )
      khh = khinit
      if( .not. allocated( khv ) ) allocate(    khv(ibc:iec,jbc:jec,kbc:kec) )
      khv = khinit

      if( .not. allocated( kmv2 ) ) allocate(    kmv2(ibc:iec,jbc:jec,kbc:kec) )
      kmv2 = kminit
      if( .not. allocated( khv2 ) ) allocate(    khv2(ibc:iec,jbc:jec,kbc:kec) )
      khv2 = khinit

      if( .not. allocated( cme ) ) allocate(    cme(ibc:iec,jbc:jec,kbc:kec) )
      cme = 0.0
      if( .not. allocated( csm ) ) allocate(    csm(ibc:iec,jbc:jec,kbc:kec) )
      csm = 0.0
      if( .not. allocated( ce1 ) ) allocate(    ce1(ibc:iec,jbc:jec,kbc:kec) )
      ce1 = 0.0
      if( .not. allocated( ce2 ) ) allocate(    ce2(ibc:iec,jbc:jec,kbc:kec) )
      ce2 = 0.0

      if( .not. allocated( tkea ) ) allocate(   tkea(ibt:iet,jbt:jet,kbt:ket) )
      tkea = 0.0
      if( .not. allocated( tke3d ) ) allocate(  tke3d(ibt:iet,jbt:jet,kbt:ket) )
      tke3d = 0.0
      if( .not. allocated( tketen ) ) allocate( tketen(ibt:iet,jbt:jet,kbt:ket) )
      tketen = 0.0

#ifdef MPI
          call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif
          if(myid.eq.0) print *,'  allocating pbl arrays ... '
#ifdef MPI
          call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif
      if( .not. allocated( thpten ) ) allocate( thpten(ibb:ieb,jbb:jeb,kbb:keb) )
      thpten = 0.0
      if( .not. allocated( qvpten ) ) allocate( qvpten(ibb:ieb,jbb:jeb,kbb:keb) )
      qvpten = 0.0
      if( .not. allocated( qcpten ) ) allocate( qcpten(ibb:ieb,jbb:jeb,kbb:keb) )
      qcpten = 0.0
      if( .not. allocated( qipten ) ) allocate( qipten(ibb:ieb,jbb:jeb,kbb:keb) )
      qipten = 0.0
      if( .not. allocated( upten ) ) allocate(  upten(ibb:ieb,jbb:jeb,kbb:keb) )
      upten = 0.0
      if( .not. allocated( vpten ) ) allocate(  vpten(ibb:ieb,jbb:jeb,kbb:keb) )
      vpten = 0.0
      if( .not. allocated( qnipten ) ) allocate( qnipten(ibb:ieb,jbb:jeb,kbb:keb) )
      qnipten = 0.0
      if( .not. allocated( qncpten ) ) allocate( qncpten(ibb:ieb,jbb:jeb,kbb:keb) )
      qncpten = 0.0
      if( .not. allocated( qrpten ) ) allocate( qrpten(ibb:ieb,jbb:jeb,kbb:keb) )
      qrpten = 0.0
      if( .not. allocated( qspten ) ) allocate( qspten(ibb:ieb,jbb:jeb,kbb:keb) )
      qspten = 0.0
      if( .not. allocated( qgpten ) ) allocate( qgpten(ibb:ieb,jbb:jeb,kbb:keb) )
      qgpten = 0.0

      !-----------------
      ! begin radiation

#ifdef MPI
          call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif
          if(myid.eq.0) print *,'  allocating rad vars ... '
#ifdef MPI
          call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif
      if( .not. allocated( swten ) ) allocate( swten(ibr:ier,jbr:jer,kbr:ker) )
      swten = 0.0
      if( .not. allocated( lwten ) ) allocate( lwten(ibr:ier,jbr:jer,kbr:ker) )
      lwten = 0.0
      if( .not. allocated( swtenc ) ) allocate( swtenc(ibr:ier,jbr:jer,kbr:ker) )
      swtenc = 0.0
      if( .not. allocated( lwtenc ) ) allocate( lwtenc(ibr:ier,jbr:jer,kbr:ker) )
      lwtenc = 0.0
      if( .not. allocated( thraten ) ) allocate( thraten(ibr:ier,jbr:jer,kbr:ker) )
      thraten = 0.0
      if( .not. allocated( cldfra ) ) allocate(cldfra(ibr:ier,jbr:jer,kbr:ker) )
      cldfra = 0.0
      if( .not. allocated( o30 ) ) allocate(   o30(ibr:ier,jbr:jer,kbr:ker) )
      o30 = 0.0
      if( .not. allocated( zir ) ) allocate(   zir(ibr:ier,jbr:jer) )
      zir = 0.0

      IF( radopt .eq. 1 )THEN
        nir = 1
        njr = 1
        nkr = nk+3
        rbufsz = n2d_radiat*nir*njr + n3d_radiat*nir*njr*nkr
      ELSE
        nir = 1
        njr = 1
        nkr = 1
        rbufsz = 1
      ENDIF

      if( .not. allocated( radsw ) ) allocate(    radsw(ni,nj) )
      radsw = 0.0
      if( .not. allocated( rnflx ) ) allocate(    rnflx(ni,nj) )
      rnflx = 0.0
      if( .not. allocated( radswnet ) ) allocate( radswnet(ni,nj) )
      radswnet = 0.0
      if( .not. allocated( radlwin ) ) allocate(  radlwin(ni,nj) )
      radlwin = 0.0
      if( .not. allocated( dsr ) ) allocate(      dsr(ni,nj) )
      dsr = 0.0
      if( .not. allocated( olr ) ) allocate(      olr(ni,nj) )
      olr = 0.0

      if( .not. allocated( rad2d ) ) allocate(    rad2d(ni,nj,nrad2d) )
      rad2d = 0.0

      if( .not. allocated( effc ) ) allocate(  effc(ibr:ier,jbr:jer,kbr:ker) )
      effc = 25.0e-6
      if( .not. allocated( effi ) ) allocate(  effi(ibr:ier,jbr:jer,kbr:ker) )
      effi = 25.0e-6
      if( .not. allocated( effs ) ) allocate(  effs(ibr:ier,jbr:jer,kbr:ker) )
      effs = 25.0e-6
      if( .not. allocated( effr ) ) allocate(  effr(ibr:ier,jbr:jer,kbr:ker) )
      effr = 25.0e-6
      if( .not. allocated( effg ) ) allocate(  effg(ibr:ier,jbr:jer,kbr:ker) )
      effg = 25.0e-6
      if( .not. allocated( effis ) ) allocate( effis(ibr:ier,jbr:jer,kbr:ker) )
      effis = 25.0e-6

      if( .not. allocated( lwupt ) ) allocate( lwupt(ibr:ier,jbr:jer) )
      lwupt = 0.0
      if( .not. allocated( lwuptc ) ) allocate( lwuptc(ibr:ier,jbr:jer) )
      lwuptc = 0.0
      if( .not. allocated( lwdnt ) ) allocate( lwdnt(ibr:ier,jbr:jer) )
      lwdnt = 0.0
      if( .not. allocated( lwdntc ) ) allocate( lwdntc(ibr:ier,jbr:jer) )
      lwdntc = 0.0
      if( .not. allocated( lwupb ) ) allocate( lwupb(ibr:ier,jbr:jer) )
      lwupb = 0.0
      if( .not. allocated( lwupbc ) ) allocate( lwupbc(ibr:ier,jbr:jer) )
      lwupbc = 0.0
      if( .not. allocated( lwdnb ) ) allocate( lwdnb(ibr:ier,jbr:jer) )
      lwdnb = 0.0
      if( .not. allocated( lwdnbc ) ) allocate( lwdnbc(ibr:ier,jbr:jer) )
      lwdnbc = 0.0

      if( .not. allocated( swupt ) ) allocate( swupt(ibr:ier,jbr:jer) )
      swupt = 0.0
      if( .not. allocated( swuptc ) ) allocate( swuptc(ibr:ier,jbr:jer) )
      swuptc = 0.0
      if( .not. allocated( swdnt ) ) allocate( swdnt(ibr:ier,jbr:jer) )
      swdnt = 0.0
      if( .not. allocated( swdntc ) ) allocate( swdntc(ibr:ier,jbr:jer) )
      swdntc = 0.0
      if( .not. allocated( swupb ) ) allocate( swupb(ibr:ier,jbr:jer) )
      swupb = 0.0
      if( .not. allocated( swupbc ) ) allocate( swupbc(ibr:ier,jbr:jer) )
      swupbc = 0.0
      if( .not. allocated( swdnb ) ) allocate( swdnb(ibr:ier,jbr:jer) )
      swdnb = 0.0
      if( .not. allocated( swdnbc ) ) allocate( swdnbc(ibr:ier,jbr:jer) )
      swdnbc = 0.0

      if( .not. allocated( lwcf ) ) allocate(   lwcf(ibr:ier,jbr:jer) )
      lwcf = 0.0
      if( .not. allocated( swcf ) ) allocate(   swcf(ibr:ier,jbr:jer) )
      swcf = 0.0
      if( .not. allocated( coszr ) ) allocate(  coszr(ibr:ier,jbr:jer) )
      coszr = 0.0

      if( .not. allocated( xice ) ) allocate( xice(ibr:ier,jbr:jer) )
      xice = 0.0
      if( .not. allocated( xsnow ) ) allocate( xsnow(ibr:ier,jbr:jer) )
      xsnow = 0.0
      if( .not. allocated( xlat ) ) allocate( xlat(ibr:ier,jbr:jer) )
      xlat = 0.0
      if( .not. allocated( xlong ) ) allocate( xlong(ibr:ier,jbr:jer) )
      xlong = 0.0
      if( .not. allocated( coszen ) ) allocate( coszen(ibr:ier,jbr:jer) )
      coszen = 0.0
      if( .not. allocated( swddir ) ) allocate( swddir(ibr:ier,jbr:jer) )
      swddir = 0.0
      if( .not. allocated( swddni ) ) allocate( swddni(ibr:ier,jbr:jer) )
      swddni = 0.0
      if( .not. allocated( swddif ) ) allocate( swddif(ibr:ier,jbr:jer) )
      swddif = 0.0
      if( .not. allocated( hrang ) ) allocate( hrang(ibr:ier,jbr:jer) )
      hrang = 0.0

      if( .not. allocated( cldfra1_flag ) ) allocate( cldfra1_flag(ibr:ier,jbr:jer,kbr:ker) )
      cldfra1_flag = 0

      if(dowr) write(outfile,*) '  rbufsz,nrad2d = ',rbufsz,nrad2d


      ! end radiation
      !-----------------

#ifdef MPI
          call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif
          if(myid.eq.0) print *,'  allocating sfc vars ... '
#ifdef MPI
          call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif
      if( .not. allocated( lu_index ) ) allocate( lu_index(ibl:iel,jbl:jel) )
      lu_index = 0
      if( .not. allocated( kpbl2d ) ) allocate(   kpbl2d(ibl:iel,jbl:jel) )
      kpbl2d = 0
      if( .not. allocated( u10 ) ) allocate(      u10(ibl:iel,jbl:jel) )
      u10 = 0.0
      if( .not. allocated( v10 ) ) allocate(      v10(ibl:iel,jbl:jel) )
      v10 = 0.0
      if( .not. allocated( s10 ) ) allocate(      s10(ibl:iel,jbl:jel) )
      s10 = 0.0
      if( .not. allocated( hfx ) ) allocate(      hfx(ibl:iel,jbl:jel) )
      hfx = 0.0
      if( .not. allocated( qfx ) ) allocate(      qfx(ibl:iel,jbl:jel) )
      qfx = 0.0
      if( .not. allocated( hpbl ) ) allocate(     hpbl(ibl:iel,jbl:jel) )
      hpbl = 100.0
      if( .not. allocated( wspd ) ) allocate(     wspd(ibl:iel,jbl:jel) )
      wspd = 0.0
      if( .not. allocated( phim ) ) allocate(     phim(ibl:iel,jbl:jel) )
      phim = 0.0
      if( .not. allocated( phih ) ) allocate(     phih(ibl:iel,jbl:jel) )
      phih = 0.0
      if( .not. allocated( psim ) ) allocate(     psim(ibl:iel,jbl:jel) )
      psim = 0.0
      if( .not. allocated( psih ) ) allocate(     psih(ibl:iel,jbl:jel) )
      psih = 0.0
      if( .not. allocated( psiq ) ) allocate(     psiq(ibl:iel,jbl:jel) )
      psiq = 0.0
      if( .not. allocated( gz1oz0 ) ) allocate(   gz1oz0(ibl:iel,jbl:jel) )
      gz1oz0 = 0.0
      if( .not. allocated( br ) ) allocate(       br(ibl:iel,jbl:jel) )
      br = 0.0
      if( .not. allocated( brcr ) ) allocate(     brcr(ibl:iel,jbl:jel) )
      brcr = 0.0
      if( .not. allocated( chs ) ) allocate(      chs(ibl:iel,jbl:jel) )
      chs = 0.0
      if( .not. allocated( chs2 ) ) allocate(     chs2(ibl:iel,jbl:jel) )
      chs2 = 0.0
      if( .not. allocated( cqs2 ) ) allocate(     cqs2(ibl:iel,jbl:jel) )
      cqs2 = 0.0
      if( .not. allocated( cpmm ) ) allocate(     cpmm(ibl:iel,jbl:jel) )
      cpmm = 0.0
      if( .not. allocated( zol ) ) allocate(      zol(ibl:iel,jbl:jel) )
      zol = 0.0
      if( .not. allocated( mavail ) ) allocate(   mavail(ibl:iel,jbl:jel) )
      mavail = 0.0
      if( .not. allocated( mol ) ) allocate(      mol(ibl:iel,jbl:jel) )
      mol = 0.0
      if( .not. allocated( rmol ) ) allocate(     rmol(ibl:iel,jbl:jel) )
      rmol = 0.0
      if( .not. allocated( regime ) ) allocate(   regime(ibl:iel,jbl:jel) )
      regime = 0.0
      if( .not. allocated( lh ) ) allocate(       lh(ibl:iel,jbl:jel) )
      lh = 0.0
      if( .not. allocated( flhc ) ) allocate(     flhc(ibl:iel,jbl:jel) )
      flhc = 0.0
      if( .not. allocated( flqc ) ) allocate(     flqc(ibl:iel,jbl:jel) )
      flqc = 0.0
      if( .not. allocated( qgh ) ) allocate(      qgh(ibl:iel,jbl:jel) )
      qgh = 0.0
      if( .not. allocated( ck ) ) allocate(       ck(ibl:iel,jbl:jel) )
      ck = 0.0
      if( .not. allocated( cka ) ) allocate(      cka(ibl:iel,jbl:jel) )
      cka = 0.0
      if( .not. allocated( cda ) ) allocate(      cda(ibl:iel,jbl:jel) )
      cda = 0.0
      if( .not. allocated( ustm ) ) allocate(     ustm(ibl:iel,jbl:jel) )
      ustm = 0.0
      if( .not. allocated( qsfc ) ) allocate(     qsfc(ibl:iel,jbl:jel) )
      qsfc = 0.0
      if( .not. allocated( t2 ) ) allocate(       t2(ibl:iel,jbl:jel) )
      t2 = 0.0
      if( .not. allocated( q2 ) ) allocate(       q2(ibl:iel,jbl:jel) )
      q2 = 0.0
      if( .not. allocated( th2 ) ) allocate(      th2(ibl:iel,jbl:jel) )
      th2 = 0.0
      if( .not. allocated( emiss ) ) allocate(    emiss(ibl:iel,jbl:jel) )
      emiss = 0.0
      if( .not. allocated( thc ) ) allocate(      thc(ibl:iel,jbl:jel) )
      thc = 0.0
      if( .not. allocated( albd ) ) allocate(     albd(ibl:iel,jbl:jel) )
      albd = 0.0
      if( .not. allocated( gsw ) ) allocate(      gsw(ibl:iel,jbl:jel) )
      gsw = 0.0
      if( .not. allocated( glw ) ) allocate(      glw(ibl:iel,jbl:jel) )
      glw = 0.0
      if( .not. allocated( chklowq ) ) allocate(  chklowq(ibl:iel,jbl:jel) )
      chklowq = 0.0
      if( .not. allocated( capg ) ) allocate(     capg(ibl:iel,jbl:jel) )
      capg = 0.0
      if( .not. allocated( snowc ) ) allocate(    snowc(ibl:iel,jbl:jel) )
      snowc = 0.0
      if( .not. allocated( snowh ) ) allocate(    snowh(ibl:iel,jbl:jel) )
      snowh = 0.0
      if( .not. allocated( qcg ) ) allocate(      qcg(ibl:iel,jbl:jel) )
      qcg = 0.0
      if( .not. allocated( dsxy ) ) allocate(     dsxy(ibl:iel,jbl:jel) )
      dsxy = 0.0
      if( .not. allocated( wstar ) ) allocate(    wstar(ibl:iel,jbl:jel) )
      wstar = 0.0
      if( .not. allocated( delta ) ) allocate(    delta(ibl:iel,jbl:jel) )
      delta = 0.0
      if( .not. allocated( prkpp ) ) allocate(    prkpp(ibl:iel,jbl:jel) )
      prkpp = 0.0
      if( .not. allocated( fm ) ) allocate(       fm(ibl:iel,jbl:jel) )
      fm = 0.0
      if( .not. allocated( fh ) ) allocate(       fh(ibl:iel,jbl:jel) )
      fh = 0.0

      ! for hwrf sfc:
      if( .not. allocated( charn ) ) allocate(    charn(ibl:iel,jbl:jel) )
      charn = 0.0185
      if( .not. allocated( msang ) ) allocate(    msang(ibl:iel,jbl:jel) )
      msang = 0.0
      if( .not. allocated( scurx ) ) allocate(    scurx(ibl:iel,jbl:jel) )
      scurx = 0.0
      if( .not. allocated( scury ) ) allocate(    scury(ibl:iel,jbl:jel) )
      scury = 0.0
      if( .not. allocated( zkmax ) ) allocate(    zkmax(ibl:iel,jbl:jel) )
      zkmax = 0.0
      if( .not. allocated( cd_out ) ) allocate(   cd_out(ibl:iel,jbl:jel) )
      cd_out = 0.0
      if( .not. allocated( ch_out ) ) allocate(   ch_out(ibl:iel,jbl:jel) )
      ch_out = 0.0
      if( .not. allocated( wscale ) ) allocate(   wscale(ibl:iel,jbl:jel) )
      wscale = 0.0
      if( .not. allocated( wscaleu ) ) allocate(  wscaleu(ibl:iel,jbl:jel) )
      wscaleu = 0.0
      if( .not. allocated( mznt ) ) allocate(     mznt(ibl:iel,jbl:jel) )
      mznt = 0.0
      if( .not. allocated( swspd ) ) allocate(    swspd(ibl:iel,jbl:jel) )
      swspd = 0.0
      if( .not. allocated( smois ) ) allocate(    smois(ibl:iel,jbl:jel) )
      smois = 0.0
      if( .not. allocated( taux ) ) allocate(     taux(ibl:iel,jbl:jel) )
      taux = 0.0
      if( .not. allocated( tauy ) ) allocate(     tauy(ibl:iel,jbl:jel) )
      tauy = 0.0
      if( .not. allocated( hpbl2d ) ) allocate(   hpbl2d(ibl:iel,jbl:jel) )
      hpbl2d = 0.0
      if( .not. allocated( evap2d ) ) allocate(   evap2d(ibl:iel,jbl:jel) )
      evap2d = 0.0
      if( .not. allocated( heat2d ) ) allocate(   heat2d(ibl:iel,jbl:jel) )
      heat2d = 0.0

      if( .not. allocated( pbl3d ) ) allocate( pbl3d(ibpbl:iepbl,jbpbl:jepbl,kbpbl:kepbl,npbl3d) )
      pbl3d = 0.0
      if( .not. allocated( pbl2d ) ) allocate( pbl2d(ibpbl:iepbl,jbpbl:jepbl,npbl2d) )
      pbl2d = 0.0
      if( .not. allocated( pbl2dint ) ) allocate( pbl2dint(ibpbl:iepbl,jbpbl:jepbl,npbl2dint) )
      pbl2dint = 0

      if( .not. allocated( tmp2d ) ) allocate( tmp2d(ibtmp:ietmp,kbtmp:ketmp,ntmp2d) )
      tmp2d = 0.0
      if( .not. allocated( tmp1d ) ) allocate( tmp1d(ibtmp:ietmp,ntmp1d) )
      tmp1d = 0.0

      ! start with very small, but non-zero, numbers:
      tem = 1.0e-6
      znt = tem
      rznt = 1.0/tem
      zntmp = tem
      ust = 1.0e-6
      ! to prevent divide-by-zeros on init for some combinations of namelist params:
      tsk  = 300.0
      psfc = 100000.0
      qsfc = 0.00001

      num_soil_layers = 5
      if( .not. allocated( slab_zs ) ) allocate(  slab_zs(num_soil_layers) )
      slab_zs = 0.0
      if( .not. allocated( slab_dzs ) ) allocate( slab_dzs(num_soil_layers) )
      slab_dzs = 0.0
      if( .not. allocated( tslb ) ) allocate(  tslb(ibl:iel,jbl:jel,num_soil_layers) )
      tslb = 0.0
      if( .not. allocated( tmn ) ) allocate(   tmn(ibl:iel,jbl:jel) )
      tmn = 0.0

      ! arrays for oml model:
      if( .not. allocated( tml ) ) allocate(   tml(ibl:iel,jbl:jel) )
      tml = 0.0
      if( .not. allocated( t0ml ) ) allocate(  t0ml(ibl:iel,jbl:jel) )
      t0ml = 0.0
      if( .not. allocated( hml ) ) allocate(   hml(ibl:iel,jbl:jel) )
      hml = 0.0
      if( .not. allocated( h0ml ) ) allocate(  h0ml(ibl:iel,jbl:jel) )
      h0ml = 0.0
      if( .not. allocated( huml ) ) allocate(  huml(ibl:iel,jbl:jel) )
      huml = 0.0
      if( .not. allocated( hvml ) ) allocate(  hvml(ibl:iel,jbl:jel) )
      hvml = 0.0
      if( .not. allocated( tmoml ) ) allocate( tmoml(ibl:iel,jbl:jel) )
      tmoml = 0.0

      if( .not. allocated( pta ) ) allocate(    pta(ibp:iep,jbp:jep,kbp:kep,npt) )
      pta = 0.0
      if( .not. allocated( pt3d ) ) allocate(   pt3d(ibp:iep,jbp:jep,kbp:kep,npt) )
      pt3d = 0.0
      if( .not. allocated( ptten ) ) allocate(  ptten(ibp:iep,jbp:jep,kbp:kep,npt) )
      ptten = 0.0

!----------------------------------------------------------------------
!  arrays for cm1ib:

      if( .not. allocated( bndy ) ) allocate( bndy(ibib:ieib,jbib:jeib,kbib:keib) )
      bndy = .false.

      if( .not. allocated( kbdy ) ) allocate( kbdy(ibib:ieib,jbib:jeib) )
      kbdy = 1

      IF( do_ib )THEN

        call   init_immersed_boundaries(                                          &
                       xh,yh,xf,yf,sigma,sigmaf,zs,zh,bndy,kbdy,out3d,            &
                       ww31(1,1,1),ww32(1,1,1),we31(1,1,1),we32(1,1,1),           &
                       ws31(1,1,1),ws32(1,1,1),wn31(1,1,1),wn32(1,1,1),reqs_p)

      ENDIF

      if( .not. allocated( hflxw ) ) allocate( hflxw(ibib:ieib,jbib:jeib,kmaxib) )
      hflxw = 0
      if( .not. allocated( hflxe ) ) allocate( hflxe(ibib:ieib,jbib:jeib,kmaxib) )
      hflxe = 0
      if( .not. allocated( hflxs ) ) allocate( hflxs(ibib:ieib,jbib:jeib,kmaxib) )
      hflxs = 0
      if( .not. allocated( hflxn ) ) allocate( hflxn(ibib:ieib,jbib:jeib,kmaxib) )
      hflxn = 0

      IF( do_ib )THEN

        call   ib_flx_init(bndy,hflxw,hflxe,hflxs,hflxn,dum1,dum2,     &
                           uw31,uw32,ue31,ue32,us31,us32,un31,un32,    &
                           vw31,vw32,ve31,ve32,vs31,vs32,vn31,vn32,reqs_u,reqs_v)

      ENDIF

!----------------------------------------------------------------------
!  allocate arrays and initialize some variables for large-scale nudging:

      IF( do_lsnudge )THEN
        lsnudge_maxlevels = max( nk , 1000 )
      ELSE
        lsnudge_maxlevels = 1
      ENDIF

      lsnudge_unit = 74
      lsnudge_count = 0
      lsnudge_time1 = 0.0
      lsnudge_time2 = 0.0
      lsnudgefac = 0.0

      if( .not. allocated( lsnudge_u ) ) allocate( lsnudge_u(lsnudge_maxlevels,lsnudge_nmax) )
      lsnudge_u = 0.0
      if( .not. allocated( lsnudge_v ) ) allocate( lsnudge_v(lsnudge_maxlevels,lsnudge_nmax) )
      lsnudge_v = 0.0
      if( .not. allocated( lsnudge_th ) ) allocate( lsnudge_th(lsnudge_maxlevels,lsnudge_nmax) )
      lsnudge_th = 0.0
      if( .not. allocated( lsnudge_qv ) ) allocate( lsnudge_qv(lsnudge_maxlevels,lsnudge_nmax) )
      lsnudge_qv = 0.0

!----------------------------------------------------------------------

#ifdef MPI
          call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif
          if(myid.eq.0) print *,'  allocating dat vars ... '
#ifdef MPI
          call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif
      if( .not. allocated( dat1 ) ) allocate( dat1(d3i,d3j) )        ! my subdomain
      dat1 = 0.0
      if( .not. allocated( dat2 ) ) allocate( dat2(d2i,d2j) )        ! entire grid (if needed)
      dat2 = 0.0
      if( .not. allocated( dat3 ) ) allocate( dat3(d3i,d3j,d3n) )    ! entire node (if needed)
      dat3 = 0.0
      if( .not. allocated( reqt ) ) allocate( reqt(d3t) )
      reqt = 0


#ifdef MPI
          call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif
          if(myid.eq.0) print *,'  allocating parcel vars ... '
#ifdef MPI
          call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif
      if( .not. allocated( pdata ) ) allocate(  pdata(nparcels,npvals) )
      pdata = 0.0
      if( .not. allocated( ploc ) ) allocate(   ploc(nparcels,  3   ) )
      ploc = 0.0

      if( .not. allocated( haildata ) ) allocate( haildata(nhailtrajs,nhvals) )
      haildata = 0.0
      if( .not. allocated( hailloc ) ) allocate(  hailloc(nhailtrajs,  13   ) )
      hailloc = 0.0

      if( .not. allocated( flag ) ) allocate( flag(ib:ie,jb:je,kb:ke) )
      flag = .false.

#ifdef MPI
          call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif
          if(myid.eq.0) print *,'  allocating misc vars ... '
#ifdef MPI
          call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif
      if( .not. allocated( cfb ) ) allocate(    cfb(ipb:ipe,jpb:jpe,kpb:kpe) )
      cfb = 0.0
      if( .not. allocated( cfa ) ) allocate(    cfa(kpb:kpe) )
      cfa = 0.0
      if( .not. allocated( cfc ) ) allocate(    cfc(kpb:kpe) )
      cfc = 0.0
      if( .not. allocated( d1 ) ) allocate(     d1(kpb:kpe) )
      d1 = 0.0
      if( .not. allocated( d2 ) ) allocate(     d2(kpb:kpe) )
      d2 = 0.0
      if( .not. allocated( pdt ) ) allocate(    pdt(ipb:ipe,jpb:jpe,kpb:kpe) )
      pdt = 0.0
      if( .not. allocated( lgbth ) ) allocate(  lgbth(ipb:ipe,jpb:jpe,kpb:kpe) )
      lgbth = 0.0
      if( .not. allocated( lgbph ) ) allocate(  lgbph(ipb:ipe,jpb:jpe,kpb:kpe) )
      lgbph = 0.0

      if( .not. allocated( pgrid1 ) ) allocate( pgrid1(pi1a:pi2a,pj1a:pj2a,pk1a:pk2a) )
      pgrid1 = 0.0

      if( .not. allocated( cgrid1 ) ) allocate( cgrid1(pi1a:pi2a,pj1a:pj2a,pk1a:pk2a) )
      cgrid1 = 0.0
      if( .not. allocated( cgrid2 ) ) allocate( cgrid2(pi1b:pi2b,pj1b:pj2b,pk1b:pk2b) )
      cgrid2 = 0.0

      if( .not. allocated( psenda ) ) allocate( psenda(pimaxa,pjmaxa,pkmaxa,pnmaxa) )
      psenda = 0.0
      if( .not. allocated( precva ) ) allocate( precva(pimaxa,pjmaxa,pkmaxa,pnmaxa) )
      precva = 0.0

      if( .not. allocated( psendb ) ) allocate( psendb(pimaxb,pjmaxb,pkmaxb,pnmaxb) )
      psendb = 0.0
      if( .not. allocated( precvb ) ) allocate( precvb(pimaxb,pjmaxb,pkmaxb,pnmaxb) )
      precvb = 0.0

      if( .not. allocated( psendc ) ) allocate( psendc(pimaxc,pjmaxc,pkmaxc,pnmaxc) )
      psendc = 0.0
      if( .not. allocated( precvc ) ) allocate( precvc(pimaxc,pjmaxc,pkmaxc,pnmaxc) )
      precvc = 0.0

      if( .not. allocated( psendd ) ) allocate( psendd(pimaxd,pjmaxd,pkmaxd,pnmaxd) )
      psendd = 0.0
      if( .not. allocated( precvd ) ) allocate( precvd(pimaxd,pjmaxd,pkmaxd,pnmaxd) )
      precvd = 0.0

      if( .not. allocated( rhsx ) ) allocate(   rhsx(pi1b:pi2b) )
      rhsx = 0.0
      if( .not. allocated( transx ) ) allocate( transx(pi1b:pi2b) )
      transx = 0.0
      if( .not. allocated( rhsy ) ) allocate(   rhsy(pj1a:pj2a) )
      rhsy = 0.0
      if( .not. allocated( transy ) ) allocate( transy(pj1a:pj2a) )
      transy = 0.0

      if( .not. allocated( cir ) ) allocate(    cir(ib:ie,jb:je) )
      cir = 0
      if( .not. allocated( crr ) ) allocate(    crr(ib:ie,jb:je) )
      crr = 0.0
      if( .not. allocated( cangle ) ) allocate( cangle(ib:ie,jb:je) )
      cangle = 0.0

      ! cm1r21.1: for cu param
      if( .not. allocated( thcuten ) ) allocate( thcuten(ibcu:iecu,jbcu:jecu,kbcu:kecu) )
      thcuten = 0.0
      if( .not. allocated( qvcuten ) ) allocate( qvcuten(ibcu:iecu,jbcu:jecu,kbcu:kecu) )
      qvcuten = 0.0
      if( .not. allocated( qccuten ) ) allocate( qccuten(ibcu:iecu,jbcu:jecu,kbcu:kecu) )
      qccuten = 0.0
      if( .not. allocated( qicuten ) ) allocate( qicuten(ibcu:iecu,jbcu:jecu,kbcu:kecu) )
      qicuten = 0.0
      if( .not. allocated( ucuten ) ) allocate(  ucuten(ibcu:iecu,jbcu:jecu,kbcu:kecu) )
      ucuten = 0.0
      if( .not. allocated( vcuten ) ) allocate(  vcuten(ibcu:iecu,jbcu:jecu,kbcu:kecu) )
      vcuten = 0.0
      if( .not. allocated( thften ) ) allocate(  thften(ibcu:iecu,jbcu:jecu,kbcu:kecu) )
      thften = 0.0
      if( .not. allocated( qvften ) ) allocate(  qvften(ibcu:iecu,jbcu:jecu,kbcu:kecu) )
      qvften = 0.0

!----------------------------------------------------------------------

      if(myid.eq.0) print *
#ifdef MPI
      call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif
      time_init3d=mytime()

      call       init3d(dt,xh,rxh,uh,ruh,xf,rxf,uf,ruf,yh,vh,rvh,yf,vf,rvf,  &
                        xfref,xhref,yfref,yhref,sigma,c1,c2,gz,zs,        &
                        arh1,arh2,arf1,arf2,rgz,gzu,gzv,rds,rdsf,         &
                        zh,mh,rmh,zf,mf,rmf,rho0s,pi0s,prs0s,             &
                        pi0,prs0,rho0,thv0,th0,rth0,qv0,                  &
                        u0,v0,qc0,qi0,rr0,rf0,rrf0,                       &
                        rain,raincv,hail,sws,svs,sps,srs,sgs,sus,shs,     &
                        thflux,qvflux,cd,ch,cq,f2d,                       &
                        dum1,dum2,dum3,dum4,dum5,dum6,dum7,dum8,divx,rho,prs, &
                        rru,ua,u3d,uten,uten1,rrv,va,v3d,vten,vten1,      &
                        rrw,wa,w3d,wten,wten1,ppi,pp3d,ppten,sten,        &
                        tha,th3d,thten,thten1,qa,q3d,qten,                &
                        kmh,kmv,khh,khv,tkea,tke3d,tketen,                &
                        pta,pt3d,ptten,pdata,haildata,out3d)

#ifdef MPI
          call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif
        time_init3d = mytime()
        if(dowr) write(outfile,*)
        if(dowr) write(outfile,*) '  time spent in subroutine init3d = ',time_init3d,' seconds '
        if(dowr) write(outfile,*)

!----------------------------------------------------------------------

      nwritea = 1

      icrs = 1
      icenter = nx/2 + 1
      jcenter = ny/2 + 1
      xcenter = minx + 0.5*(maxx-minx)
      ycenter = miny + 0.5*(maxy-miny)

      IF( doazimavg )THEN

        ! check that ddr is positive and non-zero
        if( ddr.lt.smeps )then
          if(myid.eq.0)then
          print *
          print *,'  ddr  = ',ddr
          print *
          print *,'  Stopping cm1 ... '
          print *
          endif
#ifdef MPI
          call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif
          call stopcm1
        endif

        i = 1
        do while( xhref(nx/2+i).lt.rlen .and. i.lt.nx/2 )
          i = i+1
        enddo

        icrs = i

        ! check that icrs is greater than
        if( icrs.le.0 )then
          if(myid.eq.0)then
          print *
          print *,'  icrs = ',icrs
          print *
          print *,'  Stopping cm1 ... '
          print *
          endif
#ifdef MPI
          call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif
          call stopcm1
        endif
        ! initial values for cyclone center
        ! (assumed center of domain)
        if( myid.eq.0 ) print *
        if( myid.eq.0 ) print *,'  ddr,rlen,icrs,maxr = ',ddr,rlen,icrs,xhref(nx/2+icrs)
        if( myid.eq.0 ) print *
      ENDIF

      if( .not. allocated( crref ) ) allocate( crref(icrs) )
      crref = 0.0

      do i=1,icrs
        crref(i) = xhref(nx/2+i)
      enddo

!----------------------------------------------------------------------
!  Time-averaged vars

      ! time-avg parameters are set in param.F

      if( .not. allocated( ntavg ) ) allocate( ntavg(ntim) )
      ntavg = 0
      if( .not. allocated( rtavg ) ) allocate( rtavg(ntim) )
      rtavg = 0.0

      if( .not. allocated( tavg ) ) allocate( tavg(ibta:ieta,jbta:jeta,kbta:keta,ntim,ntavr) )
      tavg = 0.0
      if( .not. allocated( timavg ) ) allocate( timavg(ibta:ieta,jbta:jeta,kbta:keta,ntavr) )
      timavg = smeps

      if( .not. allocated( sfctavg ) ) allocate( sfctavg(ibta:ieta,jbta:jeta,ntim,nsfctavr) )
      sfctavg = 0.0
      if( .not. allocated( sfctimavg ) ) allocate( sfctimavg(ibta:ieta,jbta:jeta,nsfctavr) )
      sfctimavg = smeps

!----------------------------------------------------------------------

      call       init_physics(prs0,rf0,dum1,dum2,dum3,u0,ua,v0,va,o30,   &
                             lu_index,xland,emiss,thc,albd,znt,mavail,tsk,u1,v1,s1, &
                             zh,u10,v10,wspd,pbl2dint,ust,pbl2d)

      call init_surface(num_soil_layers,dosfcflx,xh,ruh,xf,yh,rvh,yf,   &
                        lu_index,xland,tsk,slab_zs,slab_dzs,tslb, &
                        emiss,thc,albd,znt,rznt,mavail,dsxy,prs0s,prs0,   &
                        tmn,tml,t0ml,hml,h0ml,huml,hvml,tmoml,pbl2d)

      IF(irst.eq.1)THEN

        startup = .false.
        restarted = .true.
        init_pbl = 0
        init_sfcmodel = 0
        noutcheck = .true.
        writeit = .false.
        if(timestats.ge.1)then
#ifdef MPI
          call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif
          time_rstread = mytime()
        endif

        call     restart_read(nstep,srec,sirec,urec,vrec,wrec,nrec,mrec,prec,      &
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
                              rho,prs,ua,uten,va,vten,wa,ppi,tha,qa,tkea,          &
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
#ifdef MPI
        call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif
        time_rstread = mytime()
        if(dowr) write(outfile,*)
        if(dowr) write(outfile,*) '  time spent in subroutine restart_read = ',time_rstread,' seconds '
        if(dowr) write(outfile,*)
        !  In case user wants to change values on a restart:
        IF( restart_reset_frqtim )THEN 
          if( statfrq.gt.1.0e-6 ) stattim = mtime + statfrq
          if(  tapfrq.gt.1.0e-6 ) taptim  = mtime + tapfrq
          if(  rstfrq.gt.1.0e-6 ) rsttim  = mtime + rstfrq
        ENDIF
        if( prcltim.gt.0.1d60 )then
          prcltim = mtime + prclfrq
        endif
        if( radopt.ge.1 .and. radtim.gt.0.1d60 )then
          radtim = mtime
          getrad =  .true.
        endif
        if(dowr) write(outfile,*)
        if(dowr) write(outfile,*) '  Using the following: '
        if(dowr) write(outfile,*)
        if(dowr) write(outfile,*) '   stattim = ',stattim
        if(dowr) write(outfile,*) '   taptim  = ',taptim
        if(dowr) write(outfile,*) '   rsttim  = ',rsttim
        if(dowr) write(outfile,*) '   radtim  = ',radtim
        if(dowr.and.iprcl.eq.1) write(outfile,*) '   prcltim = ',prcltim
        if(dowr) write(outfile,*)
        if( interp_on_restart )then
          restart_time = mtime
          if( dot2p ) spinup_t2p = .true.
        endif
        if(dowr) write(outfile,*) '   restart_time = ',restart_time
        if(dowr) write(outfile,*) '   spinup_t2p   = ',spinup_t2p
        if(dowr) write(outfile,*)

        ! Experimental (cm1r19.8) !
        ! Insert random perturbations to potential temperature.
        ! (see subroutine in misclibs.F for more details)
!!!        call randpert(xfref,yfref,xh,yh,zh,zf,tha)

      ELSE

        restarted = .false.
        restart_prcl = .false.
        restart_hail = .false.

      ENDIF

      !  If hail trajectories are on but the restart file came from a run
      !  with ihailtraj=0, then hailtrajtim in that file is the "never"
      !  sentinel (1.0d60) that param.F assigns when ihailtraj.le.0.
      !  Reset it so hail output starts from the restart time.
      IF( ihailtraj.eq.1 .and. restarted .and. (.not.restart_hail) )THEN
        hailtrajtim = mtime
        haildata(:,hailtime) = mtime
        if(dowr) write(outfile,*)
        if(dowr) write(outfile,*) '  Starting hail trajectories from a restart'
        if(dowr) write(outfile,*) '  file that had none.  hailtrajtim = ',hailtrajtim
        if(dowr) write(outfile,*)
      ENDIF

      !  Restarting with hailstone data: continue the existing hail output
      !  file instead of starting it again at record 1 (hrec is not saved in
      !  the restart file, so it is found from the output file itself).
      IF( ihailtraj.eq.1 .and. restarted .and. restart_hail )THEN
        call hailcast_restart_hrec(hailtrajtim,hrec)
#ifdef MPI
        call MPI_BCAST(hrec,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
#endif
      ENDIF


      IF( dodomaindiag )THEN
        IF( restarted .and. diagfrq.gt.1.0e-6 )THEN
          diagtim = mtime + diagfrq
          if(dowr) write(outfile,*)
          if(dowr) write(outfile,*) '   diagtim = ',diagtim
          if(dowr) write(outfile,*)
        ELSE
          diagtim = 0.0
        ENDIF
      ELSE
        diagtim = 1.0d60
      ENDIF

      IF( doazimavg )THEN
        IF( restarted .and. azimavgfrq.gt.1.0e-6 )THEN
          azimavgtim = mtime + azimavgfrq
          if(dowr) write(outfile,*)
          if(dowr) write(outfile,*) '   azimavgtim = ',azimavgtim
          if(dowr) write(outfile,*)
        ELSE
          azimavgtim = 0.0
        ENDIF
      ELSE
        azimavgtim = 1.0d60
      ENDIF

      IF( dohifrq )THEN
        IF( restarted .and. hifrqfrq.ge.1.0e-6 )THEN
          hifrqtim = mtime + hifrqfrq
          if(dowr) write(outfile,*)
          if(dowr) write(outfile,*) '   hifrqtim = ',hifrqtim
          if(dowr) write(outfile,*)
        ELSE
          hifrqtim = 0.0
        ENDIF
      ELSE
        hifrqtim = 1.0d60
      ENDIF


      IF( run_time .gt. 0.0 )THEN
        timax = mtime + run_time
        if(dowr) write(outfile,*)
        if(dowr) write(outfile,*) '  Detected positive value for run_time. '
        if(dowr) write(outfile,*) '     Using the following values: '
        if(dowr) write(outfile,*) '        run_time = ',run_time
        if(dowr) write(outfile,*) '        timax    = ',timax
        if(dowr) write(outfile,*)
      ENDIF


!----------------------------------------------------------------------
!  final batch of allocations...

#ifdef MPI
          call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif
          if(myid.eq.0) print *
          if(myid.eq.0) print *,'      allocating turb arrays ... '
#ifdef MPI
          call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif

      if( .not. allocated( t11 ) ) allocate(  t11(ib:ie,jb:je,kb:ke) )
      t11 = 0.0
      if( .not. allocated( t12 ) ) allocate(  t12(ib:ie,jb:je,kb:ke) )
      t12 = 0.0
      if( .not. allocated( t13 ) ) allocate(  t13(ib:ie,jb:je,kb:ke) )
      t13 = 0.0
      if( .not. allocated( t22 ) ) allocate(  t22(ib:ie,jb:je,kb:ke) )
      t22 = 0.0
      if( .not. allocated( t23 ) ) allocate(  t23(ib:ie,jb:je,kb:ke) )
      t23 = 0.0
      if( .not. allocated( t33 ) ) allocate(  t33(ib:ie,jb:je,kb:ke) )
      t33 = 0.0

      if( .not. allocated( m11 ) ) allocate( m11(ibnba:ienba,jbnba:jenba,kbnba:kenba) )
      m11 = 0.0
      if( .not. allocated( m12 ) ) allocate( m12(ibnba:ienba,jbnba:jenba,kbnba:kenba) )
      m12 = 0.0
      if( .not. allocated( m13 ) ) allocate( m13(ibnba:ienba,jbnba:jenba,kbnba:kenba) )
      m13 = 0.0
      if( .not. allocated( m22 ) ) allocate( m22(ibnba:ienba,jbnba:jenba,kbnba:kenba) )
      m22 = 0.0
      if( .not. allocated( m23 ) ) allocate( m23(ibnba:ienba,jbnba:jenba,kbnba:kenba) )
      m23 = 0.0
      if( .not. allocated( m33 ) ) allocate( m33(ibnba:ienba,jbnba:jenba,kbnba:kenba) )
      m33 = 0.0

      if( .not. allocated( nm ) ) allocate(      nm(ib:ie,jb:je,kb:ke+1) )
      nm = 0.0
      if( .not. allocated( defv ) ) allocate(    defv(ib:ie,jb:je,kb:ke+1) )
      defv = 0.0
      if( .not. allocated( defh ) ) allocate(    defh(ib:ie,jb:je,kb:ke+1) )
      defh = 0.0
      if( .not. allocated( lenscl ) ) allocate(  lenscl(ib:ie,jb:je,kb:ke+1) )
      lenscl = 0.0
      if( .not. allocated( dissten ) ) allocate( dissten(ib:ie,jb:je,kb:ke+1) )
      dissten = 0.0
      if( .not. allocated( epst ) ) allocate(    epst(ib:ie,jb:je,kb:ke+1) )
      epst = 0.0

      if( .not. allocated( xkzh ) ) allocate(   xkzh(ibb:ieb,jbb:jeb,kbb:keb) )
      xkzh = 0.0
      if( .not. allocated( xkzq ) ) allocate(   xkzq(ibb:ieb,jbb:jeb,kbb:keb) )
      xkzq = 0.0
      if( .not. allocated( xkzm ) ) allocate(   xkzm(ibb:ieb,jbb:jeb,kbb:keb) )
      xkzm = 0.0

      if( .not. allocated( u2pt ) ) allocate( u2pt(ib2pt:ie2pt,jb2pt:je2pt,kb2pt:ke2pt) )
      u2pt = 0.0
      if( .not. allocated( v2pt ) ) allocate( v2pt(ib2pt:ie2pt,jb2pt:je2pt,kb2pt:ke2pt) )
      v2pt = 0.0
      if( .not. allocated( kmwk ) ) allocate( kmwk(ib2pt:ie2pt,jb2pt:je2pt,kb2pt:ke2pt) )
      kmwk = 0.0
      if( .not. allocated( ufwk ) ) allocate( ufwk(ib2pt:ie2pt,jb2pt:je2pt,kb2pt:ke2pt) )
      ufwk = 0.0
      if( .not. allocated( vfwk ) ) allocate( vfwk(ib2pt:ie2pt,jb2pt:je2pt,kb2pt:ke2pt) )
      vfwk = 0.0

!----------------------------------------------------------------------

      if( myid.eq.0 ) print *
      if( myid.eq.0 ) print *,'  allocating diagnostic arrays: '
      if( myid.eq.0 ) print *

#ifdef MPI
          call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif
          if(myid.eq.0) print *,'      tdiag ... '
#ifdef MPI
          call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif
      if( .not. allocated( tdiag ) ) allocate( tdiag(ibdt:iedt,jbdt:jedt,kbdt:kedt,ntdiag) )
      tdiag = 0.0
#ifdef MPI
          call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif
          if(myid.eq.0) print *,'      qdiag ... '
#ifdef MPI
          call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif
      if( .not. allocated( qdiag ) ) allocate( qdiag(ibdq:iedq,jbdq:jedq,kbdq:kedq,nqdiag) )
      qdiag = 0.0
#ifdef MPI
          call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif
          if(myid.eq.0) print *,'      udiag ... '
#ifdef MPI
          call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif
      if( .not. allocated( udiag ) ) allocate( udiag(ibdv:iedv,jbdv:jedv,kbdv:kedv,nudiag) )
      udiag = 0.0
#ifdef MPI
          call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif
          if(myid.eq.0) print *,'      vdiag ... '
#ifdef MPI
          call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif
      if( .not. allocated( vdiag ) ) allocate( vdiag(ibdv:iedv,jbdv:jedv,kbdv:kedv,nvdiag) )
      vdiag = 0.0
#ifdef MPI
          call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif
          if(myid.eq.0) print *,'      wdiag ... '
#ifdef MPI
          call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif
      if( .not. allocated( wdiag ) ) allocate( wdiag(ibdv:iedv,jbdv:jedv,kbdv:kedv,nwdiag) )
      wdiag = 0.0
#ifdef MPI
          call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif
          if(myid.eq.0) print *,'      kdiag ... '
#ifdef MPI
          call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif
      if( .not. allocated( kdiag ) ) allocate( kdiag(ibdk:iedk,jbdk:jedk,kbdk:kedk,nkdiag) )
      kdiag = 0.0
#ifdef MPI
          call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif
          if(myid.eq.0) print *,'      pdiag ... '
#ifdef MPI
          call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif
      if( .not. allocated( pdiag ) ) allocate( pdiag(ibdp:iedp,jbdp:jedp,kbdp:kedp,npdiag) )
      pdiag = 0.0
#ifdef MPI
          call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif
          if(myid.eq.0) print *
#ifdef MPI
          call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif


!----------------------------------------------------------------------
!  Eddy recycling:

#ifdef MPI
          call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif
          if(myid.eq.0) print *
          if(myid.eq.0) print *,'      allocating recy arrays ... '
#ifdef MPI
          call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif

      if( .not. allocated( recy_cap_we ) ) allocate( recy_cap_we(ib:ie,jb:je) )
      recy_cap_we = 0
      if( .not. allocated( recy_cap_sn ) ) allocate( recy_cap_sn(ib:ie,jb:je) )
      recy_cap_sn = 0

      if( .not. allocated( recy_inj_we ) ) allocate( recy_inj_we(ib:ie,jb:je) )
      recy_inj_we = 0
      if( .not. allocated( recy_inj_sn ) ) allocate( recy_inj_sn(ib:ie,jb:je) )
      recy_inj_sn = 0

      if( .not. allocated( recy_wsp_w ) ) allocate( recy_wsp_w(ib:ie,jb:je) )
      recy_wsp_w = 0.0
      if( .not. allocated( recy_wsp_e ) ) allocate( recy_wsp_e(ib:ie,jb:je) )
      recy_wsp_e = 0.0
      if( .not. allocated( recy_wsp_s ) ) allocate( recy_wsp_s(ib:ie,jb:je) )
      recy_wsp_s = 0.0
      if( .not. allocated( recy_wsp_n ) ) allocate( recy_wsp_n(ib:ie,jb:je) )
      recy_wsp_n = 0.0

      if( .not. allocated( recy_tau_w ) ) allocate( recy_tau_w(ib:ie,jb:je) )
      recy_tau_w = 0.0
      if( .not. allocated( recy_tau_e ) ) allocate( recy_tau_e(ib:ie,jb:je) )
      recy_tau_e = 0.0
      if( .not. allocated( recy_tau_s ) ) allocate( recy_tau_s(ib:ie,jb:je) )
      recy_tau_s = 0.0
      if( .not. allocated( recy_tau_n ) ) allocate( recy_tau_n(ib:ie,jb:je) )
      recy_tau_n = 0.0

      if( .not. allocated( recywe ) ) allocate( recywe(irecywe,jrecywe,krecy,nrecy) )
      recywe = 0.0
      if( .not. allocated( recysn ) ) allocate( recysn(irecysn,jrecysn,krecy,nrecy) )
      recysn = 0.0

!----------------------------------------------------------------------

#ifdef MPI
          call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif
          if(myid.eq.0) print *,'  ... done with major memory allocations '
          if(myid.eq.0) print *
#ifdef MPI
          call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif


!----------------------------------------------------------------------
!  Get ready for time integration:

      if(myid.eq.0) print *
#ifdef MPI
      call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif
      time_getset = mytime()

      call       getset(restarted,mass1,cloudvar,dt,   &
                        xh,rxh,arh1,arh2,uh,ruh,xf,rxf,arf1,arf2,uf,ruf,yh,vh,rvh,yf,vf,rvf, &
                        zs,gz,rgz,gzu,rgzu,gzv,rgzv,dzdx,dzdy,gx,gxu,gy,gyv,    &
                        rds,sigma,rdsf,sigmaf,zh,mh,rmh,c1,c2,zf,mf,rmf,dumk1,dumk2, &
                        pi0,th0,rho0,rr0,rf0,rrf0,prs0,thv0,qv0,ust,znt,z1,u1,v1,s1,cm0,eimaskw,eimaske,eimasks,eimaskn,u0,v0,rth0,divx,rho,rr,rf,prs, &
                        dum1,dum2,dum3,dum4,dum5,dum6,dum7,dum8,                &
                        rru,ua,u3d,uten,uten1,rrv,va,v3d,vten,vten1,rrw,wa,w3d,wten,wten1,ppi,pp3d,ppx,phi1,phi2,  &
                        tha,th3d,qa,q3d,cme,csm,ce1,ce2,                        &
                        tkea,tke3d,kmh,kmv,khh,khv,nm,defv,defh,lenscl,dissten, &
                        t11,t12,t13,t22,t23,t33,flag  ,                         &
                        pta,pt3d,pdata,haildata,tavg,timavg,pbl3d,pbl2d,pbl2dint, &
                        cfb,cfa,cfc, d1, d2,pdt,lgbth,lgbph,                    &
                        pgrid1,cgrid1,cgrid2,psenda,precva,psendb,precvb,psendc,precvc,psendd,precvd,rhsx,transx,rhsy,transy,  &
                        myi1p,myi2p,myj1p,myj2p,                                &
                        myi1pp,myi2pp,myj1pp,myj2pp,                            &
                        npcomm,nplist,npi1,npi2,npj1,npj2,npk1,npk2,            &
                        upten,vpten,ucuten,vcuten,out3d,                        &
                        bndy,kbdy,hflxw,hflxe,hflxs,hflxn,                      &
                        reqs_u,reqs_v,reqs_w,reqs_s,reqs_p,reqs_tk,             &
                        nw1,nw2,ne1,ne2,sw1,sw2,se1,se2,                        &
                        pw1,pw2,pe1,pe2,ps1,ps2,pn1,pn2,                        &
                        uw31,uw32,ue31,ue32,us31,us32,un31,un32,                &
                        vw31,vw32,ve31,ve32,vs31,vs32,vn31,vn32,                &
                        ww31,ww32,we31,we32,ws31,ws32,wn31,wn32,                &
                        sw31,sw32,se31,se32,ss31,ss32,sn31,sn32,                &
                        tkw1,tkw2,tke1,tke2,tks1,tks2,tkn1,tkn2,                &
                        kw1,kw2,ke1,ke2,ks1,ks2,kn1,kn2,                        &
                        n3w1,n3w2,n3e1,n3e2,s3w1,s3w2,s3e1,s3e2)    ! end_getset

#ifdef MPI
          call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif
        time_getset = mytime()
        if(dowr) write(outfile,*)
        if(dowr) write(outfile,*) '  time spent in subroutine getset = ',time_getset,' seconds '
        if(dowr) write(outfile,*)

!----------------------------------------------------------------------

      if( ibalance.eq.2 .and. psolver.ne.4 .and. psolver.ne.5 .and. (.not.pdcomp) )then

        if( myid.eq.0 ) print *
        if( myid.eq.0 ) print *,'  deallocating elliptic pressure arrays: '
        if( myid.eq.0 ) print *

        deallocate( cfb )
        deallocate( cfa )
        deallocate( cfc )
        deallocate( d1 )
        deallocate( d2 )
        deallocate( pdt )
        deallocate( lgbth )
        deallocate( lgbph )

        ipb = 1
        ipe = 1
        jpb = 1
        jpe = 1
        kpb = 1
        kpe = 1
        if( .not. allocated( cfb ) ) allocate(    cfb(ipb:ipe,jpb:jpe,kpb:kpe) )
        cfb = 0.0
        if( .not. allocated( cfa ) ) allocate(    cfa(kpb:kpe) )
        cfa = 0.0
        if( .not. allocated( cfc ) ) allocate(    cfc(kpb:kpe) )
        cfc = 0.0
        if( .not. allocated( d1 ) ) allocate(     d1(kpb:kpe) )
        d1 = 0.0
        if( .not. allocated( d2 ) ) allocate(     d2(kpb:kpe) )
        d2 = 0.0
        if( .not. allocated( pdt ) ) allocate(    pdt(ipb:ipe,jpb:jpe,kpb:kpe) )
        pdt = 0.0
        if( .not. allocated( lgbth ) ) allocate(  lgbth(ipb:ipe,jpb:jpe,kpb:kpe) )
        lgbth = 0.0
        if( .not. allocated( lgbph ) ) allocate(  lgbph(ipb:ipe,jpb:jpe,kpb:kpe) )
        lgbph = 0.0

        deallocate( pgrid1 )
        deallocate( cgrid1 )
        deallocate( cgrid2 )
        deallocate( psenda )
        deallocate( precva )
        deallocate( psendb )
        deallocate( precvb )
        deallocate( psendc )
        deallocate( precvc )
        deallocate( rhsx )
        deallocate( rhsy )
        deallocate( transx )
        deallocate( transy )

        pi1a = 1
        pi2a = 1
        pj1a = 1
        pj2a = 1
        pk1a = 1
        pk2a = 1

        pi1b = 1
        pi2b = 1
        pj1b = 1
        pj2b = 1
        pk1b = 1
        pk2b = 1

        pimaxa = 1
        pjmaxa = 1
        pkmaxa = 1
        pnmaxa = 1

        pimaxb = 1
        pjmaxb = 1
        pkmaxb = 1
        pnmaxb = 1

        pimaxc = 1
        pjmaxc = 1
        pkmaxc = 1
        pnmaxc = 1

        if( .not. allocated( pgrid1 ) ) allocate( pgrid1(pi1a:pi2a,pj1a:pj2a,pk1a:pk2a) )
        pgrid1 = 0.0
        if( .not. allocated( cgrid1 ) ) allocate( cgrid1(pi1a:pi2a,pj1a:pj2a,pk1a:pk2a) )
        cgrid1 = 0.0
        if( .not. allocated( cgrid2 ) ) allocate( cgrid2(pi1b:pi2b,pj1b:pj2b,pk1b:pk2b) )
        cgrid2 = 0.0
        if( .not. allocated( psenda ) ) allocate( psenda(pimaxa,pjmaxa,pkmaxa,pnmaxa) )
        psenda = 0.0
        if( .not. allocated( precva ) ) allocate( precva(pimaxa,pjmaxa,pkmaxa,pnmaxa) )
        precva = 0.0
        if( .not. allocated( psendb ) ) allocate( psendb(pimaxb,pjmaxb,pkmaxb,pnmaxb) )
        psendb = 0.0
        if( .not. allocated( precvb ) ) allocate( precvb(pimaxb,pjmaxb,pkmaxb,pnmaxb) )
        precvb = 0.0
        if( .not. allocated( psendc ) ) allocate( psendc(pimaxc,pjmaxc,pkmaxc,pnmaxc) )
        psendc = 0.0
        if( .not. allocated( precvc ) ) allocate( precvc(pimaxc,pjmaxc,pkmaxc,pnmaxc) )
        precvc = 0.0
        if( .not. allocated( rhsx ) ) allocate(   rhsx(pi1b:pi2b) )
        rhsx = 0.0
        if( .not. allocated( transx ) ) allocate( transx(pi1b:pi2b) )
        transx = 0.0
        if( .not. allocated( rhsy ) ) allocate(   rhsy(pj1a:pj2a) )
        rhsy = 0.0
        if( .not. allocated( transy ) ) allocate( transy(pj1a:pj2a) )
        transy = 0.0
      endif

!----------------------------------------------------------------------
!  get info for large-scale nudging:
!  (190724, GHB: moved this after restart_read)

      IF( do_lsnudge )THEN

        if(dowr) write(outfile,*)
        if(dowr) write(outfile,*)
        if(dowr) write(outfile,*) '  LSNUDGE is active '

        if( myid.eq.0 )then
          print *
          print *,'  do_lsnudge_u   = ',do_lsnudge_u
          print *,'  do_lsnudge_v   = ',do_lsnudge_v
          print *,'  do_lsnudge_th  = ',do_lsnudge_th
          print *,'  do_lsnudge_qv  = ',do_lsnudge_qv
          print *
        endif

        call   read_lsnudge(lsnudge_u,lsnudge_v,lsnudge_th,lsnudge_qv,sngl(mtime),zh)

        if( myid.eq.0 ) print *
        if( myid.eq.0 ) print *,'  lsnudge_tau    =  ',lsnudge_tau
        if( myid.eq.0 ) print *,'  lsnudge_start  =  ',lsnudge_start
        if( myid.eq.0 ) print *,'  lsnudge_end    =  ',lsnudge_end
        if( myid.eq.0 ) print *

      ENDIF


!----------------------------------------------------------------------
!  Get starting timestep (if adaptive timestepping is used)


    get_adapt_dt:  &
    if( adapt_dt.eq.1 )then

      dt_case:  &
      if( interp_on_restart )then
        ! adaptive dt if interpolating on restart:

        !..........  check for advection stability ..........!

        call calccflquick(dt,uh,vh,mh,u3d,v3d,w3d,reqc)
#ifdef MPI
        call mpi_wait(reqc,mpi_status_ignore,ierr)
#endif

        ! a "target" Courant number of 0.5 seems safe:

        cfl_target = 0.5

        dbldt = dt*cfl_target/max(1.0e-12,cflmax)
        dt = dbldt

        call calccflquick(dt,uh,vh,mh,u3d,v3d,w3d,reqc)
#ifdef MPI
        call mpi_wait(reqc,mpi_status_ignore,ierr)
#endif

        if( myid.eq.0 ) print *,'  Adaptive timestep for interp_on_restart: '
        if( myid.eq.0 ) print *,'    after cflcheck: cflmax,dt = ',cflmax,dt

        !..........  check for diffusion stability ..........!

        IF( cm1setup.ge.1 )THEN

          call calcksquick(dt,uh,vh,mf,kmh,kmv,khh,khv,reqk,cm0,xh,yh,zh)
#ifdef MPI
          call mpi_wait(reqk,mpi_status_ignore,ierr)
#endif

          ks_target = 0.05

          ! only do this if dt from above is too large:
          if( ksmax.gt.ks_target )then
            dbldt = dbldt*ks_target/max(1.0e-12,ksmax)
            dt = dbldt
          endif

          call calcksquick(dt,uh,vh,mf,kmh,kmv,khh,khv,reqk,cm0,xh,yh,zh)
#ifdef MPI
          call mpi_wait(reqk,mpi_status_ignore,ierr)
#endif

          if( myid.eq.0 ) print *,'    after kscheck:  ksmax,dt  = ',ksmax,dt

          ! this seems to work well: use half of dt required for stability, to begin
          dbldt = 0.5*dbldt
          dt = dbldt

          adtlast = dbldt

        ENDIF

        !..........  should be ready to go ..........!

      else  dt_case

        dt = dbldt

        if( .not. restarted )then

          adtlast  = dbldt

          ! adaptive dt if starting CM1 for first time 
          ! (ie, not a restart)

          ! get cflmax and ksmax:

          call calccflquick(dt,uh,vh,mh,u3d,v3d,w3d,reqc)
#ifdef MPI
          call mpi_wait(reqc,mpi_status_ignore,ierr)
          if(timestats.ge.1) time_cflq=time_cflq+mytime()
#endif
          if(myid.eq.0) print *,'  cflmax = ',cflmax
          if( cm1setup.ge.1 )then
            call calcksquick(dt,uh,vh,mf,kmh,kmv,khh,khv,reqk,cm0,xh,yh,zh)
#ifdef MPI
            call mpi_wait(reqk,mpi_status_ignore,ierr)
            if(timestats.ge.1) time_cflq=time_cflq+mytime()
#endif
            if(myid.eq.0) print *,'  ksmax  = ',ksmax
          endif
          stopit = .false.
          stopcondition = 0
          ndt  = 0
          adt  = 0.0
          acfl = 0.0
          call   getnewdt(ndt,dt,dtlast,adt,acfl,dbldt,                                 &
                          mtime,stattim,taptim,rsttim,prcltim,hailtrajtim,     &
                          diagtim,azimavgtim,                                   &
                          dorestart,dowriteout,dostat,doprclout,dohailout,      &
                          dotdwrite,doazimwrite,                                &
                          hifrqtim,dohifrqwrite,.true.)

        endif

      endif  dt_case

    endif  get_adapt_dt

!----------------------------------------------------------------------
!  All done with initialization and setup.  A few more odds and ends ....

      ! cm1r19:  nsound is always determined diagnostically
      call dtsmall(dt,dbldt)

!!!      if(dowr) write(outfile,*)
!!!      if(dowr) write(outfile,*) '  nsound = ',nsound
!!!      if(dowr) write(outfile,*)


      if( testcase.eq.17 )then
        ! only surface is moved for testcase=17
        umove = 0.0
        vmove = 0.0
        imove = imovesave
        umovesfc = umovesave
        vmovesfc = vmovesave
      endif


      !  Prepare for I/O:

      cmpr_output = .true.


      call       setup_stat_vars(name_stat,desc_stat,unit_stat,  &
                                 qname,qunit,budname)
      if(dowr) write(outfile,*) 'stat_out = ',stat_out
      if(dowr) write(outfile,*)
      if( .not. allocated( rstat ) ) allocate( rstat(stat_out) )
      rstat = 0.0
      noutstat = stat_out


      if( iprcl.eq.1 )then
        call       setup_parcel_vars(name_prcl,desc_prcl,unit_prcl,qname)
        if(dowr) write(outfile,*) 'prcl_out = ',prcl_out
        if(dowr) write(outfile,*)
      endif
      noutparcel = prcl_out

      if( ihailtraj.eq.1 ) THEN
        if(myid.eq.0) THEN
            print *, 'starting setup_hailcast_vars'
        endif
        call       setup_hailcast_vars(name_hail,desc_hail,unit_hail,qname)
        !if(myid.eq.0)  print *, 'after setup_hailcast_vars hail_out:', hail_out
      endif


      if( axisymm.eq.0 .and. doazimavg .or. do_adapt_move )  &
      call       getrangle(xh,yh,xcenter,ycenter,cir,crr,cangle,icrs,crref,xhref,xfref)


      if(dowr) write(outfile,*)
      if(dowr) write(outfile,*) '-------------Done with Preprocessors-----------'
      if(dowr) write(outfile,*)


!----------------------------------------------------------------------

#ifdef MPI
      call MPI_BARRIER (MPI_COMM_WORLD,ierr)

      if(myid.eq.0)then
        tstart=mpi_wtime()
      endif
#endif

      call set_time_to_zero()

!----------------------------------------------------------------------
!  Time loop

      if(timestats.ge.1)then
        steptime1 = 0.0
        steptime2 = 0.0
      endif


      nstep0 = nstep


      IF( startup )THEN
        ! starting from t=0:  skip solve and write initial conditions
        dosolve     = .false.
        dorad       = .true.
        dostat      = .true.
        dowriteout  = .true.
        doprclout   = .true.
        dohailout   = .true.
        dotdwrite   = .true.
        doazimwrite = .true.
        dohifrqwrite = .true.
        IF( imoist.eq.1 )THEN
          ! all this is copied from mp_driver
          if( getdbz .and. qd_dbz.le.0 )then
            getdbz = .false.
          endif
          getvt     = .false.
          getsed    = .false.
          getqdiags = .false.
          IF( dowriteout .and. (ptype.eq.3.or.ptype.eq.5) )THEN
            if( output_fallvel.eq.1 ) getvt     = .true.
            if( output_sedtend.eq.1 ) getsed    = .true.
            if( output_qdiags.eq.1 .and. ptype.eq.5 ) getqdiags = .true.
          ENDIF
          IF( dotdwrite .and. (ptype.eq.3.or.ptype.eq.5) )THEN
            getvt     = .true.
            getsed    = .true.
            if( ptype.eq.5 ) getqdiags = .true.
          ENDIF
          IF( doazimwrite .and. (ptype.eq.3.or.ptype.eq.5) )THEN
            getvt     = .true.
            getsed    = .true.
            if( ptype.eq.5 ) getqdiags = .true.
          ENDIF
        ENDIF
      ELSE
        ! call solve
        dosolve     = .true.
        if( getrad )then
          dorad = .true.
        else
          dorad = .false.
        endif
        dostat      = .false.
        dowriteout  = .false.
        doprclout   = .false.
        dohailout   = .false.
        dotdwrite   = .false.
        doazimwrite = .false.
        dohifrqwrite = .false.
        if( .not. restart_prcl ) doprclout = .true.
        if( ihailtraj.eq.1 .and. (.not.restart_hail) ) dohailout = .true.
      ENDIF


      IF( abs(rstfrq).lt.1.0e-10 )THEN
        ! 150820:  Write restart file if rstfrq=0 and stop
        dorestart = .true.
        dosolve = .false.
        nrst = 0
      ELSE
        dorestart = .false.
      ENDIF

      IF( restarted )THEN
        dostat = .true.
        dosolve = .false.
        if( doazimavg .and. nwritea.eq.1 ) doazimwrite = .true.
        if( dodomaindiag .and. nwritet.eq.1 ) dotdwrite = .true.
        if( dohifrq .and. nwriteh.eq.1 ) dohifrqwrite = .true.
      ENDIF

      if( interp_on_restart .or. ( restarted .and. irandprst.eq.1 ) )then

        ! nconfig is updated for interp_on_restart
        nconfig = nconfig+1
        noutcheck = .false.
        dowriteout = .true.
        taptim = mtime
        nstatout = 1
        nrec = 1
        nwrite = 1
        srec = 1
        sirec = 1
        urec = 1
        vrec = 1
        wrec = 1
        xrec = 1
        yrec = 1
        zrec = 1
        mrec = 1
        prec = 1
        nwritet = 1
        trecs = 1
        trecw = 1
        nwriteh = 1
        nwritea = 1
        arecs = 1
        arecw = 1
        nwritemv = 1
        mvrec = 1
        writeit = .true.

        if( doazimavg )then
          doazimwrite = .true.
          azimavgtim = mtime
        endif
        if( dodomaindiag )then
          dotdwrite = .true.
          diagtim = mtime
        endif
        if( dohifrq )then
          dohifrqwrite = .true.
          hifrqtim = mtime
        endif

        if(myid.eq.0) print *
        if(myid.eq.0) print *,'  write files on restart: '
        if(myid.eq.0) print *,'  dowriteout,nwrite,taptim       = ',dowriteout,nwrite,taptim
        if(myid.eq.0) print *,'  doazimwrite,nwritea,azimavgtim = ',doazimwrite,nwritea,azimavgtim
        if(myid.eq.0) print *

      endif


    !c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c!
    timeloop:  &
    DO WHILE( mtime.lt.timax )

      if( noutcheck .and. myid.eq.0 )then
        print *
        print *,'  -------------------------------------------------------------  '
        print *,'                begin noutcheck                                  '
        print *,'                writeit = ',writeit
        print *
      endif

      ifsolve:  &
      if( dosolve )then

        ! call solve, integrate model one time step:
        nstep = nstep + 1


        ! check if writeout will be called at end of timestep:
        !   (needed for some diagnostic calculations in solve)
        dowriteout = .false.
!!!    print *,'  t1,t2 = ',mtime+dbldt,taptim-dtfrac*dbldt,(mtime+dbldt)-(taptim-dtfrac*dbldt)
        if( (mtime+dbldt).ge.(taptim-dtfrac*dbldt) .or. tapfrq.lt.0.0 ) dowriteout = .true.
        IF( myid.eq.0 .and. dowriteout ) print *,'  dowriteout = ',dowriteout


        dostat = .false.
        if( (mtime+dbldt).ge.(stattim-dtfrac*dbldt) .or. statfrq.lt.0.0 )  dostat     =  .true.


        dorestart  =  .false.
        if( (mtime+dbldt).ge.(rsttim-dtfrac*dbldt) .and. rstfrq.gt.0.0 ) dorestart  =  .true.
        IF( myid.eq.0 .and. dorestart ) print *,'  dorestart = ',dorestart


        doprclout  =  .false.
        if( iprcl.eq.1 )then
          if( (mtime+dbldt).ge.(prcltim-dtfrac*dbldt) .or. prclfrq.lt.0.0 ) doprclout  =  .true.
        endif

        dohailout  =  .false.
        if( ihailtraj.eq.1 )then
          if( (mtime+dbldt).ge.(hailtrajtim-dtfrac*dbldt) .or. hailfrq.lt.0.0 ) dohailout  =  .true.
        endif


        dorad = .false.
        IF( radopt.ge.1 )THEN
          ! use time at end of timestep:
          IF( (mtime+dbldt).ge.(radtim-dtfrac*dbldt) ) dorad = .true.
          IF( myid.eq.0 .and. dorad ) print *,'  dorad = ',dorad
        ENDIF


        getdbz = .false.
        IF(output_dbz.eq.1)THEN
          if( ((mtime+dbldt).ge.(taptim-dtfrac*dbldt)) .or. tapfrq.le.0.0 )then
            getdbz = .true.
          endif
          if( dodomaindiag )then
            if( (mtime+dbldt).ge.(diagtim-dtfrac*dbldt) .or. diagfrq.le.0.0 )then
              getdbz = .true.
            endif
          endif
          if( doazimavg )then
            if( (mtime+dbldt).ge.(azimavgtim-dtfrac*dbldt) .or. azimavgfrq.le.0.0 )then
              getdbz = .true.
            endif
          endif
          if( dohifrq )then
            if( (mtime+dbldt).ge.(hifrqtim-dtfrac*dbldt) .or. hifrqfrq.le.0.0 )then
              getdbz = .true.
            endif
          endif
        ENDIF
        IF( restart_file_dbz )THEN
          if( ((mtime+dbldt).ge.(rsttim-dtfrac*dbldt)) .and. rstfrq.gt.0.0001 )then
            getdbz = .true.
          endif
        ENDIF
        IF( iprcl.eq.1 .and. prcl_dbz.eq.1 )THEN
          if( ((mtime+dbldt).ge.(prcltim-dtfrac*dbldt)) .or. prclfrq.le.0.0 )then
            getdbz = .true.
          endif
        ENDIF
        if(getdbz)then
          if(dowr) write(outfile,*) '  Getting dbz ... '
        endif


        IF( dodomaindiag )THEN
          if( (mtime+dbldt).ge.(diagtim-dtfrac*dbldt) .or. diagfrq.le.0.0 )then
            dotdwrite = .true.
          endif
        ENDIF

        IF( doazimavg )THEN
          if( (mtime+dbldt).ge.(azimavgtim-dtfrac*dbldt) .or. azimavgfrq.le.0.0 )then
            doazimwrite = .true.
          endif
        ENDIF

        IF( dohifrq )THEN
          if( (mtime+dbldt).ge.(hifrqtim-dtfrac*dbldt) .or. hifrqfrq.le.0.0 )then
            dohifrqwrite = .true.
          endif
        ENDIF

        getvt     = .false.
        getsed    = .false.
        getqdiags = .false.

        dotbud = .false.
        doqbud = .false.
        doubud = .false.
        dovbud = .false.
        dowbud = .false.

        IF( dowriteout .or. dotdwrite .or. doazimwrite )THEN
          if( output_thbudget.ge.1 .or. dodomaindiag )THEN
            dotbud = .true.
            if(myid.eq.0) print *,'  dotbud = ',dotbud
          endif
          if( output_qvbudget.ge.1 .or. dodomaindiag )THEN
            doqbud = .true.
            if(myid.eq.0) print *,'  doqbud = ',doqbud
          endif
          if( output_fallvel.eq.1 )THEN
            doqbud = .true.
            if(myid.eq.0) print *,'  doqbud (fallvel) = ',doqbud
          endif
          if( output_sedtend.ge.1 )THEN
            doqbud = .true.
            if(myid.eq.0) print *,'  doqbud (sedtend) = ',doqbud
          endif
          if( output_ubudget.ge.1 .or. dodomaindiag .or. pdcomp )THEN
            doubud = .true.
            if(myid.eq.0) print *,'  doubud = ',doubud
          endif
          if( output_vbudget.ge.1 .or. dodomaindiag .or. pdcomp )THEN
            dovbud = .true.
            if(myid.eq.0) print *,'  dovbud = ',dovbud
          endif
          if( output_wbudget.ge.1 .or. dodomaindiag .or. pdcomp )THEN
            dowbud = .true.
            if(myid.eq.0) print *,'  dowbud = ',dowbud
          endif
        ENDIF

        donudge = .false.

        if( adapt_dt.eq.1 .and. myid.eq.0 )then
          if( dt.gt.1.0e-3 )then
            write(6,122) cflmax,ksmax,dt,nsound
122         format(1x,'cflmax,ksmax,dt,nsound:',2x,f6.4,2x,f6.4,2x,f9.4,2x,i3)
          else
            write(6,123) cflmax,ksmax,dt,nsound
123         format(1x,'cflmax,ksmax,dt,nsound:',2x,f6.4,2x,f6.4,2x,ES11.4E2,2x,i3)
          endif
        endif

        if(timestats.ge.1) time_misc=time_misc+mytime()


        ! solve1: pre-RK terms
        ! (tendencies held fixed during Runge-Kutta integration and acoustic sub-stepping)
        call     solve1(nstep,num_soil_layers,                       &
                   dt,dtlast,adtlast,mtime,dbldt,                    &
                   dosfcflx,bud,bud2,qbudget,asq,bsq,                &
                   xh,rxh,arh1,arh2,uh,ruh,xf,rxf,arf1,arf2,uf,ruf,  &
                   yh,vh,rvh,yf,vf,rvf,                              &
                   rds,sigma,rdsf,sigmaf,                            &
                   tauh,taus,zh,mh,rmh,c1,c2,tauf,zf,mf,rmf,         &
                   wprof,ufrc,vfrc,thfrc,qvfrc,ug,vg,dvdr,           &
                   uavg,vavg,thavg,pavg,ulspg,vlspg,qavg,cavg,       &
                   pi0,rho0,prs0,thv0,th0,rth0,qv0,qc0,              &
                   qi0,rr0,rf0,rrf0,thrd,                            &
                   zs,gz,rgz,gzu,rgzu,gzv,rgzv,dzdx,dzdy,gx,gxu,gy,gyv, &
                   znt,zntmp,ust,thflux,qvflux,cm0,eimaskw,eimaske,eimasks,eimaskn,           &
                   radbcw,radbce,radbcs,radbcn,                      &
                   dum1,dum2,dum3,dum4,dum5,dum6,dum7,dum8,dum9,     &
                   divx,rho,rr,rf,prs,                               &
                   t11,t12,t13,t22,t23,t33,                          &
                   m11,m12,m13,m22,m23,m33,                          &
                   u0,ua,u3d,uten1,                                  &
                   v0,va,v3d,vten1,                                  &
                   rrw,wa,w3d,wten,wten1,                            &
                   ppi,pp3d,ppten,sten,sadv,                         &
                   tha,th3d,thten1,                                  &
                   qpten,qtten,qvten,qcten,qa,q3d,qten,              &
                   kmh,kmv,khh,khv,tkea,tke3d,tketen,                &
                   nm,defv,defh,lenscl,dissten,                      &
                   thpten,qvpten,qcpten,qipten,upten,vpten,qnipten,qncpten,o30,zir,  &
                   swten,lwten,thraten,                              &
                   thcuten,qvcuten,qccuten,qicuten,ucuten,vcuten,    &
                   pta,pt3d,ptten,                                   &
                   gamk,gamwall,kmw,ufw,vfw,u1b,v1b,u2pt,v2pt,ufwk,vfwk, &
                   tdiag,qdiag,udiag,vdiag,wdiag,kdiag,              &
                   out2d,out3d,                                      &
                   recy_cap_we,recy_cap_sn,recy_inj_we,recy_inj_sn,  &
                   recy_wsp_w,recy_wsp_e,recy_wsp_s,recy_wsp_n,recy_tau_w,recy_tau_e,recy_tau_s,recy_tau_n,  &
                   recywe,recysn,                                    &
                   lsnudge_u,lsnudge_v,lsnudge_th,lsnudge_qv,        &
                   timavg,                                           &
                   dowriteout,dorad,getdbz,dotdwrite,                &
                   dotbud,doqbud,doubud,dovbud,dowbud,donudge,       &
                   doazimwrite,dorestart)
        ! end_solve1


        ! solve2: RK loop and pressure solver
        ! (advection, buoyancy, pressure gradient)
        call     solve2(nstep,                                       &
                   dt,dtlast,mtime,dbldt,mass1,mass2,                &
                   dosfcflx,qmag,bud,bud2,qbudget,asq,bsq,           &
                   xh,rxh,arh1,arh2,uh,ruh,xf,rxf,arf1,arf2,uf,ruf,  &
                   yh,vh,rvh,yf,vf,rvf,                              &
                   dumk1,dumk2,dumk3,dumk4,rds,sigma,rdsf,sigmaf,    &
                   zh,mh,rmh,c1,c2,zf,mf,rmf,wprof,                  &
                   pi0,rho0,prs0,thv0,th0,rth0,qv0,qc0,              &
                   qi0,rr0,rf0,rrf0,                                 &
                   zs,gz,rgz,gzu,rgzu,gzv,rgzv,dzdx,dzdy,gx,gxu,gy,gyv,f2d,cm0, &
                   radbcw,radbce,radbcs,radbcn,dtu,dtu0,dtv,dtv0,    &
                   dum1,dum2,dum3,dum4,dum5,dum6,dum7,dum8,dum9,     &
                   divx,rho,rr,rf,prs,                               &
                   t11,t12,t13,t22,t23,t33,                          &
                   u0,rru,ua,u3d,uten,uten1,                         &
                   v0,rrv,va,v3d,vten,vten1,                         &
                   rrw,wa,w3d,wten,wten1,                            &
                   ppi,pp3d,ppten,sten,sadv,ppx,phi1,phi2,           &
                   tha,th3d,thten,thten1,thterm,                     &
                   qpten,qtten,qvten,qcten,qa,q3d,qten,              &
                   tkea,tke3d,tketen,pbl3d,                          &
                   pta,pt3d,ptten,thften,qvften,                     &
                   cfb,cfa,cfc, d1, d2,pdt,lgbth,lgbph,flag,         &
                   myi1p,myi2p,myj1p,myj2p,                          &
                   myi1pp,myi2pp,myj1pp,myj2pp,                      &
                   pgrid1,cgrid1,cgrid2,psenda,precva,psendb,precvb,psendc,precvc,psendd,precvd, &
                   rhsx,transx,rhsy,transy,                          &
                   npcomm,nplist,npi1,npi2,npj1,npj2,npk1,npk2,      &
                   reqs_u,reqs_v,reqs_w,reqs_s,reqs_p,reqs_p2,reqs_p3, &
                   reqs_x,reqs_y,reqs_z,reqs_tk,reqs_tk2,reqs_q,reqs_t, &
                   nw1,nw2,ne1,ne2,sw1,sw2,se1,se2,                  &
                   ww1,ww2,we1,we2,ws1,ws2,wn1,wn2,                  &
                   pw1,pw2,pe1,pe2,ps1,ps2,pn1,pn2,                  &
                   p2w1,p2w2,p2e1,p2e2,p2s1,p2s2,p2n1,p2n2,          &
                   p3w1,p3w2,p3e1,p3e2,p3s1,p3s2,p3n1,p3n2,          &
                   vw1,vw2,ve1,ve2,vs1,vs2,vn1,vn2,                  &
                   zw1,zw2,ze1,ze2,zs1,zs2,zn1,zn2,                  &
                   uw31,uw32,ue31,ue32,us31,us32,un31,un32,          &
                   vw31,vw32,ve31,ve32,vs31,vs32,vn31,vn32,          &
                   ww31,ww32,we31,we32,ws31,ws32,wn31,wn32,          &
                   sw31,sw32,se31,se32,ss31,ss32,sn31,sn32,          &
                   rw31,rw32,re31,re32,rs31,rs32,rn31,rn32,          &
                   qw31,qw32,qe31,qe32,qs31,qs32,qn31,qn32,          &
                   tkw1,tkw2,tke1,tke2,tks1,tks2,tkn1,tkn2,          &
                   tw1,tw2,te1,te2,ts1,ts2,tn1,tn2,                  &
                   tdiag,qdiag,udiag,vdiag,wdiag,kdiag,              &
                   out2d,out3d,                                      &
                   bndy,kbdy,hflxw,hflxe,hflxs,hflxn,                &
                   dowriteout,dorad,getdbz,dotdwrite,                &
                   dotbud,doqbud,doubud,dovbud,dowbud,               &
                   doazimwrite,dorestart)
        ! end_solve2


      IF(imoist.eq.1)THEN

        ! microphysics:

        ! cm1r20:  move calls to microphysics scheme into mp_driver.F

      if( ptype.eq.0 )then
        ! cm1r20.3:
        ! qv only:
        call pdefq(    0.0,asq(1),ruh,rvh,rmh,rho,q3d(ib,jb,kb,1))
      else
        call   mp_driver(nstep,dt,qbudget,asq,bsq,xh,ruh,xf,ruf,yh,rvh,yf,rvf,  &
                         mh,rmh,c1,c2,zh,mf,rmf,zf,rain,hail,prate,pi0,th0,rho0,prs0,qv0,   &
                         rho,prs,dum1,dum2,dum3,dum4,dum5,dum6,dum7,dum8,  &
                         w3d,ppi,pp3d,ppten,sten,tha,th3d,thten,qa,q3d,qten,   &
                         p3a,p3o,effc,effi,effs,effr,effg,effis,    &
                         tdiag,qdiag,out2d,out3d,tmp2d,tmp1d,  &
                         dowriteout,dorad,dotdwrite,doazimwrite,dorestart,  &
                         getdbz,getvt,getsed,getqdiags,dotbud,doqbud)
      endif

      ENDIF


        ! solve3: finish primary dynamical solver
        ! (finish comms, update arrays, move surface)
        call     solve3(nstep,num_soil_layers,                       &
                   dt,dtlast,mtime,dbldt,                            &
                   xh,rxh,arh1,arh2,uh,ruh,xf,rxf,arf1,arf2,uf,ruf,  &
                   yh,vh,rvh,yf,vf,rvf,                              &
                   rds,sigma,rdsf,sigmaf,dumk5,                      &
                   zh,mh,rmh,c1,c2,zf,mf,rmf,wprof,ug,vg,            &
                   pi0,rho0,prs0,thv0,th0,rth0,qv0,qc0,              &
                   qi0,rr0,rf0,rrf0,                                 &
                   zs,gz,rgz,gzu,rgzu,gzv,rgzv,dzdx,dzdy,            &
                   rain,raincv,hail,sws,svs,sps,srs,sgs,sus,shs,     &
                   dum1,dum2,dum3,dum4,dum5,dum6,dum7,dum8,          &
                   divx,rho,rr,rf,prs,                               &
                   u0,rru,ua,u3d,uten,uten1,                         &
                   v0,rrv,va,v3d,vten,vten1,                         &
                   rrw,wa,w3d,wten,wten1,                            &
                   ppi,pp3d,ppx,phi1,                                &
                   tha,th3d,qpten,qtten,qvten,qcten,qa,q3d,tkea,tke3d, &
                   pta,pt3d,pdata,haildata,u2pt,v2pt,                &
                   tsk,znt,zntmp,ust,f2d,hfx,qfx,qsfc,tslb,tmn,      &
                   reqs_u,reqs_v,reqs_w,reqs_s,reqs_p,reqs_p2,reqs_p3, &
                   reqs_x,reqs_y,reqs_z,reqs_tk,reqs_q,reqs_t,       &
                   nw1,nw2,ne1,ne2,sw1,sw2,se1,se2,                  &
                   n3w1,n3w2,n3e1,n3e2,s3w1,s3w2,s3e1,s3e2,          &
                   ww1,ww2,we1,we2,ws1,ws2,wn1,wn2,                  &
                   pw1,pw2,pe1,pe2,ps1,ps2,pn1,pn2,                  &
                   p2w1,p2w2,p2e1,p2e2,p2s1,p2s2,p2n1,p2n2,          &
                   p3w1,p3w2,p3e1,p3e2,p3s1,p3s2,p3n1,p3n2,          &
                   vw1,vw2,ve1,ve2,vs1,vs2,vn1,vn2,                  &
                   zw1,zw2,ze1,ze2,zs1,zs2,zn1,zn2,                  &
                   uw31,uw32,ue31,ue32,us31,us32,un31,un32,          &
                   vw31,vw32,ve31,ve32,vs31,vs32,vn31,vn32,          &
                   ww31,ww32,we31,we32,ws31,ws32,wn31,wn32,          &
                   sw31,sw32,se31,se32,ss31,ss32,sn31,sn32,          &
                   rw31,rw32,re31,re32,rs31,rs32,rn31,rn32,          &
                   qw31,qw32,qe31,qe32,qs31,qs32,qn31,qn32,          &
                   tkw1,tkw2,tke1,tke2,tks1,tks2,tkn1,tkn2,          &
                   kw1,kw2,ke1,ke2,ks1,ks2,kn1,kn2,                  &
                   tw1,tw2,te1,te2,ts1,ts2,tn1,tn2,                  &
                   udiag,vdiag,wdiag,                                &
                   lsnudge_u,lsnudge_v,lsnudge_th,lsnudge_qv,        &
                   dowriteout,dorad,getdbz,dotdwrite,                &
                   dotbud,doqbud,doubud,dovbud,dowbud,donudge,       &
                   doazimwrite,dorestart)
        ! end_solve3


        !--------------------------------
        !  Pressure decomposition:
        !  (a diagnostic)
        !  (NOT COMMON)

        ifpdcomp:  &
        IF( pdcomp )THEN
        IF( dowriteout .or. doazimwrite )THEN

          if(myid.eq.0) print *,'  Getting pressure diagnostics ... '

          call   pidcomp(dt,xh,rxh,arh1,arh2,uh,xf,rxf,arf1,arf2,uf,vh,vf,          &
                         gz,rgz,gzu,gzv,mh,rmh,mf,rmf,rds,rdsf,c1,c2,f2d,wprof,     &
                         pi0,th0,rth0,thv0,qv0,qc0,qi0,rho0,rr0,rf0,rrf0,u0,v0,     &
                         dum1,dum2,dum3,dum4,dum5,dum6,dum7,dum8,divx,              &
                         u3d,rru,uten,uten1,v3d,rrv,vten,vten1,w3d,rrw,wten,wten1,  &
                         rho,pp3d,th3d,q3d,udiag,vdiag,wdiag,pdiag,                 &
                         cfb,cfa,cfc, d1, d2,pdt,lgbth,lgbph,                       &
                         bndy,kbdy,hflxw,hflxe,hflxs,hflxn,                         &
                         myi1p,myi2p,myj1p,myj2p,                                   &
                         myi1pp,myi2pp,myj1pp,myj2pp,                               &
                         pgrid1,cgrid1,cgrid2,psenda,precva,psendb,precvb,psendc,precvc,psendd,precvd, &
                         rhsx,transx,rhsy,transy,npcomm,nplist,npi1,npi2,npj1,npj2,npk1,npk2)

          if(myid.eq.0) print *,'  ... finished pressure diagnostics '

          if(timestats.ge.1) time_poiss=time_poiss+mytime()

        ENDIF
        ENDIF  ifpdcomp


        !--------------------------------
        ! Step model time forward one dt
        !   (NOTE:  do not change mtime anywhere else!)

        mtime = mtime + dbldt

        domainlocx = domainlocx + dbldt*dble(umove)
        domainlocy = domainlocy + dbldt*dble(vmove)

        !--------------------------------

        if( stopit )then
          dostat      =  .true.
          dowriteout  =  .true.
          if( dodomaindiag ) dotdwrite = .true.
          if( doazimavg ) doazimwrite = .true.
          if( dohifrq ) dohifrqwrite = .true.
        endif

        !--------------------------------
        ! Time-averaging:
        ! note: ntim,rtim are set in param.F

        tavgcode:  &
        IF( dotimeavg )THEN

          ntavg(ntim) = ntavg(ntim) + 1

          rtavg(ntim) = rtavg(ntim) + dbldt

          IF( ntim.eq.5 )THEN
            ntot = ntavg(1) + ntavg(2) + ntavg(3) + ntavg(4) + ntavg(5)
            rtot = rtavg(1) + rtavg(2) + rtavg(3) + rtavg(4) + rtavg(5)
          ELSE
            ntot = 0
            qtem = 0.0
            do n=1,ntim
              ntot = ntot+ntavg(n)
              qtem = qtem+rtavg(n)
            enddo
            rtot = qtem
          ENDIF
          rdenom = 1.0d0/max( 1.0d-10 , rtot )

          printit = .true.
          if( printit .and. myid.eq.0 )then
            print *,'  getting time averages: '
!!!            do n=1,ntim
!!!              print *,'      nt,rt = ',n,ntavg(n),rtavg(n)
!!!            enddo
            print *,'    ntot,rtot = ',ntot,rtot
          endif

          k2 = min ( keta , ke2pt )

          ! 3d vars:
          do n = 1 , ntavr
            if(     n.eq.utav )then
              call tavg_plus1(rdenom,tavg(ibta,jbta,1,1,n),timavg(ibta,jbta,1,n),u3d(ib,jb,1),1,0,0,1,keta,dbldt)
            elseif( n.eq.vtav )then
              call tavg_plus1(rdenom,tavg(ibta,jbta,1,1,n),timavg(ibta,jbta,1,n),v3d(ib,jb,1),0,1,0,1,keta,dbldt)
            elseif( n.eq.wtav )then
              call tavg_plus1(rdenom,tavg(ibta,jbta,1,1,n),timavg(ibta,jbta,1,n),w3d(ib,jb,1),0,0,1,1,keta,dbldt)
            elseif( n.eq.uutav )then
              call tavg_plus1(rdenom,tavg(ibta,jbta,1,1,n),timavg(ibta,jbta,1,n),u2pt(ib2pt,jb2pt,1),0,0,0,1,k2,dbldt)
            elseif( n.eq.vvtav )then
              call tavg_plus1(rdenom,tavg(ibta,jbta,1,1,n),timavg(ibta,jbta,1,n),v2pt(ib2pt,jb2pt,1),0,0,0,1,k2,dbldt)
            elseif( n.eq.ttav )then
              call tavg_plus1(rdenom,tavg(ibta,jbta,1,1,n),timavg(ibta,jbta,1,n),th3d(ib,jb,1),0,0,0,1,keta,dbldt)
            elseif( n.eq.etav )then
              call tavg_plus1(rdenom,tavg(ibta,jbta,1,1,n),timavg(ibta,jbta,1,n),tke3d(ib,jb,1),0,0,1,1,keta,dbldt)
            elseif( n.eq.qtav )then
              call tavg_plus1(rdenom,tavg(ibta,jbta,1,1,n),timavg(ibta,jbta,1,n),q3d(ib,jb,1,nqv),0,0,0,1,keta,dbldt)
            else
              print *,'  n = ',n
              stop 32497
            endif
          enddo

          ! 2d vars:
          if( dot2p )then
          do n = 1 , nsfctavr
            if(     n.eq.  1  )then
              call tavg_plus1(rdenom,sfctavg(ibta,jbta,1,n),sfctimavg(ibta,jbta,n),ust(ib,jb),0,0,0,1,1,dbldt)
            elseif( n.eq.  2  )then
              call tavg_plus1(rdenom,sfctavg(ibta,jbta,1,n),sfctimavg(ibta,jbta,n),znt(ib,jb),0,0,0,1,1,dbldt)
            elseif( n.eq.  3  )then
              call tavg_plus1(rdenom,sfctavg(ibta,jbta,1,n),sfctimavg(ibta,jbta,n),mol(ib,jb),0,0,0,1,1,dbldt)
            else
              stop 32498
            endif
          enddo
          endif


          if( rtavg(ntim).ge.(rtim-dtfrac*dbldt) )then

            if(myid.eq.0) print *
            if(myid.eq.0) print *,'  Incrementing time-avg arrays: '
            if(myid.eq.0) print *
            if(myid.eq.0) print *,'    rtim,rtavg = ',rtim,rtavg(ntim)
            if(myid.eq.0) print *

            do n = 1 , ntim-1
              ntavg(n) = ntavg(n+1)
              rtavg(n) = rtavg(n+1)
            enddo
            ntavg(ntim) = 0
            rtavg(ntim) = 0.0

            ! 3d vars:
            do n = 1 , ntavr
              call tavg_shift(tavg(ibta,jbta,1,1,n),1,keta)
            enddo

            ! 2d vars:
            if( dot2p )then
            do n = 1 , nsfctavr
              call tavg_shift(sfctavg(ibta,jbta,1,n),1,1)
            enddo
            endif

          endif

          if(timestats.ge.1) time_tavg=time_tavg+mytime()

        ENDIF  tavgcode

        !--------------------------------

      endif  ifsolve


      !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      !cc   Identify center of TC  ccccccccccccccccccccccccccccccccccccccccccc
      !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      IF( axisymm.eq.0 .and. ( doazimavg .or. do_adapt_move ) )THEN
        IF( do_adapt_move .and. mtime.ge.adaptmovetim )THEN
        IF( dowriteout .or. dohifrqwrite .or. doazimwrite )THEN
        if( nstep.ge.1 )then

          if( .not. restarted )then
            if(myid.eq.0) print *,'  finding TC center, mtime = ',mtime
            call getpsmth(xh,yh,sigma,pp3d,dum1,xfref,yfref,icenter,jcenter,xcenter,ycenter,psmth)
            call getpcent(xh,yh,xfref,yfref,dum1,icenter,jcenter,xcenter,ycenter,psmth)
          endif

          if(myid.eq.0) print *,'  calling getrangle: '
          if(myid.eq.0) print *,'    icenter,jcenter = ',icenter,jcenter
          if(myid.eq.0) print *,'    xcenter,ycenter = ',xcenter,ycenter
          call   getrangle(xh,yh,xcenter,ycenter,cir,crr,cangle,icrs,crref,xhref,xfref)

          if(timestats.ge.1) time_misc=time_misc+mytime()
        endif
        ENDIF
        ENDIF
      ENDIF


      !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      !cc   Adaptive domain-moving option for TCs  ccccccccccccccccccccccccccc
      !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        ! Experimental (seems to work well, though)
        !
        ! Uses a smoothed field of pressure at the lowest model level to determine
        ! the position of a TC.

      IF( do_adapt_move )THEN

        IF( mtime.ge.(adaptmovetim-dtfrac*dbldt) )then

          if(myid.eq.0) print *,'  calling update_adapt_move '
          call   update_adapt_move(nstep,mtime,xh,yh,xfref,yfref,ug,vg,psmth,     &
                                   dum1,pp3d,u0,ua,u3d,v0,va,v3d,mvrec,nwritemv,  &
                                   icenter,jcenter,xcenter,ycenter,               &
                                   domainlocx,domainlocy)
          if(myid.eq.0) print *,'    umove,vmove = ',umove,vmove

          nwritemv = nwritemv+1
          adaptmovetim = adaptmovetim + adapt_move_frq

          if(timestats.ge.1) time_misc=time_misc+mytime()
        ENDIF

        if( restarted .and. myid.eq.0 ) print *,'  On restart: '
        if( restarted .and. myid.eq.0 ) print *,'    umove,vmove = ',umove,vmove

      ENDIF


      !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      !cc   radiation  ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      IF( radopt.ge.1 .and. dorad )THEN

        call     radiation_driver(mtime,radtim,dt,rbufsz,xh,yh,xf,yf,zf,rmh,c1,c2,     &
                   swten,lwten,swtenc,lwtenc,cldfra,o30,thraten,                       &
                   radsw,rnflx,radswnet,radlwin,dsr,olr,rad2d,                         &
                   effc,effi,effs,effr,effg,effis,                                     &
                   lwupt,lwuptc,lwdnt,lwdntc,lwupb,lwupbc,lwdnb,lwdnbc,                &
                   swupt,swuptc,swdnt,swdntc,swupb,swupbc,swdnb,swdnbc,                &
                   lwcf,swcf,coszr,                                                    &
                   xice,xsnow,xlat,xlong,coszen,swddir,swddni,swddif,hrang,            &
                   cldfra1_flag,dum1,dum2,dum3,dum4,dum5,dum6,dum7,dum8,dum9,sten ,thten, &
                   prs0,pi0,th0,prs,ppi,tha,rho,qa,                                    &
                   rth0s,prs0s,rho0s,tsk,albd,glw,gsw,emiss,xland,nstep,pbl3d)

        doit = .true.
        do while( doit )
          radtim = radtim+dtrad
          if( radtim.gt.mtime )then
            doit = .false.
          endif
        enddo
        if(timestats.ge.1) time_rad=time_rad+mytime()

      ENDIF


      !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      !cc   Prepare turbulence vars for next time step   cccccccccccccccccc
      !cc     (new since cm1r17)                         cccccccccccccccccc
      !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    IF( .not. noutcheck )THEN

      rtime=sngl(mtime)

      getsfc = .true.
      if( restarted ) getsfc = .false.
      if( interp_on_restart ) getsfc = .true.
      if( init_sfcmodel.eq.1 ) getsfc = .true.

      getpbl = .true.
!!!      if( restarted .and. ipbl.ge.1 ) getpbl = .false.

      update_sfc = .true.
      if( startup .or. restarted ) update_sfc = .false.

      if( restarted )then
        if( interp_on_restart )then
          dttmp = dt
        else
          dttmp = dtlast
        endif
        if(myid.eq.0) print *,'  dtlast,dt = ',dtlast,dt
        if(myid.eq.0) print *,'  dttmp for turb = ',dttmp
      else
        dttmp = dt
      endif


      call sfc_and_turb(getsfc,getpbl,nstep,dttmp,dosfcflx,cloudvar,qbudget,   &
                   avgsfcu,avgsfcv,avgsfcs,avgsfcsu,avgsfcsv,avgsfct,avgsfcq,avgsfcp, &
                   xh,rxh,arh1,arh2,uh,ruh,xf,rxf,arf1,arf2,uf,ruf,  &
                   yh,vh,rvh,yf,vf,rvf,                              &
                   rds,sigma,rdsf,sigmaf,zh,mh,rmh,c1,c2,zf,mf,rmf,  &
                   pi0s,rth0s,pi0,rho0,prs0,thv0,th0,qv0,rr0,rf0,rrf0, &
                   zs,gz,rgz,gzu,rgzu,gzv,rgzv,gx,gxu,gy,gyv,        &
                   tsk,thflux,qvflux,cd,ch,cq,z1,u1,v1,s1,th1,qv1,tlh,f2d,ustt,ut,vt,st,cm0,   &
                   dum1,dum2,dum3,dum4,dum5,dum6,dum7,dum8,dum9,     &
                   divx,rho,rr,rf,prs,                               &
                   t11,t12,t13,t22,t23,t33,                          &
                   m11,m12,m13,m22,m23,m33,                          &
                   u0,rru,u3d,uten,uten1,                            &
                   v0,rrv,v3d,vten,vten1,                            &
                   rrw,w3d,wten ,wten1,                              &
                   pp3d,ppten,sten,sadv,                             &
                   th3d,thten,thten1,thterm,q3d,                     &
                   kmh,kmv,khh,khv,cme,csm,ce1,ce2,kmv2,khv2,tkea,tke3d,       &
                   nm,defv,defh,lenscl,dissten,epst,radsw,radswnet,radlwin, &
                   thpten,qvpten,qcpten,qipten,upten,vpten,qnipten,qncpten,qrpten,qspten,qgpten, &
                   xkzh,xkzq,xkzm,pbl3d,pbl2d,pbl2dint,              &
                   lwten,swten,thraten,tdiag,raincv,pratecv,                 &
                   thcuten,qvcuten,qccuten,qicuten,ucuten,vcuten,thften,qvften, &
                   lu_index,kpbl2d,psfc,u10,v10,s10,hfx,qfx,         &
                   xland,znt,rznt,ust,stau,tst,qst,z0t,z0q,          &
                   hpbl,wspd,phim,phih,psim,psih,psiq,gz1oz0,br,brcr, &
                   CHS,CHS2,CQS2,CPMM,ZOL,MAVAIL,                    &
                   MOL,RMOL,REGIME,LH,FLHC,FLQC,QGH,                 &
                   CK,CKA,CDA,USTM,QSFC,T2,Q2,TH2,EMISS,THC,ALBD,    &
                   gsw,glw,chklowq,capg,snowc,snowh,qcg,dsxy,wstar,delta,prkpp,fm,fh,  &
                   charn,msang,scurx,scury,zkmax,cd_out,ch_out,wscale,wscaleu, &
                   zntmp,mznt,swspd,smois,taux,tauy,hpbl2d,evap2d,heat2d,  &
                   tmp2d,tmp1d,  &
                   num_soil_layers,slab_zs,slab_dzs,tslb,tmn,        &
                   tml,t0ml,hml,h0ml,huml,hvml,tmoml,                &
                   cavg,uavg,vavg,savg,gamk,gamwall,kmw,ufw,vfw,u1b,v1b,l2p,s2p,s2b,t2pm1,t2pm2,t2pm3,t2pm4,u2pt,v2pt,kmwk,ufwk,vfwk, &
                   bndy,kbdy,timavg,sfctimavg,                       &
                   reqs_u,reqs_v,reqs_w,reqs_s,reqs_p,reqs_x,reqs_y, &
                   nw1,nw2,ne1,ne2,sw1,sw2,se1,se2,                  &
                   pw1,pw2,pe1,pe2,ps1,ps2,pn1,pn2,                  &
                   vw1,vw2,ve1,ve2,vs1,vs2,vn1,vn2,                  &
                   uw31,uw32,ue31,ue32,us31,us32,un31,un32,          &
                   sw31,sw32,se31,se32,ss31,ss32,sn31,sn32,          &
                   kw1,kw2,ke1,ke2,ks1,ks2,kn1,kn2,reqk,             &
                   n3w1,n3w2,n3e1,n3e2,s3w1,s3w2,s3e1,s3e2,          &
                   flag  ,out2d,out3d,rtime,update_sfc,dotbud,dotdwrite,restarted,restart_time,spinup_t2p,  &
                   dowriteout,doazimwrite,dostat)
        ! end_sfc_and_turb


      IF( (idiff.ge.1).and.(difforder.eq.2) )THEN
        !  get stress terms for explicit diffusion scheme:
        call diff2def(uh,arh1,arh2,uf,arf1,arf2,vh,vf,mh,c1,c2,mf,ust,zntmp,u1,v1,s1,  &
                      divx,rho,rr,rf,t11,t12,t13,t22,t23,t33,u3d,v3d,w3d,dissten)
      ENDIF

#ifdef MPI
      if( cm1setup.ge.1 .and. adapt_dt.eq.1 )then
        call mpi_wait(reqk,mpi_status_ignore,ierr)
        if(timestats.ge.1) time_cflq=time_cflq+mytime()
      endif
#endif

    ENDIF

      !cccccccccccccccccccccccccccccccccccccccccccccccccccccccc!
      !cccccccccccccccccccccccccccccccccccccccccccccccccccccccc!
      !cccccccccccccccccccccccccccccccccccccccccccccccccccccccc!


        !--------------------------------
        !  Adaptive timestepping:
        !   (assumes cflmax and ksmax have already been calculated)

        IF( adapt_dt.eq.1 .and. ( dosolve .or. interp_on_restart ) )THEN

          call   getnewdt(ndt,dt,dtlast,adt,acfl,dbldt,                                 &
                          mtime,stattim,taptim,rsttim,prcltim,hailtrajtim,     &
                          diagtim,azimavgtim,                                   &
                          dorestart,dowriteout,dostat,doprclout,dohailout,      &
                          dotdwrite,doazimwrite,                                &
                          hifrqtim,dohifrqwrite,.false.)

        ENDIF


      !cccccccccccccccccccccccccccccccccccccccccccccccccccccccc!
      !cccccccccccccccccccccccccccccccccccccccccccccccccccccccc!
      !cccccccccccccccccccccccccccccccccccccccccccccccccccccccc!


      rtime=sngl(mtime)
      reset = .false.


      if( dostat .and. (.not.noutcheck) )then
        if( startup )then
          if(dowr) write(outfile,*)
          if(dowr) write(outfile,*) '  initial conditions:'
          if(dowr) write(outfile,*)
        endif
        IF(axisymm.eq.0)THEN
          ! for Cartesian grid:
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
          do k=1,nk
          do j=1,nj
          do i=1,ni
            ppten(i,j,k)=rho(i,j,k)
          enddo
          enddo
          enddo
        ELSE
          ! for axisymmetric grid:
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
          do k=1,nk
          do j=1,nj
          do i=1,ni
            ppten(i,j,k) = rho(i,j,k)*pi*(xf(i+1)**2-xf(i)**2)/(dx*dy)
          enddo
          enddo
          enddo
        ENDIF
        call     statpack(mtime,nrec,ndt,dt,dtlast,rtime,adt,acfl,cloudvar,   &
                          qname,budname,qbudget,asq,bsq,                      &
                          name_stat,desc_stat,unit_stat,                      &
                          xh,rxh,uh,ruh,xf,uf,yh,vh,rvh,vf,zh,mh,rmh,zf,mf,   &
                          zs,rgzu,rgzv,rds,sigma,rdsf,sigmaf,                 &
                          rstat,pi0,rho0,thv0,th0,qv0,u0,v0,                  &
                          dum1,dum2,dum3,dum4,dum5,ppten,prs,                 &
                          ua,va,wa,ppi,tha,qa,qten,kmh,kmv,khh,khv,tkea,      &
                          xkzh,xkzq,xkzm,pbl3d,                               &
                          pta,u10,v10,hpbl,prate,reset,nstatout,restarted)
        nrec = nrec + 1
        nstatout = nstatout + 1
        if( statfrq.gt.0.0 .and. ( .not. restarted ) )then
          doit = .true.
          do while( doit )
            stattim = stattim+statfrq
            if( stattim.gt.mtime )then
              doit = .false.
            endif
          enddo
        endif
        if(myid.eq.0) print *,'  next stattim = ',stattim
        if(timestats.ge.1) time_stat=time_stat+mytime()
      endif


      !cccccccccccccccccccccccccccccccccccccccccccccccccccccccc!
      !cccccccccccccccccccccccccccccccccccccccccccccccccccccccc!
      !cccccccccccccccccccccccccccccccccccccccccccccccccccccccc!


        if(myid.eq.0)then
          if(timeformat.eq.1)then
            write(6,110) nstep,rtime,' sec '
          elseif(timeformat.eq.2)then
            write(6,110) nstep,rtime/60.0,' min '
          elseif(timeformat.eq.3)then
            write(6,110) nstep,rtime/3600.0,' hour'
          elseif(timeformat.eq.4)then
            write(6,110) nstep,rtime/86400.0,' day '
          else
            write(6,110) nstep,rtime,' sec'
          endif
110       format(2x,i12,4x,f18.6,a5)
        endif


      !cccccccccccccccccccccccccccccccccccccccccccccccccccccccc!
      !cccccccccccccccccccccccccccccccccccccccccccccccccccccccc!
      !cccccccccccccccccccccccccccccccccccccccccccccccccccccccc!

      if( interp_on_restart .and. ( .not. noutcheck ) )then
        ! always write initial conditions for interp_on_restart
        nwrite = 1
        taptim = mtime
        dowriteout = .true.
      endif


      calling_writeout:  &
      if( dowriteout .or. noutcheck )then
      otype:  &
      IF(output_format.eq.1.or.output_format.eq.2)THEN
        nn = 1
        if(terrain_flag .and. output_interp.eq.1) nn = 2
        DO n=1,nn
          if(n.eq.1)then
            fnum = 51
            frec = srec
          else
            fnum = 71
            frec = sirec
          endif
          IF( stopit )THEN
            ! diag code does not account for stopit, so just set arrays to zero
            ! to avoid confusion
            tdiag = 0.0
            qdiag = 0.0
            udiag = 0.0
            vdiag = 0.0
            wdiag = 0.0
            kdiag = 0.0
            pdiag = 0.0
          ENDIF
          call          writeout(frec,urec,vrec,wrec,xrec,yrec,zrec,mrec,rtime,dt,fnum,nwrite,qname,qunit,    &
                        nstep,mtime,adt,ndt,domainlocx,domainlocy,                             &
                        name_output,desc_output,unit_output,grid_output,levs_output,cmpr_output, &
                        xh,xf,uf,yh,yf,vf,xfref,yfref,                                         &
                        rds,sigma,rdsf,sigmaf,zh,zf,mf,gx,gy,wprof,                            &
                        pi0,prs0,rho0,rr0,rf0,rrf0,th0,qv0,u0,v0,thv0,rth0,qc0,qi0,            &
                        zs,rgzu,rgzv,rain,raincv,hail,sws,svs,sps,srs,sgs,sus,shs,thflux,qvflux,psfc,      &
                        rxh,arh1,arh2,uh,ruh,rxf,arf1,arf2,vh,rvh,mh,rmh,rmf,rr,rf,            &
                        gz,rgz,gzu,gzv,gxu,gyv,dzdx,dzdy,c1,c2,                                &
                        cd,ch,cq,tlh,f2d,psmth,prate,ustt,cm0,eimaskw,eimaske,eimasks,eimaskn, &
                        dum1,dum2,dum3,dum4,dum5,dum6,dum7,dum8,dum9,dum10,                    &
                        t11,t12,t13,t22,t23,t33,rho,prs,divx,                                  &
                        rru,u3d,uten,uten1,rrv,v3d,vten,vten1,rrw,w3d,wten,pp3d,th3d,phi2,     &
                        sadv,thten,nm,defv,defh,lenscl,dissten,epst,                           &
                        thpten,qvpten,qcpten,qipten,upten,vpten,qnipten,qncpten,xkzh,xkzq,xkzm, &
                        lu_index,xland,mavail,tsk,tmn,tml,hml,huml,hvml,hfx,qfx,gsw,glw,tslb,  &
                        q3d,p3a,p3o,kmh,kmv,khh,khv,kmv2,khv2,cme,tke3d,swten,lwten,cldfra,o30,          &
                        radsw,rnflx,radswnet,radlwin,dsr,olr,pta,                              &
                        effc,effi,effs,effr,effg,effis,                                        &
                        lwupt,lwdnt,lwupb,lwdnb,                                               &
                        swupt,swdnt,swupb,swdnb,lwcf,swcf,coszen,                              &
                        num_soil_layers,u10,v10,s10,t2,q2,znt,ust,stau,tst,qst,z0t,z0q,u1,v1,s1,     &
                        hpbl,zol,mol,rmol,br,brcr,wscale,wscaleu,phim,phih,psim,psih,psiq,wspd,qsfc,wstar,delta,prkpp,fm,fh, &
                        mznt,taux,tauy,pbl3d,pbl2d,                                            &
                        dat1,dat2,dat3,reqt,thten1(ib,jb,1),dumk1,dumk2,                       &
                        tdiag,qdiag,udiag,vdiag,wdiag,pdiag,out2d,out3d,cir,crr,cangle,        &
                        bndy,kbdy,hflxw,hflxe,hflxs,hflxn,recy_cap_we,recy_cap_sn,recy_inj_we,recy_inj_sn,recy_wsp_w,recy_wsp_e,recy_wsp_s,recy_wsp_n,recy_tau_w,recy_tau_e,recy_tau_s,recy_tau_n,timavg,sfctimavg,kmwk, &
                        thcuten,qvcuten,qccuten,qicuten,ucuten,vcuten,                         &
                        nw1,nw2,ne1,ne2,sw1,sw2,se1,se2,myi1p,myi2p,myj1p,myj2p)
                        ! end_writeout
          if(n.eq.1)then
            srec = frec
          else
            sirec = frec
          endif
        ENDDO
      ELSE  otype
        print *,'  09832 '
        call stopcm1
      ENDIF  otype
      if( writeit )then
        nwrite=nwrite+1
        if(tapfrq.gt.0.0)then
          doit = .true.
          do while( doit )
            taptim = taptim+tapfrq
            if( taptim.gt.mtime )then
              doit = .false.
            endif
          enddo
        endif
        if(myid.eq.0) print *,'  next taptim = ',taptim
      endif
        if(timestats.ge.1) time_write=time_write+mytime()
      endif  calling_writeout


      !cccccccccccccccccccccccccccccccccccccccccccccccccccccccc!
      !cccccccccccccccccccccccccccccccccccccccccccccccccccccccc!
      !cccccccccccccccccccccccccccccccccccccccccccccccccccccccc!


      if(iprcl.eq.1)then
      IF( doprclout .and. (.not.noutcheck) )THEN
        call     parcel_interp(dt,mtime,xh,uh,ruh,xf,uf,yh,vh,rvh,yf,vf, &
                               zh,mh,rmh,zf,mf,zntmp,ust,c1,c2,        &
                               zs,sigma,sigmaf,rds,gz,                 &
                               pi0,th0,thv0,qv0,qc0,qi0,rth0,          &
                               dum1,dum2,dum3,dum4,dum5,dum6,prs,rho,  &
                               dum7,dum8,wten,wten1,                   &
                               u3d,v3d,w3d,pp3d,thten,thten1,th3d,q3d, &
                               kmh,kmv,khh,khv,tke3d,pt3d,pdata,       &
                               tdiag,qdiag,                            &
                               pw1,pw2,pe1,pe2,ps1,ps2,pn1,pn2,        &
                               nw1,nw2,ne1,ne2,sw1,sw2,se1,se2,reqs_p, &
                               tkw1,tkw2,tke1,tke2,tks1,tks2,tkn1,tkn2)
        call     parcel_write(prec,rtime,qname,name_prcl,desc_prcl,unit_prcl,pdata,ploc)
        prec = prec+1
        if(prclfrq.gt.0.0)then
          doit = .true.
          do while( doit )
            prcltim = mtime+prclfrq
            if( prcltim.gt.mtime )then
              doit = .false.
            endif
          enddo
        endif
        if(myid.eq.0) print *,'  next prcltim = ',prcltim
        if(timestats.ge.1) time_parcels=time_parcels+mytime()
      ELSE
        if( startup )then
          if(myid.eq.0) print *
          if(myid.eq.0) print *,'  NOTE:  skipping parcel_write '
          if(myid.eq.0) print *
        endif
      ENDIF
      endif

      !cccccccccccccccccccccccccccccccccccccccccccccccccccccccc!
      !cccccccccccccccccccccccccccccccccccccccccccccccccccccccc!
      !cccccccccccccccccccccccccccccccccccccccccccccccccccccccc!


      if(ihailtraj.eq.1)then
      IF( dohailout .and. (.not.noutcheck) )THEN
        call     hailcast_write(hrec,rtime,qname,name_hail,desc_hail,unit_hail,haildata,hailloc)
        hrec = hrec+1
        if(hailfrq.gt.0.0)then
          doit = .true.
          do while( doit )
            hailtrajtim = hailtrajtim+hailfrq
            if( hailtrajtim.gt.mtime )then
              doit = .false.
            endif
          enddo
        endif
        if(timestats.ge.1) time_hailtrajs=time_hailtrajs+mytime()
      ELSE
        if( startup )then
          if(myid.eq.0) print *
          if(myid.eq.0) print *,'  NOTE:  skipping hailcast_write '
          if(myid.eq.0) print *
        endif
      ENDIF
      endif


      !cccccccccccccccccccccccccccccccccccccccccccccccccccccccc!
      !cccccccccccccccccccccccccccccccccccccccccccccccccccccccc!
      !cccccccccccccccccccccccccccccccccccccccccccccccccccccccc!


!      if( restarted )then
!        nwritet = nwritet-1
!        diagtim = diagtim-diagfrq
!        dotdwrite = .true.
!      endif


    IF( dodomaindiag )THEN
      if( dotdwrite .or. noutcheck )then

        call   domaindiag(nstep,mtime,nwritet,trecs,trecw,qname,qunit,dt,dosfcflx, &
                   xh,rxh,arh1,arh2,uh,ruh,xf,rxf,arf1,arf2,uf,ruf,            &
                   yh,vh,rvh,yf,vf,rvf,                                        &
                   xfref,yfref,rds,sigma,rdsf,sigmaf,                          &
                   wprof,ufrc,vfrc,thfrc,qvfrc,ug,vg,dvdr,                     &
                   uavg,vavg,thavg,pavg,ulspg,vlspg,qavg,cavg,cloudvar,        &
                   tauh,taus,zh,mh,rmh,c1,c2,tauf,zf,mf,rmf,                   &
                   rho0s,pi0s,prs0s,rth0s,                                     &
                   pi0,rho0,prs0,thv0,th0,rth0,qv0,qc0,u0,v0,                  &
                   qi0,rr0,rf0,rrf0,                                           &
                   zs,gz,rgz,gzu,rgzu,gzv,rgzv,dzdx,dzdy,gx,gxu,gy,gyv,        &
                   tsk,znt,ust,stau,tst,qst,z0t,z0q,thflux,qvflux,             &
                   cd,ch,cq,u1,v1,s1,xland,psfc,tlh,prate,ustt,ut,vt,st,       &
                   dum1,dum2,dum3,dum4,dum5,dum6,dum7,dum8,dum9,               &
                   divx,rho,rr,rf,prs,t11,t12,t13,t22,t23,t33,                 &
                   m11,m12,m13,m22,m23,m33,                                    &
                   rru,u3d,uten,uten1,                                         &
                   rrv,v3d,vten,vten1,                                         &
                   rrw,w3d,wten,wten1,                                         &
                   pp3d,sadv,ppten ,sten,th3d,thten,thten1,thterm,phi2,        &
                   q3d,qten,kmh,kmv,khh,khv,cme,csm,ce1,ce2,tkea,tke3d,tketen, &
                   nm,defv,defh,lenscl,dissten,epst,pt3d,                      &
                   rain,raincv,hfx,qfx,u10,v10,s10,t2,q2,th2,mznt,             &
                   hpbl,wspd,zol,mol,rmol,br,brcr,phim,phih,psim,psih,         &
                   thpten,qvpten,qcpten,qipten,upten,vpten,qnipten,qncpten,xkzh,xkzq,xkzm, &
                   qsfc,o30,zir,swten,lwten,swtenc,lwtenc,cldfra,              &
                   effc,effi,effs,effr,effg,effis,                             &
                   lwupt,lwuptc,lwdnt,lwdntc,lwupb,lwupbc,lwdnb,lwdnbc,        &
                   swupt,swuptc,swdnt,swdntc,swupb,swupbc,swdnb,swdnbc,        &
                   lwcf,swcf,pbl3d,pbl2d,                                      &
                   tdiag,qdiag,udiag,vdiag,wdiag,kdiag,out2d,out3d,            &
                   gamk,gamwall,kmw,ufw,vfw,u1b,v1b,l2p,s2p,s2b,t2pm1,t2pm2,t2pm3,t2pm4,kmwk,ufwk,vfwk, &
                   lsnudge_u,lsnudge_v,lsnudge_th,lsnudge_qv,                  &
                   domainlocx,domainlocy,                                      &
                   thcuten,qvcuten,qccuten,qicuten,ucuten,vcuten,              &
                   sw31,sw32,se31,se32,ss31,ss32,sn31,sn32,flag)
        ! end_domaindiag
      if( writeit )then
        nwritet = nwritet+1
        if(diagfrq.gt.0.0)then
          doit = .true.
          do while( doit )
            diagtim = diagtim+diagfrq
            if( diagtim.gt.mtime )then
              doit = .false.
            endif
          enddo
        endif
      endif
        if(myid.eq.0) print *,'  next diagtim = ',diagtim
        if(timestats.ge.1) time_diag = time_diag+mytime()

      endif
    ENDIF


      !cccccccccccccccccccccccccccccccccccccccccccccccccccccccc!
      !cccccccccccccccccccccccccccccccccccccccccccccccccccccccc!
      !cccccccccccccccccccccccccccccccccccccccccccccccccccccccc!


    IF( dohifrq )THEN
      if( dohifrqwrite .or. noutcheck )then

        call   writeout_hifrq(                                                 &
                   nstep,mtime,nwriteh,qname,dt,dosfcflx,                      &
                   icrs,icenter,jcenter,xcenter,ycenter,                       &
                   xh,rxh,arh1,arh2,uh,ruh,xf,rxf,arf1,arf2,uf,ruf,            &
                   yh,vh,rvh,yf,vf,rvf,                                        &
                   xfref,xhref,yfref,yhref,rds,sigma,rdsf,sigmaf,              &
                   tauh,taus,zh,mh,rmh,c1,c2,tauf,zf,mf,rmf,                   &
                   rho0s,pi0s,prs0s,rth0s,                                     &
                   pi0,rho0,prs0,thv0,th0,rth0,qv0,qc0,u0,v0,                  &
                   qi0,rr0,rf0,rrf0,                                           &
                   zs,gz,rgz,gzu,rgzu,gzv,rgzv,dzdx,dzdy,gx,gxu,gy,gyv,        &
                   tsk,znt,ust,stau,thflux,qvflux,cd,ch,cq,u1,v1,s1,xland,psfc,tlh,psmth, &
                   dum1,dum2,dum3,dum4,dum5,dum6,dum7,dum8,                    &
                   divx,rho,rr,rf,prs,t11,t12,t13,t22,t23,t33,                 &
                   rru,u3d,uten,uten1,                                         &
                   rrv,v3d,vten,vten1,                                         &
                   rrw,w3d,wten,wten1,                                         &
                   pp3d,ppten,sten,th3d,sadv,thten,thten1,                     &
                   q3d,qten,kmh,kmv,khh,khv,tkea,tke3d,tketen,lenscl,          &
                   nm,defv,defh,dissten,pbl3d,                                 &
                   thpten,qvpten,qcpten,qipten,upten,vpten,xkzh,xkzq,xkzm,     &
                   rain,raincv,u10,v10,s10,br,brcr,hpbl,t2,q2,hfx,qfx,prate,cm0,eimaskw,eimaske,eimasks,eimaskn, &
                   swten,lwten,cldfra,                                         &
                   lwupt,lwuptc,lwdnt,lwdntc,lwupb,lwupbc,lwdnb,lwdnbc,        &
                   swupt,swuptc,swdnt,swdntc,swupb,swupbc,swdnb,swdnbc,        &
                   lwcf,swcf,wprof,dumk1,dumk2,                                &
                   tdiag,qdiag,udiag,vdiag,wdiag,pdiag,out2d,out3d,getdbz,     &
                   cir,crr,cangle,recy_cap_we,recy_cap_sn,recy_inj_we,recy_inj_sn,recy_wsp_w,recy_wsp_e,recy_wsp_s,recy_wsp_n,recy_tau_w,recy_tau_e,recy_tau_s,recy_tau_n,timavg,tavg,rtavg,kmwk, &
                   myi1p,myi2p,myj1p,myj2p,kbdy,bndy,                          &
                   sw31,sw32,se31,se32,ss31,ss32,sn31,sn32,flag,               &
                   dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,reqs_s)
        ! end hifrq

      if( writeit )then
        nwriteh = nwriteh+1
        if(hifrqfrq.gt.0.0)then
          doit = .true.
          do while( doit )
            hifrqtim = hifrqtim+hifrqfrq
            if( hifrqtim.gt.mtime )then
              doit = .false.
            endif
          enddo
        endif
      endif
        if(myid.eq.0) print *,'  next hifrqtim = ',hifrqtim
        if(timestats.ge.1) time_hifrq = time_hifrq+mytime()

      endif
    ENDIF


      !cccccccccccccccccccccccccccccccccccccccccccccccccccccccc!
      !cccccccccccccccccccccccccccccccccccccccccccccccccccccccc!
      !cccccccccccccccccccccccccccccccccccccccccccccccccccccccc!


!      if( restarted )then
!        nwritea = nwritea-1
!        azimavgtim = azimavgtim-azimavgfrq
!        doazimwrite = .true.
!      endif


    IF( doazimavg )THEN
      if( doazimwrite .or. noutcheck )then

        call   azimavg(nstep,mtime,nwritea,arecs,arecw,qname,dt,dosfcflx,      &
                   icrs,icenter,jcenter,xcenter,ycenter,domainlocx,domainlocy, &
                   xh,rxh,arh1,arh2,uh,ruh,xf,rxf,arf1,arf2,uf,ruf,            &
                   yh,vh,rvh,yf,vf,rvf,                                        &
                   xfref,yfref,xhref,yhref,rds,sigma,rdsf,sigmaf,              &
                   tauh,taus,zh,mh,rmh,c1,c2,tauf,zf,mf,rmf,                   &
                   rho0s,pi0s,prs0s,rth0s,                                     &
                   pi0,rho0,prs0,thv0,th0,rth0,qv0,qc0,u0,v0,                  &
                   qi0,rr0,rf0,rrf0,                                           &
                   zs,gz,rgz,gzu,rgzu,gzv,rgzv,dzdx,dzdy,gx,gxu,gy,gyv,        &
                   tsk,znt,ust,thflux,qvflux,cd,ch,cq,u1,v1,s1,xland,psfc,tlh,cm0, &
                   dum1,dum2,dum3,dum4,dum5,dum6,dum7,dum8,                    &
                   divx,rho,rf,prs,t11,t12,t13,t22,t23,t33,                    &
                   rru,u3d,uten,uten1,                                         &
                   rrv,v3d,vten,vten1,                                         &
                   rrw,w3d,wten,wten1,                                         &
                   pp3d,ppten,sten,th3d,sadv,thten,thten1,                     &
                   q3d,qten,kmh,kmv,khh,khv,cme,csm,ce1,ce2,kmv2,khv2,tkea,tke3d,tketen, &
                   nm,defv,defh,lenscl,dissten,                                &
                   thpten,qvpten,qcpten,qipten,upten,vpten,xkzh,xkzq,xkzm,     &
                   rain,raincv,u10,v10,s10,br,brcr,hpbl,prate,t2,q2,           &
                   swten,lwten,cldfra,pbl3d,                                   &
                   tdiag,qdiag,udiag,vdiag,wdiag,pdiag,out2d,out3d,getdbz,     &
                   effc,effi,effs,effr,effg,effis,                             &
                   lwupt,lwdnt,lwupb,lwdnb,                                    &
                   swupt,swdnt,swupb,swdnb,lwcf,swcf,coszen,                   &
                   cir,crr,cangle,crref,kmwk,psmth,                            &
                   thcuten,qvcuten,qccuten,qicuten,ucuten,vcuten,              &
                   dat1(1,1),dat2(1,1),dat3(1,1,1),                            &
                   sw31,sw32,se31,se32,ss31,ss32,sn31,sn32,flag)
        ! end_azimavg
      if( writeit )then
        nwritea = nwritea+1
        if(azimavgfrq.gt.0.0)then
          doit = .true.
          do while( doit )
            azimavgtim = azimavgtim+azimavgfrq
            if( azimavgtim.gt.mtime )then
              doit = .false.
            endif
          enddo
        endif
      endif
        if(myid.eq.0) print *,'  next azimavgtim = ',azimavgtim
        if(timestats.ge.1) time_azimavg = time_azimavg+mytime()

      endif
    ENDIF


      !cccccccccccccccccccccccccccccccccccccccccccccccccccccccc!
      !cccccccccccccccccccccccccccccccccccccccccccccccccccccccc!
      !cccccccccccccccccccccccccccccccccccccccccccccccccccccccc!


      IF( adapt_dt.eq.1 .and. reset .and. (.not. interp_on_restart) )THEN
        adtlast  = adt/dble(max(1,ndt))
        ndt  = 0
        adt  = 0.0
        acfl = 0.0
      ENDIF

        if( convinit.eq.1 )then
          if( mtime.gt.convtime ) convinit = 0
        endif


      !cccccccccccccccccccccccccccccccccccccccccccccccccccccccc!
      !cccccccccccccccccccccccccccccccccccccccccccccccccccccccc!
      !cccccccccccccccccccccccccccccccccccccccccccccccccccccccc!


      if( dorestart )then
        if(rstfrq.gt.0.0)then
          doit = .true.
          do while( doit )
            rsttim = rsttim+rstfrq
            if( rsttim.gt.mtime )then
              doit = .false.
            endif
          enddo
        endif
#ifdef MPI
        ! barrier for timing purposes:
        if(timestats.ge.1) call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif
        call     restart_write(nstep,srec,sirec,urec,vrec,wrec,nrec,mrec,prec,      &
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
                               rho,prs,ua,uten,va,vten,wa,ppi,tha,qa,tkea,          &
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
        if(myid.eq.0) print *,'  next rsttim = ',rsttim
        IF( nrst.eq.0 )THEN
          ! 150820:  Write restart file if rstfrq=0 and stop
          !          (useful for ensemble DA)
          if(dowr) write(outfile,*)
          if(dowr) write(outfile,*) '  Detected rstfrq = 0 '
          if(dowr) write(outfile,*) '  ... writing restart file ... '
          if(dowr) write(outfile,*)
          if(dowr) write(outfile,*) '  ... stopping ... '
          if(dowr) write(outfile,*)
#ifdef MPI
          call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif
          stop 55556
        ENDIF
        nrst = nrst+1
        ! barrier for timing purposes:
        if(timestats.ge.1)then
#ifdef MPI
          call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif
          time_restart=time_restart+mytime()
        endif
      endif


      !cccccccccccccccccccccccccccccccccccccccccccccccccccccccc!
      !cccccccccccccccccccccccccccccccccccccccccccccccccccccccc!
      !cccccccccccccccccccccccccccccccccccccccccccccccccccccccc!

        if(timestats.eq.2)then
          steptime2=time_sound+time_poiss+time_buoyan+time_turb+            &
                    time_diffu+time_microphy+time_stat+time_cflq+           &
                    time_bc+time_misc+time_integ+time_rdamp+time_divx+      &
                    time_write+time_restart+time_ttend+time_cor+time_fall+  &
                    time_satadj+time_dbz+time_sfcphys+time_parcels+         &
                    time_hailtrajs+                                         &
                    time_rad+time_pbl+time_swath+time_pdef+time_prsrho+     &
                    time_diag+time_azimavg+time_hifrq+time_ercyl+time_tavg+time_cuprm+ &
#ifdef MPI
                    time_mpu1+time_mpv1+time_mpw1+time_mpp1+                &
                    time_mpu2+time_mpv2+time_mpw2+time_mpp2+                &
                    time_mps1+time_mps3+time_mpq1+time_mptk1+                         &
                    time_mps2+time_mps4+time_mpq2+time_mptk2+time_mpb+                &
#endif
                    time_advs+time_advu+time_advv+time_advw
          write(6,157) nstep,steptime2-steptime1
157       format('    timing for time step ',i12,':',f12.4,' s')
          steptime1 = steptime2
        endif

        !--------------------------------------------------------------------

        if(iconly.eq.1)then
          if( startup .or. interp_on_restart )then
            if(dowr) write(outfile,*)
            if(dowr) write(outfile,*) '  User has requested initial conditions only'
            if(dowr) write(outfile,*) '     (iconly = 1)'
            if(dowr) write(outfile,*) '  ... stopping ... '
            if(dowr) write(outfile,*)
#ifdef MPI
            call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif
            stop 55555
          endif
        endif

        !--------------------------------------------------------------------
        ! Time step complete.  Odds and ends:

      if(stopit)then
        if(myid.eq.0)then
          print *
          if( stopcondition.eq.1 )then
            print *,' Courant number has exceeded 1.5 '
          elseif( stopcondition.eq.2 )then
            print *,' Courant number is negative '
          elseif( stopcondition.eq.3 )then
            print *,' Courant number is NaN '
          else
            print *,' unknown value of stopcondition '
          endif
          print *
          print *,' Stopping model .... '
          print *
        endif
#ifdef MPI
        call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif
        call stopcm1
      endif


      if( noutcheck )then
        increment = .false.
        myid0a:  &
        if(myid.eq.0)then
          print *
          print *,' nconfig    = ',nconfig
          print *
          print *,' noutstat   = ',noutstat,noutstat0
          print *,' nout2ds    = ',nout2ds,nout2ds0
          print *,' nout3ds    = ',nout3ds,nout3ds0
          print *,' nout3du    = ',nout3du,nout3du0
          print *,' nout3dv    = ',nout3dv,nout3dv0
          print *,' nout3dw    = ',nout3dw,nout3dw0
          print *,' nout3dxv   = ',nout3dxv,nout3dxv0
          print *,' nout3dyv   = ',nout3dyv,nout3dyv0
          print *,' nout3dzv   = ',nout3dzv,nout3dzv0
          print *,' noutdiags  = ',noutdiags,noutdiags0
          print *,' noutdiagw  = ',noutdiagw,noutdiagw0
          print *,' noutazims  = ',noutazims,noutazims0
          print *,' noutazimw  = ',noutazimw,noutazimw0
          print *,' noutparcel = ',noutparcel,noutparcel0
          print *,' nouthifrq  = ',nouthifrq,nouthifrq0

          if( (noutstat.ne.noutstat0) .or. &
              (nout2ds.ne.nout2ds0) .or. &
              (nout3ds.ne.nout3ds0) .or. &
              (nout3du.ne.nout3du0) .or. &
              (nout3dv.ne.nout3dv0) .or. &
              (nout3dw.ne.nout3dw0) .or. &
              (nout3dxv.ne.nout3dxv0) .or. &
              (nout3dyv.ne.nout3dyv0) .or. &
              (nout3dzv.ne.nout3dzv0) .or. &
              (noutdiags.ne.noutdiags0) .or. &
              (noutdiagw.ne.noutdiagw0) .or. &
              (noutazims.ne.noutazims0) .or. &
              (noutazimw.ne.noutazimw0) .or. &
              (noutparcel.ne.noutparcel0) .or. &
              (nouthifrq.ne.nouthifrq0)   )then
            increment = .true.
            nconfig = nconfig+1
          endif
          print *
          print *,' nconfig    = ',nconfig
          print *
        endif  myid0a

#ifdef MPI
        call MPI_BCAST(nconfig,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(increment,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
#endif

        if( increment )then
          nstatout = 1
          nrec = 1
          nwrite = 1
          srec = 1
          sirec = 1
          urec = 1
          vrec = 1
          wrec = 1
          xrec = 1
          yrec = 1
          zrec = 1
          mrec = 1
          prec = 1
          nwritet = 1
          trecs = 1
          trecw = 1
          nwriteh = 1
          nwritea = 1
          arecs = 1
          arecw = 1
          nwritemv = 1
          mvrec = 1
        endif

        noutcheck = .false.
        dosolve = .false.

        if( myid.eq.0 )then
          print *
          print *,'                  end noutcheck                                  '
          print *,'  -------------------------------------------------------------  '
          print *
        endif
      else
        startup = .false.
        dorestart = .false.
        interp_on_restart = .false.
        dosolve     = .true.
        dorad       = .false.
        dostat      = .false.
        dowriteout  = .false.
        doprclout   = .false.
        dotdwrite   = .false.
        doazimwrite = .false.
        dohifrqwrite = .false.
        restarted = .false.
      endif


      writeit = .true.
      startup = .false.


    ENDDO  timeloop
    !c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c!


!----------------------------------------------------------------------

#ifdef MPI
      call MPI_BARRIER (MPI_COMM_WORLD,ierr)
      if(timestats.ge.1) time_mpb=time_mpb+mytime()

      if(myid.eq.0.and.procfiles)then
        tend=mpi_wtime()
        print *
        print *,'Total time (s): ',tend-tstart,' seconds '
        print *
      endif

#endif
!----------------------------------------------------------------------

    IF(timestats.ge.1)THEN

#ifdef MPI
      ! for MPI runs without procfiles, average the timestat terms:
      IF(.not.procfiles)THEN
        sum = 0.0
        call MPI_REDUCE(time_sound   ,sum,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        time_sound = sum/float(numprocs)
        sum = 0.0
        call MPI_REDUCE(time_poiss   ,sum,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        time_poiss = sum/float(numprocs)
        sum = 0.0
        call MPI_REDUCE(time_advs    ,sum,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        time_advs = sum/float(numprocs)
        sum = 0.0
        call MPI_REDUCE(time_advu    ,sum,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        time_advu = sum/float(numprocs)
        sum = 0.0
        call MPI_REDUCE(time_advv    ,sum,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        time_advv = sum/float(numprocs)
        sum = 0.0
        call MPI_REDUCE(time_advw    ,sum,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        time_advw = sum/float(numprocs)
        sum = 0.0
        call MPI_REDUCE(time_divx    ,sum,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        time_divx = sum/float(numprocs)
        sum = 0.0
        call MPI_REDUCE(time_buoyan  ,sum,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        time_buoyan = sum/float(numprocs)
        sum = 0.0
        call MPI_REDUCE(time_turb    ,sum,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        time_turb = sum/float(numprocs)
        sum = 0.0
        call MPI_REDUCE(time_sfcphys ,sum,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        time_sfcphys = sum/float(numprocs)
        sum = 0.0
        call MPI_REDUCE(time_ttend   ,sum,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        time_ttend = sum/float(numprocs)
        sum = 0.0
        call MPI_REDUCE(time_cor     ,sum,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        time_cor = sum/float(numprocs)
        sum = 0.0
        call MPI_REDUCE(time_diffu   ,sum,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        time_diffu = sum/float(numprocs)
        sum = 0.0
        call MPI_REDUCE(time_rdamp   ,sum,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        time_rdamp = sum/float(numprocs)
        sum = 0.0
        call MPI_REDUCE(time_microphy,sum,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        time_microphy = sum/float(numprocs)
        sum = 0.0
        call MPI_REDUCE(time_satadj  ,sum,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        time_satadj = sum/float(numprocs)
        sum = 0.0
        call MPI_REDUCE(time_dbz,sum,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        time_dbz = sum/float(numprocs)
        sum = 0.0
        call MPI_REDUCE(time_fall    ,sum,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        time_fall = sum/float(numprocs)
        sum = 0.0
        call MPI_REDUCE(time_rad     ,sum,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        time_rad = sum/float(numprocs)
        sum = 0.0
        call MPI_REDUCE(time_pbl     ,sum,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        time_pbl = sum/float(numprocs)
        sum = 0.0
        call MPI_REDUCE(time_stat    ,sum,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        time_stat = sum/float(numprocs)
        sum = 0.0
        call MPI_REDUCE(time_cflq    ,sum,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        time_cflq = sum/float(numprocs)
        sum = 0.0
        call MPI_REDUCE(time_bc      ,sum,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        time_bc = sum/float(numprocs)
        sum = 0.0
        call MPI_REDUCE(time_integ   ,sum,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        time_integ = sum/float(numprocs)
        sum = 0.0
        call MPI_REDUCE(time_write   ,sum,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        time_write = sum/float(numprocs)
        sum = 0.0
        call MPI_REDUCE(time_restart ,sum,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        time_restart = sum/float(numprocs)
        sum = 0.0
        call MPI_REDUCE(time_misc    ,sum,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        time_misc = sum/float(numprocs)
        sum = 0.0
        call MPI_REDUCE(time_swath   ,sum,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        time_swath = sum/float(numprocs)
        sum = 0.0
        call MPI_REDUCE(time_pdef    ,sum,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        time_pdef = sum/float(numprocs)
        sum = 0.0
        call MPI_REDUCE(time_prsrho  ,sum,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        time_prsrho = sum/float(numprocs)
        sum = 0.0
        call MPI_REDUCE(time_diag,sum,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        time_diag = sum/float(numprocs)
        sum = 0.0
        call MPI_REDUCE(time_azimavg ,sum,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        time_azimavg = sum/float(numprocs)
        sum = 0.0
        call MPI_REDUCE(time_hifrq ,sum,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        time_hifrq = sum/float(numprocs)
        sum = 0.0
        call MPI_REDUCE(time_ercyl ,sum,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        time_ercyl = sum/float(numprocs)
        sum = 0.0
        call MPI_REDUCE(time_tavg  ,sum,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        time_tavg = sum/float(numprocs)
        sum = 0.0
        call MPI_REDUCE(time_cuprm ,sum,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        time_cuprm = sum/float(numprocs)
        sum = 0.0
        call MPI_REDUCE(time_parcels ,sum,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        time_parcels = sum/float(numprocs)
        sum = 0.0
        call MPI_REDUCE(time_hailtrajs ,sum,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        time_hailtrajs = sum/float(numprocs)
        sum = 0.0
        call MPI_REDUCE(time_mpu1    ,sum,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        time_mpu1 = sum/float(numprocs)
        sum = 0.0
        call MPI_REDUCE(time_mpv1    ,sum,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        time_mpv1 = sum/float(numprocs)
        sum = 0.0
        call MPI_REDUCE(time_mpw1    ,sum,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        time_mpw1 = sum/float(numprocs)
        sum = 0.0
        call MPI_REDUCE(time_mpp1    ,sum,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        time_mpp1 = sum/float(numprocs)
        sum = 0.0
        call MPI_REDUCE(time_mpu2    ,sum,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        time_mpu2 = sum/float(numprocs)
        sum = 0.0
        call MPI_REDUCE(time_mpv2    ,sum,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        time_mpv2 = sum/float(numprocs)
        sum = 0.0
        call MPI_REDUCE(time_mpw2    ,sum,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        time_mpw2 = sum/float(numprocs)
        sum = 0.0
        call MPI_REDUCE(time_mpp2    ,sum,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        time_mpp2 = sum/float(numprocs)
        sum = 0.0
        call MPI_REDUCE(time_mps1    ,sum,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        time_mps1 = sum/float(numprocs)
        sum = 0.0
        call MPI_REDUCE(time_mps3    ,sum,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        time_mps3 = sum/float(numprocs)
        sum = 0.0
        call MPI_REDUCE(time_mpq1    ,sum,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        time_mpq1 = sum/float(numprocs)
        sum = 0.0
        call MPI_REDUCE(time_mps2    ,sum,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        time_mps2 = sum/float(numprocs)
        sum = 0.0
        call MPI_REDUCE(time_mps4    ,sum,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        time_mps4 = sum/float(numprocs)
        sum = 0.0
        call MPI_REDUCE(time_mpq2    ,sum,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        time_mpq2 = sum/float(numprocs)
        sum = 0.0
        call MPI_REDUCE(time_mptk1   ,sum,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        time_mptk1 = sum/float(numprocs)
        sum = 0.0
        call MPI_REDUCE(time_mptk2   ,sum,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        time_mptk2 = sum/float(numprocs)
        sum = 0.0
        call MPI_REDUCE(time_mpb     ,sum,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        time_mpb = sum/float(numprocs)
      ENDIF
#endif

      time_solve=time_sound+time_poiss+time_buoyan+time_turb+             &
                  time_diffu+time_microphy+time_stat+time_cflq+           &
                  time_bc+time_misc+time_integ+time_rdamp+time_divx+      &
                  time_write+time_restart+time_ttend+time_cor+time_fall+  &
                  time_satadj+time_dbz+time_sfcphys+time_parcels+         &
                  time_hailtrajs+                                         &
                  time_rad+time_pbl+time_swath+time_pdef+time_prsrho+     &
                  time_diag+time_azimavg+time_hifrq+time_ercyl+time_tavg+time_cuprm+ &
#ifdef MPI
                  time_mpu1+time_mpv1+time_mpw1+time_mpp1+                &
                  time_mpu2+time_mpv2+time_mpw2+time_mpp2+                &
                  time_mps1+time_mps3+time_mpq1+time_mptk1+                         &
                  time_mps2+time_mps4+time_mpq2+time_mptk2+time_mpb+                &
#endif
                  time_advs+time_advu+time_advv+time_advw
      time_solve0 = time_solve

#ifdef MPI
      mp_total=time_mpu1+time_mpv1+time_mpw1+time_mpp1+                   &
               time_mpu2+time_mpv2+time_mpw2+time_mpp2+                   &
               time_mps1+time_mps3+time_mpq1+time_mptk1+                            &
               time_mps2+time_mps4+time_mpq2+time_mptk2+time_mpb
#endif

      if(dowr) write(outfile,*)
      if(dowr) write(outfile,*) 'Total time: ',sngl(time_solve),' seconds '
      if(dowr) write(outfile,*)
      time_solve=0.01*time_solve
      if(time_solve.lt.0.0001) time_solve=1.

    IF(dowr)THEN
      write(outfile,100) 'sound   ',time_sound,time_sound/time_solve
      write(outfile,100) 'poiss   ',time_poiss,time_poiss/time_solve
      write(outfile,100) 'advs    ',time_advs,time_advs/time_solve
      write(outfile,100) 'advu    ',time_advu,time_advu/time_solve
      write(outfile,100) 'advv    ',time_advv,time_advv/time_solve
      write(outfile,100) 'advw    ',time_advw,time_advw/time_solve
      write(outfile,100) 'divx    ',time_divx,time_divx/time_solve
      write(outfile,100) 'buoyan  ',time_buoyan,time_buoyan/time_solve
      write(outfile,100) 'turb    ',time_turb,time_turb/time_solve
      write(outfile,100) 'sfcphys ',time_sfcphys,time_sfcphys/time_solve
      write(outfile,100) 'ttend   ',time_ttend,time_ttend/time_solve
      write(outfile,100) 'cor     ',time_cor,time_cor/time_solve
      write(outfile,100) 'diffu   ',time_diffu,time_diffu/time_solve
      write(outfile,100) 'rdamp   ',time_rdamp,time_rdamp/time_solve
      write(outfile,100) 'microphy',time_microphy,time_microphy/time_solve
      write(outfile,100) 'satadj  ',time_satadj,time_satadj/time_solve
      write(outfile,100) 'dbz     ',time_dbz,time_dbz/time_solve
      write(outfile,100) 'fallout ',time_fall,time_fall/time_solve
      write(outfile,100) 'radiatio',time_rad,time_rad/time_solve
      write(outfile,100) 'pbl     ',time_pbl,time_pbl/time_solve
      write(outfile,100) 'stat    ',time_stat,time_stat/time_solve
      write(outfile,100) 'cflq    ',time_cflq,time_cflq/time_solve
      write(outfile,100) 'bc      ',time_bc,time_bc/time_solve
      write(outfile,100) 'integ   ',time_integ,time_integ/time_solve
      write(outfile,100) 'write   ',time_write,time_write/time_solve
      write(outfile,100) 'restart ',time_restart,time_restart/time_solve
      write(outfile,100) 'misc    ',time_misc,time_misc/time_solve
      write(outfile,100) 'swaths  ',time_swath,time_swath/time_solve
      write(outfile,100) 'pdef    ',time_pdef,time_pdef/time_solve
      write(outfile,100) 'prsrho  ',time_prsrho,time_prsrho/time_solve
      write(outfile,100) 'parcels ',time_parcels,time_parcels/time_solve
      write(outfile,100) 'hailtraj',time_hailtrajs,time_hailtrajs/time_solve
      if( dodomaindiag )  &
      write(outfile,100) 'domndiag',time_diag,time_diag/time_solve
      if( doazimavg )  &
      write(outfile,100) 'azimavg ',time_azimavg,time_azimavg/time_solve
      if( dohifrq )  &
      write(outfile,100) 'hifrq   ',time_hifrq,time_hifrq/time_solve
      if( do_recycle )  &
      write(outfile,100) 'recycle ',time_ercyl,time_ercyl/time_solve
      if( dotimeavg )  &
      write(outfile,100) 'timeavg ',time_tavg,time_tavg/time_solve
      if( cuparam.ge.1 )  &
      write(outfile,100) 'cuparam ',time_cuprm,time_cuprm/time_solve
#ifdef MPI
      write(outfile,100) 'mp_total',mp_total,mp_total/time_solve
      write(outfile,*)
      write(outfile,100) 'mpu1    ',time_mpu1,time_mpu1/time_solve
      write(outfile,100) 'mpu2    ',time_mpu2,time_mpu2/time_solve
      write(outfile,100) 'mpv1    ',time_mpv1,time_mpv1/time_solve
      write(outfile,100) 'mpv2    ',time_mpv2,time_mpv2/time_solve
      write(outfile,100) 'mpw1    ',time_mpw1,time_mpw1/time_solve
      write(outfile,100) 'mpw2    ',time_mpw2,time_mpw2/time_solve
      write(outfile,100) 'mpp1    ',time_mpp1,time_mpp1/time_solve
      write(outfile,100) 'mpp2    ',time_mpp2,time_mpp2/time_solve
      write(outfile,100) 'mps1    ',time_mps1,time_mps1/time_solve
      write(outfile,100) 'mps2    ',time_mps2,time_mps2/time_solve
      write(outfile,100) 'mps3    ',time_mps3,time_mps3/time_solve
      write(outfile,100) 'mps4    ',time_mps4,time_mps4/time_solve
      write(outfile,100) 'mpq1    ',time_mpq1,time_mpq1/time_solve
      write(outfile,100) 'mpq2    ',time_mpq2,time_mpq2/time_solve
      write(outfile,100) 'mptk1   ',time_mptk1,time_mptk1/time_solve
      write(outfile,100) 'mptk2   ',time_mptk2,time_mptk2/time_solve
      write(outfile,100) 'mpb     ',time_mpb,time_mpb/time_solve
#endif
      write(outfile,*)
    ENDIF

100   format(3x,a8,' :  ',f10.2,2x,f6.2,'%')

      if(myid.eq.0) print *
      time_solve0 = max( time_solve0 , 1.0e-20 )
      if(myid.eq.0) write(outfile,101) (nx*ny*nz)*float(nstep-nstep0)/((time_solve0/3600.0)*float(numprocs*nomp))
101   format('  MN = ',es12.4e2)
      if(myid.eq.0) print *
      if(myid.eq.0) print *,'  approximate core-hours = ',sngl((time_solve0/3600.0)*float(numprocs*nomp))
      if(myid.eq.0) print *
#ifdef MPI
      call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif

    ENDIF

!  End time loop
!----------------------------------------------------------------------

      close(unit=51)
      close(unit=52)
      close(unit=53)
      close(unit=54)
      close(unit=60)

!----------------------------------------------------------------------

#ifdef MPI
      if(myid.eq.0)then
        print *,'Program terminated normally:  myid=',myid
        do n=1,(numprocs-1)
          call MPI_IRECV(doit,1,mpi_logical,n,999,MPI_COMM_WORLD,reqs,ierr)
          call MPI_WAIT(reqs,status,ierr)
          print *,'Program terminated normally:  myid=',n
        enddo
      else
        doit = .true.
        call MPI_ISEND(doit,1,mpi_logical,0,999,MPI_COMM_WORLD,reqs,ierr)
        call MPI_WAIT(reqs,status,ierr)
      endif
      call MPI_FINALIZE(ierror=rc)
#else
      print *,'Program terminated normally'
#endif

      stop

8778  print *
      print *,'  8778: error opening namelist.input '
      print *,'    ... stopping cm1 ... '
      print *
      call stopcm1

8800  print *
      print *,'  8800: error reading param0 section of namelist.input '
      print *,'    ... stopping cm1 ... '
      call stopcm1

8802  print *
      print *,'  8802: error reading param2 section of namelist.input '
      print *,'    ... stopping cm1 ... '
      call stopcm1

8808  print *
      print *,'  8808: error reading param8 section of namelist.input '
      print *,'    ... stopping cm1 ... '
      call stopcm1

      end program cm1

