cJuly 22, 2013 Jx and |B| and bfiled lines in 300 x 100 x 100 (300 < x 600)
c  y-z plane Jx with Byz
c  larger number of particles, az is remmoved for thin jet 
c        September 7, 2002
c  Figure is rotated 90 degree cw with trimmed edges
c  contour: Ne: dense, Ni: densi, Jy: ajy, Eke: denqe, Eki: denqi
c           Btot: btot, Ey: aey, Etot: etot
c           GAMMA: vgam
c  arrows:  Bzx, emyx, emyz, Ezx; elyx, elyz
c           Vezx: fvlex, fvley, Vizx; fvlix, fvliy
c
c Contours with color for number, energy, flux, velocity, density
c with magnetic and electric fields
c  
      parameter (Nproc=81)
c     parameter (Npx=80,Npy=2,Npz=2)
      parameter (mxs=645,mys=131,mzs=131)

c     parameter (nptl=4600000)
c     parameter (mb=nptl,mj=3500000)
      parameter (nptl=2000000)
      parameter (mb=nptl,mj=1800000)


      parameter (mg = mb*Nproc,mjg = mj*Nproc)
c     integer*8,parameter:: mg=mb*Nproc
c     integer*8,parameter:: mjg=mj*Nproc


c     parameter (imax=32767)
c  attention 90 dgree rotation
c     parameter (mxx=mxs-12,myy=mys-12,mzz=mzs-12)
c New ken
      parameter (mxx=200,myy=mys-5,mzz=mzs-5)
      parameter( pi=3.141592 )
c     parameter( imax=200, jmax=200, kmax=200 )
c     parameter( lmax=200, mmax=200, nmax=200 )
      parameter( imax=300, jmax=100, kmax=100 )
      parameter( lmax=300, mmax=100, nmax=100 )
      parameter( nq=12, ns=7, nv=2, nnh=10 )
c
      parameter( ncol=40, nrow=20 )
c ken
      parameter(llmax=300,mmmax=100,nnmax=100)


c attention
c     parameter (mxr1=52,   myr1=27,   mzr1=27)
c     parameter (mxr=53,  myr=28,   mzr=28)
c     parameter (mxq=28,   myq=53,   mzq=53)
c     parameter (mxr1=50,   myr1=25,   mzr1=25)
c     parameter (mxr=50,  myr=25,   mzr=25)
c     parameter (mxq=50,   myq=25,   mzq=25)
      parameter (mxr1=40,   myr1=20,   mzr1=20)
      parameter (mxr=40,  myr=20,   mzr=20)
      parameter (mxq=40,  myq=20,   mzq=20)
c     real*8 me,mi,mi2
      dimension rhi(mxs,mys,mzs),rhe(mxs,mys,mzs)
      dimension eki(mxs,mys,mzs),eke(mxs,mys,mzs)
      dimension fex(mxs,mys,mzs),fey(mxs,mys,mzs),fez(mxs,mys,mzs)
      dimension fix(mxs,mys,mzs),fiy(mxs,mys,mzs),fiz(mxs,mys,mzs)
      dimension rhij(mxs,mys,mzs),rhej(mxs,mys,mzs)
      dimension ekij(mxs,mys,mzs),ekej(mxs,mys,mzs)
      dimension fexj(mxs,mys,mzs),feyj(mxs,mys,mzs),fezj(mxs,mys,mzs)
      dimension fixj(mxs,mys,mzs),fiyj(mxs,mys,mzs),fizj(mxs,mys,mzs)

c  for 3DRPIC for VisIt
      real*4 avsdna(llmax,nnmax,mmmax),avsbxa(llmax,nnmax,mmmax),
     &       avsbya(llmax,nnmax,mmmax),avsbza(llmax,nnmax,mmmax),
     &       avspra(llmax,nnmax,mmmax),avsvxa(llmax,nnmax,mmmax),
     &       avsvya(llmax,nnmax,mmmax),avsvza(llmax,nnmax,mmmax),
     &       avsxxa(llmax,nnmax,mmmax),avsxya(llmax,nnmax,mmmax),
     &       avsxza(llmax,nnmax,mmmax)


c     dimension ixx(mb),iyy(mb),izz(mb)
c  ambient ions and electrons positions and velocities arrays
c     real*4 xi(mb),yi(mb),zi(mb),ui(mb),vi(mb),wi(mb)
c     real*4 xe(mb),ye(mb),ze(mb),ue(mb),ve(mb),we(mb)

c  jet ions and electrons positions and velocities arrays
c     real*4 xij(mj),yij(mj),zij(mj),uij(mj),vij(mj),wij(mj)
c     real*4 xej(mj),yej(mj),zej(mj),uej(mj),vej(mj),wej(mj)

c  global field arrays
c     real bxx,byy,bzz,exx,eyy,ezz
      real*4 bxx(mxs,mys,mzs),byy(mxs,mys,mzs),bzz(mxs,mys,mzs)
      real*4 exx(mxs,mys,mzs),eyy(mxs,mys,mzs),ezz(mxs,mys,mzs)
c  global particle arrays
c     real xig,yig,zig,uig,vig,wig
c     real xeg,yeg,zeg,ueg,veg,weg
c     real xijg,yijg,zijg,uijg,vijg,wijg
c     real xejg,yejg,zejg,uejg,vejg,wejg
c     real*4 xig(mg),yig(mg),zig(mg),uig(mg),vig(mg),wig(mg)
c     real*4 xeg(mg),yeg(mg),zeg(mg),ueg(mg),veg(mg),weg(mg)
c     real*4 xijg(mjg),yijg(mjg),zijg(mjg)
c     real*4 uijg(mjg),vijg(mjg),wijg(mjg)
c     real*4 xejg(mjg),yejg(mjg),zejg(mjg)
c     real*4 uejg(mjg),vejg(mjg),wejg(mjg)

      real, allocatable, dimension(:) :: xig,yig,zig,uig,vig,wig
      real, allocatable, dimension(:) :: xeg,yeg,zeg,ueg,veg,weg
      real, allocatable, dimension(:) :: xijg,yijg,zijg
      real, allocatable, dimension(:) :: uijg,vijg,wijg
      real, allocatable, dimension(:) :: xejg,yejg,zejg
      real, allocatable, dimension(:) :: uejg,vejg,wejg


      real*4 
     &   c4,DT4,qi4,qe4,qmi4,qme4,vithml4,vethml4,vijet4,vejet4,
     &   vithmj4,vethmj4,refli4,refle4,rselect4,xj04,yj04,zj04,
     &   b0x4,b0y4,e0z4
c
c  attention myy mzz
      dimension dense(myy,mzz),densi(myy,mzz),ajy(myy,mzz),
     1          btot(myy,mzz),etot(myy,mzz),vgam(myy,mzz),
     1          denqe(myy,mzz),denqi(myy,mzz),aey(myy,mzz)
      dimension fvlex(mxr,mzr),fvley(mxr,mzr),fvlix(mxr,mzr),
     1   fvliy(mxr,mzr),fvqex(mxq,mzq),fvqey(mxq,mzq),
     1   fvqix(mxq,mzq),fvqiy(mxq,mzq)
c     dimension argx(mxs),argy(mys)
c     dimension xeb(mbeam), yeb(mbeam),
c    1          zeb(mbeam), vzb(mbeam),
c    2          xeg0(mback), yeg0(mback),
c    3          zeg0(mback), vzg(mback)
c
      dimension rho(mxs,mys),rho1(mxs,mzs)
      dimension gridr(4),gridm(4),gride(4),bnum(2),bnu1(2) 
      dimension xd(3),yd(3)                                             
      data gridr/0.10,0.85,0.100,0.85/
      data gride/0.10,0.75,0.25,0.7/,bnum/0.030,0.030/
      data gridm/0.25,0.95,0.25,0.5478/,bnu1/0.030,0.030/
      data xd/11.,12.,15./, yd/1.0,2.0,3.0/                             
c                                                                       
c     dimension x(nptl),y(nptl),z(nptl),
c    1          u(nptl),v(nptl),w(nptl)
c     dimension     mcont(mbeam),kcont(mback),
c    1              mcon2(mbeam),kcon2(mback),
c    1              mcone(mlpar),kconi(mlpar)
      COMMON /VEC1/   ASH        ,EXT        ,ICTRFG     ,ILAB       ,
     +                IOFFD      ,IOFFM      ,ISX        ,ISY        ,
     +                RMN        ,RMX        ,SIDE       ,SIZE       ,
     +                XLT        ,YBT        ,ZMN        ,ZMX
      common /srfip1/ ifr    ,istp   ,irots  ,idrx ,idry,               
     1                idrz   ,iupper ,iskirt ,ncla ,theta,              
     2                hskirt ,chi    ,clo    ,cinc ,ispval              
c------------------------------------------------------------
      character num*3,num1,num2,num3
      character step*6,step1,step2,step3,step4,step5
      character st01,st02,st03,st0*3,st*4,hyph
      character dchar1*6,dchar2*4,file1*10,file2*8
      character strpd*6
      character message2*6
      character dchar5*10
      character file8*19
      character dot,vk*3,vvt*4

      character message*28                                              
      character message1*22
      character*1 ittmp(21)
!     logical, intent(in) :: ascii ! if T writes ascii version
      logical :: ascii ! if T writes ascii version
!     character(LEN=20),intent(in) :: filename ! name of VTK file
      character(LEN=35):: filename ! name of VTK file
      open(6,file='dwriteot',status='unknown',form='formatted')

c     call opngks
c     call gsclip (0)

c specify "nst" - which data need to be processed
c     open(unit=2, file='nstn', status='old')

c      read(2,*) nst
      
c      nst = 5
       do ins=25,15,-1

       nst = ins
      
       c = 1.0

       strpd = 'partd_'
       hyph='_'

       st01 = char(int(nst/100.)+48)
       st02 = char(int((nst-int(nst/100.)*100)/10.)+48)
       st03 = char(int(nst-int(nst/10.)*10)+48)
       st0 = st01//st02//st03
       st = hyph//st0

       write(6,*) 'before reading'  

       allocate(xig(mg),yig(mg),zig(mg),uig(mg),vig(mg),wig(mg))
       allocate(xeg(mg),yeg(mg),zeg(mg),ueg(mg),veg(mg),weg(mg))
       allocate(xijg(mjg),yijg(mjg),zijg(mjg))
       allocate(uijg(mjg),vijg(mjg),wijg(mjg))
       allocate(xejg(mjg),yejg(mjg),zejg(mjg))
       allocate(uejg(mjg),vejg(mjg),wejg(mjg))


c  read from file
        dchar1='bfield'
        file1= dchar1//st

c for visitdat 
        dchar5='vJxBlLn03n'
        dot='.'
        vk='vtk'
        vvt=dot//vk
        file8=dchar5//st//vvt   
c       file8=dchar5//st


        open(15,file=file1,form='unformatted')
c       open(15,file=bfield_007,form='unformatted')
        read(15) bxx,byy,bzz
        close(15)
       write(6,*) 'after reading 15'  
       write(6,*)'bxx(300,35,35),byy(300,35,35),bzz(300,35,35)=',
     1 bxx(300,35,35),byy(300,35,35),bzz(300,35,35)

        dchar1='efield'
        file1= dchar1//st

        open(16,file=file1,form='unformatted')
c       open(16,file=efield_007,form='unformatted')
        read(16) exx,eyy,ezz
        close(16)
       write(6,*)'exx(300,35,35),eyy(300,35,35),ezz(300,35,35)=',
     1 exx(300,35,35),eyy(300,35,35),ezz(300,35,35)

       write(6,*) 'after reading 16'  
        dchar2='ions'
        file2= dchar2//st

cheung
        open(17,file=file2,form='unformatted')
c       open(17,file=ions_007,form='unformatted')
        read(17) iong
c       read(17) xig(1:iong),yig(1:iong),zig(1:iong)
c       read(17) uig(1:iong),vig(1:iong),wig(1:iong)
        read(17) xig(1:iong)
        read(17) yig(1:iong)
        read(17) zig(1:iong)
        read(17) uig(1:iong)
        read(17) vig(1:iong)
        read(17) wig(1:iong)

        close(17)
       write(6,*) 'after reading 17,iong=',iong

        dchar2='lecs'
        file2= dchar2//st
        
cheung
        open(18,file=file2,form='unformatted')
c       open(18,file=lecs_007,form='unformatted')
        read(18) lecg
c       read(18) xeg(1:lecg),yeg(1:lecg),zeg(1:lecg)
c       read(18) ueg(1:lecg),veg(1:lecg),weg(1:lecg)
        read(18) xeg(1:lecg)
        read(18) yeg(1:lecg)
        read(18) zeg(1:lecg)
        read(18) ueg(1:lecg)
        read(18) veg(1:lecg)
        read(18) weg(1:lecg)

        close(18)
       write(6,*) 'after reading 18,lecg=',lecg
        
        dchar2='ionj'
        file2= dchar2//st

cheung
        open(19,file=file2,form='unformatted')
c       open(19,file=ionj_007,form='unformatted')
        read(19) ionjg
c       read(19) xijg(1:ionjg),yijg(1:ionjg),zijg(1:ionjg)
c       read(19) uijg(1:ionjg),vijg(1:ionjg),wijg(1:ionjg)
        read(19) xijg(1:ionjg)
        read(19) yijg(1:ionjg)
        read(19) zijg(1:ionjg)
        read(19) uijg(1:ionjg)
        read(19) vijg(1:ionjg)
        read(19) wijg(1:ionjg)

        close(19) 
       write(6,*) 'after reading 19,ionjg=',ionjg

        dchar2='lecj'
        file2= dchar2//st

cheung
        open(20,file=file2,form='unformatted')
c       open(20,file=lecj_007,form='unformatted')
        read(20) lecjg
c       read(20) xejg(1:lecjg),yejg(1:lecjg),zejg(1:lecjg)
c       read(20) uejg(1:lecjg),vejg(1:lecjg),wejg(1:lecjg)
        read(20) xejg(1:lecjg)
        read(20) yejg(1:lecjg)
        read(20) zejg(1:lecjg)
        read(20) uejg(1:lecjg)
        read(20) vejg(1:lecjg)
        read(20) wejg(1:lecjg)

        close(20)
       write(6,*) 'after reading 20,lecjg=',lecjg

        dchar1='kenout'
        file1= dchar1//st

cheung
c       open(20,file=file2,form='unformatted')

c       message='kenout'

        open(21,file=file1,form='unformatted')
c       message2='kenout'
c       open(21,file=message2,form='unformatted')
        read(21) nstep
     &  ,c,DT,qi,qe,qmi,qme,vithml,vethml,vijet,vejet,
c    & vithmj,vethmj,refli,refle,rselect,xj0,yj0,zj0,
c    &  vithmj,vethmj,refli,refle,rselect,xj0,yj0,zj0,b0x,
     &  vithmj,vethmj
c new Jacek
c    &  b0y,e0z
        close(21)
       write(6,*) 'after reading 21,nstep=',nstep


      time = 0.1*(nstep -1)

c     write(6,*) ions, lecs,maxptl,maxhlf
c     write(6,*) 'ionm,ionj,lecm,lecj= ',ionm,ionj,lecm,lecj
      je=47
c     ke=48
c     ie=70
      in=500
      jn=100
      kn=100
      imid=in
      jmid=jn
      kmid=kn
c
      b0z=b0z*0.992
c
c   shifting the location
c
c      ishift = vejet*DT*nstep -350
       ishift = 0         
c     if(ishift.lt.0) ishift = 0
c     if((ishift+408).ge.4000) ishift = 4000-408
      write(6,*)'DT,vejet,nstep,ishift=',DT,vejet,nstep,ishift

c     xliml=150.
c     xlimr=290.
c     xliml=0.
c     xlimr=150.
      do 8 k=1,mzs
      do 8 j=1,mys
      do 8 i=1,mxs
      rhi(i,j,k)=0.
      rhe(i,j,k)=0.
      fex(i,j,k)=0.
      fey(i,j,k)=0.
      fez(i,j,k)=0.
      fix(i,j,k)=0.
      fiy(i,j,k)=0.
      fiz(i,j,k)=0.
      rhij(i,j,k)=0.
      rhej(i,j,k)=0.
      fexj(i,j,k)=0.
      feyj(i,j,k)=0.
      fezj(i,j,k)=0.
      fixj(i,j,k)=0.
      fiyj(i,j,k)=0.
      fizj(i,j,k)=0.
    8 continue

C  Volume weighting:
c for ions
      call Vol_Weighting(mxs,mys,mzs,iong,rhi,eki,uig,vig,wig,
     &                         xig,yig,zig,fix,fiy,fiz)
c
       write(6,*)'after background ions'

c for electrons
      call Vol_Weighting(mxs,mys,mzs,lecg,rhe,eke,ueg,veg,weg,
     &                         xeg,yeg,zeg,fex,fey,fez)
       write(6,*)'after background electrons'
        write(6,*)' non-zero rhe, zero rhej'
c     write(6,*)'rhe(350,25,25),rhej(350,25,25)=',
c    1 rhe(350,25,25),rhej(350,25,25)

C  Volume weighting:
c for jet ions
      call Vol_Weighting(mxs,mys,mzs,ionjg,rhij,ekij,uijg,vijg,wijg,
     &                         xijg,yijg,zijg,fixj,fiyj,fizj)
      write(6,*)'after jet ions'

c
c for jet electrons
      call Vol_Weighting(mxs,mys,mzs,lecjg,rhej,ekej,uejg,vejg,wejg,
     &                         xejg,yejg,zejg,fexj,feyj,fezj)

cheung
      write(6,*)'after jet electrons'

        write(6,*)' non-zero rhe, non-zero rhej'
      write(6,*)'rhe(350,25,25),rhej(350,25,25)=',
     1 rhe(350,25,25),rhej(350,25,25)
c     write(6,*)"s cx dx cy dy cz dz", s,cx,dx,cy,dy,cz,dz
c
c  attention the sustem is original size  
c
        do 534 m=1,mmax
         do 534 n=1,nmax
          do 534 l=1,lmax
           avsxxa(l,n,m) = float(l)
           avsxya(l,n,m) = float(n)
           avsxza(l,n,m) = float(m)
  534   continue
        do 535 m=1,mmax
         do 535 n=1,nmax
          do 535 l=1,lmax
           avsxxa(l,n,m) = avsxxa(l,1,1)
           avsxya(l,n,m) = avsxya(1,n,1)
           avsxza(l,n,m) = avsxza(1,1,m)
  535   continue

c      dchar5//st
       ascii=.true.
c        filename='visitdata'
         filename=file8


c attention z --->x  x --> z
c     do 11 j=1,mzz
c        ii=j+2
c     do 11 i=1,myy
c        jj=i+2
c     
c     densi(i,j)=rhi(in,jj,ii)+rhij(in,jj,ii)
c     dense(i,j)=rhe(in,jj,ii)+rhej(in,jj,ii)
c     btot(i,j)=sqrt(
c    1  (0.5*(bxx(in,jj,ii) +bxx(jj,jj-1,ii-1)))**2
c    1 +(0.5*(byy(in,jj,ii) +byy(jj-1,jj,ii-1)))**2
c    1 +(0.5*(bzz(in,jj,ii) +bzz(jj-1,jj-1,ii)))**2)
c     etot(i,j)=sqrt(
c    1  (0.5*(exx(in,jj,ii) +exx(jj-1,jj,ii)))**2
c    1 +(0.5*(eyy(in,jj,ii) +eyy(jj,jj-1,ii)))**2
c    1 +(0.5*(ezz(in,jj,ii) +ezz(jj,jj,ii-1)))**2)
c attention   ajy = Jx
c     ajy(i,j)=fiy(in,jj,ii)-fey(j,jj,ii)
c    ajy = By(x,z)
c     ajy(i,j)= 0.5*(
c    1 +0.5*(byy(jj,jn,ii) +byy(jj-1,jn,ii))
c    1 +0.5*(byy(jj,jn,ii-1) +byy(jj-1,jn,ii-1)))
c     vtemp= 1.0 -(fex(in,jj,ii)**2+fey(in,jj,ii)**2
c    1      +fez(in,jj,ii)**2)/(c*rhe(in,jj,ii))**2
c    1           -(fexj(in,jj,ii)**2+feyj(in,jj,ii)**2
c    1      +fezj(in,jj,ii)**2)/(c*rhej(in,jj,ii))**2
c     if(vtemp.le.0.0) vtemp=0.0001
c     vgam(i,j)=1.0/sqrt(vtemp)
c  11 continue

         umin = 0.0
         umax = 200.
         vmin = 0.0
         vmax = 200.
         tmin = 0.0
         tmax = 200.
         write(6,*)  umin, umax, vmin, vmax, tmin, tmax
        do 234 m=1,mmax
               mm = m +15
         do 234 n=1,nmax
               nn = n +15
          do 234 l=1,lmax
               ll = l +319
           avsdna(l,n,m) =sqrt(
     1  (0.5*(bxx(ll,nn,mm) +bxx(ll,nn-1,mm-1)))**2
     1 +(0.5*(byy(ll,nn,mm) +byy(ll-1,nn,mm-1)))**2
     1 +(0.5*(bzz(ll,nn,mm) +bzz(ll-1,nn-1,mm)))**2)

           avsvxa(l,n,m) = fix(ll,nn,mm)-fex(ll,nn,mm)
     1        +fixj(ll,nn,mm)-fexj(ll,nn,mm)
           avsvya(l,n,m) = fiy(ll,nn,mm)-fey(ll,nn,mm)
     1        +fiyj(ll,nn,mm)-feyj(ll,nn,mm)
           avsvza(l,n,m) =  fiz(ll,nn,mm)-fez(ll,nn,mm)
     1        +fizj(ll,nn,mm)-fezj(ll,nn,mm)
           avsbxa(l,n,m) = 0.5*(bxx(ll,nn,mm) +bxx(ll,nn-1,mm-1))
           avsbya(l,n,m) = 0.5*(byy(ll,nn,mm) +byy(ll-1,nn,mm-1))
           avsbza(l,n,m) = 0.5*(bzz(ll,nn,mm) +bzz(ll-1,nn-1,mm))
           avspra(l,n,m) =rhe(ll,nn,mm)+rhej(ll,nn,mm)
c          avspra(l,n,m) =sqrt(
c    1  (0.5*(exx(ll,nn,mm) +exx(ll-1,nn,mm)))**2
c    1 +(0.5*(eyy(ll,nn,mm) +eyy(ll,nn-1,mm)))**2
c    1 +(0.5*(ezz(ll,nn,mm) +ezz(ll,nn,mm-1)))**2)
 
  234   continue

c  for AVS

        write(6,*) 'after calculating data'

c       open(unit=401,file='avsro',status='unknown')
c       open(unit=402,file='avsvx',status='unknown')
c       open(unit=403,file='avsvy',status='unknown')
c       open(unit=404,file='avsvz',status='unknown')
c       open(unit=405,file='avsbx',status='unknown')
c       open(unit=406,file='avsby',status='unknown')
c       open(unit=407,file='avsbz',status='unknown')
c       open(unit=408,file='avspr',status='unknown')
c       open(unit=409,file='avsc1',status='unknown')
c       open(unit=410,file='avsc2',status='unknown')
c       open(unit=411,file='avsc3',status='unknown')
c       open(unit=412,file='avstm',status='unknown')

        write(6,*) 'after defining fine names'

c       do 334 m=1,mmax
c        do 334 n=1,nmax
c         do 334 l=1,lmax
c          write(401,*) avsdna(l,n,m)
c          write(402,*) avsvxa(l,n,m)
c          write(403,*) avsvya(l,n,m)
c          write(404,*) avsvza(l,n,m)
c          write(405,*) avsbxa(l,n,m)
c          write(406,*) avsbya(l,n,m)
c          write(407,*) avsbza(l,n,m)
c          write(408,*) avspra(l,n,m)
c 334   continue

c          write(409,*) avsxxa(l,1,1)
  342   continue

c       do 344 n=1,nmax
           write(410,*) avsxya(1,n,1)
  344   continue

c       do 343 m=1,mmax
           write(411,*) avsxza(1,1,m)
  343   continue


c        write(412,*) ' time =',time
c        write(412,*)  umin, umax, vmin, vmax, tmin, tmax


        write(6,*) 'after writing data in the file'

c      write(6,*)'dense(100,35),densi(100,35),ajy(100,35)=',
c    1  dense(100,35),densi(100,35),ajy(100,35)
c
c  for visit

      call dump_vtk_3d(ascii,filename,llmax,nnmax,mmmax,avsxxa,avsxya,
     & avsxza,avsdna,avspra,avsvxa,avsvya,avsvza,avsbxa,avsbya,avsbza)


c 110 format('JX      (EFLU) T=',f8.1)
c 110 format('BY      (MAG)  T=',f8.1)
c 111 format('ELE DEN (MAG)  T=',f8.1)
c 112 format('ELE DEN (EFLU) T=',f8.1)
c 112 format('ELE DEN (ELE)  T=',f8.1)
c 113 format('ION DEN (MAG)  T=',f8.1)
c 114 format('BTOT    (MAG)  T=',f8.1)
c 115 format('ETOT    (ELE)  T=',f8.1)
c 116 format('GAMMAE  (EFLU) T=',f8.1)
c 117 format('EY      (ELE)  T=',f8.1)
c 118 format('EKIN    (EFLU) T=',f8.1)
c 119 format('IKIN    (IFLU) T=',f8.1)
c 214 format('CUR DEN (MAG)  T=',f8.1)
c 215 format('CUR DEN (ELE)  T=',f8.1)
c 275 format('CUR DEN (FLUX) T=',f8.1)
c 216 format('CUR DEN (CUR)  T=',f8.1)
c 223 format('ELE ENE (MAG)  T=',f8.1)
c 224 format('ION ENE (MAG)  T=',f8.1)
c 323 format('ELE ENE (FLUX) T=',f8.1)
c 324 format('ION ENE (FLUX) T=',f8.1)
c  25 format('X-Y SPACE (ION) T=',F6.0,'$')                            
c  26 format('X-Y SPACE (ELE) T=',F6.0,'$')                            
c  27 format('Z-VZ  SPACE     T=',F6.0,'$')                            
c                                                                       
c     go to 4
c     call clsgks

c     mv visitdata vJxBlLn03n20.vtk
c     mv vistidat file8
c     mv vistidat dchar5//st
      end do
      stop
      end
c
      subroutine arwmax(kxr,kyr,fvex,fvey,fvix,fviy,
     1            flo1,hi1,flo2,hi2)
      dimension fvex(kxr,kyr),fvey(kxr,kyr), 
     1          fvix(kxr,kyr),fviy(kxr,kyr)
      flo1=1.0e10
      hi1=-1.0
      flo2=1.0e10
      hi2=-1.0
c
      do 31 i = 1, kxr
      do 31 j = 1, kyr
      exy=sqrt(fvex(i,j)**2 +fvey(i,j)**2)
      bxy=sqrt(fvix(i,j)**2 +fviy(i,j)**2)
      if(exy.le.flo1) flo1=exy
      if(exy.ge.hi1)  hi1=exy
      if(bxy.le.flo2) flo2=bxy
      if(bxy.ge.hi2)  hi2=bxy
   31 continue
c
      return
      end
C ---------------------------------------------------------------------
C
      subroutine Vol_Weighting(mxs,mys,mzs,lecg,rhe,eke,ueg,veg,weg,
     &                         xeg,yeg,zeg,fex,fey,fez)
      implicit none

      integer mxs,mys,mzs,mg,lecg
      real*4 xeg(lecg),yeg(lecg),zeg(lecg)
      real*4 ueg(lecg),veg(lecg),weg(lecg)
      real rhe(mxs,mys,mzs),eke(mxs,mys,mzs)
      real fex(mxs,mys,mzs),fey(mxs,mys,mzs),fez(mxs,mys,mzs)
      real*4 dx,dy,dz,cx,cy,cz
      real*4 sl,sr,sn,s     
      integer i,j,k,l,m,n,n0

C  Volume weighting:

c     mhlec = maxhlf +lecs
      do 3 n0=1, lecg
c     if(n0.gt.mhlec.and.n0.lt.lecm) go to 3

      i=xeg(n0)
      dx=xeg(n0)-i
      cx=1.-dx
      j=yeg(n0)
      dy=yeg(n0)-j
      cy=1.-dy
      k=zeg(n0)
      dz=zeg(n0)-k
      cz=1.-dz
C  Smoothing with the (.25,.5,.25) profile in each dimension:
      sl=.5
      do 187 l=-1,1
      sl=.75-sl
      sr=.5
      do 187 m=-1,1
      sr=.75-sr
      sn=.5
      do 187 n=-1,1
      sn=.75-sn
      s=sl*sr*sn
      rhe(i+l  ,j+m  ,k+n  )=rhe(i+l  ,j+m  ,k+n  )+s*cx*cy*cz
      rhe(i+l+1,j+m  ,k+n  )=rhe(i+l+1,j+m  ,k+n  )+s*dx*cy*cz
      rhe(i+l  ,j+m+1,k+n  )=rhe(i+l  ,j+m+1,k+n  )+s*cx*dy*cz
      rhe(i+l  ,j+m  ,k+n+1)=rhe(i+l  ,j+m  ,k+n+1)+s*cx*cy*dz
      rhe(i+l  ,j+m+1,k+n+1)=rhe(i+l  ,j+m+1,k+n+1)+s*cx*dy*dz
      rhe(i+l+1,j+m  ,k+n+1)=rhe(i+l+1,j+m  ,k+n+1)+s*dx*cy*dz
      rhe(i+l+1,j+m+1,k+n  )=rhe(i+l+1,j+m+1,k+n  )+s*dx*dy*cz
      rhe(i+l+1,j+m+1,k+n+1)=rhe(i+l+1,j+m+1,k+n+1)+s*dx*dy*dz

c
c     eke(i+l  ,j+m  ,k+n  )=eke(i+l  ,j+m  ,k+n  )
c    1 +(ueg(n0)*ueg(n0)+veg(n0)*veg(n0)+weg(n0)*weg(n0))*s*cx*cy*cz
c     eke(i+l+1,j+m  ,k+n  )=eke(i+l+1,j+m  ,k+n  )
c    1 +(ueg(n0)*ueg(n0)+veg(n0)*veg(n0)+weg(n0)*weg(n0))*s*dx*cy*cz
c     eke(i+l  ,j+m+1,k+n  )=eke(i+l  ,j+m+1,k+n  )
c    1 +(ueg(n0)*ueg(n0)+veg(n0)*veg(n0)+weg(n0)*weg(n0))*s*cx*dy*cz
c     eke(i+l  ,j+m  ,k+n+1)=eke(i+l  ,j+m  ,k+n+1)
c    1 +(ueg(n0)*ueg(n0)+veg(n0)*veg(n0)+weg(n0)*weg(n0))*s*cx*cy*dz
c     eke(i+l  ,j+m+1,k+n+1)=eke(i+l  ,j+m+1,k+n+1)
c    1 +(ueg(n0)*ueg(n0)+veg(n0)*veg(n0)+weg(n0)*weg(n0))*s*cx*dy*dz
c     eke(i+l+1,j+m  ,k+n+1)=eke(i+l+1,j+m  ,k+n+1)
c    1 +(ueg(n0)*ueg(n0)+veg(n0)*veg(n0)+weg(n0)*weg(n0))*s*dx*cy*dz
c     eke(i+l+1,j+m+1,k+n  )=eke(i+l+1,j+m+1,k+n  )
c    1 +(ueg(n0)*ueg(n0)+veg(n0)*veg(n0)+weg(n0)*weg(n0))*s*dx*dy*cz
c     eke(i+l+1,j+m+1,k+n+1)=eke(i+l+1,j+m+1,k+n+1)
c    1 +(ueg(n0)*ueg(n0)+veg(n0)*veg(n0)+weg(n0)*weg(n0))*s*dx*dy*dz
c
      fex(i+l  ,j+m  ,k+n  )=fex(i+l  ,j+m  ,k+n  )
     1        +ueg(n0)*s*cx*cy*cz
      fex(i+l+1,j+m  ,k+n  )=fex(i+l+1,j+m  ,k+n  )
     1        +ueg(n0)*s*dx*cy*cz
      fex(i+l  ,j+m+1,k+n  )=fex(i+l  ,j+m+1,k+n  )
     1        +ueg(n0)*s*cx*dy*cz
      fex(i+l  ,j+m  ,k+n+1)=fex(i+l  ,j+m  ,k+n+1)
     1        +ueg(n0)*s*cx*cy*dz
      fex(i+l  ,j+m+1,k+n+1)=fex(i+l  ,j+m+1,k+n+1)
     1        +ueg(n0)*s*cx*dy*dz
      fex(i+l+1,j+m  ,k+n+1)=fex(i+l+1,j+m  ,k+n+1)
     1        +ueg(n0)*s*dx*cy*dz
      fex(i+l+1,j+m+1,k+n  )=fex(i+l+1,j+m+1,k+n  )
     1        +ueg(n0)*s*dx*dy*cz
      fex(i+l+1,j+m+1,k+n+1)=fex(i+l+1,j+m+1,k+n+1)
     1        +ueg(n0)*s*dx*dy*dz
c
      fey(i+l  ,j+m  ,k+n  )=fey(i+l  ,j+m  ,k+n  )
     1        +veg(n0)*s*cx*cy*cz
      fey(i+l+1,j+m  ,k+n  )=fey(i+l+1,j+m  ,k+n  )
     1        +veg(n0)*s*dx*cy*cz
      fey(i+l  ,j+m+1,k+n  )=fey(i+l  ,j+m+1,k+n  )
     1        +veg(n0)*s*cx*dy*cz
      fey(i+l  ,j+m  ,k+n+1)=fey(i+l  ,j+m  ,k+n+1)
     1        +veg(n0)*s*cx*cy*dz
      fey(i+l  ,j+m+1,k+n+1)=fey(i+l  ,j+m+1,k+n+1)
     1        +veg(n0)*s*cx*dy*dz
      fey(i+l+1,j+m  ,k+n+1)=fey(i+l+1,j+m  ,k+n+1)
     1        +veg(n0)*s*dx*cy*dz
      fey(i+l+1,j+m+1,k+n  )=fey(i+l+1,j+m+1,k+n  )
     1        +veg(n0)*s*dx*dy*cz
      fey(i+l+1,j+m+1,k+n+1)=fey(i+l+1,j+m+1,k+n+1)
     1        +veg(n0)*s*dx*dy*dz
c
      fez(i+l  ,j+m  ,k+n  )=fez(i+l  ,j+m  ,k+n  )
     1        +weg(n0)*s*cx*cy*cz
      fez(i+l+1,j+m  ,k+n  )=fez(i+l+1,j+m  ,k+n  )
     1        +weg(n0)*s*dx*cy*cz
      fez(i+l  ,j+m+1,k+n  )=fez(i+l  ,j+m+1,k+n  )
     1        +weg(n0)*s*cx*dy*cz
      fez(i+l  ,j+m  ,k+n+1)=fez(i+l  ,j+m  ,k+n+1)
     1        +weg(n0)*s*cx*cy*dz
      fez(i+l  ,j+m+1,k+n+1)=fez(i+l  ,j+m+1,k+n+1)
     1        +weg(n0)*s*cx*dy*dz
      fez(i+l+1,j+m  ,k+n+1)=fez(i+l+1,j+m  ,k+n+1)
     1        +weg(n0)*s*dx*cy*dz
      fez(i+l+1,j+m+1,k+n  )=fez(i+l+1,j+m+1,k+n  )
     1        +weg(n0)*s*dx*dy*cz
      fez(i+l+1,j+m+1,k+n+1)=fez(i+l+1,j+m+1,k+n+1)
     1        +weg(n0)*s*dx*dy*dz
  187 continue
    3 continue
      return
      end subroutine Vol_Weighting

!        
! Given a filename of VTK file and data arrays, this routine writes the
! the data in vkt format. This is based on dump_vtk routine.
!  
! Note: time of dump is not returned (this will be implemented later).
!
!
! Writes stuructured_grid data  (For 3-D case)
!
! Created on 07-jun-2011 :: R. Kurosawa

       subroutine dump_vtk_3d(ascii,filename,n1,n2,n3,x,y,z,
     & rho, temp, vx, vy, vz, hx, hy, hz)

!    !
      implicit NONE
!    !
      logical, intent(in) :: ascii ! if T writes ascii version
      character(LEN=35),intent(in) :: filename ! name of VTK file
!     character(LEN=*), intent(in) :: filename ! name of VTK file
      integer, intent(in) :: n1 ! first dimension data arrays
      integer, intent(in) :: n2 ! second dimension of data arrays
      integer, intent(in) :: n3 ! thrid dimension of data arrays
!    ! Coordinates
      real*4, intent(in) :: x(n1,n2,n3) ! coordinates
      real*4, intent(in) :: y(n1,n2,n3) ! coordinates
      real*4, intent(in) :: z(n1,n2,n3) ! coordinates
!    !
      real*4, intent(in) :: rho(n1,n2,n3)    ! Density
      real*4, intent(in) :: temp(n1,n2,n3)   ! temperature
      real*4, intent(in) :: Vx(n1,n2,n3)     ! 1st component of velocity
      real*4, intent(in) :: Vy(n1,n2,n3)     ! 2nd component of velocity
      real*4, intent(in) :: Vz(n1,n2,n3)     ! 3rd component of velocity
      real*4, intent(in) :: Hx(n1,n2,n3)     ! 1st component of magnetic field vec
      real*4, intent(in) :: Hy(n1,n2,n3)     ! 2nd component of magnetic field vec
      real*4, intent(in) :: Hz(n1,n2,n3)     ! 3rd component of magnetic field vec

!    !
!    !
      integer :: i, j, k
      integer, parameter :: LUOUT=61
      character(LEN=1), parameter :: newline = ACHAR(10) ! newline symbol!
      real*8 :: ri
      real*4 :: dum_r
!    !
!    ! Now the data must be cartesian coordinates!
!    ! The gird points does not need to be regular.
!    !


!    ! Opening the outputfile
      open(unit=LUOUT, file=filename, status="replace",
     &     form="formatted")


!    ! There are five basic parts to the VTK file format.  
!    ! 1. Write file version and identifier (Header)
      write(LUOUT, "(a)") "# vtk DataFile Version 2.0"

!    ! 2. Tilte 
      write(LUOUT, "(a)") "RPIC Simulation Result"


!    ! 3. Data type (ASCII or BINARY)
      if (ascii) then
        write(LUOUT, "(a)") "ASCII"
      else
        write(LUOUT, "(a)") "BINARY"
      end if

!    ! 4.  Dataset structure 
10    format(a11, 3(1x, i6))
11    format(a, 2x, i9, 2x, a5)
      write(LUOUT, "(a)") "DATASET STRUCTURED_GRID"
      write(LUOUT, 10)    "DIMENSIONS ", n1, n2, n3
      write(LUOUT, 11)    "POINTS", n1*n2*n3, "float"

       write(6,*) 'after header'

!    ! -- now dumps coordnates
      if (.not. ascii) then
        close(LUOUT)
!       ! reopen it as unformatted and append
        open(unit=LUOUT, file=filename, status="old",
     &   form="unformatted", position="append")
c    &   form="unformatted", position="append", access="stream")  
      end if
!    ! --- now dumps data
      do k = 1, n3
        do j = 1, n2
           do i = 1, n1
              if (ascii) then
                 write(LUOUT,*) x(i,j,k),y(i,j,k),z(i,j,k)
              else
                 write(LUOUT) x(i,j,k), y(i,j,k), z(i,j,k)
              end if
           end do
        end do
      end do

       write(6,*) 'after coordinates'
!    ! 5 . Data 
      if (.not. ascii) then
        close(LUOUT)
!       ! reopen it as formatted and append
        open(unit=LUOUT, file=filename, status="old",
     & form="formatted", position="append")
      end if

12    format(a, a11, 1x, i12)
!    ! You need to add newline character for safty after writing data in binary
!     write(LUOUT, "(a)") 
       write(LUOUT, 12) newline, "POINT_DATA ", n1*n2*n3

       write(6,*) 'after new line'
!    !
!    ! Write density
!    !
c     write(LUOUT, "(a)") "SCALARS density float"
      write(LUOUT, "(a)") "SCALARS total_Bfield float"
      write(LUOUT, "(a)") "LOOKUP_TABLE default"

!    ! --- now dumps data
      if (.not. ascii) then
        close(LUOUT)
!       ! reopen it as unformatted and append
        open(unit=LUOUT, file=filename, status="old",
     &   form="unformatted", position="append")
c    &   form="unformatted", position="append", access="stream")   
      end if
      if (ascii) then
        do k = 1, n3
           do j = 1, n2
              do i = 1, n1
                 dum_r = MAX(rho(i,j,k), 1e-30) ! for safety
                 write(LUOUT,*) dum_r
              end do
           end do
        end do
      else
        do k =1, n3
           do j = 1, n2
              do i = 1, n1
                 dum_r = MAX(rho(i,j,k), 1e-30) ! for safety
                 write(LUOUT) dum_r
              end do
           end do
        end do
      end if


       write(6,*) 'after density'

!    !
!    ! Write Jx current
!    !
c     write(LUOUT, "(a)") "SCALARS density float"
      write(LUOUT, "(a)") "SCALARS jx_current float"
      write(LUOUT, "(a)") "LOOKUP_TABLE default"

!    ! --- now dumps data
      if (.not. ascii) then
        close(LUOUT)
!       ! reopen it as unformatted and append
        open(unit=LUOUT, file=filename, status="old",
     &   form="unformatted", position="append")
c    &   form="unformatted", position="append", access="stream")
      end if
      if (ascii) then
        do k = 1, n3
           do j = 1, n2
              do i = 1, n1
c                dum_r = MAX(Vx(i,j,k), 1e-30) ! for safety
c                write(LUOUT,*) dum_r
                  write(LUOUT,*) Vx(i,j,k)
              end do
           end do
        end do
      else
        do k =1, n3
           do j = 1, n2
              do i = 1, n1
c                dum_r = MAX(Vx(i,j,k), 1e-30) ! for safety
c                write(LUOUT) dum_r
                  write(LUOUT,*) Vx(i,j,k)
              end do
           end do
        end do
      end if


       write(6,*) 'after jx_current'


!    !
!    ! Write temperature
!    !
      if (.not. ascii) then
        close(LUOUT)
!       ! reopen it as formatted and append
        open(unit=LUOUT, file=filename, status="old",
     &  form="formatted", position="append")
      end if

c     write(LUOUT, "(a)") "SCALARS temperature float"
      write(LUOUT, "(a)") "SCALARS total_Efield float"
      write(LUOUT, "(a)") "LOOKUP_TABLE default"

!    ! --- now dumps data
      if (.not. ascii) then
        close(LUOUT)
!       ! reopen it as unformatted and append
        open(unit=LUOUT, file=filename, status="old",
     &       form="unformatted", position="append")
c    &       form="unformatted", position="append", access="stream")   
      end if
      if (ascii) then
        do k = 1, n3
           do j = 1, n2
              do i = 1, n1
                 dum_r = MAX(temp(i,j,k), 1e-30) ! for safety
                 write(LUOUT,*) dum_r
              end do
           end do
        end do
      else
        do k =1, n3
           do j = 1, n2
              do i = 1, n1
                 dum_r = MAX(temp(i,j,k), 1e-30) ! for safety
                 write(LUOUT) dum_r
              end do
           end do
        end do
      end if


       write(6,*) 'after temprature'




! Write veleeity
!
       if (.not. ascii) then
         close(LUOUT)
! reopen it as formatted and append
         open(unit=LUOUT, file=filename, status="old",
     & form="formatted", position="append")
       end if

c      write(LUOUT, "(a, a)") newline, "VECTORS velocity float"
       write(LUOUT, "(a, a)") newline, "VECTORS xcurrent float"

       if (.not. ascii) then
         close(LUOUT)
         ! reopen it as unformatted and append
         open(unit=LUOUT, file=filename, status="old",
     &    form="unformatted", position="append")
c    &    form="unformatted", position="append", access="stream")   
       end if
! --- now dumps data
       if (ascii) then
         do k = 1, n3
            do j = 1, n2
               do i = 1, n1
                  write(LUOUT,*) Vx(i,j,k), Vy(i,j,k), Vz(i,j,k)
               end do
            end do
         end do
       else
         do k = 1, n3
            do j = 1, n2
               do i = 1, n1
                  write(LUOUT) Vx(i,j,k), Vy(i,j,k), Vz(i,j,k)
               end do
            end do
         end do
       end if

  
       write(6,*) 'after velocity'

!
! Write B field vectees
!
       if (.not. ascii) then
         close(LUOUT)
         ! reopen it as formatted and append
         open(unit=LUOUT, file=filename, status="old",
     &  form="formatted", position="append")
       end if

       write(LUOUT, "(a, a)") newline, "VECTORS B-field float"

       if (.not. ascii) then
         close(LUOUT)
         ! reopen it as unformatted and append
         open(unit=LUOUT, file=filename, status="old",
     &   form="unformatted", position="append")
c    &   form="unformatted", position="append", access="stream")   
       end if
      ! --- now dumps data
       if (ascii) then
         do k = 1, n3
            do j = 1, n2
               do i = 1, n1
                  write(LUOUT,*) Hx(i,j,k), Hy(i,j,k), Hz(i,j,k)
               end do
            end do
         end do
       else
         do k = 1, n3
            do j = 1, n2
               do i = 1, n1
                  write(LUOUT) Hx(i,j,k), Hy(i,j,k), Hz(i,j,k)
               end do
            end do
         end do
       end if

       write(6,*) 'end of the data'

       close(LUOUT)
       return 
       end subroutine dump_vtk_3d
