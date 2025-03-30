      program plot

      parameter (Nproc=81) 
      parameter (Npx=1,Npy=9,Npz=9)

      parameter (mx=645,my=131,mz=131)
      
      parameter (nFx=(mx-5)/Npx,mFx=nFx+5)
      parameter (nFy=(my-5)/Npy,mFy=nFy+5)
      parameter (nFz=(mz-5)/Npz,mFz=nFz+5)

c     parameter (nptl=3000000)
c     parameter (nptl=2500000)
c     parameter (mb=nptl,mj=2500000)
c     parameter (mb=nptl,mj=1700000)

c     parameter (nptl=4600000)
c     parameter (mb=nptl,mj=3500000)

c     parameter (nptl=5800000)
c     parameter (mb=nptl,mj=5500000)
cJsm
      parameter (nptl=2000000)
      parameter (mb=nptl,mj=1800000)



c global particle array size: depends on number of particles in a single 
c processor domain 
      parameter (mg = mb*Nproc,mjg = mj*Nproc)
c     integer*8,parameter:: mg=mb*Nproc
c     integer*8,parameter:: mjg=mj*Nproc

            
      parameter (imax=32767)
      
      integer coords(0:(Nproc-1),3)
      
c  temporary arrays in 16-bit integers for data dump      
      integer*2 irho,ifx,ify,ifz,ixx,iyy,izz
      
      dimension ifx(mFx,mFy,mFz),ify(mFx,mFy,mFz),ifz(mFx,mFy,mFz)
      dimension irho(mFx,mFy,mFz)    
      dimension ixx(mb),iyy(mb),izz(mb)  
      integer FBDRx,FBDRy,FBDRz
      integer FBDLx,FBDLy,FBDLz
      integer FBDRxe,FBDLxe
      integer FBDRxp,FBDLxp
      integer FBD_BRx,FBD_BRy,FBD_BRz 
      integer FBD_BLx,FBD_BLy,FBD_BLz
      integer FBD_ERx,FBD_ERy,FBD_ERz
      integer FBD_ELx,FBD_ELy,FBD_ELz

      integer ions,lecs,ionj,lecj,njskip
      
      double precision bx,by,bz,ex,ey,ez
      double precision flxi,flyi,flzi,rhoi
      double precision flxe,flye,flze,rhoe
      double precision flxij,flyij,flzij,rhoij
      double precision flxej,flyej,flzej,rhoej
      real*8 me,mi
c new ken
      double precision PBLeft,PBRght,PBFrnt,PBRear,PBBot,PBTop,
     1     umax,umin,vmax,vmin,wmax,wmin,
c    1 c,DT,qi,qe,mi,me,qmi,qme,vithml,vethml,vijet,vejet,
     1 c,DT,qi,qe,qmi,qme,vithml,vethml,vijet,vejet,
     &         vithmj,vethmj,refli,refle,rselect,xj0,yj0,zj0,b0x,
c new Jacek
     &          b0y,e0z,
     &         dlxj,dlyj,dlzj,
c     &         mc,mrl,mrh,mc2,mrh2,mcol,mrow,isis,
     &         GBLeft,GBRght,DHDx,DHDy,DHDz,
c    &         FBD_BLx,FBD_BRx,
c    &         FBD_BLy,FBD_BRy,FBD_BLz,FBD_BRz,FBD_ELx,FBD_ERx,
c    &         FBD_ELy,FBD_ERy,FBD_ELz,FBD_ERz,FBDLx,FBDLy,FBDLz,
c    &         FBDRx,FBDRy,FBDRz,FBDLxe,FBDRxe,FBDLxp,FBDRxp,
     &         PVLeft,PVRght,sm1,sm2,sm3
c new Ken   
      real*4      
     &   c4,DT4,qi4,qe4,qmi4,qme4,vithml4,vethml4,vijet4,vejet4,
     &   vithmj4,vethmj4,refli4,refle4,rselect4,xj04,yj04,zj04,
     &   b0x4,b0y4,e0z4,
     &         dlxj4,dlyj4,dlzj4

      dimension sm1(-1:1,-1:1,-1:1),sm2(-1:1,-1:1,-1:1) 
      dimension sm3(-1:1,-1:1,-1:1) 
      integer   mc,mrl,mrh,mc2,mrh2,mcol,mrow,isis,
     &  nsmooth,nfilt,nstep
      dimension nfilt(16)

      
      dimension bx(mFx,mFy,mFz),by(mFx,mFy,mFz),bz(mFx,mFy,mFz)
      dimension ex(mFx,mFy,mFz),ey(mFx,mFy,mFz),ez(mFx,mFy,mFz)
      
      dimension flxi(mFx,mFy,mFz),flyi(mFx,mFy,mFz),flzi(mFx,mFy,mFz) 
      dimension rhoi(mFx,mFy,mFz)
      dimension flxe(mFx,mFy,mFz),flye(mFx,mFy,mFz),flze(mFx,mFy,mFz) 
      dimension rhoe(mFx,mFy,mFz)
      dimension flxij(mFx,mFy,mFz),flyij(mFx,mFy,mFz),flzij(mFx,mFy,mFz) 
      dimension rhoij(mFx,mFy,mFz)
      dimension flxej(mFx,mFy,mFz),flyej(mFx,mFy,mFz),flzej(mFx,mFy,mFz) 
      dimension rhoej(mFx,mFy,mFz)

c     double precision xi,yi,zi,ui,vi,wi,xe,ye,ze,ue,ve,we
c     double precision xij,yij,zij,uij,vij,wij,xej,yej,zej
c  ambient ions and electrons positions and velocities arrays       
      real*8 xi(mb),yi(mb),zi(mb),ui(mb),vi(mb),wi(mb)
      real*8 xe(mb),ye(mb),ze(mb),ue(mb),ve(mb),we(mb)

c  jet ions and electrons positions and velocities arrays      
      real*8 xij(mj),yij(mj),zij(mj),uij(mj),vij(mj),wij(mj)
      real*8 xej(mj),yej(mj),zej(mj),uej(mj),vej(mj),wej(mj)

c  global field arrays  
      real*4 bxx,byy,bzz,exx,eyy,ezz     
      dimension bxx(mx,my,mz),byy(mx,my,mz),bzz(mx,my,mz)
      dimension exx(mx,my,mz),eyy(mx,my,mz),ezz(mx,my,mz) 

c  global particle arrays
      real*4 xig,yig,zig,uig,vig,wig
      real*4 xeg,yeg,zeg,ueg,veg,weg
      real*4 xijg,yijg,zijg,uijg,vijg,wijg
      real*4 xejg,yejg,zejg,uejg,vejg,wejg
      dimension xig(mg),yig(mg),zig(mg),uig(mg),vig(mg),wig(mg)
      dimension xeg(mg),yeg(mg),zeg(mg),ueg(mg),veg(mg),weg(mg)
      dimension xijg(mjg),yijg(mjg),zijg(mjg)
      dimension uijg(mjg),vijg(mjg),wijg(mjg)
      dimension xejg(mjg),yejg(mjg),zejg(mjg)
      dimension uejg(mjg),vejg(mjg),wejg(mjg)
            
      character num*3,num1,num2,num3
      character step*6,step1,step2,step3,step4,step5
      character st01,st02,st03,st0*3,st*4,hyph
      character dchar1*6,dchar2*4,file1*10,file2*8
      character strpd*6
      character message*6
      
      double precision bxmax,bxmin,bymax,bymin,bzmax,bzmin
      double precision exmax,exmin,eymax,eymin,ezmax,ezmin
      double precision rhomax,rhomin
      double precision fxmax,fxmin,fymax,fymin,fzmax,fzmin
c     real*8 mi,me

c assign coordinates to process ids       
       open(3,file="topology")
       do i = 0,Nproc-1
        read(3,*) ii,coords(i,1),coords(i,2),coords(i,3)      
       end do
       close(3)

c specify "nst" - which data need to be processed
c      nst = 7
c      open(unit=2, file='nstn', status='old')

c      read(2,*) nst

       c = 1.0

       do in=1,25
       
       nst = in

        if(nst.eq.3) go to 1578
        if(nst.eq.5) go to 1578
        if(nst.eq.7) go to 1578
        if(nst.eq.9) go to 1578
        if(nst.eq.10) go to 1578
        if(nst.eq.11) go to 1578
        if(nst.eq.13) go to 1578
        if(nst.eq.15) go to 1578
        if(nst.eq.17) go to 1578
        if(nst.eq.19) go to 1578
        if(nst.eq.20) go to 1578
        if(nst.eq.21) go to 1578

c        go to 2578

 1578  continue
       strpd = 'partd_'
       hyph='_'

       st01 = char(int(nst/100.)+48)
       st02 = char(int((nst-int(nst/100.)*100)/10.)+48)
       st03 = char(int(nst-int(nst/10.)*10)+48)	 
       st0 = st01//st02//st03
       st = hyph//st0     

       do k = 1,mz
        do j = 1,my
	 do i = 1,mx
	  bxx(i,j,k)=0.0  
	  byy(i,j,k)=0.0
	  bzz(i,j,k)=0.0
          exx(i,j,k)=0.0
	  eyy(i,j,k)=0.0
	  ezz(i,j,k)=0.0
         end do
	end do
       end do

       iong = 0
       lecg = 0
       ionjg = 0
       lecjg = 0
  
c  first calculate quantities and currents at ALL grid points, then trimm edges       
       do np = 0,Nproc-1       
        print *, np
       
        num1 = char(int(np/100.)+48)
        num2 = char(int((np-int(np/100.)*100)/10.)+48)
        num3 = char(int(np-int(np/10.)*10)+48)
        num = num1//num2//num3

        open(7,file=strpd//num//st,form='unformatted')

        write(*,*) 'before read 7'
	
        read(7) c
        read(7) ions,lecs,ionj,lecj
        write(*,*) 'c, ions,lecs,ionj,lecj =', 
     1   c, ions,lecs,ionj,lecj
        read(7) PBLeft,PBRght,PBFrnt,PBRear,PBBot,PBTop
        read(7) bxmax,bxmin,bymax,bymin,bzmax,bzmin
        read(7) ifx,ify,ifz
        call F_reconv(bx,by,bz,ifx,ify,ifz,mFx,mFy,mFz,imax,
     &                  bxmax,bxmin,bymax,bymin,bzmax,bzmin)
        read(7) exmax,exmin,eymax,eymin,ezmax,ezmin
        read(7) ifx,ify,ifz
        call F_reconv(ex,ey,ez,ifx,ify,ifz,mFx,mFy,mFz,imax,
     &                  exmax,exmin,eymax,eymin,ezmax,ezmin)
c  read ambient ions       
       read(7) (ixx(i),i=1,ions),(iyy(i),i=1,ions),(izz(i),i=1,ions)
       call X_reconv(ions,xi,yi,zi,ixx,iyy,izz,mb,mb,imax,
     &            PBLeft,PBRght,PBFrnt,PBRear,PBBot,PBTop)
       read(7) umax,umin,vmax,vmin,wmax,wmin
       read(7) (ixx(i),i=1,ions),(iyy(i),i=1,ions),(izz(i),i=1,ions)
c      call V_reconv(ions,ui,vi,wi,ixx,iyy,izz,mb,mb,imax,
c    &                     umax,umin,vmax,vmin,wmax,wmin)
cNewJacek changes V_convert into P_convert everywhere below
       call P_reconv(ions,ui,vi,wi,ixx,iyy,izz,mb,mb,imax,
     &                      umax,umin,vmax,vmin,wmax,wmin,c)
cPconv
c  read ambient electrons       
       read(7) (ixx(i),i=1,lecs),(iyy(i),i=1,lecs),(izz(i),i=1,lecs)
       call X_reconv(lecs,xe,ye,ze,ixx,iyy,izz,mb,mb,imax,
     &            PBLeft,PBRght,PBFrnt,PBRear,PBBot,PBTop)       
       read(7) umax,umin,vmax,vmin,wmax,wmin
       read(7) (ixx(i),i=1,lecs),(iyy(i),i=1,lecs),(izz(i),i=1,lecs)
c      call V_reconv(lecs,ue,ve,we,ixx,iyy,izz,mb,mb,imax,
c    &                      umax,umin,vmax,vmin,wmax,wmin)       
       call P_reconv(lecs,ue,ve,we,ixx,iyy,izz,mb,mb,imax,
     &                      umax,umin,vmax,vmin,wmax,wmin,c)
c  read CR ions       
c      write(6,*)'ionj=',ionj
       if (ionj.gt.0) then
        read(7) (ixx(i),i=1,ionj),(iyy(i),i=1,ionj),(izz(i),i=1,ionj)
        call X_reconv(ionj,xij,yij,zij,ixx,iyy,izz,mj,mb,imax,
     &                PBLeft,PBRght,PBFrnt,PBRear,PBBot,PBTop)
        read(7) umax,umin,vmax,vmin,wmax,wmin
        read(7) (ixx(i),i=1,ionj),(iyy(i),i=1,ionj),(izz(i),i=1,ionj)
c       call V_reconv(ionj,uij,vij,wij,ixx,iyy,izz,mj,mb,imax,
c    &                          umax,umin,vmax,vmin,wmax,wmin)
        call P_reconv(ionj,uij,vij,wij,ixx,iyy,izz,mj,mb,imax,
     &                          umax,umin,vmax,vmin,wmax,wmin,c)
       end if

c  read CR electrons
       if (lecj.gt.0) then     
        read(7) (ixx(i),i=1,lecj),(iyy(i),i=1,lecj),(izz(i),i=1,lecj)
        call X_reconv(lecj,xej,yej,zej,ixx,iyy,izz,mj,mb,imax,
     &                PBLeft,PBRght,PBFrnt,PBRear,PBBot,PBTop)       
        read(7) umax,umin,vmax,vmin,wmax,wmin
        read(7) (ixx(i),i=1,lecj),(iyy(i),i=1,lecj),(izz(i),i=1,lecj)
c       call V_reconv(lecj,uej,vej,wej,ixx,iyy,izz,mj,mb,imax,
c    &                          umax,umin,vmax,vmin,wmax,wmin)
        call P_reconv(lecj,uej,vej,wej,ixx,iyy,izz,mj,mb,imax,
     &                          umax,umin,vmax,vmin,wmax,wmin,c)
       end if
       
       read(7) c,DT,qi,qe,mi,me,qmi,qme,vithml,vethml,vijet,vejet,
     &         vithmj,vethmj,refli,refle,rselect,xj0,yj0,zj0,b0x,
c new Jacek
     &          b0y,e0z,
     &         dlxj,dlyj,dlzj,lyj1,lzj1,
     &         mc,mrl,mrh,mc2,mrh2,mcol,mrow,isis    
c NewKen
c    &          ,njskip,zz1,zz2
     &          ,njskip

       
       read(7) GBLeft,GBRght,DHDx,DHDy,DHDz,FBD_BLx,FBD_BRx,
     &         FBD_BLy,FBD_BRy,FBD_BLz,FBD_BRz,FBD_ELx,FBD_ERx,
     &         FBD_ELy,FBD_ERy,FBD_ELz,FBD_ERz,FBDLx,FBDLy,FBDLz,
     &         FBDRx,FBDRy,FBDRz,FBDLxe,FBDRxe,FBDLxp,FBDRxp,
     &         PVLeft,PVRght,nsmooth,sm1,sm2,sm3,nfilt,nstep 
       
       
       close(7)

       i0=coords(np,1)*nFx
       j0=coords(np,2)*nFy
       k0=coords(np,3)*nFz 

       ibmin=3
       ibmax=nFx+2
       jbmin=3
       jbmax=nFy+2
       kbmin=3
       kbmax=nFz+2
 
       if (coords(np,1).eq.0) ibmin=1
       if (coords(np,1).eq.(Npx-1))ibmax=mFx
       if (coords(np,2).eq.0) jbmin=1
       if (coords(np,2).eq.(Npy-1))jbmax=mFy
       if (coords(np,3).eq.0) kbmin=1
       if (coords(np,3).eq.(Npz-1))kbmax=mFz
       
       do kb = kbmin,kbmax
        do jb = jbmin,jbmax
	 do ib = ibmin,ibmax
   	  i = ib+i0
          j = jb+j0
	  k = kb+k0
	  bxx(i,j,k)=bx(ib,jb,kb)/c  
	  byy(i,j,k)=by(ib,jb,kb)/c
	  bzz(i,j,k)=bz(ib,jb,kb)/c 
          exx(i,j,k)=ex(ib,jb,kb)
	  eyy(i,j,k)=ey(ib,jb,kb)
	  ezz(i,j,k)=ez(ib,jb,kb)
         end do
	end do
       end do

       write(6,*)'bxx(300,35,35),byy(300,35,35),bzz(300,35,35)=',
     1 bxx(300,35,35),byy(300,35,35),bzz(300,35,35)

       write(6,*)'exx(300,35,35),eyy(300,35,35),ezz(300,35,35)=',
     1 exx(300,35,35),eyy(300,35,35),ezz(300,35,35)

       
       do i = 1,ions
        iong = iong+1
        xig(iong)=xi(i)
        yig(iong)=yi(i)
        zig(iong)=zi(i)
        uig(iong)=ui(i)
        vig(iong)=vi(i)
        wig(iong)=wi(i)
       end do

       do i = 1,lecs
        lecg = lecg+1
        xeg(lecg)=xe(i)
        yeg(lecg)=ye(i)
        zeg(lecg)=ze(i)
        ueg(lecg)=ue(i)
        veg(lecg)=ve(i)
        weg(lecg)=we(i)
       end do

       do i = 1,ionj
        ionjg = ionjg+1
        xijg(ionjg)=xij(i)
        yijg(ionjg)=yij(i)
        zijg(ionjg)=zij(i)
        uijg(ionjg)=uij(i)
        vijg(ionjg)=vij(i)
        wijg(ionjg)=wij(i)
       end do

       do i = 1,lecj
        lecjg = lecjg+1
        xejg(lecjg)=xej(i)
        yejg(lecjg)=yej(i)
        zejg(lecjg)=zej(i)
        uejg(lecjg)=uej(i)
        vejg(lecjg)=vej(i)
        wejg(lecjg)=wej(i)
       end do

      end do
      write(6,*)'xejg(7),uejg(7),weg(7),weg(7),ueg(7)=',
     1   xejg(7),uejg(7),weg(7),weg(7),ueg(7)
                         
c  write to file 
        dchar1='bfield'
        file1= dchar1//st

        open(15,file=file1,form='unformatted')
        write(15) bxx,byy,bzz
        close(15)

        dchar1='efield'
        file1= dchar1//st

        open(16,file=file1,form='unformatted')
        write(16) exx,eyy,ezz
        close(16)

        dchar2='ions'
        file2= dchar2//st

cheung
        open(17,file=file2,form='unformatted')
        write(17) iong
c       write(17) xig(1:iong),yig(1:iong),zig(1:iong)
c       write(17) uig(1:iong),vig(1:iong),wig(1:iong)
        write(17) xig(1:iong)
        write(17) yig(1:iong)
        write(17) zig(1:iong)
        write(17) uig(1:iong)
        write(17) vig(1:iong)
        write(17) wig(1:iong)
        close(17)

        dchar2='lecs'
        file2= dchar2//st

cheung
        open(18,file=file2,form='unformatted')
        write(18) lecg
c       write(18) xeg(1:lecg),yeg(1:lecg),zeg(1:lecg)
c       write(18) ueg(1:lecg),veg(1:lecg),weg(1:lecg)
        write(18) xeg(1:lecg)
        write(18) yeg(1:lecg)
        write(18) zeg(1:lecg)
        write(18) ueg(1:lecg)
        write(18) veg(1:lecg)
        write(18) weg(1:lecg)
        close(18)

        dchar2='ionj'
        file2= dchar2//st

cheung
        open(19,file=file2,form='unformatted')
        write(19) ionjg
c       write(19) xijg(1:ionjg),yijg(1:ionjg),zijg(1:ionjg)
c       write(19) uijg(1:ionjg),vijg(1:ionjg),wijg(1:ionjg)
        write(19) xijg(1:ionjg)
        write(19) yijg(1:ionjg)
        write(19) zijg(1:ionjg)
        write(19) uijg(1:ionjg)
        write(19) vijg(1:ionjg)
        write(19) wijg(1:ionjg)
        close(19)

        dchar2='lecj'
        file2= dchar2//st

cheung
        open(20,file=file2,form='unformatted')
        write(20) lecjg
c       write(20) xejg(1:lecjg),yejg(1:lecjg),zejg(1:lecjg)
c       write(20) uejg(1:lecjg),vejg(1:lecjg),wejg(1:lecjg)
        write(20) xejg(1:lecjg)
        write(20) yejg(1:lecjg)
        write(20) zejg(1:lecjg)
        write(20) uejg(1:lecjg)
        write(20) vejg(1:lecjg)
        write(20) wejg(1:lecjg)
        close(20)

c    &   c4,DT4,qi4,qe4,qmi4,qme4,vithml4,vethml4,vijet4,vejet4,
c    &   vithmj4,vethmj4,refli4,refle4,rselect4,xj04,yj04,zj04,
c    &   b0x4,b0y4,e0z4,
c    &         dlxj4,dlyj4,dlzj4
        c4    = c
        DT4   = DT
        qi4   = qi
        qe4   = qe
        qmi4  = qmi
        qme4  = qme
        vithml4 = vithml
        vethml4 = vethml
        vijet4  = vijet
        vejet4  = vejet 
        vithmj4 = vithmj
        vethmj4 = vethmj
        refli4  = refli
        refle4  = refle
       rselect4 = rselect
        xj04    = xj0
        yj04    = yj0
        zj04    = zj0
        b0x4    = b0x
        b0y4    = b0y
        e0z4    = e0z
c    &         dlxj4,dlyj4,dlzj4

        dchar1='kenout'
        file1= dchar1//st

cheung
c       open(20,file=file2,form='unformatted')

c       message='kenout'

        open(21,file=file1,form='unformatted')
c       write(21) nstep,c,DT,qi,qe,mi,me,qmi,qme,vithml,vethml,vijet,vejet,
c    &         vithmj,vethmj,refli,refle,rselect,xj0,yj0,zj0,b0x,
c new Jacek
c    &          b0y,e0z 
c New Ken
        write(21) nstep,
     &   c4,DT4,qi4,qe4,qmi4,qme4,vithml4,vethml4,vijet4,vejet4,
     &   vithmj4,vethmj4,refli4,refle4,rselect4,xj04,yj04,zj04,
     &   b0x4,b0y4,e0z4

c,DT,qi,qe,mi,me,qmi,qme,vithml,vethml,vijet,vejet,
c    &         vithmj,vethmj,refli,refle,rselect,xj0,yj0,zj0,b0x,
c new Jacek
c    &          b0y,e0z 
        close(21)
 2578  continue
       enddo

       end
c *****************************************************************      
      subroutine F_reconv(fx,fy,fz,ifx,ify,ifz,mFx,mFy,mFz,imax,
     &                      fxmax,fxmin,fymax,fymin,fzmax,fzmin)    
            
      double precision fxmax,fxmin,fymax,fymin,fzmax,fzmin
      double precision fx,fy,fz
      integer*2 ifx,ify,ifz
      integer   mFx,mFy,mFz
      dimension fx(mFx,mFy,mFz),fy(mFx,mFy,mFz),fz(mFx,mFy,mFz)
      dimension ifx(mFx,mFy,mFz),ify(mFx,mFy,mFz),ifz(mFx,mFy,mFz)      
      
      do k = 1,mFz
       do j = 1,mFy
	do i = 1,mFx	 
	 fx(i,j,k)=dfloat(ifx(i,j,k))*(fxmax-fxmin)/dfloat(imax) + fxmin
	 fy(i,j,k)=dfloat(ify(i,j,k))*(fymax-fymin)/dfloat(imax) + fymin
	 fz(i,j,k)=dfloat(ifz(i,j,k))*(fzmax-fzmin)/dfloat(imax) + fzmin
        end do
       end do
      end do
           
      return
      end       
c *****************************************************************      
      subroutine X_reconv(ipar,x,y,z,ix,iy,iz,mh,mb,imax,
     &                 PBLeft,PBRght,PBFrnt,PBRear,PBBot,PBTop)    
      
      integer*2 ix,iy,iz
      double precision PBLeft,PBRght,PBFrnt,PBRear,PBBot,PBTop
c     dimension x(mh),y(mh),z(mh)
      real*8 x(mh),y(mh),z(mh)
      dimension ix(mb),iy(mb),iz(mb)
     
      do i = 1,ipar
       x(i) = dfloat(ix(i))*(PBRght-PBLeft)/dfloat(imax) + PBLeft
       y(i) = dfloat(iy(i))*(PBRear-PBFrnt)/dfloat(imax) + PBFrnt
       z(i) = dfloat(iz(i))*(PBTop-PBBot)/dfloat(imax) + PBBot 
      end do
                       
      return
      end                  
c *****************************************************************      
      subroutine V_reconv(ipar,u,v,w,iu,iv,iw,mh,mb,imax,
     &                     umax,umin,vmax,vmin,wmax,wmin)    
      
      integer*2 iu,iv,iw
      double precision  umax,umin,vmax,vmin,wmax,wmin
c     dimension u(mh),v(mh),w(mh)
      real*8 u(mh),v(mh),w(mh)
      dimension iu(mb),iv(mb),iw(mb)
                              
      do i = 1,ipar
       u(i) = dfloat(iu(i))*(umax-umin)/dfloat(imax) + umin
       v(i) = dfloat(iv(i))*(vmax-vmin)/dfloat(imax) + vmin
       w(i) = dfloat(iw(i))*(wmax-wmin)/dfloat(imax) + wmin
      end do
                       
      return
      end      
c *****************************************************************
      subroutine P_convert(ipar,u,v,w,iu,iv,iw,mh,mb,imax,
     &                    umax,umin,vmax,vmin,wmax,wmin,c)

      integer*2 iu,iv,iw
     
      double precision  umax,umin,vmax,vmin,wmax,wmin,c
      real*8 u(mh),v(mh),w(mh)
      real*8 ptemp(mh)
      dimension iu(mb),iv(mb),iw(mb)
      
      umax=-1000.0
      umin= 1000.0
      vmax=-1000.0
      vmin= 1000.0
      wmax=-1000.0
      wmin= 1000.0

      do i = 1,ipar
       gamma=c/dsqrt(c**2-u(i)**2-v(i)**2-w(i)**2)
       ptemp(i)=u(i)*gamma

       umax=max(umax,ptemp(i))
       umin=min(umin,ptemp(i))
      end do

      do i = 1,ipar
       iu(i) = imax*(ptemp(i)-umin)/(umax-umin)
      end do

      do i = 1,ipar
       gamma=c/dsqrt(c**2-u(i)**2-v(i)**2-w(i)**2)
       ptemp(i)=v(i)*gamma

       vmax=max(vmax,ptemp(i))
       vmin=min(vmin,ptemp(i))
      end do

      do i = 1,ipar
       iv(i) = imax*(ptemp(i)-vmin)/(vmax-vmin)
      end do

      do i = 1,ipar
       gamma=c/dsqrt(c**2-u(i)**2-v(i)**2-w(i)**2)
       ptemp(i)=w(i)*gamma
      
       wmax=max(wmax,ptemp(i))
       wmin=min(wmin,ptemp(i))
      end do 
      
      do i = 1,ipar
       iw(i) = imax*(ptemp(i)-wmin)/(wmax-wmin)
      end do
       
      return
      end
       
c *****************************************************************
      subroutine P_reconv(ipar,u,v,w,iu,iv,iw,mh,mb,imax,
     &                   umax,umin,vmax,vmin,wmax,wmin,c)
       
      integer*2 iu,iv,iw
      double precision  umax,umin,vmax,vmin,wmax,wmin,c

      real*8 u(mh),v(mh),w(mh)
      dimension iu(mb),iv(mb),iw(mb)

      do i = 1,ipar
       u(i) = dfloat(iu(i))*(umax-umin)/dfloat(imax) + umin
       v(i) = dfloat(iv(i))*(vmax-vmin)/dfloat(imax) + vmin
       w(i) = dfloat(iw(i))*(wmax-wmin)/dfloat(imax) + wmin

       gammac=dsqrt(u(i)*u(i)+v(i)*v(i)+w(i)*w(i)+c**2)

       u(i)=u(i)*c/gammac
       v(i)=v(i)*c/gammac
       w(i)=w(i)*c/gammac
      end do

      return
      end

