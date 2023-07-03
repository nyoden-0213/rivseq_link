      program calc_rivseq_link
! ================================================
!* PURPOSE: Calculate River Sequence considering Canal Link Connectivity
!
! (C) Naho Yoden & Dai Yamazaki (UTokyo)  July 2023
!
! Licensed under the Apache License, Version 2.0 (the "License");
!   You may not use this file except in compliance with the License.
!   You may obtain a copy of the License at: http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software distributed under the License is 
!  distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. 
! See the License for the specific language governing permissions and limitations under the License.
! ================================================
      implicit none
! ==============================
! calculation type
      character*256       ::  buf
      character*256       ::  type                    !! 'bin' for binary, 'cdf' for netCDF
      data                    type /'bin'/
! input files
      character*256       ::  diminfo                 !! dimention info file
      data                    diminfo /'testdata/diminfo_test-1deg.txt'/
! dimention
      integer             ::  ix, iy, jx, jy
      integer             ::  nx, ny                  !! river map dimention
      integer             ::  nxin, nyin              !! input map dimention
      integer             ::  nflp                    !! floodplain layer
! river network map
      integer,allocatable ::  nextx(:,:), nexty(:,:)  !! downstream xy
      real,allocatable    ::  lon(:), lat(:)          !! longitude, latidude [deg]
      real                ::  west,east,north,south   !! domain boundary [deg]
      real                ::  gsize                   !! grid size [deg]
      integer,allocatable ::  upst(:,:)             !! # of upstreams
      integer,allocatable ::  upnow(:,:)             !! count upstream
      integer,allocatable ::  xseq(:), yseq(:)       !! calculation sequence from upstream
      integer             ::  iseq, jseq, nseqpre, nseqnow, again
! input matrix
      integer             ::  inpn                    !! max number of input grid for one river grid
! variables
      real,allocatable    ::  rivseqmod(:,:)
      real,allocatable    ::  rivseqmod2(:)
! files
      character*256       ::  cnextxy        !! river network map, grid area
      integer             ::  ios
      character*256       ::  cinpmat                 !! input matrix
      character*256       ::  crivseqmod                 !! river sequence
      parameter              (crivseqmod='testdata/rivseqmod.bin')
! undef
      integer             ::  imis                !! integer undefined
! value
      real                ::  rmis                !! real    undefined
! value
      parameter              (imis = -9999)
      parameter              (rmis = 1.e+20)
! link canal
      real                ::  r1tmp(15)
      integer,allocatable ::  i2llnk_odx(:,:,:), i2llnk_ody(:,:,:), i2llnk_id(:,:,:)
      integer             ::  i,j,i0rec
      integer             ::  kx, ky
      integer             ::  m0rec
      parameter              (m0rec=5)
      character*256       ::  ccanalinfo              !!canal info
      parameter              (ccanalinfo='testdata/canal_info.csv')
      character*256       ::  dummy
! ================================================
      print *, 'CALC_OUTCLM - calculate annual mean discharge from runoff climatology'
      print *, 'calc_outclm: read parameters from arguments'

      call getarg(1,buf)
       if( buf/='' ) read(buf,*) type
      call getarg(2,buf)
       if( buf/='' ) read(buf,'(a128)') diminfo
      print *, 'TYPE=',    trim(type)
      print *, 'DIMINFO=', trim(diminfo)
! ===============================
! read river network and input file dimentions
      if( type=='cdf')then
        print *, 'calculation for netCDF map'
      else
        type='bin'
        print *, 'calculation for binary map'
      endif

      open(11,file=diminfo,form='formatted')
      read(11,*) nx
      read(11,*) ny
      read(11,*) nflp
      read(11,*) nxin
      read(11,*) nyin
      read(11,*) inpn
      read(11,'(a)') cinpmat
      read(11,*) west
      read(11,*) east
      read(11,*) north
      read(11,*) south
      close(11)

      print *, trim(cinpmat)
      print *, nx, ny, nxin, nyin, inpn
      allocate(nextx(nx,ny),nexty(nx,ny))
      allocate(i2llnk_odx(m0rec,nx,ny),i2llnk_ody(m0rec,nx,ny),i2llnk_id(m0rec,nx,ny))
! ===========================================
      allocate(lon(nx),lat(ny))
      gsize=(east-west)/real(nx)
      do ix=1,nx
        lon(ix)=west+(real(ix)-0.5)*gsize
      enddo
      do iy=1,ny
        lat(iy)=north-(real(iy)-0.5)*gsize
      enddo
! =========
      print *, 'calc_outclm: read nextxy.bin'
      cnextxy='testdata/nextxy.bin'
      open(11,file=cnextxy,form='unformatted',access='direct',recl=4*nx*ny,status='old',iostat=ios)
      read(11,rec=1) nextx
      read(11,rec=2) nexty
      close(11)
! ==========
! read link canal info
      i2llnk_odx=0
      i2llnk_ody=0
      i2llnk_id=0
      open(15,file=ccanalinfo,status='old')
      read (15, '()')
      do i=0,100,1
        read (15, *, iostat=ios) dummy,dummy,(r1tmp(j),j=1,15)
        if(ios<0)exit
        if(int(r1tmp(12)).eq.-999)cycle
        do i0rec=1,m0rec
          if(i2llnk_odx(i0rec,int((r1tmp(12)-65)*12)+1,int((38-r1tmp(13))*12)+1).ne.0)cycle
          i2llnk_odx(i0rec,int((r1tmp(12)-65)*12)+1,int((38-r1tmp(13))*12)+1)=int((r1tmp(14)-65)*12)+1
          i2llnk_ody(i0rec,int((r1tmp(12)-65)*12)+1,int((38-r1tmp(13))*12)+1)=int((38-r1tmp(15))*12)+1
          i2llnk_id(i0rec,int((r1tmp(12)-65)*12)+1,int((38-r1tmp(13))*12)+1)=int(r1tmp(1))
          exit
        end do
      end do
      close (15) 
! ==========
      print *, 'calc_outclm: calculate river sequence'
      allocate(upst(nx,ny),upnow(nx,ny),xseq(nx*ny),yseq(nx*ny),rivseqmod(nx,ny),rivseqmod2(nx*ny))
      upst(:,:)=0
      upnow(:,:)=0
      xseq(:)=-9999
      yseq(:)=-9999
      rivseqmod(:,:)=0
! count number of upstreams
      do iy=1, ny
        do ix=1, nx
          if( nextx(ix,iy)>0 )then
            jx=nextx(ix,iy)
            jy=nexty(ix,iy)
            upst(jx,jy)=upst(jx,jy)+1
            do i0rec=1,m0rec
              if( i2llnk_odx(i0rec,ix,iy)>0 )then
                kx=i2llnk_odx(i0rec,ix,iy)
                ky=i2llnk_ody(i0rec,ix,iy)
                upst(kx,ky)=upst(kx,ky)+1
                print *, 'link', i2llnk_id(i0rec,ix,iy), 'is', ix, iy, 'to', kx, ky
              endif
            end do
          elseif( nextx(ix,iy)==-9999 )then
            upst(ix,iy)=-9999
          endif
        end do
      end do
      print *, 'upst',upst(103,64)
! find topmost grid, and register to xseq & yseq
      nseqpre=0
      iseq=nseqpre
      do iy=1, ny
        do ix=1, nx
          if( upst(ix,iy)==0 )then
            iseq=iseq+1
            xseq(iseq)=ix   !! sort from upstream to downstream
            yseq(iseq)=iy
            rivseqmod(ix,iy)=1
          endif
        end do
      end do
      nseqnow=iseq

! find the next downstream, register to xseq & yseq when upst=upnow
      print *, 'main calculation'
      again=1
      do while( again==1 )
        again=0
        jseq=nseqnow
        do iseq=nseqpre+1, nseqnow          
          ix=xseq(iseq)
          iy=yseq(iseq)    
          if( nextx(ix,iy)>0 )then
            jx=nextx(ix,iy)
            jy=nexty(ix,iy)
            upnow(jx,jy)=upnow(jx,jy)+1
            rivseqmod(jx,jy)=max(rivseqmod(jx,jy),rivseqmod(ix,iy)+1)
            if( upnow(jx,jy)==upst(jx,jy) )then  !! if all upstream are registered, then target downstream can be registered
              again=1
              jseq=jseq+1
              xseq(jseq)=jx
              yseq(jseq)=jy
            endif
            do i0rec=1,m0rec
              if( i2llnk_odx(i0rec,ix,iy)>0 )then
                kx=i2llnk_odx(i0rec,ix,iy)
                ky=i2llnk_ody(i0rec,ix,iy)
                upnow(kx,ky)=upnow(kx,ky)+1
                rivseqmod(kx,ky)=max(rivseqmod(kx,ky),rivseqmod(ix,iy)+1)
                print *, 'link', i2llnk_id(i0rec,ix,iy), 'is', ix, iy, 'to', kx, ky
                if( upnow(kx,ky)==upst(kx,ky) )then  !! if all upstream are registered, then target downstream can be registered
                  again=1
                  jseq=jseq+1
                  xseq(jseq)=kx
                  yseq(jseq)=ky
                endif
              endif
            end do
            
          endif
        end do
        nseqpre=nseqnow
        nseqnow=jseq
      end do

      do iy=1, ny
        do ix=1, nx
          if(nextx(ix,iy)>0)then
            if(rivseqmod(ix,iy)>=rivseqmod(nextx(ix,iy),nexty(ix,iy)))then
              print *, 'error: next grid of (', ix, ',', iy, ') has smaller river sequence'
            endif
          endif
        end do
      end do
      do iy=1, ny
        do ix=1, nx
          do i0rec=1,m0rec
            if( i2llnk_odx(i0rec,ix,iy)>0 )then
              if(rivseqmod(ix,iy)>=rivseqmod(i2llnk_odx(i0rec,ix,iy),i2llnk_ody(i0rec,ix,iy)))then
                print *, 'error: link canal (', ix, ',', iy, ') to (', i2llnk_odx(i0rec,ix,iy), ',', i2llnk_ody(i0rec,ix,iy), &
                &') have river sequence', rivseqmod(ix,iy), 'and', rivseqmod(i2llnk_odx(i0rec,ix,iy),i2llnk_ody(i0rec,ix,iy))
              endif
            endif
          end do
        end do
      end do
! write river sequence
      do iy=1, ny
        do ix=1, nx
          rivseqmod2(ix+iy*nx)=rivseqmod(ix,iy)
        end do
      end do
      open(11,file=crivseqmod,form='unformatted',access='direct',recl=4*nx*ny)
      write(11,rec=1) rivseqmod2
      close(11)
!!================================================
      end program calc_rivseq_link