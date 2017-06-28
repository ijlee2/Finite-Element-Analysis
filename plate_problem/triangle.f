!     ------------------------------------------------------------------
!
!     triangle.f
!
!
!     Created on Mon Nov 26 14:57:33 2007
!     Copyright (c) 2007 MyCompany. All rights reserved.
!
!
!     ------------------------------------------------------------------

      program triangle
      integer maxnodes,maxelem
      parameter (maxnodes=2000,maxelem=3000)
      real*8 nodedata(maxnodes,2),EE,vv
      real*8 kglob(2*maxnodes,2*maxnodes),fval(2*maxnodes)
      real*8 uval(2*maxnodes),f(2*maxnodes)
      real*8 u(maxnodes),v(maxnodes),B(maxelem,3,6)
      real*8 strn(maxelem,3),strs(maxelem,3)
      integer fnode(2*maxnodes),unode(2*maxnodes)
      integer fdof(2*maxnodes),udof(2*maxnodes)
      integer eledata(maxelem,3),gcon(maxnodes,2),nnodes,nelem
      integer nforces,ndisp,i
      
      call readinp(maxnodes,maxelem,nnodes,nelem,nforces,ndisp,
     1             nodedata,eledata,gcon,fnode,fdof,
     2             fval,unode,udof,uval,EE,vv)
     
      call stiffness(maxelem,nelem,maxnodes,eledata,nodedata,
     1               gcon,kglob,EE,vv,B)
     
      call applyf(maxnodes,nnodes,nforces,fnode,fdof,fval,gcon,f)

      call applydisp(maxnodes,nnodes,ndisp,unode,udof,uval,gcon,
     1               kglob,f)

      call gaussj(kglob,2*nnodes,2*maxnodes,f,1)

      do 10 i=1,nnodes
        u(i)=f(2*i-1)
        v(i)=f(2*i)
 10   continue
     
      call stress(maxelem,maxnodes,nelem,eledata,nodedata,u,v,
     1            strn,strs,EE,vv,B)
     
      call structure(maxelem,maxnodes,nelem,eledata,nodedata,u,v)
     
!      print*,'Nodal displacements, node#, u, v'
!      do 20 i=1,nnodes
!        print*,i,u(i),v(i)
! 20   continue
!      print*,' '
      
      stop
      end
      
      subroutine readinp(maxnodes,maxelem,nnodes,nelem,nforces,ndisp,
     1                   nodedata,eledata,gcon,fnode,fdof,
     2                   fval,unode,udof,uval,EE,vv)
      real*8 nodedata(maxnodes,2),EE,vv,fval(2*maxnodes)
      real*8 uval(2*maxnodes)
      integer eledata(maxelem,3),gcon(maxnodes,2),nnodes,nelem
      integer maxnodes,maxelem,i,num,fnode(2*maxnodes)
      integer fdof(2*maxnodes),unode(2*maxnodes),udof(2*maxnodes)
      integer nforces,ndisp
      
      open(unit=12,file='nodesR')
      open(unit=13,file='elementsR')
      open(unit=14,file='forcesR')
      open(unit=15,file='displacementsR')
     
      read(12,*) nnodes
      do 10 i=1,nnodes
        read(12,*) num,nodedata(i,1),nodedata(i,2)
        gcon(i,1)=2*i-1
        gcon(i,2)=2*i
 10   continue
 
      read(13,*) nelem,EE,vv
      do 20 i=1,nelem
        read(13,*) num,eledata(i,1),eledata(i,2),eledata(i,3)
 20   continue
 
      read(14,*) nforces
      do 30 i=1,nforces
        read(14,*) fnode(i),fdof(i),fval(i)
 30   continue
 
      read(15,*) ndisp
      do 40 i=1,ndisp
        read(15,*) unode(i),udof(i),uval(i)
 40   continue
 
      return
      end
      
      subroutine stiffness(maxelem,nelem,maxnodes,eledata,nodedata,
     1                     gcon,kglob,EE,vv,B)
      real*8 nodedata(maxnodes,2),EE,vv,x1,y1,x2,y2,x3,y3,AA
      real*8 kele(6,6),kglob(2*maxnodes,2*maxnodes),B(maxelem,3,6)
      real*8 a0,a1,a2,b0,b1,b2,c0,c1,c2,CC(3,3)
      integer i,j,maxelem,maxnodes,nelem,eledata(maxelem,3)
      integer inode,idof,jnode,jdof,gcon(maxnodes,2),nodei,nodej
      integer nodek,ii,jj,kk,ll
      
      CC(1,1)=EE/(1.d0-vv*vv)
      CC(1,2)=vv*CC(1,1)
      CC(1,3)=0.d0
      CC(2,1)=CC(1,2)
      CC(2,2)=CC(1,1)
      CC(2,3)=0.d0
      CC(3,1)=0.d0
      CC(3,2)=0.d0
      CC(3,3)=EE/2.d0/(1.d0+vv)
      
      do 10 i=1,2*maxnodes
        do 11 j=1,2*maxnodes
          kglob(i,j)=0.d0
 11     continue
 10   continue
 
      do 20 i=1,nelem
        nodei=eledata(i,1)
        nodej=eledata(i,2)
        nodek=eledata(i,3)
        x1=nodedata(nodei,1)
        x2=nodedata(nodej,1)
        x3=nodedata(nodek,1)
        y1=nodedata(nodei,2)
        y2=nodedata(nodej,2)
        y3=nodedata(nodek,2)
        AA=x1*y2-x2*y1+x2*y3-x3*y2+x3*y1-x1*y3
        a0=(x2*y3-x3*y2)/AA
        a1=(y2-y3)/AA
        a2=(x3-x2)/AA
        b0=(x3*y1-x1*y3)/AA
        b1=(y3-y1)/AA
        b2=(x1-x3)/AA
        c0=(x1*y2-x2*y1)/AA
        c1=(y1-y2)/AA
        c2=(x2-x1)/AA
        B(i,1,1)=a1
        B(i,1,2)=0.d0
        B(i,1,3)=b1
        B(i,1,4)=0.d0
        B(i,1,5)=c1
        B(i,1,6)=0.d0
        B(i,2,1)=0.d0
        B(i,2,2)=a2
        B(i,2,3)=0.d0
        B(i,2,4)=b2
        B(i,2,5)=0.d0
        B(i,2,6)=c2
        B(i,3,1)=a2
        B(i,3,2)=a1
        B(i,3,3)=b2
        B(i,3,4)=b1
        B(i,3,5)=c2
        B(i,3,6)=c1

        do 30 ii=1,6
          do 31 jj=1,6
            kele(ii,jj)=0.d0
            do 32 kk=1,3
              do 33 ll=1,3
                kele(ii,jj)=kele(ii,jj)+
     1                      AA/2.d0*B(i,kk,ii)*CC(kk,ll)*B(i,ll,jj)
 33           continue
 32         continue
 31       continue
 30     continue

        do 21 inode=1,3
          nodei=eledata(i,inode)
          do 22 idof=1,2
            do 23 jnode=1,3
              nodej=eledata(i,jnode)
              do 24 jdof=1,2
                kglob(gcon(nodei,idof),gcon(nodej,jdof))=
     1            kglob(gcon(nodei,idof),gcon(nodej,jdof))+
     2            kele((inode-1)*2+idof,(jnode-1)*2+jdof)
 24           continue
 23         continue
 22       continue
 21     continue
 20   continue
        
      return
      end

      subroutine applyf(maxnodes,nnodes,nforces,fnode,fdof,fval,gcon,f)
      real*8 fval(2*maxnodes),f(2*maxnodes)
      integer i,fnode(2*maxnodes),fdof(2*maxnodes),nnodes,dof,maxnodes
      integer nforces,gcon(maxnodes,2)
      
      do 10 i=1,2*nnodes
        f(i)=0.d0
 10   continue
      
      do 20 i=1,nforces
        dof=gcon(fnode(i),fdof(i))
        f(dof)=fval(i)
 20   continue
 
      return
      end
      
      subroutine applydisp(maxnodes,nnodes,ndisp,unode,udof,uval,gcon,
     1                     kglob,f)
      real*8 kglob(2*maxnodes,2*maxnodes),uval(2*maxnodes),kpenalty
      real*8 f(2*maxnodes),maxk
      integer maxnodes,ndisp,gcon(maxnodes,2)
      integer i,unode(2*maxnodes),udof(2*maxnodes),nnodes,dof
      
      maxk=0.d0
      do 10 i=1,2*nnodes
        if (kglob(i,i).gt.maxk) maxk=kglob(i,i)
 10   continue
      kpenalty=1.d6*maxk
      
      do 20 i=1,ndisp
        dof=gcon(unode(i),udof(i))
        kglob(dof,dof)=kglob(dof,dof)+kpenalty
        f(dof)=kpenalty*uval(i)
 20   continue
 
      return
      end
 
      subroutine structure(maxelem,maxnodes,nelem,eledata,nodedata,u,v)
      integer eledata(maxelem,3),maxelem,maxnodes,nelem,i,nodei,nodej
      integer nodek
      real*8 nodedata(maxnodes,2),u(maxnodes),v(maxnodes)
      
      open(unit=21,file='undeformed')
      open(unit=22,file='deformed')
      
      do 10 i=1,nelem
        nodei=eledata(i,1)
        nodej=eledata(i,2)
        nodek=eledata(i,3)
        write(21,*)nodedata(nodei,1),nodedata(nodei,2)
        write(21,*)nodedata(nodej,1),nodedata(nodej,2)
        write(21,*)nodedata(nodek,1),nodedata(nodek,2)
        write(21,*)nodedata(nodei,1),nodedata(nodei,2)
        write(21,*)' '
        write(22,*)nodedata(nodei,1)+u(nodei),
     1             nodedata(nodei,2)+v(nodei)
        write(22,*)nodedata(nodej,1)+u(nodej),
     1             nodedata(nodej,2)+v(nodej)
        write(22,*)nodedata(nodek,1)+u(nodek),
     1             nodedata(nodek,2)+v(nodek)
        write(22,*)nodedata(nodei,1)+u(nodei),
     1             nodedata(nodei,2)+v(nodei)
        write(22,*)' '
 10   continue
 
      return
      end
      
      subroutine stress(maxelem,maxnodes,nelem,eledata,nodedata,u,v,
     1                  strn,strs,EE,vv,B)
      integer eledata(maxelem,3),maxelem,maxnodes,nelem,i,j,k
      real*8 nodedata(maxnodes,2),u(maxnodes),v(maxnodes)
      real*8 strn(maxelem,3),strs(maxelem,3),EE,vv,CC(3,3)
      real*8 B(maxelem,3,6),uu(6),x,y
      
      open(unit=21,file='strain')
      open(unit=22,file='stress')
      
      CC(1,1)=EE/(1.d0-vv*vv)
      CC(1,2)=vv*CC(1,1)
      CC(1,3)=0.d0
      CC(2,1)=CC(1,2)
      CC(2,2)=CC(1,1)
      CC(2,3)=0.d0
      CC(3,1)=0.d0
      CC(3,2)=0.d0
      CC(3,3)=EE/2.d0/(1.d0+vv)
      
      do 100 i=1,nelem
        uu(1)=u(eledata(i,1))
        uu(2)=v(eledata(i,1))
        uu(3)=u(eledata(i,2))
        uu(4)=v(eledata(i,2))
        uu(5)=u(eledata(i,3))
        uu(6)=v(eledata(i,3))
        do 110 j=1,3
          strn(i,j)=0.d0
          do 120 k=1,6
            strn(i,j)=strn(i,j)+B(i,j,k)*uu(k)
 120      continue
 110    continue
        do 130 j=1,3
          strs(i,j)=0.d0
          do 140 k=1,3
            strs(i,j)=strs(i,j)+CC(j,k)*strn(i,k)
 140      continue
 130    continue
        x=(nodedata(eledata(i,1),1)+nodedata(eledata(i,2),1)+
     1     nodedata(eledata(i,3),1))/3.d0
        y=(nodedata(eledata(i,1),2)+nodedata(eledata(i,2),2)+
     1     nodedata(eledata(i,3),2))/3.d0
        if ((nodedata(eledata(i,1),1).lt.1.d-6.and.
     1       nodedata(eledata(i,2),1).lt.1.d-6).or.
     2      (nodedata(eledata(i,1),1).lt.1.d-6.and.
     1       nodedata(eledata(i,3),1).lt.1.d-6).or.
     2      (nodedata(eledata(i,2),1).lt.1.d-6.and.
     1       nodedata(eledata(i,3),1).lt.1.d-6)) then
          write(21,*)y,strs(i,1),strs(i,2)
        endif
        if ((nodedata(eledata(i,1),2).lt.1.d-6.and.
     1       nodedata(eledata(i,2),2).lt.1.d-6).or.
     2      (nodedata(eledata(i,1),2).lt.1.d-6.and.
     1       nodedata(eledata(i,3),2).lt.1.d-6).or.
     2      (nodedata(eledata(i,2),2).lt.1.d-6.and.
     1       nodedata(eledata(i,3),2).lt.1.d-6)) then
          write(22,*)x,strs(i,1),strs(i,2)
        endif
 100  continue
 
      return
      end
      
      SUBROUTINE gaussj(a,n,np,b,m)
      INTEGER m,n,np,NMAX
      REAL*8 a(np,np),b(np)
      PARAMETER (NMAX=10000)
      INTEGER i,icol,irow,j,k,l,ll,indxc(NMAX),indxr(NMAX),ipiv(NMAX)
      REAL*8 big,dum,pivinv
      do 11 j=1,n
        ipiv(j)=0
11    continue
      do 22 i=1,n
        big=0.d0
        do 13 j=1,n
          if(ipiv(j).ne.1)then
            do 12 k=1,n
              if (ipiv(k).eq.0) then
                if (dabs(a(j,k)).ge.big)then
                  big=dabs(a(j,k))
                  irow=j
                  icol=k
                endif
              else if (ipiv(k).gt.1) then
                pause 'singular matrix in gaussj'
              endif
12          continue
          endif
13      continue
        ipiv(icol)=ipiv(icol)+1
        if (irow.ne.icol) then
          do 14 l=1,n
            dum=a(irow,l)
            a(irow,l)=a(icol,l)
            a(icol,l)=dum
14        continue
          do 15 l=1,m
            dum=b(irow)
            b(irow)=b(icol)
            b(icol)=dum
15        continue
        endif
        indxr(i)=irow
        indxc(i)=icol
        if (a(icol,icol).eq.0.d0) pause 'singular matrix in gaussj'
        pivinv=1.d0/a(icol,icol)
        a(icol,icol)=1.d0
        do 16 l=1,n
          a(icol,l)=a(icol,l)*pivinv
16      continue
        do 17 l=1,m
          b(icol)=b(icol)*pivinv
17      continue
        do 21 ll=1,n
          if(ll.ne.icol)then
            dum=a(ll,icol)
            a(ll,icol)=0.d0
            do 18 l=1,n
              a(ll,l)=a(ll,l)-a(icol,l)*dum
18          continue
            do 19 l=1,m
              b(ll)=b(ll)-b(icol)*dum
19          continue
          endif
21      continue
22    continue
      do 24 l=n,1,-1
        if(indxr(l).ne.indxc(l))then
          do 23 k=1,n
            dum=a(k,indxr(l))
            a(k,indxr(l))=a(k,indxc(l))
            a(k,indxc(l))=dum
23        continue
        endif
24    continue
      return
      END

