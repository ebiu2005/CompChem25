      implicit real*8(a-h,o-z)

      dimension x(500),xnew(500)
      character*80 conf
      character*2 aa(100)

c   input

      open(unit=5,file='rotate.in')

      read(5,*) conf
      read(5,*) natom
      read(5,*) th,ax,ay,az

      close(5)
      
c   store initial conf

      open(unit=7,file=conf)

      do i=1,natom

      ii=3*(i-1)

      read(7,*) aa(i),xx,y,z

      x(1+ii)=xx
      x(2+ii)=y
      x(3+ii)=z

      enddo

      close(7)

c  rotate the structure

      do i=1,natom

      ii=3*(i-1)
    
      xx=x(1+ii)
      y=x(2+ii)
      z=x(3+ii)
    
      call rotaxis(th,ax,ay,az,xx,y,z)

      write(*,*) aa(i),xx,y,z

      enddo

      stop
      end

c-----------------------------------------------------------------

      subroutine rotaxis(th,ax,ay,az,xx,y,z)

      implicit real*8(a-h,o,z)

      dimension ee(3),a(3)

C rotate dh into dp around the axis ee = dh x dp (made it a unit vector)

c      ee(1)=dh(2)*dp(3)-dh(3)*dp(2)
c      ee(2)=-(dh(1)*dp(3)-dh(3)*dp(1))
c      ee(3)=dh(1)*dp(2)-dh(2)*dp(1)

      ee(1)=ax
      ee(2)=ay
      ee(3)=az

c   make sure it is a unit vector

      eee=ee(1)*ee(1)+ee(2)*ee(2)+ee(3)*ee(3)
      eee=dsqrt(eee)

      do i=1,3
      ee(i)=ee(i)/eee
      enddo

C apply formulas for R

      ct=dcos(th)
      dt=1.d0-dcos(th)
      st=dsin(th)

      R(1,1)=ct+ee(1)*ee(1)*dt
      R(1,2)=ee(1)*ee(2)*dt-ee(3)*st
      R(1,3)=ee(1)*ee(3)*dt+ee(2)*st
      R(2,1)=ee(2)*ee(1)*dt+ee(3)*st
      R(2,2)=ct+ee(2)*ee(2)*dt
      R(2,3)=ee(2)*ee(3)*dt-ee(1)*st
      R(3,1)=ee(3)*ee(1)*dt-ee(2)*st
      R(3,2)=ee(3)*ee(2)*dt+ee(1)*st
      R(3,3)=ct+ee(3)*ee(3)*dt

C     rotation

      a(1)=xx
      a(2)=y
      a(3)=z

      aa(1)=0.d0
      aa(2)=0.d0
      aa(3)=0.d0

      do i=1,3
      
      do j=1,3
      aa(i)=aa(i)+R(i,j)*a(j)
      enddo

      xx=aa(1)
      y=aa(2)
      z=aa(3)

      return
      end


