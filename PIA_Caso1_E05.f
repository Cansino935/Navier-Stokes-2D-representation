      program PIA_Caso1_E05
       !Este programa calcula la sol. num. de una e.d.p. (Uyy=Cte.)
       !Con u=u(y)

      implicit none
      real*8 dy, dl, h,la
      integer*8 i, j, k, n
      real*8, allocatable, dimension (:,:):: a,c !nxn
      real*8, allocatable, dimension (:,:):: b1,d1,d2  !nx1
      
!****************************************************************
      !Fijacion de parametros
      write(*,*)"******************************************************"
      write(*,*)"Fijacion de parametros:"
      write(*,*)"Este programa calcula la sol. num. de una e.d.p.
     &(Uyy=Cte.), con u=u(y)"
      write(*,*)"Donde la e.d.p proviene del planteamiento del caso: 1"
      write(*,*)"******************************************************"
      write(*,*)"Escriba el numero de subintervalos (N) entero:"
      read(*,*) n
      write(*,*)"Escriba la longuitud en y (h)"
      read(*,*) h
      write(*,*)"Escriba el termino cte. lamda (lamda<0):"
      read(*,*) la
      !write(*,*)"Si N>8, escriba 1:"
      !read(*,*)k
      write(*,*)"******************************************************"

!****************************************************************
      !apertura de archivos y escritura de script
      OPEN( 7, FILE= "Vectorsol.txt")
      OPEN( 8, FILE="PIA_Caso1.gpl")

      write(8,*)"set zeroaxis"
      write(8,*)"set xlabel ""y"" font ""Times-Roman, 10"""
      write(8,*)"set ylabel ""u(y)"" font ""Times-Roman, 10"""
      write(8,*)"set xrange [-1:",h+1,"]"
      write(8,*)"set yrange [-1:",h+1,"]"
      write(8,*)"set title ""Pia Caso:1"""
      write(8,*)"unset key"
      write(8,*)"set xtics font ""Times-Roman, 10"""
      write(8,*)"set ytics font ""Times-Roman, 10"""
      write(8,*)"h(x)=",-0.5*la*(h**2),"*( (x/",h,")- (x**2/",h**2,") )"
      write(8,*)"plot ""Vectorsol.txt"""
      write(8,*)"replot h(x)"
!****************************************************************
      !incremento en y
      dy=h/n
      !fijacion de elementos
      allocate(a(n,n))
      allocate(b1(n,1))
      allocate(c(n,n))
      allocate (d1(n,1))
      allocate (d2(n,1))
!*****************************************************************
      !calculo de matriz a
      call matriz1(a,n)
      write (*,200) !Escritura de a

       Do i=1,n
       write (*,201) (a(i,j),j=1,n)
       end do
      write(*,*)"******************************************************"

!*****************************************************************
      !calculo de matriz inv
      call matriz3(a,c,n)
      write (*,202)

       Do i = 1,n      !impresion de matriz inversa
        write (*,201)  (c(i,j),j=1,n)
       end do
      write(*,*)"******************************************************"

!*****************************************************************
      !vector 1
      call matriz2(b1,n,la,dy)
      write(*,*)"Vector 1:"
        Do i = 1,n      !impresion de matriz inversa
        write (*,201) (b1(i,j),j=1,1)
        end do
      write(*,*)"******************************************************"
      
!*****************************************************************
      !vector inc. en y
      call matriz5 (d2,dy,n)
      write(*,*)"Inc. en y:"
        Do i = 1,n      !impresion de incrementos en y
        write (*,201) (d2(i,j),j=1,1)
        end do
      write(*,*)"******************************************************"

!*****************************************************************
      !vector sol.
      call matriz4 (c,b1,d1,n)
      write(*,*)"Vector sol.:"
        Do i = 1,n      !impresion del vector sol-
        write (*,201) (d1(i,j),j=1,1)
        write(7,*)  (d2(i,j),j=1,1), (d1(i,j),j=1,1)
        end do
      write(*,*)"******************************************************"

!*****************************************************************
      close(7)
      close(8)
200    format (' Matriz A:')
201    format (6f12.6)
202    format (/,' Inversa matriz A^{-1}:')

      call system("PIA_Caso1.gpl")
      
       write(*,*)"Aqui acaba el programa:)"
       Pause
       Stop
       Continue
       end
       
!***********************************************************************
      subroutine matriz1(a,n)
      integer*8 i, j, k, n
      double precision a(n,n)

      Do i=1, n    !Llenado de a
      Do j=1, n

      if (i.ge.2) then
      if (i.eq.j) then  !coeficientes
      a(i,j)= -2
      a(i-1,j) = 1
      if (i.lt.n) then
      a(i+1,j)= 1
      end if
      end if
      else
      a(i,j) = 0        !los demas elementos
      end if
      a(1,1)= 1      !condiciones de borde u(0)
      a(1,2)= 0
      a(2,1)= 1
      a(n,n)= 1      !condiciones de borde u(h)
      a(n,n-1)= 0
      end do
      end do

      end subroutine
!***********************************************************************
      subroutine matriz2(b1,n,la,dy)
      real*8 dy, dl, h,la
      integer*8 i, j, k, n
      double precision b1(n,1)

      Do i=1,n
      Do j=1,1

      if(i.eq.1) then !condicion de borde
      b1(i,j)= 0

      else
      b1(i,j)= la*(dy**2)
      end if

      b1(n,1)=0   !condicion de borde
      end do
      end do
      end subroutine
!***********************************************************************
       subroutine matriz3(a,c,n)
       implicit none
       integer n
       double precision a(n,n), c(n,n)
       double precision L(n,n), U(n,n), b(n), d(n), x(n)
       double precision coeff
       integer i, j, k

       L=0.0
       U=0.0
       b=0.0

       do k=1, n-1
        do i=k+1,n
         coeff=a(i,k)/a(k,k)
           L(i,k) = coeff
        do j=k+1,n
         a(i,j) = a(i,j)-coeff*a(k,j)
       end do
       end do
       end do


       do i=1,n
        L(i,i) = 1.0
       end do

        do j=1,n
        do i=1,j
         U(i,j) = a(i,j)
        end do
        end do

        do k=1,n
          b(k)=1.0
          d(1) = b(1)

        do i=2,n
        d(i)=b(i)
        do j=1,i-1
         d(i) = d(i) - L(i,j)*d(j)
        end do
        end do

        x(n)=d(n)/U(n,n)
       do i = n-1,1,-1
        x(i) = d(i)
       do j=n,i+1,-1
        x(i)=x(i)-U(i,j)*x(j)
        end do
        x(i) = x(i)/u(i,i)
        end do

         do i=1,n
         c(i,k) = x(i)
         end do
         b(k)=0.0
        end do
        end subroutine
!***********************************************************************
      subroutine matriz4 (c,b1,d1,n)
      integer*8 i, j, k, n
      double precision c(n,n), b1(n,1), d1(n,1)

      !calculo del vector sol.  i.e. multiplicacion de c*b1
      Do i=1, n
      Do j=1, 1

      d1(i,j)= 0

      Do k=1, n
      d1(i,j)=  d1(i,j)+ ( c(i,k)*b1(k,j) )
      end do

      end do
      end do
      end subroutine
!***********************************************************************
      subroutine matriz5 (d2,dy,n)
      integer*8 i, j, k, n
      real*8 dy
      double precision d2(n,1)

      Do i=1, n  !calculo de incrementos
      Do j=1, 1

      d2(i,j)=i*dy

      end do
      end do
      end subroutine
!***********************************************************************
