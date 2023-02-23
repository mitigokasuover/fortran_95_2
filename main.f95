PROGRAM Runge_Kutt
REAL(4) :: K, Kappa
DATA t_max/5./, p2/7./, Beta/1./, Kappa/6./, x0/.1/    !Див. завдання
! На екран
WRITE (*,*)' *** Project Runge-Kutt ***'
WRITE (*,*)' ****** Input  maximum steps i_max:'
READ (*,*) i_max
h= t_max/ (i_max-1) !Крок
OPEN(UNIT=1, FILE='rezult.txt', STATUS='replace') !До файлу
WRITE (1,*)' Модель реакції'
WRITE (1,*)' в''яжучих у водному середовищі - метод Рунге-Кутта.'
WRITE (1,*)'  h=', h, '       t_max=', t_max, '             i_max=', i_max
WRITE(1,*) '  p2=', p2, '   Beta=', Beta, '    Kappa=', Kappa,'   x0=',x0
WRITE(1,*) '_________________________________________________'
WRITE (1,*)'    i             t               K(t)    (p2-Beta/Kappa*x)      x(t)'
! Те саме на екран
WRITE (*,*)'  h=', h, '       t_max=', t_max, '             i_max=', i_max
WRITE(*,*) '  p2=', p2, '   Beta=', Beta, '    Kappa=', Kappa,'   x0=',x0
WRITE(*,*) '_________________________________________________'
WRITE (*,*)'    i             t               K(t)    (p2-Beta/Kappa*x)      x(t)'
x=x0 ! Початкові умови 
DO i=1, i_max  ! Розрахунок по вузлах
t=t+h !Поточний час
CALL KOEF (t, h, x, p2, Beta, Kappa, A, B, C, D) ! Розрахунок коеф. A, B, C, D
Dx =(A+2*B+2*C+D)/6   ! ?xi = 1/6·(Ai + 2Bi + 2Ci + Di)
x=x+Dx       ! xi+1 = xi+?xi
aK=K(t)  ! Подивимось K(t)
! Остача Р1    (p1 - Alfa/Kappa * x)·= 1   - завжди надлишок
w=p2-Beta/Kappa*x	! Остача P2
WRITE(1,*) i, t, aK, w, x ! результати до файлу
WRITE(*,*) i, t, aK, w, x ! результати до файлу
ENDDO ! Закінчення циклу

CLOSE(UNIT=1) !Файл збережений і закритий
STOP
END PROGRAM Runge_Kutt

! Підпрограма розрахунку коеф. A, B, C, D
SUBROUTINE KOEF(t, h, x, p2, Beta, Kappa, A, B, C, D)
REAL(4) :: Kappa
    A=h*f(t, x, p2, Beta, Kappa)       ! Ai = h · f(ti, xi)
    B=h*f(t+h/2, x+A/2, p2, Beta, Kappa)      ! Bi  = h ·  (ti + h/2, xi+Ai/2)
    C=h*f(t+h/2, x+B/2, p2, Beta, Kappa)      ! Ci = h ·  (ti + h/2, xi + Bi/2)
    D=h*f(t+ h/2, x+C, p2, Beta, Kappa)        ! Di = h ·  (ti + h/2, xi + Ci)
RETURN
END SUBROUTINE KOEF

!Права частина ДР
! f(t, x)=K(t) · (p1 - Alfa/Kappa * x)·(p2 - Beta/Kappa * x)
! (p1 - Alfa/Kappa * x)·= 1
FUNCTION f(t, x, p2, Beta, Kappa) 
    REAL(4) :: K, Kappa
       f=K(t)*(p2-Beta/Kappa*x)
END FUNCTION f

             REAL FUNCTION K(t) ! Див. завдання
                      K=12*t+3.2  ! Див. завдання К(t)
            END FUNCTION K
