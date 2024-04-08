module math

implicit none

contains

function math_axisAngleToR(axis,omega) result(math_axisAngleToR1)

  implicit none
  
  real, dimension(3), intent(in) :: axis
  real, intent(in) :: omega
  real, dimension(3) :: n
  real :: norm,s,c,c1
  
  real, dimension(3,3), parameter :: &
  I3 = real(reshape([&
    1, 0, 0, &
    0, 1, 0, &
    0, 0, 1  &
    ],shape(I3)))                                                                      !< 3x3 Identity

  real, dimension(3,3) :: math_axisAngleToR1
  
  norm = norm2(axis)
  wellDefined: if (norm > 1.0e-8) then
    n = axis/norm                                                                                    ! normalize axis to be sure
  
    s = sin(omega)
    c = cos(omega)
    c1 = 1.0 - c
  
    math_axisAngleToR1(1,1) =  c + c1*n(1)**2.0
    math_axisAngleToR1(1,2) =  c1*n(1)*n(2) - s*n(3)
    math_axisAngleToR1(1,3) =  c1*n(1)*n(3) + s*n(2)
                              
    math_axisAngleToR1(2,1) =  c1*n(1)*n(2) + s*n(3)
    math_axisAngleToR1(2,2) =  c + c1*n(2)**2.0
    math_axisAngleToR1(2,3) =  c1*n(2)*n(3) - s*n(1)
                              
    math_axisAngleToR1(3,1) =  c1*n(1)*n(3) - s*n(2)
    math_axisAngleToR1(3,2) =  c1*n(2)*n(3) + s*n(1)
    math_axisAngleToR1(3,3) =  c + c1*n(3)**2.0
  else wellDefined
    math_axisAngleToR1 = I3
  endif wellDefined
  
end function

end module math


program corresponcence_matrix
  use math
  implicit none

    integer, dimension(4) :: &
        active = [6,6,6,6], &                 !< number of active twin systems
        potential = [6,6,6,6]                 !< all the potential twin systems

    real, dimension(3) :: &
        direction, normal

    real, dimension(3,24) :: normal_vector, direction_vector

    real, dimension(3,3,24) :: SchmidMatrix, corresponcenceMatrix

    real, dimension(24) :: characteristicShear

    real :: cOverA = 1.6235

    real :: pi = 3.14159274 

    real, dimension(8,24) :: &
    system = reshape(real([&
    ! <-10.1>{10.2} systems, shear = (3-(c/a)^2)/(sqrt(3) c/a)
    ! tension in Co, Mg, Zr, Ti, and Be; compression in Cd and Zn
      -1,  0,  1,  1,     1,  0, -1,  2, & !
       0, -1,  1,  1,     0,  1, -1,  2, &
       1, -1,  0,  1,    -1,  1,  0,  2, &
       1,  0, -1,  1,    -1,  0,  1,  2, &
       0,  1, -1,  1,     0, -1,  1,  2, &
      -1,  1,  0,  1,     1, -1,  0,  2, &
    ! <11.6>{-1-1.1} systems, shear = 1/(c/a)
    ! tension in Co, Re, and Zr
      -1, -1,  2,  6,     1,  1, -2,  1, &
       1, -2,  1,  6,    -1,  2, -1,  1, &
       2, -1, -1,  6,    -2,  1,  1,  1, &
       1,  1, -2,  6,    -1, -1,  2,  1, &
      -1,  2, -1,  6,     1, -2,  1,  1, &
      -2,  1,  1,  6,     2, -1, -1,  1, &
    ! <10.-2>{10.1} systems, shear = (4(c/a)^2-9)/(4 sqrt(3) c/a)
    ! compression in Mg
       1,  0, -1, -2,     1,  0, -1,  1, &
       0,  1, -1, -2,     0,  1, -1,  1, &
      -1,  1,  0, -2,    -1,  1,  0,  1, &
      -1,  0,  1, -2,    -1,  0,  1,  1, &
       0, -1,  1, -2,     0, -1,  1,  1, &
       1, -1,  0, -2,     1, -1,  0,  1, &
    ! <11.-3>{11.2} systems, shear = 2((c/a)^2-2)/(3 c/a)
    ! compression in Ti and Zr
       1,  1, -2, -3,     1,  1, -2,  2, &
      -1,  2, -1, -3,    -1,  2, -1,  2, &
      -2,  1,  1, -3,    -2,  1,  1,  2, &
      -1, -1,  2, -3,    -1, -1,  2,  2, &
       1, -2,  1, -3,     1, -2,  1,  2, &
       2, -1, -1, -3,     2, -1, -1,  2  &
      ]),shape(system))

    real, dimension(3,3), parameter :: &
    I3 = real(reshape([&
      1, 0, 0, &
      0, 1, 0, &
      0, 0, 1  &
      ],shape(I3)))                                                                      !< 3x3 Identity

    integer :: &
    a, &                                                                        !< index of active system
    p, &                                                                        !< index in potential system matrix
    f, &                                                                        !< index of my family
    s, f1, s1, e1, i, j, k                                                                      !< index of my system in current family


    
    
    a = 0
    do f = 1, size(active,1)                                                    !< Loops 1 to 4 for hP
        do s = 1, active(f)                                                     !< 1 to 6 two times
            
            a = a + 1
            p = sum(potential(1:f-1))+s                                         !< 1 to 6 and 13 to 18
            
            direction = [ system(1,p)*1.5, &
            (system(1,p)+2.0*system(2,p))*sqrt(0.75), &
             system(4,p)*cOverA ]                                               ! direction [uvtw]->[3u/2 (u+2v)*sqrt(3)/2 w*(p/a)])

            normal    = [ system(5,p), &
            (system(5,p)+2.0*system(6,p))/sqrt(3.0), &
            system(8,p)/cOverA ]                                              ! plane (hkil)->(h (h+2k)/sqrt(3) l/(p/a))
            

            
            normal_vector(1:3,a) = normal /norm2(normal)
            direction_vector(1:3,a) = direction / norm2(direction)
            
        end do
    end do

    do f1 = 1,size(active,1)
      s1 = sum(active(:f1-1)) + 1
      e1 = sum(active(:f1))
      select case(f1)
      case (1)
        characteristicShear(s1:e1) = (3.0-cOverA**2)/sqrt(3.0)/cOverA
      case (2)
        characteristicShear(s1:e1) = 1.0/cOverA
      case (3)
        characteristicShear(s1:e1) = (4.0*cOverA**2-9.0)/sqrt(48.0)/cOverA
      case (4)
        characteristicShear(s1:e1) = 2.0*(cOverA**2-2.0)/3.0/cOverA
      end select
    enddo
    
    write(6,*)'characteristic shear, 1,  0, -1,  1,    -1,  0,  1,  2, = ',characteristicShear(4)
    write(6,*)'characteristic shear, -1, -1,  2,  6,     1,  1, -2,  1, =',characteristicShear(7)
    write(6,*)'characteristic shear, 1,  0, -1, -2,     1,  0, -1,  1, =',characteristicShear(13)
    write(6,*)'characteristic shear, 1,  1, -2, -3,     1,  1, -2,  2, =',characteristicShear(19)


    do i = 1, sum(active)
      forall(j=1:3, k=1:3) SchmidMatrix(j,k,i) = direction_vector(j,i) * normal_vector(k,i)
    enddo

    do i = 1, sum(active)
      corresponcenceMatrix(1:3,1:3,i) = matmul(math_axisAngleToR(normal_vector(1:3,i),pi),&
                                               I3+characteristicShear(i)*SchmidMatrix(1:3,1:3,i))
      
    enddo
    
    write(6,*)'correspondence matrix for 1,  0, -1,  1,    -1,  0,  1,  2 =====',corresponcenceMatrix(1:3,1:3,4)
    write(6,*)'ooooooooooooooooooooooooooooooooooooooooooooooo'
    write(6,*)'correspondence matrix for -1, -1,  2,  6,     1,  1, -2,  1 =====',corresponcenceMatrix(1:3,1:3,7)
    write(6,*)'ooooooooooooooooooooooooooooooooooooooooooooooo'
    write(6,*)'correspondence matrix for 1,  0, -1, -2,     1,  0, -1,  1 =====',corresponcenceMatrix(1:3,1:3,13)
    write(6,*)'ooooooooooooooooooooooooooooooooooooooooooooooo'
    write(6,*)'correspondence matrix for 1,  1, -2, -3,     1,  1, -2,  2 =====',corresponcenceMatrix(1:3,1:3,19)


end program corresponcence_matrix

