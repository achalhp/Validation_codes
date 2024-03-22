program corresponcence_matrix
    implicit none
    integer, dimension(4) :: &
        active = [6,0,6,0], &
        potential = [6,6,6,6]

    real, dimension(3) :: &
        direction, normal

    real, dimension(3,12) :: normal_vector

    real :: cOverA = 1.6235

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

    integer :: &
    a, &                                                                        !< index of active system
    p, &                                                                        !< index in potential system matrix
    f, &                                                                        !< index of my family
    s                                                                           !< index of my system in current family
    
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
            

            
            normal_vector(1:3,a) = normal   /norm2(normal)
            write(6,*)'normal vector', normal_vector
        end do
    end do

    do f = 1,size(active,1)
    enddo

end program corresponcence_matrix