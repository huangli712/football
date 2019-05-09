
   program ball
     implicit none

! global parameters
     integer, parameter :: nsite = 60      ! number of sites
     integer, parameter :: nedge = 90      ! number of edges
     !
     integer, parameter :: nsweep = 100000 ! number of Monte Carlo sweeps
     integer, parameter :: nwarm = 1000    ! number of warm up steps
     integer, parameter :: naver = 100     ! interval for successive measurements
     !
! local variables
     integer :: ip      ! loop index for Monte Carlo sweeps
     integer :: is      ! loop index for site
     integer :: i, j, k ! point to the three neighbours
     integer :: naccu   ! number of successful measurements
     integer :: spin    ! spin in the current site
     integer :: spin_i  ! spin in the 1st neighbour site
     integer :: spin_j  ! spin in the 2nd neighbour site
     integer :: spin_k  ! spin in the 3rd neighbour site
     integer :: nseed   ! number of seeds
     !
     integer :: ising(nsite)     ! ising fields
     integer :: neigh(nsite,3)   ! list of neighbours
     !
     real(8) :: ene              ! energy for the current configuration
     real(8) :: ediff            ! energy difference for two configuration
     real(8) :: etot             ! total energy 
     real(8) :: lnz              ! ln Z
     real(8) :: p                ! transition probability
     real(8) :: r                ! random number
     !
     integer, allocatable :: seed(:) ! seeds for random number generator

! init seed
     call random_seed()
     call random_seed(size=nseed)
     allocate(seed(nseed))

! get default seed
     call random_seed(get=seed)

! update seed
     do i=1,nseed 
        call system_clock(j)
        seed(i)=i * j + 2019
     enddo
     call random_seed(put=seed)

! init ising field 
     call init_ising(nsite,ising)

! init neighbours
     call init_neigh(nsite,neigh)

! try to calculate energy for the initial configuration
     call init_energy(nsite, ising, neigh, ene)

! reset some key observables
     etot = 0.0d0
     lnz = 0.0d0
     naccu = 0

     do ip=1,nsweep                   ! Monte Carlo sweeps
         do is=1,nsite                ! loop over all sites
             !
             spin = ising(is)         ! get current spin
             !
             i = neigh(is,1)          ! get three neighbours
             j = neigh(is,2) 
             k = neigh(is,3)
             !
             spin_i = ising(i)        ! get spins for the three neighbours
             spin_j = ising(j)
             spin_k = ising(k)
             !
             ediff = 2.0d0 * spin * ( spin_i + spin_j + spin_k ) ! calculate energy difference
             !
             p = exp(-ediff)          ! calculate transition probability
             !
             call random_number(r)    ! Metropolis algorithm
             if ( min(1.0d0,p) > r ) then ! accept
                 ising(is) = -spin    ! flip the spin, update the configuration
                 ene = ene + ediff    ! update the energy
             endif
             !
             if ( ip > nwarm .and. mod(ip, naver) == 0 ) then
                 etot = etot + ene    ! try to measure the energy and lnZ
                 lnz = lnz + exp(-ene)
                 naccu = naccu + 1    ! do not update the counter
             endif
             !
         enddo
     enddo
     print *, 'energy:', etot / naccu, 'lnZ/N:', log( lnz / naccu ) / nsite
!<     print *, 'final ising fields:'
!<     print *, ising

   end program ball

! init the ising fields
   subroutine init_ising(nsite, ising)
     implicit none

     integer, intent(in)  :: nsite
     integer, intent(out) :: ising(nsite)

     integer :: i
     real(8) :: r

     do i=1,nsite
         call random_number(r)
         if ( r > 0.5d0 ) then
             ising(i) = +1
         else
             ising(i) = -1
         endif
     enddo

     return
   end subroutine init_ising

! init the neighbours
! now the implement is very ugly and dirty. but it works
! it tooks me half and hour to write it. i am too stupid.
   subroutine init_neigh(nsite,neigh)
     implicit none

     integer, intent(in)  :: nsite
     integer, intent(out) :: neigh(nsite,3)

     integer :: i, j, k, ind1, ind2
     logical :: find
     integer :: hist(nsite)

! check the size of array
     if ( nsite /= 60 ) STOP

     neigh(1,1) = 2
     neigh(1,2) = 5
     neigh(1,3) = 6

     neigh(2,1) = 1
     neigh(2,2) = 3
     neigh(2,3) = 7

     neigh(3,1) = 2
     neigh(3,2) = 4
     neigh(3,3) = 8

     neigh(4,1) = 3
     neigh(4,2) = 5
     neigh(4,3) = 9

     neigh(5,1) = 1
     neigh(5,2) = 4
     neigh(5,3) = 10

     neigh(6,1) = 1
     neigh(6,2) = 11
     neigh(6,3) = 12

     neigh(7,1) = 2
     neigh(7,2) = 13
     neigh(7,3) = 14

     neigh(8,1) = 3
     neigh(8,2) = 15
     neigh(8,3) = 16

     neigh(9,1) = 4
     neigh(9,2) = 17
     neigh(9,3) = 18

     neigh(10,1) = 5
     neigh(10,2) = 19
     neigh(10,3) = 20

     neigh(11,1) = 6
     neigh(11,2) = 20
     neigh(11,3) = 21

     neigh(12,1) = 6
     neigh(12,2) = 22
     neigh(12,3) = 13

     neigh(13,1) = 7
     neigh(13,2) = 12
     neigh(13,3) = 23

     neigh(14,1) = 7
     neigh(14,2) = 24
     neigh(14,3) = 15

     neigh(15,1) = 8
     neigh(15,2) = 14
     neigh(15,3) = 25

     neigh(16,1) = 8
     neigh(16,2) = 26
     neigh(16,3) = 17

     neigh(17,1) = 9
     neigh(17,2) = 16
     neigh(17,3) = 27

     neigh(18,1) = 9
     neigh(18,2) = 28
     neigh(18,3) = 19

     neigh(19,1) = 10
     neigh(19,2) = 18
     neigh(19,3) = 29

     neigh(20,1) = 10
     neigh(20,2) = 11
     neigh(20,3) = 30

     neigh(21,1) = 11
     neigh(21,2) = 31
     neigh(21,3) = 22

     neigh(22,1) = 12
     neigh(22,2) = 21
     neigh(22,3) = 32

     neigh(23,1) = 13
     neigh(23,2) = 33
     neigh(23,3) = 24

     neigh(24,1) = 14
     neigh(24,2) = 23
     neigh(24,3) = 34

     neigh(25,1) = 15
     neigh(25,2) = 26
     neigh(25,3) = 35

     neigh(26,1) = 16
     neigh(26,2) = 25
     neigh(26,3) = 36

     neigh(27,1) = 17
     neigh(27,2) = 37
     neigh(27,3) = 28

     neigh(28,1) = 27
     neigh(28,2) = 18
     neigh(28,3) = 38

     neigh(29,1) = 19
     neigh(29,2) = 39
     neigh(29,3) = 30

     neigh(30,1) = 20
     neigh(30,2) = 29
     neigh(30,3) = 40

     neigh(31,1) = 21
     neigh(31,2) = 40
     neigh(31,3) = 41

     neigh(32,1) = 22
     neigh(32,2) = 42
     neigh(32,3) = 33

     neigh(33,1) = 32
     neigh(33,2) = 23
     neigh(33,3) = 43

     neigh(34,1) = 24
     neigh(34,2) = 44
     neigh(34,3) = 35

     neigh(35,1) = 34
     neigh(35,2) = 25
     neigh(35,3) = 45

     neigh(36,1) = 26
     neigh(36,2) = 37
     neigh(36,3) = 46

     neigh(37,1) = 36
     neigh(37,2) = 27
     neigh(37,3) = 47

     neigh(38,1) = 28
     neigh(38,2) = 39
     neigh(38,3) = 48

     neigh(39,1) = 29
     neigh(39,2) = 38
     neigh(39,3) = 49

     neigh(40,1) = 30
     neigh(40,2) = 50
     neigh(40,3) = 31

     neigh(41,1) = 31
     neigh(41,2) = 42
     neigh(41,3) = 51

     neigh(42,1) = 41
     neigh(42,2) = 32
     neigh(42,3) = 52

     neigh(43,1) = 33
     neigh(43,2) = 52
     neigh(43,3) = 44

     neigh(44,1) = 43
     neigh(44,2) = 34
     neigh(44,3) = 53

     neigh(45,1) = 35
     neigh(45,2) = 53
     neigh(45,3) = 46

     neigh(46,1) = 45
     neigh(46,2) = 36
     neigh(46,3) = 54

     neigh(47,1) = 37
     neigh(47,2) = 54
     neigh(47,3) = 48

     neigh(48,1) = 38
     neigh(48,2) = 47
     neigh(48,3) = 55

     neigh(49,1) = 39
     neigh(49,2) = 55
     neigh(49,3) = 50

     neigh(50,1) = 49
     neigh(50,2) = 40
     neigh(50,3) = 51

     neigh(51,1) = 50
     neigh(51,2) = 41
     neigh(51,3) = 56

     neigh(52,1) = 42
     neigh(52,2) = 57
     neigh(52,3) = 43

     neigh(53,1) = 44
     neigh(53,2) = 45
     neigh(53,3) = 58

     neigh(54,1) = 46
     neigh(54,2) = 47
     neigh(54,3) = 59

     neigh(55,1) = 48
     neigh(55,2) = 49
     neigh(55,3) = 60

     neigh(56,1) = 57
     neigh(56,2) = 51
     neigh(56,3) = 60

     neigh(57,1) = 52
     neigh(57,2) = 56
     neigh(57,3) = 58

     neigh(58,1) = 53
     neigh(58,2) = 57
     neigh(58,3) = 59

     neigh(59,1) = 58
     neigh(59,2) = 54
     neigh(59,3) = 60

     neigh(60,1) = 59
     neigh(60,2) = 55
     neigh(60,3) = 56

! check A
! for a given site A, suppose that one of its neighbours is B
! A must be one of neighbours of B
     do i=1,nsite
         do j=1,3
             ind1 = neigh(i,j) ! get its neighbour
             find = .false.
             do k=1,3          ! loop over the three neighbours of the target neighbour
                               ! for the target neighbour, there must be one neighbour
                               ! that is the current site
                 ind2 = neigh(ind1,k)
                 if ( ind2 == i ) then
                     find = .true. 
                 else          ! yes, we find the current site
                     continue
                 endif
             enddo
             if ( find .eqv. .false. ) then
                 print *, i, j, neigh(i,j), ' failed'
                 STOP          ! fail to find the current site, we have to check it
!< DEBUG CODE
!<------------------------------------------------------------------------
!<             else
!<                 print *, i, j, neigh(i,j), ' pass '
!<------------------------------------------------------------------------
             endif
         enddo
     enddo

! check B
! for each site, it should be counted three times
     hist = 0
     do i=1,nsite
         do j=1,3
             ind1 = neigh(i,j)
             hist(ind1) = hist(ind1) + 1
         enddo
     enddo
!< DEBUG CODE
!<------------------------------------------------------------------------
!<     do i=1,nsite
!<         if ( hist(i) == 3 ) then
!<             print *, 'site: ', i, ' count:', hist(i), ' pass'
!<         else
!<             print *, 'site: ', i, ' count:', hist(i), ' failed'
!<             STOP
!<         endif
!<     enddo
!<------------------------------------------------------------------------

     return
   end subroutine init_neigh

! try to calculate the energy for the initial configuration
   subroutine init_energy(nsite, ising, neigh, ene)
     implicit none

     integer, intent(in)  :: nsite
     integer, intent(out) :: ising(nsite)
     integer, intent(in)  :: neigh(nsite,3)
     real(8), intent(out) :: ene

     integer :: i, j   ! loop index
     integer :: ind1   ! index for site A
     integer :: ind2   ! index for site B
     integer :: spin1  ! spin for site A
     integer :: spin2  ! spin for site B
     integer :: edge   ! number of edges, it should be 90

     edge = 0
     ene = 0.0d0

! go through all sites 
     do i=1,nsite
         ind1 = i
         spin1 = ising(ind1)
! go through all three neighbours
         do j=1,3
             ind2 = neigh(i,j)
             spin2 = ising(ind2)
             ! 
             ene = ene - spin1 * spin2 / 2.0d0
             edge = edge + 1
         enddo
     enddo
!<     print *, 'number of edges:', edge/2

     return
   end subroutine init_energy
