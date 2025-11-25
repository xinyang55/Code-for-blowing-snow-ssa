! The following code is derived from the p-TOMCAT model, which is used to calculate 
! blowing-snow-related sea salt aerosol production flux at any given conditions.
! (see the work by Yang et al., 2019 and 2008).
! 
! This single Fortran script can be run independently, but you need to set up 
! temperature, wind speed, RH, and also tell which hemisphere (NH or SH) it is for.   
! To run the code: mpif90 -r8 -extend_source -O2 -convert big_endian Blowingsnow_ssa.f 
!
! Note: To use it in models, you will need an interface for the subroutine.
!
! Xin Yang: January 2022 
!
! A wind-dependent blowing-snow size distribution (for alpha and beta values) is added, 
! based on Ranjithkumar et al. (2025), Elementa: Science of the Anthropocene, 
! 13, 27, 10.1525/elementa.2024.00006.
!
! A bug in the calculation water content of each particle (wmass) is fixed 
! Xin Yang: 26/02/2024
! 
! A snow salinity distribution dataset for the NH is added, based on Macfarlane, et al. (2021),
! PANGAEA, https://doi.org/10.1594/PANGAEA.935934. 
!  
! Xin Yang: 25/11/2025
! 
      integer mich, ish
      real fcoe 
      real drynacl(1000,24), flux(1000,24) ! array to hold ssa size and flux
      real dryr(21), dryrup(21), dryrlow(21) ! dry radius of SSA, upper and lower size
      real fssa(21), fssa_m(21)  ! final SSA production flux in bins for number (particles/m^2/s) and mass (kg/m^2/s)
      real fa(24) ! Snow salinity distribution frequency (in fraction), corresponding to sa(24) 
      real sa(24) ! snow salinity in psu (in 24 bins) 
      real gambin(19) ! Look-up table for gamma value 
! Choose a mechanism for testing:  
! mech =1: original mechanism in Yang et al., 2008; namely dm/dt = constant (Yang et al., 2019);
! mech =2 (not recommended): alternative mechanism, namely dm/dt is proportional to r (see Yang et al.,2019)
        mich=1 !  or =2
        if (mich.ne.1.and.mich.ne.2) then
          print *, 'Wrong mechanism, model stopped!' 
          stop
        endif 
        if (mich.eq.1) then
         print *, 'Mechanism 1: dm/dt =constant'
        endif 
        if (mich.eq.2) then
         print *, 'Mechanism 2: dm/dt =r '
        endif
         print *, '****************'

! Choose which hemisphere is for: NH or SH, as snow salinity datasets differ for them.
        ish=1 ! NH
!        ish=2 ! SH          
        if (ish.ne.1.and.ish.ne.2) then
          print *, 'Wrong choice for ish (ish=1 for NH and =2 for SH), model stopped!' 
          stop
        endif 

! You need to set up the meteorology factors at surface layer:
! at 10 m: 
         t3d=263.0 ! temperature in Kelvin 
         RH=90.0 ! relative humidity in %
         u10=12.0 ! 10 m wind speed in m/s


         z10=10.0 ! height in m 
         zz0=5.6E-5 ! roughness height in m 
         u_star=u10*0.4/alog(z10/zz0) ! friction velocity u_star in m/s
!        print *, 'Temperature(K), RH(%), wind speed (m/s), friction velocity(m/s)'
         print *, t3d, RH, u10,u_star
         print *, '*********************'

         pi=3.1415926
         den=2.165 ! density of dry NaCl (g/cm^3)
         denice=0.9167 ! density of ice (g/cm^3)        

! --- bulk blowing snow sublimation flux calculation ----
! based on Dery and Yau (1999, 2001) parameterisations
! rvv: individual gas constant for water vapour (=461.5 J/kg/K) 
         rvv=461.5
! roo: ice density in kg/m^3
         roo=denice*1000 ! = 916.7 
! heat: latent heat (=2834-2839 J/g from T=0C to -30C)
         heat=2835.*1000.0  !XY: can be improved to allow T-dependent
! est: saturated water vapour pressure over water/ice in Pa
         est=100.*6.112*exp(17.67*(t3d-273.15)/((t3d-273.15)+243.5))
! fcoe and dcoe: coefficients, see pp103 in Rogers and Yau (1989): A Short Course 
! in Cloud Physics, Third Edition, Pergamon Press.
         fcoe=2.2E-2
         dcoe=2.0E-5
! Fk(in m s/kg): conductivity term associated with sublimation
         fk=(heat/rvv/t3d-1.0)*heat*roo/fcoe/t3d
! Fd (in m s/kg): diffusion term associated with sublimaiton
! 
         fd=roo*rvv*t3d/dcoe/est
! coefficients
         qa0=3.78407E-1
         qa1=-8.64089E-2
         qa2=-1.60570E-2
         qa3=7.25516E-4
         qa4=-1.25650E-1
         qa5=2.48430E-2
         qa6=-9.56871E-4
         qa7=1.24600E-2
         qa8=1.56862E-3
         qa9=-2.93002E-4

! note: thi must be in units of -1E-12m^2/s (Rogers and Yau 1989)
         thi=(RH-100.0)/100.0/2.0/(fk+fd)/(-1E-12)

! wind speed (m/s) threshold for blowing snow 
         ut0=6.975+0.0033*(t3d-273.15+27.27)**2.0

         if (u10.gt.ut0) then
            qb0=0.385*(1.0-6.975/u10)**2.59/u_star ! 
            qbsalt=0.385*(1.0-ut0/u10)**2.59/u_star !
            qnormal=qb0/qbsalt
         else
            qnormal=0.0
         endif

! calculate water sublimation flux from blowing snow
! coefficients
         qa0=3.78407E-1
         qa1=-8.64089E-2
         qa2=-1.60570E-2
         qa3=7.25516E-4
         qa4=-1.25650E-1
         qa5=2.48430E-2
         qa6=-9.56871E-4
         qa7=1.24600E-2
         qa8=1.56862E-3
         qa9=-2.93002E-4

         qs0=qa0+qa1*thi+qa2*thi**2+qa3*thi**3+qa4*u10
         qs1=qa5*thi*u10+qa6*thi**2*u10
         qs2=qa7*u10**2.0+qa8*thi*u10**2.0+qa9*u10**3.0
! qsx: sublimation flux in mm/day 
         qsx=(qs0+qs1+qs2)/qnormal
         if (qsx.lt.0) qf=0.  ! force a negative value to zero

! convert above sublimation flux to kg/m^2/s  
         qf=qsx/24.0/60.0/60.0/10.0/100.*1E4 

          print *, 'BS sublimation flux:'
          print *, '  mm/day,    kg/m^2/s'
          print *, qsx, qf
          print *, '*********************'

! --- end of blowing snow sublimation flux calculation ---

! --- snow salinity ---
! salinity (psu) bins ranging from 0.002 to 80 psu
         sa(1)=0.002
         sa(2)=0.004
         sa(3)=0.006
         sa(4)=0.008
         sa(5)=0.01
         sa(6)=0.02
         sa(7)=0.04
         sa(8)=0.06
         sa(9)=0.08
         sa(10)=0.1
         sa(11)=0.2
         sa(12)=0.4
         sa(13)=0.6
         sa(14)=0.8
         sa(15)=1.0
         sa(16)=2.0
         sa(17)=4.0
         sa(18)=6.0
         sa(19)=8.0
         sa(20)= 10.0  ! 10.0
         sa(21)= 20.0  ! 20.0
         sa(22)= 40.0  ! 40 psu
         sa(23)= 60.0  ! 60 psu
         sa(24)= 80.0  ! 80 psu0
! 
! Find the right snow salinity frequency distribution
! in the NH: salinity frequency distribution from the Arctic (Macfarlane, et al., 2021)
      If (ish.eq.1) then 
           fa(           1 )=  0.0075
           fa(           2 )=  0.0251
           fa(           3 )=  0.0176
           fa(           4 )=  0.0201
           fa(           5 )=  0.0478
           fa(           6 )=  0.0251
           fa(           7 )=  0.0126
           fa(           8 )=  0.005
           fa(           9 )=  0.0075
           fa(          10 )=  0.3040
           fa(          11 )=  0.0704
           fa(          12 )=  0.0276
           fa(          13 )=  0.0251
           fa(          14 )=  0.0302
           fa(          15 )=  0.0905
           fa(          16 )=  0.0603
           fa(          17 )=  0.0327
           fa(          18 )=  0.0402
           fa(          19 )=  0.0226
           fa(          20 )=  0.0879
           fa(          21 )=  0.0302
           fa(          22 )=  0.0050
           fa(          23 )=  0.0050
           fa(          24 )=  0.0
       endif

! in the SH: salinity frequency distribution from the Weddell Sea (Frey et al., 2020)
       if (ish.eq.2) then 
           fa(           1 )=  1.4084507E-02
           fa(           2 )=  1.8779343E-02
           fa(           3 )=  5.6338027E-02
           fa(           4 )=  5.6338027E-02
           fa(           5 )=  0.1126761
           fa(           6 )=  0.1314554
           fa(           7 )=  0.1032864
           fa(           8 )=  7.5117372E-02
           fa(           9 )=  0.1126761
           fa(          10 )=  8.4507041E-02
           fa(          11 )=  6.1032865E-02
           fa(          12 )=  9.3896715E-03
           fa(          13 )=  1.4084507E-02
           fa(          14 )=  3.2863848E-02
           fa(          15 )=  1.8779343E-02
           fa(          16 )=  3.2863848E-02
           fa(          17 )=  1.4084507E-02
           fa(          18 )=  1.4084507E-02
           fa(          19 )=  0.0000000E+00
           fa(          20 )=  1.4084507E-02
           fa(          21 )=  2.3474179E-02
           fa(          22 )=  0.0000000E+00
           fa(          23 )=  0.0000000E+00
           fa(          24 )=  0.0000000E+00
       endif
! --- end of snow salinity ---

! BS particle size distribution following a two-parameter gamma probability density function (Schmidt 1982)
! Set up a look-up table for gamma at different alpha values
gambin(1)=0.9314 ! alfa=1.8
gambin(2)=0.9618 ! alfa=1.9
gambin(3)=1.0000 ! alfa=2
gambin(4)=1.0465 ! alfa=2.1
gambin(5)=1.1018 ! alfa=2.2
gambin(6)=1.1667 ! alfa=2.3
gambin(7)=1.2422 ! alfa=2.4
gambin(8)=1.3293 ! alfa=2.5
gambin(9)=1.4296 ! alfa=2.6
gambin(10)=1.5447 ! alfa=2.7
gambin(11)=1.6765 ! alfa=2.8
gambin(12)=1.8274 ! alfa=2.9
gambin(13)=2.0000 ! alfa=3.0
gambin(14)=2.1976 ! alfa=3.1
gambin(15)=2.4240 ! alfa=3.2
gambin(16)=2.6834 ! alfa=3.3
gambin(17)=2.9812 ! alfa=3.4
gambin(18)=3.3234 ! alfa=3.5
gambin(19)=3.7170 ! alfa=3.6


! --- calculate blowings snow SSA production ---
!
! maximum blowing snow particle size (diameter in micrometres)
         req=1000.0
! BS bin size interval (in micrometres)
         ddx=1.0
! total BS bin number 
         mk=int(req/ddx) 
       
       alfa=11.58*exp(-0.426*u10)+1.911    ! alfa always > 1.911
       beta=37.97*alog(0.38*u10)   
       dia=alfa*beta     !  diameter in micrometres
! get gamma values using the above look-up table (with 0.1 interval)
       do kk=1,19
           alfa1=1.75+0.1*(kk-1)
           alfa2=alfa1+0.1
       if (alfa.ge.alfa1.and.alfa.lt.alfa2) then
          gam=gambin(kk)
       else 
!          print *, 'warning: alfa is too large due to small winds (<5m/s), set to a default value',kk, alfa 
          gam=3.7  ! a default value 
       endif
       enddo
! --- end of gamma function ---

! --- snowage  ---
! in NH, snowage =1 day (Yang et al., 2024), or =3 days (Huang and Jeagle, 2016); in SH, snowage=1.5 days
!       if (ish.eq.1) then snowage=1 
!       if (ish.eq.2) then snowage=1.5 
          snowage=0.0*24  ! set =0 for a maximum SSA production (Yang et al., 2019)
         
          snowx=1./(1.038+0.03758*snowage-0.00014349*snowage**2.+1.911315E-7*snowage**3)
           if (snowx.le.0.22) snowx=0.22
! --- end of snowage  ---

! ---- 1st loop for blowing snow size bins ---
! To get a sum of "sumd" for later normalisation calculation (see Yang et al., 2019)

            sumd=0.0

         do lk=1,mk
             dd=lk*ddx  ! size of BS particles
! mechanism =1 or =2
             if (mich .eq. 1) then ! dm/dt is constant 
               dmdt=1. 
               sumd=sumd+dmdt*ddx*dd**(alfa-1)*exp(-1*dd/beta)/(beta**alfa)/gam
             else ! mich=2, dm/dt is proportional to r
               dmdt=dd
               sumd=sumd+dmdt*ddx*dd**(alfa-1)*exp(-1*dd/beta)/(beta**alfa)/gam
             endif
         enddo
! ---- end of 1st loop for blowing snow size bins ---


! --- Re-loop for blowing snow size bins ---
         do lk=1,mk

             dd=lk*ddx  ! diameter in micrometer (=1E-4 cm, or =1E-6 m)
! get water mass (in kg) of each BS particle
                wmass=4.0/3.0*PI*(dd/2.0)**3*1E-12*denice/1000.0
! now get the normalised size distribution "df" for later BS sublimation flux allocation

             if (mich.eq.1) then ! dm/dt is constant
                dmdt=1.
                df=dmdt*ddx*dd**(alfa-1)*exp(-1*dd/beta)/(beta**alfa)/gam/sumd
             else ! =2, dm/dt is proportional to r
                dmdt=dd
                df=dmdt*ddx*dd**(alfa-1+1)*exp(-1*dd/beta)/(beta**alfa)/gam/sumd
             endif

! ---- then loop snow salinity bins ----
            do iii=1,24

! Derive BS particle production flux (particles/m^2/s) from the allocated sublimation flux(mass). 
! This is the key flux to be used later for SSA flux derivation
                flux(lk,iii)=snowx*qf*df*fa(iii)/wmass           
! Get the corresponding naked SSA size (dry radius in micrometre) under a given BS size and under a given salinity
! (within the double integrations), which will be used later for re-grouping them to the 21 dry SSA bins.
                drynacl(lk,iii)=((dd/2.)**3.*sa(iii)/1000.0*denice/den)**(1./3.)

            enddo
! --- end of loop snow salinity ---

         enddo  ! end of do lk=1,mk
! --- end of Re-loop for blowing snow size bins ---
 
! Set 21 dry SSA bins (in radius from 0.01 to 10 um) for grouping (bin size can be manualluy changed if needed)
              do j=1,21
                  dryr(j)=10.0**(0.15*(j-1.0)-2.0)
              enddo

! set each bin's boundaries (upper and lower values) 
              do j=1,21

                 if (j.eq.1) then
                    dryrup(1)=0.5*(dryr(1)+dryr(2))
                    dryrlow(1)=0.5*dryr(1)
                 elseif (j.eq.21) then
                    dryrup(21)=dryr(21)+0.5*dryr(21)
                    dryrlow(21)=0.5*(dryr(20)+dryr(21))
                 else
                    dryrup(j)=(dryr(j)+dryr(j+1))/2.0
                    dryrlow(j)=(dryr(j)+dryr(j-1))/2.0
                 endif
              enddo

! SSA number emission flux (partile/m^2/s) 
         do j=1,21
             fssa(j)=0.0 ! initiate
             fssa_m(j)=0.0 
         enddo   
        print *, 'SSA ID, Diameter (µm), Lower D(µm), Upper D(µm), Nu-flux(particles/m^2/s), Mass-flux(kg/m^2/s)'

! --- loop SSA size bins  ---
         do j=1,21 ! dry SSA size bin
! --- loop blowing snow particle bins ---
            do lk =1,  mk
! --- loop snow salinity bins ----
               do iii=1,24
                    if (drynacl(lk,iii).ge.dryrlow(j).and.drynacl(lk,iii).lt.dryrup(j)) then 
                        fssa(j)=fssa(j)+flux(lk,iii) 
                        fssa_m(j)=fssa_m(j)+flux(lk,iii)*4./3.*pi*dryr(j)**3*1e-12*den/1000 ! mass flux in kg/m^2/s 
                       
                    endif
               enddo
! --- end of loop snow salinity bins ---
            enddo
! --- end of loop blowing snow particle bins ---
                        print *, 99, j, 2*dryr(j), 2*dryrlow(j), 2*dryrup(j), fssa(j), fssa_m(j) 
99    format (i6, 1x, 5(f9.4,1x))
         enddo
! --- end of loop SSA size bin for grouping ---

! Now fssa(j) contains ssa number emission flux (particles/m^2/s) from each of the 21 size bins, 
! same for fssa_m(j), which contains mass emission flux (kg/m^2/s) from each of the 21 size bins, 
! You need to use your cutoff size, e.g. for diameter=0-10 micrometres to get what you need for use. 

      end

