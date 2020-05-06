#ifdef NAVIERSTOKES
      subroutine navierstokes_init
      use param_mod
      use navierstokes_mod
      use chemkin_m
      implicit none

!     Parameters
      allocate(mu(xmn:xpn,tmn:tpn));  mu = 0.0_wp  ! viscosity
      allocate(lam(xmn:xpn,tmn:tpn)); lam = 0.0_wp ! heat diff coeff
      allocate(diff(xmn:xpn,tmn:tpn,v_Y1:v_Yn)); diff = 0.0_wp ! mass diff coeffs
      allocate(cpmix(xmn:xpn,tmn:tpn)); cpmix = 0.0_wp
      allocate(MW(xmn:xpn,tmn:tpn)); MW = 0.0_wp
      allocate(MWinv(xmn:xpn)); MWinv = 0.0_wp
      allocate(rr_r(1:nx,1:n_spec)); rr_r = 0.0_wp
      
!     RHS
      allocate(rhs_cont(1:nx)); rhs_cont = 0.0_wp
      allocate(rhs_mom(1:nx));  rhs_mom  = 0.0_wp
      allocate(rhs_ene(1:nx));  rhs_ene  = 0.0_wp
      allocate(rhs_Y(1:nx,v_Y1:v_Yn));  rhs_Y = 0.0_wp

!     Derivatives
      allocate(dUdx(1:nx)); dUdx = 0.0_wp
      allocate(d2Udx2(1:nx)); d2Udx2 = 0.0_wp
      allocate(dRHOdx(1:nx)); dRHOdx = 0.0_wp
      allocate(dEdx(1:nx)); dEdx = 0.0_wp
      allocate(dPdx(1:nx)); dPdx = 0.0_wp
      allocate(d2Pdx2(1:nx)); d2Pdx2 = 0.0_wp
      allocate(dTdx(1:nx)); dTdx = 0.0_wp
      allocate(d2Tdx2(1:nx)); d2Tdx2 = 0.0_wp
      allocate(dqdx(1:nx)); dqdx = 0.0_wp
      allocate(dmudx(1:nx)); dmudx = 0.0_wp
      allocate(dlambdadx(1:nx)); dlambdadx = 0.0_wp
      allocate(dYdx(1:nx,v_Y1:v_Yn)); dYdx = 0.0_wp
      allocate(dJdx(1:nx,v_Y1:v_Yn)); dJdx = 0.0_wp
      
      end subroutine

      subroutine navierstokes_ic
      use param_mod
      use navierstokes_mod
      use chemkin_m
      use thermchem_m
      use reference_m
      use transport_m
      implicit none

      real(wp),dimension(xmn:xpn) :: rho, vel, ene, P, Temp
      real(wp) :: min_P, max_P
      real(wp) :: min_T, max_T
      real(wp) :: min_rho, max_rho
      real(wp) :: hmix
      real(wp) :: h
      integer :: i,v,v_
      
      integer :: i_shock
      real(wp) :: RHO_mean, RHO_diff
      real(wp) :: U_mean, U_diff
      real(wp) :: P_mean, P_diff
      real(wp) :: T_mean, T_diff
#ifdef TRAVELING
      real(wp),dimension(1:nx) :: wave
#endif
      
!     Ensures mass fractions are at least 1.0e-20
      call massfrac_set_min

!     Ensures sum of mass fractions are unity
      call massfrac_unity
      
!     Copy and nondimensionalize variables -- waves already computed
      vel  = u(:,0,v_vel) / a_ref
      rho  = u(:,0,v_rho) / rho_ref
      P    = u(:,0,v_P) / p_ref     ! or P_atm?
      Temp = u(:,0,v_T) / t_ref     ! or t_o?
#ifdef DATAINPUT
      ene = u(:,0,v_ene) / a_ref / a_ref
#endif

      
!     Mixed molecular weights
      call calc_inv_avg_mol_wt_r(u(1:nx,0,v_Y1:v_Yn),MWinv(1:nx),MW(1:nx,0),1,nx) ! dim of MW = [kg/mol]
      mixMW(:,1,1) = MW(1:nx,0) ! Store current MW in thermochem_m for transport computations

!     Initializes data from input parameters
!     MUST skip this section if input data is used
#ifndef DATAINPUT

!     Compute remaining values: T, P, rho
      min_P = minval(P)
      max_P = maxval(P)
      min_T = minval(Temp)
      max_T = maxval(Temp)
      min_rho = minval(rho)
      max_rho = maxval(rho)

!     Initial values -- ideal gas law
      if (min_rho.eq.0.0_wp .and. max_rho.eq.0.0_wp) then
         rho(1:nx) = MW(1:nx,0) * P(1:nx) / Ru / Temp(1:nx)
      else if (min_P.eq.0.0_wp .and. max_P.eq.0.0_wp) then
         P(1:nx) = rho(1:nx) * Ru * Temp(1:nx) / MW(1:nx,0)
      else if (min_T.eq.0.0_wp .and. max_T.eq.0.0_wp) then
         Temp(1:nx) = P(1:nx) * MW(1:nx,0) / Ru / rho(1:nx)
      else
         print *, 'RHO, P, and T overspecified in input file'
         call abort
      end if

#ifdef TRAVELING
      ! Compute wave perturbation
      wave(1:nx) = alpha * exp(-(beta * (x-pi)/twopi)**2) / a_ref

#ifdef PERTURBSPECIES
      do v = v_Y1,v_Yn
         print *, perturb_var
         print *, var_names
         if (var_names(v).eq.perturb_var) then
            u(1:nx,0,v) = u(1:nx,0,v) + wave(1:nx)
            print *, v
         end if
      end do
      print *, maxval(u(1:nx,0,v_Y1))
      call massfrac_set_min
      call massfrac_unity
#else      
      ! Compute initial conditions
      vel(1:nx) = vel(1:nx) + wave(1:nx)
      rho(1:nx) = rho(1:nx) + wave(1:nx)
      P(1:nx) = P(1:nx) + wave(1:nx)
      Temp(1:nx) = P(1:nx) * MW(1:nx,0) / Ru / rho(1:nx)
#endif
#endif

!     Compute initial energy
      do i = 1,nx
         hmix = mixEnth(u(i,0,v_Y1:v_Yn),Temp(i))               ! nondim mix enthalpy
         ene(i) = 0.5_wp * vel(i)*vel(i) - P(i) / RHO(i) + hmix ! nondim total energy
      end do
#endif
      
!     Compute mixed heat capacity
      do i = 1,nx
         cpmix(i,0) = mixcp(u(i,0,v_Y1:v_Yn),Temp(i))   ! nondim
      end do

!     Store values -- all nondim
      u(1:nx,0,v_rho) = rho(1:nx)
      u(1:nx,0,v_vel) = vel(1:nx)
      u(1:nx,0,v_ene) = ene(1:nx)
      u(1:nx,0,v_P)   = P(1:nx)
      u(1:nx,0,v_T)   = Temp(1:nx)
      
!     Compute transport -- visc, lambda, diff
      call computeCoefficients(u(1:nx,0,v_P), u(1:nx,0,v_T),
     +     u(1:nx,0,v_Y1:v_Yn), u(1:nx,0,v_rho))

!     Store transport
      lam(1:nx,0) = lambda(:,1,1)
      mu(1:nx,0)  = viscosity(:,1,1)
      do v = v_Y1,v_Yn
         v_ = v - v_Y1 + 1
         diff(1:nx,0,v) = ds_mixavg(1:nx,1,1,v_)
      end do

!     Set boundary conditions
      call updateBCs(0)

      end subroutine
      
      subroutine navierstokes_step_rk
      use param_mod
      use navierstokes_mod
      use chem
      use deriv
      use reference_m
      use chemkin_m
      use thermchem_m
      use transport_m
      implicit none

      integer :: i,ig,v,v_, rkstep
      real(wp) :: dtime = 0.0_wp
      integer,parameter :: nl = 0 ! always 0 for RK schemes
      real(wp),dimension(xmn:xpn) :: lam_, mu_, MW_
      real(wp),dimension(xmn:xpn,v_Y1:v_Yn) :: diff_
      real(wp) :: dt_           ! nondim dt
      
!     Variables for diffusion/chemistry stuff
      real(wp) :: di, ddidx, Vi, dVidx, Xi, dXidx, d2Xidx2, d2Ydx2, ddiffdx, J
      real(wp) :: dMWdx, d2MWdx2

!     Temporary variables for easy reading
      real(wp) :: Y_, rho_, D_, U_, P_

!     Variables for storage
      real(wp),dimension(xmn:xpn) :: rho, vel, ene, P, Temp
      real(wp),dimension(xmn:xpn,v_Y1:v_Yn) :: Y
      
!     Store data to local variables.
      rho  = u(:,0,v_rho)
      vel  = u(:,0,v_vel)
      ene  = u(:,0,v_ene)
      P    = u(:,0,v_P)
      Temp = u(:,0,v_T)
      lam_ = lam(:,0)
      mu_  = mu(:,0)
      MW_  = MW(:,0)
      do v = v_Y1,v_Yn
         Y(:,v) = u(:,0,v)
         diff_(:,v) = diff(:,0,v)
      end do

!     RK Time Integration
!     -------------------
      rk_loop : do rkstep = 1,rk%nsteps
      rk%step = rkstep

!     Compute derivatives
      call computeDerivatives(nl)
      
!     Compute source terms
      call computeSourceTerms(nl)

!     Nondim timestep
      dt_ = dt / time_ref

!     Rescale RHS
      rhs_cont = rhs_cont * rk%k(rk%step,1)
      rhs_mom = rhs_mom * rk%k(rk%step,1)
      rhs_ene = rhs_ene * rk%k(rk%step,1)
      rhs_Y = rhs_Y * rk%k(rk%step,1)

!     Update RHS
      call computeConRHS(nl,dt_)
      call computeMomRHS(nl,dt_)
      call computeEneRHS(nl,dt_)
      call computeYspRHS(nl,dt_)
      
!     Update variables
      u(1:nx,nl,v_rho) = u(1:nx,nl,v_rho) + rhs_cont(1:nx) * rk%k(rk%step,2)
      u(1:nx,nl,v_vel) = u(1:nx,nl,v_vel) + rhs_mom(1:nx) * rk%k(rk%step,2)
      u(1:nx,nl,v_ene) = u(1:nx,nl,v_ene) + rhs_ene(1:nx) * rk%k(rk%step,2)
      u(1:nx,nl,v_Y1:v_Yn) = u(1:nx,nl,v_Y1:v_Yn) + rhs_Y(1:nx,v_Y1:v_Yn) * rk%k(rk%step,2)
      call massfrac_set_min
      call massfrac_unity
      
!     Update properties
      call updatePropertiesWithS3D(nl)
      
!     BCs
      call updateBCs(nl)

!     Update time
      dtime = rk%k(rk%step,1)*dtime + dt
      t = t + rk%k(rk%step,2)*dtime

      enddo rk_loop

!     Store timesteps
!     ---------------

!     Timestep n+1
      do i = xmn,xpn
         u(i,1,v_rho) = u(i,0,v_rho)
         u(i,1,v_vel) = u(i,0,v_vel)
         u(i,1,v_ene) = u(i,0,v_ene)
         u(i,1,v_P) = u(i,0,v_P)
         u(i,1,v_T) = u(i,0,v_T)
         lam(i,1) = lam(i,0)
         mu(i,1) = mu(i,0)
         MW(i,1) = MW(i,0)
      end do
      do v = v_Y1,v_Yn
         do i = xmn,xpn
            u(i,1,v) = u(i,0,v)
            diff(i,1,v) = diff(i,0,v)
         end do
      end do

!     Timestep n
      u(:,0,v_rho) = rho
      u(:,0,v_vel) = vel
      u(:,0,v_ene) = ene
      u(:,0,v_P) = P
      u(:,0,v_T) = Temp
      lam(:,0) = lam_
      mu(:,0) = mu_
      MW(:,0) = MW_
      do v = v_Y1,v_Yn
         u(:,0,v) = Y(:,v)
         diff(:,0,v) = diff_(:,v)
      end do

!     Update timesteps
      do ig = tmn,0
         u(:,ig,v_rho) = u(:,ig+1,v_rho)
         u(:,ig,v_vel) = u(:,ig+1,v_vel)
         u(:,ig,v_ene) = u(:,ig+1,v_ene)
         u(:,ig,v_P) = u(:,ig+1,v_P)
         u(:,ig,v_T) = u(:,ig+1,v_T)
         lam(:,ig) = lam(:,ig+1)
         mu(:,ig) = mu(:,ig+1)
         MW(:,ig) = MW(:,ig+1)
         do v = v_Y1,v_Yn
            u(:,ig,v) = u(:,ig+1,v)
            diff(:,ig,v) = diff(:,ig+1,v)
         end do
      end do

      end subroutine

      subroutine navierstokes_step_ab
#ifdef AB
      use param_mod
      use navierstokes_mod
      use chem
      use deriv
      use reference_m
      use chemkin_m
      use thermchem_m
      use transport_m
      implicit none

      integer :: i,ig,v,v_, abstep
      real(wp) :: dtime = 0.0_wp
      integer :: nl
      real(wp),dimension(xmn:xpn) :: rho, vel, ene, P, Temp
      real(wp),dimension(xmn:xpn) :: lam_, mu_, MW_
      real(wp),dimension(xmn:xpn,v_Y1:v_Yn) :: Y, diff_
      real(wp) :: dt_           ! nondim dt
      
!     Variables for diffusion/chemistry stuff
      real(wp) :: di, ddidx, Vi, dVidx, Xi, dXidx, d2Xidx2, d2Ydx2, ddiffdx, J
      real(wp) :: dMWdx, d2MWdx2

!     Temporary variables for easy reading
      real(wp) :: Y_, rho_, D_, U_, P_

!     Reset RHS to 0
      rhs_cont = R0
      rhs_mom  = R0
      rhs_ene  = R0
      rhs_Y    = R0

!     AB step
      ab_loop : do abstep = 1,ab%nsteps
      ab%step = abstep

!     AB scheme setup
      nl  = 1 - abstep
      dt_ = dt * ab%k(abstep) / time_ref
      
!     Compute derivatives
      call computeDerivatives(nl)

!     Compute source terms
      call computeSourceTerms(nl)

!     Update RHS
      call computeConRHS(nl,dt_)
      call computeMomRHS(nl,dt_)
      call computeEneRHS(nl,dt_)
      call computeYspRHS(nl,dt_)
      
      enddo ab_loop
      
!     Update variables
!     ----------------

      t = t + dt
      
!     Current timestep
      u(1:nx,1,v_rho) = u(1:nx,0,v_rho) + rhs_cont(1:nx)
      u(1:nx,1,v_vel) = u(1:nx,0,v_vel) + rhs_mom(1:nx)
      u(1:nx,1,v_ene) = u(1:nx,0,v_ene) + rhs_ene(1:nx)
      u(1:nx,1,v_Y1:v_Yn) = u(1:nx,0,v_Y1:v_Yn) + rhs_Y(1:nx,v_Y1:v_Yn)
      call massfrac_set_min
      call massfrac_unity
      
!     Update timesteps
      do ig = tmn,0
         u(:,ig,v_rho) = u(:,ig+1,v_rho)
         u(:,ig,v_vel) = u(:,ig+1,v_vel)
         u(:,ig,v_ene) = u(:,ig+1,v_ene)
         u(:,ig,v_P) = u(:,ig+1,v_P)
         u(:,ig,v_T) = u(:,ig+1,v_T)
         lam(:,ig) = lam(:,ig+1)
         mu(:,ig) = mu(:,ig+1)
         MW(:,ig) = MW(:,ig+1)
         do v = v_Y1,v_Yn
            u(:,ig,v) = u(:,ig+1,v)
         end do
      end do

!     Update properties
      nl = 0
      call updatePropertiesWithS3D(nl)
      
!     BCs
      call updateBCs(nl)

      if (maxval(u(:,0,v_T)).gt.3000/T_ref) print *, minval(u(:,0,v_T))*T_ref, maxval(u(:,0,v_T))*T_ref
#endif
      end subroutine

      subroutine navierstokes_dimensional
      use navierstokes_mod
      use reference_m
      implicit none

      u(:,0,v_rho) = u(:,0,v_rho) * rho_ref
      u(:,0,v_vel) = u(:,0,v_vel) * a_ref
      u(:,0,v_ene) = u(:,0,v_ene) * a_ref * a_ref
      u(:,0,v_T) = u(:,0,v_T) * t_ref
      u(:,0,v_P) = u(:,0,v_P) * p_ref
      
      end subroutine navierstokes_dimensional

      subroutine navierstokes_nondimensional
      use navierstokes_mod
      use reference_m
      implicit none
      
      u(:,0,v_rho) = u(:,0,v_rho) / rho_ref
      u(:,0,v_vel) = u(:,0,v_vel) / a_ref
      u(:,0,v_ene) = u(:,0,v_ene) / a_ref / a_ref
      u(:,0,v_T) = u(:,0,v_T) / t_ref
      u(:,0,v_P) = u(:,0,v_P) / p_ref
      
      end subroutine navierstokes_nondimensional

      subroutine massfrac_set_min
      use param_mod
      use navierstokes_mod
      implicit none
      integer :: i,v
      real(wp) :: min_ = 1.0e-60
      real(wp) :: test

      do v = v_Y1,v_Yn
         test = minval(u(:,0,v))
         if (test.lt.min_) then
            do i = 1,nx
               if (u(i,0,v).lt.min_) u(i,0,v) = min_
            end do
         end if
      end do
      
      end subroutine massfrac_set_min

      subroutine massfrac_unity
      use param_mod
      use navierstokes_mod
      implicit none
      integer :: i,v

      u(:,0,v_Yn) = 1.0_wp
      do v = v_Y1,v_Yn-1
         u(:,0,v_Yn) = u(:,0,v_Yn) - u(:,0,v)
      end do

      end subroutine massfrac_unity

      subroutine computeDerivatives(nl)
      use param_mod
      use navierstokes_mod
      use chem
      use deriv
      use reference_m
      use chemkin_m
      use thermchem_m
      use transport_m
      implicit none
      integer,intent(in) :: nl
      integer :: i, v, v_
      real(wp) :: Y_, rho_, D_, Xi, Vi, J
      real(wp) :: ddiffdx, d2Ydx2, dMWdx, d2MWdx2, dXidx, d2Xidx2, dVidx
      
      do i = 1,nx
         ! Gradients
         dUdx(i)   = ddx(u(:,:,v_vel),i,nl)
         dRHOdx(i) = ddx(u(:,:,v_rho),i,nl)
         dPdx(i)   = ddx(u(:,:,v_P),i,nl)
         dTdx(i)   = ddx(u(:,:,v_T),i,nl)
         dEdx(i)   = ddx(u(:,:,v_ene),i,nl)
         dmudx(i)  = ddx(mu,i,nl)
         dlambdadx(i) = ddx(lam,i,nl)
         
         ! Curvatures
         d2Udx2(i) = d2dx2(u(:,:,v_vel),i,nl)
         d2Pdx2(i) = d2dx2(u(:,:,v_P),i,nl)
         d2Tdx2(i) = d2dx2(u(:,:,v_T),i,nl)

         ! Nonreacting derivative of heat flux
         dqdx(i) = - dlambdadx(i) * dTdx(i) - lam(i,nl) * d2Tdx2(i)
      end do

      ! Compute chemistry derivatives
      do v = v_Y1,v_Yn
         v_ = v - v_Y1 + 1 ! 1 to n_spec
         do i = 1,nx
            ! Temporary local variables for easy reading
            Y_ = u(i,nl,v)
            rho_ = u(i,nl,v_rho)
            D_ = diff(i,nl,v)
            ddiffdx = ddx(diff(:,:,v),i,nl)
            dYdx(i,v) = ddx(u(:,:,v),i,nl) ! Stored for later use
            d2Ydx2 = d2dx2(u(:,:,v),i,nl)
            dMWdx = ddx(MW,i,nl)

            ! Mass fraction and derivatives
            Xi = Y_ * MW(i,nl) * molwt_c(v_)
            dXidx = molwt_c(v_) * (MW(i,nl) * dYdx(i,v) + Y_ * dMWdx)
            d2Xidx2 = molwt_c(v_) * (MW(i,nl) * d2Ydx2
     +           + Y_ * d2MWdx2 + 2 * dMWdx * dYdx(i,v))

            ! Species mass diffusion velocity and derivative -- neglects both Soret and diffusion
            Vi = - D_ * dXidx / Xi
            dVidx = D_ * dXidx * dXidx / Xi / Xi 
     +           - D_  * d2Xidx2 / Xi
     +           - ddiffdx * dXidx / Xi
            
            ! mass diffusion flux
            J = rho_ * Y_ * Vi
            dJdx(i,v) = dRhodx(i) * Y_ * Vi
     +           + rho_ * dYdx(i,v) * Vi
     +           + rho_ * Y_ * dVidx

            ! Chemistry contribution to heat flux
            dqdx(i) = dqdx(i) + specEnth(v_,u(i,nl,v_T)) * dJdx(i,v)
     +           + specCp(v_,u(i,nl,v_T)) * dTdx(i) * J
         end do
      end do
      
      end subroutine computeDerivatives

!     Computes chemical source terms for species transport equations
      subroutine computeSourceTerms(nl)
      use param_mod
      use navierstokes_mod
      use chem
      use deriv
      use reference_m
      use chemkin_m
      use thermchem_m
      use transport_m
      implicit none
      integer,intent(in) :: nl
      integer :: i
      
!     Compute source terms
      do i = 1,nx
         call getrates(u(i,nl,v_P)*p_ref_cgs, u(i,nl,v_T)*t_ref,
     +        u(i,nl,v_Y1:v_Yn), ickwrk, rckwrk, rr_r(i,:))
         ! [1]   = [mol/cm^3/s]*[1000 g/kg]*[kg/mol] / [g/cm^3/s]
         rr_r(i,:) = rr_r(i,:) * 1.0e3_wp * molwt(:) / rr_ref_cgs
      end do

      end subroutine computeSourceTerms

!     Continuity equation
      subroutine computeConRHS(nl,dt_)
      use param_mod
      use navierstokes_mod
      use chem
      use deriv
      use reference_m
      use chemkin_m
      use thermchem_m
      use transport_m
      implicit none
      integer,intent(in) :: nl
      real(wp),intent(in) :: dt_
      integer :: i
      real(wp) :: rho_, U_
      
      do i = 1,nx
         !Temporary variables
         rho_ = u(i,nl,v_rho)
         U_ = u(i,nl,v_vel)

         ! Compute RHS
         rhs_cont(i) = rhs_cont(i) - dt_ * rho_ * dUdx(i)
         rhs_cont(i) = rhs_cont(i) - dt_ * U_ * dRHOdx(i)
      end do

      end subroutine computeConRHS

!     Momentum equations
      subroutine computeMomRHS(nl,dt_)
      use param_mod
      use navierstokes_mod
      use chem
      use deriv
      use reference_m
      use chemkin_m
      use thermchem_m
      use transport_m
      implicit none
      integer,intent(in) :: nl
      real(wp),intent(in) :: dt_
      integer :: i
      real(wp) :: rho_, U_

      do i = 1,nx
         !Temporary variables
         rho_ = u(i,nl,v_rho)
         U_ = u(i,nl,v_vel)
         
         ! Compute RHS
         rhs_mom(i) = rhs_mom(i) - dt_ * U_ * dUdx(i)   ! Convection
         rhs_mom(i) = rhs_mom(i) - dt_ * dPdx(i) / rho_ ! Pressure force
         rhs_mom(i) = rhs_mom(i) + dt_ * mu(i,nl) * d2Udx2(i) / rho_ ! diff 1
         rhs_mom(i) = rhs_mom(i) + dt_ * dmudx(i) * dUdx(i) / rho_ ! diff 2
      end do

      end subroutine computeMomRHS

!     Energy equation
      subroutine computeEneRHS(nl,dt_)
      use param_mod
      use navierstokes_mod
      use chem
      use deriv
      use reference_m
      use chemkin_m
      use thermchem_m
      use transport_m
      implicit none
      integer,intent(in) :: nl
      real(wp),intent(in) :: dt_
      integer :: i
      real(wp) :: rho_, U_, P_

      do i = 1,nx
         !Temporary variables
         rho_ = u(i,nl,v_rho)
         U_ = u(i,nl,v_vel)
         P_ = u(i,nl,v_P)

         ! Compute RHS
         rhs_ene(i) = rhs_ene(i) - dt_ * U_ * dEdx(i) ! convection
         rhs_ene(i) = rhs_ene(i) - dt_ * U_ * dPdx(i) / rho_                ! pressure 1
         rhs_ene(i) = rhs_ene(i) - dt_ * P_ * dUdx(i) / rho_                ! pressure 2
         rhs_ene(i) = rhs_ene(i) + dt_ * U_ * dmudx(i) * dUdx(i) / rho_ ! diff 1
         rhs_ene(i) = rhs_ene(i) + dt_ * mu(i,nl) * dUdx(i) * dUdx(i) / rho_ ! diff 2
         rhs_ene(i) = rhs_ene(i) + dt_ * mu(i,nl) * U_ * d2Udx2(i) / rho_ ! diff 3
         rhs_ene(i) = rhs_ene(i) - dt_ * dqdx(i) / rho_
      end do

      end subroutine computeEneRHS

!     Chemical species equations
      subroutine computeYspRHS(nl,dt_)
      use param_mod
      use navierstokes_mod
      use chem
      use deriv
      use reference_m
      use chemkin_m
      use thermchem_m
      use transport_m
      implicit none
      integer,intent(in) :: nl
      real(wp),intent(in) :: dt_
      integer :: i
      real(wp) :: rho_, U_
      integer :: v, v_
      
      do v = v_Y1,v_Yn
         v_ = v - v_Y1 + 1
         do i = 1,nx
            ! Temporary variables
            rho_ = u(i,nl,v_rho)
            U_ = u(i,nl,v_vel)

            ! Compute RHS
            rhs_Y(i,v) = rhs_Y(i,v) - dt_ * U_ * dYdx(i,v) ! Convection
            rhs_Y(i,v) = rhs_Y(i,v) - dt_ * dJdx(i,v) / rho_ ! Diffusion
            rhs_Y(i,v) = rhs_Y(i,v) + dt_ * rr_r(i,v_) / rho_ ! Source
         end do
      end do

      end subroutine computeYspRHS

!     Update transport properties
      subroutine updatePropertiesWithS3D(nl)
      use param_mod
      use navierstokes_mod
      use chem
      use deriv
      use reference_m
      use chemkin_m
      use thermchem_m
      use transport_m
      implicit none
      integer,intent(in) :: nl
      integer :: v, v_

      call calc_inv_avg_mol_wt_r(u(1:nx,nl,v_Y1:v_Yn),MWinv(1:nx),MW(1:nx,nl),1,nx)
      mixMW(:,1,1) = MW(1:nx,nl) ! Store current MW in thermochem_m for transport computations
      call calc_temp_r(u(1,nl,v_T),u(1,nl,v_ene),u(1,nl,v_vel),
     +     u(1:nx,nl,v_Y1:v_Yn),cpmix(1:nx,nl),MWinv(1:nx),1,nx)
      call calc_press_r(u(1,nl,v_P),u(1,nl,v_rho),u(1,nl,v_T),MWinv(1),1,nx)
      call computeCoefficients(u(1:nx,nl,v_P), u(1:nx,nl,v_T),
     +     u(1:nx,nl,v_Y1:v_Yn), u(1:nx,nl,v_rho)*rho_ref_cgs)

      lam(1:nx,nl) = lambda(:,1,1)
      mu(1:nx,nl)  = viscosity(:,1,1)
      do v = v_Y1,v_Yn
         v_ = v - v_Y1 + 1
         diff(1:nx,nl,v) = ds_mixavg(1:nx,1,1,v_)
      end do

      end subroutine updatePropertiesWithS3D

!     Update boundary conditions
      subroutine updateBCs(nl)
      use param_mod
      use navierstokes_mod
      use chem
      use deriv
      use reference_m
      use chemkin_m
      use thermchem_m
      use transport_m
      implicit none
      integer,intent(in) :: nl
      integer :: v

      call bc_array(u(:,:,v_rho),nl)
      call bc_array(u(:,:,v_vel),nl)
      call bc_array(u(:,:,v_ene),nl)
      call bc_array(u(:,:,v_P),nl)
      call bc_array(u(:,:,v_T),nl)
      call bc_array(lam,nl)
      call bc_array(mu,nl)
      call bc_array(MW,nl)
      do v = v_Y1,v_Yn
         call bc_array(u(:,:,v),nl)
         call bc_array(diff(:,:,v),nl)
      end do      
      
      end subroutine updateBCs
      
#endif
