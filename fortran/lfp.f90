!> LFP half-cell model: Li foil | separator | porous LiFePO4 cathode.
!>
!> A self-contained program: edit the input file, run, and it writes Time_Voltage.txt.
!> The equations are described in docs/model.md; parameters in docs/parameters.md.
!> The linear block system of each time step is solved with bandsolver (Newman's BAND).
!>
!>   usage:  lfp_f [input.nml]        (default input file: lfp.nml)
!>
!> mode = 'faithful' reproduces the original research code, including the defects
!> listed in docs/deviations.md; mode = 'corrected' applies the fixes.
!>
!> SPDX-License-Identifier: BSD-3-Clause
program lfp
    use, intrinsic :: iso_fortran_env, only: dp => real64, sp => real32, error_unit
    use, intrinsic :: ieee_arithmetic, only: ieee_is_nan, ieee_value, ieee_quiet_nan
    use bandsolver_kernel, only: band_solve, BAND_OK, PIVOT_LEGACY
    implicit none

    integer, parameter :: NV = 4                  ! unknowns per node
    integer, parameter :: IC = 1, IP1 = 2, IP2 = 3, ICS = 4

    ! ---------------- parameters (defaults = original research code) ----------------
    real(dp) :: L_cath = 24.0e-4_dp, L_sep = 25.0e-4_dp
    integer  :: nj = 101, sep_node = 22
    real(dp) :: eps = 0.5_dp, eps_AM = 0.8_dp, eps_sep = 0.39_dp, tau_sep = 4.0_dp, bruggeman = -0.5_dp
    real(dp) :: D = 2.0e-6_dp, t_plus = 0.25_dp, c_bulk = 1.0e-3_dp, z_plus = 1.0_dp, z_minus = -1.0_dp
    real(dp) :: sigma = 3.0e-3_dp, M = 125.759_dp, rho = 3.6_dp, Q_th = 0.170_dp, R_p = 200.0e-7_dp
    real(dp) :: k_rxn = -1.0_dp                    ! < 0: use the default for the chosen mode
    real(dp) :: alpha_a = 0.5_dp, alpha_c = 0.5_dp, k_Li = 1.0e-6_dp, c_Li_ref = 1.0e-3_dp
    real(dp) :: R = 8.314_dp, T = 298.0_dp, F = 96485.0_dp
    real(dp) :: C_rate = 1.0_dp, phi1_init = 3.6_dp, phi2_init = 0.0_dp, cs_init = 1.0e-5_dp
    real(dp) :: t_max = 36000.0_dp
    integer  :: n_steps = 36000
    real(dp) :: V_min = 2.5_dp, V_max = 4.2_dp        ! corrected-mode cutoffs (D-3)
    real(dp) :: fd_step = 1.0e-6_dp
    real(dp) :: newton_tol = 1.0e-10_dp                ! corrected-mode Newton tolerance (D-7)
    integer  :: newton_max_iter = 25
    character(len=16)  :: mode = 'faithful'
    character(len=256) :: file = 'Time_Voltage.txt'

    namelist /cell/ L_cath, L_sep, nj, sep_node, eps, eps_AM, eps_sep, tau_sep, bruggeman
    namelist /electrolyte/ D, t_plus, c_bulk, z_plus, z_minus
    namelist /active/ sigma, M, rho, Q_th, R_p, k_rxn, alpha_a, alpha_c, k_Li, c_Li_ref
    namelist /constants/ R, T, F
    namelist /operation/ C_rate, phi1_init, phi2_init, cs_init, t_max, n_steps, V_min, V_max
    namelist /numerics/ fd_step, newton_tol, newton_max_iter, mode
    namelist /output/ file

    ! ---------------- derived quantities and mode switches ----------------
    logical  :: faithful
    real(dp) :: lit36, ocp_c(7), spec_a, tortuosity, i_spec, i_app, eps_sep_face, phi1_sign
    logical  :: full_current
    real(dp) :: dcat_s, dan_s, ucat_s, uan_s, dcat_c, dan_c, ucat_c, uan_c

    ! ---------------- mesh and state ----------------
    integer :: s                                  ! interface node index (1-based = sep_node)
    real(dp), allocatable :: dx(:), aW(:), aE(:), bW(:), bE(:)
    real(dp), allocatable :: c(:,:), dc(:,:), A(:,:,:), B(:,:,:), Dm(:,:,:), G(:,:)

    ! ---------------- time loop ----------------
    integer  :: it, status, ounit, nsolve
    real(dp) :: time, dt, mAhg, write_every, last_write, v, h
    logical  :: want_row, header, ok, stopped
    character(len=1) :: state
    character(len=64) :: exit_reason
    character(len=256) :: input_file

    ! ---------------- read input ----------------
    input_file = 'lfp.nml'
    if (command_argument_count() >= 1) call get_command_argument(1, input_file)
    call read_input(trim(input_file))
    call setup()

    open(newunit=ounit, file=trim(file), status='replace', action='write')

    allocate(c(NV,nj), dc(NV,nj), A(NV,NV,nj), B(NV,NV,nj), Dm(NV,NV,nj), G(NV,nj))
    c(IC,:) = c_bulk
    c(IP1,:) = phi1_init
    c(IP2,:) = phi2_init
    c(ICS,:) = cs_init
    dc = 0.0_dp
    nsolve = 0
    time = 0.0_dp
    mAhg = 0.0_dp
    dt = t_max/real(n_steps, dp)
    last_write = 0.0_dp
    write_every = t_max/n_steps/200
    state = 'D'
    if (C_rate < 0) state = 'C'
    exit_reason = 'max_steps'

    do it = 1, n_steps
        want_row = .false.
        header = .false.
        if (it == 1) then
            want_row = .true.
            header = .true.
        else if ((time - last_write)/3600 >= write_every) then
            want_row = .true.
            if (faithful) then
                last_write = real(int(time - dt), dp)   ! an integer in the original (D-4)
            else
                last_write = time
            end if
        else if (it >= n_steps) then
            want_row = .true.
        else if (c(IP1,nj) >= 99.0_dp .and. state == 'C') then
            call write_row(.false.)
            exit_reason = 'end_of_charge'
            exit
        else if (ieee_is_nan(dc(IC,1))) then
            call write_row(.false.)
            exit_reason = 'nan'
            exit
        else if (time >= 99.0_dp*3600.0_dp) then
            call write_row(.false.)
            exit_reason = 'max_time'
            exit
        end if
        if (.not. faithful) then
            v = cell_voltage()
            if (v <= V_min .or. v >= V_max) then
                call write_row(header)
                exit_reason = merge('cutoff_low ', 'cutoff_high', v <= V_min)
                exit
            end if
        end if
        if (want_row) call write_row(header)

        if (state == 'D') then
            mAhg = mAhg + 1000.0_dp*i_spec*dt/3600.0_dp
        else if (state == 'C') then
            mAhg = mAhg - 1000.0_dp*i_spec*dt/3600.0_dp
        end if

        if (faithful) then
            call assemble(dt)
            call band_solve(NV, nj, A, B, Dm, G, dc, status, pivot=PIVOT_LEGACY)
            if (status /= BAND_OK) dc = ieee_value(1.0_dp, ieee_quiet_nan)
            c = c + dc
        else
            call advance(dt, h, stopped, ok)
            if (.not. ok) then
                exit_reason = 'solver_fail'
                exit
            end if
            if (stopped) then
                if (state == 'D') then
                    mAhg = mAhg - 1000.0_dp*i_spec*(dt - h)/3600.0_dp
                else
                    mAhg = mAhg + 1000.0_dp*i_spec*(dt - h)/3600.0_dp
                end if
                time = time + h
                v = cell_voltage()
                call write_row(.false.)
                exit_reason = merge('cutoff_low ', 'cutoff_high', v <= V_min)
                nsolve = it
                exit
            end if
        end if
        nsolve = it

        if (state == 'R') then
            dt = dt*1.0001_dp
        else
            dt = t_max/real(n_steps, dp)
        end if
        time = time + dt
    end do

    close(ounit)
    write(*,'(A,G0,A,I0,A)') trim(mode)//' run, C-rate ', C_rate, ': exit '//trim(exit_reason)//' after ', &
        nsolve, ' steps; wrote '//trim(file)

contains

    ! =============================== input ===============================
    subroutine read_input(path)
        character(len=*), intent(in) :: path
        integer :: u, ios
        logical :: exists
        inquire(file=path, exist=exists)
        if (.not. exists) then
            write(error_unit,'(A)') 'input file not found: '//path
            error stop 2
        end if
        open(newunit=u, file=path, status='old', action='read')
        rewind(u); read(u, nml=cell, iostat=ios);        call check(ios, 'cell')
        rewind(u); read(u, nml=electrolyte, iostat=ios); call check(ios, 'electrolyte')
        rewind(u); read(u, nml=active, iostat=ios);      call check(ios, 'active')
        rewind(u); read(u, nml=constants, iostat=ios);   call check(ios, 'constants')
        rewind(u); read(u, nml=operation, iostat=ios);   call check(ios, 'operation')
        rewind(u); read(u, nml=numerics, iostat=ios);    call check(ios, 'numerics')
        rewind(u); read(u, nml=output, iostat=ios);      call check(ios, 'output')
        close(u)
    end subroutine read_input

    subroutine check(ios, group)
        integer, intent(in) :: ios
        character(len=*), intent(in) :: group
        if (ios > 0) then
            write(error_unit,'(A)') 'error reading namelist group &'//group
            error stop 2
        end if
        ! ios < 0: group absent, keep defaults
    end subroutine check

    pure real(dp) function r32(x)
        !! Round to single precision (reproduces the original's single-precision literals).
        real(dp), intent(in) :: x
        r32 = real(real(x, sp), dp)
    end function r32

    ! =============================== setup ===============================
    subroutine setup()
        real(dp), parameter :: OCP_FIT(7) = [3.114559_dp, 4.438792_dp, 71.7352_dp, 70.85337_dp, &
                                             4.240252_dp, 68.5605_dp, 67.730082_dp]
        ! 10**0.966 as the original evaluated it: single precision, constant-folded
        real(dp), parameter :: K_FAITHFUL = 1.0e-8_dp*real(10.0_sp**0.966_sp, dp)
        real(dp) :: h_sep, h_cat, t_an, d_cat, d_an, u_cat, u_an
        integer :: j

        select case (trim(mode))
        case ('faithful')
            faithful = .true.
        case ('corrected')
            faithful = .false.
        case default
            write(error_unit,'(A)') 'mode must be ''faithful'' or ''corrected'''
            error stop 2
        end select

        if (faithful) then
            R = r32(R); c_bulk = r32(c_bulk); Q_th = r32(Q_th); M = r32(M); rho = r32(rho)
            phi1_init = r32(phi1_init); eps_sep = r32(eps_sep); eps_AM = r32(eps_AM); C_rate = r32(C_rate)
            if (k_rxn < 0) k_rxn = K_FAITHFUL
            lit36 = r32(3.6_dp)
            ocp_c = [(r32(OCP_FIT(j)), j = 1, 7)]
            eps_sep_face = eps          ! D-2
            phi1_sign = -1.0_dp         ! D-1
            full_current = .false.      ! D-11
        else
            if (k_rxn < 0) k_rxn = 1.0e-8_dp*10.0_dp**0.966_dp
            lit36 = 3.6_dp
            ocp_c = OCP_FIT
            eps_sep_face = eps_sep
            phi1_sign = 1.0_dp
            full_current = .true.
        end if

        spec_a = 3*eps_AM/R_p
        tortuosity = eps**bruggeman
        i_spec = Q_th*C_rate
        i_app = i_spec*L_cath*eps_AM*rho

        ! ion diffusivities and mobilities, then effective values per region
        t_an = 1.0_dp - t_plus
        d_cat = D*(1.0_dp + (t_an/t_plus))/(2.0_dp*t_an/t_plus)
        d_an = d_cat*t_an/t_plus
        u_cat = d_cat/(R*T)
        u_an = d_an/(R*T)
        dcat_s = d_cat/tau_sep;    dan_s = d_an/tau_sep;    ucat_s = u_cat/tau_sep;    uan_s = u_an/tau_sep
        dcat_c = d_cat/tortuosity; dan_c = d_an/tortuosity; ucat_c = u_cat/tortuosity; uan_c = u_an/tortuosity

        s = sep_node
        allocate(dx(nj), aW(nj), aE(nj), bW(nj), bE(nj))
        h_sep = L_sep/real(sep_node - 2, dp)
        h_cat = L_cath/real(nj - sep_node - 1, dp)
        dx = 0.0_dp
        dx(2:s-1) = h_sep
        dx(s+1:nj-1) = h_cat
        aW = 0.0_dp; bW = 0.0_dp; aE = 0.0_dp; bE = 0.0_dp
        do j = 2, nj
            aW(j) = dx(j-1)/(dx(j-1) + dx(j))
            bW(j) = 2.0_dp/(dx(j-1) + dx(j))
        end do
        do j = 1, nj - 1
            aE(j) = dx(j)/(dx(j+1) + dx(j))
            bE(j) = 2.0_dp/(dx(j) + dx(j+1))
        end do
    end subroutine setup

    ! =============================== kinetics ===============================
    real(dp) function cs_max()
        cs_max = (rho/M)*M*Q_th*1000.0_dp*lit36/F
    end function cs_max

    real(dp) function ocp(cs)
        real(dp), intent(in) :: cs
        real(dp) :: th
        th = (cs/(rho/M))/(M*Q_th*1000.0_dp*lit36/F)
        ocp = ocp_c(1) + ocp_c(2)*atan(-(ocp_c(3)*th) + ocp_c(4)) - ocp_c(5)*atan(-(ocp_c(6)*th) + ocp_c(7))
    end function ocp

    real(dp) function rate(cc, cs, p1, p2)
        !! Butler-Volmer current density per interfacial area [A/cm2], anodic positive.
        real(dp), intent(in) :: cc, cs, p1, p2
        real(dp) :: eta, i0
        eta = p1 - p2 - ocp(cs)
        i0 = F*k_rxn*(cc**alpha_a)*((cs_max() - cs)**alpha_a)*(cs**alpha_c)
        if (faithful) i0 = r32(i0)     ! D-4: implicitly single precision in the original
        rate = i0*(exp(alpha_a*F*eta/(R*T)) - exp(-(alpha_c*F*eta/(R*T))))
    end function rate

    subroutine rate_derivs(cc, cs, p1, p2, i, di)
        !! Rate and finite-difference derivatives w.r.t. (c, phi1, phi2, cs) (D-6).
        real(dp), intent(in) :: cc, cs, p1, p2
        real(dp), intent(out) :: i, di(NV)
        real(dp) :: h
        if (.not. faithful) then
            call rate_derivs_exact(cc, cs, p1, p2, i, di)
            return
        end if
        h = fd_step
        i = rate(cc, cs, p1, p2)
        if (cc <= h) then
            di(IC) = (rate(cc + h, cs, p1, p2) - i)/h
        else
            di(IC) = (rate(cc + h, cs, p1, p2) - rate(cc - h, cs, p1, p2))/(2.0_dp*h)
        end if
        if (cs <= h) then
            di(ICS) = (rate(cc, cs + h, p1, p2) - i)/h
        else
            di(ICS) = (rate(cc, cs + h, p1, p2) - rate(cc, cs - h, p1, p2))/(2.0_dp*h)
        end if
        di(IP1) = (rate(cc, cs, p1 + h, p2) - rate(cc, cs, p1 - h, p2))/(2.0_dp*h)
        di(IP2) = (rate(cc, cs, p1, p2 + h) - rate(cc, cs, p1, p2 - h))/(2.0_dp*h)
    end subroutine rate_derivs

    real(dp) function ocp_slope(cs)
        !! dU/dcs [V cm3/mol]
        real(dp), intent(in) :: cs
        real(dp) :: th, dth, du
        th = (cs/(rho/M))/(M*Q_th*1000.0_dp*lit36/F)
        dth = 1.0_dp/((rho/M)*(M*Q_th*1000.0_dp*lit36/F))
        du = ocp_c(2)*(-ocp_c(3))/(1.0_dp + (-(ocp_c(3)*th) + ocp_c(4))**2) &
           - ocp_c(5)*(-ocp_c(6))/(1.0_dp + (-(ocp_c(6)*th) + ocp_c(7))**2)
        ocp_slope = du*dth
    end function ocp_slope

    subroutine rate_derivs_exact(cc, cs, p1, p2, i, di)
        !! Rate and exact derivatives w.r.t. (c, phi1, phi2, cs); corrected mode (fixes D-6).
        real(dp), intent(in) :: cc, cs, p1, p2
        real(dp), intent(out) :: i, di(NV)
        real(dp) :: rt, aa, bb, eta, i0, ea, ec, di_deta
        rt = R*T
        aa = alpha_a*F/rt
        bb = alpha_c*F/rt
        eta = p1 - p2 - ocp(cs)
        i0 = F*k_rxn*(cc**alpha_a)*((cs_max() - cs)**alpha_a)*(cs**alpha_c)
        ea = exp(aa*eta)
        ec = exp(-bb*eta)
        i = i0*(ea - ec)
        di_deta = i0*(aa*ea + bb*ec)
        di(IC) = alpha_a*i/cc
        di(ICS) = i*(-alpha_a/(cs_max() - cs) + alpha_c/cs) - di_deta*ocp_slope(cs)
        di(IP1) = di_deta
        di(IP2) = -di_deta
    end subroutine rate_derivs_exact

    ! =============================== corrected-mode time step ===============================
    subroutine time_terms(dt, Tt)
        !! Storage coefficients: the time-derivative part of each row is Tt*(c - c_old).
        real(dp), intent(in) :: dt
        real(dp), intent(out) :: Tt(NV,nj)
        Tt = 0.0_dp
        Tt(ICS,1) = -(eps_AM/dt)
        Tt(IC,2:s-1) = -(eps_sep/dt*dx(2:s-1))
        Tt(ICS,2:s-1) = -((1.0_dp - eps_sep)/dt)
        Tt(IC,s+1:nj-1) = -((eps/dt)*dx(s+1:nj-1))
        Tt(ICS,s:nj) = -(eps_AM/dt)
    end subroutine time_terms

    real(dp) function bounded_step(dcc)
        !! Largest step <= 1 keeping 0 < c and 0 < cs < cs_max (at most 90 % of the way to a bound).
        real(dp), intent(in) :: dcc(NV,nj)
        real(dp), parameter :: keep = 0.9_dp
        real(dp) :: csm
        integer :: j
        csm = cs_max()
        bounded_step = 1.0_dp
        do j = 1, nj
            if (dcc(IC,j) < 0) bounded_step = min(bounded_step, keep*(c(IC,j) - 0.0_dp)/(-dcc(IC,j)))
            if (dcc(ICS,j) < 0) bounded_step = min(bounded_step, keep*(c(ICS,j) - 0.0_dp)/(-dcc(ICS,j)))
            if (dcc(ICS,j) > 0) bounded_step = min(bounded_step, keep*(csm - c(ICS,j))/dcc(ICS,j))
        end do
        bounded_step = max(bounded_step, 0.0_dp)
    end function bounded_step

    subroutine newton_step(h, ok)
        !! One backward-Euler step of length h, solved with Newton's method; c is updated on success.
        real(dp), intent(in) :: h
        logical, intent(out) :: ok
        real(dp) :: c_old(NV,nj), Tt(NV,nj), scale(NV), lam
        integer :: k, st
        c_old = c
        call time_terms(h, Tt)
        scale = [c_bulk, 1.0_dp, 1.0_dp, cs_max()]
        ok = .false.
        do k = 1, newton_max_iter
            call assemble(h)
            G = G - Tt*(c - c_old)
            call band_solve(NV, nj, A, B, Dm, G, dc, st)
            if (st /= BAND_OK) exit
            lam = bounded_step(dc)
            c = c + lam*dc
            if (lam == 1.0_dp .and. maxval(abs(dc)/spread(scale, 2, nj)) <= newton_tol) then
                ok = .true.
                return
            end if
        end do
        c = c_old
    end subroutine newton_step

    subroutine advance(dt, t_done, stopped, ok)
        !! Advance by dt, halving the sub-step on Newton failure and stopping at a voltage cutoff.
        real(dp), intent(in) :: dt
        real(dp), intent(out) :: t_done
        logical, intent(out) :: stopped, ok
        real(dp), parameter :: min_dt = 1.0e-6_dp
        real(dp) :: hh, vv
        logical :: good
        t_done = 0.0_dp
        hh = dt
        stopped = .false.
        ok = .true.
        do while (t_done < dt)
            hh = min(hh, dt - t_done)
            call newton_step(hh, good)
            if (.not. good) then
                if (hh/2 < min_dt) then
                    ok = .false.
                    return
                end if
                hh = hh/2
                cycle
            end if
            t_done = t_done + hh
            vv = cell_voltage()
            if (vv <= V_min .or. vv >= V_max) then
                stopped = .true.
                return
            end if
        end do
    end subroutine advance

    ! =============================== assembly ===============================
    subroutine face_coeffs(e, dcat, ucat, dan, uan, cface, gphi2, dd, ff)
        !! Cation-flux (row 1) and ionic-current (row 3) coefficients of one face.
        real(dp), intent(in) :: e, dcat, ucat, dan, uan, cface, gphi2
        real(dp), intent(inout) :: dd(NV,NV), ff(NV,NV)
        real(dp) :: k
        k = z_plus**2*ucat + z_minus**2*uan
        dd(IC,IC)   = -(e*dcat)
        ff(IC,IC)   = -(e*z_plus*ucat*F*gphi2)
        dd(IC,IP2)  = -(e*z_plus*ucat*F*cface)
        dd(IP2,IC)  = -(e*F*(z_plus*dcat + z_minus*dan))
        ff(IP2,IC)  = -(e*F**2*k*gphi2)
        dd(IP2,IP2) = -(e*F**2*k*cface)
    end subroutine face_coeffs

    subroutine assemble(dt)
        real(dp), intent(in) :: dt
        real(dp) :: dW(NV,NV), dE(NV,NV), fW(NV,NV), fE(NV,NV), rj(NV,NV), gg(NV)
        real(dp) :: cW(NV), cE(NV), gW(NV), gE(NV), i, di(NV)
        integer :: j, k

        A = 0.0_dp; B = 0.0_dp; Dm = 0.0_dp; G = 0.0_dp
        do j = 1, nj
            dW = 0.0_dp; dE = 0.0_dp; fW = 0.0_dp; fE = 0.0_dp; rj = 0.0_dp; gg = 0.0_dp
            cW = 0.0_dp; cE = 0.0_dp; gW = 0.0_dp; gE = 0.0_dp
            if (j > 1) then
                cW = aW(j)*c(:,j) + (1.0_dp - aW(j))*c(:,j-1)
                gW = bW(j)*(c(:,j) - c(:,j-1))
            end if
            if (j < nj) then
                cE = aE(j)*c(:,j+1) + (1.0_dp - aE(j))*c(:,j)
                gE = bE(j)*(c(:,j+1) - c(:,j))
            end if
            call rate_derivs(c(IC,j), c(ICS,j), c(IP1,j), c(IP2,j), i, di)

            if (j == 1) then
                ! ---- Li-foil face ----
                dE(IC,IC)  = -(eps_sep_face*dcat_s)
                fE(IC,IC)  = -(eps_sep_face*z_plus*ucat_s*F*gE(IP2))
                dE(IC,IP2) = -(eps_sep_face*z_plus*ucat_s*F*cE(IC))
                gg(IC) = -i_app/F + (dE(IC,IC)*gE(IC) + fE(IC,IC)*cE(IC))
                rj(ICS,ICS) = 0.0_dp - 1.0_dp*eps_AM/dt
                dE(IP1,IP1) = -(1.0_dp - eps)*sigma
                gg(IP1) = phi1_sign*dE(IP1,IP1)*gE(IP1)
                rj(IP2,IP2) = 1.0_dp
                gg(IP2) = 0.0_dp - c(IP2,j)

            else if (j < s) then
                ! ---- separator interior ----
                call face_coeffs(eps_sep, dcat_s, ucat_s, dan_s, uan_s, cW(IC), gW(IP2), dW, fW)
                call face_coeffs(eps_sep, dcat_s, ucat_s, dan_s, uan_s, cE(IC), gE(IP2), dE, fE)
                rj(IC,IC) = -(eps_sep/dt*dx(j))
                gg(IC) = 0.0_dp - (fW(IC,IC)*cW(IC) + dW(IC,IC)*gW(IC)) + (fE(IC,IC)*cE(IC) + dE(IC,IC)*gE(IC))
                rj(ICS,ICS) = -(1.0_dp*(1.0_dp - eps_sep)/dt)
                dW(IP1,IP1) = -(1.0_dp - eps_sep)*sigma
                dE(IP1,IP1) = -(1.0_dp - eps_sep)*sigma
                gg(IP1) = 0.0_dp - (fW(IP1,IP1)*cW(IP1) + dW(IP1,IP1)*gW(IP1)) &
                                 + (fE(IP1,IP1)*cE(IP1) + dE(IP1,IP1)*gE(IP1))
                gg(IP2) = 0.0_dp - current(fW, dW, cW, gW) + current(fE, dE, cE, gE)

            else if (j == s) then
                ! ---- separator/cathode interface ----
                call face_coeffs(eps_sep_face, dcat_s, ucat_s, dan_s, uan_s, cW(IC), gW(IP2), dW, fW)
                call face_coeffs(eps, dcat_c, ucat_c, dan_c, uan_c, cE(IC), gE(IP2), dE, fE)
                gg(IC) = 0.0_dp - (fW(IC,IC)*cW(IC) + dW(IC,IC)*gW(IC)) + (fE(IC,IC)*cE(IC) + dE(IC,IC)*gE(IC))
                dW(IP1,IP1) = -(1.0_dp - eps)*sigma
                dE(IP1,IP1) = -(1.0_dp - eps)*sigma
                gg(IP1) = 0.0_dp - (fW(IP1,IP1)*cW(IP1) + dW(IP1,IP1)*gW(IP1)) &
                                 + (fE(IP1,IP1)*cE(IP1) + dE(IP1,IP1)*gE(IP1))
                gg(IP2) = 0.0_dp - current(fW, dW, cW, gW) + current(fE, dE, cE, gE)
                call solid_row(i, di, dt, rj, gg)

            else if (j < nj) then
                ! ---- cathode interior ----
                call face_coeffs(eps, dcat_c, ucat_c, dan_c, uan_c, cW(IC), gW(IP2), dW, fW)
                call face_coeffs(eps, dcat_c, ucat_c, dan_c, uan_c, cE(IC), gE(IP2), dE, fE)
                do k = 1, NV
                    rj(IC,k) = (spec_a*di(k)/F)*dx(j)
                    rj(IP1,k) = -((spec_a*di(k))*dx(j))
                    rj(IP2,k) = (spec_a*di(k))*dx(j)
                end do
                rj(IC,IC) = (spec_a*di(IC)/F)*dx(j) - (eps/dt)*dx(j)
                gg(IC) = -((spec_a*i/F)*dx(j)) - (fW(IC,IC)*cW(IC) + dW(IC,IC)*gW(IC)) &
                                               + (fE(IC,IC)*cE(IC) + dE(IC,IC)*gE(IC))
                dW(IP1,IP1) = -(1.0_dp - eps)*sigma
                dE(IP1,IP1) = -(1.0_dp - eps)*sigma
                gg(IP1) = (spec_a*i)*dx(j) - (fW(IP1,IP1)*cW(IP1) + dW(IP1,IP1)*gW(IP1)) &
                                           + (fE(IP1,IP1)*cE(IP1) + dE(IP1,IP1)*gE(IP1))
                gg(IP2) = -((spec_a*i)*dx(j)) - current(fW, dW, cW, gW) + current(fE, dE, cE, gE)
                call solid_row(i, di, dt, rj, gg)

            else
                ! ---- current collector ----
                call face_coeffs(eps, dcat_c, ucat_c, dan_c, uan_c, cW(IC), gW(IP2), dW, fW)
                gg(IC) = 0.0_dp - dW(IC,IC)*gW(IC) - fW(IC,IC)*cW(IC)
                dW(IP1,IP1) = -(1.0_dp - eps)*sigma
                gg(IP1) = i_app - dW(IP1,IP1)*gW(IP1)
                gg(IP2) = 0.0_dp - dW(IP2,IC)*gW(IC) - fW(IP2,IC)*cW(IC)
                call solid_row(i, di, dt, rj, gg)
            end if

            ! control-volume blocks: A dc(j-1) + B dc(j) + D dc(j+1) = G
            if (j == 1) then
                B(:,:,j) = rj - (1.0_dp - aE(j))*fE + bE(j)*dE
                Dm(:,:,j) = -(aE(j)*fE) - bE(j)*dE
            else if (j == nj) then
                A(:,:,j) = (1.0_dp - aW(j))*fW - bW(j)*dW
                B(:,:,j) = rj + bW(j)*dW + aW(j)*fW
            else
                A(:,:,j) = (1.0_dp - aW(j))*fW - bW(j)*dW
                B(:,:,j) = rj + bW(j)*dW + aW(j)*fW - (1.0_dp - aE(j))*fE + bE(j)*dE
                Dm(:,:,j) = -(aE(j)*fE) - bE(j)*dE
            end if
            G(:,j) = gg
        end do
    end subroutine assemble

    real(dp) function current(ff, dd, cf, gf)
        !! Ionic-current face value used in the residual (D-11: faithful mode omits diffusion).
        real(dp), intent(in) :: ff(NV,NV), dd(NV,NV), cf(NV), gf(NV)
        current = ff(IP2,IP2)*cf(IP2) + dd(IP2,IP2)*gf(IP2)
        if (full_current) current = current + dd(IP2,IC)*gf(IC)
    end function current

    subroutine solid_row(i, di, dt, rj, gg)
        !! eps_AM dcs/dt = -a i_n / F  (uniform particles)
        real(dp), intent(in) :: i, di(NV), dt
        real(dp), intent(inout) :: rj(NV,NV), gg(NV)
        integer :: k
        do k = 1, NV
            rj(ICS,k) = -(spec_a*di(k)/F)
        end do
        rj(ICS,ICS) = -(spec_a*di(ICS)/F) - 1.0_dp*eps_AM/dt
        gg(ICS) = +(spec_a*i/F)
    end subroutine solid_row

    ! =============================== output ===============================
    real(dp) function li_eta()
        !! Overpotential of the lithium counter electrode (output only).
        real(dp) :: i0_li, alpha
        i0_li = F*k_Li*(c(IC,1)**0.5_dp)*(c_Li_ref**0.5_dp)
        alpha = 0.5_dp
        if (state == 'C') then
            li_eta = 0.5_dp*log(i_app/i0_li)/(alpha*F/(R*T))
        else if (state == 'D') then
            li_eta = -(0.5_dp*log(i_app/i0_li))/(alpha*F/(R*T))
        else
            li_eta = 0.0_dp
        end if
    end function li_eta

    real(dp) function cell_voltage()
        cell_voltage = c(IP1,nj) + li_eta()
    end function cell_voltage

    subroutine write_row(header)
        logical, intent(in) :: header
        real(dp) :: equiv, i0_li, eta
        if (header) then
            write(ounit,'(A5,1X,2(A12,1X),20(A15,1X))') 'State', 'Time', 'Voltage', 'Equivalence', 'Anode_Eta', &
                'anode_exchange_c', 'Edge_c0'
            write(ounit,'(A5,1X,2(A12,1X),20(A15,1X))') 'CDR', 'hours', 'Volts', 'electron_equivs', 'mV', &
                'mA/cm2', 'mol/cm3'
        end if
        equiv = mAhg*M*lit36/F
        i0_li = F*k_Li*(c(IC,1)**0.5_dp)*(c_Li_ref**0.5_dp)
        eta = li_eta()
        write(ounit,'(A5,1X,2(F12.5,1X),20(ES15.5,1X))') state, time/real(3600, dp), c(IP1,nj) + eta, equiv, &
            eta*1.0e3_dp, i0_li*1.0e3_dp, c(IC,1)
    end subroutine write_row

end program lfp
