subroutine sparsity(modes, nth_tier, index_minus, index_plus, d_operators, hamil, rho_0, rho_sparsity, &
                    N_nonzeros, dimen_el, dimen_hierarchy, number_index_plus, truncation_tier, number_mode, number_el_state)
        implicit none
        integer, intent(in) :: dimen_el, dimen_hierarchy, number_index_plus
        integer, intent(in) :: truncation_tier, number_mode, number_el_state
        integer, intent(in), dimension(0:number_mode-1, 4) :: modes
        integer, intent(in), dimension(0:dimen_hierarchy-1) :: nth_tier
        integer, intent(in), dimension(0:dimen_hierarchy-1, truncation_tier, 4) :: index_minus
        integer, intent(in), dimension(0:number_index_plus-1, 0:number_mode-1, 3) :: index_plus
        logical, intent(in), dimension(0:1,number_el_state, dimen_el, dimen_el) :: d_operators
        logical, intent(in), dimension(dimen_el, dimen_el) :: hamil
        logical, intent(in), dimension(dimen_el, dimen_el) :: rho_0
        integer, intent(out), dimension(dimen_el, dimen_el, 0:dimen_hierarchy-1) :: rho_sparsity
        integer, intent(out)  :: N_nonzeros 
        logical, dimension(dimen_el, dimen_el, 0:dimen_hierarchy-1) :: rho_temp
        integer :: i_hierarchy,ni,nj

        rho_temp = .False.
        rho_temp(:,:,0) = rho_0(:,:)
        call one_step_propagation(modes, nth_tier, index_minus, index_plus, d_operators, hamil, &
                                    rho_temp, dimen_el, dimen_hierarchy, number_index_plus, truncation_tier, &
                                    number_mode, number_el_state)
        N_nonzeros = 0
        rho_sparsity = -1
        !open(100, file='ADM_logic.txt')
        do i_hierarchy = 0, dimen_hierarchy - 1
        !write(100,*) i_hierarchy, ((rho_temp(ni,nj,i_hierarchy), nj=1,dimen_el), ni=1,dimen_el)
           do ni = 1, dimen_el
              do nj = 1, dimen_el
                 if (rho_temp(ni,nj,i_hierarchy)) then
                    N_nonzeros = N_nonzeros +1
                    rho_sparsity(ni,nj,i_hierarchy) = N_nonzeros
                 endif
              enddo
           enddo
           !write(200,*) i_hierarchy, ((rho_sparsity(ni,nj,i_hierarchy), nj=1,dimen_el), ni=1,dimen_el)
        enddo
        !close(100) 


end subroutine       

subroutine one_step_propagation(modes, nth_tier, index_minus, index_plus, d_operators, hamil, rho_in, &
                                dimen_el, dimen_hierarchy, number_index_plus, truncation_tier, number_mode, &
                                number_el_state)
        implicit none
        integer, intent(in) :: dimen_el, dimen_hierarchy, number_index_plus
        integer, intent(in) :: truncation_tier, number_mode, number_el_state
        integer, intent(in), dimension(0:number_mode-1, 4) :: modes
        integer, intent(in), dimension(0:dimen_hierarchy-1) :: nth_tier
        integer, intent(in), dimension(0:dimen_hierarchy-1, truncation_tier, 4) :: index_minus
        integer, intent(in), dimension(0:number_index_plus-1, 0:number_mode-1, 3) :: index_plus
        logical, intent(in), dimension(0:1,number_el_state, dimen_el, dimen_el) :: d_operators
        logical, intent(in), dimension(dimen_el, dimen_el) :: hamil
        logical, intent(inout), dimension(dimen_el, dimen_el, 0:dimen_hierarchy-1) :: rho_in
        logical, dimension(dimen_el, dimen_el, 0:dimen_hierarchy-1) :: rho_temp
        logical, dimension(dimen_el, dimen_el, 0:dimen_hierarchy-1) :: drho
        integer :: nth_order, i_hierarchy, i_hierarchy_minus, i_hierarchy_plus, tier
        integer :: i_mode, i_tier, i_elec, i_sign, i_sym_on

        rho_temp = rho_in
        do nth_order = 1, 4
           drho=.False.
           do i_hierarchy = 0, dimen_hierarchy - 1
              drho(:, :, i_hierarchy) = (drho(:, :, i_hierarchy) .or. ( MATMUL(hamil, rho_temp(:, :, i_hierarchy) ) &
                                                                   .or. MATMUL(rho_temp(:, :, i_hierarchy), hamil ) ) )
              tier = nth_tier(i_hierarchy)
              if ( tier > 0 ) then
                 drho(:, :, i_hierarchy) = (drho(:, :, i_hierarchy) .or. rho_temp(:, :, i_hierarchy))

                 do i_tier = 1, tier
                    i_mode = index_minus(i_hierarchy, i_tier, 1)
                    i_elec=modes(i_mode,2)
                    i_sign=modes(i_mode,4)
                    i_hierarchy_minus = index_minus(i_hierarchy, i_tier, 2)
                    i_sym_on = index_minus(i_hierarchy, i_tier, 4)
                    if ( i_sym_on == 0 ) then
                    drho(:, :, i_hierarchy) = (drho(:, :, i_hierarchy) .or. (&
                             MATMUL( d_operators(i_sign, i_elec+1, :, :), rho_temp(:, :, i_hierarchy_minus) ) .or. &
                             MATMUL( rho_temp(:, :, i_hierarchy_minus), d_operators(i_sign, i_elec+1, :, :) ) ) )
                    elseif ( i_sym_on == 1 ) then
                    drho(:, :, i_hierarchy) = (drho(:, :, i_hierarchy) .or. (&
                             MATMUL( d_operators(i_sign, i_elec+1, :, :), transpose(rho_temp(:, :, i_hierarchy_minus)) ) .or. &
                             MATMUL( transpose(rho_temp(:, :, i_hierarchy_minus)), d_operators(i_sign, i_elec+1, :, :) ) ) )
                    endif

                 enddo !i_tier
              endif ! tier>0

              if ( tier < truncation_tier ) then
                 do i_mode = 0, number_mode - 1
                    i_elec = modes(i_mode, 2)
                    i_sign = modes(i_mode, 4)
                    i_hierarchy_plus = index_plus(i_hierarchy, i_mode, 1)
                    i_sym_on = index_plus(i_hierarchy, i_mode, 3)
                    if ( i_hierarchy_plus .ne. -1 ) then
                       i_sign = 1 - i_sign
                       if ( i_sym_on == 0 ) then
                       drho(:, :, i_hierarchy) = (drho(:, :, i_hierarchy) .or. (&
                               MATMUL(d_operators(i_sign, i_elec+1, :, :), rho_temp(:, :, i_hierarchy_plus) ) .or. &
                               MATMUL(rho_temp(:, :, i_hierarchy_plus), d_operators(i_sign, i_elec+1, :, :) ) ) )
                       elseif ( i_sym_on == 1 ) then
                       drho(:, :, i_hierarchy) = (drho(:, :, i_hierarchy) .or. (&
                               MATMUL(d_operators(i_sign, i_elec+1, :, :), transpose(rho_temp(:, :, i_hierarchy_plus)) ) .or. &
                               MATMUL(transpose(rho_temp(:, :, i_hierarchy_plus)), d_operators(i_sign, i_elec+1, :, :) ) ) )
                       endif
                    endif
                 enddo
              endif

           enddo !i_hierarchy

        rho_in = (rho_in .or. drho)
        rho_temp = drho
        enddo
end subroutine

subroutine sparse_matrix_elements_a(modes, nth_tier, index_minus, index_plus, d_operators, hamil, &
                                rho_sparsity, rho_nonzeros, number_pair, number_pair_hamil_l, number_pair_hamil_r, &
                                number_nonzeros, dimen_el, dimen_hierarchy, &
                                number_index_plus, truncation_tier,  number_mode, number_el_state)
        implicit none
        integer, intent(in) :: number_nonzeros, dimen_el, dimen_hierarchy, number_index_plus
        integer, intent(in) :: truncation_tier, number_mode, number_el_state
        integer, intent(in), dimension(0:number_mode-1, 4) :: modes
        integer, intent(in), dimension(0:dimen_hierarchy-1) :: nth_tier
        integer, intent(in), dimension(0:dimen_hierarchy-1, truncation_tier, 4) :: index_minus
        integer, intent(in), dimension(0:number_index_plus-1, 0:number_mode-1, 3) :: index_plus
        integer, intent(in), dimension(0:1,number_el_state, dimen_el, dimen_el) :: d_operators
        complex*16, intent(in), dimension(dimen_el, dimen_el) :: hamil
        integer, intent(in), dimension(dimen_el, dimen_el, 0:dimen_hierarchy-1) :: rho_sparsity
        integer, intent(out), dimension(number_nonzeros, 3) :: rho_nonzeros
        integer, intent(out) :: number_pair
        integer, intent(out) :: number_pair_hamil_l, number_pair_hamil_r
        integer :: i_hierarchy, i_hierarchy_minus, i_hierarchy_plus, tier
        integer :: i_mode, i_tier, i_elec, i_sign, i_sym_on
        integer :: op_value
        integer :: ni, nj, nk, i_count 
        integer :: i_a,i_b

        i_count = 0
        do i_hierarchy = 0, dimen_hierarchy - 1
           do ni = 1, dimen_el
              do nj = 1, dimen_el
                 if (rho_sparsity(ni,nj,i_hierarchy) /= -1) then
                    i_count = i_count + 1
                    rho_nonzeros(i_count,1) = i_hierarchy
                    rho_nonzeros(i_count,2) = ni
                    rho_nonzeros(i_count,3) = nj
                    !write(500,*) i_count, (rho_nonzeros(i_count, nk), nk = 1,3)
                 endif
              enddo
           enddo
        enddo

        number_pair = 0
        number_pair_hamil_l = 0
        number_pair_hamil_r = 0
        do i_a = 1, number_nonzeros
           i_hierarchy = rho_nonzeros(i_a,1)
           tier = nth_tier(i_hierarchy)
           ni = rho_nonzeros(i_a,2)
           nj = rho_nonzeros(i_a,3)

           if ( tier > 0 ) then
              do i_tier = 1, tier
                 i_mode = index_minus(i_hierarchy, i_tier, 1)
                 i_elec=modes(i_mode,2) + 1
                 i_sign=modes(i_mode,4)
                 i_hierarchy_minus = index_minus(i_hierarchy, i_tier, 2)
                 i_sym_on = index_minus(i_hierarchy, i_tier, 4)

                 if ( i_sym_on == 0 ) then
                    do nk = 1, dimen_el
                       i_b = rho_sparsity(nk,nj,i_hierarchy_minus)
                       op_value = d_operators(i_sign, i_elec, ni, nk)
                       if ( (i_b .ne. -1) .and. (op_value .ne. 0) ) number_pair = number_pair + 1

                       i_b = rho_sparsity(ni,nk,i_hierarchy_minus)
                       op_value = d_operators(i_sign, i_elec, nk, nj)
                       if ( (i_b .ne. -1) .and. (op_value .ne. 0) ) number_pair = number_pair + 1
                    enddo
                 elseif ( i_sym_on == 1 ) then
                    do nk = 1, dimen_el
                       i_b = rho_sparsity(nj,nk,i_hierarchy_minus)
                       op_value = d_operators(i_sign, i_elec, ni, nk)
                       if ( (i_b .ne. -1) .and. (op_value .ne. 0) ) number_pair = number_pair + 1

                       i_b = rho_sparsity(nk,ni,i_hierarchy_minus)
                       op_value = d_operators(i_sign, i_elec, nk, nj) 
                       if ( (i_b .ne. -1) .and. (op_value .ne. 0) ) number_pair = number_pair + 1
                    enddo
                 endif
              enddo !i_tier

           endif ! tier>0

           do nk = 1, dimen_el
              i_b = rho_sparsity(nk,nj,i_hierarchy)
              if ( (i_b .ne. -1) .and. (hamil(ni,nk) .ne. 0.) ) number_pair_hamil_l = number_pair_hamil_l + 1

              i_b = rho_sparsity(ni,nk,i_hierarchy)
              if ( (i_b .ne. -1) .and. (hamil(nk,nj) .ne. 0.) ) number_pair_hamil_r = number_pair_hamil_r + 1
           enddo

          if ( tier < truncation_tier ) then
             do i_mode = 0, number_mode - 1
                i_elec = modes(i_mode, 2) + 1
                i_sign = 1 - modes(i_mode, 4)
                i_hierarchy_plus = index_plus(i_hierarchy, i_mode, 1)
                i_sym_on = index_plus(i_hierarchy, i_mode, 3)

                if ( i_hierarchy_plus .ne. -1 ) then

                   if ( i_sym_on == 0 ) then
                      do nk = 1, dimen_el
                         i_b = rho_sparsity(nk,nj,i_hierarchy_plus)
                         op_value = d_operators(i_sign, i_elec, ni, nk)
                         if ( (i_b .ne. -1) .and. (op_value .ne. 0) ) number_pair = number_pair + 1

                         i_b = rho_sparsity(ni,nk,i_hierarchy_plus)
                         op_value = d_operators(i_sign, i_elec, nk, nj) 
                         if ( (i_b .ne. -1) .and. (op_value .ne. 0) ) number_pair = number_pair + 1
                      enddo
                   elseif ( i_sym_on == 1 ) then
                      do nk = 1, dimen_el
                         i_b = rho_sparsity(nj,nk,i_hierarchy_plus)
                         op_value = d_operators(i_sign, i_elec, ni, nk)
                         if ( (i_b .ne. -1) .and. (op_value .ne. 0) ) number_pair = number_pair + 1

                         i_b = rho_sparsity(nk,ni,i_hierarchy_plus)
                         op_value = d_operators(i_sign, i_elec, nk, nj) 
                         if ( (i_b .ne. -1) .and. (op_value .ne. 0) ) number_pair = number_pair + 1
                      enddo
                   endif
                endif

             enddo
          endif

        enddo !i_a

end subroutine

subroutine sparse_matrix_elements_b(modes, nth_tier, indices, index_minus, index_plus, gammas, etas, d_operators, hamil, &
                                rho_sparsity, rho_nonzeros, matrix, matrix_hamil_l, matrix_hamil_r, matrix_elem, &
                                matrix_gammas, number_pair, number_pair_hamil_l, number_pair_hamil_r, number_nonzeros, &
                                dimen_el, dimen_hierarchy, number_index_plus, truncation_tier,  number_mode, number_lead, &
                                number_el_state, number_Pade_poles, sign_pm)
        implicit none
        integer, intent(in) :: number_pair, number_pair_hamil_l, number_pair_hamil_r
        integer, intent(in) :: number_nonzeros, dimen_el, dimen_hierarchy, number_index_plus
        integer, intent(in) :: truncation_tier, number_mode
        integer, intent(in) :: number_lead, number_el_state, number_Pade_poles, sign_pm
        integer, intent(in), dimension(0:number_mode-1, 4) :: modes
        integer, intent(in), dimension(0:dimen_hierarchy-1) :: nth_tier
        integer, intent(in), dimension(0:dimen_hierarchy-1, truncation_tier) :: indices
        integer, intent(in), dimension(0:dimen_hierarchy-1, truncation_tier, 4) :: index_minus
        integer, intent(in), dimension(0:number_index_plus-1, 0:number_mode-1, 3) :: index_plus
        complex*16, intent(in), dimension(0:number_lead-1, 0:number_Pade_poles, 0:sign_pm-1) :: gammas
        complex*16, intent(in), dimension(0:number_lead-1, 0:number_Pade_poles) :: etas
        integer, intent(in), dimension(0:1,number_el_state, dimen_el, dimen_el) :: d_operators
        complex*16, intent(in), dimension(dimen_el, dimen_el) :: hamil
        integer, intent(in), dimension(dimen_el, dimen_el, 0:dimen_hierarchy-1) :: rho_sparsity
        integer, intent(in), dimension(number_nonzeros, 3) :: rho_nonzeros
        integer, intent(out), dimension(number_pair, 6) :: matrix
        integer, intent(out), dimension(number_pair_hamil_l, 4) :: matrix_hamil_l
        integer, intent(out), dimension(number_pair_hamil_r, 4) :: matrix_hamil_r
        complex*16, intent(out), dimension(number_pair) :: matrix_elem
        complex*16, intent(out), dimension(number_nonzeros) :: matrix_gammas
        complex*16, parameter :: ci = (0.d0, 1.d0)
        integer :: i_hierarchy, i_hierarchy_minus, i_hierarchy_plus, tier
        integer :: i_mode, i_tier, i_lead, i_elec, i_pole, i_sign, i_sym_on, i_permutation
        complex*16 :: gamma_sum
        integer :: op_value
        real*8 :: sym_button, permutation_button
        integer :: ni, nj, nk, i_count,  i_count_hamil_l, i_count_hamil_r
        integer :: i_a,i_b

        matrix_gammas = 0.d0
        matrix(:,1) = -1
        matrix(:,2) = -1
        matrix(:,3) = 0
        matrix(:,4) = 0
        matrix_hamil_l = 0
        matrix_hamil_r = 0
        matrix_elem = 0.d0
        i_count_hamil_l = 0
        i_count_hamil_r = 0
        i_count = 0
        do i_a = 1, number_nonzeros
           i_hierarchy = rho_nonzeros(i_a,1)
           tier = nth_tier(i_hierarchy)
           ni = rho_nonzeros(i_a,2)
           nj = rho_nonzeros(i_a,3)

           if ( tier > 0 ) then
              do i_tier = 1, tier
                 i_mode = index_minus(i_hierarchy, i_tier, 1)
                 i_lead=modes(i_mode,1)
                 i_elec=modes(i_mode,2)
                 i_pole=modes(i_mode,3)
                 i_sign=modes(i_mode,4)
                 i_hierarchy_minus = index_minus(i_hierarchy, i_tier, 2)
                 i_permutation=index_minus(i_hierarchy,i_tier,3)
                 i_sym_on = index_minus(i_hierarchy, i_tier, 4)

                 permutation_button=(-1.d0)**i_permutation
                 if ( i_sym_on == 0 ) then
                    do nk = 1, dimen_el
                       i_b = rho_sparsity(nk,nj,i_hierarchy_minus)
                       op_value = d_operators(i_sign, i_elec+1, ni, nk)
                       if ( (i_b .ne. -1) .and. (op_value .ne. 0) ) then
                          i_count = i_count + 1
                          matrix(i_count,1) =  i_a
                          matrix(i_count,2) =  i_b
                          matrix(i_count,5) =  i_lead + 1
                          matrix(i_count,6) =  i_elec + 1
                          matrix_elem(i_count) =  - ci * (-1)**(tier-i_tier) * permutation_button * &
                                                       etas(i_lead, i_pole) * op_value  
                       endif

                       i_b = rho_sparsity(ni,nk,i_hierarchy_minus)
                       op_value = d_operators(i_sign, i_elec+1, nk, nj)
                       if ( (i_b .ne. -1) .and. (op_value .ne. 0) ) then
                          i_count = i_count + 1
                          matrix(i_count,1) =  i_a
                          matrix(i_count,2) =  i_b
                          matrix(i_count,4) =  1
                          matrix(i_count,5) =  i_lead + 1
                          matrix(i_count,6) =  i_elec + 1
                          matrix_elem(i_count) =  ci * (-1)**(2*tier-i_tier-1)  * permutation_button * &
                                                       conjg( etas(i_lead, i_pole) ) * op_value
                       endif
                    enddo
                 elseif ( i_sym_on == 1 ) then
                    sym_button = (-1.d0)**floor( (tier-1.d0)/2.d0 )
                    do nk = 1, dimen_el
                       i_b = rho_sparsity(nj,nk,i_hierarchy_minus)
                       op_value = d_operators(i_sign, i_elec+1, ni, nk)
                       if ( (i_b .ne. -1) .and. (op_value .ne. 0) ) then
                          i_count = i_count + 1
                          matrix(i_count,1) =  i_a
                          matrix(i_count,2) =  i_b
                          matrix(i_count,3) =  i_sym_on
                          matrix(i_count,5) =  i_lead + 1
                          matrix(i_count,6) =  i_elec + 1
                          matrix_elem(i_count) =  - ci *(-1)**(tier-i_tier) * permutation_button * &
                                                       etas(i_lead, i_pole) * sym_button * op_value 
                       endif

                       i_b = rho_sparsity(nk,ni,i_hierarchy_minus)
                       op_value = d_operators(i_sign, i_elec+1, nk, nj) 
                       if ( (i_b .ne. -1) .and. (op_value .ne. 0) ) then
                          i_count = i_count + 1
                          matrix(i_count,1) =  i_a
                          matrix(i_count,2) =  i_b
                          matrix(i_count,3) =  i_sym_on
                          matrix(i_count,4) =  1 
                          matrix(i_count,5) =  i_lead + 1
                          matrix(i_count,6) =  i_elec + 1
                          matrix_elem(i_count) =  ci * (-1)**(2*tier-i_tier-1)  * permutation_button * &
                                                       sym_button * conjg( etas(i_lead, i_pole) ) * op_value 
                       endif
                    enddo
                 endif
              enddo !i_tier

              gamma_sum = 0.d0
              do i_tier = 1, tier
                 i_mode = indices(i_hierarchy, i_tier)
                 i_lead = modes(i_mode, 1)
                 i_elec = modes(i_mode, 2)
                 i_pole = modes(i_mode, 3)
                 i_sign = modes(i_mode, 4)
                 gamma_sum = gamma_sum + gammas(i_lead, i_pole, i_sign) 
              enddo 
              matrix_gammas(i_a) = - gamma_sum
           endif ! tier>0

           do nk = 1, dimen_el
              i_b = rho_sparsity(nk,nj,i_hierarchy)
              if ( (i_b .ne. -1) .and. (hamil(ni,nk) .ne. 0.) ) then
                 i_count_hamil_l = i_count_hamil_l + 1
                 matrix_hamil_l(i_count_hamil_l,1) =  i_a
                 matrix_hamil_l(i_count_hamil_l,2) =  i_b
                 matrix_hamil_l(i_count_hamil_l,3) =  ni
                 matrix_hamil_l(i_count_hamil_l,4) =  nk
              endif

              i_b = rho_sparsity(ni,nk,i_hierarchy)
              if ( (i_b .ne. -1) .and. (hamil(nk,nj) .ne. 0.) )then
                 i_count_hamil_r = i_count_hamil_r + 1
                 matrix_hamil_r(i_count_hamil_r,1) =  i_a
                 matrix_hamil_r(i_count_hamil_r,2) =  i_b
                 matrix_hamil_r(i_count_hamil_r,3) =  nk
                 matrix_hamil_r(i_count_hamil_r,4) =  nj
              endif
           enddo

          if ( tier < truncation_tier ) then
             do i_mode = 0, number_mode - 1
                i_lead = modes(i_mode,1)
                i_elec = modes(i_mode, 2)
                i_sign = 1 - modes(i_mode, 4)
                i_hierarchy_plus = index_plus(i_hierarchy, i_mode, 1)
                i_permutation = index_plus(i_hierarchy, i_mode, 2)
                i_sym_on = index_plus(i_hierarchy, i_mode, 3)

                permutation_button = (-1.d0)**i_permutation

                if ( i_hierarchy_plus .ne. -1 ) then

                   if ( i_sym_on == 0 ) then
                      do nk = 1, dimen_el
                         i_b = rho_sparsity(nk,nj,i_hierarchy_plus)
                         op_value = d_operators(i_sign, i_elec+1, ni, nk)
                         if ( (i_b .ne. -1) .and. (op_value .ne. 0) ) then
                            i_count = i_count + 1
                            matrix(i_count,1) =  i_a
                            matrix(i_count,2) =  i_b
                            matrix(i_count,5) =  i_lead + 1
                            matrix(i_count,6) =  i_elec + 1
                            matrix_elem(i_count) =  - ci * permutation_button * op_value 
                         endif

                         i_b = rho_sparsity(ni,nk,i_hierarchy_plus)
                         op_value = d_operators(i_sign, i_elec+1, nk, nj) 
                         if ( (i_b .ne. -1) .and. (op_value .ne. 0) ) then
                            i_count = i_count + 1
                            matrix(i_count,1) =  i_a
                            matrix(i_count,2) =  i_b
                            matrix(i_count,4) =  1
                            matrix(i_count,5) =  i_lead + 1
                            matrix(i_count,6) =  i_elec + 1
                            matrix_elem(i_count) =  - ci * permutation_button * (-1)**(tier+1)* op_value 
                         endif
                      enddo
                   elseif ( i_sym_on == 1 ) then
                      sym_button = (-1.d0)**floor( (tier+1.d0)/2.d0 )
                      do nk = 1, dimen_el
                         i_b = rho_sparsity(nj,nk,i_hierarchy_plus)
                         op_value = d_operators(i_sign, i_elec+1, ni, nk)
                         if ( (i_b .ne. -1) .and. (op_value .ne. 0) ) then
                            i_count = i_count + 1
                            matrix(i_count,1) =  i_a
                            matrix(i_count,2) =  i_b
                            matrix(i_count,3) =  i_sym_on
                            matrix(i_count,5) =  i_lead + 1
                            matrix(i_count,6) =  i_elec + 1
                            matrix_elem(i_count) =  - ci * permutation_button * sym_button * op_value 
                         endif

                         i_b = rho_sparsity(nk,ni,i_hierarchy_plus)
                         op_value = d_operators(i_sign, i_elec+1, nk, nj) 
                         if ( (i_b .ne. -1) .and. (op_value .ne. 0) ) then
                            i_count = i_count + 1
                            matrix(i_count,1) =  i_a
                            matrix(i_count,2) =  i_b
                            matrix(i_count,3) =  i_sym_on
                            matrix(i_count,4) =  1
                            matrix(i_count,5) =  i_lead + 1
                            matrix(i_count,6) =  i_elec + 1
                            matrix_elem(i_count) =  - ci * permutation_button * (-1)**(tier+1)*  sym_button *  op_value 
                         endif
                      enddo
                   endif
                endif

             enddo
          endif

        enddo !i_a

end subroutine


subroutine sparse_matrix_elements_a_wbl(modes, nth_tier, index_minus, index_plus, d_operators, hamil, &
                                rho_sparsity, rho_nonzeros, number_pair, number_pair_hamil_l, number_pair_hamil_r, &
                                number_pair_wbl, number_nonzeros, dimen_el, dimen_hierarchy, &
                                number_index_plus, truncation_tier,  number_mode, number_el_state)
        implicit none
        integer, intent(in) :: number_nonzeros, dimen_el, dimen_hierarchy, number_index_plus
        integer, intent(in) :: truncation_tier, number_mode, number_el_state
        integer, intent(in), dimension(0:number_mode-1, 4) :: modes
        integer, intent(in), dimension(0:dimen_hierarchy-1) :: nth_tier
        integer, intent(in), dimension(0:dimen_hierarchy-1, truncation_tier, 4) :: index_minus
        integer, intent(in), dimension(0:number_index_plus-1, 0:number_mode-1, 3) :: index_plus
        integer, intent(in), dimension(0:1,number_el_state, dimen_el, dimen_el) :: d_operators
        complex*16, intent(in), dimension(dimen_el, dimen_el) :: hamil
        integer, intent(in), dimension(dimen_el, dimen_el, 0:dimen_hierarchy-1) :: rho_sparsity
        integer, intent(out), dimension(number_nonzeros, 3) :: rho_nonzeros
        integer, intent(out) :: number_pair
        integer, intent(out) :: number_pair_hamil_l, number_pair_hamil_r
        integer, intent(out) :: number_pair_wbl
        integer :: i_hierarchy, i_hierarchy_minus, i_hierarchy_plus, tier
        integer :: i_mode, i_tier, i_elec, i_sign, i_sym_on
        integer :: op_value
        integer :: ni, nj, nk, nl, i_count 
        integer :: i_a,i_b

        i_count = 0
        do i_hierarchy = 0, dimen_hierarchy - 1
           do ni = 1, dimen_el
              do nj = 1, dimen_el
                 if (rho_sparsity(ni,nj,i_hierarchy) /= -1) then
                    i_count = i_count + 1
                    rho_nonzeros(i_count,1) = i_hierarchy
                    rho_nonzeros(i_count,2) = ni
                    rho_nonzeros(i_count,3) = nj
                    !write(500,*) i_count, (rho_nonzeros(i_count, nk), nk = 1,3)
                 endif
              enddo
           enddo
        enddo

        number_pair = 0
        number_pair_wbl = 0
        number_pair_hamil_l = 0
        number_pair_hamil_r = 0
        do i_a = 1, number_nonzeros
           i_hierarchy = rho_nonzeros(i_a,1)
           tier = nth_tier(i_hierarchy)
           ni = rho_nonzeros(i_a,2)
           nj = rho_nonzeros(i_a,3)

           if ( tier > 0 ) then
              do i_tier = 1, tier
                 i_mode = index_minus(i_hierarchy, i_tier, 1)
                 i_elec=modes(i_mode,2) + 1
                 i_sign=modes(i_mode,4)
                 i_hierarchy_minus = index_minus(i_hierarchy, i_tier, 2)
                 i_sym_on = index_minus(i_hierarchy, i_tier, 4)

                 if ( i_sym_on == 0 ) then
                    do nk = 1, dimen_el
                       i_b = rho_sparsity(nk,nj,i_hierarchy_minus)
                       op_value = d_operators(i_sign, i_elec, ni, nk)
                       if ( (i_b .ne. -1) .and. (op_value .ne. 0) ) number_pair = number_pair + 1

                       i_b = rho_sparsity(ni,nk,i_hierarchy_minus)
                       op_value = d_operators(i_sign, i_elec, nk, nj)
                       if ( (i_b .ne. -1) .and. (op_value .ne. 0) ) number_pair = number_pair + 1
                    enddo
                 elseif ( i_sym_on == 1 ) then
                    do nk = 1, dimen_el
                       i_b = rho_sparsity(nj,nk,i_hierarchy_minus)
                       op_value = d_operators(i_sign, i_elec, ni, nk)
                       if ( (i_b .ne. -1) .and. (op_value .ne. 0) ) number_pair = number_pair + 1

                       i_b = rho_sparsity(nk,ni,i_hierarchy_minus)
                       op_value = d_operators(i_sign, i_elec, nk, nj) 
                       if ( (i_b .ne. -1) .and. (op_value .ne. 0) ) number_pair = number_pair + 1
                    enddo
                 endif
              enddo !i_tier

           endif ! tier>0

           do nk = 1, dimen_el
              i_b = rho_sparsity(nk,nj,i_hierarchy)
              if ( (i_b .ne. -1) .and. (hamil(ni,nk) .ne. 0.) ) number_pair_hamil_l = number_pair_hamil_l + 1

              i_b = rho_sparsity(ni,nk,i_hierarchy)
              if ( (i_b .ne. -1) .and. (hamil(nk,nj) .ne. 0.) ) number_pair_hamil_r = number_pair_hamil_r + 1
           enddo

           do nk = 1, dimen_el
           do nl = 1, dimen_el
              i_b = rho_sparsity(nk,nl,i_hierarchy)
              if ( i_b .ne. -1 ) then
                 do i_elec=1,number_el_state
                    op_value = d_operators(0, i_elec, ni, nk) * d_operators(1, i_elec, nl, nj) 
                    if ( op_value .ne. 0 ) number_pair_wbl = number_pair_wbl + 1

                    op_value = d_operators(1, i_elec, ni, nk) * d_operators(0, i_elec, nl, nj) 
                    if ( op_value .ne. 0 ) number_pair_wbl = number_pair_wbl + 1
                 enddo
              endif
           enddo
           enddo

          if ( tier < truncation_tier ) then
             do i_mode = 0, number_mode - 1
                i_elec = modes(i_mode, 2) + 1
                i_sign = 1 - modes(i_mode, 4)
                i_hierarchy_plus = index_plus(i_hierarchy, i_mode, 1)
                i_sym_on = index_plus(i_hierarchy, i_mode, 3)

                if ( i_hierarchy_plus .ne. -1 ) then

                   if ( i_sym_on == 0 ) then
                      do nk = 1, dimen_el
                         i_b = rho_sparsity(nk,nj,i_hierarchy_plus)
                         op_value = d_operators(i_sign, i_elec, ni, nk)
                         if ( (i_b .ne. -1) .and. (op_value .ne. 0) ) number_pair = number_pair + 1

                         i_b = rho_sparsity(ni,nk,i_hierarchy_plus)
                         op_value = d_operators(i_sign, i_elec, nk, nj) 
                         if ( (i_b .ne. -1) .and. (op_value .ne. 0) ) number_pair = number_pair + 1
                      enddo
                   elseif ( i_sym_on == 1 ) then
                      do nk = 1, dimen_el
                         i_b = rho_sparsity(nj,nk,i_hierarchy_plus)
                         op_value = d_operators(i_sign, i_elec, ni, nk)
                         if ( (i_b .ne. -1) .and. (op_value .ne. 0) ) number_pair = number_pair + 1

                         i_b = rho_sparsity(nk,ni,i_hierarchy_plus)
                         op_value = d_operators(i_sign, i_elec, nk, nj) 
                         if ( (i_b .ne. -1) .and. (op_value .ne. 0) ) number_pair = number_pair + 1
                      enddo
                   endif
                endif

             enddo
          endif

        enddo !i_a

end subroutine

subroutine sparse_matrix_elements_b_wbl(modes, nth_tier, indices, index_minus, index_plus, gammas, etas, d_operators, hamil, &
                                rho_sparsity, rho_nonzeros, matrix, matrix_hamil_l, matrix_hamil_r, matrix_wbl, matrix_elem, &
                                matrix_gammas, number_pair, number_pair_hamil_l, number_pair_hamil_r, number_pair_wbl, &
                                number_nonzeros, dimen_el, dimen_hierarchy, number_index_plus, truncation_tier,  number_mode,&
                                number_lead, number_el_state, number_Pade_poles, sign_pm)
        implicit none
        integer, intent(in) :: number_pair, number_pair_hamil_l, number_pair_hamil_r, number_pair_wbl
        integer, intent(in) :: number_nonzeros, dimen_el, dimen_hierarchy, number_index_plus
        integer, intent(in) :: truncation_tier, number_mode
        integer, intent(in) :: number_lead, number_el_state, number_Pade_poles, sign_pm
        integer, intent(in), dimension(0:number_mode-1, 4) :: modes
        integer, intent(in), dimension(0:dimen_hierarchy-1) :: nth_tier
        integer, intent(in), dimension(0:dimen_hierarchy-1, truncation_tier) :: indices
        integer, intent(in), dimension(0:dimen_hierarchy-1, truncation_tier, 4) :: index_minus
        integer, intent(in), dimension(0:number_index_plus-1, 0:number_mode-1, 3) :: index_plus
        complex*16, intent(in), dimension(0:number_lead-1, 0:number_Pade_poles-1, 0:sign_pm-1) :: gammas
        complex*16, intent(in), dimension(0:number_lead-1, 0:number_Pade_poles-1) :: etas
        integer, intent(in), dimension(0:1,number_el_state, dimen_el, dimen_el) :: d_operators
        complex*16, intent(in), dimension(dimen_el, dimen_el) :: hamil
        integer, intent(in), dimension(dimen_el, dimen_el, 0:dimen_hierarchy-1) :: rho_sparsity
        integer, intent(in), dimension(number_nonzeros, 3) :: rho_nonzeros
        integer, intent(out), dimension(number_pair, 6) :: matrix
        integer, intent(out), dimension(number_pair_hamil_l, 4) :: matrix_hamil_l
        integer, intent(out), dimension(number_pair_hamil_r, 4) :: matrix_hamil_r
        integer, intent(out), dimension(number_pair_wbl, 5) :: matrix_wbl
        complex*16, intent(out), dimension(number_pair) :: matrix_elem
        complex*16, intent(out), dimension(number_nonzeros) :: matrix_gammas
        complex*16, parameter :: ci = (0.d0, 1.d0)
        integer :: i_hierarchy, i_hierarchy_minus, i_hierarchy_plus, tier
        integer :: i_mode, i_tier, i_lead, i_elec, i_pole, i_sign, i_sym_on, i_permutation
        complex*16 :: gamma_sum
        integer :: op_value
        real*8 :: sym_button, permutation_button
        integer :: ni, nj, nk, nl, i_count,  i_count_hamil_l, i_count_hamil_r, i_count_wbl   
        integer :: i_a,i_b

        matrix_gammas = 0.d0
        matrix(:,1) = -1
        matrix(:,2) = -1
        matrix(:,3) = 0
        matrix(:,4) = 0
        matrix_hamil_l = 0
        matrix_hamil_r = 0
        matrix_wbl = 0
        matrix_elem = 0.d0
        i_count_hamil_l = 0
        i_count_hamil_r = 0
        i_count_wbl = 0
        i_count = 0
        do i_a = 1, number_nonzeros
           i_hierarchy = rho_nonzeros(i_a,1)
           tier = nth_tier(i_hierarchy)
           ni = rho_nonzeros(i_a,2)
           nj = rho_nonzeros(i_a,3)

           if ( tier > 0 ) then
              do i_tier = 1, tier
                 i_mode = index_minus(i_hierarchy, i_tier, 1)
                 i_lead=modes(i_mode,1)
                 i_elec=modes(i_mode,2)
                 i_pole=modes(i_mode,3)
                 i_sign=modes(i_mode,4)
                 i_hierarchy_minus = index_minus(i_hierarchy, i_tier, 2)
                 i_permutation=index_minus(i_hierarchy,i_tier,3)
                 i_sym_on = index_minus(i_hierarchy, i_tier, 4)

                 permutation_button=(-1.d0)**i_permutation
                 if ( i_sym_on == 0 ) then
                    do nk = 1, dimen_el
                       i_b = rho_sparsity(nk,nj,i_hierarchy_minus)
                       op_value = d_operators(i_sign, i_elec+1, ni, nk)
                       if ( (i_b .ne. -1) .and. (op_value .ne. 0) ) then
                          i_count = i_count + 1
                          matrix(i_count,1) =  i_a
                          matrix(i_count,2) =  i_b
                          matrix(i_count,5) =  i_lead + 1
                          matrix(i_count,6) =  i_elec + 1
                          matrix_elem(i_count) =  - ci*(-1)**(tier-i_tier) * permutation_button * &
                                                       etas(i_lead, i_pole) * op_value  
                       endif

                       i_b = rho_sparsity(ni,nk,i_hierarchy_minus)
                       op_value = d_operators(i_sign, i_elec+1, nk, nj)
                       if ( (i_b .ne. -1) .and. (op_value .ne. 0) ) then
                          i_count = i_count + 1
                          matrix(i_count,1) =  i_a
                          matrix(i_count,2) =  i_b
                          matrix(i_count,4) =  1
                          matrix(i_count,5) =  i_lead + 1
                          matrix(i_count,6) =  i_elec + 1
                          matrix_elem(i_count) =  ci* (-1)**(2*tier-i_tier-1)  * permutation_button * &
                                                       conjg( etas(i_lead, i_pole) ) * op_value
                       endif
                    enddo
                 elseif ( i_sym_on == 1 ) then
                    sym_button = (-1.d0)**floor( (tier-1.d0)/2.d0 )
                    do nk = 1, dimen_el
                       i_b = rho_sparsity(nj,nk,i_hierarchy_minus)
                       op_value = d_operators(i_sign, i_elec+1, ni, nk)
                       if ( (i_b .ne. -1) .and. (op_value .ne. 0) ) then
                          i_count = i_count + 1
                          matrix(i_count,1) =  i_a
                          matrix(i_count,2) =  i_b
                          matrix(i_count,3) =  i_sym_on
                          matrix(i_count,5) =  i_lead + 1
                          matrix(i_count,6) =  i_elec + 1
                          matrix_elem(i_count) =  - ci*(-1)**(tier-i_tier) * permutation_button * &
                                                       etas(i_lead, i_pole) * sym_button * op_value 
                       endif

                       i_b = rho_sparsity(nk,ni,i_hierarchy_minus)
                       op_value = d_operators(i_sign, i_elec+1, nk, nj) 
                       if ( (i_b .ne. -1) .and. (op_value .ne. 0) ) then
                          i_count = i_count + 1
                          matrix(i_count,1) =  i_a
                          matrix(i_count,2) =  i_b
                          matrix(i_count,3) =  i_sym_on
                          matrix(i_count,4) =  1 
                          matrix(i_count,5) =  i_lead + 1
                          matrix(i_count,6) =  i_elec + 1
                          matrix_elem(i_count) =  ci* (-1)**(2*tier-i_tier-1)  * permutation_button * &
                                                       sym_button * conjg( etas(i_lead, i_pole) ) * op_value 
                       endif
                    enddo
                 endif
              enddo !i_tier

              gamma_sum = 0.d0
              do i_tier = 1, tier
                 i_mode = indices(i_hierarchy, i_tier)
                 i_lead = modes(i_mode, 1)
                 i_elec = modes(i_mode, 2)
                 i_pole = modes(i_mode, 3)
                 i_sign = modes(i_mode, 4)
                 gamma_sum = gamma_sum + gammas(i_lead, i_pole, i_sign) 
              enddo 
              matrix_gammas(i_a) = - gamma_sum
           endif ! tier>0

           do nk = 1, dimen_el
              i_b = rho_sparsity(nk,nj,i_hierarchy)
              if ( (i_b .ne. -1) .and. (hamil(ni,nk) .ne. 0.) ) then
                 i_count_hamil_l = i_count_hamil_l + 1
                 matrix_hamil_l(i_count_hamil_l,1) =  i_a
                 matrix_hamil_l(i_count_hamil_l,2) =  i_b
                 matrix_hamil_l(i_count_hamil_l,3) =  ni
                 matrix_hamil_l(i_count_hamil_l,4) =  nk
              endif

              i_b = rho_sparsity(ni,nk,i_hierarchy)
              if ( (i_b .ne. -1) .and. (hamil(nk,nj) .ne. 0.) )then
                 i_count_hamil_r = i_count_hamil_r + 1
                 matrix_hamil_r(i_count_hamil_r,1) =  i_a
                 matrix_hamil_r(i_count_hamil_r,2) =  i_b
                 matrix_hamil_r(i_count_hamil_r,3) =  nk
                 matrix_hamil_r(i_count_hamil_r,4) =  nj
              endif
           enddo

           do nk = 1, dimen_el
           do nl = 1, dimen_el
              i_b = rho_sparsity(nk,nl,i_hierarchy)
              if ( i_b .ne. -1 ) then
                 do i_elec=1,number_el_state
                    op_value = d_operators(0, i_elec, ni, nk) * d_operators(1, i_elec, nl, nj) 
                    if ( op_value .ne. 0 ) then
                       i_count_wbl = i_count_wbl + 1
                       matrix_wbl(i_count_wbl,1) =  i_a
                       matrix_wbl(i_count_wbl,2) =  i_b
                       matrix_wbl(i_count_wbl,3) =  i_elec
                       matrix_wbl(i_count_wbl,4) =  tier
                       matrix_wbl(i_count_wbl,5) =  op_value
                    endif

                    op_value = d_operators(1, i_elec, ni, nk) * d_operators(0, i_elec, nl, nj) 
                    if ( op_value .ne. 0 ) then
                       i_count_wbl = i_count_wbl + 1
                       matrix_wbl(i_count_wbl,1) =  i_a
                       matrix_wbl(i_count_wbl,2) =  i_b
                       matrix_wbl(i_count_wbl,3) =  i_elec
                       matrix_wbl(i_count_wbl,4) =  tier
                       matrix_wbl(i_count_wbl,5) =  op_value
                    endif
                 enddo
              endif
           enddo
           enddo

          if ( tier < truncation_tier ) then
             do i_mode = 0, number_mode - 1
                i_lead = modes(i_mode, 1)
                i_elec = modes(i_mode, 2)
                i_sign = 1 - modes(i_mode, 4)
                i_hierarchy_plus = index_plus(i_hierarchy, i_mode, 1)
                i_permutation = index_plus(i_hierarchy, i_mode, 2)
                i_sym_on = index_plus(i_hierarchy, i_mode, 3)

                permutation_button = (-1.d0)**i_permutation

                if ( i_hierarchy_plus .ne. -1 ) then

                   if ( i_sym_on == 0 ) then
                      do nk = 1, dimen_el
                         i_b = rho_sparsity(nk,nj,i_hierarchy_plus)
                         op_value = d_operators(i_sign, i_elec+1, ni, nk)
                         if ( (i_b .ne. -1) .and. (op_value .ne. 0) ) then
                            i_count = i_count + 1
                            matrix(i_count,1) =  i_a
                            matrix(i_count,2) =  i_b
                            matrix(i_count,5) =  i_lead + 1
                            matrix(i_count,6) =  i_elec + 1
                            matrix_elem(i_count) =  - ci * permutation_button * op_value 
                         endif

                         i_b = rho_sparsity(ni,nk,i_hierarchy_plus)
                         op_value = d_operators(i_sign, i_elec+1, nk, nj) 
                         if ( (i_b .ne. -1) .and. (op_value .ne. 0) ) then
                            i_count = i_count + 1
                            matrix(i_count,1) =  i_a
                            matrix(i_count,2) =  i_b
                            matrix(i_count,4) =  1
                            matrix(i_count,5) =  i_lead + 1
                            matrix(i_count,6) =  i_elec + 1
                            matrix_elem(i_count) =  - ci * permutation_button * (-1)**(tier+1)* op_value 
                         endif
                      enddo
                   elseif ( i_sym_on == 1 ) then
                      sym_button = (-1.d0)**floor( (tier+1.d0)/2.d0 )
                      do nk = 1, dimen_el
                         i_b = rho_sparsity(nj,nk,i_hierarchy_plus)
                         op_value = d_operators(i_sign, i_elec+1, ni, nk)
                         if ( (i_b .ne. -1) .and. (op_value .ne. 0) ) then
                            i_count = i_count + 1
                            matrix(i_count,1) =  i_a
                            matrix(i_count,2) =  i_b
                            matrix(i_count,3) =  i_sym_on
                            matrix(i_count,5) =  i_lead + 1
                            matrix(i_count,6) =  i_elec + 1
                            matrix_elem(i_count) =  - ci * permutation_button * sym_button * op_value 
                         endif

                         i_b = rho_sparsity(nk,ni,i_hierarchy_plus)
                         op_value = d_operators(i_sign, i_elec+1, nk, nj) 
                         if ( (i_b .ne. -1) .and. (op_value .ne. 0) ) then
                            i_count = i_count + 1
                            matrix(i_count,1) =  i_a
                            matrix(i_count,2) =  i_b
                            matrix(i_count,3) =  i_sym_on
                            matrix(i_count,4) =  1
                            matrix(i_count,5) =  i_lead + 1
                            matrix(i_count,6) =  i_elec + 1
                            matrix_elem(i_count) =  - ci * permutation_button * (-1)**(tier+1)*  sym_button *  op_value 
                         endif
                      enddo
                   endif
                endif

             enddo
          endif

        enddo !i_a
        
end subroutine


