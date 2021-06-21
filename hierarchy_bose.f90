subroutine aw_number(N_hierarchy,truncation_tier,index_number)
    implicit none
    integer, intent(in) :: truncation_tier,index_number
    integer, intent(out) :: N_hierarchy
    integer :: i
     
    N_hierarchy=1
    do i=1,truncation_tier
        N_hierarchy=N_hierarchy*(index_number+i)
    enddo

    do i=1,truncation_tier
        N_hierarchy=N_hierarchy/i
    enddo
    N_hierarchy=N_hierarchy-1
end
!================================================================      
module index_setting
    implicit none
    integer :: index_number
    integer, dimension(:,:), allocatable,save ::   jn
    integer, dimension(:,:), allocatable,save ::   Npj
    integer, dimension(:,:), allocatable,save ::   Nmj
    integer, dimension(:), allocatable :: jtmp
    !-----------------------------------------------------------------
contains
    recursive subroutine setindex(Nh0,i1,ja,nj)       
    implicit none
    integer ::  Nh0,i1,j1,ja,jb,nj,mj,i2 
    
    if ( i1 < index_number ) then
        i2 = i1 + 1
        mj = nj - ja
        do jb = 0, mj
            jtmp(i1) = jb
            call setindex(Nh0,i2,jb,mj)
        enddo
    endif
    
    if ( i1 == index_number ) then
        do jb = 0, nj-ja
            jtmp(i1) = jb
            do j1 = 1, index_number
                jn(j1,Nh0) = jtmp(j1)
            enddo
            Nh0 = Nh0 + 1
        enddo
    endif
    end subroutine setindex
!--------------------------------------------------------
    recursive subroutine setNpNm(Nh1,Nh2,i1,jt)         
    implicit none
    integer Nh1,Nh2,jt,i1,i2
    
    if ( i1 < (index_number + 1) ) then
        i2 = i1 + 1
        if ( i1 == jt ) then
            call setNpNm(Nh1,Nh2,i2,jt)
        elseif ( jn(i1,Nh1) == jn(i1,Nh2) ) then
            call setNpNm(Nh1,Nh2,i2,jt)
        endif
    else
        if ( jn(jt,Nh1) == (jn(jt,Nh2)-1) ) Npj(jt,Nh1) = Nh2
        if ( jn(jt,Nh1) == (jn(jt,Nh2)+1) ) Nmj(jt,Nh1) = Nh2
    endif
    end subroutine setnpnm
!--------------------------------------------------------
end module index_setting
!================================================================      
subroutine set_aw(indices,indices_plus,indices_minus,N_index,truncation_tier,N_hierarchy)
    use index_setting
    implicit none
    integer  :: i1,ja,nj,jt,Nh1,Nh2,Nh0
    integer, intent(in) :: N_index,truncation_tier,N_hierarchy
    integer, dimension(N_index,0:N_hierarchy), intent(out) :: indices, indices_plus, indices_minus
    integer :: i,j
    
    allocate( jn(N_index,0:N_hierarchy) )
    allocate( Npj(N_index,0:N_hierarchy) )
    allocate( Nmj(N_index,0:N_hierarchy) )
    allocate( jtmp(N_index) )
    index_number = N_index
    
    jn = 0
    Nh0 = 0
    i1 = 1
    ja = 0
    nj = truncation_tier
    call setindex(Nh0,i1,ja,nj)
    indices = jn
    
    Npj = -1
    Nmj = -1
    do jt = 1, index_number
        do Nh1 = 0, N_hierarchy
            do Nh2 = 0, N_hierarchy
                call setNpNm(Nh1,Nh2,i1,jt)
            enddo
        enddo
    enddo
    indices_plus = Npj
    indices_minus = Nmj
    
    deallocate( jn )
    deallocate( Npj )
    deallocate( Nmj )
    deallocate( jtmp )
end subroutine set_aw
