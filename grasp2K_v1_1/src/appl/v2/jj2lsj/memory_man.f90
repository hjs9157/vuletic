MODULE memory_man
   USE vast_kind_param
   INTERFACE alloc
      MODULE PROCEDURE alloc_0i, alloc_1i, alloc_2i, alloc_3i,     &
                       alloc_1r, alloc_2r,                         &
                       alloc_0b, alloc_1b, alloc_2b, alloc_3b
   END INTERFACE

   INTERFACE ralloc
      MODULE PROCEDURE ralloc_0i, ralloc_1i, ralloc_2i, ralloc_3i, &
                       ralloc_1r, ralloc_2r,                       &
                       ralloc_0b, ralloc_1b, ralloc_2b, ralloc_3b
   END INTERFACE ralloc

   INTERFACE dalloc
      MODULE PROCEDURE dalloc_1i, dalloc_2i, dalloc_3i,            &
                       dalloc_1r, dalloc_2r,                       &
                       dalloc_1b, dalloc_2b, dalloc_3b
   END INTERFACE

   CONTAINS
   SUBROUTINE alloc_0i(p, first, last, var, sub ) 
     IMPLICIT NONE
     INTEGER, POINTER, DIMENSION(:) :: p
     INTEGER, INTENT(IN) :: first, last
     CHARACTER(LEN=*), INTENT(IN) :: var, sub
     INTEGER :: error

     IF (ASSOCIATED(p)) then
        Write(6,'(T10,A,A,A,A,A)') 'Variable ',trim(var),' in subroutine ', &
                                trim(sub),  ' already allocated'
        STOP
     ELSE
        allocate(p(first:last),STAT=error)
        IF (error/= 0) then
          Write(6,'(T10,A,A,A,A)') 'Unable to allocate ',trim(var), &
                        ' in subroutine ', trim(sub)
          Write(6, '(T10,A,2I15)') 'Number of elements requested: ',first, last
        STOP
        END IF
     END IF
   END SUBROUTINE alloc_0i

   SUBROUTINE alloc_1i(p, n, var, sub) 
     IMPLICIT NONE
     INTEGER, POINTER, DIMENSION(:) :: p
     INTEGER, INTENT(IN) :: n
     CHARACTER(LEN=*), INTENT(IN) :: var, sub
     INTEGER :: error

     IF (ASSOCIATED(p)) then
        Write(6,'(T10,A,A,A,A,A)') 'Variable ',trim(var),' in subroutine ', &
                                trim(sub),  ' already allocated'
        STOP
     ELSE
        allocate(p(n),STAT=error)
        IF (error/= 0) then
          Write(6,'(T10,A,A,A,A)') 'Unable to allocate ',trim(var), &
                        ' in subroutine ', trim(sub)
          Write(6, '(T10,A,I15)') 'Number of elements requested: ',n
        STOP
        END IF
     END IF
   END SUBROUTINE alloc_1i

   SUBROUTINE alloc_2i(p, n1, n2, var, sub) 
     IMPLICIT NONE
     INTEGER, POINTER, DIMENSION(:,:) :: p
     INTEGER, INTENT(IN) :: n1, n2
     CHARACTER(LEN=*), INTENT(IN) :: var, sub
     INTEGER :: error

     IF (ASSOCIATED(p)) then
        Write(6,'(T10,A,A,A,A,A)') 'Variable ',trim(var),' in subroutine ', &
                                trim(sub),  ' already allocated'
        STOP
     ELSE
        allocate(p(n1,n2),STAT=error)
        IF (error /= 0) then
          Write(6,'(T10,A,A,A,A)') 'Unable to allocate ',trim(var), &
                        ' in subroutine ', trim(sub)
          Write(6, '(T10,A,I15,I15)') 'Number of elements requested: ',n1,n2
        STOP
        END IF
     END IF
   END SUBROUTINE alloc_2i

   SUBROUTINE alloc_3i(p, n1, n2, n3, var, sub)
     IMPLICIT NONE
     INTEGER, POINTER, DIMENSION(:,:,:) :: p
     INTEGER, INTENT(IN) :: n1, n2, n3
     CHARACTER(LEN=*), INTENT(IN) :: var, sub
     INTEGER :: error

     IF (ASSOCIATED(p)) then
        Write(6,'(T10,A,A,A,A,A)') 'Variable ',trim(var),' in subroutine ', &
                                trim(sub),  ' already allocated'
        STOP
     ELSE
        allocate(p(n1,n2,n3),STAT=error)
        IF (error /= 0) then
          Write(6,'(T10,A,A,A,A)') 'Unable to allocate ',trim(var), &
                        ' in subroutine ', trim(sub)
          Write(6, '(T10,A,3I15)') 'Number of elements requested: ',n1,n2,n3
        STOP
        END IF
     END IF
   END SUBROUTINE alloc_3i

   
   SUBROUTINE alloc_1r(p, n, var, sub) 
     IMPLICIT NONE
     real(double), POINTER, DIMENSION(:) :: p
     INTEGER, INTENT(IN) :: n
     CHARACTER(LEN=*), INTENT(IN) :: var, sub
     INTEGER :: error

     IF (ASSOCIATED(p)) then
        Write(6,'(T10,A,A,A,A,A)') 'Variable ',trim(var),' in subroutine ', &
                                trim(sub),  ' already allocated'
        STOP
     ELSE
        allocate(p(n),STAT=error)
        IF (error /= 0) then
          Write(6,'(T10,A,A,A,A)') 'Unable to allocate ',trim(var), &
                        ' in subroutine ', trim(sub)
          Write(6, '(T10,A,I15)') 'Number of elements requested: ',n
        STOP
        END IF
     END IF
   END SUBROUTINE alloc_1r

   SUBROUTINE alloc_2r(p, n1, n2, var, sub) 
     IMPLICIT NONE
     REAL(double), POINTER, DIMENSION(:,:) :: p
     INTEGER, INTENT(IN) :: n1, n2
     CHARACTER(LEN=*), INTENT(IN) :: var, sub
     INTEGER :: error

     IF (ASSOCIATED(p)) then
        Write(6,'(T10,A,A,A,A,A)') 'Variable ',trim(var),' in subroutine ', &
                                trim(sub),  ' already allocated'
        STOP
     ELSE
        allocate(p(n1,n2),STAT=error)
        IF (error /= 0) then
          Write(6,'(T10,A,A,A,A)') 'Unable to allocate ',trim(var), &
                        ' in subroutine ', trim(sub)
          Write(6, '(T10,A,I15,I15)') 'Number of elements requested: ',n1, n2
        STOP
        END IF
     END IF
   END SUBROUTINE alloc_2r

   SUBROUTINE alloc_0b(p, first, last, var, sub ) 
     IMPLICIT NONE
     INTEGER(BYTE), POINTER, DIMENSION(:) :: p
     INTEGER, INTENT(IN) :: first, last
     CHARACTER(LEN=*), INTENT(IN) :: var, sub
     INTEGER :: error

     IF (ASSOCIATED(p)) then
        Write(6,'(T10,A,A,A,A,A)') 'Variable ',trim(var),' in subroutine ', &
                                trim(sub),  ' already allocated'
        STOP
     ELSE
        allocate(p(first:last),STAT=error)
        IF (error/= 0) then
          Write(6,'(T10,A,A,A,A)') 'Unable to allocate ',trim(var), &
                        ' in subroutine ', trim(sub)
          Write(6, '(T10,A,2I15)') 'Number of elements requested: ',first, last
        STOP
        END IF
     END IF
   END SUBROUTINE alloc_0b

   SUBROUTINE alloc_1b(p, n, var, sub) 
     IMPLICIT NONE
     INTEGER(BYTE), POINTER, DIMENSION(:) :: p
     INTEGER, INTENT(IN) :: n
     CHARACTER(LEN=*), INTENT(IN) :: var, sub
     INTEGER :: error

     IF (ASSOCIATED(p)) then
        Write(6,'(T10,A,A,A,A,A)') 'Variable ',trim(var),' in subroutine ', &
                                trim(sub),  ' already allocated'
        STOP
     ELSE
        allocate(p(n),STAT=error)
        IF (error/= 0) then
          Write(6,'(T10,A,A,A,A)') 'Unable to allocate ',trim(var), &
                        ' in subroutine ', trim(sub)
          Write(6, '(T10,A,I15)') 'Number of elements requested: ',n
        STOP
        END IF
     END IF
   END SUBROUTINE alloc_1b

   SUBROUTINE alloc_2b(p, n1, n2, var, sub) 
     IMPLICIT NONE
     INTEGER(BYTE), POINTER, DIMENSION(:,:) :: p
     INTEGER, INTENT(IN) :: n1, n2
     CHARACTER(LEN=*), INTENT(IN) :: var, sub
     INTEGER :: error

     IF (ASSOCIATED(p)) then
        Write(6,'(T10,A,A,A,A,A)') 'Variable ',trim(var),' in subroutine ', &
                                trim(sub),  ' already allocated'
        STOP
     ELSE
        allocate(p(n1,n2),STAT=error)
        IF (error /= 0) then
          Write(6,'(T10,A,A,A,A)') 'Unable to allocate ',trim(var), &
                        ' in subroutine ', trim(sub)
          Write(6, '(T10,A,I15,I15)') 'Number of elements requested: ',n1,n2
        STOP
        END IF
     END IF
   END SUBROUTINE alloc_2b

   SUBROUTINE alloc_3b(p, n1, n2, n3, var, sub)
     IMPLICIT NONE
     INTEGER(BYTE), POINTER, DIMENSION(:,:,:) :: p
     INTEGER, INTENT(IN) :: n1, n2, n3
     CHARACTER(LEN=*), INTENT(IN) :: var, sub
     INTEGER :: error

     IF (ASSOCIATED(p)) then
        Write(6,'(T10,A,A,A,A,A)') 'Variable ',trim(var),' in subroutine ', &
                                trim(sub),  ' already allocated'
        STOP
     ELSE
        allocate(p(n1,n2,n3),STAT=error)
        IF (error /= 0) then
          Write(6,'(T10,A,A,A,A)') 'Unable to allocate ',trim(var), &
                        ' in subroutine ', trim(sub)
          Write(6, '(T10,A,3I15)') 'Number of elements requested: ',n1,n2,n3
        STOP
        END IF
     END IF
   END SUBROUTINE alloc_3b


   SUBROUTINE dalloc_1i(p, var, sub)
     IMPLICIT NONE
     INTEGER, POINTER, DIMENSION(:) :: p
     CHARACTER(LEN=*) :: var, sub
     IF (ASSOCIATED(p)) THEN
        DEALLOCate (p)
     ELSE
        Write(6,'(T10,A,A,A,A,A)') ' DALLOC:  Variable ',trim(var), &
                   ' in subroutine ', trim(sub),  ' not allocated'
     END IF
   END SUBROUTINE dalloc_1i

   SUBROUTINE dalloc_2i(p, var, sub)
     IMPLICIT NONE
     INTEGER, POINTER, DIMENSION(:,:) :: p
     CHARACTER(LEN=*) :: var, sub
     IF (ASSOCIATED(p)) THEN
        DEALLOCate (p)
     ELSE
        Write(6,'(T10,A,A,A,A,A)') ' DALLOC:  Variable ',trim(var), &
                   ' in subroutine ', trim(sub),  ' not allocated'
     END IF
   END SUBROUTINE dalloc_2i

   SUBROUTINE dalloc_3i(p, var, sub)
     IMPLICIT NONE
     INTEGER, POINTER, DIMENSION(:,:,:) :: p
     CHARACTER(LEN=*) :: var, sub
     IF (ASSOCIATED(p)) THEN
        DEALLOCate (p)
     ELSE
        Write(6,'(T10,A,A,A,A,A)') ' DALLOC:  Variable ',trim(var), &
                   ' in subroutine ', trim(sub),  ' not allocated'
     END IF
   END SUBROUTINE dalloc_3i

   SUBROUTINE dalloc_1r(p, var, sub)
     IMPLICIT NONE
     REAL(DOUBLE), POINTER, DIMENSION(:) :: p
     CHARACTER(LEN=*) :: var, sub
     IF (ASSOCIATED(p)) THEN
        DEALLOCate (p)
     ELSE
        Write(6,'(T10,A,A,A,A,A)') ' DALLOC:  Variable ',trim(var), &
                   ' in subroutine ', trim(sub),  ' not allocated'
     END IF
   END SUBROUTINE dalloc_1r

   SUBROUTINE dalloc_2r(p, var, sub)
     IMPLICIT NONE
     REAL(DOUBLE), POINTER, DIMENSION(:,:) :: p
     CHARACTER(LEN=*) :: var, sub
     IF (ASSOCIATED(p)) THEN
        DEALLOCate (p)
     ELSE
        Write(6,'(T10,A,A,A,A,A)') ' DALLOC:  Variable ',trim(var), &
                   ' in subroutine ', trim(sub),  ' not allocated'
     END IF
   END SUBROUTINE dalloc_2r

   SUBROUTINE dalloc_1b(p, var, sub)
     IMPLICIT NONE
     INTEGER(BYTE), POINTER, DIMENSION(:) :: p
     CHARACTER(LEN=*) :: var, sub
     IF (ASSOCIATED(p)) THEN
        DEALLOCate (p)
     ELSE
        Write(6,'(T10,A,A,A,A,A)') ' DALLOC:  Variable ',trim(var), &
                   ' in subroutine ', trim(sub),  ' not allocated'
     END IF
   END SUBROUTINE dalloc_1b

   SUBROUTINE dalloc_2b(p, var, sub)
     IMPLICIT NONE
     INTEGER(BYTE), POINTER, DIMENSION(:,:) :: p
     CHARACTER(LEN=*) :: var, sub
     IF (ASSOCIATED(p)) THEN
        DEALLOCate (p)
     ELSE
        Write(6,'(T10,A,A,A,A,A)') ' DALLOC:  Variable ',trim(var), &
                   ' in subroutine ', trim(sub),  ' not allocated'
     END IF
   END SUBROUTINE dalloc_2b

   SUBROUTINE dalloc_3b(p, var, sub)
     IMPLICIT NONE
     INTEGER(BYTE), POINTER, DIMENSION(:,:,:) :: p
     CHARACTER(LEN=*) :: var, sub
     IF (ASSOCIATED(p)) THEN
        DEALLOCate (p)
     ELSE
        Write(6,'(T10,A,A,A,A,A)') ' DALLOC:  Variable ',trim(var), &
                   ' in subroutine ', trim(sub),  ' not allocated'
     END IF
   END SUBROUTINE dalloc_3b


   SUBROUTINE ralloc_0i(p, first, last, var, sub)
     IMPLICIT NONE
     INTEGER, POINTER, DIMENSION(:) :: p, pnew
     INTEGER, INTENT(IN) :: first, last
     CHARACTER(LEN=*), INTENT(IN) :: var, sub
     INTEGER :: error, nold
     IF (.NOT.  ASSOCIATED(p)) then
        Write(6,'(T10,A,A,A,A,A)') 'Variable ',trim(var),' in subroutine ', &
                                trim(sub),  ' not already allocated'
        STOP
     ELSE
        allocate(pnew(first:last),STAT=error)
        IF (error/= 0) then
          Write(6,'(T10,A,A,A,A)') 'Unable to reallocate ',trim(var), &
                        ' in subroutine ', trim(sub)
          Write(6, '(T10,A,I15)') 'Elements requested for range: ', &
                            first, ' to ', last
        STOP
        END IF
        nold = MIN(SIZE(p), last+1) - 1
        pnew(first:nold) = p(first:nold)
        DEALLOCATE(p)
        p => pnew
     END IF
   END SUBROUTINE ralloc_0i


   SUBROUTINE ralloc_1i(p, n, var, sub) 
     IMPLICIT NONE
     INTEGER, POINTER, DIMENSION(:) :: p, pnew
     INTEGER, INTENT(IN) :: n
     CHARACTER(LEN=*), INTENT(IN) :: var, sub
     INTEGER :: error, nold
     IF (.NOT.  ASSOCIATED(p)) then
        Write(6,'(T10,A,A,A,A,A)') 'Variable ',trim(var),' in subroutine ', &
                                trim(sub),  ' not already allocated'
        STOP
     ELSE
        allocate(pnew(n),STAT=error)
        IF (error/= 0) then
          Write(6,'(T10,A,A,A,A)') 'Unable to reallocate ',trim(var), &
                        ' in subroutine ', trim(sub)
          Write(6, '(T10,A,I15)') 'Number of elements requested: ',n
        STOP
        END IF
        nold = MIN(SIZE(p),n)
        pnew(1:nold) = p(1:nold)
        DEALLOCATE(p)
        p => pnew
     END IF
   END SUBROUTINE ralloc_1i

   SUBROUTINE ralloc_2i(p, n1, n2, var, sub) 
     IMPLICIT NONE
     INTEGER, POINTER, DIMENSION(:,:) :: p, pnew
     INTEGER, INTENT(IN) :: n1, n2
     CHARACTER(LEN=*), INTENT(IN) :: var, sub
     INTEGER :: error, n1_old, n2_old

     IF (.NOT. ASSOCIATED(p)) then
        Write(6,'(T10,A,A,A,A,A)') 'Variable ',trim(var),' in subroutine ', &
                                trim(sub),  ' not already allocated'
        STOP
     ELSE
        allocate(pnew(n1,n2),STAT=error)
        IF (error /= 0) then
          Write(6,'(T10,A,A,A,A)') 'Unable to allocate ',trim(var), &
                        ' in subroutine ', trim(sub)
          Write(6, '(T10,A,I15,I15)') 'Number of elements requested: ',n1,n2
        STOP
        END IF
        n1_old = MIN(SIZE(p, DIM=1),n1)
        n2_old = MIN(SIZE(p, DIM=2),n2)
        pnew(1:n1_old, 1:n2_old) = p(1:n1_old, 1:n2_old)
        DEALLOCATE(p)
        p => pnew
     END IF
   END SUBROUTINE ralloc_2i
   
   SUBROUTINE ralloc_3i(p, n1, n2, n3, var, sub) 
     IMPLICIT NONE
     INTEGER, POINTER, DIMENSION(:,:,:) :: p, pnew
     INTEGER, INTENT(IN) :: n1, n2, n3
     CHARACTER(LEN=*), INTENT(IN) :: var, sub
     INTEGER :: error, n1_old, n2_old, n3_old

     IF (.NOT. ASSOCIATED(p)) then
        Write(6,'(T10,A,A,A,A,A)') 'Variable ',trim(var),' in subroutine ', &
                                trim(sub),  ' not already allocated'
        STOP
     ELSE
        allocate(pnew(n1,n2,n3),STAT=error)
        IF (error /= 0) then
          Write(6,'(T10,A,A,A,A)') 'Unable to allocate ',trim(var), &
                        ' in subroutine ', trim(sub)
          Write(6, '(T10,A,I15,I15)') 'Number of elements requested: ',n1,n2,n3
        STOP
        END IF
        n1_old = MIN(SIZE(p, DIM=1),n1)
        n2_old = MIN(SIZE(p, DIM=2),n2)
        n3_old = MIN(SIZE(p, DIM=3),n3)
        pnew(1:n1_old, 1:n2_old, 1:n3_old) = p(1:n1_old, 1:n2_old, 1:n3_old)
        DEALLOCATE(p)
        p => pnew
     END IF
   END SUBROUTINE ralloc_3i
   
   SUBROUTINE ralloc_1r(p, n, var, sub) 
     IMPLICIT NONE
     real(double), POINTER, DIMENSION(:) :: p, pnew
     INTEGER, INTENT(IN) :: n
     CHARACTER(LEN=*), INTENT(IN) :: var, sub
     INTEGER :: error, nold

     IF (.NOT. ASSOCIATED(p)) then
        Write(6,'(T10,A,A,A,A,A)') 'Variable ',trim(var),' in subroutine ', &
                                trim(sub),  ' not already allocated'
        STOP
     ELSE
        allocate(pnew(n),STAT=error)
        IF (error /= 0) then
          Write(6,'(T10,A,A,A,A)') 'Unable to allocate ',trim(var), &
                        ' in subroutine ', trim(sub)
          Write(6, '(T10,A,I15)') 'Number of elements requested: ',n
        STOP
        END IF
        nold = MIN(SIZE(p),n)
        pnew(1:nold) = p(1:nold)
        DEALLOCATE(p)
        p => pnew
     END IF
   END SUBROUTINE ralloc_1r

   SUBROUTINE ralloc_2r(p, n1, n2, var, sub) 
     IMPLICIT NONE
     REAL(double), POINTER, DIMENSION(:,:) :: p, pnew
     INTEGER, INTENT(IN) :: n1, n2
     CHARACTER(LEN=*), INTENT(IN) :: var, sub
     INTEGER :: error, n1_old, n2_old

     IF (.NOT. ASSOCIATED(p)) then
        Write(6,'(T10,A,A,A,A,A)') 'Variable ',trim(var),' in subroutine ', &
                                trim(sub),  ' not already allocated'
        STOP
     ELSE
        allocate(pnew(n1,n2),STAT=error)
        IF (error /= 0) then
          Write(6,'(T10,A,A,A,A)') 'Unable to allocate ',trim(var), &
                        ' in subroutine ', trim(sub)
          Write(6, '(T10,A,I15,I15)') 'Number of elements requested: ',n1, n2
        STOP
        END IF
        n1_old = MIN(SIZE(p, DIM=1),n1)
        n2_old = MIN(SIZE(p, DIM=2),n2)
        pnew(1:n1_old, 1:n2_old) = p(1:n1_old, 1:n2_old)
        DEALLOCATE(p)
        p => pnew 
     END IF
   END SUBROUTINE ralloc_2r

   SUBROUTINE ralloc_0b(p, first, last, var, sub)
     IMPLICIT NONE
     INTEGER(BYTE), POINTER, DIMENSION(:) :: p, pnew
     INTEGER, INTENT(IN) :: first, last
     CHARACTER(LEN=*), INTENT(IN) :: var, sub
     INTEGER :: error, nold
     IF (.NOT.  ASSOCIATED(p)) then
        Write(6,'(T10,A,A,A,A,A)') 'Variable ',trim(var),' in subroutine ', &
                                trim(sub),  ' not already allocated'
        STOP
     ELSE
        allocate(pnew(first:last),STAT=error)
        IF (error/= 0) then
          Write(6,'(T10,A,A,A,A)') 'Unable to reallocate ',trim(var), &
                        ' in subroutine ', trim(sub)
          Write(6, '(T10,A,I15)') 'Elements requested for range: ', &
                            first, ' to ', last
        STOP
        END IF
        nold = MIN(SIZE(p), last+1) - 1
        pnew(first:nold) = p(first:nold)
        DEALLOCATE(p)
        p => pnew
     END IF
   END SUBROUTINE ralloc_0b

   SUBROUTINE ralloc_1b(p, n, var, sub) 
     IMPLICIT NONE
     INTEGER(BYTE), POINTER, DIMENSION(:) :: p, pnew
     INTEGER, INTENT(IN) :: n
     CHARACTER(LEN=*), INTENT(IN) :: var, sub
     INTEGER :: error, nold
     IF (.NOT.  ASSOCIATED(p)) then
        Write(6,'(T10,A,A,A,A,A)') 'Variable ',trim(var),' in subroutine ', &
                                trim(sub),  ' not already allocated'
        STOP
     ELSE
        allocate(pnew(n),STAT=error)
        IF (error/= 0) then
          Write(6,'(T10,A,A,A,A)') 'Unable to reallocate ',trim(var), &
                        ' in subroutine ', trim(sub)
          Write(6, '(T10,A,I15)') 'Number of elements requested: ',n
        STOP
        END IF
        nold = MIN(SIZE(p),n)
        pnew(1:nold) = p(1:nold)
        DEALLOCATE(p)
        p => pnew
     END IF
   END SUBROUTINE ralloc_1b

   SUBROUTINE ralloc_2b(p, n1, n2, var, sub) 
     IMPLICIT NONE
     INTEGER(BYTE), POINTER, DIMENSION(:,:) :: p, pnew
     INTEGER, INTENT(IN) :: n1, n2
     CHARACTER(LEN=*), INTENT(IN) :: var, sub
     INTEGER :: error, n1_old, n2_old

     IF (.NOT. ASSOCIATED(p)) then
        Write(6,'(T10,A,A,A,A,A)') 'Variable ',trim(var),' in subroutine ', &
                                trim(sub),  ' not already allocated'
        STOP
     ELSE
        allocate(pnew(n1,n2),STAT=error)
        IF (error /= 0) then
          Write(6,'(T10,A,A,A,A)') 'Unable to allocate ',trim(var), &
                        ' in subroutine ', trim(sub)
          Write(6, '(T10,A,I15,I15)') 'Number of elements requested: ',n1,n2
        STOP
        END IF
        n1_old = MIN(SIZE(p, DIM=1),n1)
        n2_old = MIN(SIZE(p, DIM=2),n2)
        pnew(1:n1_old, 1:n2_old) = p(1:n1_old, 1:n2_old)
        DEALLOCATE(p)
        p => pnew
     END IF
   END SUBROUTINE ralloc_2b
   
   SUBROUTINE ralloc_3b(p, n1, n2, n3, var, sub) 
     IMPLICIT NONE
     INTEGER(BYTE), POINTER, DIMENSION(:,:,:) :: p, pnew
     INTEGER, INTENT(IN) :: n1, n2, n3
     CHARACTER(LEN=*), INTENT(IN) :: var, sub
     INTEGER :: error, n1_old, n2_old, n3_old, i
!     INTEGER :: error, n1_old, n2_old, n3_old

     IF (.NOT. ASSOCIATED(p)) then
        Write(6,'(T10,A,A,A,A,A)') 'Variable ',trim(var),' in subroutine ', &
                                trim(sub),  ' not already allocated'
        STOP
     ELSE
        allocate(pnew(n1,n2,n3),STAT=error)
        IF (error /= 0) then
          Write(6,'(T10,A,A,A,A)') 'Unable to allocate ',trim(var), &
                        ' in subroutine ', trim(sub)
          Write(6, '(T10,A,I15,I15)') 'Number of elements requested: ',n1,n2,n3
        STOP
        END IF
        n1_old = MIN(SIZE(p, DIM=1),n1)
        n2_old = MIN(SIZE(p, DIM=2),n2)
        n3_old = MIN(SIZE(p, DIM=3),n3)
        do i = 1,n2_old
           pnew(1:n1_old, i, 1:n3_old) = p(1:n1_old, i, 1:n3_old)
        end do
!        pnew(1:n1_old, 1:n2_old, 1:n3_old) = p(1:n1_old, 1:n2_old, 1:n3_old)
        DEALLOCATE(p)
        p => pnew
     END IF
   END SUBROUTINE ralloc_3b
   
END MODULE memory_man

