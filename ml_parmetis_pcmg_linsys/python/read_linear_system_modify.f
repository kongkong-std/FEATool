      program main
      implicit none
      integer(kind=8) num_row, num_col, num_nzeros
      integer(kind=8),dimension(:),allocatable::row_index, col_index
      double precision,dimension(:),allocatable:: matrix, rhs, sol

      integer(kind=8) call_times

      integer(kind=8) ii, jj

      call_times = 0
      call read_linear_system_summary(call_times, 
     1    num_row, num_col, num_nzeros)

      allocate(row_index(num_nzeros) )
      allocate(col_index(num_nzeros) )
      allocate(matrix(num_nzeros) )
      allocate(rhs(num_row) )
      allocate(sol(num_row) )

      call READ_LINEAR_SYSTEM(call_times, num_row, num_col, num_nzeros, 
     1    row_index, col_index, matrix, rhs, sol)

C      do ii=1, num_nzeros
C        write(*,*) row_index(ii), col_index(ii), matrix(ii) 
C      enddo 

      end program

! ==========================================================================================

      subroutine read_linear_system_summary(call_times, 
     1    num_row, num_col, num_nzeros)
      implicit none
      integer(kind=8),intent(in)::call_times
      integer(kind=8),intent(out):: num_row, num_col, num_nzeros

      character(len=128) file_name

      write(file_name, 001)"Matrix", call_times
  001 format(A,"_",I5.5,".txt")

      open(unit=999, file=trim(file_name), form="unformatted")
      read(unit=999) num_row
      read(unit=999) num_col
      read(unit=999) num_nzeros
      close(unit=999)
      end subroutine

! ==========================================================================================

      SUBROUTINE READ_LINEAR_SYSTEM(CALL_TIMES,
     1    NUM_ROW, NUM_COL, NUM_NZEROS, 
     2    ROW_INDEX, COL_INDEX, MATRIX, RHS, SOL)

      IMPLICIT NONE
      INTEGER(kind=8),INTENT(IN):: CALL_TIMES
      INTEGER(kind=8),INTENT(INOUT):: NUM_ROW, NUM_COL, NUM_NZEROS

      INTEGER(kind=8),INTENT(OUT):: ROW_INDEX(num_nzeros)
      INTEGER(kind=8),INTENT(OUT):: COL_INDEX(num_nzeros)
      DOUBLE PRECISION,INTENT(OUT)::MATRIX(num_nzeros)
      DOUBLE PRECISION,INTENT(OUT)::RHS(num_row)
      DOUBLE PRECISION,INTENT(OUT)::SOL(num_row)

      CHARACTER(LEN=128) FILE_NAME

      WRITE(FILE_NAME,001)"Matrix", CALL_TIMES
 001  FORMAT(A,"_",I5.5,".txt")

      OPEN(UNIT=999,FILE=TRIM(FILE_NAME), FORM="unformatted")
      READ(UNIT=999)NUM_ROW
      READ(UNIT=999)NUM_COL
      READ(UNIT=999)NUM_NZEROS

      READ(UNIT=999)ROW_INDEX(1:NUM_NZEROS)
      READ(UNIT=999)COL_INDEX(1:NUM_NZEROS)
      READ(UNIT=999)MATRIX(1:NUM_NZEROS)
      CLOSE(UNIT=999)

      WRITE(FILE_NAME,001)"RHS", CALL_TIMES
      OPEN(UNIT=999,FILE=TRIM(FILE_NAME), FORM="unformatted")
      READ(UNIT=999) NUM_ROW
      READ(UNIT=999) RHS(1:NUM_ROW)
      CLOSE(UNIT=999)

      WRITE(FILE_NAME,001)"Solution", CALL_TIMES
      OPEN(UNIT=999,FILE=TRIM(FILE_NAME), FORM="unformatted")
      READ(UNIT=999) NUM_ROW
      READ(UNIT=999) SOL(1:NUM_ROW)
      CLOSE(UNIT=999)

      RETURN
      END SUBROUTINE
