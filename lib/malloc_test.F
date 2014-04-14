c
c     --- Program to test malloc models ---
c
      program test
c
c     --- declare array & pointer at top level
c
      integer ia(1)
      pointer(iptr, ia)
c
c     --- allocate array in subroutine
c
      call make(iptr, ia)
      write(6,*) ' top level: assign values to array'
      do 1, i=1,10
	ia(i) = i
    1 continue
c
c     --- reallocate array in subroutine
c
      call make2(iptr, ia)
      write(6,*) ' print original array values'
      write(6,*) (ia(i), i=1,10)
      write(6,*) ' top level: assign values to longer array'
      do 3, i=1,20
	ia(i) = i * 10
    3 continue
      write(6,*) ' print array values'
      write(6,*) (ia(i), i=1,20)
      end
c
c--------------------------------------------
c     --- allocate array; experiment whether 
c         pointer needs to be declared as such
c
      subroutine make(iptr, ia)
      integer ia(1)
c
c     -- HP, SGI R8000, AXP: can't have pointer() in
c        this subroutine. SGI R4400 - pointer() ok.
c
c     pointer(iptr, ia)
      integer iptr
      integer  cmalloc
      external cmalloc

      write(6,*) ' call cmalloc'
      iptr = cmalloc(10 * 4)
      if (iptr.eq.0) stop 'iptr 0'
c     write(6,*) ' cmallocd 40'
c
c     -- HP: crashes if assignment is done here, ok above
c
c     write(6,*) ' subr: assign a value to array'
c     ia(1) = 0
      return
      end
c
c--------------------------------------------
c     --- reallocate array; experiment whether 
c         pointer needs to be declared as such
c
      subroutine make2(iptr, ia)
      integer ia(1)
c     pointer(iptr, ia)
      integer iptr
      integer  crealloc
      external crealloc

c     write(6,*) ' call crealloc'
      iptr = crealloc(iptr, 2*10 * 4)
      if (iptr.eq.0) stop 'crealloc iptr 0'
c     write(6,*) ' creallocd 80'
      return
      end
c--------------------------------------------
