subroutine READ_POLY_GRID(grid)

    use MOD_OPENFOAM_GRID

    implicit none

    type(mesh) :: grid

    !!!! @todo cambiare interfaccia introducendo casi


    select case (grid%reader)

       case('openfoam')
           ! TODO: 
           ! - handle parallelism
           write(*,*) "reading..."
          call READ_MESH_POLY_OPENFOAM(grid) 

     ! LEGACY routines:
     ! they are actually implemented in AeroSPEED
     ! Please check AeroSPEED for the actual implementation, but currently NOT SUPPORTED
     ! case('intersections')
     !     ! this files are generated offline
     !     ! they are directly read and then discarded during the intersections
     !     ! should ease the memory consumption
     !     nb_elem_fluid = 0
     !     
     !     open(20, file='FILES_MESH_FLUID/info.txt')
     !     read(20,*) nb_elem_fluid
     !     read(20,*) hexa_fluid
     !     read(20,*) polyhedra_fluid
     !     read(20,*) tetra_fluid
     !     close(20)
     !     write(*,*) "elements will be read during the intersections"

     ! case('ensight')
     !     call READ_MESH_FLUID_POLYHEDRA_ENSIGHT() !=> 3D


   end select


end subroutine
