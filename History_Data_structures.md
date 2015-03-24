# History and Data structures #

The development of this software depend also on the understanding of the developers on the simulation output and on the tools available.
The library "polymer\_lib.h" was born to analyze trajectories obtained with Hoomd-Blue.  Such software create "dcd" files.

To read and analyze these files we were initially converting them with an home made script that was using VMD.

Such script (in tcl/tk) is the following:
```
set file [open output.dat w]
set rings [atomselect 0 "all"]
set nfram [molinfo 0 get numframes]
set natoms [$rings num]
puts $file "Num_of_atoms $natoms"
puts $file "Num_of_frames $nfram"
for {set j 0} {$j < $nfram } { incr j }  {
  $rings frame $j
  $rings update
  set coords [$rings get {x y z}  ]
  puts $file "$coords"
}
close $file
```

If we have the script in the file "dcd\_to\_coord.tcl", the configuration file of the simulation in "conf.xml" and the trajectory in "traj.dcd" we can call vmd as following:

vmd -hoomd conf.xml traj.dcd -dispdev text -eofexit <dcd\_to\_coord.tcl

With that we were obtaining a much readable file called "output.dat".
Using this file and C++ routine we can obtain the number of atoms and frames automatically using the routine "frames\_atoms\_size".
With these numbers we can allocate a one dimensional array to store the data in the same C++ routine.

It is then possible to read the file and store it in the 1D array (the routine "load\_matrix" do that).
The file conserve the index of the beads and the number of frames such that, to find the coordinates of the bead "k"  of the frame "j" with a total number of beads "N" (considering 3 coordinates for each bead) one have to point to the array location:
j\*N\*3 + k\*3 (to get the X coordinate), j\*N\*3 + k\*3+1 (Z), j\*N\*3 + k\*3+2(Z).


Usually, in polymer\_lib.h, the pointer `*data` is used.
Such pointer is the one dimensional array that contain the trajectory of the simulation.