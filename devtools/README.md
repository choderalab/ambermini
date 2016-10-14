# Documentation for updating AmberMini

AmberMini is a stripped down version of the AmberTools resources.
Its not meant to be a standalone package but we feel that having the tools 
available a a conda-installable package helps spread usability.

This page is suppose to help people when updating AmberMini as its not exactly 
a smooth or easy process. This not perfect or complete, but its a start.

This process relies on a tool called Meld to easily compare differences 
and merge between files. In reality, any diff program will work, but you 
have to make selective choices on every difference.

Keep in mind, the structure of AmberTools may change with subsequent versions, 
so use your judgement on these files. Not every file is needed is the problem, 
compare the versions, see if new dependencies are needed. 

| -------------------------|-------------------------------------------------|
| The Most Critical Thing! |There are several files (*mostly* `Makefile`s)\  |
|                          |that have custom functions to make AmberMini\    |
|                          |work on Windows and as a standalone! These\      |
|                          |create private compiled binaries needed for the\ |
|                          |standalone functionality.                        |
|--------------------------|--------------------------------------------------|
|                          | DO NOT DELETE THESE FUNCTIONS!\                  |
|                          | Incorperate them!                                |
                          
                          

We compare files 
with the AmberMini on the right and AmberTools on the left (e.g. 
AmberMini <-> AmberTools).

1. Pull down and untar the AmberTools tarbal into its own folder
2. Compare the following *folders* and files therein (These are the 
   easy ones)
 * `arpack`   <-> `AmberTools/src/arpack`
 * `blas`     <-> `AmberTools/src/blas`
 * `cifparse` <-> `AmberTools/src/cifparse`
 * `include`  <-> `AmberTools/src/include`
 * `lapack`   <-> `AmberTools/src/lapack`
 * `lib`      <-> `AmberTools/src/lib`
 * `sff`      <-> `AmberTools/src/sff`
 * `sqm`      <-> `AmberTools/src/sqm`
3. Compare the following *folder* making absolutley sure you preserve 
   the private compile calls in the `Makefile`. If need be, you may 
   have to add new functions as AmberTools are developed
 * `antechamber` <-> `AmberTools/src/antechamber`
4. Compare the following entries in the `share/amber/dat` folder. This 
   will be the most time consuming step to get right, especially as the 
   FFs get developed.
 * `share/amber/dat/antechamber` <-> `dat/antechamber`
 * `share/amber/dat/dgdb`        <-> `dat/dgdb`
 * `share/amber/dat/leap`        <-> `dat/leap`
 * `share/amber/dat/reslib`      <-> `dat/reslib`
 * `share/amber/dat/slko`        <-> `dat/slko`
5. Single file you need to manually edit and specify new programs to 
   compile as needed. I have not figure out which files need the `_pvt` 
   invoke.
 * `Makefile`
6. Folders I have yet to figure out ¯\\\_(ツ)\_/¯
 * `share/amber/dat/charmmff_in_amber`         <-> ???
 * `share/amber/dat/contrib`                   <-> ???
 * `share/amber/dat/leap/pol_solvent_database` <-> ???
 * `share/amber/dat/pixmap`                    <-> ???
 * `share/amber/dat/solvents`                  <-> ???
 