# Documentation for updating AmberMini

AmberMini is a stripped down version of the AmberTools resources.
Its not meant to be a standalone package but we feel that having the tools 
available a a conda-installable package helps spread usability.

This page is suppose to help people when updating AmberMini as its not exactly 
a smooth or easy process. This not perfect or complete, but its a start.

This process relies on a tool called Meld to easily compare differences 
and merge between files. In reality, any diff program will work, but you 
have to make selective choices on every difference.

Keep in mind, the structure of AmberTools may change with each version, 
so use your judgement on these files. Not every file is needed, 
compare the versions, see if new dependencies are needed. 

You may also 
find it helpful to compare a 3-way diff between the version of 
AmberTools(X), AmberTools(X+1), and AmberMini to see all the changes 
and make sure you dont remove files AmberMini needs that AmberTools does 
not (X is current version AmberMini is built against and you are 
upgrading from).

There are a rules-of-thumb you can think of when choosing files in the 
3-way diff:

| In AmberTools(X) | In AmberTools(X+1) | In AmberMini | Suggestion |
| :---: | :---: | :---: | :--- |
|       |       |   X   | Its probably a critical AmberMini File, check it, updated by hand as needed |
|   X   |       |   X   | Probably removed in AmberTools X+1, probably remove |
|       |   X   |   X   | Not sure under what condition this would happen |
|   X   |   X   |   X   | File updated, incorporate new changes, watch for custom lines |
|   X   |       |       | Nothing to do
|       |   X   |       | Possibly a new file to add, check the folders its in, and see if there are any calls to it from existing files in AmberMini
|   X   |   X   |       | Nothing to do, probably. There is a chance changes to existing files may now require the file now, be aware if there are errors.

## Critical Notes!

There are several files (*mostly* `Makefile`s) that have custom 
functions to make AmberMini work on Windows and as a standalone! These 
create private compiled binaries needed for the standalone 
functionality.

**DO NOT OVERWRITE THESE FUNCTIONS!**

Incorporate them!

### Several library files are BINARIES somehow

There are several files, especially in the `dat` folder which got flagged 
as binaries somewhere in the life of AmberTools. If you are using a `diff` 
program which filters out binaries, you **must turn this off** to see 
the complete file list. Meld is one such tool, disable this filter.
                          
# File list

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
 * `paramfit` <-> `AmberTools/src/paramfit`
 * `sff`      <-> `AmberTools/src/sff`
 * `sqm`      <-> `AmberTools/src/sqm`
3. Compare the following *folder* making absolutley sure you preserve 
   the private compile calls in the `Makefile`. If need be, you may 
   have to add new functions as AmberTools are developed
 * `antechamber` <-> `AmberTools/src/antechamber`
4. **Keep** the following executables in `antechamber`: 
 * `antechamber`
 * `antechamber.bat`
 * `parmchk`
 * `parmchk.bat`
 * `parmchk2`
 * `parmchk2.bat`
5. Compare the following entries in the `share/amber/dat` folder. This 
   will be the most time consuming step to get right, especially as the 
   FFs get developed.
 * `share/amber/dat/antechamber` <-> `dat/antechamber`
 * `share/amber/dat/dgdb`        <-> `dat/dgdb`
 * `share/amber/dat/leap`        <-> `dat/leap`
 * `share/amber/dat/reslib`      <-> `dat/reslib`
 * `share/amber/dat/slko`        <-> `dat/slko`
6. Single file you need to manually edit and specify new programs to 
   compile as needed. I have not figure out which files need the `_pvt` 
   invoke.
 * `Makefile`
7. Folders I have yet to figure out ¯\\\_(ツ)\\\_/¯
 * `share/amber/dat/charmmff_in_amber`         <-> ???
 * `share/amber/dat/contrib`                   <-> ???
 * `share/amber/dat/leap/pol_solvent_database` <-> ???
 * `share/amber/dat/pixmap`                    <-> ???
 * `share/amber/dat/solvents`                  <-> ???

 
# Creating an executable which needs `$AMBERHOME`

Certain AmberTools programs need the `$AMBERHOME` environment variable 
when compiled. In that case, you need to create a special version of that 
program ending in `_pvt` along with shell scripts to trick the program 
into pointing at the AmberMini's `$AMBERHOME`. The instructions here are 
by no means a catch-all, and more meant as a guideline.

These instructions use `newProg` as the program name that would normally 
be added by AmberTools which we will make a custom version of.

1. In the AmberMini head `Makefile`, add the following:
 * Add to the `PROGS` list `newProg$(SFX)` and `newProg_pvt$(WRAPPER_SFX)`
2. If this is a standalone program (like `antechamber`), add `newProg` 
   as an entry.
3. If this program is a dependency of another (like `parmchk`), you will 
   want to work with that source program's `Makefile`
4. Change to the `newProg` source code directory (this may be a dependant directory):
5. In the `Makefile` which compiles `newProg`:
 * Change the entry name of the compiling directive from `newProg` to `newProg_pvt`
   preserving all other parts of the directive name. e.g.
   ```Makefile
   newProg_pvt$(SFX): 	$(OBJECTS)
	  $(CC) $(CFLAGS) $(AMBERCFLAGS) -o newProg_pvt$(SFX) $(OBJECTS) \
		  $(LDFLAGS) $(AMBERLDFLAGS) $(LM) 
   ```
 * Add the following `if...else` directives: 
   ```Makefile
   ifeq ($(OS),Windows_NT)
   $(BINDIR)/newProg$(SFX):
	  cp -L newProg.bat $(BINDIR)
   else
   $(BINDIR)/paramfit$(SFX):
	  cp -L newProg $(BINDIR)
   endif
   ```
 * Make sure there is a directive that when called compiles the `newProg_pvt`
   directive and calls the `newProg` directive in the `if...else` set.
6. Still in the `newProg` source code directory, create an executable 
   file (`chmod +x`) called `newProg` with contents:
   ```Shell
   #!/bin/sh

   # Shell script to wrap newProg so it knows where to find the necessary data
   # files

   amberhome=`dirname $0`

   # Point AMBERHOME to the location that has dat/leap
   export AMBERHOME="$amberhome/.."

   `dirname $0`/newProg_pvt $@
   ```
 7. Still in the `newProg` source code directory, create an executable 
   file (`chmod +x`) called `newProg.bat` with contents:
   ```Batchfile
   :: Point AMBERHOME to the location that has dat/leap
   set AMBERHOME=%~dp0\..
   %~dp0\newProg_pvt %*
   ```
 8. If correct, the top level Makefile will know to expect the `_pvt` 
    program along with a wrapper program based on OS. Then the compiling 
    `Makefile` will compile the program as normal under the `_pvt` name 
    and the the correct wrapper script sets an `$AMBERHOME` variable for 
    the program before calling the compiled program. When the program is 
    called, it should point to the wrapper, which then points to the 
    program.