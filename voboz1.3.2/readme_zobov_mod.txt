#### Instructions to run ZOBOV in a non-cubic set of points and without periodic boundary condtions. #####

The input file should be a ascii-format file with information of the position of the particles to run the tessellation. The file should have the following structure:

- first line: 2 (A,B) integers, one indicating the number of total points to be tessellated (e.g. real galaxies of a survey + mocks in edges and holes) and the second is the number of real points.
- next B lines: the position of the A number of points which correspond to galaxies.
- next (A - B) lines: the positions of the mocks.


1.- run vozisol executable file in voboz1.3.2/bin/ folder.
	
	voboz1.3.2/bin/vozisol arg1 arg2 arg3

- arg1: file prefix (for a .txt file).
- arg2: boxsize (needs to be higher than the maximun cordinate distance).
- arg3: y or n if the output files are need in ascii format.

The outputs are 2 files that contains the volumes (*.vol) and the adjacency of the voronoi cells (*.adj), which are need to run the next step.

2.- run jozov to join cells and form voids.

	voboz1.3.2/bin/jozov arg1 arg2 arg3 arg4

- arg1: v (to find voids) or c (to find clusters).
- arg2: file prefix (same as in 1.-) .
- arg3: Density threshold (0 for no threshold).
- arg4: Border-particle density (0 to run this version).

3 files are generated in this step:

- *.void contains the information of different zones that forms a void. First line is the number of voids. The next N lines (with N the number of voids) has a different number of columns: first column is the number of the void; second column could be 0 or 1; third column is the value that define the ZOBOV voids (if higher is more probable to be a void). If the second column is 1 there will be 3 more columns with the same information.
- *.zone: This files has the number of the void that every particle belong. The first line is the number of particle and voids. The next M lines (where M is the number of particles) is an integer with the number of the void label.
- *.txt contains multiple information of the voids found by ZOBOV, the columns are:
- first line: M particles, N vloidsters (after wathershed)
- next M lines: Void# FileVoid# CoreParticle CoreDens ZoneVol Zone#Part Void#Zones VoidVol Void#Part VoidDensContrast VoidProb



