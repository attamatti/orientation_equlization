WHAT THE SCRIPTS DO:

*******************
part_angdist_eq.py
*******************

1) Divides the particles into a user specified number of bins based on their orientations -
        I haven't figured out the best method to determine what this number should be yet
        In the initial tests 240K particle dataset 3000 gave a better reults than 1000 
        
2) It divides the bins into 'good' and 'bad' based on how many paticles each has
        if a bin's particle count is greater than a user specified number of standard deviations above the mean it is labelled 'bad'
        Inital runs were tested with 0 (using the mean value) values larger than 0 didn't seem to have very much effect , negative values (0.25) should work, but make sure the target particle number isn't negative   

3) Then divides each bad bin into n/10 sub-bins
        (n = total number of particles in the original bin)

3) It then begins throwing away particles from the most populated sub-bin(s) until the target number of particles is reached for that bin 
        Target number of particle is the same user specified number of standard deviations above the mean for all bins as above
        It chooses which particles to throw away based on their _rlnMaxValueProbDistribution values
        This should throw away the 'worst' particles keeping the better ones (need to check if this is the right way to be doing this)
        If every sub-bin has reached 1 particle it stops, even if the target number has not been reached
        Varying the numbe rof sub-bins doesn't seem to have much effect

4) It then repeats the process for each of the bad bins amd writes out a new star file

USAGE:

python part_angdist_eq.py <particles starfile> <number of bins> <number of standard deviations>

IE:

python part_angdist_eq.py run_class01_data.star 500 1

OUTPUTS:
plt.png - a histogram of the user specified bins with # of particles in each
BBnnn-fixed.png - for each bad bin a histogram of particle counts of the sub-bins, blue are the particles that were kept and yellow were those that were thrownaway
filtered_starfile.star - the new starfile.


*******************
analyze_bildfile.py
*******************

compares the bild files from different relion reconstructions

USAGE:
python analyze_bildfile.py <file1> <file2> <filen>

input as many files as desired. if they are a series of tests of different paramerters put them in order so the output graph will also be in order

OUTPUT:
bildcomp.png - histogram of the lengths of the bild file cylinders
    an evenly distributed set of views will have the mean centered around 5-10
    a poorly distributed set of views will have the most of the values clustered around 0-1

still thinking about better ways to visualize this...
