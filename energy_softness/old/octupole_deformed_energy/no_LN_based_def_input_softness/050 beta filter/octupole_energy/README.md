0.5 set contains complementary set of 0.35 cutoff set:

In 0.35 set, the octupole masstable's ground state is determined by the lowest energy among all initial configurations, with an additional condition that the resulting beta2 and beta3 is <= 0.35.
Later on, I discovered the 0.35 cut is removing a lot of valid minimum of beta2 between 0.35 and 0.4/0.5. So I decided to increase the cutoff to 0.5 and investigate the octupole energy/softness.

There's two different types in 0.5 calculation:
1. Those determined to not have octupole deformation (beta3<0.01) because of the cut on beta2<=0.35, eg. nucleus with -100MeV at beta2/3 = 0.36,0.1, with -99MeV at beta2/3 = 0.1,0.0; if we remove all beta2>0.35, this will remove the first nucleus.
 If we increased the cutoff to beta2<=0.5, these nuclei "become" beta3 deformed. So we need to calculate their octupole energy and softness as well, by using a grid of beta2 = -0.5~0.5, 0.05 step; beta3 = 0 (impose reflection symmetry)

2. Those found to be octupole deformed in 0.35 set, we complement the previous calculation (of 0.35 set) that used beta2 = -0.4~0.4, 0.05 step; beta3 = 0 (impose reflection symmetry), with 4 new grid points: beta2 = -0.5, -0.45, 0.45, 0.5; beta3 = 0, so that we don't have to waste resource to calculate the less deformationed configurations again.

After the calculations are completed, we do the following to find the true  beta2 (beta3=0 )minimum:
1. If a nucleus is both in 0.35 set and 0.5 set, we compare the 2 minimas, and select the lower ground state. The reduced data set for control group already has a cut on 0.5, so we don't have to worry about beta2/3>0.5 mixing in.

2. If a nucleus is only in 0.5 set, we simply use this ground state.

After we determined the beta2 (beta3=0 )minimum, we create a set of hfbtho_MASSTABLE.dat using these ground state beta2 values, and beta3 = 0.01, the resulting calculation should be 1 for each nucleus, so the total calculation should only be a few.


