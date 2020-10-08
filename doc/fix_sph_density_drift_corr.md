[LIGGGHTS(R)-TUG WWW Site](http://www.cfdem.com),
[LIGGGHTS(R)-TUG Commands](Section_commands.html#comm)

fix sph/density/drift/corr command
===============
* * *
Syntax
---------------------

```
fix ID group-ID style args

```

* ID, group-ID are documented in [fix](fix.html) command
* style = sph/density/drift/corr
* args = list of arguments

> _density_ value = &#961;<sub>0</sub>
>> &#961;<sub>0</sub> = nominal density of the fluid

> _coeff_ value = C
>> C = coefficient

> _after_ value = step
>> step = time step after which the correction is applied

* * *
Examples
---------------------
```
fix driftcorr all sph/density/drift/corr density 1000 coeff 10 after 5000
```

* * *
LIGGGHTS(R)-TUG vs. LIGGGHTS(R)-TUG Info
---------------------
This command is not available in LIGGGHTS(R)-TUG.

* * *
Description
---------------------
This is an SPH density correction applicable together with the
[fix sph/density/continuity](fix_sph_density_continuity.html)
for liquids in completely filled geometries.

A drawback of the [fix sph/density/continuity](fix_sph_density_continuity.html)
is that numerical errors due to the time integration of the density can accumulate.
In SPH simulations of liquids (weakly compressible) in completely filled geometries
(where the entire mass of the fluid and the volume are strictly constant), this
can cause unphysical values of the average density (i.e., a drift
of the average density versus time towards an unphysical level).

In such cases, this fix can be applied and adds a uniform correction value to
the time derivative of the density, which is proportional to the
difference of the nominal fluid density &#961;<sub>0</sub>
(i.e., the ratio of entire mass and entire volume) and the
actual density averaged over all particles &#961;<sub>av</sub>:

![Eq1](Eqs/fix_sph_density_drift_corr_eq1.jpg)

This correction is added to the time derivative of the density (obtained from
[fix sph/density/continuity](fix_sph_density_continuity.html)) for each
particle, thus, density gradients are not affected.

The argument _after_ specifies a time step, after which the correction will be
applied. Before this step, the correction will not be applied. This can be required
to allow an expansion of the fluid initially in order to fill the geometry completely.

* * *
Restart, fix_modify, output, run start/stop, minimize info
---------------------
No information is written to [binary restart files](restart.html).
None of the [fix_modify](fix_modify.html) options are relevant for this fix.
No parameter of this fix can be used with the start/stop keywords of the
[run](run.html) command. This fix is not invoked during
[energy minimization](minimize.html).

* * *
Restrictions
---------------------
none

* * *
Related commands
---------------------
[fix sph/density/continuity](fix_sph_density_continuity.html)

* * *
Default
---------------------
none
