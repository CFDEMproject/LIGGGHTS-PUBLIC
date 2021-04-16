[LIGGGHTS(R)-TUG WWW Site](http://www.cfdem.com),
[LIGGGHTS(R)-TUG Commands](Section_commands.html#comm)

fix wall/sph/general/gap  command
===============
* * *
Syntax
---------------------

```
fix ID group-ID style wallstyle wallstyleargs model_keywords model_values

```

* ID, group-ID are documented in [fix](fix.html) command
* style = wall/sph/general/gap
* wallstyle = _mesh_
* wallstyle args for wallstyle _mesh_ = _n_meshes_ and _meshes_

> _n_meshes_ value = nm
>> nm = # of meshes (see [fix mesh/surface](fix_mesh_surface.html)) to use
for the wall (positive integer)

> _meshes_ values = meshlist
>> meshlist =  ID(s) of the mesh(es) (see
[fix mesh/surface](fix_mesh_surface.html)) to be used.
These must be defined before

* model_keyword/model_value pairs

> _r0_ value = r<sub>0</sub>
>> r<sub>0</sub> = interaction range of the repulsive force, for consistency
with the boundary conditions r<sub>0</sub> should be a half particle distance
(e.g., compared to the initial particle distance)

> _D_ value = D
>> D = repulsion force coefficient (has the unit of a force), specifically
required to avoid wall penetration due to static fluid pressure, D can be
estimated from the expected fluid pressure P and the smoothing length h by
approximately 10 * P * h<sup>2</sup>

> _gap_ value = g
>> g = maximum gap distance, for SPH particles in gaps between two meshes tighter
than the specified value g the gap treatment will be applied.

> _vwall_ values = vx vy vz
>> vx, vy, vz = components of the wall velocity (optional)

* * *
Examples
---------------------
```
fix sphMeshWall1 all wall/sph/general/gap mesh n_meshes 1 meshes cad1 r0 0.0005 D 0.1 gap 0
fix sphMeshWall2 all wall/sph/general/gap mesh n_meshes 2 meshes cad2a cad2b r0 0.0005 D 0.1 gap 0
fix sphMeshWall3 all wall/sph/general/gap mesh n_meshes 1 meshes cad3 r0 0.0005 D 0.1 gap 0.0006
```

* * *
LIGGGHTS(R)-TUG vs. LIGGGHTS(R)-TUG Info
---------------------
This command is not available in LIGGGHTS(R)-TUG.

* * *
Description
---------------------
Includes all features of [fix wall/sph/general](fix_wall_sph_general.md).

Additionally, this fix applies a specific treatment to SPH particles within
thight gaps between two parallel walls, where the resolution is not sufficient
to obtain correct velocity profiles. This is based on the analytical solution
of the Newtonian, laminar flow between two parallel walls, and yields particle
velocities within unresolved gaps, which correspond to the flow rate
through the gap. This means, in the most extreme case of only 1 layer of SPH particles
between the walls, the obtained particle velocity will be the flow rate equivalent
average velocity. More details and the applied equations are given in
[Eitzlmayr et al. (2014)](#Eitzlmayr2014).

The conditions for the application of this gap treatment are:

* the considered SPH particle is in contact with two meshes (defined by
[fix mesh/surface](fix_mesh_surface.html)), for each of them a
separate fix wall/sph/general/gap is used (e.g. the fixes sphMeshWall1 and
sphMeshWall3 in the above 'Examples'). For each of these fixes a parameter g
can be specified by the keyword _gap_. If the values of g are different,
the order of these fixes in the input script is relevant (see below).

* the distances of the considered SPH particle to the contact points of both meshes are
r<sub>1</sub> and r<sub>2</sub>, respectively. For the application of the gap treatment
the sum r<sub>1</sub> + r<sub>2</sub> (i.e., the distance of the parallel walls) has to
be smaller than the parameter g specified for the fix wall/sph/general/gap listed in
the input script secondly (in this example fixMeshWall3, g = 0.0006).

This means for the fixes defined in the above 'Examples', the gap treatment will not
be applied between the meshes cad1 and cad2a-cad2b due to g = 0, but between the meshes
cad1/cad3 as well as cad2a-cad2b/cad3 the gap treatment will be applied if the gap
distance is smaller than the specified g = 0.0006.

If these conditions are not fulfilled, the single wall interaction in the same way
as for [fix wall/sph/general](fix_wall_sph_general.md) will be applied.

If more than one mesh is specified in a single fix wall/sph/general/gap
(e.g., cad2a and cad2b in the above 'Examples'), for these meshes only one wall contact
will be calculated, and the gap treatment will not be applied between these meshes.

* * *
Restart, fix_modify, output, run start/stop, minimize info
---------------------
No information is written to [binary restart files](restart.html).
None of the [fix_modify](fix_modify.html) options are relevant for this fix.
No parameter of this fix can be used with the start/stop keywords of the
[run](run.html) command. This fix is not invoked during
[energy minimization](minimize.html).

This fix generates the following fixes:

* fix property/atom with the ID 'wallContact_id', where 'id'
is the ID of this fix wall/sph/general. wallContact_id stores a vector including
contact data for each particle, namely the following 7 scalars: delta, deltax,
deltay, deltaz, vx, vy, vz. deltax, delty, deltaz are the components of the vector
from the SPH particle to the wall contact point (which is the closest point of the wall).
delta is the magnitude of this vector (the distance from the wall). vx, vy, vz are
the components of the wall velocity at the contact point.

* fix property/atom with the ID 'F_id', where 'id' is the ID of
this fix wall/sph/general. F_id stores for each particle the three components of
the wall force exerted on this particle. This can be used to evaluate the integral
force and torque exerted on the wall by using [compute/reduce](compute_reduce.html).

* fix property/atom with the ID 'wallCount'. This fix is generated only once, also if
fix wall/sph/general/gap is used more than once. wallCount stores for each particle,
how much fix wall/sph/general/gap the particle was in contact with
(wallCount > 1 is a precondition for the application of the gap treatment).

* fix property/atom with the ID 'usedGapmodel'. This fix is also generated only once
and stores a scalar which indicates if the gap treatment was applied to the
considered particle (usedGapmodel = 1) or not (usedGapmodel = 0).

* * *
Restrictions
---------------------
When using style _mesh_, you have to use the style _bin_ for
the [neighbor command](neighbor.html).

Style _mesh_ can not be used in conjunction with triclinic simulation boxes.

If you use a particle insertion command, you have to keep a minium distance
of r<sub>0</sub> between the wall and the insertion region, otherwise newely inserted
SPH particles can be too close to the wall, which can cause instabilities
since the wall repulsion is inversely proportional to the distance from the wall.

If [fix sph/density/summation](fix_sph_density_summation.html), or
[fix sph/velgrad](fix_sph_velgrad.md) are used, the fix wall/sph/general
command has to be used subsequently in the input script, otherwise the
required wall boundary contribution are not accounted for these fixes.

The total number of used fix wall/sph/general/gap is limited to three.

This fix also needs a [fix sph/integrity](fix_sph_integrity.md).

* * *
Related commands
---------------------
[fix wall/sph/general](fix_wall_sph_general.md)
[fix mesh/surface](fix_mesh_surface.html), [fix_move_mesh](fix_move_mesh.html),
[pair_style sph/artVisc/tensCorr](pair_sph_artvisc_tenscorr.md),
[pair_style sph/morris/tensCorr](pair_sph_morris_tenscorr.md),
[fix sph/density/summation](fix_sph_density_summation.html),
[fix sph/density/continuity](fix_sph_density_continuity.html),
[fix sph/velgrad](fix_sph_velgrad.md)
[fix sph/integrity](fix_sph_integrity.md)

* * *
Default
---------------------
none

* * *
<a name="Eitzlmayr2014"/>
**(Eitzlmayr et al., 2014)** A. Eitzlmayr, G. Koscher, J. Khinast,
Simulation of Mixing in Completely and Partially Filled Co-Rotating
Twin-Screw Extruders via SPH. Paper Draft, 2014.
