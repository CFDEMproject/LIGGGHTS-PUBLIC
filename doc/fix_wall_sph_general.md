[LIGGGHTS(R)-TUG WWW Site](http://www.cfdem.com),
[LIGGGHTS(R)-TUG Commands](Section_commands.html#comm)

fix wall/sph/general  command
===============
* * *
Syntax
---------------------

```
fix ID group-ID style wallstyle wallstyleargs model_keywords model_values

```

* ID, group-ID are documented in [fix](fix.html) command
* style = wall/sph/general
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

> _vwall_ values = vx vy vz
>> vx, vy, vz = components of the wall velocity (optional)

* * *
Examples
---------------------
```
fix sphMeshWall all wall/sph/general mesh n_meshes 1 meshes cad1 r0 0.0005 D 0.1 vwall 0 0.1 0
fix sphMeshWalls all wall/sph/general mesh n_meshes 2 meshes cad1 cad2 r0 0.0005 D 0.1
```

* * *
LIGGGHTS(R)-TUG vs. LIGGGHTS(R)-TUG Info
---------------------
This command is not available in LIGGGHTS(R)-TUG.

* * *
Description
---------------------
Bounds the SPH simulation domain with a solid wall, enforcing the no-slip condition.
All particles in the group interact with the wall when they are close enough
(i.e., closer than the interaction range of the used smoothing kernel - note that this
is different from r<sub>0</sub>, a typical value is 2h).

Details about the applied boundary method are provided in
[Eitzlmayr et al. (2014)](#Eitzlmayr2014), specifically Eqs. 15 - 31 and 36 - 44.
This fix calculates the boundary contributions of a solid wall in a generalized
form by fitted polynomials (Eqs. 25, 26, 38, 43, 44) and adds the required
contributions to the SPH neighbor summations, namely in
[fix sph/density/summation](fix_sph_density_summation.html) (Eq. 37),
[fix sph/density/continuity](fix_sph_density_continuity.html) (Eq. 17),
[pair_style sph/artVisc/tensCorr](pair_sph_artvisc_tenscorr.md) (Eq. 42),
[pair_style sph/morris/tensCorr](pair_sph_morris_tenscorr.md) (Eq. 22),
[fix sph/velgrad](fix_sph_velgrad.md) (analogous). Moreover, a repulsion
force is applied if a particle is closer to the wall than the specified
r<sub>0</sub> (Eq. 31).

For wallstyle mesh, the meshlist must be a list of IDs of valid fixes
 of type fix [fix mesh/surface](fix_mesh_surface.html), defining the mesh to be used
for the wall. Triangle-particle neighbor lists are built to efficiently track
particle-triangle contacts. Particle-tri neighbor list build is triggered if any
particle moves more than half the skin distance or (in case of moving mesh) if the
mesh itself moves more than half the skin distance since the last build.
 A warning is generated if a dangerous particle-tri neigh list build
is detected (e.g. if particles are inserted too close to a wall,
see section 'Restrictions').

For a moving wall with style mesh, use [fix move/mesh](fix_move_mesh.html).

If more than one mesh is specified (e.g., cad1 and cad2 in the above 'Examples'),
for these meshes only 1 wall interaction is calculated for each particle, based
on the closest distance.

The optional keyword _vwall_ is required if a mesh should represent a moving
wall, but the movement of the mesh with [fix move/mesh](fix_move_mesh.html)
is not possible due to periodic boundaries. Then, the wall velocity specified
with _vwall_ is applied at the surface of the static mesh.

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
delta is the magnitude of this vector (the distance from the wall). vx, vy, vz are the
components of the wall velocity at the contact point.

* fix property/atom with the ID 'F_id', where 'id' is the ID of
this fix wall/sph/general. F_id stores for each particle the three components of
the wall force exerted on this particle. This can be used to evaluate the integral
force and torque exerted on the wall by using [compute/reduce](compute_reduce.html).

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

* * *
Related commands
---------------------
[fix wall/sph/general/gap](fix_wall_sph_general_gap.md)
[fix mesh/surface](fix_mesh_surface.html), [fix_move_mesh](fix_move_mesh.html),
[pair_style sph/artVisc/tensCorr](pair_sph_artvisc_tenscorr.md),
[pair_style sph/morris/tensCorr](pair_sph_morris_tenscorr.md),
[fix sph/density/summation](fix_sph_density_summation.html),
[fix sph/density/continuity](fix_sph_density_continuity.html),
[fix sph/velgrad](fix_sph_velgrad.md)

* * *
Default
---------------------
none

* * *
<a name="Eitzlmayr2014"/>
**(Eitzlmayr et al., 2014)** A. Eitzlmayr, G. Koscher, J. Khinast,
A novel method for modeling of complex wall geometries in SPH.
Comp. Phys. Comm. 185 (2014) 2436-2448.
