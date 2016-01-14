The POEMS - LIGGGHTS Interface
================================

Preface
-------------
the fix 'fix_poems' contains a point to an object of type 'Workspace', the latter being the top-level class of the POEMS library [Anderson et al](Anderson_2007_POEMS_ParallelizableOpen-Source.pdf). The key functions (within the "Workspace" class) called from LIGGGHTS are (i) MakeSystem, (ii) LobattoOne (during LIGGGHTS' 'initial_integrate' step), and (iii) LobattoTwo (during LIGGGHTS' 'final_integrate' step). 'Lobatto' is a time integration scheme based on the work of [Chun et al.](Chun_MBOND_MultibodyMethodForLont-TimeMolecularDynamics.pdf).

Chun et al. (2000) MBO(N)D: a multibody method for long-time molecular dynamics simulations. J Comput Chem 21(3):159–184.

The class "System"
------------------
This is the main class holding the information of a connected system of rigid bodies. The key data hold is (i) a list of "Body" called 'bodies', and (ii) a list of "Joint" called 'joints'.

The class "Solver"
------------------
This is the main class holding the current state of the system, and advancing the system in time. The only implemented solver at this point is "OnSolver". This solver contains a 'bodyarray' that acts as the temporal object for the solution procedure, and contains a matrix of objects of class 'OnBody'. The 'OnBody' class is hence of central importance for the solution procedure.

The class "Body"
------------------
Key class for describing and solving the problem. It cotains all key information of the body (position, speed, etc.), as well as lists with (i) 'Joint' objects called 'joints, and (ii) 'Point' objects called 'points'.

The class "OnBody"
------------------
Key class for describing and solving the problem. The class contains a pointer to a body object, which is of derived type 'InertialFrame'.

The class "Joint"
------------------
This is a key class for hosting joint information, and has a number of derived classes.

The function 'MakeSystem'
----------------
MakeSystem is used to generate and maintain the description and topology of the multibody system and get the system ready for solving by choosing the desired solver. 

### Check List of bodies and Joints objects. 
Go through and check the input chain list, get the bodies, joints topological information, i.e., number of elements on chain, the head elements of chain, location of joints and reference joints, get 
SystemProcessor.processArray
### Set and prepare system information for solving 
After checking the input chain information, if node exist, set the chain system information for LAMMPS using function in system class: Create_System_LAMMPS  
If isolated unconected bodies found, prepare system for solving using MakeDegenerateSystem function in system class. 
###Function Create_System_LAMMPS
setting up inertial frame, gravity and origin;defining the two handles positions and looping over each body recursive kinematics; creating bodies and setting joints with associated properties. 
### Choose solver type 
One can choose solver type in MakeSystem, i.e., choosing ONsolver or DCAsolver 


The function 'LobattoOne'
--------------

### Assemble Matrix with Forces and Torques
Using a mapping defined in the 'system' data element of the SysData struct, a 6xnumbodies matrix is fille with (i) torques and (ii) forces.

### Solve
This is the most important step in the LobattoOne function, and consist of the below detailed

#### Solver Setup
In this initial step the current system to be solved is handed over to the solver object, and then the function "CreateModel" is executed. The "CreateModel" function is the core function, and is documented in higher [Anderson et al](Anderson_2007_POEMS_ParallelizableOpen-Source.pdf). Essentially, the solver calls 

#### Solver Execution
Then, the selected solver is called TWICE (!), and the state (+1st and 2nd derivative) is then retrieved from the solver object, and assigned to local variables. Then, the solvers is called THREE times (!), and the acceleration is used to perform a half-step and update the velocities. Finally, the position is updated with a full step.

### ForwardKinematics
The kinematics are then forwarded (internal function of the 'system' data element.

### Update the positions, orientation, velocities, and rotation rates of the LIGGGHTS data
In this section the above information is set (and handed over to LIGGGHTS, since the variables are pointers in the function calls).

The function 'LobattoTwo'
--------------
Updating the velocities and forces by the last 1/2 step, using the same methods and functions as 'LobattoOne' does
