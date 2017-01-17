particle_mesh nGridPoints 10

#particle_data number_particles 10
particle_data number_particles 10  
#coupling json myCoupling
coupling liggghts

#Heat properties
model propertiesThermo heatThermalConductivity_solid 
model propertiesThermo heatCapacity_solid 
model propertiesThermo heatDensity_solid 


#Equations
modelEqn 1DSpherical  heat     BC0 1  BC1 2	# 0 Neumann (fixed flux); 1 Dirichlet (fixed value); 2 Convective

control outputTimeStep 0.01
control timeStep 5e-6
#control run 0.000025 init yes
