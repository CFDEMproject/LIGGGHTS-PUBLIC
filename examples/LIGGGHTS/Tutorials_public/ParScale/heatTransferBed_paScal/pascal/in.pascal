particle_mesh nGridPoints 10

particle_data number_particles 800 #must be number of particles + 1!
#coupling json myCoupling
coupling liggghts

#Heat properties
model propertiesThermo heatThermalConductivity_solid 
model propertiesThermo heatCapacity_solid 
model propertiesThermo heatDensity_solid 
model propertiesThermo heatTransferCoeff


#Equations
modelEqn 1DSpherical  heat     BC0 1  BC1 2

control outputTimeStep 0.01
control timeStep 1e-5
#control run 0.000025 init yes
