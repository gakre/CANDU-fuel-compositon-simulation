######################################
# Garrett Akre
# OpenMC project
# Models and stimulates neutron flux, energy level and KEff in CANDU reactor
#
#######################################
import openmc
import math
import numpy as np
import matplotlib.pyplot as plt
import openmc.deplete
import os


#Define Materials


#defines enriched fuel
'''
fuel = openmc.Material(name='fuel')
fuel.add_nuclide('U235',0.03)
fuel.add_nuclide('U238',0.97)
fuel.add_nuclide('O16',2)
fuel.set_density('g/cm3', 10.97)
'''

'''
#defines natural fuel
fuel = openmc.Material(name='fuel')
fuel.add_nuclide('U235',0.007)
fuel.add_nuclide('U238',0.993)
fuel.add_nuclide('O16',2)
fuel.set_density('g/cm3', 10.97)
'''

'''
#defines thorium 5%
fuel = openmc.Material(name='fuel')
fuel.add_nuclide('Th232', 1*.95)
fuel.add_nuclide('Pu239', 1*.05)
fuel.add_element('O',2)
fuel.set_density('g/cm3', 10.97)

'''

#defines thorium 1%
fuel = openmc.Material(name='fuel')
fuel.add_nuclide('Th232', 1*.99)
fuel.add_nuclide('Pu239', 1*.01)
fuel.add_element('O',2)
fuel.set_density('g/cm3', 10.97)

#defines heavy water
Heavy_water = openmc.Material(name='heavy water')
Heavy_water.add_nuclide('H2', 2.0)
Heavy_water.add_nuclide('O16', 1.0)
Heavy_water.set_density('g/cm3', 1.07801)


#defiens Zirch 4
Zirc4 = openmc.Material(name='Zirch4')
Zirc4.add_nuclide('Zr90',.985*.5515)
Zirc4.add_nuclide('Zr91',.985*.112)
Zirc4.add_nuclide('Zr92',.985*.171)
Zirc4.add_nuclide('Zr94',.985*.174)
Zirc4.add_nuclide('Zr96',.985*.028)
Zirc4.add_nuclide('Sn112',.014)
Zirc4.add_nuclide('O16',.0012)
Zirc4.add_nuclide('Fe56',.0020)
Zirc4.add_nuclide('Cr52',.001)
Zirc4.set_density('g/cm3',6.56)

#define zirch 2.5
Zirc2_5=openmc.Material(name='Zirch 2.5')
Zirc2_5.add_nuclide('Zr90',.975*.55)
Zirc2_5.add_nuclide('Zr91',.975*.112)
Zirc2_5.add_nuclide('Zr92',.975*.171)
Zirc2_5.add_nuclide('Zr94',.975*.174)
Zirc2_5.add_nuclide('Zr96',.975*.028)
Zirc2_5.add_element('Nb',.025)
Zirc2_5.set_density('g/cm3',.644)

#Defines CO2
CO2=openmc.Material(name="carbon dioxide")
CO2.add_nuclide('C12',1.00)
CO2.add_nuclide('O16',2.00)
CO2.set_density('g/cm3',1.56)

#creates steel
Steel= openmc.Material(name = 'stainless steel')
Steel.add_nuclide('C12',0.02)
Steel.add_nuclide('Mn55',1)
Steel.add_nuclide('Mo98',2)
Steel.add_nuclide('Cr52',17)
Steel.add_nuclide('Ni58',11)
Steel.add_nuclide('Fe56',68.98)

#exports materials
mats = openmc.Materials((fuel, Heavy_water,Zirc4,Zirc2_5,CO2,Steel))
mats.export_to_xml()

#regions for fuel rods
assembly_surface = openmc.ZCylinder(r=5.1689,)
r_pin = openmc.ZCylinder(r=1.2243/2)
c_pin = openmc.ZCylinder(r=1.3081/2)
f_pins =[[r_pin,c_pin]]

#regieons for fuel assembly
press_surface =openmc.ZCylinder(r=11.2064/2,)
gas_surface = openmc.ZCylinder(r=12.90/2,)
calan_surface = openmc.ZCylinder(r=13.179/2,boundary_type='reflective')

#region for fuel assembly end
front_ref = openmc.ZPlane(z0=25, boundary_type='reflective')
end_ref = openmc.ZPlane(z0=-25, boundary_type='reflective')

#region for outer calandria
#reflect = openmc.ZCylinder(r=337.8, boundary_type='reflective')
#Calandria_outer =openmc.ZCylinder(r=379.7, boundary_type='vacuum')

#region for reactor front ends
#front_ref =openmc.ZPlane(z0=297.2,boundary_type='reflective')
#fronr_calan =openmc.ZPlane(z0=297.2+(379.7-337.8),boundary_type='reflective')

#regions for reactor end:
#end_ref =openmc.ZPlane(z0=-297.2,boundary_type='reflective')
#end_calan =openmc.ZPlane(z0=-297.2-(379.7-337.8),boundary_type='vacuum')


#creates first ring
theta=0
r= 2.9769/2
while theta < (2* math.pi):
    f_pin = openmc.ZCylinder(r=1.2243/2, x0=r*math.cos(theta), y0=r*math.sin(theta) )
    c_pin = openmc.ZCylinder(r=1.3081/2,x0=r*math.cos(theta), y0=r*math.sin(theta))
    f_pins.append([f_pin,c_pin])
    theta += 2*math.pi/6

#creates second ring
theta=math.pi/12
r= 5.7506/2
while theta < (2* math.pi):

    f_pin = openmc.ZCylinder(r=1.2243/2, x0=r*math.cos(theta), y0=r*math.sin(theta))
    c_pin = openmc.ZCylinder(r=1.3081/2,x0=r*math.cos(theta), y0=r*math.sin(theta))
    f_pins.append([f_pin,c_pin])
    theta += 2*math.pi/12

#creates outer ring
theta=0
r= 8.6614/2
while theta < (2* math.pi):

    f_pin = openmc.ZCylinder(r=1.2243/2, x0=r*math.cos(theta), y0=r*math.sin(theta))
    c_pin = openmc.ZCylinder(r=1.3081/2,x0=r*math.cos(theta), y0=r*math.sin(theta))
    f_pins.append([f_pin,c_pin])
    theta += 2*math.pi/18


#creates a list to store all cells for one fuel assembly
fcells =[]

#creates a fuel cell for each pin we have saved.
for i in f_pins:
    cell = openmc.Cell(fill=fuel,region=-i[0] &-front_ref &+end_ref)
    fcells.append(cell)
    cell = openmc.Cell(fill=Zirc4,region= -i[1] & +i[0]&-front_ref &+end_ref)
    fcells.append(cell)


water_cell=openmc.Cell(fill=Heavy_water,region=-assembly_surface & +f_pins[0][1] &-front_ref &+end_ref)
pressure_cell=openmc.Cell(fill=Zirc2_5,region=+assembly_surface & -press_surface &-front_ref &+end_ref)
gas_cell=openmc.Cell(fill=CO2,region=+press_surface &-gas_surface &-front_ref &+end_ref)
calan_cell=openmc.Cell(fill=Zirc4, region=+gas_surface & -calan_surface &-front_ref &+end_ref)
blank_cell =openmc.Cell(fill = Heavy_water)
#mod_cell = openmc.Cell(fill = Heavy_water, region = +calan_surface &  +openmc.XPlane(-14.2875) & -openmc.XPlane(14.2875)&  +openmc.YPlane(-14.2875) & -openmc.YPlane(14.2875) &-front_ref &+end_ref)


#creates a list of cells for 1 fuel pin.
fcells.append(water_cell)
fcells.append(pressure_cell)
fcells.append(gas_cell)
fcells.append(calan_cell)
#fcells.append(mod_cell)
#print(fcells)

#uncomment the fallowing to run a simulation on the entire reactor
""" 
#pin_universe = openmc.Universe(cells=(fuel_cells,water_cell))
modr= openmc.Universe(cells=[blank_cell])
fuce = openmc.Universe(cells=fcells)
core = openmc.RectLattice()
core.lower_left = (-314.325,-314.325)
core.pitch=(28.575,28.575)
core.center = (0,0)
core.background = Heavy_water


#sets empty space to heavy water
core.outer = modr


#map the calandria fuel pins
core.universes =[
    [modr,modr,modr,modr,modr,modr,modr,modr,fuce,fuce,fuce,fuce,fuce,fuce,modr,modr,modr,modr,modr,modr,modr,modr,],
    [modr,modr,modr,modr,modr,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,modr,modr,modr,modr,modr,],
    [modr,modr,modr,modr,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,modr,modr,modr,modr,],
    [modr,modr,modr,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,modr,modr,modr,],
    [modr,modr,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,modr,modr,],            
    [modr,modr,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,modr,modr,],
    [modr,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,modr,],
    [modr,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,modr,],
    [fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,],
    [fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,],
    [fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,],
    [fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,],
    [fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,],
    [fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,],
    [modr,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,modr,],
    [modr,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,modr,],
    [modr,modr,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,modr,modr,],            
    [modr,modr,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,modr,modr,],
    [modr,modr,modr,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,modr,modr,modr,],
    [modr,modr,modr,modr,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,modr,modr,modr,modr,],
    [modr,modr,modr,modr,modr,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,fuce,modr,modr,modr,modr,modr,],
    [modr,modr,modr,modr,modr,modr,modr,modr,fuce,fuce,fuce,fuce,fuce,fuce,modr,modr,modr,modr,modr,modr,modr,modr,]
    ]


#applies lattice.
core_cell = openmc.Cell()
core_cell.region = +openmc.ZPlane(-297.2) & -openmc.ZPlane(297.2) &-reflect    
core_cell.fill=core
core_cell.translation=(0,0,0)


#fills empty space with heavy water.
moderator_cell = openmc.Cell(fill=Heavy_water, region = -reflect &-front_ref &+end_ref)


#creates outer calandria shell. the working boundry of our project
calandria_cell = openmc.Cell(fill = Steel,region =+reflect & -Calandria_outer &-front_ref &+end_ref)
calandria_front = openmc.Cell(fill = Steel, region = -Calandria_outer & +front_ref & -fronr_calan)
calandria_end = openmc.Cell(fill = Steel, region = -Calandria_outer & -end_ref & +end_calan)
"""

#exports geometry
#geom = openmc.Geometry([core_cell,moderator_cell,calandria_cell,calandria_front,calandria_end]) # for full reactor

geom = openmc.Geometry(fcells) #for only one cell
geom.export_to_xml()


# Run Tally simulations:


settings = openmc.Settings()
settings.batches = 20
settings.inactive = 3
settings.particles = 1000
settings.export_to_xml()


mesh = openmc.RegularMesh()
mesh.dimension = (300,300)
mesh.lower_left = (-8,-8)
mesh.upper_right = (8,8)


# Create a mesh filter that can be used in a tally
mesh_filter = openmc.MeshFilter(mesh)


# Now use the mesh filter in a tally and indicate what scores are desired
mesh_tally = openmc.Tally(name="Mesh tally")
mesh_tally.filters = [mesh_filter]
mesh_tally.scores = ['fission']

spectrum_tally = openmc.Tally(name='spectrum')
spectrum_tally.filters = [openmc.EnergyFilter(np.geomspace(1e-5, 20.e6, 501)), openmc.ParticleFilter('neutron')]
spectrum_tally.scores = ['flux']



# Instantiate a Tallies collection and export to XML
tallies = openmc.Tallies([mesh_tally, spectrum_tally])
tallies.export_to_xml()
openmc.run()

with openmc.StatePoint('statepoint.20.h5') as sp:
    radial_tally = sp.get_tally(name='Mesh tally')
    spectral_tally = sp.get_tally(name='spectrum')
    energy_source = sp.source['E']
    #heat_source = sp.source['heating']
    #print(f"heat source {all_scores}")
       
data21 = radial_tally.get_reshaped_data(expand_dims=True).squeeze()


# Plot radial flux spectrum(xy)
fig, ax = plt.subplots()
pos = ax.imshow(data21, origin='lower')
cbar = plt.colorbar(pos)
cbar.set_label('Flux [neutron-cm/source]')
ax.set_xlabel('x [cm]')
ax.set_ylabel('y [cm]')
fig.savefig('radial_2D')

# Compare flux spectra
openmc_spectra = spectral_tally.mean.ravel()

energies = spectral_tally.filters[0].values
fig, ax1 = plt.subplots()
ax1.stairs(openmc_spectra, energies, label='Enriched Uranium CANDU')
ax1.grid(True)
ax1.set_title("Distribution of flux energy levels across reactor")
#ax1.set_xlim(1e-3, 0.5)
ax1.set_ylabel('Flux [neutron-cm/source]', fontsize=10)
ax1.legend()
ax1.set_xlabel('Energy [eV]')
fig.savefig('spectrum')


fig2, ax2 = plt.subplots()
# Create log-spaced energy bins from 1 keV to 10 MeV
energy_bins = np.logspace(3,7)

# Calculate pdf for source energies
probability, bin_edges = np.histogram(energy_source, energy_bins, density=True)


# Make sure integrating the PDF gives us unity
print(sum(probability*np.diff(energy_bins)))
log_bin = []


print(f"mean neutron energy {np.mean(energy_source)}")
print(f"mean neutron flux/source {np.mean(data21)}")
# Plot source energy PDF
ax2.set_title("Distribution of neutron energy Natural enrichment")
ax2.grid(True)
ax2.semilogx(energy_bins[:-1], probability*np.diff(energy_bins), color='none')
ax2.step(energy_bins[:-1], probability*np.diff(energy_bins))
ax2.set_xlabel('Energy (eV)')


ax2.set_ylabel('Probability/eV')
fig2.savefig('energy')



#run depletion calculations

fuel.volume = (2*50*math.pi * (1.2243/2)**2*37) #calculates the volume of all fuel elements
model = openmc.Model(geometry=geom, settings=settings)
chain_file = 'chain8.xml' # see openmc data sets for thermal VIII decay chain see openmc libraries for more info
op = openmc.deplete.CoupledOperator(model, chain_file)

power = (153.0e3)
timesteps = [30,30,30,30,30,30,30,30,30,30]  # 12 month simulation
#openmc.deplete.CECMIntegrator(op, timesteps, power, timestep_units='d').integrate()





'''
results = openmc.deplete.Results("depletion_results.h5")
time, keff = results.get_keff()
print(keff)
days =[]
for i in time:
    days.append(i/86400)
print(days)
k=[]
error =[]
for i in keff:
    k.append(i[0])
    error.append(i[1])


fig3,ax3=plt.subplots()
ax3.errorbar(days,k,error,label="natural U")
ax3.grid(True)




ax3.set_title("Keffecive over time considering depletion")
ax3.set_ylabel('Keff', fontsize=10)
ax3.legend()
ax3.set_xlabel('time [days]')
fig3.savefig('keff_overtime_dep')
np.savetxt("keeff_en.txt", np.column_stack((k,error)))
'''
