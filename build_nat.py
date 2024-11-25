###################################################################
# Garrett Akre
# Open MC CANDU fuel flux/e ergy comparision
#build_nat.p 
#
# simulates a candu reactor core and loads it with natural uranium
#
#
# Questions? email gakre04@gamil.com
####################################################################

import openmc
import math

#Define Materials

#defines fuels
fuel = openmc.Material(name='fuel')
fuel.add_nuclide('U238', 1.0)
fuel.set_density('g/cm3', 10.0)

#defines heavy water
Heavy_water = openmc.Material(name='heavy water')
Heavy_water.add_nuclide('H2', 1.996)
Heavy_water.add_nuclide('O16', 1.0)
Heavy_water.add_nuclide('H1',0.004)
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
Zirc2_5.add_nuclide('Zr90',.975*.5515)
Zirc2_5.add_nuclide('Zr91',.975*.112)
Zirc2_5.add_nuclide('Zr92',.975*.171)
Zirc2_5.add_nuclide('Zr94',.975*.174)
Zirc2_5.add_nuclide('Zr96',.975*.028)
Zirc2_5.add_nuclide('Nb93',2.5)
Zirc2_5.set_density('g/cm3',6.44)

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
calan_surface = openmc.ZCylinder(r=13.179/2)


#region for outer calandria
reflect = openmc.ZCylinder(r=337.8, boundary_type='reflective')
Calandria_outer =openmc.ZCylinder(r=379.7, boundary_type='reflective')

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

    f_pin = openmc.ZCylinder(r=1.2243/2, x0=r*math.cos(theta), y0=r*math.sin(theta) )

    c_pin = openmc.ZCylinder(r=1.3081/2,x0=r*math.cos(theta), y0=r*math.sin(theta))
    f_pins.append([f_pin,c_pin])


    theta += 2*math.pi/12

#creates outer ring
theta=0
r= 8.6614/2
while theta < (2* math.pi):


    f_pin = openmc.ZCylinder(r=1.2243/2, x0=r*math.cos(theta), y0=r*math.sin(theta) )
    c_pin = openmc.ZCylinder(r=1.3081/2,x0=r*math.cos(theta), y0=r*math.sin(theta))
    f_pins.append([f_pin,c_pin])

    theta += 2*math.pi/18

#creates a list to store all cells for one fuel assembly
fcells =[]

#creates a fuel cell for each pin we have saved.
for i in f_pins:
    cell = openmc.Cell(fill=fuel,region=-i[0])
    fcells.append(cell)
    cell = openmc.Cell(fill=Zirc4,region= -i[1] & +i[0])
    fcells.append(cell)

water_cell=openmc.Cell(fill=Heavy_water,region=-assembly_surface & +f_pins[0][1]) 
pressure_cell=openmc.Cell(fill=Zirc2_5,region=+assembly_surface & -press_surface)
gas_cell=openmc.Cell(fill=CO2,region=+press_surface &-gas_surface)
calan_cell=openmc.Cell(fill=Zirc4, region=+gas_surface & -calan_surface)
blank_cell =openmc.Cell(fill = Heavy_water)
mod_cell = openmc.Cell(fill = Heavy_water, region = +calan_surface &  +openmc.XPlane(-14.2875) & -openmc.XPlane(14.2875)&  +openmc.YPlane(-14.2875) & -openmc.YPlane(14.2875)      )

#creates a list of cells for 1 fuel pin.
fcells.append(water_cell)
fcells.append(pressure_cell)
fcells.append(gas_cell)
fcells.append(calan_cell)
fcells.append(mod_cell)
print(fcells)


#pin_universe = openmc.Universe(cells=(fuel_cell,water_cell))
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
core_cell.region = +openmc.ZPlane(-10) & -openmc.ZPlane(10) &-reflect    
core_cell.fill=core
core_cell.translation=(0,0,0)


#fills empty space with heavy water.
moderator_cell = openmc.Cell(fill=Heavy_water, region = -reflect)

#creates outer calandria shell. the working boundry of our project
calandria_cell = openmc.Cell(fill = Steel,region =+reflect & -Calandria_outer)

geom = openmc.Geometry([core_cell,moderator_cell,calandria_cell])
geom.export_to_xml()
