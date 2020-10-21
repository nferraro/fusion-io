import fpy
import numpy as np

mysim = fpy.sim_data('C1.h5')
print('the available fields are stored in the sim_data.available_fields attribute:')
print(mysim.available_fields)


print('the available series are stored in the sim_data.available_traces')
print(mysim.available_fields)

print('use the interactive  python interpreter (we prefer Ipython) to explore the available series.')
print('most all functions have docstrings available')

print('access a certain field by specifying the field that you want an the time')

mymagneticfield1 =  mysim.get_field('magnetic field', time=0)

#evaluate a field:
fieldValue = mymagneticfield1.evaluate((0,0,0))

#Use python's hyper-efficient array casting to evaluate your function:
FieldMagnitudeGrid= [[np.sqrt(np.sum(mymagneticfield1.evaluate((r, z, 0)))) #value of vector
                            for r in np.linspace(2,4, 100)] for z in np.linspace(-1,1, 100)]
                                                                            #on a 100x100 grid.
#Ooh yeah, that's some efficient python!

#Plotting is now as easy as plt.imshow(FieldMagnitudeGrid)!

print('Information on all functions can be found in docstrings. Type \'mymagneticfield.evaluate?\' to get information on the function call. \n Most all functions have docstrings, and if they do not, this should be corrected.')

print('we also provide a few shorthands to easier access fields. You can for example use \'B\' to get a magnetic field by calling mysim.getfield(\'B\', time=2). available shorthands are described in the docstring of sim_data.field')

print('to remove a sim_data from memory you can explicitly delete it, or it will delete itself
        when the last reference to it is removed. deleting the sim_data object whilst there are
        still fields around can result in unexpected behavior!')

#Access series and diagnostics by:
mydiagnostic = mysim.get_diagnostic('slice times')
#plot by plt.plot(mydiagnostic.diagnostic, mydiagnostic.x_axis) if you use matplotlib

#access time traces:
mytimetrace =  mysim.get_time_traces('bharmonics')
