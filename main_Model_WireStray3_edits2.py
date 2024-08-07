import numpy as np
import scipy.io as sio
import nibabel as nib
import matplotlib.pyplot as plt

"""main_Model_WireStray3
Author: Aryan Lund
Modified: Warren Boschen, 08062024

Returns:
B_leadZ_NAME.mat and B_leadZ_NAME.nii.gz

TO-DO:
- Why does including the u0/4pi Biot-Savart factor produce incorrect results?
- Translation and scaling issues with the input data (I think it's the input anyway).
- Data fabrication and verification. Seems that the test data created here is found exclusively in the anterior right inferior region.
    Could this be referring to the voxel location in index order and not the true voxel location?
"""

# Function to calculate magnetic field from a current-carrying wire
def calculate_magnetic_field(wire_pos, current, voxel_size, bounds, im_size, nslices):
    X = np.linspace(bounds[0], bounds[1], im_size)
    Y = np.linspace(bounds[0], bounds[1], im_size)
    Z = np.linspace(bounds[2], bounds[3], nslices) 
    X, Y, Z = np.meshgrid(X, Y, Z)
    
    Bx, By, Bz = np.zeros(X.shape), np.zeros(Y.shape), np.zeros(Z.shape)
    
    for idx, pos in enumerate(wire_pos[:-1]):
        wire_diff = wire_pos[idx+1] - wire_pos[idx]
        wireX = wire_diff[0]
        wireY = wire_diff[1]
        wireZ = wire_diff[2]
        Rx = X - pos[0]
        Ry = Y - pos[1]
        Rz = Z - pos[2]
        R = np.sqrt(Rx**2 + Ry**2 + Rz**2)
        R[R == 0] = np.inf  # Avoid division by zero

        # 10^-7 refers to u0/4pi
        # 10^3 is converting voxels [mm] to match units of u0/4pi [T-m/A]
        Bx += 10**-7 * current * (wireY*Rz - wireZ*Ry) / (R**3) * 10**3
        By += 10**-7 * current * (wireZ*Rx - wireX*Rz) / (R**3) * 10**3
        Bz += 10**-7 * current * (wireX*Ry - wireY*Rx) / (R**3) * 10**3

    return Bx, By, Bz

# Initial setup and data loading
F3 = sio.loadmat('sub004_F3.mat')['smCenter']
F4 = sio.loadmat('sub004_F4.mat')['smCenter']
OZ = sio.loadmat('sub004_OZ.mat')['smCenter']

# Imaging parameters
im_size = 224
nslices = 208
fov = 224
fovz = 208
bounds = [0, fov-1, 0, nslices-1]
voxel_size = [fov/im_size, fov/im_size, fovz/nslices]

# Fabricate test data
test_data = np.zeros((100,3))
test_data[:,0] = fov/2
test_data[:,1] = np.linspace(1, fov, 100)
test_data[:,2] = fovz/2

# Current values for the wires
current_pos = 1.5*(10**-3)  # positive wire current
current_neg = -1.5*(10**-3)  # negative wire current

# Calculate magnetic field for test data
Bx_test_wire, By_test_wire, Bz_test_wire = calculate_magnetic_field(test_data, current_pos, voxel_size, bounds, im_size, nslices)

# Load expected theoretical magnetic field (this should be provided or calculated separately)
expected_data = sio.loadmat('expected_magnetic_field.mat')
B_x_expected = expected_data['B_x_expected']
B_y_expected = expected_data['B_y_expected']
B_z_expected = expected_data['B_z_expected']

# Calculate the magnetic field for the provided wire positions
Bx_calculated, By_calculated, Bz_calculated = calculate_magnetic_field(test_data, current_pos, voxel_size, bounds, im_size, nslices)

# Compute the differences between calculated and expected fields
difference_Bx = Bx_calculated - B_x_expected
difference_By = By_calculated - B_y_expected
difference_Bz = Bz_calculated - B_z_expected

# Save the results
B_leadZ = Bz_test_wire
sio.savemat('B_leadZ_test.mat', {'B_leadZ': B_leadZ})

MPRAGE = nib.load('MPRAGE_Transverse.nii')
nii_BleadZ = nib.Nifti1Image(B_leadZ, affine=MPRAGE.affine)
nii_BleadZ.header['xyzt_units'] = 10
nib.save(nii_BleadZ, 'B_leadZ_test.nii.gz')

print('Magnetic field saved as .nii and .mat')

# Visualization of differences (optional)
plt.figure(figsize=(12, 4))
plt.subplot(1, 3, 1)
plt.imshow(difference_Bx[..., 0], cmap='jet')
plt.title('Difference in Bx')

plt.subplot(1, 3, 2)
plt.imshow(difference_By[..., 0], cmap='jet')
plt.title('Difference in By')

plt.subplot(1, 3, 3)
plt.imshow(difference_Bz[..., 0], cmap='jet')
plt.title('Difference in Bz')

plt.show()

# Notes:
# Loading Expected Magnetic Field:
# Ensure expected_magnetic_field.mat contains B_x_expected, B_y_expected, and B_z_expected arrays.
# The dimensions of these arrays should match the calculated fields.
# Visualization:
# Visualization of the differences is optional but useful for debugging and analysis.
# Data:
# Replace the fabricated test data and current values with actual data if available.
# Adjust the imaging parameters and bounds as per the actual dataset used.
# Debugging:
# If the u0/4pi prefactor is causing issues, consider analyzing the input data and scaling factors more closely.






