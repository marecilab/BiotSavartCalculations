import numpy as np
import scipy.io as sio
import nibabel as nib
import matplotlib.pyplot as plt

"""main_Model_WireStray2
Author: Aryan Lund
Modified: Warren Boschen, 31072024

Returns:
B_leadZ_NAME.mat and B_leadZ_NAME.nii.gz

TO-DO:
- Why does including the u0/4pi Biot-Savart factor produce incorrect results?
- Translation and scaling issues with the input data (I think it's the input anyway).
"""

# Function to calculate magnetic field from a current-carrying wire
def calculate_magnetic_field(wire_pos, current, voxel_size, bounds, im_size, nslices):
    X = np.linspace(bounds[0], bounds[1], im_size)
    Y = np.linspace(bounds[0], bounds[1], im_size)
    Z = np.linspace(bounds[2], bounds[2] + voxel_size[2] * nslices, nslices)
    X, Y, Z = np.meshgrid(X, Y, Z)
    
    Bx, By, Bz = np.zeros(X.shape), np.zeros(Y.shape), np.zeros(Z.shape)

    for idx, pos in enumerate(wire_pos[:-1]):
        wire_diff = abs(wire_pos[idx+1] - wire_pos[idx])
        wireX = wire_diff[0] # x-distance from next centroid to current centroid
        wireY = wire_diff[0] # y-distance from next centroid to current centroid
        wireZ = wire_diff[0] # z-distance from next centroid to current centroid
        Rx = X - pos[0] # x-distance from voxel to centroid
        Ry = Y - pos[1] # y-distance from voxel to centroid
        Rz = Z - pos[2] # z-distance from voxel to centroid
        R = np.sqrt(Rx**2 + Ry**2 + Rz**2)
        R[R == 0] = np.inf  # Avoid division by zero

        Bx += current * (wireY*Rz - wireZ*Ry) / (R**3) # currently excluding prefactor u0/4pi
        By += current * (wireZ*Rx - wireX*Rz) / (R**3)
        Bz += current * (wireX*Ry - wireY*Rx) / (R**3)

    return Bx, By, Bz

# Initial setup and data loading
F3 = sio.loadmat('sub004_F3.mat')['smCenter'] * 0.001
F4 = sio.loadmat('sub004_F4.mat')['smCenter'] * 0.001
OZ = sio.loadmat('sub004_OZ.mat')['smCenter'] * 0.001

# Imaging parameters
im_size = 224
nslices = 208
fov = 224
fovz = 208
bounds = [-fov/2, fov/2, nslices/2]
voxel_size = [fov/im_size, fov/im_size, fovz/nslices]

# Reference points for current-carrying wires
# P3_pos = np.array([
#     [-0.145, 0, 0.075],
#     [-0.145, 0, 0.175],
#     [-0.029, 0, 0.225]
# ])

# P3_neg = np.array([
#     [0.145, 0, 0.075],
#     [0.145, 0, 0.175],
#     [0.029, 0, 0.225]
# ])

# Current values for the wires
current_pos = 1.5*(10**-3)  # positive wire current
current_neg = -1.5*(10**-3)  # negative wire current

# Calculate magnetic field for both wires
Bx_pos_wire, By_pos_wire, Bz_pos_wire = calculate_magnetic_field(F3, current_pos, voxel_size, bounds, im_size, nslices)
Bx_neg_wire, By_neg_wire, Bz_neg_wire = calculate_magnetic_field(F4, current_neg, voxel_size, bounds, im_size, nslices)

# Sum the magnetic fields from both wires
Bx = Bx_pos_wire + Bx_neg_wire
By = By_pos_wire + By_neg_wire
Bz = Bz_pos_wire + Bz_neg_wire

# Save results
# B_leadZ = Bz[:, :, 0]

B_leadZ = Bz

# Save as .mat file
sio.savemat('B_leadZ_test.mat', {'B_leadZ': B_leadZ})

# Save as .nii file
MPRAGE = nib.load('MPRAGE_Transverse.nii')
nii_BleadZ = nib.Nifti1Image(B_leadZ, affine=MPRAGE.affine)
nii_BleadZ.header['xyzt_units'] = 10
nib.save(nii_BleadZ, 'B_leadZ_test.nii.gz')

print('Magnetic field saved as .nii and .mat')

# Displaying the results (optional)
# plt.figure('color', 'white')
# plt.imshow(B_leadZ, cmap='jet')
# plt.axis('off')
# plt.show()

# Explanation
# Library Imports:

# numpy: Numerical operations.
# scipy.io: Handling .mat files.
# nibabel: Handling .nii files.
# matplotlib: Visualization.
# Function to Calculate Magnetic Field:

# calculate_magnetic_field: Computes the magnetic field components (Bx, By, Bz) from a given wire position and current.
# Initial Setup and Data Loading:

# Loads and scales the data from .mat files.
# Imaging Parameters:

# Defines parameters such as image size, field of view, voxel size, etc.
# Reference Points for Wires:

# Defines positions of current-carrying wires.
# Magnetic Field Calculation:

# Calculates the magnetic field for both positive and negative wires.
# Sums the magnetic fields from both wires.
# Save Results:

# Saves the resulting magnetic field as .mat and .nii files.
# Display Results (Optional):

# Visualizes the magnetic field using matplotlib.

