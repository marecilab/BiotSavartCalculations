import numpy as np
import scipy.io as sio
import nibabel as nib
import matplotlib.pyplot as plt

# Constants
MU_0 = 4 * np.pi * 10**-7  # Magnetic constant (T·m/A)

# Function to calculate magnetic field from a current-carrying wire
def calculate_magnetic_field(wire_pos, current, voxel_size, bounds, im_size, nslices):
    X = np.linspace(bounds[0], bounds[1], im_size)
    Y = np.linspace(bounds[0], bounds[1], im_size)
    Z = np.linspace(bounds[2], bounds[2] + voxel_size[2] * nslices, nslices)
    X, Y, Z = np.meshgrid(X, Y, Z)
    
    Bx, By, Bz = np.zeros(X.shape), np.zeros(Y.shape), np.zeros(Z.shape)
    
    for idx, pos in enumerate(wire_pos[:-1]):
        wire_diff = wire_pos[idx+1] - wire_pos[idx]
        wireX, wireY, wireZ = wire_diff[0], wire_diff[1], wire_diff[2]
        Rx, Ry, Rz = X - pos[0], Y - pos[1], Z - pos[2]
        R = np.sqrt(Rx**2 + Ry**2 + Rz**2)
        R[R == 0] = np.inf  # Avoid division by zero
        
        Bx += MU_0 / (4 * np.pi) * current * (wireY*Rz - wireZ*Ry) / (R**3)
        By += MU_0 / (4 * np.pi) * current * (wireZ*Rx - wireX*Rz) / (R**3)
        Bz += MU_0 / (4 * np.pi) * current * (wireX*Ry - wireY*Rx) / (R**3)

    return Bx, By, Bz

# Initial setup and data loading
F3 = sio.loadmat('sub004_F3.mat')['smCenter'] * 0.001
F4 = sio.loadmat('sub004_F4.mat')['smCenter'] * 0.001

# Imaging parameters
im_size = 224
nslices = 208
fov = 224
fovz = 208
bounds = [-fov/2, fov/2, -fovz/2]
voxel_size = [fov/im_size, fov/im_size, fovz/nslices]

# Current values for the wires
current_pos = 1.5 * 10**-3  # positive wire current (A)
current_neg = -1.5 * 10**-3  # negative wire current (A)

# Calculate magnetic field for both wires
Bx_pos_wire, By_pos_wire, Bz_pos_wire = calculate_magnetic_field(F3, current_pos, voxel_size, bounds, im_size, nslices)
Bx_neg_wire, By_neg_wire, Bz_neg_wire = calculate_magnetic_field(F4, current_neg, voxel_size, bounds, im_size, nslices)

# Sum the magnetic fields from both wires
Bx = Bx_pos_wire + Bx_neg_wire
By = By_pos_wire + By_neg_wire
Bz = Bz_pos_wire + Bz_neg_wire

# Save results
B_leadZ = Bz

# Save as .mat file
sio.savemat('B_leadZ_test.mat', {'B_leadZ': B_leadZ})

# Save as .nii file
MPRAGE = nib.load('MPRAGE_Transverse.nii')
nii_BleadZ = nib.Nifti1Image(B_leadZ, affine=MPRAGE.affine)
nii_BleadZ.header['xyzt_units'] = 10
nib.save(nii_BleadZ, 'B_leadZ_test.nii.gz')

print('Magnetic field saved as .nii and .mat')

# Optional visualization
plt.figure('color', 'white')
plt.imshow(B_leadZ[:, :, int(nslices/2)], cmap='jet')
plt.axis('off')
plt.show()
