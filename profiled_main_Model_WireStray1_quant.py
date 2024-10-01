
import cProfile
import pstats

def main():
    # Your existing code content
    import numpy as np
    import scipy.io as sio
    import nibabel as nib
    import matplotlib.pyplot as plt
    from matplotlib import bezier as bez

    """main_Model_WireStray_quant
    Authors: Aryan Lund and Warren Boschen

    Returns:
    B_leadZ_NAME.mat and B_leadZ_NAME.nii.gz

    TO-DO:
    - Why does including the u0/4pi Biot-Savart factor produce incorrect results?
    - Data fabrication and verification.
    """

    # Function to calculate magnetic field from a current-carrying wire
    def calculate_magnetic_field(wire_pos, current, voxel_size, bounds, im_size, nslices):
        X = np.linspace(bounds[0], bounds[1], im_size)
        Y = np.linspace(bounds[0], bounds[1], im_size)
        Z = np.linspace(bounds[2], bounds[3], nslices) # Strangely written. Why?
        X, Y, Z = np.meshgrid(X, Y, Z)
        
        Bx, By, Bz = np.zeros(X.shape), np.zeros(Y.shape), np.zeros(Z.shape)
        
        
        if current < 0:
            wire_pos = np.flipud(wire_pos) # Necessary for experimental centroid data. Not necessary in fabricated data.

        for idx, pos in enumerate(wire_pos[:-1]):
            # w/n NIFTI file, x refers to left/right, y refers to anterior/posterior, and z refers to superior/inferior
            wire_diff = wire_pos[idx+1] - wire_pos[idx]
            wireX = wire_diff[0] # x-distance from next centroid to current centroid
            wireY = wire_diff[1] # y-distance from next centroid to current centroid
            wireZ = wire_diff[2] # z-distance from next centroid to current centroid
            Rx = X - pos[0] # x-distance from voxel to centroid
            Ry = Y - pos[1] # y-distance from voxel to centroid
            Rz = Z - pos[2] # z-distance from voxel to centroid
            R = np.sqrt(Rx**2 + Ry**2 + Rz**2)
            R[R == 0] = np.inf  # Avoid division by zero

            # 10^-7 refers to u0/4pi
            # 10^3 is converting voxels [mm] to match units of u0/4pi [T-m/A]
            Bx += 10**-7 * current * (wireY*Rz - wireZ*Ry) / (R**3) * 10**3
            By += 10**-7 * current * (wireZ*Rx - wireX*Rz) / (R**3) * 10**3
            Bz += 10**-7 * current * (wireX*Ry - wireY*Rx) / (R**3) * 10**3

        return Bx, By, Bz

    # Initial setup and data loading
    F3 = sio.loadmat('sub004_F3.mat')['smCenter'] # * 0.001 # No need to multiply by scaling factor since the centroid list is voxelwise.
    F4 = sio.loadmat('sub004_F4.mat')['smCenter'] # * 0.001
    OZ = sio.loadmat('sub004_OZ.mat')['smCenter'] # * 0.001

    # Imaging parameters
    im_size = 224
    nslices = 208
    fov = 224
    fovz = 208
    bounds = [0, fov-1, 0, nslices-1]
    voxel_size = [fov/im_size, fov/im_size, fovz/nslices]

    # Fabricate test data
    # test_data = np.zeros((200,3))
    # test_data[:,0] = fov/2
    # test_data[:,1] = np.linspace(1, fov, 200)
    # test_data[:,2] = fovz/2

    # Construct a 1000-point Bezier curve from the provided data
    bez_curve1 = bez.BezierSegment(F3)
    params = np.arange(0, 1, 0.002)
    bez_points1 = np.zeros((len(params), 3))
    for t in range(len(params)):
        bez_points1[t][:] = np.asarray(bez_curve1.point_at_t(params[t]))

    bez_curve2 = bez.BezierSegment(F4)
    params = np.arange(0, 1, 0.001)
    bez_points2 = np.zeros((len(params), 3))
    for t in range(len(params)):
        bez_points2[t][:] = np.asarray(bez_curve2.point_at_t(params[t]))

    # Current values for the wires
    current_pos = 1.5*(10**-3)  # positive wire current
    current_neg = -1.5*(10**-3)  # negative wire current

    # Calculate magnetic field for both wires
    Bx_pos_wire, By_pos_wire, Bz_pos_wire = calculate_magnetic_field(bez_points1, current_pos, voxel_size, bounds, im_size, nslices)
    Bx_neg_wire, By_neg_wire, Bz_neg_wire = calculate_magnetic_field(bez_points2, current_neg, voxel_size, bounds, im_size, nslices)

    # Calculate magnetic field for test data
    # Bx_test_wire, By_test_wire, Bz_test_wire = calculate_magnetic_field(bez_points, current_pos, voxel_size, bounds, im_size, nslices)

    # Sum the magnetic fields from both wires
    Bx = Bx_pos_wire + Bx_neg_wire
    By = By_pos_wire + By_neg_wire
    Bz = Bz_pos_wire + Bz_neg_wire

    # Save results
    B_leadZ = Bz
    # B_leadZ = Bz_test_wire

    # Save as .mat file
    sio.savemat('B_leadZ_F3F4_quant.mat', {'B_leadZ': B_leadZ})

    # Save as .nii file
    MPRAGE = nib.load('MPRAGE_Transverse.nii')
    nii_BleadZ = nib.Nifti1Image(B_leadZ, affine=MPRAGE.affine)
    nii_BleadZ.header['xyzt_units'] = 10
    nib.save(nii_BleadZ, 'B_leadZ_F3F4_quant.nii.gz')

    print('Magnetic field saved as .nii and .mat')

    # def plot_vector_field(Bx, By, Bz, X, Y, Z):
    #     import numpy as np
    #     import matplotlib.pyplot as plt
    #     from mpl_toolkits.mplot3d import Axes3D

    #     # Create a 3D figure
    #     fig = plt.figure()
    #     ax = fig.add_subplot(111, projection='3d')

    #     # Plot the vector field using quiver
    #     ax.quiver(X, Y, Z, Bx, By, Bz, length=0.1, normalize=True)

    #     # Set labels and title
    #     ax.set_xlabel('X')
    #     ax.set_ylabel('Y')
    #     ax.set_zlabel('Z')
    #     ax.set_title('3D Vector Field of Magnetic Field Components Bx, By, Bz')

    #     # Show the plot
    #     plt.show()

    # # Call the function to plot the vector field
    # plot_vector_field(Bx, By, Bz, X, Y, Z)

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

if __name__ == "__main__":
    # Profiling the main function

    profiler = cProfile.Profile()
    profiler.enable()
    
    main()
    
    profiler.disable()
    # Save profiling results to a file
    with open("profile_results.txt", "w") as f:
        ps = pstats.Stats(profiler, stream=f)
        ps.sort_stats('cumtime')  # Sorting by cumulative time
        ps.print_stats()




