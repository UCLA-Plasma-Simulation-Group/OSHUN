from oshun1d_ld__and_collisions import *
import pylab
import numpy




sim_metadata = pickle.load( open( "./%s/oshun_metadata.p" % ('.',), "rb" ) )
sim_metadata_cfg = pickle.load( open( "./%s/oshun_metadata__cfg.p" % ('.',), "rb" ) )

output_step 				= sim_metadata.output_interval__time_steps
time_step_size 		 	= sim_metadata.dt   	
output_time_step			= float(output_step) * time_step_size

box_size__x 				= sim_metadata_cfg.system_size_max__global[0]
num_cells__x 				= sim_metadata_cfg.num_cells__global[0]
cell_size__x 				= box_size__x / float( num_cells__x )

# also read in the plasma paramters...
density						= sim_metadata_cfg.species_config[0].density_profile.val	
vth_e							= sim_metadata_cfg.species_config[0].temp_profile.val
wp_e							= np.sqrt( density )

# grab the location and freq of drivers...
x1 = sim_metadata_cfg.wave_driver.x_for_w1
x2 = sim_metadata_cfg.wave_driver.x_for_w2
w1 = sim_metadata_cfg.wave_driver.spatial_echo_w1
w2 = sim_metadata_cfg.wave_driver.spatial_echo_w2

x_global_index_for_w1 =  sim_metadata_cfg.wave_driver.x_global_index_for_w1
x_global_index_for_w2 =  sim_metadata_cfg.wave_driver.x_global_index_for_w2

seperation_distance = np.fabs( x2 - x1)
echo_distance = seperation_distance * w2/(w2-w1)
x_of_echo = x1 + echo_distance

print "w1=%f and is at %f, w2=%f and is at %f and we Looking for echo at x=%f \n\n" % (w1, x1, w2, x2, x_of_echo)
print x_global_index_for_w1
print x_global_index_for_w2
print ""





# first load the saved E-field data....
(e_field_data, freq_axis, k_axis, t_axis, x_axis) = load_packed_E_data(output_time_step = output_time_step, limit_data_to_before_time=60.0, cell_size__x=cell_size__x,  box_size__x = box_size__x)


clf()
imshow( e_field_data , extent=[x_axis[0], x_axis[-1], t_axis[0], t_axis[-1]], aspect='auto', origin='lower', cmap='Purples')
title('Time versus Electric Field')
xlabel('Distance (skindepths)')
ylabel("Time ( $w_p^{-1}$ )")
colorbar()
xlim([0.0,6.0])
savefig('time_vs_e_field.png', dpi=300)


e_slice_index = int(x_of_echo/cell_size__x)
w1_index = int(x1/cell_size__x)

clf()
plot(t_axis,e_field_data[:,e_slice_index])
#plot(t_axis,e_field_data[:,e_slice_index-1])
#plot(t_axis,e_field_data[:,e_slice_index+1])
xlim([10.0, 55.0])
ylim([-6e-10, 6e-10])
title('Electric Field at x=%0.2f versus time\n' % x_of_echo )
xlabel("Time ( $w_p^{-1}$ )")
ylabel('Electric Field Amplitude')
savefig('e_field_spatial_echo.png')

clf()
plot(t_axis,e_field_data[:,e_slice_index], label="Echo at x=%0.2f" % x_of_echo)
plot(t_axis,e_field_data[:,int(x1/cell_size__x)], label="Driver w1=%0.2f at x=%0.2f" % (w1, x1) )
plot(t_axis,e_field_data[:,int(x2/cell_size__x)], label="Driver w2=%0.2f at x=%0.2f" % (w2, x2) )
legend()
xlim([10.0, 55.0])
ylim([-2e-9, 2e-9])
title("Waveforms (at fixed location) versus Time\n")
xlabel("Time ( $w_p^{-1}$ )")
ylabel('Electric Field Amplitude')
savefig('e_field_spatial_echo__w_drivers.png')


freq_data_echo = numpy.fft.fftshift( np.absolute( numpy.fft.fft(e_field_data[:,e_slice_index] ) ) )
clf()
plot( freq_axis, freq_data_echo)
xlim([-6.0, 6.0])
title('FFT (Temporal) of Echo Waveform')
xlabel('Frequency (units of $w_p$)')
ylabel('FFT Magnitude')
savefig('fft_time__echo.png')
xlim([-2.0, 2.0])
savefig('fft_time__echo__zoom.png')