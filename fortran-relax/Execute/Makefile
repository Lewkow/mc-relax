#	Makefile

SRCS = ../Modules/PHYSICS_MODULE.f90 \
       ../Modules/RAND_SEED_MODULE.f90 \
       ../Modules/INPUTS_MODULE.f90 \
       ../Modules/TABLES_MODULE.f90 \
       ../Modules/MPI_INFO_MODULE.f90 \
       ../Modules/PLANET_MODULE.f90 \
       ../Modules/ENA_MODULE.f90 \
       ../Modules/LISM_MODULE.f90 \
       ../Modules/LISM_CLICK_MODULE.f90 \
       ../Modules/UNIVERSAL_MODULE.f90 \
       ../Modules/SECONDARY_MODULE.f90 \
       ../Modules/FLUX_MAPPING.f90 \
       ../Modules/ESCAPE_TRANS.f90 \
       ../Modules/CLICK_DATA_MODULE.f90 \
       ../Routines/input_reader.f90 \
       ../Routines/rand_angle.f90 \
       ../Routines/HS_rand_angle.f90 \
       ../Routines/HS_tcs.f90 \
       ../Routines/rand_angle_CO2.f90 \
       ../Routines/rand_mfp.f90 \
       ../Routines/read_table.f90 \
       ../Routines/cubic_spline.f90 \
       ../Routines/new_energy.f90 \
       ../Routines/output_writer.f90 \
       ../Routines/plfg.f90 \
       ../Routines/experiments.f90 \
       ../Routines/average_scattering_angle.f90 \
       ../Routines/write_energy_rand_angle.f90 \
       ../Routines/energy_to_velocity.f90 \
       ../Routines/relax_distance_uniform_gas.f90 \
       ../Routines/relax_distance_uniform_gas_no_array.f90 \
       ../Routines/frame_converter.f90 \
			 ../Routines/He_tcs.f90 \
       ../Routines/HH_tcs.f90 \
       ../Routines/Hp_H_tcs.f90 \
       ../Routines/get_tcs.f90 \
       ../Routines/fill_phase_spaces.f90 \
       ../Routines/particle_transport.f90 \
       ../Routines/ion_transport.f90 \
       ../Routines/LISM_density.f90 \
       ../Routines/LISM_temp.f90 \
       ../Routines/LISM_ENA_displacement.f90 \
       ../Routines/LISM_ENA_start.f90 \
       ../Routines/rand_ENA_init_energy.f90 \
       ../Routines/Mars_ENA_start.f90 \
       ../Routines/rand_init_velocity.f90 \
       ../Routines/rand_init_position.f90 \
       ../Routines/rand_init_energy.f90 \
       ../Routines/histogram.f90 \
       ../Routines/vec_histogram.f90 \
       ../Routines/get_click_histograms.f90 \
       ../Routines/energy_loss.f90 \
       ../Routines/electron_density.f90 \
       ../Routines/distance_finder.f90 \
       ../Routines/call_tests.f90 \
       ../Routines/H_CO2_tcs.f90 \
       ../Routines/CX_cross_sections.f90 \
       ../Routines/mars_density.f90 \
       ../Routines/mars_gravity.f90 \
       ../Routines/secondary_hot.f90 \
       ../Routines/ENA_production.f90 \
       ../Routines/chebyshev.f90 \
       ../Routines/erff.f90 \
       ../Routines/universal_rand_angle.f90 \
       ../Routines/universal_tcs.f90 \
       ../Routines/universal_dcs.f90 \
       ../Routines/universal_core.f90 \
       ../Routines/universal_min_angle.f90 \
       ../Routines/write_run_info.f90 \
       ../Routines/planet_3d.f90 \
       ../Routines/planet_onestep_3d.f90 \
       ../Routines/planet_transport.f90 \
       ../Routines/planet_onestep_transport.f90 \
       ../Routines/planet_onestep_ion_transport.f90 \
       ../Routines/reflection.f90 \
       ../Routines/relax.f90 

OBJS = $(SRCS:.f90=.o) 

#Compile and link flags
CC	    = mpif90 
#CC      = /opt/openmpi/bin/mpif90
CFLAGS  = -O3 -fno-range-check -fbounds-check
DEBUG   =
LDFLAGS =  
EXE     = RELAX 

#Link and executable create
$(EXE) : $(OBJS) 
	$(CC) $(DEBUG) $(LDFLAGS) -o $@ $(OBJS) 

# Compile using a pattern rule 
%.o : %.f90 $(HDRS)   
	$(CC) -c $(CFLAGS) $(DEBUG) -o $@ $<

clean:
	rm -f  ../Routines/*.o ../Modules/*.o core *.mod

clobber:
	rm -f  ../Data/*.dat ../Modules/*.o ../Routines/*.o core *.mod $(EXE) 

#tar:
#	tar -cf data.tar ./data
#	tar -cf mie.tar $(SRCS) $(HDRS) $(DATA) Makefile 
