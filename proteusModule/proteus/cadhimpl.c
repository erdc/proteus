#include "adh_revision.h"
#include "share_declare.h"
/*#define ADH_SHARE_INCLUDED*/
#include "cadhimpl.h"
#include "assert.h"

void readBC(ADH_Input* self, char* filebase)
{
  /*adh*/
  self->iout = 0;			/* output counter */
  self->ixy = 0;			/* counter for the XY series */
  self->irn = 2;			/* run name counter */
  self->interval = 0;		/* interval var */
  self->hot_count = 0;		/* keeps track of when a hot start file needs to be printed */
  self->fp_in = NULL;		/* input file */
  self->fp_sflx_in = NULL;	/* input file */
  self->fp_sflx_out = NULL;	/* output surface heat flux */
#ifdef _ADH_XT3
  self->stdbuf = NULL;
#endif
  
  /* set values of myid and npes, and start up the message passing */
  debug_initialize();
  
  /*messg_initialize(argc, argv);*/
  /* mpi already initialized */
  int i_processor = 0;		/* Loop Counter over processors */
#ifdef _MPI
  int ierr_code = MPI_ERR_UNKNOWN;	/* the error code from an mpi call */
  
  /*ierr_code = MPI_Init(&argc, &argv);
    if(ierr_code != MPI_SUCCESS)
    messg_err(ierr_code);*/
  
  ierr_code = MPI_Comm_dup(MPI_COMM_WORLD, &ADH_COMM);
  if(ierr_code != MPI_SUCCESS)
    messg_err(ierr_code);
  
  msg_status = (MPI_Status *) NULL;
  msg_request = (MPI_Request *) NULL;
#endif
  
  npes = messg_comm_size();
  myid = messg_comm_rank();
  msg_starttime = messg_wtime();
  
  /* allocates the message arrays */
  send_msg = (MESSG_BUFFER *) tl_alloc(sizeof(MESSG_BUFFER), npes);
  recv_msg = (MESSG_BUFFER *) tl_alloc(sizeof(MESSG_BUFFER), npes);
  send_edge_msg = (MESSG_BUFFER *) tl_alloc(sizeof(MESSG_BUFFER), npes);
  recv_edge_msg = (MESSG_BUFFER *) tl_alloc(sizeof(MESSG_BUFFER), npes);
  send_key = (MESSG_KEY *) tl_alloc(sizeof(MESSG_KEY), npes);
  recv_init = (int *)tl_alloc(sizeof(int), npes);
  nsend = (int *)tl_alloc(sizeof(int), npes);
  nrecv = (int *)tl_alloc(sizeof(int), npes);
  nsend_edge = (int *)tl_alloc(sizeof(int), npes);
  nrecv_edge = (int *)tl_alloc(sizeof(int), npes);
  send_edge_key = (EDGE_LIST_ITEM ***) tl_alloc(sizeof(EDGE_LIST_ITEM **), npes);
  recv_edge_key = (EDGE_LIST_ITEM ***) tl_alloc(sizeof(EDGE_LIST_ITEM **), npes);
  for(i_processor = 0; i_processor < npes; i_processor++)
    {
      messg_init_buff(send_msg + i_processor, i_processor);
      messg_init_buff(recv_msg + i_processor, i_processor);
      messg_init_buff(send_edge_msg + i_processor, i_processor);
      messg_init_buff(recv_edge_msg + i_processor, i_processor);
      send_key[i_processor].size = 0;
      send_key[i_processor].key = NULL;
      recv_init[i_processor] = UNSET_INT;
      nsend[i_processor] = 0;
      nrecv[i_processor] = 0;
    }

  nmsg_counter = 0;
  nmsg_status = 0;
  nmsg_request = 0;
  max_nmsg_status = 0;
  max_nmsg_request = 0;
  /* */
  solv_initialize();
  
  self->fp_in = tl_efopen(filebase, "", UNSET_INT, ".adh", UNSET_INT, "rb");
  
#ifdef _ADH_XT3
  /* create buffer for standard output to improve performance on XT3 */
  /* self->stdbuf=malloc(51200); */
  /* setvbuf(stdout,self->stdbuf,_IOFBF,51200); */
  self->stdbuf = malloc(10240);
  setvbuf(stdout, self->stdbuf, _IOFBF, 10240);
#endif
  
  /* reads the input file */
  printf("Read binary input file\n");
  init_binin(self->fp_in);
  
  /* initialize boat info if required */
  if(BOAT)
    boat_init(filebase);
  
  /* sets up the linked lists */
  tl_list_setup();
  /* after reading in all the elements, find the connections from the */
  /* 1d elements to the 2d elements ( the 3d to 2d comes about in the boundary read */
  /* this can't be before the tl_list_setup  */
  elem1d_find_elem2d_init();

  elem1d_outward_nrml();
  
  fp_surf_temp = NULL;
  printf("open file pointers\n");
  if(HEAT && RAY_TRACING && myid <= 0)
    {
      self->fp_sflx_out =
	tl_efopen(proj_name, "_surface_fluxes", UNSET_INT, ".dat", UNSET_INT, "w");
      if(!SOCKETS)
	{
	  self->fp_sflx_in =
	    tl_efopen(proj_name, "_FacetEnergy", UNSET_INT, ".dat", UNSET_INT, "r");
	}
    }
}
/**A modified version of messg_finalize to allow PyADH to handle MPI*/
void messg_finalize_proteus(void)
{
  double elapsed_time = DZERO;	/* runtimes for the model */
  int i = 0;			/* Loop Counter */
#ifdef _MPI
  int ierr_code = MPI_ERR_UNKNOWN;	/* error code from an mpi call */
#endif

  messg_barrier();
  printf("**********************************************************************\n");
  elapsed_time = messg_wtime();
  elapsed_time -= msg_starttime;
  elapsed_time = messg_dmax(elapsed_time);
  if(myid == 0)
    printf("The elapsed time in seconds is: %14.6e\n", elapsed_time);

  if(recv_edge_key != NULL)
    recv_edge_key =
      (EDGE_LIST_ITEM ***) tl_free(sizeof(EDGE_LIST_ITEM **), npes, recv_edge_key);
  if(send_edge_key != NULL)
    send_edge_key =
      (EDGE_LIST_ITEM ***) tl_free(sizeof(EDGE_LIST_ITEM **), npes, send_edge_key);
  if(nrecv_edge != NULL)
    nrecv_edge = (int *)tl_free(sizeof(int), npes, nrecv_edge);
  if(nsend_edge != NULL)
    nsend_edge = (int *)tl_free(sizeof(int), npes, nsend_edge);
  if(nrecv != NULL)
    nrecv = (int *)tl_free(sizeof(int), npes, nrecv);
  if(nsend != NULL)
    nsend = (int *)tl_free(sizeof(int), npes, nsend);
  if(recv_init != NULL)
    recv_init = (int *)tl_free(sizeof(int), npes, recv_init);
  if(send_key != NULL)
    {
      for(i = 0; i < npes; i++)
	{
	  if(send_key[i].size != 0)
	    send_key[i].key =
	      (int *)tl_free(sizeof(int), send_key[i].size, send_key[i].key);
	}
      send_key = (MESSG_KEY *) tl_free(sizeof(MESSG_KEY), npes, send_key);
    }
  if(recv_edge_msg != NULL)
    {
      for(i = 0; i < npes; i++)
	{
	  messg_free(recv_edge_msg + i);
	}
      recv_edge_msg = (MESSG_BUFFER *) tl_free(sizeof(MESSG_BUFFER), npes, recv_edge_msg);
    }
  if(send_edge_msg != NULL)
    {
      for(i = 0; i < npes; i++)
	{
	  messg_free(send_edge_msg + i);
	}
      send_edge_msg = (MESSG_BUFFER *) tl_free(sizeof(MESSG_BUFFER), npes, send_edge_msg);
    }
  if(recv_msg != NULL)
    {
      for(i = 0; i < npes; i++)
	{
	  messg_free(recv_msg + i);
	}
      recv_msg = (MESSG_BUFFER *) tl_free(sizeof(MESSG_BUFFER), npes, recv_msg);
    }
  if(send_msg != NULL)
    {
      for(i = 0; i < npes; i++)
	{
	  messg_free(send_msg + i);
	}
      send_msg = (MESSG_BUFFER *) tl_free(sizeof(MESSG_BUFFER), npes, send_msg);
    }

#ifdef _MPI
  if((max_nmsg_request > 0) || (msg_request != NULL))
    {
      msg_request =
	(MPI_Request *) tl_free(sizeof(MPI_Request), max_nmsg_request, msg_request);
    }
  if((max_nmsg_status > 0) || (msg_status != NULL))
    {
      msg_status = (MPI_Status *) tl_free(sizeof(MPI_Status), max_nmsg_status, msg_status);
    }

  ierr_code = MPI_Comm_free(&ADH_COMM);
  if(ierr_code != MPI_SUCCESS)
    messg_err(ierr_code);

  /*  ierr_code = MPI_Finalize();
  if(ierr_code != MPI_SUCCESS)
  messg_err(ierr_code);*/
#endif
  nmsg_counter = 0;
  nmsg_status = 0;
  nmsg_request = 0;
  max_nmsg_status = 0;
  max_nmsg_request = 0;
}

int ADH_NumericalSolution_step(ADH_NumericalSolution* self)
{
  /* evaluate met data at beginning of time step */
  if(HEAT && METEOROLOGY)
    tc_eval_met_data(t_prev);
  /* loops over the XY series and sets the values */
  for(self->ixy = 0; self->ixy < nseries; self->ixy++)
    if(xy_sers[self->ixy].type == TIME_SERIES)
      tc_set_value(self->ixy, t_prev, t_prev + dt);
    else if(xy_sers[self->ixy].type == WIND_SERIES)
      {
	self->interval = tc_find_interval(self->ixy, t_prev);
	xy_sers[self->ixy].x_stress = tc_eval_series(self->ixy, self->interval, t_prev);
	xy_sers[self->ixy].y_stress = tc_eval_series_2(self->ixy, self->interval, t_prev);
      }
  if(HOT_START)
    {
      if(self->hot_count == hot_step)
	{
	  ps_print_hot(t_prev);
	  self->hot_count = 0;
	}
      else
	self->hot_count++;
    }

  return 0;
}

int ADH_NumericalSolution_stepTaken(ADH_NumericalSolution* self)
{

  if(STEADY_STATE_HYDRO && !NEW_STEADY_HYDRO)
    {
      NEW_STEADY_HYDRO = tc_new_interval(sth_sers, t_prev, t_prev + dt);
    }

  if(ADAPTION)
    err_main();

  /* prints the flux line data */
  if(SW2_FLUX && SW2_TRANSPORT)
    print_sw2_concentration_flux(t_prev + dt);
  if(SW2_FLUX)
    print_sw2_flux(t_prev + dt);

  /* prints the solution if necessary */

  if(fabs(t_prev + dt - tfinal) < SMALL8 ||
     (STEADY_STATE && old_spatial_residual < tol_nonlin))
    ioutput_flag = YES;

  if(ioutput_flag == YES)
    {
      if(t_prev + dt == tfinal)
	{
	  ps_print(tfinal);
	  if(HEAT)
	    {
	      if(myid <= 0)
		printf("Write final surface temperatures.\n");
	      print_surf_ts(fp_surf_temp, temper, tfinal);
	    }
	}
      else if(old_spatial_residual < tol_nonlin && STEADY_STATE)
	ps_print(t_prev + dt);
      else
	while((xy_sers[output_sers].size > self->iout) &&
	      (xy_sers[output_sers].sers[self->iout].time <= t_prev + dt))
	  {
	    ps_print(xy_sers[output_sers].sers[self->iout].time);
	    self->iout++;
	    /* write the adapted meshes */

		/*mwf stacy says to blow this away print_2dmesh not an adaptive routine set PRN_ADPT flag instead*/
/* 	    if(ADAPTION) */
/* 	      { */
/* 		if(nelem3d > 0) */
/* 		  print_3dmesh(); */
/* 		if(nelem2d > 0) */
/* 		  { */
/* 		    if(nelem3d == 0) */
/* 		      print_2dmesh(); */
/* 		    if(nelem3d > 0 && RAY_TRACING) */
/* 		      print_2dmesh(); */
/* 		  } */
/* 	      } */

	    /* read surface heat info */
	    if(HEAT)
	      {
		/* write the computed surface temperatures */

		if(myid <= 0)
		  printf("Writing temperatures time %f\n", t_prev + dt);
		print_surf_ts(fp_surf_temp, temper, t_prev + dt);

		if(RAY_TRACING)
		  {
		    if(myid <= 0)
		      printf("Reading heat fluxes at time %f\n", t_prev + dt);
		    read_surface_heat_flux(self->fp_sflx_in, self->fp_sflx_out);
		  }
	      }
	  }
    }

  /**err_main();  JNT moved up to remove -3 in error file **/
  /* adapts the mesh */
  if(ADAPTION)
    {
      adpt_main();

      /* update elemental quantities for GW after adaption */
      if(GW_FLOW)
	fe_gw_elem_sat_vel();
    }

  if(SW3_FLOW && ADPT_FLAG == TRUE)
    get_column();

  return tc_end() == NO;
}

int ADH_NumericalSolution_init(ADH_NumericalSolution* self, int argc, char* argv[])
{

  /*adh*/
  int ii,ie;			/* counters */
  self->iout = 0;			/* output counter */
  self->ixy = 0;			/* counter for the XY series */
  self->irn = 2;			/* run name counter */
  self->interval = 0;		/* interval var */
  self->hot_count = 0;		/* keeps track of when a hot start file needs to be printed */
  self->fp_in = NULL;		/* input file */
  self->fp_sflx_in = NULL;	/* input file */
  self->fp_sflx_out = NULL;	/* output surface heat flux */
#ifdef _ADH_XT3
  self->stdbuf = NULL;
#endif

  /* set values of myid and npes, and start up the message passing */
  debug_initialize();

  /*messg_initialize(argc, argv);*/
  /* mpi already initialized */
  int i_processor = 0;		/* Loop Counter over processors */
#ifdef _MPI
  int ierr_code = MPI_ERR_UNKNOWN;	/* the error code from an mpi call */

  /*ierr_code = MPI_Init(&argc, &argv);
  if(ierr_code != MPI_SUCCESS)
  messg_err(ierr_code);*/

  ierr_code = MPI_Comm_dup(MPI_COMM_WORLD, &ADH_COMM);
  if(ierr_code != MPI_SUCCESS)
    messg_err(ierr_code);

  msg_status = (MPI_Status *) NULL;
  msg_request = (MPI_Request *) NULL;
#endif

  npes = messg_comm_size();
  myid = messg_comm_rank();
  msg_starttime = messg_wtime();

  /* allocates the message arrays */
  send_msg = (MESSG_BUFFER *) tl_alloc(sizeof(MESSG_BUFFER), npes);
  recv_msg = (MESSG_BUFFER *) tl_alloc(sizeof(MESSG_BUFFER), npes);
  send_edge_msg = (MESSG_BUFFER *) tl_alloc(sizeof(MESSG_BUFFER), npes);
  recv_edge_msg = (MESSG_BUFFER *) tl_alloc(sizeof(MESSG_BUFFER), npes);
  send_key = (MESSG_KEY *) tl_alloc(sizeof(MESSG_KEY), npes);
  recv_init = (int *)tl_alloc(sizeof(int), npes);
  nsend = (int *)tl_alloc(sizeof(int), npes);
  nrecv = (int *)tl_alloc(sizeof(int), npes);
  nsend_edge = (int *)tl_alloc(sizeof(int), npes);
  nrecv_edge = (int *)tl_alloc(sizeof(int), npes);
  send_edge_key = (EDGE_LIST_ITEM ***) tl_alloc(sizeof(EDGE_LIST_ITEM **), npes);
  recv_edge_key = (EDGE_LIST_ITEM ***) tl_alloc(sizeof(EDGE_LIST_ITEM **), npes);
  for(i_processor = 0; i_processor < npes; i_processor++)
    {
      messg_init_buff(send_msg + i_processor, i_processor);
      messg_init_buff(recv_msg + i_processor, i_processor);
      messg_init_buff(send_edge_msg + i_processor, i_processor);
      messg_init_buff(recv_edge_msg + i_processor, i_processor);
      send_key[i_processor].size = 0;
      send_key[i_processor].key = NULL;
      recv_init[i_processor] = UNSET_INT;
      nsend[i_processor] = 0;
      nrecv[i_processor] = 0;
    }

  nmsg_counter = 0;
  nmsg_status = 0;
  nmsg_request = 0;
  max_nmsg_status = 0;
  max_nmsg_request = 0;
  /* */
  solv_initialize();

  /* checks usage */
  if(argc >= 2)
    {
      /* version check */
      if(strncmp(argv[1], "-v", 2) == AGREE)
	{
	  printf("ADH svn revision #%s\n", ADHREV);
	  /*messg_finalize();*/
	  messg_finalize_proteus();
	  exit(0);
	}
      else
	{			/* opens the input file */
	  self->fp_in = tl_efopen(argv[1], "", UNSET_INT, ".adh", UNSET_INT, "rb");
	  strcpy(proj_name, argv[1]);	/* save the project name for later */
	  strcpy(run_name, "");
	  for(self->irn = 2; self->irn < argc; self->irn++)
	    {
	      strcat(proj_name, "_");
	      strcat(run_name, " ");
	      strcat(run_name, argv[self->irn]);
	      strcat(proj_name, argv[self->irn]);
	    }
	}
    }
  else
    {
      /* bad call */
      fprintf(stderr, "Incorrect usage: \n    adh file_base optional_run_name\n");
      fprintf(stderr, "  or,\n    adh -v\n");
      /*messg_finalize();*/
      messg_finalize_proteus();
      exit(0);
    }

#ifdef _ADH_XT3
  /* create buffer for standard output to improve performance on XT3 */
  /* self->stdbuf=malloc(51200); */
  /* setvbuf(stdout,self->stdbuf,_IOFBF,51200); */
  self->stdbuf = malloc(10240);
  setvbuf(stdout, self->stdbuf, _IOFBF, 10240);
#endif

  /* reads the input file */
  printf("Read binary input file\n");
  init_binin(self->fp_in);

  /* initialize boat info if required */
  if(BOAT)
    boat_init(argv[1]);

  /* sets up the linked lists */
  tl_list_setup();
  /* after reading in all the elements, find the connections from the */
  /* 1d elements to the 2d elements ( the 3d to 2d comes about in the boundary read */
  /* this can't be before the tl_list_setup  */
  elem1d_find_elem2d_init();

  elem1d_outward_nrml();

  fp_surf_temp = NULL;
  printf("open file pointers\n");
  if(HEAT && RAY_TRACING && myid <= 0)
    {
      self->fp_sflx_out =
	tl_efopen(proj_name, "_surface_fluxes", UNSET_INT, ".dat", UNSET_INT, "w");
      if(!SOCKETS)
	{
	  self->fp_sflx_in =
	    tl_efopen(proj_name, "_FacetEnergy", UNSET_INT, ".dat", UNSET_INT, "r");
	}
    }

  if(myid <= 0)
    {
      /* open the file pointers for the gms-formatted files */
      print_xms_pointers();
    }

  /* partitions the problem */
  if(myid <= 0)
    printf("Initial partition of the mesh.\n");
  partition_init();

  /* cleans up the mesh */
  partition_cleanup();

  /* renumber the nodes and elements */
  if(myid <= 0)
    printf("Renumber nodes and elements.\n");
  node_renumber();
  elem3d_renumber();
  elem2d_renumber();
  elem1d_renumber();
  /* sets up the communications */
  comm_set_keys();

  /* update the ghost information */
  comm_update_GLOBAL_NODE(0, node_pair);

  /* initialize elemental quantities for GW */
  if(GW_FLOW)
    {
      fe_gw_elem_sat_vel();

      /* have to initialize the old saturations here */
      for(ii = 0; ii < NDPRELM; ii++)
	{
	  for(ie = 0; ie < nelem3d; ie++)
	    {
	      old_saturation[ie][ii] = saturation[ie][ii];
	    }
	}
    }
  /* prints the initial conditions */
  ps_print(tinitial);

  /* calculate the initial value for each XY series and time step size */
  tc_init();

  /* prints the solution initially if necessary */
  if(ioutput_flag == YES)
    while((xy_sers[output_sers].size > self->iout) &&
	  (xy_sers[output_sers].sers[self->iout].time < t_prev))
      self->iout++;

  /* starts the timer */
  /* messg_start_time(); */

#ifdef _METIS
  partition_main();
  global_update();
#endif

  /* reassign global node and element numbers */
  renumber_globals();

  /* print the initial 2d and 3d meshes (even if adaption is off) */
  if(nelem3d > 0)
    print_3dmesh();
  if(nelem2d > 0)
    {
      /*mwf stacy says to skip this because print_2dmesh not an adaptive routine set PRN_ADPT flag instead*/
      if(nelem3d == 0 && !ADAPTION)
	print_2dmesh();
      if(nelem3d > 0 && RAY_TRACING)
	print_2dmesh();
    }
  if(SW3_FLOW)
    get_column();

  if(HEAT)
    {
      if(myid <= 0)
	/* printf("Writing initial temperatures and reading initial surface fluxes.\n"); */
	printf("Writing initial temperatures and reading initial surface fluxes.\n");
      print_surf_ts(fp_surf_temp, temper, t_prev);

      if(RAY_TRACING)
	read_surface_heat_flux(self->fp_sflx_in, self->fp_sflx_out);
    }

  return 0;
}

int ADH_NumericalSolution_dealloc(ADH_NumericalSolution* self)
{
  /* Close Out All Pieces/Parts */
  tl_free_global();
  /* closes out the communications */
  messg_barrier();
  /* messg_elapsed_time(); */
  /*messg_finalize();*/
  messg_finalize_proteus();
  solv_finalize();
  debug_finalize();

  /*cek try resetting all these as in their declarations*/
vwork_size = 0;
flx_dt = 0.0;		/* dt used for flux calculations */
prof_rows = NULL;	/* the rows of the profile matrix */
prof_cols = NULL;	/* the columns of the profile matrix */
coarse_inv = NULL;	/* the inverse of the coarse matrix */
vr = NULL;		/* the restricted vector */
my_vr = NULL;		/* the restricted vector for my processor */
coarse_result = NULL;	/* the result of the coarse matrix vector product */
latitude = 32.97;	/* latitude, default is Yuma, AZ */
longitude = -114.26;	/* longitude, default is Yuma, AZ */
std_stdt = 1.;		/* initial dt for steady-state calculations */
sth_sers = UNSET_INT;	/* series controlling the quasi-steady-state updates */
sth_maxi = 20;		/* maximum number of steady-state iterations for each quasi-steady-state update */
sth_stdt = 1.;		/* initial dt for quasi-steady-state calculations */
temp_dt = 0.0;
original_dt = 1.0;
spatial_residual = 0.0;
old_spatial_residual = 0.0;
c_mu = 0.09;
c_mu_root = 0.3;
c_epsilon_1 = 1.44;
c_epsilon_2 = 1.92;
von_karman = 0.418;
julian_day_start = UNSET_INT;	/* the julian day at the start of the simulation */
julian_day_end = UNSET_INT;	/* the julian day at the end of the simulation */
num_prints = 0;		/* number of snapshots printed */
fetch_split = 72;		/* the number of angles to use in nodal fetch search */
prof_size = 0;		/* the size of the profile matrix */
prs_adj_node = UNSET_INT;	/* the pressure adjustment node number */
prs_adj_series = UNSET_INT;	/* the pressure adjustment series number */
iprint_step = 0;		/* the current print step */
coarse_size = 0;		/* the size of the coarse problem */
coarse_vect_size = 0;	/* the size of a coarse vector */
my_coarse_vect_size = 0;	/* the size of my part of a coarse vector */
iunref_flag = NO;		/* flag to determine whether unrefinement took place */
NS_NFS = FALSE;		/* if yes then velocity grad printedout, NFS flag */
SW2_NFS = FALSE;		/*if true then velocity grad printout in SW2 GSavant 10/10/2006 */
ADAPTION = FALSE;		/* if TRUE then permit spatial adaption */
PRN_ADPT = FALSE;		/* if TRUE then write adapted mesh and data files - currently only SW2 */
STEADY_STATE = FALSE;	/* if TRUE then abandon time accuracy and push to long-time solution */
GW_FLOW = FALSE;		/* if TRUE then run groundwater */
GW_TRANSPORT = FALSE;	/* if TRUE then run transport with groundwater */
HEAT = FALSE;		/* if TRUE then run heat transport (w/GW only for now) */
OL_FLOW = FALSE;		/* if TRUE then run diffusive wave overland flow */
OL_TRANSPORT = FALSE;	/* if TRUE then run transport with diffusive wave overland flow */
CH_FLOW = FALSE;		/* if TRUE then run channel flow */
NS_FLOW = FALSE;		/* if TRUE then run navier stokes */
NS_TRANSPORT = FALSE;	/* if TRUE then run transport with navier stokes */
RAY_TRACING = FALSE;	/* if TRUE then couple to ray tracing code for heat bcs */
SW2_FLOW = FALSE;		/* if TRUE then run 2D shallow water */
SW2_TRANSPORT = FALSE;	/* if TRUE then run transport with 2D shallow water */
SW2_FLUX = FALSE;		/* if TRUE then print SW2 flux output file */
SW3_FLOW = FALSE;		/* if TRUE then run 3D shallow water */
SW3_TRANSPORT = FALSE;	/* if TRUE then run transport with 2D shallow water */
USER_PHYSICS = FALSE;	/* if TRUE then run user-defined physics */
MOVING_GRID = FALSE;	/* if TRUE then the grid is allowed to move */
refresh_pc = YES;		/* refresh_pc set to YES means that the preconditioner is updated */
WAVE = FALSE;		/* the flag indicating short wave effects are included */
BOAT = FALSE;		/* the flag indicating boats are included (SW2) */
BOAT_ENT = FALSE;		/* the flag indicating vessel entrainment is included (SW2) */
CONVEYANCE = FALSE;		/* the flag indicating total conveyance is needed (SW2) */
METEOROLOGY = FALSE;	/* the flag indicating met data will be read */
SOCKETS = FALSE;		/* the flag indicating sockets rather than files will be used */
lin_fail = NO;		/* the flag representing a linear failure */
STEADY_STATE_HYDRO = FALSE;	/* the flag indicating hydro will be done in quasi-steady-state for sw2 */
NEW_STEADY_HYDRO = TRUE;	/* the flag indicating a new hydro solution is needed for quasi-steady-state */
iac_flag = NO;		/* flag for timestep interpolation accuracy checking */
time_zone = 7;		/* time zone, default is Yuma, AZ */
windwavemod = OFF;		/* flag to determine which model to use to calculate combined wave-current shear stress */
SEDIMENT = FALSE;		/* sediment calculations done or not */
SED_HOT = FALSE;		/* sediment hotstarted or not */
HOT_START = FALSE;		/* all output data for last time step converted to hotstart */
ICE = FALSE;		/* indicates that ice is included (SW2) */
INS = FALSE;		/* indicates that ice is designated by a circular string */
baroclinic_code = 0;	/* if the baroclinic code = 0  no density effects */
ADPT_FLAG = 0;
wcdt_cnt = 0, wcdt_size = 0;
units_flag = UNSET_INT;
timestart = 0;

  return 0;
}

int ADH_NumericalSolution_calculateSolution(ADH_NumericalSolution* self)
{
  /* runs the time steps */
  if(myid <= 0)
    printf("Begin the master time step loop.\n");
  do
    {
      /* repeats the time step until it should keep it */
      /* evaluate met data at beginning of time step */
      if(HEAT && METEOROLOGY)
  	tc_eval_met_data(t_prev);
      do
  	{
  	  /* loops over the XY series and sets the values */
  	  for(self->ixy = 0; self->ixy < nseries; self->ixy++)
  	    if(xy_sers[self->ixy].type == TIME_SERIES)
  	      tc_set_value(self->ixy, t_prev, t_prev + dt);
  	    else if(xy_sers[self->ixy].type == WIND_SERIES)
  	      {
  		self->interval = tc_find_interval(self->ixy, t_prev);
  		xy_sers[self->ixy].x_stress = tc_eval_series(self->ixy, self->interval, t_prev);
  		xy_sers[self->ixy].y_stress = tc_eval_series_2(self->ixy, self->interval, t_prev);
  	      }
  	  if(HOT_START)
  	    {
  	      if(self->hot_count == hot_step)
  		{
  		  ps_print_hot(t_prev);
  		  self->hot_count = 0;
  		}
  	      else
  		self->hot_count++;
  	    }
  	}
      while(fe_main() == NO);

      if(STEADY_STATE_HYDRO && !NEW_STEADY_HYDRO)
  	{
  	  NEW_STEADY_HYDRO = tc_new_interval(sth_sers, t_prev, t_prev + dt);
  	}

      if(ADAPTION)
  	err_main();

      /* prints the flux line data */
      if(SW2_FLUX && SW2_TRANSPORT)
  	print_sw2_concentration_flux(t_prev + dt);
      if(SW2_FLUX)
  	print_sw2_flux(t_prev + dt);

      /* prints the solution if necessary */

      if(fabs(t_prev + dt - tfinal) < SMALL8 ||
  	 (STEADY_STATE && old_spatial_residual < tol_nonlin))
  	ioutput_flag = YES;

      if(ioutput_flag == YES)
  	{
  	  if(t_prev + dt == tfinal)
  	    {
  	      ps_print(tfinal);
  	      if(HEAT)
  		{
  		  if(myid <= 0)
  		    printf("Write final surface temperatures.\n");
  		  print_surf_ts(fp_surf_temp, temper, tfinal);
  		}
  	    }
  	  else if(old_spatial_residual < tol_nonlin && STEADY_STATE)
  	    ps_print(t_prev + dt);
  	  else
  	    while((xy_sers[output_sers].size > self->iout) &&
  		  (xy_sers[output_sers].sers[self->iout].time <= t_prev + dt))
  	      {
  		ps_print(xy_sers[output_sers].sers[self->iout].time);
  		self->iout++;
  		/* write the adapted meshes */
		
		/*mwf stacy says to skip this if adaption because print_2dmesh not an adaptive routine set PRN_ADPT flag instead*/
/*   		if(ADAPTION) */
/*   		  { */
/*   		    if(nelem3d > 0) */
/*   		      print_3dmesh(); */
/*   		    if(nelem2d > 0) */
/*   		      { */
/*   			if(nelem3d == 0) */
/*   			  print_2dmesh(); */
/*   			if(nelem3d > 0 && RAY_TRACING) */
/*   			  print_2dmesh(); */
/*   		      } */
/*   		  } */

  		/* read surface heat info */
  		if(HEAT)
  		  {
  		    /* write the computed surface temperatures */

  		    if(myid <= 0)
  		      printf("Writing temperatures time %f\n", t_prev + dt);
  		    print_surf_ts(fp_surf_temp, temper, t_prev + dt);

  		    if(RAY_TRACING)
  		      {
  			if(myid <= 0)
  			  printf("Reading heat fluxes at time %f\n", t_prev + dt);
  			read_surface_heat_flux(self->fp_sflx_in, self->fp_sflx_out);
  		      }
  		  }
  	      }
  	}

      /**err_main();  JNT moved up to remove -3 in error file **/
      /* adapts the mesh */
      if(ADAPTION)
  	{
  	  adpt_main();

  	  /* update elemental quantities for GW after adaption */
  	  if(GW_FLOW)
  	    fe_gw_elem_sat_vel();
  	}

      if(SW3_FLOW && ADPT_FLAG == TRUE)
  	get_column();
    }
  while(tc_end() == NO);

  if(HEAT && RAY_TRACING && !SOCKETS && myid <= 0)
    fclose(self->fp_sflx_in);

  if(HOT_START)
    ps_print_hot(tfinal);

  if(myid <= 0)
    {
      printf("Printing xms trailers\n");
      print_xms_trailers();
    }
  return 0;
}
/** A modifed version of the reallocation portion of fe_main */
int realloc_fe_main(ADH_OneLevelTransport* fe)
{
  int i,j;
  /* int i, j, k, ie, n;		/\* loop counter *\/ */
  /* int itrn;			/\* loop counter over the transported quantities *\/ */
  /* int nsys;			/\* the number of systems being solved *\/ */
  /* int nsys_sq;			/\* the number of systems being solved squared *\/ */
  int isize_prev = 0;		/* the previous isize */
  /* int sed_iteration;		/\* the number of iterations done to balance the fractions in active layer *\/ */
  /* int sed_number;		/\* given the number of clay or silt what is its number in sediment *\/ */
  /* double ald_check;		/\* maximum change in the active layer distribution since the last iteration *\/ */
  /* double ald_max = DONE; */
  /* double t_sed;			/\* the real time for the sediment current sediment step *\/ */
  /* double flowrate_x;		/\* the flowrate for the hydro results to be kept for the sed timesteps *\/ */
  /* double flowrate_y;		/\* the flowrate for the hydro results to be kept for the sed timesteps *\/ */
  /* double ald_diff; */
  /* int max_node; */
  /* double sed_it_tol; */
  /* double bed_thick_old; */
  /* double bed_thick_new; */
  /* double esrh;			/\* the equivalent sand roughness height *\/ */
  /* double blsg; */
  /* double albd; */
  /* double asbd; */
  /* VECT2D tvel;			/\* the local transverse velocity components *\/ */
  /* int count = 0;		/\* counter for sed time step *\/ */
  /* int istr;			/\* counter for ceq boundary condition *\/ */
/* #ifdef _MESSG */
/*   double sed_time_start; */
/*   double sed_time_end; */
/*   double sed_time_total; */
/* #endif */
  /* double qus_depth_factor = 0.1; */
  /* double qus_max_change; */
  /* double sdtemp; */
  /* double temp_dpl; */
  /* double ust; */
  /* double blv; */
  /* double blvelmag; */
  /* allocates memory if necessary */
  if(nnode > fe->isize)
    {

      isize_prev = fe->isize;
      fe->isize = nnode;
      fe->bc_mask =
	(int *)tl_realloc(sizeof(int), fe->isize * max_nsys, isize_prev * max_nsys, fe->bc_mask);

      if(NS_FLOW && transport.tke >= 0)
	{
	  fe->nodal_production =
	    (double *)tl_realloc(sizeof(double), fe->isize, isize_prev, fe->nodal_production);
	  fe->nodal_viscosity =
	    (double *)tl_realloc(sizeof(double), fe->isize, isize_prev, fe->nodal_viscosity);
	  fe->nodal_volume =
	    (double *)tl_realloc(sizeof(double), fe->isize, isize_prev, fe->nodal_volume);
	}

      if(SW2_FLOW || NS_FLOW)
	{
	  if(SEDIMENT)
	    {
	      if(isize_prev == 0)
		{
		  fe->old_sed_ald = (double *)tl_alloc(sizeof(double), nsed);
		  fe->fraction = (double *)tl_alloc(sizeof(double), nsed);
		  fe->old_fraction = (double *)tl_alloc(sizeof(double), nsed);
		  fe->ald = (double *)tl_alloc(sizeof(double), nsed);
		  fe->asd = (double *)tl_alloc(sizeof(double), nsed);
		  fe->blt = (double *)tl_alloc(sizeof(double), number_bed_layers);
		  fe->bl_bd = (double *)tl_alloc(sizeof(double), number_bed_layers);
		  fe->bl_ces = (double *)tl_alloc(sizeof(double), number_bed_layers);
		  fe->bl_erc = (double *)tl_alloc(sizeof(double), number_bed_layers);
		  fe->bl_ere = (double *)tl_alloc(sizeof(double), number_bed_layers);
		  fe->bld = (double **)tl_alloc(sizeof(double *), number_bed_layers);
		  fe->bl_por = (double *)tl_alloc(sizeof(double), number_bed_layers);
		  fe->bl_sg = (double *)tl_alloc(sizeof(double), number_bed_layers);
		  for(j = 0; j < number_bed_layers; j++)
		    {
		      fe->bld[j] = (double *)tl_alloc(sizeof(double), nsed);
		    }
		  if(ncti > 1)
		    {
		      fe->cbp = (double **)tl_alloc(sizeof(double *), ncti);
		      for(j = 0; j < ncti; j++)
			{
			  fe->cbp[j] = (double *)tl_alloc(sizeof(double), 5);
			}
		    }
		  fe->old_sed_active_layer_distribution =
		    (double **)tl_alloc(sizeof(double *), nsed);
		  fe->old_sed_active_stratum_distribution =
		    (double **)tl_alloc(sizeof(double *), nsed);
		  fe->old_it_ald = (double **)tl_alloc(sizeof(double *), nsed);
		  fe->old_it_rhs = (double **)tl_alloc(sizeof(double *), nsed + 1);
		  fe->old_it_concentration = (double **)tl_alloc(sizeof(double *), nsed);
		}
	      fe->bed_load_3d =
		(VECT **) tl_realloc(sizeof(VECT *), fe->isize, isize_prev, fe->bed_load_3d);
	      fe->source_coef =
		(double **)tl_realloc(sizeof(double *), fe->isize, isize_prev, fe->source_coef);
	      fe->bl_eq_coef =
		(double **)tl_realloc(sizeof(double *), fe->isize, isize_prev, fe->bl_eq_coef);
	      fe->prob_dep =
		(double **)tl_realloc(sizeof(double *), fe->isize, isize_prev, fe->prob_dep);
	      fe->ustar = (double *)tl_realloc(sizeof(double), fe->isize, isize_prev, fe->ustar);
	      fe->eq_bed_load =
		(VECT2D **) tl_realloc(sizeof(VECT2D *), fe->isize, isize_prev, fe->eq_bed_load);
	      fe->bedload_eq_con =
		(double **)tl_realloc(sizeof(double *), fe->isize, isize_prev, fe->bedload_eq_con);
	      fe->bl_thick =
		(double **)tl_realloc(sizeof(double *), fe->isize, isize_prev, fe->bl_thick);
	      
	      fe->bl_slc = (double **)tl_realloc(sizeof(double *), fe->isize, isize_prev, fe->bl_slc);
	      fe->bl_vel = (VECT2D **) tl_realloc(sizeof(VECT2D *), fe->isize, isize_prev, fe->bl_vel);
	      fe->bed_character =
		(int *)tl_realloc(sizeof(int), fe->isize, isize_prev, fe->bed_character);
	      fe->aldkf = (int *)tl_realloc(sizeof(int), fe->isize, isize_prev, fe->aldkf);
	      fe->asdkf = (int *)tl_realloc(sizeof(int), fe->isize, isize_prev, fe->asdkf);
	      fe->sum_snd_ald =
		(double *)tl_realloc(sizeof(double), fe->isize, isize_prev, fe->sum_snd_ald);
	      fe->sum_nodal_source =
		(double *)tl_realloc(sizeof(double), fe->isize, isize_prev, fe->sum_nodal_source);
	      fe->sum_snd_con =
		(double *)tl_realloc(sizeof(double), fe->isize, isize_prev, fe->sum_snd_con);
	      fe->sum_blt = (double *)tl_realloc(sizeof(double), fe->isize, isize_prev, fe->sum_blt);
	      fe->bed_flag = (int *)tl_realloc(sizeof(int), fe->isize, isize_prev, fe->bed_flag);
	      fe->old_sed_active_stratum_ceiling =
		(double *)tl_realloc(sizeof(double), fe->isize, isize_prev,
				     fe->old_sed_active_stratum_ceiling);
	      fe->old_sed_displacement =
		(double *)tl_realloc(sizeof(double), fe->isize, isize_prev,
				     fe->old_sed_displacement);
	      fe->critical_shear_vel =
		(double **)tl_realloc(sizeof(double *), fe->isize, isize_prev,
				      fe->critical_shear_vel);
	      fe->nodal_grad =
		(VECT2D *) tl_realloc(sizeof(VECT2D), fe->isize, isize_prev, fe->nodal_grad);
	      fe->elev_disp =
		(double *)tl_realloc(sizeof(double), fe->isize, isize_prev, fe->elev_disp);
	      fe->AL_thickness =
		(double *)tl_realloc(sizeof(double), fe->isize, isize_prev, fe->AL_thickness);
	      fe->old_sed_active_layer_porosity =
		(double *)tl_realloc(sizeof(double), fe->isize, isize_prev,
				     fe->old_sed_active_layer_porosity);
	      fe->old_sed_active_layer_critical_erosion_shear =
		(double *)tl_realloc(sizeof(double), fe->isize, isize_prev,
				     fe->old_sed_active_layer_critical_erosion_shear);
	      fe->old_sed_active_layer_erosion_rate_constant =
		(double *)tl_realloc(sizeof(double), fe->isize, isize_prev,
				     fe->old_sed_active_layer_erosion_rate_constant);
	      fe->old_sed_active_layer_erosion_rate_exponent =
		(double *)tl_realloc(sizeof(double), fe->isize, isize_prev,
				     fe->old_sed_active_layer_erosion_rate_exponent);
	      fe->sum_sed_decay =
		(double *)tl_realloc(sizeof(double), fe->isize, isize_prev, fe->sum_sed_decay);
	      fe->limiter = (int *)tl_realloc(sizeof(int), fe->isize, isize_prev, fe->limiter);
	      fe->limiter2 = (int *)tl_realloc(sizeof(int), fe->isize, isize_prev, fe->limiter2);
	      fe->conc_converge =
		(int *)tl_realloc(sizeof(int), fe->isize, isize_prev, fe->conc_converge);

	      for(i = 0; i < nsed; i++)
		{
		  fe->old_sed_active_layer_distribution[i] =
		    (double *)tl_realloc(sizeof(double), fe->isize, isize_prev,
					 fe->old_sed_active_layer_distribution[i]);
		  fe->old_sed_active_stratum_distribution[i] =
		    (double *)tl_realloc(sizeof(double), fe->isize, isize_prev,
					 fe->old_sed_active_stratum_distribution[i]);
		  fe->old_it_ald[i] =
		    (double *)tl_realloc(sizeof(double), fe->isize, isize_prev, fe->old_it_ald[i]);
		  fe->old_it_concentration[i] =
		    (double *)tl_realloc(sizeof(double), fe->isize, isize_prev,
					 fe->old_it_concentration[i]);
		}
	      for(i = 0; i < nsed + 1; i++)
		{
		  fe->old_it_rhs[i] =
		    (double *)tl_realloc(sizeof(double), fe->isize, isize_prev, fe->old_it_rhs[i]);
		}

	      for(i = isize_prev; i < fe->isize; i++)
		{
		  fe->bed_load_3d[i] = (VECT *) tl_alloc(sizeof(VECT), nsed);
		  fe->eq_bed_load[i] = (VECT2D *) tl_alloc(sizeof(VECT2D), nsed);
		  fe->source_coef[i] = (double *)tl_alloc(sizeof(double), nsed);
		  fe->bl_eq_coef[i] = (double *)tl_alloc(sizeof(double), nsed);
		  fe->prob_dep[i] = (double *)tl_alloc(sizeof(double), nsed);
		  fe->critical_shear_vel[i] = (double *)tl_alloc(sizeof(double), nsed);
		  fe->bedload_eq_con[i] = (double *)tl_alloc(sizeof(double), nsed);
		  fe->bl_thick[i] = (double *)tl_alloc(sizeof(double), nsed);
		  fe->bl_slc[i] = (double *)tl_alloc(sizeof(double), nsed);
		  fe->bl_vel[i] = (VECT2D *) tl_alloc(sizeof(VECT2D), nsed);
		}
	    }
	}
      /* general allocations */
      fe->residual =
	(double *)tl_realloc(sizeof(double), fe->isize * max_nsys, isize_prev * max_nsys,
			     fe->residual);
      fe->sol =
	(double *)tl_realloc(sizeof(double), fe->isize * max_nsys, isize_prev * max_nsys, fe->sol);
      /* initialize sol */
      for(i = isize_prev * max_nsys; i < fe->isize * max_nsys; i++)
	{
	  fe->sol[i] = 0.0;
	}

      fe->scale_vect =
	(double *)tl_realloc(sizeof(double), fe->isize * max_nsys, isize_prev * max_nsys,
			     fe->scale_vect);
      fe->diagonal =
	(double *)tl_realloc(sizeof(double), fe->isize * max_nsys_sq, isize_prev * max_nsys_sq,
			     fe->diagonal);
      fe->matrix = (SPARSE_VECT *) tl_realloc(sizeof(SPARSE_VECT), fe->isize, isize_prev, fe->matrix);
      for(i = isize_prev; i < fe->isize; i++)
	{
	  fe->matrix[i].value = (double *)tl_alloc(sizeof(double), SPV_BLOCK * max_nsys_sq);
	  fe->matrix[i].index = (int *)tl_alloc(sizeof(int), SPV_BLOCK);
	  fe->matrix[i].max_size = SPV_BLOCK;
	  fe->matrix[i].size = 0;
	}
    }
  return 0;
}

extern int ADH_OneLevelTransport_init(ADH_OneLevelTransport* self)
{
  
  /* sets the default to a good time step */
  t_fail_flag = NO;

  /* updates the block numbers for the nodes */
#ifdef _MESSG
  comm_update_int(node_block, 1);
#endif

  /* sets the number of systems being solved */
  /* ordered by increasing number of dofs to get max size (seh) */
  max_nsys = 0;
  max_nsys_sq = 0;
  if(GW_FLOW || HEAT || OL_FLOW)
    {
      max_nsys = 1;
      max_nsys_sq = 1;
    }
  if(SW2_FLOW)
    {
      max_nsys = 3;
      max_nsys_sq = 9;
    }
  if(SW3_FLOW)
    {
      max_nsys = 3;
      max_nsys_sq = 9;
      /* these are made at their largest size which is actually the 2d sw equations which have 3 dof per node */
    }
  if(NS_FLOW)
    {
      max_nsys = 4;
      max_nsys_sq = 16;
    }
  /*see if can force refresh variable to make sure solv_blk_set is called to allocate preconditioner before first call to solve*/
  refresh_pc = YES; 
  realloc_fe_main(self);
  return 0;
}

int ADH_OneLevelTransport_dealloc(ADH_OneLevelTransport* self)
{
  int i;

  if(self->matrix != NULL)
    {
      for(i = 0; i < self->isize; i++)
	{
	  self->matrix[i].value =
	    (double *)tl_free(sizeof(double), self->matrix[i].max_size * max_nsys_sq,
			      self->matrix[i].value);
	  self->matrix[i].index =
	    (int *)tl_free(sizeof(int), self->matrix[i].max_size, self->matrix[i].index);
	}
      self->matrix = (SPARSE_VECT *) tl_free(sizeof(SPARSE_VECT), self->isize, self->matrix);
    }

  if(self->diagonal != NULL)
    self->diagonal = (double *)tl_free(sizeof(double), self->isize * max_nsys_sq, self->diagonal);

  if(self->bc_mask != NULL)
    self->bc_mask = (int *)tl_free(sizeof(int), self->isize * max_nsys, self->bc_mask);
  if(self->nodal_volume != NULL)
    self->nodal_volume = (double *)tl_free(sizeof(double), self->isize, self->nodal_volume);
  if(self->nodal_production != NULL)
    self->nodal_production = (double *)tl_free(sizeof(double), self->isize, self->nodal_production);
  if(self->nodal_viscosity != NULL)
    self->nodal_viscosity = (double *)tl_free(sizeof(double), self->isize, self->nodal_viscosity);

  if(self->residual != NULL)
    self->residual = (double *)tl_free(sizeof(double), self->isize * max_nsys, self->residual);
  if(self->sol != NULL)
    self->sol = (double *)tl_free(sizeof(double), self->isize * max_nsys, self->sol);
  if(self->scale_vect != NULL)
    self->scale_vect = (double *)tl_free(sizeof(double), self->isize * max_nsys, self->scale_vect);

  if(self->eq_bed_load != NULL)
    {
      for(i = 0; i < self->isize; i++)
	{
	  if(self->eq_bed_load[i] != NULL)
	    self->eq_bed_load[i] = (VECT2D *) tl_free(sizeof(VECT2D), nsand, self->eq_bed_load[i]);
	}
      self->eq_bed_load = (VECT2D * *)tl_free(sizeof(VECT2D *), self->isize, self->eq_bed_load);
    }

  if(self->bed_load_3d != NULL)
    {
      for(i = 0; i < self->isize; i++)
	{
	  if(self->bed_load_3d[i] != NULL)
	    self->bed_load_3d[i] = (VECT *) tl_free(sizeof(VECT), nsand, self->bed_load_3d[i]);
	}
      self->bed_load_3d = (VECT * *)tl_free(sizeof(VECT *), self->isize, self->bed_load_3d);
    }

  if(self->ustar != NULL)
    self->ustar = (double *)tl_free(sizeof(double), self->isize, self->ustar);

  if(self->ald != NULL)
    self->ald = (double *)tl_free(sizeof(double), nsed, self->ald);

  if(self->asd != NULL)
    self->asd = (double *)tl_free(sizeof(double), nsed, self->asd);

  if(self->blt != NULL)
    self->blt = (double *)tl_free(sizeof(double), number_bed_layers, self->blt);

  if(self->bl_por != NULL)
    self->bl_por = (double *)tl_free(sizeof(double), number_bed_layers, self->bl_por);

  if(self->bl_bd != NULL)
    self->bl_bd = (double *)tl_free(sizeof(double), number_bed_layers, self->bl_bd);

  if(self->bl_sg != NULL)
    self->bl_sg = (double *)tl_free(sizeof(double), number_bed_layers, self->bl_sg);

  if(self->bl_ces != NULL)
    self->bl_ces = (double *)tl_free(sizeof(double), number_bed_layers, self->bl_ces);

  if(self->bl_erc != NULL)
    self->bl_erc = (double *)tl_free(sizeof(double), number_bed_layers, self->bl_erc);

  if(self->bl_ere != NULL)
    self->bl_ere = (double *)tl_free(sizeof(double), number_bed_layers,self->bl_ere);

  if(self->bld != NULL)
    {
      for(i = 0; i < number_bed_layers; i++)
	{
	  if(self->bld[i] != NULL)
	    self->bld[i] = (double *)tl_free(sizeof(double), nsed, self->bld[i]);
	}
      self->bld = (double **)tl_free(sizeof(double *), number_bed_layers, self->bld);
    }
  if(self->cbp != NULL)
    {
      for(i = 0; i < ncti; i++)
	{
	  if(self->cbp[i] != NULL)
	    self->cbp[i] = (double *)tl_free(sizeof(double), 5, self->cbp[i]);
	}
      self->cbp = (double **)tl_free(sizeof(double *), ncti, self->cbp);
    }

  if(self->sum_snd_ald != NULL)
    self->sum_snd_ald = (double *)tl_free(sizeof(double), self->isize, self->sum_snd_ald);

  if(self->bed_character != NULL)
    self->bed_character = (int *)tl_free(sizeof(int), self->isize, self->bed_character);

  if(self->aldkf != NULL)
    self->aldkf = (int *)tl_free(sizeof(int), self->isize, self->aldkf);

  if(self->asdkf != NULL)
    self->asdkf = (int *)tl_free(sizeof(int), self->isize, self->asdkf);

  if(self->sum_nodal_source != NULL)
    self->sum_nodal_source = (double *)tl_free(sizeof(double), self->isize, self->sum_nodal_source);

  if(self->sum_snd_con != NULL)
    self->sum_snd_con = (double *)tl_free(sizeof(double), self->isize, self->sum_snd_con);

  if(self->source_coef != NULL)
    {
      for(i = 0; i < self->isize; i++)
	{
	  if(self->source_coef[i] != NULL)
	    self->source_coef[i] = (double *)tl_free(sizeof(double), nsed, self->source_coef[i]);
	}
      self->source_coef = (double **)tl_free(sizeof(double *), self->isize, self->source_coef);
    }

  if(self->prob_dep != NULL)
    {
      for(i = 0; i < self->isize; i++)
	{
	  if(self->prob_dep[i] != NULL)
	    self->prob_dep[i] = (double *)tl_free(sizeof(double), nsed, self->prob_dep[i]);
	}
      self->prob_dep = (double **)tl_free(sizeof(double *), self->isize, self->prob_dep);
    }

  if(self->bedload_eq_con != NULL)
    {
      for(i = 0; i < self->isize; i++)
	{
	  if(self->bedload_eq_con[i] != NULL)
	    self->bedload_eq_con[i] = (double *)tl_free(sizeof(double), nsed, self->bedload_eq_con[i]);
	}
      self->bedload_eq_con = (double **)tl_free(sizeof(double *), self->isize, self->bedload_eq_con);
    }

  if(self->bl_thick != NULL)
    {
      for(i = 0; i < self->isize; i++)
	{
	  if(self->bl_thick[i] != NULL)
	    self->bl_thick[i] = (double *)tl_free(sizeof(double), nsed, self->bl_thick[i]);
	}
      self->bl_thick = (double **)tl_free(sizeof(double *), self->isize, self->bl_thick);
    }

  if(self->bl_slc != NULL)
    {
      for(i = 0; i < self->isize; i++)
	{
	  if(self->bl_slc[i] != NULL)
	    self->bl_slc[i] = (double *)tl_free(sizeof(double), nsed, self->bl_slc[i]);
	}
      self->bl_slc = (double **)tl_free(sizeof(double *), self->isize, self->bl_slc);
    }
  if(self->critical_shear_vel != NULL)
    {
      for(i = 0; i < self->isize; i++)
	{
	  if(self->critical_shear_vel[i] != NULL)
	    self->critical_shear_vel[i] =
	      (double *)tl_free(sizeof(double), nsed, self->critical_shear_vel[i]);
	}
      self->critical_shear_vel = (double **)tl_free(sizeof(double *), self->isize, self->critical_shear_vel);
    }

  if(self->bl_vel != NULL)
    {
      for(i = 0; i < self->isize; i++)
	{
	  if(self->bl_vel[i] != NULL)
	    self->bl_vel[i] = (VECT2D *) tl_free(sizeof(VECT2D), nsand, self->bl_vel[i]);
	}
      self->bl_vel = (VECT2D **) tl_free(sizeof(VECT2D *), self->isize, self->bl_vel);
    }

  if(self->sum_blt != NULL)
    self->sum_blt = (double *)tl_free(sizeof(double), self->isize, self->sum_blt);

  if(self->bed_flag != NULL)
    self->bed_flag = (int *)tl_free(sizeof(int), self->isize, self->bed_flag);

  if(self->limiter != NULL)
    self->limiter = (int *)tl_free(sizeof(int), self->isize, self->limiter);

  if(self->limiter2 != NULL)
    self->limiter2 = (int *)tl_free(sizeof(int), self->isize, self->limiter2);

  if(self->conc_converge != NULL)
    self->conc_converge = (int *)tl_free(sizeof(int), self->isize, self->conc_converge);

  if(self->old_sed_displacement != NULL)
    self->old_sed_displacement = (double *)tl_free(sizeof(double), self->isize, self->old_sed_displacement);

  if(self->old_sed_active_stratum_ceiling != NULL)
    self->old_sed_active_stratum_ceiling =
      (double *)tl_free(sizeof(double), self->isize, self->old_sed_active_stratum_ceiling);

  if(self->old_sed_active_layer_distribution != NULL)
    {
      for(i = 0; i < nsed; i++)
	{
	  if(self->old_sed_active_layer_distribution[i] != NULL)
	    self->old_sed_active_layer_distribution[i] =
	      (double *)tl_free(sizeof(double), self->isize,
				self->old_sed_active_layer_distribution[i]);
	}
      self->old_sed_active_layer_distribution =
	(double **)tl_free(sizeof(double *), nsed, self->old_sed_active_layer_distribution);
    }

  if(self->old_sed_active_stratum_distribution != NULL)
    {
      for(i = 0; i < nsed; i++)
	{
	  if(self->old_sed_active_stratum_distribution[i] != NULL)
	    self->old_sed_active_stratum_distribution[i] =
	      (double *)tl_free(sizeof(double), self->isize,
				self->old_sed_active_stratum_distribution[i]);
	}
      self->old_sed_active_stratum_distribution =
	(double **)tl_free(sizeof(double *), nsed, self->old_sed_active_stratum_distribution);
    }

  if(self->old_it_ald != NULL)
    {
      for(i = 0; i < nsed; i++)
	{
	  if(self->old_it_ald[i] != NULL)
	    self->old_it_ald[i] = (double *)tl_free(sizeof(double), self->isize, self->old_it_ald[i]);
	}
      self->old_it_ald = (double **)tl_free(sizeof(double *), nsed, self->old_it_ald);
    }
  if(self->old_it_rhs != NULL)
    {
      for(i = 0; i < nsed + 1; i++)
	{
	  if(self->old_it_rhs[i] != NULL)
	    self->old_it_rhs[i] = (double *)tl_free(sizeof(double), self->isize, self->old_it_rhs[i]);
	}
      self->old_it_rhs = (double **)tl_free(sizeof(double *), nsed + 1, self->old_it_rhs);
    }

  if(self->old_it_concentration != NULL)
    {
      for(i = 0; i < nsed; i++)
	{
	  if(self->old_it_concentration[i] != NULL)
	    self->old_it_concentration[i] =
	      (double *)tl_free(sizeof(double), self->isize, self->old_it_concentration[i]);
	}
      self->old_it_concentration =
	(double **)tl_free(sizeof(double *), nsed,self->old_it_concentration);
    }

  if(self->nodal_grad != NULL)
    self->nodal_grad = (VECT2D *) tl_free(sizeof(VECT2D), self->isize, self->nodal_grad);

  if(self->elev_disp != NULL)
    self->elev_disp = (double *)tl_free(sizeof(double), self->isize, self->elev_disp);

  if(self->fraction != NULL)
    self->fraction = (double *)tl_free(sizeof(double), nsed, self->fraction);

  if(self->old_fraction != NULL)
    self->old_fraction = (double *)tl_free(sizeof(double), nsed, self->old_fraction);

  if(self->AL_thickness != NULL)
    self->AL_thickness = (double *)tl_free(sizeof(double), self->isize, self->AL_thickness);

  if(self->old_sed_active_layer_porosity != NULL)
    self->old_sed_active_layer_porosity =
      (double *)tl_free(sizeof(double), self->isize, self->old_sed_active_layer_porosity);

  if(self->old_sed_active_layer_critical_erosion_shear != NULL)
    self->old_sed_active_layer_critical_erosion_shear =
      (double *)tl_free(sizeof(double), self->isize, self->old_sed_active_layer_critical_erosion_shear);

  if(self->old_sed_active_layer_erosion_rate_constant != NULL)
    self->old_sed_active_layer_erosion_rate_constant =
      (double *)tl_free(sizeof(double), self->isize, self->old_sed_active_layer_erosion_rate_constant);

  if(self->old_sed_active_layer_erosion_rate_exponent != NULL)
    self->old_sed_active_layer_erosion_rate_exponent =
      (double *)tl_free(sizeof(double), self->isize, self->old_sed_active_layer_erosion_rate_exponent);

  if(self->sum_sed_decay != NULL)
    self->sum_sed_decay = (double *)tl_free(sizeof(double), self->isize, self->sum_sed_decay);
  return 0;
}

/** Call the ADH Newton solves for all the models. This will kick out
    with a failure for any model so the caller needs handle nonlinear
    solver failures */
int ADH_OneLevelTransport_solve(ADH_OneLevelTransport* self)
{
  int i, j, k, ie, n;		/* loop counter */
  int itrn;			/* loop counter over the transported quantities */
  int nsys;			/* the number of systems being solved */
  int nsys_sq;			/* the number of systems being solved squared */
  /* int isize_prev = 0;		/\* the previous isize *\/ */
  int sed_iteration;		/* the number of iterations done to balance the fractions in active layer */
  int sed_number;		/* given the number of clay or silt what is its number in sediment */
  double ald_check;		/* maximum change in the active layer distribution since the last iteration */
  double ald_max = DONE;
  double t_sed;			/* the real time for the sediment current sediment step */
  double flowrate_x;		/* the flowrate for the hydro results to be kept for the sed timesteps */
  double flowrate_y;		/* the flowrate for the hydro results to be kept for the sed timesteps */
  double ald_diff;
  int max_node;
  double sed_it_tol;
  double bed_thick_old;
  double bed_thick_new;
  double esrh;			/* the equivalent sand roughness height */
  double blsg;
  double albd;
  double asbd;
  VECT2D tvel;			/* the local transverse velocity components */
  int count = 0;		/* counter for sed time step */
  int istr;			/* counter for ceq boundary condition */
#ifdef _MESSG
  double sed_time_start;
  double sed_time_end;
  double sed_time_total;
#endif
  double qus_depth_factor = 0.1;
  double qus_max_change;
  /* double sdtemp; */
  /* double temp_dpl; */
  /* double ust; */
  /* double blv; */
  double blvelmag;
  /* reset the matrix */
  realloc_fe_main(self);
  /* ************************** */
  /* SOLVES GROUNDWATER PROBLEM */
  if(GW_FLOW)
    {

      nsys = 1;
      nsys_sq = 1;

      /* solves the groundwater flow problem - returns from fe_main with a NO if 
         the nonlinear iterations failed */
      if(myid <= 0)
	printf("\n time step: groundwater flow. (time = %f, %3.0f %% ), dt = %15.6e\n",
	       t_prev, 100. * t_prev / tfinal, dt);
      if(fe_newton
	 (self->bc_mask, self->matrix, self->diagonal, self->scale_vect, self->residual, self->sol, nnode, my_nnode, nsys,
	  nsys_sq, 0, fe_gw_init, fe_gw_update, fe_gw_resid, fe_gw_load, fe_gw_inc) == NO)
	return 1;

      /* update the elemental velocities and saturations for use by the transport functions */
      fe_gw_elem_sat_vel();

      /* solves the transport problems - returns from fe_main with a NO if 
         the nonlinear iterations failed */
      if(myid <= 0 && ntransport > 0)
	{
	  printf(" \n");
	  printf("time step: GW TRANSPORT \n");
	}
      for(itrn = 0; itrn < ntransport; itrn++)
	{
	  if(myid <= 0)
	    printf("\n Transporting Constituent %d\n", itrn + 1);
	  if(fe_newton
	     (self->bc_mask, self->matrix, self->diagonal, self->scale_vect, self->residual, self->sol, nnode, my_nnode, nsys,
	      nsys_sq, itrn, fe_trns_init, fe_trns_update, fe_trns_resid, fe_trns_load,
	      fe_trns_inc) == NO)
	    return 1;

	  /* operator split reactions and sorption */
	}
    }
  /* ************************** */
  /* SOLVES HEAT TRANSPORT PROBLEM */
  if(HEAT)
    {

      nsys = 1;
      nsys_sq = 1;

      /* solves the heat transport problem - returns from fe_main with a NO if
         the nonlinear iterations failed */
      if(myid <= 0)
	printf("\n time step: heat transport. (time = %f, %3.0f %% ), dt = %15.6e\n",
	       t_prev, 100. * t_prev / tfinal, dt);
      if(fe_newton
	 (self->bc_mask, self->matrix, self->diagonal, self->scale_vect, self->residual, self->sol, nnode, my_nnode, nsys,
	  nsys_sq, 0, fe_heat_init, fe_heat_update, fe_heat_resid, fe_heat_load,
	  fe_heat_inc) == NO)
	return 1;
    }
  /* ************************** */
  /* SOLVES DIFFUSIVE WAVE OVERLAND FLOW PROBLEM */
  if(OL_FLOW)
    {

      nsys = 1;
      nsys_sq = 1;

      /* solves the groundwater flow problem - returns from fe_main with a NO if
         the nonlinear iterations failed */
      if(myid <= 0)
	printf("\n time step: dwave overland flow. (time = %f, %3.0f %% ), dt = %15.6e\n",
	       t_prev, 100. * t_prev / tfinal, dt);

      /* solves the channel flow problem */
      /*if(fe_newton(self->bc_mask, self->matrix, self->diagonal, self->scale_vect, self->residual, self->sol,
         nnode, my_nnode, nsys, nsys_sq, 0,
         fe_ch_init, fe_ch_update, fe_ch_resid, fe_ch_load, fe_ch_inc) == NO) */

      /* copies the solutions to the old solutions */
      /*for(i=0;i<nnode;i++)
         old_ch_head[i] = ch_head[i]; */
      /* the "ch" routines have been commented out, but will be used eventually */
      /* BR per RCB, 6/04 */

      /* solves the diffusive wave equation for overland flow */
      /* returns from fe_main with a NO if the nonlinear iterations failed */

      if(fe_newton
	 (self->bc_mask, self->matrix, self->diagonal, self->scale_vect, self->residual, self->sol, nnode, my_nnode, nsys,
	  nsys_sq, 0, fe_ol_init, fe_ol_update, fe_ol_resid, fe_ol_load, fe_ol_inc) == NO)
	return 1;

      /* solves the transport problems - returns from fe_main with a NO if 
         the nonlinear iterations failed */
      if(myid <= 0 && OL_TRANSPORT)
	{
	  printf(" \n");
	  printf("time step: OL TRANSPORT \n");
	}
      for(itrn = 0; itrn < ntransport; itrn++)
	if(fe_newton
	   (self->bc_mask, self->matrix, self->diagonal, self->scale_vect, self->residual, self->sol, nnode, my_nnode, nsys,
	    nsys_sq, itrn, fe_trns_init, fe_trns_update, fe_trns_resid, fe_trns_load,
	    fe_trns_inc) == NO)
	  return 1;

      /* operator split reactions and sorption */

    }
  /* ******************************** */
  /* SOLVES THE NAVIER STOKES PROBLEM */
  if(NS_FLOW)
    {
      /* solves the navier stokes flow problem - returns to fe_main with a NO if 
         the nonlinear iterations failed */

      nsys = 4;
      nsys_sq = 16;
      refresh_pc = YES;

      if(myid <= 0)
	printf("\n time step: Navier Stokes flow. (time = %f, %3.0f %% ), dt = %15.6e\n",
	       t_prev, 100. * t_prev / tfinal, dt);

      if(fe_newton
	 (self->bc_mask, self->matrix, self->diagonal, self->scale_vect, self->residual, self->sol, nnode, my_nnode, nsys,
	  nsys_sq, 0, fe_ns_init, fe_ns_update, fe_ns_resid, fe_ns_load, fe_ns_inc) == NO)
	return 1;

      fe_press_adjust();

      if(MOVING_GRID)
	{
	  if(myid <= 0)
	    printf("move the grid\n");

	  nsys = 1;
	  nsys_sq = 1;

	  if(fe_newton
	     (self->bc_mask, self->matrix, self->diagonal, self->scale_vect, self->residual, self->sol, nnode, my_nnode, nsys,
	      nsys_sq, 0, fe_mg_init, fe_mg_update, fe_mg_resid, fe_mg_load,
	      fe_mg_inc) == NO)
	    return 1;

	  nsys = 4;
	  nsys_sq = 16;
	}
      if(NS_NFS)
	{
	  if(myid <= 0)
	    printf("calculate the velocity gradients\n");

	  nsys = 1;
	  nsys_sq = 1;
	  refresh_pc = YES;

	  if(fe_newton
	     (self->bc_mask, self->matrix, self->diagonal, self->scale_vect, self->residual, self->sol, nnode, my_nnode, nsys,
	      nsys_sq, 0, fe_ns_velgrad_init, fe_ns_u_grad_update, fe_ns_dudx_resid,
	      fe_ns_velgrad_load, fe_ns_dudx_inc) == NO)
	    return 1;
	  refresh_pc = NO;
	  if(fe_newton
	     (self->bc_mask, self->matrix, self->diagonal, self->scale_vect, self->residual, self->sol, nnode, my_nnode, nsys,
	      nsys_sq, 0, fe_ns_velgrad_init, fe_ns_u_grad_update, fe_ns_dudy_resid,
	      fe_ns_velgrad_load, fe_ns_dudy_inc) == NO)
	    return 1;
	  refresh_pc = NO;
	  if(fe_newton
	     (self->bc_mask, self->matrix, self->diagonal, self->scale_vect, self->residual, self->sol, nnode, my_nnode, nsys,
	      nsys_sq, 0, fe_ns_velgrad_init, fe_ns_u_grad_update, fe_ns_dudz_resid,
	      fe_ns_velgrad_load, fe_ns_dudz_inc) == NO)
	    return 1;
	  if(fe_newton
	     (self->bc_mask, self->matrix, self->diagonal, self->scale_vect, self->residual, self->sol, nnode, my_nnode, nsys,
	      nsys_sq, 0, fe_ns_velgrad_init, fe_ns_v_grad_update, fe_ns_dvdx_resid,
	      fe_ns_velgrad_load, fe_ns_dvdx_inc) == NO)
	    return 1;
	  if(fe_newton
	     (self->bc_mask, self->matrix, self->diagonal, self->scale_vect, self->residual, self->sol, nnode, my_nnode, nsys,
	      nsys_sq, 0, fe_ns_velgrad_init, fe_ns_v_grad_update, fe_ns_dvdy_resid,
	      fe_ns_velgrad_load, fe_ns_dvdy_inc) == NO)
	    return 1;
	  if(fe_newton
	     (self->bc_mask, self->matrix, self->diagonal, self->scale_vect, self->residual, self->sol, nnode, my_nnode, nsys,
	      nsys_sq, 0, fe_ns_velgrad_init, fe_ns_v_grad_update, fe_ns_dvdz_resid,
	      fe_ns_velgrad_load, fe_ns_dvdz_inc) == NO)
	    return 1;
	  if(fe_newton
	     (self->bc_mask, self->matrix, self->diagonal, self->scale_vect, self->residual, self->sol, nnode, my_nnode, nsys,
	      nsys_sq, 0, fe_ns_velgrad_init, fe_ns_w_grad_update, fe_ns_dwdx_resid,
	      fe_ns_velgrad_load, fe_ns_dwdx_inc) == NO)
	    return 1;
	  if(fe_newton
	     (self->bc_mask, self->matrix, self->diagonal, self->scale_vect, self->residual, self->sol, nnode, my_nnode, nsys,
	      nsys_sq, 0, fe_ns_velgrad_init, fe_ns_w_grad_update, fe_ns_dwdy_resid,
	      fe_ns_velgrad_load, fe_ns_dwdy_inc) == NO)
	    return 1;
	  if(fe_newton
	     (self->bc_mask, self->matrix, self->diagonal, self->scale_vect, self->residual, self->sol, nnode, my_nnode, nsys,
	      nsys_sq, 0, fe_ns_velgrad_init, fe_ns_w_grad_update, fe_ns_dwdz_resid,
	      fe_ns_velgrad_load, fe_ns_dwdz_inc) == NO)
	    return 1;

	  refresh_pc = YES;
	  nsys = 4;
	  nsys_sq = 16;
	}

      /* solves the transport problems - returns from fe_main with a NO if 
         the nonlinear iterations failed */
      nsys = 1;
      nsys_sq = 1;

      if(myid <= 0)
	{
	  printf(" \n");
	  printf("time step: NS TRANSPORT \n");
	}
      /* this section is solving transport for all constituents that are not sediment */
      for(itrn = 0; itrn < ntransport; itrn++)
	{
	  if(con[itrn].type != SND && con[itrn].type != CLA)
	    {
	      if(fe_newton
		 (self->bc_mask, self->matrix, self->diagonal, self->scale_vect, self->residual, self->sol, nnode, my_nnode,
		  nsys, nsys_sq, itrn, fe_nstrns_init, fe_nstrns_update, fe_nstrns_resid,
		  fe_nstrns_load, fe_nstrns_inc) == NO)
		return 1;
	    }
	}
      /* solves runge-kutta for the reactions */
      if(transport.tke >= 0)
	tr_3d_reaction(nnode, my_nnode, self->nodal_production, self->nodal_volume, self->nodal_viscosity,
		       concentration[transport.tke], concentration[transport.tds]);

      /*  JNT 7-21-03 ****************************** */
      /* calculate the bed change */
      if(SEDIMENT)
	{
	  for(k = 0; k < nsed; k++)
	    {
	      for(i = 0; i < nnode; i++)
		{
		  active_layer_distribution[k][i] = old_active_layer_distribution[k][i];
		}
	    }

	  /* zero the arrays */

	  for(i = 0; i < nnode; i++)
	    {
	      self->sum_blt[i] = 0.0;
	      nodal_source[i] = 0.;
	      nodal_decay_coefficient[i] = 0.;
	      for(j = 0; j < nsed; j++)
		{
		  self->source_coef[i][j] = 0.0;
		  self->prob_dep[i][j] = 0.0;
		}
	      self->bed_character[i] = 0;
	      self->bed_flag[i] = 0;
	      for(j = 0; j < number_bed_layers; j++)
		{
		  self->sum_blt[i] += bed_layer_thickness[j][i];
		}
	    }

	  sed_iteration = -1;
	  do
	    {
	      sed_iteration++;
	      printf("sed_iteration = %d\n", sed_iteration);
	      for(i = 0; i < nnode; i++)
		{
		  self->sum_snd_ald[i] = 0.0;
		  for(k = 0; k < nsand; k++)
		    {
		      self->sum_snd_ald[i] += active_layer_distribution[k][i];
		    }
		  if(self->sum_snd_ald[i] > 0.9)
		    {
		      self->bed_character[i] = 1;
		      /* 1 for a sand bed, 2 for a clay bed */
		    }
		  else
		    {
		      self->bed_character[i] = 2;
		    }
		}

	      if(nsed > 0)
		{
		  sed_snd_source_decay_ns(self->source_coef, rouse_coef, active_layer_thickness);
		}
	      for(i = 0; i < nnode; i++)
		{
		  self->sum_nodal_source[i] = 0.0;
		}

	      for(itrn = 0; itrn < ntransport; itrn++)
		{
		  for(i = 0; i < nnode; i++)
		    {
		      if(self->bed_character[i] == 1)
			{
			  if(con[itrn].type == SND)
			    {
			      k = transport_to_sand[itrn];
			      nodal_source[i] =
				-sand_settling_velocity[k] * sed_specific_gravity[k] *
				active_layer_distribution[k][i] * self->source_coef[i][k] /
				con[itrn].property[0];
			      self->sum_nodal_source[i] += nodal_source[i];
			      nodal_decay_coefficient[i] =
				sand_settling_velocity[k] * rouse_coef[k][i];
			    }

			  /* now that the source and decay have been determined for the sand, get the fraction for the clays */
			  if(con[itrn].type == CLA)
			    {
			      j = transport_to_clay[itrn];
			      sed_number = j + nsand + nsilt;
			      nodal_source[i] =
				self->sum_nodal_source[i] *
				(active_layer_distribution[sed_number][i] / self->sum_snd_ald[i]);
			      nodal_decay_coefficient[i] =
				rouse_coef[sed_number][i] * clay_settling_velocity[j];
			    }
			}
		      if(self->bed_character[i] == 2)
			{
			  if(con[itrn].type == CLA)
			    {
			      k = transport_to_clay[itrn];
			      sed_number = nsand + nsilt + k;
			      if((self->sum_blt[i] + old_displacement[i]) <= 0.0)
				{
				  self->source_coef[i][sed_number] = 0.0;
				}
			      nodal_source[i] = -self->source_coef[i][sed_number];
			      self->sum_nodal_source[i] += nodal_source[i];
			      nodal_decay_coefficient[i] =
				rouse_coef[sed_number][i] * clay_settling_velocity[k];
			    }
			  if(con[itrn].type == SND)
			    {
			      j = transport_to_sand[itrn];
			      sed_number = j;
			      nodal_source[i] =
				self->sum_nodal_source[i] *
				(active_layer_distribution[sed_number][i] /
				 (1.0 - self->sum_snd_ald[i]));
			      nodal_decay_coefficient[i] =
				rouse_coef[sed_number][i] * sand_settling_velocity[j];
			    }
			}
		    }

		  if(fe_newton
		     (self->bc_mask, self->matrix, self->diagonal, self->scale_vect, self->residual, self->sol, nnode, my_nnode,
		      nsys, nsys_sq, itrn, fe_nstrns_init, fe_nstrns_update,
		      fe_nstrns_resid, fe_nstrns_load, fe_nstrns_inc) == NO)
		    return 1;
		}
	      /* sed_snd_init_ns(eq_bed_load, self->source_coef, rouse_coef, ustar,
	         self->bed_character, active_layer_thickness);

	         fe_sed_activechange_ns(eq_bed_load, self->source_coef, rouse_coef, ustar,
	         self->bed_character); */

	      for(i = 0; i < nnode; i++)
		{
		  self->aldkf[i] = 0;
		  self->asdkf[i] = 0;
		}

	      sed_active_layer_distribution_balance(self->aldkf, self->asdkf);
	    }
	  while(ald_max > 1.0e-8 && sed_iteration < 4);

	}

      /******************************************************************************************************/

      nsys = 4;
      nsys_sq = 16;

      /* copies the solutions to the old solutions since the 
         time step has made it to the end and is therefore to be kept */
      for(i = 0; i < nnode; i++)
	{
	  old_vel[i].x = vel[i].x;
	  old_vel[i].y = vel[i].y;
	  old_vel[i].z = vel[i].z;
	  old_prs[i] = prs[i];
	  for(itrn = 0; itrn < ntransport; itrn++)
	    {
	      old_concentration[itrn][i] = concentration[itrn][i];
	    }
	}
      if(MOVING_GRID)
	{

      /************************JNT 1-23-04********************************************/
	  if(SEDIMENT)
	    {
	      /* if I pass the border node values for the active_stratum_ceiling to the other processors */
	      /* then each of them should be able to calculate the bed_thickness and bed distribution on */
	      /* their own.  this is why the sed_bed_layer_update loops on nnode rather than my_nnode.   */

#ifdef _MESSG
	      comm_update_double(active_stratum_ceiling, ONE);
#endif

	      for(i = 0; i < nnode; i++)
		{
		  for(k = 0; k < nsed; k++)
		    {
		      self->ald[k] = active_layer_distribution[k][i];
		      self->asd[k] = active_stratum_distribution[k][i];
		      for(j = 0; j < number_bed_layers; j++)
			{
			 self-> bld[j][k] = bed_layer_distribution[j][k][i];
			}
		    }
		  for(j = 0; j < number_bed_layers; j++)
		    {
		      self->blt[j] = bed_layer_thickness[j][i];
		      self->bl_por[j] = bed_layer_bulk_density[j][i];
		      self->bl_ces[j] = bed_layer_critical_erosion_shear[j][i];
		      self->bl_erc[j] = bed_layer_erosion_rate_constant[j][i];
		      self->bl_ere[j] = bed_layer_erosion_rate_exponent[j][i];
		    }
		  blsg = 0.0;
		  for(k = 0; k < nsed; k++)
		    {
		      blsg += self->ald[k] * sed_specific_gravity[k];
		    }
		  albd =
		    density * (active_layer_porosity[i] +
			       (1. - active_layer_porosity[i]) * blsg);
		  blsg = 0.0;
		  for(k = 0; k < nsed; k++)
		    {
		      blsg += self->ald[k] * sed_specific_gravity[k];
		    }
		  asbd =
		    density * (active_stratum_porosity[i] +
			       (1. - active_stratum_porosity[i]) * blsg);
		  sed_bed_layer_update(old_displacement[i], displacement[i],
				       active_stratum_ceiling[i],
				       self->old_sed_active_stratum_ceiling[i], self->bed_character[i],
				       number_bed_layers, nsed, self->bld, self->blt, self->ald, self->asd, albd,
				       active_layer_critical_erosion_shear[i],
				       active_layer_erosion_rate_constant[i],
				       active_layer_erosion_rate_exponent[i], asbd,
				       active_stratum_critical_erosion_shear[i],
				       active_stratum_erosion_rate_constant[i],
				       active_stratum_erosion_rate_exponent[i], self->bl_por,
				       self->bl_ces, self->bl_erc, self->bl_ere);
		  for(k = 0; k < nsed; k++)
		    {
		      active_layer_distribution[k][i] = self->ald[k];
		      active_stratum_distribution[k][i] = self->asd[k];
		      for(j = 0; j < number_bed_layers; j++)
			{
			  bed_layer_distribution[j][k][i] = self->bld[j][k];
			}
		    }
		  for(j = 0; j < number_bed_layers; j++)
		    {
		      bed_layer_thickness[j][i] = self->blt[j];
		      bed_layer_bulk_density[j][i] = self->bl_por[j];
		      bed_layer_critical_erosion_shear[j][i] = self->bl_ces[j];
		      bed_layer_erosion_rate_constant[j][i] = self->bl_erc[j];
		      bed_layer_erosion_rate_exponent[j][i] = self->bl_ere[j];
		    }
		}

	      for(i = 0; i < nnode; i++)
		{
		  old_active_stratum_ceiling[i] = active_stratum_ceiling[i];
		}
	      for(k = 0; k < nsed; k++)
		{
		  for(i = 0; i < nnode; i++)
		    {
		      old_active_layer_distribution[k][i] = active_layer_distribution[k][i];
		    }
		}
	    }
      /*****************************************************************/
	  for(i = 0; i < nnode; i++)
	    {
	      old_displacement[i] = displacement[i];
	      old_grid_speed[i] = grid_speed[i];
	    }

	}
    }
  if(SW2_FLOW)
    {
      /* solves the overland flow problem - returns from fe_main with a NO if 
         the nonlinear iterations failed */
      /*mwf hack see if setting refresh_pc here matters, doesn't seem to*/
      refresh_pc = YES;
      

      if(myid <= 0)
	{
	  printf("\n\n time step: 2D Hydrodynamics. (time = %f, %3.0f %% ), dt = %15.6e\n",
		 t_prev, 100. * t_prev / tfinal, dt);
	}

      nsys = 3;
      nsys_sq = 9;
      if(BOAT)
	boat_move_fields((t_prev + dt), node, nnode, gravity, prs);
      if(ICE)
	for(i = 0; i < nstring; i++)
	  if(str_values[i].ice_string == STR_ICE)
	    fe_sw2_move_ice(i);

      if(STEADY_STATE_HYDRO)
	{
	  if(NEW_STEADY_HYDRO)
	    {
	      if(myid <= 0)
		printf("NEW QUASI-STEADY-STATE SOLUTION AT TIME: %f\n", t_prev);
	      if(fe_newton_std
		 (self->bc_mask, self->matrix, self->diagonal, self->scale_vect, self->residual, self->sol, nnode, my_nnode,
		  nsys, nsys_sq, 0, fe_sw2_init, fe_sw2_update, fe_sw2_resid, fe_sw2_load,
		  fe_sw2_inc) == NO)
		tl_error("Steady State Hydro Not Reached For Quasi-Steady-State Run.\n");
	    }
	}
      else
	{
	  if(CONVEYANCE)
	    fe_sw2_conveyance();

	  if(fe_newton
	     (self->bc_mask, self->matrix, self->diagonal, self->scale_vect, self->residual, self->sol, nnode, my_nnode, nsys,
	      nsys_sq, 0, fe_sw2_init, fe_sw2_update, fe_sw2_resid, fe_sw2_load,
	      fe_sw2_inc) == NO)
	    return 1;
	}

      if(BOAT_ENT)
	boat_get_stresses((t_prev + dt), node, nnode, gravity, fe_sw2_centroid_values,
			  ol_vel, ol_head, vessel_stress);

      /* solves the transport problems - returns from fe_main with a NO if 
         the nonlinear iterations failed */

      nsys = 1;
      nsys_sq = 1;
      if(ntransport > 0)
	{
	  if(myid <= 0)
	    printf("\ntime step: 2D TRANSPORT \n");
	  for(i = 0; i < nnode; i++)
	    {
	      nodal_source[i] = 0.;
	      nodal_decay_coefficient[i] = 0.;
	    }
	}

      for(itrn = 0; itrn < ntransport; itrn++)
	{
	  if(con[itrn].type != SND && con[itrn].type != CLA)
	    {
      /***Set dt_sed = dt since it isn't incremented as with sediment  Can this be done in the time routine???***/
	      dt_sed = dt;
	      printf("dt = %f   dt_sed = %f\n", dt, dt_sed);
	      if(fe_newton
		 (self->bc_mask, self->matrix, self->diagonal, self->scale_vect, self->residual, self->sol, nnode, my_nnode,
		  nsys, nsys_sq, itrn, fe_sw2trns_init, fe_sw2trns_update, fe_sw2trns_resid,
		  fe_sw2trns_load, fe_sw2trns_inc) == NO)
		return 1;
	    }
	}
      if(SW2_NFS)		/* velocity gradients for SW2 GS 12/13/2006 */
	{
	  if(myid <= 0)
	    nsys = 1;
	  nsys_sq = 1;
	  refresh_pc = YES;
	  if(fe_newton
	     (self->bc_mask, self->matrix, self->diagonal, self->scale_vect, self->residual, self->sol, nnode, my_nnode, nsys,
	      nsys_sq, 0, fe_sw2_velgrad_init, fe_sw2_u_grad_update, fe_sw2_dudx_resid,
	      fe_sw2_velgrad_load, fe_sw2_dudx_inc) == NO)
	    return 1;
	  refresh_pc = NO;
	  if(fe_newton
	     (self->bc_mask, self->matrix, self->diagonal, self->scale_vect, self->residual, self->sol, nnode, my_nnode, nsys,
	      nsys_sq, 0, fe_sw2_velgrad_init, fe_sw2_u_grad_update, fe_sw2_dudy_resid,
	      fe_sw2_velgrad_load, fe_sw2_dudy_inc) == NO)
	    return 1;
	  refresh_pc = NO;
	  if(fe_newton
	     (self->bc_mask, self->matrix, self->diagonal, self->scale_vect, self->residual, self->sol, nnode, my_nnode, nsys,
	      nsys_sq, 0, fe_sw2_velgrad_init, fe_sw2_v_grad_update, fe_sw2_dvdx_resid,
	      fe_sw2_velgrad_load, fe_sw2_dvdx_inc) == NO)
	    return 1;
	  refresh_pc = NO;
	  if(fe_newton
	     (self->bc_mask, self->matrix, self->diagonal, self->scale_vect, self->residual, self->sol, nnode, my_nnode, nsys,
	      nsys_sq, 0, fe_sw2_velgrad_init, fe_sw2_v_grad_update, fe_sw2_dvdy_resid,
	      fe_sw2_velgrad_load, fe_sw2_dvdy_inc) == NO)
	    return 1;

	  refresh_pc = YES;

	  nsys = 3;
	  nsys_sq = 9;
	  /*for (i = 0; i < nnode; i++)
	     { 
	     printf(" dudx is %f\n", u_grad[i].x);
	     printf("dudy is %f\n", u_grad[i].y);
	     printf("dvdux is %f\n",v_grad[i].x);
	     printf("dvdy us %f\n", v_grad[i].y);
	     } */

	  /* Rate of strain tennsor GS 12/27/2006 */

	  /* for ( i = 0; i < nnode; i++)
	     {
	     vel2dg[i].xx = u_grad[i].x;
	     printf("Def tensor x is %f\n", vel2dg[i].xx);
	     vel2dg[i].xy = 0.5 * (u_grad[i].y + v_grad[i].x);
	     vel2dg[i].yy = v_grad[i].y;
	     vort[i] = (v_grad[i].x - u_grad[i].y);
	     } */

	}

      if(SEDIMENT)
	{
	  for(i = 0; i < nnode; i++)
	    {
	      old_sed_ol_head[i] = old_ol_head[i];
	      self->old_sed_displacement[i] = old_displacement[i];
	      self->old_sed_active_stratum_ceiling[i] = old_active_stratum_ceiling[i];
	      for(n = 0; n < number_bed_layers; n++)
		{
		  bed_layer_thickness[n][i] = old_bed_layer_thickness[n][i];
		}
	    }
	  for(k = 0; k < nsed; k++)
	    {
	      for(i = 0; i < nnode; i++)
		{
		  if(k < nsand)
		    {
		      j = sand_to_transport[k];
		    }
		  else
		    {
		      j = clay_to_transport[k - nsand];
		    }
		  active_layer_distribution[k][i] = old_active_layer_distribution[k][i];
		  self->old_sed_active_layer_distribution[k][i] =
		    old_active_layer_distribution[k][i];
		  old_sed_concentration[k][i] = old_concentration[j][i];
		  old_sed_bedload_con[k][i] = old_bedload_con[j][i];
		  old_sed_bedload_thick[k][i] = old_bedload_thick[j][i];
		  active_stratum_distribution[k][i] = old_active_stratum_distribution[k][i];
		  self->old_sed_active_stratum_distribution[k][i] =
		    old_active_stratum_distribution[k][i];
		  for(n = 0; n < number_bed_layers; n++)
		    {
		      bed_layer_distribution[n][k][i] = old_bed_layer_distribution[n][k][i];
		    }
		}
	    }
	  for(i = 0; i < nnode; i++)
	    {
	      active_layer_porosity[i] = old_active_layer_porosity[i];
	      self->old_sed_active_layer_porosity[i] = old_active_layer_porosity[i];
	      active_layer_critical_erosion_shear[i] =
		old_active_layer_critical_erosion_shear[i];
	      self->old_sed_active_layer_critical_erosion_shear[i] =
		old_active_layer_critical_erosion_shear[i];
	      active_layer_erosion_rate_constant[i] =
		old_active_layer_erosion_rate_constant[i];
	      self->old_sed_active_layer_erosion_rate_constant[i] =
		old_active_layer_erosion_rate_constant[i];
	      active_layer_erosion_rate_exponent[i] =
		old_active_layer_erosion_rate_exponent[i];
	      self->old_sed_active_layer_erosion_rate_exponent[i] =
		old_active_layer_erosion_rate_exponent[i];
	    }

	  /* calculate the bed change */
	  /* Calculate the transport timestep size and loop within the input dt */
	  if(fabs(dt - dt_max) < SMALL8)
	    {
	      dt_sed = dt_sed_max;
	    }
	  else
	    {
	      dt_sed = dt;
	      do
		{
		  count += 1;
		  if(dt_sed <= dt_sed_max && dt_sed <= dt)
		    {
		      dt_sed = dt_sed;
		    }
		  else
		    {
		      dt_sed = dt_sed / 2;
		    }
		}
	      while(count < 10);
	    }

	  t_sed = t_prev;
	  do
	    {
	      t_sed = t_sed + dt_sed;
	      printf("\n sed_time = %e\n", t_sed);
#ifdef _MESSG
	      sed_time_start = MPI_Wtime();
#endif
	      /* zero the arrays */
	      for(i = 0; i < nnode; i++)
		{
		  self->limiter[i] = 0;
		  self->sum_blt[i] = 0.0;
		  self->sum_snd_ald[i] = 0.0;
		  nodal_source[i] = 0.;
		  nodal_decay_coefficient[i] = 0.;
		  for(j = 0; j < nsed; j++)
		    {
		      /* really only want to zero these first timestep */
		      self->source_coef[i][j] = 0.0;
		      self->prob_dep[i][j] = 0.0;
		      self->old_sed_ald[j] = self->old_sed_active_layer_distribution[j][i];
		    }
		  self->bed_character[i] = 0;
		  for(j = 0; j < number_bed_layers; j++)
		    {
		      self->sum_blt[i] += bed_layer_thickness[j][i];
		    }
		  for(k = 0; k < nsand; k++)
		    {
		      self->sum_snd_ald[i] += active_layer_distribution[k][i];
		    }
		  if(nclay > 0)
		    {		/* only distinguish when clays are present JNT 8-27 */
		      if(self->sum_snd_ald[i] > 0.9)
			{
			  self->bed_character[i] = 1;	/* 1 for a sand bed, 2 for a clay bed */
			}
		      else
			{
			  self->bed_character[i] = 2;
			}
		    }
		  else
		    {
		      self->bed_character[i] = 1;
		    }
		}
	      sed_iteration = 0;
/*** SEDIMENT ITERATION LOOP ***/

	      do
		{
		  sed_iteration++;
		  printf("sed_iteration = %d\n", sed_iteration);
		  ald_check = 0.0;
		  sed_it_tol = 1.0;

		  for(i = 0; i < nnode; i++)
		    {
		      self->limiter2[i] = 0;
		      self->sum_snd_ald[i] = 0.0;
		      for(k = 0; k < nsand; k++)
			{
			  self->sum_snd_ald[i] += self->old_sed_active_layer_distribution[k][i];
			}

		      /* calculate the nodal bed elevation */
		      self->elev_disp[i] = bed_elev[i] + self->old_sed_displacement[i];
		    }
		  /* calculate the nodal bed gradient to be used later, this is a VECT2D[nnode] */
		  vt_2d_grad_nodal_avg(self->elev_disp, self->nodal_grad);
      /*** this update is in ***/
#ifdef _MESSG
		  comm_update_VECT2D(self->nodal_grad);
#endif

#ifdef _MESSG
		  comm_update_double(active_layer_thickness, ONE);
#endif

		  for(i = 0; i < nnode; i++)
		    {
		      self->sum_nodal_source[i] = 0.0;
		    }

		  sed_source_decay(self->source_coef, rouse_coef, self->prob_dep, self->critical_shear_vel,
				   self->nodal_grad, self->old_sed_displacement, self->ustar,
				   active_layer_critical_erosion_shear,
				   active_layer_erosion_rate_constant,
				   active_layer_erosion_rate_exponent, sed_iteration);

		  sed_bedload_source_decay(self->bl_vel, self->bedload_eq_con, self->bl_thick, self->bl_slc, self->ustar,
					   self->bed_character, self->critical_shear_vel, self->nodal_grad,
					   self->old_sed_active_layer_distribution,
					   self->old_sed_displacement,
					   active_layer_critical_erosion_shear,
					   sed_iteration);

		  for(i = 0; i < nnode; i++)
		    {
		      for(itrn = 0; itrn < ntransport; itrn++)
			{
			  if(con[itrn].type == SND || con[itrn].type == CLA)
			    {
			      if(con[itrn].type == SND)
				{
				  j = transport_to_sand[itrn];
				  sed_number = j;
				  bedload_thick[itrn][i] = self->bl_thick[i][j];
				  bedload_vel[itrn][i].x = self->bl_vel[i][j].x;
				  bedload_vel[itrn][i].y = self->bl_vel[i][j].y;
				}
			      if(con[itrn].type == CLA)
				{
				  j = transport_to_clay[itrn];
				  sed_number = j + nsand + nsilt;
				  bedload_thick[itrn][i] = self->bl_thick[i][j];
				  bedload_vel[itrn][i].x = self->bl_vel[i][j].x;
				  bedload_vel[itrn][i].y = self->bl_vel[i][j].y;
				}
			    }
			}
		    }

#ifdef _MESSG
		  for(itrn = 0; itrn < ntransport; itrn++)
		    comm_update_VECT2D(bedload_vel[itrn]);
		  for(itrn = 0; itrn < ntransport; itrn++)
		    comm_update_double(bedload_thick[itrn], ONE);
		  for(k = 0; k < nsed; k++)
		    comm_update_double(rouse_coef[k], ONE);
#endif

		  for(i = 0; i < nnode; i++)
		    {
		      for(k = 0; k < nsand; k++)
			{
			  j = sand_to_transport[k];
			  if(node_flags[i] > NORMAL)
			    {
			      istr = node_flags[i];
    /***SET THE EQUILIBRIUM CONCENTRATION BOUNDARY CONDITION***/
			      if(str_values[istr].trans[j].bc_flag == BCT_CEQ)
				{
				  concentration[j][i] =
				    active_layer_distribution[k][i] / rouse_coef[k][i] *
				    self->source_coef[i][k] / con[j].property[0];
				  bedload_con[j][i] =
				    active_layer_distribution[k][i] * self->bedload_eq_con[i][k];
				}
			    }
			}
		    }
		  for(i = 0; i < nnode; i++)
		    {
		      self->sum_sed_decay[i] = 0.0;
		    }
		  for(itrn = 0; itrn < ntransport; itrn++)
		    {
		      if(((con[itrn].type == SND || con[itrn].type == CLA) &&
			  itrn < (ntransport - 1)))
			{
			  for(i = 0; i < nnode; i++)
			    {
			      if(con[itrn].type == SND)
				{
				  k = transport_to_sand[itrn];
				  self->sum_sed_decay[i] +=
				    sand_settling_velocity[k] * rouse_coef[k][i] *
				    concentration[itrn][i];
				}
			      if(con[itrn].type == CLA)
				{
				  j = transport_to_clay[itrn];
				  sed_number = j + nsand + nsilt;
				  self->sum_sed_decay[i] +=
				    rouse_coef[sed_number][i] * self->prob_dep[i][sed_number] *
				    clay_settling_velocity[j] * concentration[itrn][i];
				}
			    }
			}
		    }

		  for(itrn = 0; itrn < ntransport; itrn++)
		    {
		      if(((con[itrn].type == SND || con[itrn].type == CLA) &&
			  itrn < (ntransport - 1)))
			{

			  for(i = 0; i < nnode; i++)
			    {
			      if(self->bed_character[i] == 1)
				{
				  if(con[itrn].type == SND)
				    {
				      k = transport_to_sand[itrn];
				      sed_number = k;
				      nodal_source[i] =
					-sand_settling_velocity[k] *
					active_layer_distribution[k][i] *
					self->source_coef[i][k] / con[itrn].property[0];
				      self->sum_nodal_source[i] +=
					-sand_settling_velocity[k] *
					self->old_sed_active_layer_distribution[k][i] *
					self->source_coef[i][k] / con[itrn].property[0];
				      nodal_decay_coefficient[i] =
					sand_settling_velocity[k] * rouse_coef[k][i];
				    }

				  /* now that the source and decay have been determined for the sand, get the fraction for the clays */
				  if(con[itrn].type == CLA)
				    {
				      j = transport_to_clay[itrn];
				      sed_number = j + nsand + nsilt;
				      nodal_source[i] =
					self->sum_nodal_source[i] *
					active_layer_distribution[sed_number][i] /
					(MAX(self->sum_snd_ald[i], SMALL));
				      nodal_decay_coefficient[i] =
					rouse_coef[sed_number][i] *
					self->prob_dep[i][sed_number] * clay_settling_velocity[j];
				      /* if clay exists, then largest clay is solid boundary */
				      if(sed_number == nsed - 1)
					{
					  nodal_source[i] = 0.0;
					  nodal_decay_coefficient[i] = 0.0;
					}
				    }
				}
			      if(self->bed_character[i] == 2)
				{
				  if(con[itrn].type == SND)
				    {
				      if(transport_to_sand[itrn] == 0)
					{
					  for(k = 0; k < nclay - 1; k++)
					    {
					      sed_number = nsand + nsilt + k;
					      self->sum_nodal_source[i] +=
						-self->old_sed_active_layer_distribution
						[sed_number][i] *
						self->source_coef[i][sed_number] /
						con[itrn].property[0];
					    }
					}
				    }
				  if(con[itrn].type == CLA)
				    {
				      k = transport_to_clay[itrn];
				      sed_number = nsand + nsilt + k;
				      nodal_source[i] =
					-active_layer_distribution[sed_number][i] *
					self->source_coef[i][sed_number] / con[itrn].property[0];
				      nodal_decay_coefficient[i] =
					rouse_coef[sed_number][i] *
					self->prob_dep[i][sed_number] * clay_settling_velocity[k];
				      /* if clay exists, then largest clay is solid boundary */
				      if(sed_number == nsed - 1)
					{
					  nodal_source[i] = 0.0;
					  nodal_decay_coefficient[i] = 0.0;
					}
				    }
				  if(con[itrn].type == SND)
				    {
				      j = transport_to_sand[itrn];
				      sed_number = j;
				      if(self->source_coef[i][j] > SMALL)
					{
					  nodal_source[i] =
					    -self->source_coef[i][nsand] *
					    active_layer_distribution[sed_number][i] /
					    con[itrn].property[0];
					}
				      else
					{
					  nodal_source[i] = 0.0;
					}
				      nodal_decay_coefficient[i] =
					rouse_coef[sed_number][i] *
					sand_settling_velocity[j];
				    }
				}
			      if(active_layer_erosion_rate_exponent[i] < 0.5 && nsfssi > 1)
				{
				  if(con[itrn].type == SND)
				    {
				      j = transport_to_sand[itrn];
				      sed_number = j;
				      nodal_source[i] =
					-active_layer_distribution[sed_number][i] *
					self->source_coef[i][sed_number] / con[itrn].property[0];
				    }
				  if(con[itrn].type == CLA)
				    {
				      k = transport_to_clay[itrn];
				      sed_number = nsand + nsilt + k;
				      nodal_source[i] =
					-active_layer_distribution[sed_number][i] *
					self->source_coef[i][sed_number] / con[itrn].property[0];
				    }
				}
			      if(ol_head[i] < 0.0)
				{
				  nodal_source[i] = 0.0;
				  nodal_decay_coefficient[i] = 0.0;
				}
			    }
			  esrh = 0.02;
			  for(i = 0; i < nnode; i++)
			    {
			      sed_vorticity_velocity_components(i, &tvel);
			      sed_transport_correction_factors(esrh, ol_head[i], ol_vel[i],
							       tvel,
							       rouse_coef[sed_number][i],
							       &mfcf[sed_number][i],
							       &vcf[sed_number][i]);
			    }
			  if(fe_newton
			     (self->bc_mask, self->matrix, self->diagonal, self->scale_vect, self->residual, self->sol, nnode,
			      my_nnode, nsys, nsys_sq, itrn, fe_sw2trns_init,
			      fe_sw2trns_update, fe_sw2trns_resid, fe_sw2trns_load,
			      fe_sw2trns_inc) == NO)
			    return 1;

			}

		      for(i = 0; i < nnode; i++)
			{

			  if(con[itrn].type == SND)
			    {
			      k = transport_to_sand[itrn];
			      sed_number = k;
			      blvelmag = 0.0;
			      blvelmag =
				pow(bedload_vel[itrn][i].x,
				    2.0) + pow(bedload_vel[itrn][i].y, 2.0);
			      if(blvelmag > 0.0)
				blvelmag = sqrt(blvelmag);
			      blvelmag = MAX(blvelmag, bedload_thick[itrn][i]);
			      nodal_source[i] =
				-active_layer_distribution[sed_number][i] *
				self->bedload_eq_con[i][sed_number] * blvelmag *
				bedload_thick[itrn][i] * self->bl_slc[i][sed_number];
			      nodal_decay_coefficient[i] =
				blvelmag * bedload_thick[itrn][i] * self->bl_slc[i][sed_number];
			    }

			  if(con[itrn].type == CLA)
			    {
			      if(nsand > 0)
				{
				  k = transport_to_clay[itrn];
				  sed_number = nsand + nsilt + k;
				  blvelmag = 0.0;
				  blvelmag =
				    pow(bedload_vel[itrn][i].x,
					2.0) + pow(bedload_vel[itrn][i].y, 2.0);
				  if(blvelmag > 0.0)
				    blvelmag = sqrt(blvelmag);
				  blvelmag = MAX(blvelmag, bedload_thick[itrn][i]);
				  nodal_source[i] =
				    -active_layer_distribution[sed_number][i] *
				    self->bedload_eq_con[i][sed_number] * blvelmag *
				    bedload_thick[itrn][i] * self->bl_slc[i][0];
				  nodal_decay_coefficient[i] =
				    blvelmag * bedload_thick[itrn][i] * self->bl_slc[i][0];
				}
			      else
				{
				  nodal_source[i] = 0.0;
				  nodal_decay_coefficient[i] = 0.0;
				}
			    }

			}

		      if(fe_newton
			 (self->bc_mask, self->matrix, self->diagonal, self->scale_vect, self->residual, self->sol, nnode,
			  my_nnode, nsys, nsys_sq, itrn, fe_sw2trns_bl_init,
			  fe_sw2trns_bl_update, fe_sw2trns_bl_resid, fe_sw2trns_bl_load,
			  fe_sw2trns_bl_inc) == NO)
			{
			  return 1;
			}
		    }

#ifdef _MESSG

		  for(itrn = 0; itrn < ntransport; itrn++)
		    comm_update_double(concentration[itrn], ONE);
		  for(itrn = 0; itrn < ntransport; itrn++)
		    comm_update_double(bedload_con[itrn], ONE);
		  for(k = 0; k < nsed; k++)
		    comm_update_double(self->old_it_concentration[k], ONE);

#endif

		  /*calculate bedload vector sum for each node */

		  for(i = 0; i < nnode; i++)
		    {
		      bedload_vector[i].x = 0.0;
		      bedload_vector[i].y = 0.0;
		      for(itrn = 0; itrn < ntransport; itrn++)
			{
			  bedload_vector[i].x +=
			    bedload_vel[itrn][i].x * bedload_con[itrn][i] *
			    bedload_thick[itrn][i] * density;
			  bedload_vector[i].y +=
			    bedload_vel[itrn][i].y * bedload_con[itrn][i] *
			    bedload_thick[itrn][i] * density;
			}
		    }

#ifdef _MESSG
		  comm_update_VECT2D(bedload_vector);
#endif

		  for(i = 0; i < nnode; i++)
		    {
		      for(k = 0; k < nsed; k++)
			{
			  self->old_it_ald[k][i] = active_layer_distribution[k][i];
			}
		      /* calculate the nodal bed elevation */
		      self->elev_disp[i] = bed_elev[i] + self->old_sed_displacement[i];
		    }
		  /* calculate the nodal bed gradient to be used later, this is a VECT2D[nnode] */
		  vt_2d_grad_nodal_avg(self->elev_disp, self->nodal_grad);

      /*** this update is in ***/
#ifdef _MESSG
		  comm_update_VECT2D(self->nodal_grad);
#endif

		  sed_source_decay(self->source_coef, rouse_coef, self->prob_dep, self->critical_shear_vel,
				   self->nodal_grad, self->old_sed_displacement, self->ustar,
				   active_layer_critical_erosion_shear,
				   active_layer_erosion_rate_constant,
				   active_layer_erosion_rate_exponent, sed_iteration);

		  sed_bedload_source_decay(self->bl_vel, self->bedload_eq_con, self->bl_thick, self->bl_slc, self->ustar,
					   self->bed_character, self->critical_shear_vel, self->nodal_grad,
					   self->old_sed_active_layer_distribution,
					   self->old_sed_displacement,
					   active_layer_critical_erosion_shear,
					   sed_iteration);

		  for(i = 0; i < nnode; i++)
		    {
		      for(itrn = 0; itrn < ntransport; itrn++)
			{
			  if(con[itrn].type == SND || con[itrn].type == CLA)
			    {
			      if(con[itrn].type == SND)
				{
				  j = transport_to_sand[itrn];
				  sed_number = j;
				  bedload_thick[itrn][i] = self->bl_thick[i][j];
				  bedload_vel[itrn][i].x = self->bl_vel[i][j].x;
				  bedload_vel[itrn][i].y = self->bl_vel[i][j].y;
				}
			      if(con[itrn].type == CLA)
				{
				  j = transport_to_clay[itrn];
				  sed_number = j + nsand + nsilt;
				  bedload_thick[itrn][i] = self->bl_thick[i][j];
				  bedload_vel[itrn][i].x = self->bl_vel[i][j].x;
				  bedload_vel[itrn][i].y = self->bl_vel[i][j].y;
				}
			    }
			}
		    }

#ifdef _MESSG
		  for(itrn = 0; itrn < ntransport; itrn++)
		    comm_update_VECT2D(bedload_vel[itrn]);
		  for(itrn = 0; itrn < ntransport; itrn++)
		    comm_update_double(bedload_thick[itrn], ONE);
		  for(k = 0; k < nsed; k++)
		    comm_update_double(rouse_coef[k], ONE);
#endif

		  fe_bedchange_activechange(bedload_vel, self->bedload_eq_con, self->bl_slc,
					    self->source_coef, rouse_coef, self->prob_dep,
					    self->bed_character, self->old_sed_displacement, self->ustar,
					    self->old_sed_active_stratum_ceiling,
					    self->old_sed_active_layer_distribution,
					    self->critical_shear_vel,
					    self->old_sed_active_stratum_distribution,
					    active_layer_porosity,
					    active_layer_critical_erosion_shear,
					    active_layer_erosion_rate_constant,
					    active_layer_erosion_rate_exponent,
					    self->old_sed_active_layer_porosity,
					    self->old_sed_active_layer_critical_erosion_shear,
					    self->old_sed_active_layer_erosion_rate_constant,
					    self->old_sed_active_layer_erosion_rate_exponent,
					    active_stratum_porosity,
					    active_stratum_critical_erosion_shear,
					    active_stratum_erosion_rate_constant,
					    active_stratum_erosion_rate_exponent,
					    self->old_it_rhs, self->limiter, self->limiter2, sed_iteration);

#ifdef _MESSG
		  comm_update_double(displacement, ONE);
		  comm_update_double(active_layer_porosity, ONE);
		  comm_update_double(active_layer_critical_erosion_shear, ONE);
		  comm_update_double(active_layer_erosion_rate_constant, ONE);
		  comm_update_double(active_layer_erosion_rate_exponent, ONE);
		  comm_update_double(active_stratum_porosity, ONE);
		  comm_update_double(active_stratum_critical_erosion_shear, ONE);
		  comm_update_double(active_stratum_erosion_rate_constant, ONE);
		  comm_update_double(active_stratum_erosion_rate_exponent, ONE);
		  for(k = 0; k < nsed; k++)
		    {
		      comm_update_double(active_layer_distribution[k], ONE);
		      comm_update_double(active_stratum_distribution[k], ONE);
		    }
#endif

		  for(i = 0; i < nnode; i++)
		    {
		      ald_diff = 0.0;
		      for(k = 0; k < nsed; k++)
			{
			  ald_diff =
			    MAX(fabs(active_layer_distribution[k][i] - self->old_it_ald[k][i]),
				ald_diff);
			}
		      if(ald_diff >= ald_check)
			{
			  ald_check = ald_diff;
			  max_node = orig_nd_number[i];
			}
		    }

		  ald_max = messg_dmax(ald_check);

		  printf("self->ald_check = %15.6e  at node %d\n", ald_check, max_node);
		  printf("self->ald_max_ = %15.6e, sed_iteration=%d\n", ald_max, sed_iteration);

		}
	      while(sed_iteration < 10);

#ifdef _MESSG
	      comm_update_double(displacement, ONE);
	      comm_update_double(active_layer_porosity, ONE);
	      comm_update_double(active_layer_critical_erosion_shear, ONE);
	      comm_update_double(active_layer_erosion_rate_constant, ONE);
	      comm_update_double(active_layer_erosion_rate_exponent, ONE);
	      comm_update_double(active_stratum_porosity, ONE);
	      comm_update_double(active_stratum_critical_erosion_shear, ONE);
	      comm_update_double(active_stratum_erosion_rate_constant, ONE);
	      comm_update_double(active_stratum_erosion_rate_exponent, ONE);
	      for(k = 0; k < nsed; k++)
		{
		  comm_update_double(active_layer_distribution[k], ONE);
		  comm_update_double(active_stratum_distribution[k], ONE);
		}
#endif

      /*** EXIT WITH A NO IF CONVERGENCE NOT REACHED IN 50 SED_ITERATIONS***/
	      if(ald_max > sed_it_tol)
		{
		  if(myid <= 0)
		    printf("Bed convergence FAILURE!\n");
		  tc_scale();
		  fe_sw2_init(0, self->bc_mask);
		  for(k = 0; k < ntransport; k++)
		    {
		      fe_sw2bed_init(k, self->bc_mask);
		      for(i = 0; i < nnode; i++)
			{
			  displacement[i] = self->old_sed_displacement[i];
			  for(k = 0; k < nsed; k++)
			    {
			      active_layer_distribution[k][i] =
				self->old_sed_active_layer_distribution[k][i];
			      active_stratum_distribution[k][i] =
				self->old_sed_active_stratum_distribution[k][i];
			    }
			}
		    }
		  return 1;
		}

	      /* constituent reactions using Runge-Kutta Technique */
	      /* BTR, 7/04 */
	      if(ntransport > 0)
		fe_rxn_rk5(concentration, dt, self->bc_mask);

	      /* This is used when there is sedimentation that changes the bed elevation */
	      /* if I pass the border node values for the active_stratum_ceiling to the other processors */
	      /* then each of them should be able to calculate the bed_thickness and bed distribution on */
	      /* their own.  this is why the sed_bed_layer_update loops on nnode rather than my_nnode.   */

#ifdef _MESSG
	      comm_update_double(active_layer_thickness, ONE);
#endif
	      for(i = 0; i < nnode; i++)
		{
		  active_stratum_ceiling[i] =
		    l_reference * displacement[i] - active_layer_thickness[i];
		}
#ifdef _MESSG
	      comm_update_double(active_stratum_ceiling, ONE);
#endif

	      for(i = 0; i < nnode; i++)
		{
		  bed_thick_old = 0.0;
		  bed_thick_new = 0.0;
		  for(k = 0; k < nsed; k++)
		    {
		      self->ald[k] = active_layer_distribution[k][i];
		      self->asd[k] = active_stratum_distribution[k][i];
		      for(j = 0; j < number_bed_layers; j++)
			{
			  self->bld[j][k] = bed_layer_distribution[j][k][i];
			}
		    }
		  for(k = 0; k < 5; k++)
		    {
		      for(j = 0; j < ncti; j++)
			{
			  self->cbp[j][k] = consolidation_bed_properties[k][j][i];
			}
		    }
		  for(j = 0; j < number_bed_layers; j++)
		    {
		      self->blt[j] = bed_layer_thickness[j][i];
		      blsg = 0.0;
		      for(k = 0; k < nsed; k++)
			{
			  blsg += self->bld[j][k] * sed_specific_gravity[k];
			}
		      self->bl_por[j] =
			(bed_layer_bulk_density[j][i] -
			 density * blsg) / (density * (1. - blsg));
		      self->bl_ces[j] = bed_layer_critical_erosion_shear[j][i];
		      self->bl_erc[j] = bed_layer_erosion_rate_constant[j][i];
		      self->bl_ere[j] = bed_layer_erosion_rate_exponent[j][i];
		      bed_thick_old += bed_layer_thickness[j][i];
		    }

		  /* self->old_sed_displacement = displacement - delta_displ */
		  sed_bed_layer_update(self->old_sed_displacement[i], displacement[i],
				       active_stratum_ceiling[i],
				       self->old_sed_active_stratum_ceiling[i], self->bed_character[i],
				       number_bed_layers, nsed, self->bld, self->blt, self->ald, self->asd,
				       active_layer_porosity[i],
				       active_layer_critical_erosion_shear[i],
				       active_layer_erosion_rate_constant[i],
				       active_layer_erosion_rate_exponent[i],
				       active_stratum_porosity[i],
				       active_stratum_critical_erosion_shear[i],
				       active_stratum_erosion_rate_constant[i],
				       active_stratum_erosion_rate_exponent[i], self->bl_por,
				       self->bl_ces, self->bl_erc, self->bl_ere);

		  for(j = 0; j < number_bed_layers; j++)
		    {
		      self->bl_sg[j] = 0.0;
		      for(k = 0; k < nsed; k++)
			{
			  self->bl_sg[j] += self->bld[j][k] * sed_specific_gravity[k];
			}
		      self->bl_bd[j] = density * (self->bl_por[j] + (1. - self->bl_por[j]) * self->bl_sg[j]);
		    }

		  sed_consolidate(number_bed_layers, ncti, dt_sed, self->cbp, self->blt, self->bl_bd, self->bl_sg,
				  self->bl_por, self->bl_ces, self->bl_erc, self->bl_ere);

		  for(k = 0; k < nsed; k++)
		    {
		      active_layer_distribution[k][i] = self->ald[k];
		      active_stratum_distribution[k][i] = self->asd[k];
		      for(j = 0; j < number_bed_layers; j++)
			{
			  bed_layer_distribution[j][k][i] = self->bld[j][k];
			}
		    }
		  for(j = 0; j < number_bed_layers; j++)
		    {
		      bed_layer_thickness[j][i] = self->blt[j];
		      bed_layer_bulk_density[j][i] = self->bl_bd[j];
		      bed_layer_critical_erosion_shear[j][i] = self->bl_ces[j];
		      bed_layer_erosion_rate_constant[j][i] = self->bl_erc[j];
		      bed_layer_erosion_rate_exponent[j][i] = self->bl_ere[j];
		      bed_thick_new += self->blt[j];
		    }
		  /* NOW UPDATE BED ONCE MORE TO CORRECT FOR CONSOLIDATION AND SMALL EROSIONS OF THE SOLID BOUNDARY CAUSED BY BEDLOAD */
		  displacement[i] = self->old_sed_displacement[i] + bed_thick_new - bed_thick_old;
		  active_stratum_ceiling[i] =
		    l_reference * displacement[i] - active_layer_thickness[i];
		}
#ifdef _MESSG
	      comm_update_double(displacement, ONE);
	      comm_update_double(active_stratum_ceiling, ONE);
#endif
	      /* copies the transport dependent solutions to the old solutions since the 
	         time step has made it to the end and is therefore to be kept */
	      for(i = 0; i < nnode; i++)
		{
		  /*Calculate flowrate from hydro results...then redefine head and velocity after displacement */
		  flowrate_x = ol_head[i] * ol_vel[i].x;
		  flowrate_y = ol_head[i] * ol_vel[i].y;
		  self->old_sed_active_stratum_ceiling[i] = active_stratum_ceiling[i];

		  if(STEADY_STATE_HYDRO)
		    qus_delta_displacement[i] += displacement[i] - self->old_sed_displacement[i];

		  old_sed_ol_head[i] = ol_head[i];
		  self->old_sed_displacement[i] = displacement[i];

		  /* update the velocity to maintain constant flowrate */
		  /* want to maintain the direction of the flow */
		  /* WAS DRYING_LOWER_LIMIT INSTEAD OF 0.01 */
		  /* if(ol_head[i] < drying_lower_limit)
		     {
		     ol_vel[i].x = (MAX(ol_head[i],0.0) / drying_lower_limit) * flowrate_x / drying_lower_limit;
		     ol_vel[i].y = (MAX(ol_head[i],0.0) / drying_lower_limit) * flowrate_y / drying_lower_limit;
		     }
		     else
		     {
		     ol_vel[i].x = flowrate_x / ol_head[i];
		     ol_vel[i].y = flowrate_y / ol_head[i];
		     } */

		  /* update the remaining sediment variables for this sed timestep */
		  for(k = 0; k < nsed; k++)
		    {
		      if(k < nsand)
			{
			  j = sand_to_transport[k];
			}
		      else
			{
			  j = clay_to_transport[k - nsand];
			}
		      self->old_sed_active_layer_distribution[k][i] =
			active_layer_distribution[k][i];
		      self->old_sed_active_stratum_distribution[k][i] =
			active_stratum_distribution[k][i];

		      older_concentration[j][i] = old_sed_concentration[k][i];
		      old_sed_concentration[k][i] = concentration[j][i];
		      old_sed_bedload_con[k][i] = bedload_con[j][i];
		      old_sed_bedload_thick[k][i] = bedload_thick[j][i];
		    }
		  self->old_sed_active_layer_porosity[i] = active_layer_porosity[i];
		  self->old_sed_active_layer_critical_erosion_shear[i] =
		    active_layer_critical_erosion_shear[i];
		  self->old_sed_active_layer_erosion_rate_constant[i] =
		    active_layer_erosion_rate_constant[i];
		  self->old_sed_active_layer_erosion_rate_exponent[i] =
		    active_layer_erosion_rate_exponent[i];
		}

#ifdef _MESSG
      /*** this update is in ***/
	      comm_update_double(ol_head, ONE);
	      comm_update_double(old_sed_ol_head, ONE);
	      comm_update_double(self->old_sed_active_layer_porosity, ONE);
	      comm_update_double(self->old_sed_active_layer_critical_erosion_shear, ONE);
	      comm_update_double(self->old_sed_active_layer_erosion_rate_constant, ONE);
	      comm_update_double(self->old_sed_active_layer_erosion_rate_exponent, ONE);
	      comm_update_VECT2D(ol_vel);
	      for(k = 0; k < nsed; k++)
		comm_update_double(old_sed_concentration[k], ONE);
	      for(k = 0; k < nsed; k++)
		comm_update_double(old_sed_bedload_con[k], ONE);
	      for(k = 0; k < nsed; k++)
		comm_update_double(old_sed_bedload_thick[k], ONE);

	      sed_time_end = MPI_Wtime();
	      sed_time_total = sed_time_end - sed_time_start;
	      if(myid == 0)
		printf("sed time = %e\n", sed_time_total);
#endif
	      if(STEADY_STATE_HYDRO)
		{
		  for(i = 0; i < nnode; i++)
		    {
		      qus_max_change =
			MAX(qus_depth_factor * (ol_head[i] - qus_delta_displacement[i]),
			    drying_lower_limit);
		      if(qus_delta_displacement[i] > qus_max_change)
			NEW_STEADY_HYDRO = YES;
		    }
		  NEW_STEADY_HYDRO = messg_imax(NEW_STEADY_HYDRO);

		  if(NEW_STEADY_HYDRO == YES)
		    {
		      for(i = 0; i < nnode; i++)
			qus_delta_displacement[i] = 0.0;
		    }
		}

	    }
	  while(t_sed < t_prev + dt);

	}

      /* copies the solutions to the old solutions since the 
         time step has made it to the end and is therefore to be kept */
      for(i = 0; i < nnode; i++)
	{
	  /* dss >>> */
	  older_ol_head[i] = old_ol_head[i];
	  older_ol_vel[i].x = old_ol_vel[i].x;
	  older_ol_vel[i].y = old_ol_vel[i].y;
	  /* dss <<< */
	  old_ol_head[i] = ol_head[i];
	  old_ol_vel[i].x = ol_vel[i].x;
	  old_ol_vel[i].y = ol_vel[i].y;
	  for(itrn = 0; itrn < ntransport; itrn++)
	    {
	      older_concentration[itrn][i] = old_concentration[itrn][i];
	      old_concentration[itrn][i] = concentration[itrn][i];
	    }
	}
      if(SEDIMENT == YES)
	{
	  for(i = 0; i < nnode; i++)
	    {
	      /* update ol_head such that the water doesn't rise when the displacement changes JNT 4-1-04 */
	      delta_displ[i] = displacement[i] - old_displacement[i];
	      /* ol_head[i] -= delta_displ[i]; */
	      old_ol_head[i] = ol_head[i];
	      old_active_stratum_ceiling[i] = active_stratum_ceiling[i];
	      old_active_layer_porosity[i] = active_layer_porosity[i];
	      old_active_layer_critical_erosion_shear[i] =
		active_layer_critical_erosion_shear[i];
	      old_active_layer_erosion_rate_constant[i] =
		active_layer_erosion_rate_constant[i];
	      old_active_layer_erosion_rate_exponent[i] =
		active_layer_erosion_rate_exponent[i];
	      old_active_stratum_porosity[i] = active_stratum_porosity[i];
	      old_active_stratum_critical_erosion_shear[i] =
		active_stratum_critical_erosion_shear[i];
	      old_active_stratum_erosion_rate_constant[i] =
		active_stratum_erosion_rate_constant[i];
	      old_active_stratum_erosion_rate_exponent[i] =
		active_stratum_erosion_rate_exponent[i];
	      old_displacement[i] = displacement[i];
	      for(k = 0; k < nsed; k++)
		{
		  old_active_layer_distribution[k][i] = active_layer_distribution[k][i];
		  old_active_stratum_distribution[k][i] = active_stratum_distribution[k][i];
		  for(n = 0; n < number_bed_layers; n++)
		    old_bed_layer_distribution[n][k][i] = bed_layer_distribution[n][k][i];
		}
	      for(itrn = 0; itrn < ntransport; itrn++)
		{
		  old_bedload_con[itrn][i] = bedload_con[itrn][i];
		  old_bedload_thick[itrn][i] = bedload_thick[itrn][i];
		}
	      for(n = 0; n < number_bed_layers; n++)
		{
		  old_bed_layer_thickness[n][i] = bed_layer_thickness[n][i];
		  old_bed_layer_bulk_density[n][i] = bed_layer_bulk_density[n][i];
		  old_bed_layer_critical_erosion_shear[n][i] =
		    bed_layer_critical_erosion_shear[n][i];
		  old_bed_layer_erosion_rate_constant[n][i] =
		    bed_layer_erosion_rate_constant[n][i];
		  old_bed_layer_erosion_rate_exponent[n][i] =
		    bed_layer_erosion_rate_exponent[n][i];
		}
	    }
	}

      /*  END OF SW2 FLOW SECTION */

    }
  if(SW3_FLOW)
    {
      /* solves the 3D shallow water flow problem - returns from fe_main with a NO if 
         the nonlinear iterations failed */
      if(myid <= 0)
	{
	  printf
	    ("\n time step: 2D Shallow Water Flow. (time = %f, %3.0f %% ), dt = %15.6e\n",
	     t_prev, 100. * t_prev / tfinal, dt);
	}

      /* get the density based on salinity and/or temperature */
      switch (baroclinic_code)
	{
	case 1:
	  tl_density_calculator_metric(NULL, 1., concentration[transport.sal],
				       con[transport.sal].property[0], nnode, dnst_3d);
	  break;
	case 10:
	  tl_density_calculator_metric(concentration[transport.tmp],
				       con[transport.tmp].property[0], NULL, 1., nnode,
				       dnst_3d);
	  break;
	case 11:
	  tl_density_calculator_metric(concentration[transport.tmp],
				       con[transport.tmp].property[0],
				       concentration[transport.sal],
				       con[transport.sal].property[0], nnode, dnst_3d);
	  break;
	}
      if(baroclinic_code > 0)
	{
	  comm_update_double(dnst_3d, 1);
	  integ_density(dnst_3d, dnst_2d);
	  comm_update_double(dnst_2d, 1);
	  for(i = 0; i < nnode; i++)
	    {
	      dnst_2d[i] = dnst_2d[i] / density - 1.;
	    }
	}
      get_bottom_vel(bed_vel);
      comm_update_VECT2D(bed_vel);
      get_bottom_elev(bed_elev);
      comm_update_double(bed_elev, 1);
      integ_velocity2(vel, uu_integral, vv_integral, uv_integral);
      comm_update_double(uu_integral, 1);
      comm_update_double(vv_integral, 1);
      comm_update_double(uv_integral, 1);

      for(i = 0; i < nnode; i++)
	{
	  if(ol_head[i] > 0.)
	    {
	      uu_integral[i] -= ol_head[i] * ol_vel[i].x * ol_vel[i].x;
	      uv_integral[i] -= ol_head[i] * ol_vel[i].x * ol_vel[i].y;
	      vv_integral[i] -= ol_head[i] * ol_vel[i].y * ol_vel[i].y;
	    }
	  else
	    {
	      uu_integral[i] = 0.;
	      uv_integral[i] = 0.;
	      vv_integral[i] = 0.;
	    }
	}
      nsys = 3;
      nsys_sq = 9;

      if(CONVEYANCE)
	fe_sw2_conveyance();

      if(fe_newton
	 (self->bc_mask, self->matrix, self->diagonal, self->scale_vect, self->residual, self->sol, nnode, my_nnode, nsys,
	  nsys_sq, 0, fe_sw2_init, fe_sw2_update, fe_sw2_resid, fe_sw2_load,
	  fe_sw2_inc) == NO)
	return 1;

      /* calculate the vertical grid location */
      if(myid <= 0)
	{
	  printf(" \n");
	  printf("GRID MOVEMENT \n");
	}
      nsys = 1;
      nsys_sq = 1;
      l_reference = 1.;		/* this is put in temporarily since the shallow water calculations */
      u_reference = 1.;		/* are not non-dimensionalized, but the Navier Stokes solver is.   */
      /* And the grid movement was written for the Navier Stokes Solver. */

      fe_sw3_head_disp();
      if(fe_newton
	 (self->bc_mask, self->matrix, self->diagonal, self->scale_vect, self->residual, self->sol, nnode, my_nnode, nsys,
	  nsys_sq, 0, fe_mg_init, fe_mg_update, fe_mg_resid, fe_mg_load, fe_mg_inc) == NO)
	return 1;

      /* distribute the depth average velocity over the 3D nodes */
      if(myid <= 0)
	{
	  printf("  \n");
	  printf("DISTRIBUTE OL_VEL CALCULATION \n");
	}
      refresh_pc = YES;

      /* calculate the hydrostatic pressure */

      if(myid <= 0)
	{
	  printf("  \n");
	  printf("PRESSURE CALCULATION \n");
	}

      nsys = 1;
      nsys_sq = 1;

      get_surface_double(ol_head);
      comm_update_double(ol_head, 1);
      get_surface_vect2d(ol_vel);
      comm_update_VECT2D(ol_vel);
      if(baroclinic_code > 0)
	{
	  integ_density_press(prs);
	}
      else
	{
	  integ_press(prs);
	}
      comm_update_double(prs, 1);

      /* calculate the horizontal velocity (momentum)and vertical velocity (continuity) */
      if(myid <= 0)
	{
	  printf(" \n");
	  printf("VEL_HAT CALCULATION \n");
	}
      refresh_pc = YES;

      if(fe_newton
	 (self->bc_mask, self->matrix, self->diagonal, self->scale_vect, self->residual, self->sol, nnode, my_nnode, 2, 4, 0,
	  fe_sw3_frhvel_init, fe_sw3_frhvel_update, fe_sw3_frhvel_resid, fe_sw3_frhvel_load,
	  fe_sw3_frhvel_inc) == NO)
	return 1;
      refresh_pc = YES;
      if(fe_newton
	 (self->bc_mask, self->matrix, self->diagonal, self->scale_vect, self->residual, self->sol, nnode, my_nnode, 1, 1, 0,
	  fe_sw3_frvvel_init, fe_sw3_frvvel_update, fe_sw3_frvvel_resid, fe_sw3_frvvel_load,
	  fe_sw3_frvvel_inc) == NO)
	return 1;

      for(i = 0; i < nnode; i++)
	{
	  vel_cor[i].x = vel_hat[i].x;
	  vel_cor[i].y = vel_hat[i].y;
	  vel_cor[i].z = vel_hat[i].z;
	}
      if(fe_newton
	 (self->bc_mask, self->matrix, self->diagonal, self->scale_vect, self->residual, self->sol, nnode, my_nnode, 1, 1, 0,
	  fe_sw3_wcbc_init, fe_sw3_wcbc_update, fe_sw3_wcbc_resid, fe_sw3_wcbc_load,
	  fe_sw3_wcbc_inc) == NO)
	return 1;

      /* calculate the correction velocity */
      if(myid <= 0)
	{
	  printf(" \n");
	  printf("U CORRECTION CALCULATION \n");
	}

      nsys = 1;
      nsys_sq = 1;
      if(myid <= 0)
	{
	  printf(" \n");
	  printf("V CORRECTION CALCULATION \n");
	}
      if(myid <= 0)
	{
	  printf(" \n");
	  printf("W CORRECTION CALCULATION \n");
	}
      if(myid <= 0)
	{
	  printf(" \n");
	  printf("Calculate the potential to conserve mass\n");
	}

      if(fe_newton
	 (self->bc_mask, self->matrix, self->diagonal, self->scale_vect, self->residual, self->sol, nnode, my_nnode, 1, 1, 0,
	  fe_sw3_frwspt_init, fe_sw3_frwspt_update, fe_sw3_frwspt_resid, fe_sw3_frwspt_load,
	  fe_sw3_frwspt_inc) == NO)
	return 1;

      if(myid <= 0)
	{
	  printf(" \n");
	  printf("Calculate the final x-component of velocity\n");
	}

      if(fe_newton
	 (self->bc_mask, self->matrix, self->diagonal, self->scale_vect, self->residual, self->sol, nnode, my_nnode, 1, 1, 0,
	  fe_sw3_frwsptgrad_init, fe_sw3_frvel_update, fe_sw3_frdwsptdx_resid,
	  fe_sw3_frwsptgrad_load, fe_sw3_frdwsptdx_inc) == NO)
	return 1;

      if(myid <= 0)
	{
	  printf(" \n");
	  printf("Calculate the final y-component of velocity\n");
	}

      if(fe_newton
	 (self->bc_mask, self->matrix, self->diagonal, self->scale_vect, self->residual, self->sol, nnode, my_nnode, 1, 1, 0,
	  fe_sw3_frwsptgrad_init, fe_sw3_frvel_update, fe_sw3_frdwsptdy_resid,
	  fe_sw3_frwsptgrad_load, fe_sw3_frdwsptdy_inc) == NO)
	return 1;

      if(myid <= 0)
	{
	  printf(" \n");
	  printf("Calculate the final z-component of velocity\n");
	}

      if(fe_newton
	 (self->bc_mask, self->matrix, self->diagonal, self->scale_vect, self->residual, self->sol, nnode, my_nnode, 1, 1, 0,
	  fe_sw3_frwsptgrad_init, fe_sw3_frvel_update, fe_sw3_frdwsptdz_resid,
	  fe_sw3_frwsptgrad_load, fe_sw3_frdwsptdz_inc) == NO)
	return 1;

      /* calculate the transport of constituents */
      if(myid <= 0)
	{
	  printf("\n");
	  printf(" TRANSPORT CONSTITUENTS \n");
	}
      refresh_pc = YES;
      nsys = 1;
      nsys_sq = 1;

      for(itrn = 0; itrn < ntransport; itrn++)
	{
	  if(fe_newton
	     (self->bc_mask, self->matrix, self->diagonal, self->scale_vect, self->residual, self->sol, nnode, my_nnode, nsys,
	      nsys_sq, itrn, fe_nstrns_init, fe_nstrns_update, fe_sw3trns_resid,
	      fe_sw3trns_load, fe_nstrns_inc) == NO)
	    return 1;
	}

      /* copies the solutions to the old solutions since the 
         time step has made it to the end and is therefore to be kept */

      for(i = 0; i < nnode; i++)
	{
	  older_ol_head[i] = old_ol_head[i];
	  older_ol_vel[i].x = old_ol_vel[i].x;
	  older_ol_vel[i].y = old_ol_vel[i].y;
	  old_ol_head[i] = ol_head[i];
	  old_ol_vel[i].x = ol_vel[i].x;
	  old_ol_vel[i].y = ol_vel[i].y;
	  old_vel[i].x = vel[i].x;
	  old_vel[i].y = vel[i].y;
	  old_vel[i].z = vel[i].z;
	  old_prs[i] = prs[i];
	  old_displacement[i] = displacement[i];
	  old_grid_speed[i] = grid_speed[i];
	  old_ws_pot[i] = ws_pot[i];
	  for(itrn = 0; itrn < ntransport; itrn++)
	    old_concentration[itrn][i] = concentration[itrn][i];
	}
    }

  /* if it makes it to here then the calculations were good */
  /* must wait til the end to update in case of failure in second or third
     equation set */
  if(GW_FLOW)
    {
      /* copies the solutions to the old solutions since the 
         time step has made it to the end and is therefore to be kept */
      for(i = 0; i < nnode; i++)
	{
	  old_head[i] = head[i];
	}
      for(i = 0; i < nnode; i++)
	{
	  for(itrn = 0; itrn < ntransport; itrn++)
	    {
	      old_concentration[itrn][i] = concentration[itrn][i];
	      old_sp_concentration[itrn][i] = sp_concentration[itrn][i];
	    }
	}
      for(ie = 0; ie < nelem3d; ie++)
	{
	  for(i = 0; i < NDPRELM; i++)
	    {
	      old_saturation[ie][i] = saturation[ie][i];
	    }
	}
    }
  if(HEAT)
    {
      /* copies the solution to the old solution */
      for(i = 0; i < nnode; i++)
	{
	  old_temper[i] = temper[i];
	}
    }
  if(OL_FLOW)
    {
      /* copies the solution to the old solution */
      for(i = 0; i < nnode; i++)
	{
	  old_ol_head[i] = ol_head[i];
	}
      for(i = 0; i < nnode; i++)
	{
	  for(itrn = 0; itrn < ntransport; itrn++)
	    {
	      old_concentration[itrn][i] = concentration[itrn][i];
	      old_sp_concentration[itrn][i] = sp_concentration[itrn][i];
	    }
	}
    }
  return 0;
 }

/** 
    solution info
**/
int get_ADH_NumberOfComponents() {return max_nsys;}
int get_ADH_NumberOfDegreesOfFreedom() {return ndof;}
double get_ADH_newTimeLevel() { return t_prev+dt;}
double get_ADH_prevTimeLevel() { return t_prev;}
double get_ADH_dt() { return dt;}
double get_ADH_tinitial() { return tinitial;}
void* get_ADH_solutionComponent(int ci)
{
  if(GW_FLOW)
    return (void*)head;
  if(HEAT)
    return (void*)temper;
  if(SW2_FLOW)
    {
      if (ci == 2)
	return (void*)ol_head;
      return (void*) ol_vel;
    }
  if(NS_FLOW)
    {
      if (ci == 3)
	return (void*)prs;
      return (void*) vel;
    }
  /*todo need to fill in other models*/
  assert(0);
  return 0;
}
int get_ADH_solutionComponentDimension(int ci)
{
  if(GW_FLOW)
    return 1;
  if(HEAT)
    return 1;
  if(SW2_FLOW)
    {
      if (ci == 2)
	return 1;
      return 2;
    }
  if(NS_FLOW)
    {
      if (ci == 3)
	return 1;
      return 3;
    }
  /*todo fill in other models*/
  assert(0);
  return -1;
}
int get_ADH_firstSolutionComponentForVectorComponent(int ci)
{
  if(SW2_FLOW)
    {
      /*storage order (u,v,h)*/
      return 0;
    }
  if(NS_FLOW)
    {
      /*storage order (u,v,w,p)*/
      return 0;
    }
  /*todo fill in other models*/
  return 0;
}

int get_ADH_ADAPTION_flag(void)
{
  if(ADAPTION)
    return 1;
  return 0;
}
/** 
    mesh info
**/
int get_ADH_nElements_global()
{
  if (nelem3d > 0)
    return nelem3d;
  if (nelem2d > 0)
    return nelem2d;
  if (nelem1d > 0)
    return nelem1d;
  return -1;
}
int get_ADH_nNodes_global()
{
  /*todo: make sure if number of nodes on local processor mesh is
    nnode and my_nnode is the number owned or if nnode is the global
    number of nodes*/
  return nnode;
}
int get_ADH_nSpace_global()
{
  if (nelem3d > 0)
    return 3;
  if (nelem2d > 0)
    return 2;
  if (nelem1d > 0)
    return 1;
  return -1;
}
int get_ADH_nNodes_element()
{
  if (nelem3d > 0)
    return 4;
  if (nelem2d > 0)
    return 3;
  if (nelem1d > 0)
    return 2;
  return -1;
}

int get_nodeArray(double * nodeArray)
{
  /*
    load adh nodes into a double* array for proteus
    nodeArray is logically nNodes_global x 3
   */
  int nN;
  for (nN = 0; nN < nnode; nN++)
    {
      nodeArray[nN*3+0] = node[nN].x;
      nodeArray[nN*3+1] = node[nN].y;
      nodeArray[nN*3+2] = node[nN].z;

    }
  return 0;
}

int get_elementNodesArray(int* elementNodesArray)
{
  /*
    load adh element-node connectivity table into a int* array for proteus
    elementNodesArray is logically nElements_global x (nNodes_element)
   */
  int eN,nN,nNodes_element;
  if (nelem3d > 0)
    {
      nNodes_element = 4;
      for (eN=0; eN < nelem3d; eN++)
	for (nN=0; nN < nNodes_element; nN++)
	  elementNodesArray[eN*nNodes_element+nN] = elem3d[eN].nodes[nN];
    }
  else if (nelem2d > 0)
    {
      nNodes_element = 3;
      for (eN=0; eN < nelem2d; eN++)
	for (nN=0; nN < nNodes_element; nN++)
	  elementNodesArray[eN*nNodes_element+nN] = elem2d[eN].nodes[nN];

    }
  else
    {
      nNodes_element = 2;
      for (eN=0; eN < nelem1d; eN++)
	for (nN=0; nN < nNodes_element; nN++)
	  elementNodesArray[eN*nNodes_element+nN] = elem1d[eN].nodes[nN];

    }
  return 0;
}

int get_elementMaterialTypes(int* elementMaterialTypes)
{
  /* copy adh element material identifier to proteus
     may be better to go with elem3d[eN].mat instead
  */
  int eN;
  if (nelem3d > 0)
    {
      for (eN=0; eN < nelem3d; eN++)
	elementMaterialTypes[eN] = elem3d_flags[eN];
    } 
  else if (nelem2d > 0)
    {
      for (eN=0; eN < nelem2d; eN++)
	{
	  /*looks like for trunk this should be elem2d_mat*/
	  elementMaterialTypes[eN] = elem2d_mat[eN];
	}
    } 
  else
    {
      assert(0);
    }
  return 0;
}
int get_nodeMaterialTypes(int* nodeMaterialTypes)
{
  /* copy adh node identifier to proteus*/
  int nN;
  for (nN=0; nN < nnode; nN++)
    nodeMaterialTypes[nN] = node_flags[nN];
   
  return 0;
}

int ADH_PETSc_Interface_init(ADH_PETSc_Interface* self, ADH_OneLevelTransport* cadh_transport)
{
  printf("trying to init ADH_PETSC %d %d \n",npes,nnode);
  /* solve_petsc_init */
  int i;
  /* allocate storage for the integer mapping */
  printf("all gather\n");
  self->nnodeLocal =
    (int *)tl_alloc(sizeof(int), npes);
  printf("all gather\n");
  self->globalNodeNumberStart =
    (int *)tl_alloc(sizeof(int), npes);
  printf("all gather\n");
  self->nnodeLocal_old =
    (int *)tl_alloc(sizeof(int), npes);
  printf("all gather\n");
  self->globalNodeNumber =
    (int *)tl_alloc(sizeof(int), nnode);
  printf("all gather\n");
  /* get the number of local (non-ghost) nodes on each proc and generate local (processor) to global mapping */
  MPI_Allgather(&my_nnode, 1, MPI_INT, self->nnodeLocal, 1, MPI_INT,PETSC_COMM_WORLD);
  self->globalNodeNumberStart[0] = 0;
  for(i = 1; i < npes; i++)
    {
      self->nnodeLocal_old[i] = self->nnodeLocal[i];
      self->globalNodeNumberStart[i] =
	self->globalNodeNumberStart[i - 1] + self->nnodeLocal[i - 1];
    }
  for(i = 0; i < nnode; i++)
    self->globalNodeNumber[i] =
      self->globalNodeNumberStart[node_pair[i].sd] + node_pair[i].rnode;
  printf("trying to Create mapping \n");
  ISLocalToGlobalMappingCreate(PETSC_COMM_WORLD, nnode,
			       self->globalNodeNumber,
			       &self->localToGlobalMapping);
  self->mappingHasChanged = FALSE;
  /* solver, matrices, and vectors for each type of physics */
  for (i=0;i<10;i++)
    {
      self->L[i] = NULL;
      self->du[i] = NULL;
      self->r[i] = NULL;
      self->ksp[i] = NULL;
    }
  printf("trying to update storage \n");
  ADH_PETSc_updateStorage(self,cadh_transport);
  return 0;
}

void ADH_PETSc_matrix_resize(ADH_PETSc_Interface* self,
			     ADH_OneLevelTransport* cadh_transport,
			     KSP* ksp, 
			     Mat* matrix,
			     int nsys, 
			     Vec* sol, 
			     Vec* residual)
{
  /* estimate of  the nonzero structure, used for preallocating sparse matrix, need to check this */
  int nNonZeroDiagonalBands = 3 * nsys;
  int nNonZeroOffDiagonalBands = 24 * nsys;
  int allocatedNewMatrix = NO;
  if(self->mappingHasChanged && matrix != NULL)
    {
      printf("mapping has changed = %i", self->mappingHasChanged);
      ADH_PETSc_matrix_destroy(ksp, matrix, sol, residual);
    }
  /* allocate the matrix if necessary */
  if(matrix == NULL)
    {
      allocatedNewMatrix = YES;
#ifdef _MESSG
      VecCreateMPIWithArray(PETSC_COMM_WORLD, my_nnode * nsys, PETSC_DETERMINE, cadh_transport->sol, sol);
      VecCreateMPIWithArray(PETSC_COMM_WORLD, my_nnode * nsys, PETSC_DETERMINE, cadh_transport->residual, residual);
      /*VecCreateMPI(PETSC_COMM_WORLD, my_nnode * nsys, PETSC_DETERMINE, &matrix->petsc_diag);*/
#else
      VecCreateSeqWithArray(my_nnode * nsys, cadh_transport->sol, sol);
      VecCreateSeqWithArray(my_nnode * nsys, cadh_transport->residual, residual);
      /*VecCreateSeq(my_nnode * nsys, &matrix->petsc_diag);*/
#endif
      if(nsys == 1)
	{
#ifdef _MESSG
	  printf("Allocating petsc MatCreateMPIAIJ matrix \n");
	  printf("number of nodes %i \n", my_nnode);
	  printf("number of rows %i \n", my_nnode * nsys);
	  printf("number of columns %i \n", my_nnode * nsys);
	  MatCreateMPIAIJ(PETSC_COMM_WORLD, my_nnode,
			  my_nnode, PETSC_DETERMINE, PETSC_DETERMINE,
			  nNonZeroDiagonalBands, PETSC_NULL,
			  nNonZeroOffDiagonalBands, PETSC_NULL,
			  matrix);
#else
	  printf("Allocating petsc MatCreateSeqAIJ matrix \n");
	  printf("number of nodes %i \n", my_nnode);
	  printf("number of rows %i \n", my_nnode * nsys);
	  printf("number of columns %i \n", my_nnode * nsys);
	  MatCreateSeqAIJ(PETSC_COMM_WORLD, nnode, nnode,
			  nNonZeroDiagonalBands + nNonZeroOffDiagonalBands,
			  PETSC_NULL, matrix);
#endif
	}
      else
	{
#ifdef _MESSG
	  printf("Allocating petsc MatCreateMPIBAIJ matrix \n");
	  printf
	    ("Petsc Matrix block size (number of components per node) %i \n",
	     nsys);
	  printf("number of nodes %i \n", my_nnode);
	  printf("number of rows %i \n", my_nnode * nsys);
	  printf("number of columns %i \n", my_nnode * nsys);
	  MatCreateMPIBAIJ(PETSC_COMM_WORLD,
			   nsys,
			   my_nnode * nsys,
			   my_nnode * nsys,
			   PETSC_DETERMINE, PETSC_DETERMINE,
			   nNonZeroDiagonalBands * nsys, PETSC_NULL,
			   nNonZeroOffDiagonalBands * nsys, PETSC_NULL,
			   matrix);
#else
	  printf("Allocating petsc MatCreateSeqBAIJ  matrix");
	  printf("block size %i", nsys);
	  printf("number of rows %i", my_nnode * nsys);
	  printf("number of columns %i", my_nnode * nsys);
	  MatCreateSeqBAIJ(PETSC_COMM_WORLD, nsys,
			   my_nnode * nsys, my_nnode * nsys,
			   nNonZeroDiagonalBands * nsys +
			   nNonZeroOffDiagonalBands * nsys, PETSC_NULL,
			   matrix);
#endif
	}
      MatSetFromOptions(*matrix);
      MatSetOption(*matrix, MAT_IGNORE_OFF_PROC_ENTRIES,PETSC_TRUE);
      if(nsys == 1)
	MatSetLocalToGlobalMapping(*matrix, self->localToGlobalMapping);
      else
	MatSetLocalToGlobalMappingBlock(*matrix,
					self->localToGlobalMapping);
      /* build krylov solver */
      KSPCreate(PETSC_COMM_WORLD, ksp);
      KSPSetFromOptions(*ksp);
      /* draw some simple line graphs of the residual and newton correction */
      /* PetscViewerDrawOpen(PETSC_COMM_WORLD, PETSC_NULL, residual->name, */
      /* 			  PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE, */
      /* 			  PETSC_DECIDE, &residual->petscViewer); */
      /* PetscViewerDrawOpen(PETSC_COMM_WORLD, PETSC_NULL, sol->name, */
      /* 			  PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE, */
      /* 			  PETSC_DECIDE, &sol->petscViewer); */
    }
}

void ADH_PETSc_matrix_destroy(KSP* ksp, Mat *matrix, Vec* sol, Vec* residual)
{
  if(matrix != NULL)
    {
      printf("Deallocating PETSc matrix\n %p", matrix);
      MatDestroy(*matrix);
      printf("Finished with PETSc matrix \n");
      printf("Deallocating PETSc vectors\n");
      VecDestroy(*sol);
      VecDestroy(*residual);
      printf("Finished with PETSc vectors \n");
/*       PetscViewerDestroy(*sol); */
/*       PetscViewerDestroy(*residual); */
      KSPDestroy(*ksp);/* cek todo seems like this should go first*/
    }
}

int ADH_PETSc_updateStorage(ADH_PETSc_Interface* self,
			    ADH_OneLevelTransport* cadh_transport)
{
  int nsys;
  if (NS_FLOW)
    {
      nsys = 4;
      ADH_PETSc_matrix_resize(self,cadh_transport,&self->ksp[NS_FLAG],&self->L[NS_FLAG],nsys,&self->du[NS_FLAG],&self->r[NS_FLAG]);
    }
  if (GW_FLOW)
    {
      nsys=1;
      ADH_PETSc_matrix_resize(self,cadh_transport,&self->ksp[GW_FLAG],&self->L[GW_FLAG],nsys,&self->du[GW_FLAG],&self->r[GW_FLAG]);
    }
  if (ntransport > 0)
    {
      nsys = 1;
      ADH_PETSc_matrix_resize(self,cadh_transport,&self->ksp[USR_FLAG],&self->L[USR_FLAG],nsys,&self->du[USR_FLAG],&self->r[USR_FLAG]);
    }
  if (HEAT)
    {
      nsys=1; 
      ADH_PETSc_matrix_resize(self,cadh_transport,&self->ksp[HT_FLAG],&self->L[HT_FLAG],nsys,&self->du[HT_FLAG],&self->r[HT_FLAG]);
    }
  if (OL_FLOW)
    {
      nsys=1;
      ADH_PETSc_matrix_resize(self,cadh_transport,&self->ksp[OL_FLAG],&self->L[OL_FLAG],nsys,&self->du[OL_FLAG],&self->r[OL_FLAG]);
    }
  if (SW2_FLOW)
    {
      nsys=3;
      ADH_PETSc_matrix_resize(self,cadh_transport,&self->ksp[SW2_FLAG],&self->L[SW2_FLAG],nsys,&self->du[SW2_FLAG],&self->r[SW2_FLAG]);
    }
  if (SW3_FLOW)
    {
      nsys=3;
      ADH_PETSc_matrix_resize(self,cadh_transport,&self->ksp[SW3_FLAG],&self->L[SW3_FLAG],nsys,&self->du[SW3_FLAG],&self->r[SW3_FLAG]);
    }
  if (SEDIMENT)
    {
      nsys=1;
      ADH_PETSc_matrix_resize(self,cadh_transport,&self->ksp[SED_FLAG],&self->L[SED_FLAG],nsys,&self->du[SED_FLAG],&self->r[SED_FLAG]);
    }
  if (MOVING_GRID)
    {
      nsys=1;
      ADH_PETSc_matrix_resize(self,cadh_transport,&self->ksp[0],&self->L[0],nsys,&self->du[0],&self->r[0]);
    }
  return 0;
}

int ADH_PETSc_destroyStorage(ADH_PETSc_Interface* self)
{
  int nsys;
  if (NS_FLOW)
    {
      nsys = 4;
      ADH_PETSc_matrix_destroy(&self->ksp[NS_FLAG],&self->L[NS_FLAG],&self->du[NS_FLAG],&self->r[NS_FLAG]);
    }
  if (GW_FLOW)
    {
      nsys=1;
      ADH_PETSc_matrix_destroy(&self->ksp[GW_FLAG],&self->L[GW_FLAG],&self->du[GW_FLAG],&self->r[GW_FLAG]);
    }
  if (ntransport > 0)
    {
      nsys = 1;
      ADH_PETSc_matrix_destroy(&self->ksp[USR_FLAG],&self->L[USR_FLAG],&self->du[USR_FLAG],&self->r[USR_FLAG]);
    }
  if (HEAT)
    {
      nsys=1; 
      ADH_PETSc_matrix_destroy(&self->ksp[HT_FLAG],&self->L[HT_FLAG],&self->du[HT_FLAG],&self->r[HT_FLAG]);
    }
  if (OL_FLOW)
    {
      nsys=1;
      ADH_PETSc_matrix_destroy(&self->ksp[OL_FLAG],&self->L[OL_FLAG],&self->du[OL_FLAG],&self->r[OL_FLAG]);
    }
  if (SW2_FLOW)
    {
      nsys=3;
      ADH_PETSc_matrix_destroy(&self->ksp[SW2_FLAG],&self->L[SW2_FLAG],&self->du[SW2_FLAG],&self->r[SW2_FLAG]);
    }
  if (SW3_FLOW)
    {
      nsys=3;
      ADH_PETSc_matrix_destroy(&self->ksp[SW3_FLAG],&self->L[SW3_FLAG],&self->du[SW3_FLAG],&self->r[SW3_FLAG]);
    }
  if (SEDIMENT)
    {
      nsys=1;
      ADH_PETSc_matrix_destroy(&self->ksp[SED_FLAG],&self->L[SED_FLAG],&self->du[SED_FLAG],&self->r[SED_FLAG]);
    }
  if (MOVING_GRID)
    {
      nsys=1;
      ADH_PETSc_matrix_destroy(&self->ksp[0],&self->L[0],&self->du[0],&self->r[0]);
    }
  return 0;
}

int ADH_PETSc_Interface_update(ADH_PETSc_Interface* self,
                               ADH_OneLevelTransport* cadh_transport)
{
  /* petsc_mapping_update */
  int i;
  self->mappingHasChanged = FALSE;
  /* save the old nnode on each processor */
  for(i = 0; i < npes; i++)
    self->nnodeLocal_old[i] = self->nnodeLocal[i];
  /* get the number of local (non-ghost) nodes on each processor */
  MPI_Allgather(&my_nnode, 1, MPI_INT, self->nnodeLocal, 1, MPI_INT,
		PETSC_COMM_WORLD);
  /* check to see if the number of nodes has changed on any processor */
  for(i = 0; i < npes; i++)
    {
      if(self->nnodeLocal[i] != self->nnodeLocal_old[i])
	self->mappingHasChanged = TRUE;
    }
  if(self->mappingHasChanged)
    {
      ISLocalToGlobalMappingDestroy(self->localToGlobalMapping);
      if(nnode > self->nnodeLocal_old[myid])
	self->globalNodeNumber =
	  (int *)tl_realloc(sizeof(int), nnode, self->nnodeLocal_old[myid],
			    (void *)self->globalNodeNumber);
      self->globalNodeNumberStart[0] = 0;
      for(i = 1; i < npes; i++)
	self->globalNodeNumberStart[i] =
	  self->globalNodeNumberStart[i - 1] + self->nnodeLocal[i - 1];
      for(i = 0; i < nnode; i++)
	self->globalNodeNumber[i] =
	  self->globalNodeNumberStart[node_pair[i].sd] + node_pair[i].rnode;
      ISLocalToGlobalMappingCreate(PETSC_COMM_WORLD, nnode,
				   self->globalNodeNumber,
				   &self->localToGlobalMapping);
    }
  /* solver, matrices, and vectors for each type of physics */
  ADH_PETSc_updateStorage(self,cadh_transport);
  return 0;
}

int ADH_PETSc_Interface_dealloc(ADH_PETSc_Interface* self)
{
  /* solve_petsc_*/
  printf("trying to dealloc ADH_PETSC \n");
  ADH_PETSc_destroyStorage(self);
  ISLocalToGlobalMappingDestroy(self->localToGlobalMapping);
  tl_free(sizeof(int), my_nnode, self->globalNodeNumber);
  tl_free(sizeof(int), npes, self->nnodeLocal);
  tl_free(sizeof(int), npes, self->globalNodeNumberStart);
  tl_free(sizeof(int), npes, self->nnodeLocal_old);
  return 0;
}
/*cek another broken up fe_main to allow calling just the resid and load without fe_newton*/

/** Call the ADH Newton solves for all the models. This will kick out
    with a failure for any model so the caller needs handle nonlinear
    solver failures */
int ADH_OneLevelTransport_updateStorage(ADH_OneLevelTransport* self)
{
 /*  int i, j, k, ie, n;		/\* loop counter *\/ */
/*   int itrn;			/\* loop counter over the transported quantities *\/ */
/*   int nsys;			/\* the number of systems being solved *\/ */
/*   int nsys_sq;			/\* the number of systems being solved squared *\/ */
/*   /\* int isize_prev = 0;		/\\* the previous isize *\\/ *\/ */
/*   int sed_iteration;		/\* the number of iterations done to balance the fractions in active layer *\/ */
/*   int sed_number;		/\* given the number of clay or silt what is its number in sediment *\/ */
/*   double ald_check;		/\* maximum change in the active layer distribution since the last iteration *\/ */
/*   double ald_max = DONE; */
/*   double t_sed;			/\* the real time for the sediment current sediment step *\/ */
/*   double flowrate_x;		/\* the flowrate for the hydro results to be kept for the sed timesteps *\/ */
/*   double flowrate_y;		/\* the flowrate for the hydro results to be kept for the sed timesteps *\/ */
/*   double ald_diff; */
/*   int max_node; */
/*   double sed_it_tol; */
/*   double bed_thick_old; */
/*   double bed_thick_new; */
/*   double esrh;			/\* the equivalent sand roughness height *\/ */
/*   double blsg; */
/*   double albd; */
/*   double asbd; */
/*   VECT2D tvel;			/\* the local transverse velocity components *\/ */
/*   int count = 0;		/\* counter for sed time step *\/ */
/*   int istr;			/\* counter for ceq boundary condition *\/ */
/* #ifdef _MESSG */
/*   double sed_time_start; */
/*   double sed_time_end; */
/*   double sed_time_total; */
/* #endif */
/*   double qus_depth_factor = 0.1; */
/*   double qus_max_change; */
/*   /\* double sdtemp; *\/ */
/*   /\* double temp_dpl; *\/ */
/*   /\* double ust; *\/ */
/*   /\* double blv; *\/ */
/*   double blvelmag; */
  /* reset the matrix */
  realloc_fe_main(self);
  return 0;
}

