/******************************************************************************/
/*** file:   starspotmodule.c                                               ***/
/*** author: X. Bonfils                                                     ***/
/***         Centro de Astronomia e Astrofisica da Universidade de Lisboa   ***/
/***                                                                        ***/
/*** version: 0.1 2006/08/29                                                ***/
/******************************************************************************/

#include <Python.h>
#include <numpy/arrayobject.h>
#include <starspot.h>


static PyObject *phasm_gauss_bis(PyObject *self, PyObject *args)
{
    /* Declarations */
    PyArrayObject *aa, *bb, *cc;
    PyArrayObject *dd,*ee,*ff;
    int  nstep=100, margin=5, len_depth=nstep-2*margin+1; // Bisectors properties
    int i, n1, dim_mod[1], dim_bis[1];
    double *depth, *bis, *mod, *para;
    double *vrad_ccf, *ccf, *err_ccf;
    
    depth    = (double *)malloc(sizeof(double)*len_depth);
    bis      = (double *)malloc(sizeof(double)*len_depth);
    
    // Initialization...
    for (i=0; i<(nstep-2*margin+1); i++) depth[i] = (double )(i+margin)/nstep;
    
    /* Parsing arguments */
    if (!PyArg_ParseTuple(args, "O!O!O!",&PyArray_Type, &aa, &PyArray_Type, &bb, &PyArray_Type, &cc))
    return NULL;
    
    /* Allocations */
    vrad_ccf = (double *) (aa->data + 0*aa->strides[0]);
    ccf      = (double *) (bb->data + 0*bb->strides[0]);
    err_ccf  = (double *) (cc->data + 0*cc->strides[0]);
    n1 = (int)aa->dimensions[0];
    
    mod  = (double *)malloc(sizeof(double)*n1);
    para = (double *)malloc(sizeof(double)*9);
    
    /* Calculate Gaussian fit on the CCF, as well as the bisector of the CCF */
    gauss_bis(vrad_ccf,ccf,err_ccf,n1,mod,para,depth,bis,len_depth);

    dim_mod[0] = n1;
    dim_bis[0] = len_depth;
    dd = (PyArrayObject *) PyArray_FromDims(1, dim_mod, PyArray_DOUBLE);
    ee = (PyArrayObject *) PyArray_FromDims(1, dim_bis, PyArray_DOUBLE);
    ff = (PyArrayObject *) PyArray_FromDims(1, dim_bis, PyArray_DOUBLE);
    for (i = 0; i < len_depth; i++)
    {
        *(double *) (dd->data + i*dd->strides[0]) = mod[i];
        *(double *) (ee->data + i*ee->strides[0]) = depth[i];
        *(double *) (ff->data + i*ff->strides[0]) = bis[i];
    }
    
    free(depth); free(bis); free(mod);
    
    // return parameter of the Gaussian fit and the bisector calculation
    // model,cte,ampli,vrad,fwhm,sig_cte,sig_ampli,sig_vrad_sig_fwhm,span,bis,depth
    return Py_BuildValue("NdddddddddNN", dd, para[0], para[1], para[2], para[3], para[4], para[5], para[6], para[7], para[8],ff,ee);
}




static PyObject *phasm_itot(PyObject *self, PyObject *args)
{
  /* Declarations */
  PyArrayObject *f_star2,*a, *b;
  double v, i, limba1, limba2, modif_bis_quad, modif_bis_lin, modif_bis_cte, v_interval;
  int grid, n_v, n, j;
  double *f_star,*vrad_ccf, *intensity_ccf, sum_star;
  
  /* Parsing arguments */
  if (!PyArg_ParseTuple(args, "dddddddiO!O!dii", &v, &i, &limba1, &limba2, &modif_bis_quad, &modif_bis_lin, &modif_bis_cte, &grid,
			&PyArray_Type, &a, &PyArray_Type, &b, &v_interval, &n_v, &n)) 
    return NULL;

  /* Allocations */
  f_star2 = (PyArrayObject *) PyArray_FromDims(1, &n_v, PyArray_DOUBLE);
  f_star  = (double *)f_star2->data;
  vrad_ccf = (double *) (a->data + 0*a->strides[0]);
  intensity_ccf = (double *) (b->data + 0*b->strides[0]);

  /* Initialisations */
  for (j=0; j<n_v; j++) f_star[j] = 0;
  sum_star = 0;

  /* Total intensity */
  itot(v,i,limba1,limba2, modif_bis_quad, modif_bis_lin, modif_bis_cte,grid,vrad_ccf,intensity_ccf,v_interval,n_v,n,f_star,&sum_star);

  return Py_BuildValue("Nd", f_star2, sum_star);
}

static PyObject *phasm_starmap(PyObject *self, PyObject *args)
{
  /* Declarations */
  PyArrayObject *f_map2, *v_map2;
  double v, i, limba1, limba2;
  int grid, j, k, dimensions[2];
  double **Fmap, **Vmap;
  
  /* Parsing arguments */
  if (!PyArg_ParseTuple(args, "ddddi", &v, &i, &limba1, &limba2, &grid)) 
    return NULL;

  /* Allocations */
  Fmap = (double **)malloc(sizeof(double *)*grid);
  Vmap = (double **)malloc(sizeof(double *)*grid);
  for (j=0; j<grid; j++) Fmap[j] = (double *)malloc(sizeof(double)*grid);
  for (j=0; j<grid; j++) Vmap[j] = (double *)malloc(sizeof(double)*grid);
  
  /* Initialisations */
  for (j=0; j<grid; j++) for (k=0; k<grid; k++) Fmap[j][k] = 0.;
  for (j=0; j<grid; j++) for (k=0; k<grid; k++) Vmap[j][k] = 0.;

  /* Star Maps (Intensity & Velocity Maps) */
  starmap(v, i, limba1, limba2, grid, Fmap, Vmap);

  /* Creation of the returned python object (Numeric array)*/
  dimensions[0] = grid; dimensions[1] = grid;
	
  f_map2 = (PyArrayObject *) PyArray_FromDims(2, dimensions, PyArray_DOUBLE);
  v_map2 = (PyArrayObject *) PyArray_FromDims(2, dimensions, PyArray_DOUBLE);

  for (j=0; j<dimensions[0]; j++) for (k=0; k<dimensions[1]; k++)
    *(double *)(f_map2->data + j*f_map2->strides[0] + k*f_map2->strides[1]) = Fmap[j][k];
  for (j=0; j<dimensions[0]; j++) for (k=0; k<dimensions[1]; k++)
    *(double *)(v_map2->data + j*v_map2->strides[0] + k*v_map2->strides[1]) = Vmap[j][k];

  for (j=0; j<grid; j++) free(Fmap[j]); free(Fmap);
  for (j=0; j<grid; j++) free(Vmap[j]); free(Vmap);

  return Py_BuildValue("NN", f_map2, v_map2);
}

static PyObject *phasm_spot_init(PyObject *self, PyObject *args)
{
  /* Declarations */
  PyArrayObject *xyz2;
  double s, longitude, latitude, inclination;
  int nrho, i, j;
  double **xyz;
  
  /* Parsing arguments */
  if (!PyArg_ParseTuple(args, "ddddi", &s, &longitude, &latitude, &inclination, 
			&nrho))
    return NULL;
 
  /* Allocations & Init.*/
  xyz    = (double **)malloc(sizeof(double *)*nrho);
  for (j=0; j<nrho; j++) xyz[j] = (double *)malloc(sizeof(double)*3);

  /* Init Spot */
  spot_init(s, longitude, latitude, inclination, nrho, xyz);

  /* Creation of the returned python object (Numeric array)                    */
  int dimensions[2]={nrho,3};
  xyz2 = (PyArrayObject *) PyArray_FromDims(2, dimensions, PyArray_DOUBLE);
  for (i=0; i<dimensions[0]; i++) for (j=0; j<dimensions[1]; j++)
    *(double *)(xyz2->data + i*xyz2->strides[0] + j*xyz2->strides[1]) = xyz[i][j];

  for (j=0; j<nrho; j++) free(xyz[j]); free(xyz);  
  
  return Py_BuildValue("N", xyz2);
}

static PyObject *phasm_spot_phase(PyObject *self, PyObject *args)
{
  /* Declarations */
  PyArrayObject *xyz_in, *xyz_out;
  double inclination, phase;
  int nrho, i, j;
  double **xyz, **xyz2;
  
  /* Parsing arguments */
  if (!PyArg_ParseTuple(args, "O!did", &PyArray_Type, &xyz_in, &inclination, 
			&nrho, &phase))
    return NULL;
 
  /* Allocations & Init.*/
  xyz    = (double **)malloc(sizeof(double *)*nrho);
  xyz2   = (double **)malloc(sizeof(double *)*nrho);
  for (j=0; j<nrho; j++) xyz[j]  = (double *)malloc(sizeof(double)*3);
  for (j=0; j<nrho; j++) xyz2[j] = (double *)malloc(sizeof(double)*3);
  for (i=0; i<nrho; i++) for (j=0; j<3; j++) 
    xyz[i][j] = *(double *)(xyz_in->data 
			    + i*xyz_in->strides[0] + j*xyz_in->strides[1]);

  /* Init Spot */
  spot_phase(xyz, inclination, nrho, phase, xyz2);

  /* Creation of the returned python object (Numeric array)                    */
  int dimensions[2]={nrho,3};
  xyz_out = (PyArrayObject *) PyArray_FromDims(2, dimensions, PyArray_DOUBLE);
  for (i=0; i<dimensions[0]; i++) for (j=0; j<dimensions[1]; j++)
    *(double *)(xyz_out->data + i*xyz_out->strides[0] + j*xyz_out->strides[1]) = xyz2[i][j];
    
  for (j=0; j<nrho; j++) free(xyz[j]); free(xyz);
  for (j=0; j<nrho; j++) free(xyz2[j]); free(xyz2);
  return Py_BuildValue("N", xyz_out);
}

static PyObject *phasm_spot_area(PyObject *self, PyObject *args)
{
  /* Declarations */
  PyArrayObject *xlylzl_in;
  int nrho, i, j, vis, grid;
  int iminy, iminz, imaxy, imaxz;
  double **xlylzl;
  
  /* Parsing arguments */
  if (!PyArg_ParseTuple(args, "O!ii", &PyArray_Type, &xlylzl_in, &nrho, &grid))
    return NULL;
 
  /* Allocations & Init.*/
  xlylzl    = (double **)malloc(sizeof(double *)*nrho);
  for (j=0; j<nrho; j++) xlylzl[j]  = (double *)malloc(sizeof(double)*3);
  for (i=0; i<nrho; i++) for (j=0; j<3; j++) 
    xlylzl[i][j] = *(double *)(xlylzl_in->data 
			    + i*xlylzl_in->strides[0] + j*xlylzl_in->strides[1]);

  /* Spot Area */
  vis = spot_area(xlylzl, nrho, grid, &iminy, &iminz, &imaxy, &imaxz);

  for (j=0; j<nrho; j++) free(xlylzl[j]); free(xlylzl);
  return Py_BuildValue("iiiii", vis, iminy, iminz, imaxy, imaxz);
}

static PyObject *phasm_spot_scan(PyObject *self, PyObject *args)
{
  /* Declarations */
  PyArrayObject *a, *b, *c, *f_spot2, *f_spot3, *f_spot4;
  double v, i, limba1, limba2, modif_bis_quad, modif_bis_lin, modif_bis_cte, v_interval;
  int grid, n_v, n;
  double s, longitude, phase, latitude;
  int iminy, iminz, imaxy, imaxz, magn_feature_type,T_star,T_diff_spot; 
  double *vrad_ccf, *intensity_ccf, *intensity_ccf_spot, *f_spot_flux,*f_spot_bconv,*f_spot_tot, sum_spot;

  /* Parsing arguments */
  if (!PyArg_ParseTuple(args, "dddddddiO!O!O!diiddddiiiiiii", &v, &i, &limba1, &limba2, &modif_bis_quad, &modif_bis_lin, &modif_bis_cte, &grid, 
      &PyArray_Type, &a,&PyArray_Type, &b,&PyArray_Type, &c, &v_interval, &n_v, &n, &s, &longitude, 
      &phase, &latitude, &iminy, &iminz, &imaxy, &imaxz,
      &magn_feature_type, &T_star, &T_diff_spot))
    return NULL;
 
  /* Allocations */
  vrad_ccf = (double *) (a->data + 0*a->strides[0]);
  intensity_ccf = (double *) (b->data + 0*b->strides[0]);
  intensity_ccf_spot = (double *) (c->data + 0*c->strides[0]);
  f_spot2 = (PyArrayObject *) PyArray_FromDims(1, &n_v, PyArray_DOUBLE);
  f_spot_flux = (double *)f_spot2->data;
  f_spot3 = (PyArrayObject *) PyArray_FromDims(1, &n_v, PyArray_DOUBLE);
  f_spot_bconv = (double *)f_spot3->data;
  f_spot4 = (PyArrayObject *) PyArray_FromDims(1, &n_v, PyArray_DOUBLE);
  f_spot_tot = (double *)f_spot4->data;

  sum_spot = 0;

  /* Scan Spot */
  spot_scan(v, i, limba1, limba2, modif_bis_quad, modif_bis_lin, modif_bis_cte, grid, vrad_ccf, intensity_ccf, intensity_ccf_spot, v_interval, n_v, n, s, longitude, 
            phase, latitude, iminy, iminz, imaxy, imaxz, f_spot_flux, f_spot_bconv,f_spot_tot, 
            &sum_spot,magn_feature_type,T_star,T_diff_spot);

  /*return Py_BuildValue("Nd", f_spot2, sum_spot);*/
  return Py_BuildValue("NNNd", f_spot2, f_spot3,f_spot4, sum_spot);
}

static PyObject *phasm_spot_inverse_rotation(PyObject *self, PyObject *args)
{
  /* Declarations */
  PyArrayObject *xyz_in, *xiyizi2;
  double longitude, latitude, inclination, phase;
  int i, dimensions[1];
  double *xyz, *xiyizi;
  
  /* Parsing arguments */
  if (!PyArg_ParseTuple(args, "O!dddd", &PyArray_Type, &xyz_in, &longitude, 
      &latitude, &inclination, &phase))
    return NULL;
 
  /* Allocations & Init.*/
  dimensions[0] = 3;
  xiyizi2 = (PyArrayObject *) PyArray_FromDims(1, dimensions, PyArray_DOUBLE);
  xiyizi  = (double *)xiyizi2->data;
  xyz    = (double *)malloc(sizeof(double *)*3);
  for (i=0; i<3; i++) xyz[i] = *(double *)(xyz_in->data + i*xyz_in->strides[0]);
  for (i=0; i<3; i++) xiyizi[i] = 0.;
   
  /* Init Spot */
  spot_inverse_rotation(xyz, longitude, latitude, inclination, phase, xiyizi);

  free(xyz);
  return Py_BuildValue("N", xiyizi2);
}

static PyObject *phasm_spot_scan_npsi(PyObject *self, PyObject *args)
{
  /* Declarations */
  PyArrayObject *a,*b,*c, *xyz_in, *psi_in, *f_spot2, *f_spot3, *f_spot4, *sum_spot2;
  int i,j;   int dimensions[2];
  int nrho, npsi, grid, n_v, n, magn_feature_type,T_star,T_diff_spot;
  double v, inclination, limba1, limba2, modif_bis_quad, modif_bis_lin, modif_bis_cte, *vrad_ccf, *intensity_ccf, *intensity_ccf_spot, v_interval, s, longitude;
  double latitude;
  double **xyz, *psi, **f_spot_flux, **f_spot_bconv, **f_spot_tot, *sum_spot;
  
  /* Parsing arguments */
  if (!PyArg_ParseTuple(args, "O!iO!idddddddiO!O!O!diidddiii", &PyArray_Type, &xyz_in, 
        &nrho, &PyArray_Type, &psi_in, &npsi, &v, &inclination, &limba1, &limba2, &modif_bis_quad, &modif_bis_lin, &modif_bis_cte, &grid,
		&PyArray_Type, &a, &PyArray_Type, &b, &PyArray_Type, &c, &v_interval, &n_v, &n, &s, &longitude, &latitude,
	    &magn_feature_type,&T_star,&T_diff_spot))
    return NULL;
     
  /* Allocations */
  vrad_ccf = (double *) (a->data + 0*a->strides[0]);
  intensity_ccf = (double *) (b->data + 0*b->strides[0]);
  intensity_ccf_spot = (double *) (c->data + 0*c->strides[0]);
    
  sum_spot2 = (PyArrayObject *) PyArray_FromDims(1, &npsi, PyArray_DOUBLE);
  sum_spot  = (double *)sum_spot2->data;
  f_spot_flux    = (double **)malloc(sizeof(double *)*npsi);
  f_spot_bconv    = (double **)malloc(sizeof(double *)*npsi);
  f_spot_tot    = (double **)malloc(sizeof(double *)*npsi);
  for (j=0; j<npsi; j++) f_spot_flux[j] = (double *)malloc(sizeof(double)*n_v);
  for (j=0; j<npsi; j++) f_spot_bconv[j] = (double *)malloc(sizeof(double)*n_v);
  for (j=0; j<npsi; j++) f_spot_tot[j] = (double *)malloc(sizeof(double)*n_v);
  xyz    = (double **)malloc(sizeof(double *)*nrho);
  for (j=0; j<nrho; j++) xyz[j]  = (double *)malloc(sizeof(double)*3);
  for (i=0; i<nrho; i++) for (j=0; j<3; j++) 
    xyz[i][j] = *(double *)(xyz_in->data 
			    + i*xyz_in->strides[0] + j*xyz_in->strides[1]);
  psi = (double *)malloc(sizeof(double)*npsi);
  for (i=0; i<npsi; i++) psi[i] = *(double *)(psi_in->data + i*psi_in->strides[0]);

  /* Initialisations */
  for (i=0; i<npsi; i++) for (j=0; j<n_v; j++) f_spot_flux[i][j] = 0.;
  for (i=0; i<npsi; i++) for (j=0; j<n_v; j++) f_spot_bconv[i][j] = 0.;
  for (i=0; i<npsi; i++) for (j=0; j<n_v; j++) f_spot_tot[i][j] = 0.;
  for (i=0; i<npsi; i++) sum_spot[i] = 0.;

  /* Scan Spot */
  spot_scan_npsi(xyz, nrho, psi, npsi, v, inclination, limba1, limba2, modif_bis_quad, modif_bis_lin, modif_bis_cte, grid, 
				 vrad_ccf, intensity_ccf, intensity_ccf_spot, v_interval, n_v, n, s, longitude, latitude, 
	      		 f_spot_flux, f_spot_bconv, f_spot_tot, sum_spot, magn_feature_type,T_star,T_diff_spot);

  /* Creation of the returned python object (Numeric array)*/
  dimensions[0] = npsi; dimensions[1] = n_v;
  f_spot2 = (PyArrayObject *) PyArray_FromDims(2, dimensions, PyArray_DOUBLE);
  f_spot3 = (PyArrayObject *) PyArray_FromDims(2, dimensions, PyArray_DOUBLE);
  f_spot4 = (PyArrayObject *) PyArray_FromDims(2, dimensions, PyArray_DOUBLE);
  for (j=0; j<dimensions[0]; j++) for (i=0; i<dimensions[1]; i++)
  {
    *(double *)(f_spot2->data + j*f_spot2->strides[0] + i*f_spot2->strides[1]) = f_spot_flux[j][i];
    *(double *)(f_spot3->data + j*f_spot3->strides[0] + i*f_spot3->strides[1]) = f_spot_bconv[j][i];
    *(double *)(f_spot4->data + j*f_spot4->strides[0] + i*f_spot4->strides[1]) = f_spot_tot[j][i];
  }
  for (j=0; j<npsi; j++) free(f_spot_flux[j]);
  for (j=0; j<npsi; j++) free(f_spot_bconv[j]);
  for (j=0; j<npsi; j++) free(f_spot_tot[j]);
  for (j=0; j<nrho; j++) free(xyz[j]);
  free(f_spot_flux); free(f_spot_bconv); free(f_spot_tot); free(xyz); free(psi);
  
  return Py_BuildValue("NNNN", f_spot2, f_spot3, f_spot4, sum_spot2);
}

static PyMethodDef PhasmMethods[] = {
    {"gauss_bis", phasm_gauss_bis, METH_VARARGS,
        "Fit of a Gaussian on the CCF and calculate the bisector"},
    {"itot", phasm_itot, METH_VARARGS, 
     "Calculation of total intensity with no spot"},
    {"starmap", phasm_starmap, METH_VARARGS, 
     "Returns intensity & velocity maps of the star"},
    {"spot_init", phasm_spot_init, METH_VARARGS, "Init spot"},
    {"spot_phase", phasm_spot_phase, METH_VARARGS, "Phase spot"},
    {"spot_area", phasm_spot_area, METH_VARARGS, 
     "Determine a smaller area to locate the spot"},
    {"spot_scan", phasm_spot_scan, METH_VARARGS, 
     "Scan spot non-contribution to the flux and velocity field of the star"},
    {"spot_inverse_rotation", phasm_spot_inverse_rotation, METH_VARARGS, 
     "Inverse rotaion of xyz to check whether it belongs to the spot or not"},
    {"spot_scan_npsi", phasm_spot_scan_npsi, METH_VARARGS, 
     "Same as spot_scan for an array of phases psi"},
    {NULL, NULL, 0, NULL}
};

void initstarspot(void){
    (void) Py_InitModule("starspot", PhasmMethods);
    import_array();
}
