/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Contains basic functions for problem initialization.

  The init.c file collects most of the user-supplied functions useful
  for problem configuration.
  It is automatically searched for by the makefile.

  \author A. Mignone (mignone@ph.unito.it)
  \date   Sept 10, 2012

  * Edited for Molecular Cloud Cooling/ Jet interaction *
    		Version 1.0 (7/3)
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

static void   GetJetValues (double x1, double x2, double x3, double *vj);

/* ********************************************************************* */
void Init (double *v, double x1, double x2, double x3)
/*!
 * The Init() function can be used to assign initial conditions as
 * as a function of spatial position.
 *
 * \param [out] v   a pointer to a vector of primitive variables
 * \param [in] x1   coordinate point in the 1st dimension
 * \param [in] x2   coordinate point in the 2nd dimension
 * \param [in] x3   coordinate point in the 3rdt dimension
 *
 * The meaning of x1, x2 and x3 depends on the geometry:
 * \f[ \begin{array}{cccl}
 *    x_1  & x_2    & x_3  & \mathrm{Geometry}    \\ \noalign{\medskip}
 *     \hline
 *    x    &   y    &  z   & \mathrm{Cartesian}   \\ \noalign{\medskip}
 *    R    &   z    &  -   & \mathrm{cylindrical} \\ \noalign{\medskip}
 *    R    & \phi   &  z   & \mathrm{polar}       \\ \noalign{\medskip}
 *    r    & \theta & \phi & \mathrm{spherical}
 *    \end{array}
 *  \f]
 *
 * Variable names are accessed by means of an index v[nv], where
 * nv = RHO is density, nv = PRS is pressure, nv = (VX1, VX2, VX3) are
 * the three components of velocity, and so forth.
 *
 *********************************************************************** */
{
      g_gamma=5.0/3.0;


      /* -- ISM values ---- */
      v[RHO] = g_inputParam[rho_ism];
      EXPAND(v[VX1] = 0.0; ,
             v[VX2] = 0.0; ,
             v[VX3] = 0.0; )
      #if HAVE_ENERGY
      v[PRS] = g_inputParam[pressure_gas];
      #endif
      v[TRC] = 0.0;

      /* one cloud */
      double rs,A,B,r0,z0;
      A=10.0;
      B=0.002;
      r0=2.0;
      z0=3.0;

      rs=sqrt((x1-r0)*(x1-r0)+(x2-z0)*(x2-z0));
      if (rs<g_inputParam[radius_cloud]) {
        // Plummer Profile
        v[RHO]=A/(B+pow(rs,2.3));
      }


      #if PHYSICS == MHD || PHYSICS == RMHD

       v[BX1] = 0.0;
       v[BX2] = 0.0;
       v[BX3] = 0.0;

       v[AX1] = 0.0;
       v[AX2] = 0.0;
       v[AX3] = 0.0;

      #endif
}

/* ********************************************************************* */
void InitDomain (Data *d, Grid *grid)
/*!
 * Assign initial condition by looping over the computational domain.
 * Called after the usual Init() function to assign initial conditions
 * on primitive variables.
 * Value assigned here will overwrite those prescribed during Init().
 *
 *
 *********************************************************************** */
{
}

/* ********************************************************************* */
void Analysis (const Data *d, Grid *grid)
/*!
 *  Perform runtime data analysis.
 *
 * \param [in] d the PLUTO Data structure
 * \param [in] grid   pointer to array of Grid structures
 *
 *********************************************************************** */
{

}
#if PHYSICS == MHD
/* ********************************************************************* */
void BackgroundField (double x1, double x2, double x3, double *B0)
/*!
 * Define the component of a static, curl-free background
 * magnetic field.
 *
 * \param [in] x1  position in the 1st coordinate direction \f$x_1\f$
 * \param [in] x2  position in the 2nd coordinate direction \f$x_2\f$
 * \param [in] x3  position in the 3rd coordinate direction \f$x_3\f$
 * \param [out] B0 array containing the vector componens of the background
 *                 magnetic field
 *********************************************************************** */
{
   B0[0] = 0.0;
   B0[1] = 0.0;
   B0[2] = 0.0;
}
#endif

/* ********************************************************************* */
void UserDefBoundary (const Data *d, RBox *box, int side, Grid *grid)

/*!
 *  Assign user-defined boundary conditions.
 *
 * \param [in,out] d  pointer to the PLUTO data structure containing
 *                    cell-centered primitive quantities (d->Vc) and
 *                    staggered magnetic fields (d->Vs, when used) to
 *                    be filled.
 * \param [in] box    pointer to a RBox structure containing the lower
 *                    and upper indices of the ghost zone-centers/nodes
 *                    or edges at which data values should be assigned.
 * \param [in] side   specifies the boundary side where ghost zones need
 *                    to be filled. It can assume the following
 *                    pre-definite values: X1_BEG, X1_END,
 *                                         X2_BEG, X2_END,
 *                                         X3_BEG, X3_END.
 *                    The special value side == 0 is used to control
 *                    a region inside the computational domain.
 * \param [in] grid  pointer to an array of Grid structures.
 *
 *********************************************************************** */
{
  int   i, j, k, nv;
  double  *x1, *x2, *x3;
  double B0,Rj,Pj,r0,rstar,beta;

  //x1 = grid[IDIR].x;
  //x2 = grid[JDIR].x;
  //x3 = grid[KDIR].x;

  x1 = grid->xgc[IDIR];  /* -- array pointer to x1 coordinate -- */
  x2 = grid->xgc[JDIR];  /* -- array pointer to x2 coordinate -- */
  x3 = grid->xgc[KDIR];  /* -- array pointer to x3 coordinate -- */

  if (side == 0) {    /* -- check solution inside domain -- */
    DOM_LOOP(k,j,i){};
  }

  if (side == X1_BEG){  /* -- X1_BEG boundary -- */
    if (box->vpos == CENTER) {
      BOX_LOOP(box,k,j,i){ }
    }else if (box->vpos == X1FACE){
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X2FACE){
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X3FACE){
      BOX_LOOP(box,k,j,i){  }
    }
  }

  if (side == X1_END){  /* -- X1_END boundary -- */
    if (box->vpos == CENTER) {
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X1FACE){
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X2FACE){
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X3FACE){
      BOX_LOOP(box,k,j,i){  }
    }
  }

    if (side == X2_BEG){  /* -- X2_BEG boundary -- */
        if (box->vpos == CENTER) {
            BOX_LOOP(box,k,j,i){
            if (g_inputParam[gJET]==1) { // Boundary conditions IN THE JET
                if (fabs(x1[i])<=g_inputParam[jet_window]) {
                    //if (fabs(x1[i]) < 1.e-10) x1[i] = 1.e-10;
                    d->Vc[RHO][k][j][i] = g_inputParam[rho_jet];
                    d->Vc[PRS][k][j][i] = g_inputParam[pressure_jet_thermal];
                    d->Vc[VX1][k][j][i] = 0.0;
                    d->Vc[VX2][k][j][i] = sqrt(1.0-(1.0/(g_inputParam[lorentz_jet]*g_inputParam[lorentz_jet])));
                    d->Vc[BX1][k][j][i] = 0.0;
                    if (g_inputParam[gBphi]==1) {
                        Rj=g_inputParam[jet_window];
                        r0=0.2*Rj;
                        rstar=0.04*Rj;
                        Pj=g_inputParam[pressure_jet_thermal];
                        beta=0.05;
                        B0=sqrt(2.*Pj)/sqrt(beta) *(1.+pow((rstar/r0),2))/(rstar/r0);
                        //printf('B0 = ',B0);
                        d->Vc[BX2][k][j][i] = g_inputParam[lorentz_jet ]*B0*(x2[j]/r0)/(1.+pow((x2[j]/r0),2));
                    } else {
                        d->Vc[BX2][k][j][i] = 0.0;
                    }
                    //sqrt(2.*g_inputParam[pressure_jet_thermal])/sqrt(0.05) *(1.+pow(0.04*g_inputParam[jet_window]/+0.2*g_inputParam[jet_window],2))/(0.04*g_inputParam[jet_window]/+0.2*g_inputParam[jet_window]);
                    // d->Vc[BX3][k][j][i] = 0.0;
                    // d->Vc[AX1][k][j][i] = 0.0;
                    // d->Vc[AX2][k][j][i] = 0.0;
                    // d->Vc[AX3][k][j][i] = 0.0;
                } else {
                    VAR_LOOP(nv) d->Vc[nv][k][j][i] = d->Vc[nv][k][JBEG][i];
                }
            } else {
                VAR_LOOP(nv) d->Vc[nv][k][j][i] = d->Vc[nv][k][JBEG][i];
            }
    }

    }else if (box->vpos == X1FACE){
      // #ifdef STAGGERED_MHD
      //  x1 = grid->xr[IDIR];
      //  BOX_LOOP(box,k,j,i){
      //    vout[BX1] = -d->Vs[BX1s][k][2*JBEG - j - 1][i];
      //    d->Vs[BX1s][k][j][i] =    vout[BX1]
      //                           - (vout[BX1] - vjet[BX1])*Profile(fabs(x1[i]), BX1);
      //  }
      // #endif
    }else if (box->vpos == X2FACE){
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X3FACE){
      BOX_LOOP(box,k,j,i){  }
    }
  }

  if (side == X2_END){  /* -- X2_END boundary -- */
    if (box->vpos == CENTER) {
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X1FACE){
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X2FACE){
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X3FACE){
      BOX_LOOP(box,k,j,i){  }
    }
  }

  if (side == X3_BEG){  /* -- X3_BEG boundary -- */
    if (box->vpos == CENTER) {
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X1FACE){
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X2FACE){
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X3FACE){
      BOX_LOOP(box,k,j,i){  }
    }
  }

  if (side == X3_END){  /* -- X3_END boundary -- */
    if (box->vpos == CENTER) {
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X1FACE){
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X2FACE){
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X3FACE){
      BOX_LOOP(box,k,j,i){  }
    }
  }
}

#if BODY_FORCE != NO
/* ********************************************************************* */
void BodyForceVector(double *v, double *g, double x1, double x2, double x3)
/*!
 * Prescribe the acceleration vector as a function of the coordinates
 * and the vector of primitive variables *v.
 *
 * \param [in] v  pointer to a cell-centered vector of primitive
 *                variables
 * \param [out] g acceleration vector
 * \param [in] x1  position in the 1st coordinate direction \f$x_1\f$
 * \param [in] x2  position in the 2nd coordinate direction \f$x_2\f$
 * \param [in] x3  position in the 3rd coordinate direction \f$x_3\f$
 *
 *********************************************************************** */
{
	double M,rs,G,A;
	rs=sqrt(x1*x1+x2*x2);
	G=1.114e-13;
	A=g_inputParam[pl_rho1]*7.39198;
	if (rs>g_inputParam[radius_cloud]) {
	  M=A*pow(g_inputParam[radius_cloud],1.7);
	  g[IDIR] = -G*M*x1/rs/rs/rs;
	  g[JDIR] = -G*M*x2/rs/rs/rs;
	  g[KDIR] = 0.0;
	}
	else {
	  M=A*pow(rs,1.7);
	  g[IDIR]= -G*M*x1/rs/rs/rs;
	  g[JDIR]= -G*M*x2/rs/rs/rs;
	  g[KDIR]= 0.0;

	}
}
/* ********************************************************************* */
double BodyForcePotential(double x1, double x2, double x3)
/*!
 * Return the gravitational potential as function of the coordinates.
 *
 * \param [in] x1  position in the 1st coordinate direction \f$x_1\f$
 * \param [in] x2  position in the 2nd coordinate direction \f$x_2\f$
 * \param [in] x3  position in the 3rd coordinate direction \f$x_3\f$
 *
 * \return The body force potential \f$ \Phi(x_1,x_2,x_3) \f$.
 *
 *********************************************************************** */
{
    	double M,A,rs,G;
	rs=sqrt(x1*x1+x2*x2);
	G=1.114e-13;

	A=g_inputParam[pl_rho1]*7.39198;  //Power-Law Density Profile Integration
	//A=(4.0/3.0)*3.14156*g_inputParam[rho_cloud] //Box Density Profile Integration

	if (rs>g_inputParam[radius_cloud]) {
	  M=A*pow(g_inputParam[radius_cloud],1.7); //Power-Law Density Profile Integration
	  //M=A*pow(g_inputParam[radius_cloud],3); //Box Density Profile Integration
	  return -G*M/rs;
	}
	else {
	  M=A*pow(rs,1.7); //Power-Law Density Profile Integration
	  //M=A*pow(rs,3); //Box Density Profile Integration
	  return -G*M/rs;
	}
}
#endif

/* **************************************************************** */
void GetJetValues (double x1, double x2, double x3, double *vj)
/*
 *
 *
 ****************************************************************** */
{
  static int  first_call = 1;
  double lor,r,r0;

  r   = x1;
  r0 = 0.1* g_inputParam[jet_window];
  lor = g_inputParam[lorentz_jet];
  if (fabs(r) < 1.e-10) r = 1.e-10;

  vj[RHO] = g_inputParam[rho_jet];

  EXPAND(vj[VX1] = 0.0;                        ,
         vj[VX2] = sqrt(1.0 - 1.0/(lor*lor));  , /* 3-vel */
         vj[VX3] = 0.0;)

  vj[PRS] =  g_inputParam[pressure_jet_thermal];
  #if GEOMETRY == CYLINDRICAL
   EXPAND(vj[iBR]   = 0.0;                                 ,
          vj[iBZ]   = 0.0;  ,
          vj[iBPHI] = 0.0;)
   vj[AX1] = 0.0;
   vj[AX2] = 0.0;
   vj[AX3] = 0.0; //0.5*r*vj[iBZ];
  #endif

}
