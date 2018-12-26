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
static double Profile (double, double);

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

      rs=sqrt((x1-g_inputParam[cloudx])*(x1-g_inputParam[cloudx])+(x2-g_inputParam[cloudy])*(x2-g_inputParam[cloudy]));
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

void Analysis (const Data *d, Grid *grid)
{
}
#if PHYSICS == MHD
/* ********************************************************************* */
void BackgroundField (double x1, double x2, double x3, double *B0)
{
   B0[0] = 0.0;
   B0[1] = 0.0;
   B0[2] = 0.0;
}
#endif

/* ********************************************************************* */
void UserDefBoundary (const Data *d, RBox *box, int side, Grid *grid)
{
  int   i, j, k, nv;
  double  *x1, *x2, *x3;
  double B0,B,P,Rj,Pj,r0,rstar,beta;

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

                    d->Vc[VX1][k][j][i] = 0.0;
                    d->Vc[VX2][k][j][i] = sqrt(1.0-(1.0/(g_inputParam[lorentz_jet]*g_inputParam[lorentz_jet])));
                    d->Vc[BX1][k][j][i] = 0.0;
                    d->Vc[BX2][k][j][i] = 0.0;
                    Pj=g_inputParam[pressure_jet_thermal];
                    if (g_inputParam[gBphi]==1) {
                        Rj=g_inputParam[jet_window];
                        r0=g_inputParam[r0toR]*Rj;
                        B0=sqrt(2.*Pj*(1.0+(Rj/r0)*(Rj/r0)));
                        B=-1.0*g_inputParam[lorentz_jet ]*B0*(x1[i]/r0)/(1.0+pow((x1[i]/r0),2));
                        //printf('B0 = ',B0);
                        d->Vc[BX3][k][j][i] = B;
                        P=B0*B0/pow((2.0*(1.0+(x1[i]/r0)*(x1[i]/r0))),2); //+B*B/(2.*g_inputParam[lorentz_jet])*(1.0 -  MIN(x1[i]*x1[i],1.0));
                        if (P<1e-13) P=1e-13;
                        d->Vc[PRS][k][j][i] = P;
                    } else {
                        d->Vc[BX3][k][j][i] = 0.0;
                        d->Vc[PRS][k][j][i] = Pj;
                    }
                } else {
                    VAR_LOOP(nv) d->Vc[nv][k][j][i] = d->Vc[nv][k][JBEG][i];
                }
            } else {
                VAR_LOOP(nv) d->Vc[nv][k][j][i] = d->Vc[nv][k][JBEG][i];
            }
    }

    }else if (box->vpos == X1FACE){
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

/* ********************************************************************* */
double Profile (double R, double nv)
/*
 *
 *
 *********************************************************************** */
{
  double R4 = R*R*R*R, R8 = R4*R4;

#if PHYSICS == MHD && COMPONENTS == 3    /* Start with a smoother profile */
  if (g_time < 0.1)  return 1.0/cosh(R4); /* with toroidal magnetic fields */
#endif
  return 1.0/cosh(R8);
}
