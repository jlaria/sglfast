#include <iostream>
#include <math.h>

using namespace std;

extern "C" {

  int single = 0;

void VCalc(int *riskSetInd, int *riskSet, int *nDeath, int *nrow, int *ncol, double *eta, double *V)
{

  int ndeath = nDeath[0];
  double* C = NULL;
  C = new double[ndeath];
  double w;
  double r;
  int n = nrow[0];
  double* expEta = NULL;
  expEta = new double[n];

  for(int i = 0; i < n; i++)
    {
      expEta[i] = exp(eta[i]);
    }

  for(int i = 0; i < ndeath; i++)
    {
      C[i] = 0;
	  for(int j = (riskSetInd[i]-1); j < (riskSetInd[i+1]-1); j++)
	    {
	      C[i] = C[i] + expEta[j];
	    }
    }


  V[ndeath - 1] = C[ndeath - 1];
  
  for(int i = 2; i <= ndeath; i++) 
    {
      V[ndeath - i] = V[ndeath - i + 1] + C[ndeath - i];
    }
  delete [] C;
  delete [] expEta;
}



void gradCalc(int *riskSetInd, int *riskSet, int *numDeath, int *status, int *nDeath, int *nrow, int *ncol, double *eta, double *ldot, double *V)
{

  int ndeath = nDeath[0];
  double* C = NULL;
  C = new double[ndeath];
  double w;
  double r;
  int n = nrow[0];
  double* expEta = NULL;
  expEta = new double[n];

  for(int i = 0; i < n; i++)
    {
      expEta[i] = exp(eta[i]);
    }

  for(int i = 0; i < ndeath; i++)
    {
      C[i] = 0;
	  for(int j = (riskSetInd[i]-1); j < (riskSetInd[i+1]-1); j++)
	    {
	      C[i] = C[i] + expEta[j];
	    }
    }


  V[ndeath - 1] = C[ndeath - 1];
  
  for(int i = 2; i <= ndeath; i++) 
    {
      V[ndeath - i] = V[ndeath - i + 1] + C[ndeath - i];
    }

  for(int k = 0; k < n; k++)
    {
      w = 0;
      for(int i = 0; i < riskSet[k]; i++)
	{
	  w = w + expEta[k] * numDeath[i]/V[i];
	}
      ldot[k] = -(status[k] - w)/nrow[0];    /* OR MAYBE NOT?? */
    }

  delete [] C;
  delete [] expEta;
}

  // deathInd is indices of all deaths, totDeath is length of deathInd

double negLogLikelihoodCalc(int *riskSetInd, int *riskSet, int *numDeath, int *nDeath, int *nrow, int *ncol, double *eta, double* V, int *deathInd, int *totDeath)
{
  double first = 0;
  double second = 0;
  for(int i = 0; i < totDeath[0]; i++)
    {
      first = first + eta[deathInd[i]-1];
    }


  for(int i = 0; i < nDeath[0]; i++)
    {
      second = second + numDeath[i] * log(V[i]);
    }
  return (-first + second)/nrow[0];  /* OR MAYBE NOT? */
}



  void Solver(double *X, int* index, int *nrow, int *ncol, int *numGroup, double *beta, int *rangeGroupInd, int *groupLen, double *lambda1, double *lambda2, int *innerIter, double *thresh, double *ldot, double *nullBeta, double *gamma, double *V, int *deathInd, int *totDeath, int *riskSetInd, int *riskSet, int *numDeath, int *status, int *nDeath, double *eta, int* betaIsZero, int& groupChange, int* isActive, int* useGroup, double *step)
{
  double *theta = new double[ncol[0]];
  double *thetaNew = new double[ncol[0]];
  int startInd = 0;
  double zeroCheck = 0;
  double check = 0;
  int count = 0;
  double t = step[0];
  double diff = 1;
  double norm = 0;
  double uOp = 0;
  double Lnew = 0;
  double Lold = 0;
  double sqNormG = 0;
  double iProd = 0;
  double *etaNew = NULL;
  etaNew = new double[nrow[0]];
  double *etaNull = NULL;
  etaNull = new double[nrow[0]];
  int reset = 20;

  for(int i = 0; i < numGroup[0]; i++)
    {
      if(useGroup[i] == 1)
	{
      startInd = rangeGroupInd[i];
      
      // Setting up null gradient calc to check if group is 0
      for(int k = 0; k < nrow[0]; k++)
	{
	  etaNull[k] = eta[k];
	  for(int j = startInd; j < rangeGroupInd[i] + groupLen[i]; j++)
	    {
	      etaNull[k] = etaNull[k] - X[k + nrow[0] * j] * beta[j]; 
	    }
 	}

      // Calculating Null Gradient
      gradCalc(riskSetInd, riskSet, numDeath, status, nDeath, nrow, ncol, etaNull, ldot, V);

      double *grad = NULL;
      grad = new double[groupLen[i]];

      for(int j = 0; j < groupLen[i]; j++)
	{
	  grad[j] = 0;
	  for(int k = 0; k < nrow[0]; k++)
	    {
	      grad[j] = grad[j] + X[k + nrow[0] * (j + rangeGroupInd[i])] * ldot[k];
	    }
	  if(grad[j] < lambda1[0] && grad[j] > -lambda1[0])
	    {
	      grad[j] = 0;
	    }
	  if(grad[j] > lambda1[0])
	    {
	      grad[j] = grad[j] - lambda1[0];
	    }
	  if(grad[j] < -lambda1[0])
	    {
	      grad[j] = grad[j] + lambda1[0];
	    }
	  if(pow(grad[j],2) == pow(lambda1[0],2))
	    {
	      grad[j] = 0;
	    }

	}
      
      zeroCheck = 0;
      for(int j = 0; j < groupLen[i]; j++)
	{
	  zeroCheck = zeroCheck + pow(grad[j],2);
	}

      if(zeroCheck <= pow(lambda2[0],2)* groupLen[i])  //Or not?
	{
	  if(betaIsZero[i] == 0)
	    {
	      for(int k = 0; k < nrow[0]; k++)
		{
		  for(int j = rangeGroupInd[i]; j < rangeGroupInd[i] + groupLen[i]; j++)
		    {
		      eta[k] = eta[k] - X[k + nrow[0] * j] * beta[j];
		    }
		}
	    }
	  betaIsZero[i] = 1;
	  for(int j = 0; j < groupLen[i]; j++)
	    {
	      beta[j + rangeGroupInd[i]] = 0;
	    }
	}
      else
	{
	  if(isActive[i] == 0)
	    {
	      groupChange = 1;
	    }
	  isActive[i] = 1;

	  for(int k = 0; k < ncol[0]; k++)
	    {
	      theta[k] = beta[k];
	    }

	  betaIsZero[i] = 0;
	  double *z = NULL;
	  z = new double[groupLen[i]];
	  double *U = NULL;
	  U = new double[groupLen[i]];
	  double *G = NULL;
	  G = new double[groupLen[i]];
	  double *betaNew = NULL;
	  betaNew = new double[ncol[0]];

	  count = 0;
	  check = 1;
	  


	  while(count <= innerIter[0] && check > thresh[0])
	    {

	      count++;

	      gradCalc(riskSetInd, riskSet, numDeath, status, nDeath, nrow, ncol, eta, ldot, V);

	      for(int j = 0; j < groupLen[i]; j++)
		{		  
		  grad[j] = 0;
		  for(int k = 0; k < nrow[0]; k++)
		    {
		      grad[j] = grad[j] + X[k + nrow[0] * (j + rangeGroupInd[i])] * ldot[k];
		    }

		}
	      
	      diff = -1;
	      //	      t = 0.5;
	      	      
	      Lold = negLogLikelihoodCalc(riskSetInd, riskSet, numDeath, nDeath, nrow, ncol, eta, V, deathInd, totDeath);

	      // Back-tracking

	      while(diff < 0)
		{
		  for(int j = 0; j < groupLen[i]; j++)
		    {
		      z[j] = beta[j + rangeGroupInd[i]] - t * grad[j];
		      if(z[j] < lambda1[0] * t && z[j] > -lambda1[0] * t)
			{
			  z[j] = 0;
			}
		      if(z[j] >= lambda1[0] * t)
			{
			  z[j] = z[j] - lambda1[0] * t;
			}
		      if(z[j] <= -lambda1[0] * t)
			{
			  z[j] = z[j] + lambda1[0] * t;
			}
		    }
		  
		  norm = 0;
		  for(int j = 0; j < groupLen[i]; j++)
		    {
		      norm = norm + pow(z[j],2);
		    }
		  norm = sqrt(norm);
		  if(norm != 0){
		    uOp = (1 - lambda2[0]*sqrt(double(groupLen[i]))*t/norm);   //Or not?
		  }
		  else{uOp = 0;}

		  if(uOp < 0)
		    {
		      uOp = 0;
		    }

		  for(int j = 0; j < groupLen[i]; j++)
		    {
		      U[j] = uOp*z[j];
		      G[j] = 1/t *(beta[j + rangeGroupInd[i]] - U[j]);
		      
		    }

		  // Setting up betaNew and etaNew in direction of Grad for descent step

		  for(int j = 0; j < rangeGroupInd[i]; j++)
		    {
		      thetaNew[j] = beta[j];
		    }
		  for(int j = rangeGroupInd[i]; j < rangeGroupInd[i] + groupLen[i]; j++)
		    {
		      thetaNew[j] = beta[j] - t * G[j];
		    }
		  
		  for(int j = rangeGroupInd[i] + groupLen[i]; j < ncol[0]; j++)
		    {
		      thetaNew[j] = beta[j];
		    }
		  for(int k = 0; k < nrow[0]; k++)
		    {
		      etaNew[k] = eta[k];
			for(int j = 0; j < groupLen[i]; j++)
			  {
			    etaNew[k] = etaNew[k] - t*G[j] * X[k + nrow[0]*(rangeGroupInd[i] + j)];
			  }
		    }

		  VCalc(riskSetInd, riskSet, nDeath, nrow, ncol, etaNew, V);

		  Lnew = negLogLikelihoodCalc(riskSetInd, riskSet, numDeath, nDeath, nrow, ncol, etaNew, V, deathInd, totDeath);
		    
		  sqNormG = 0;
		  iProd = 0;

		  for(int j = 0; j < groupLen[i]; j++)
		    {
		      sqNormG = sqNormG + pow(G[j],2);
		      iProd = iProd + grad[j] * G[j];
		    }
		  
		  diff = Lold - Lnew - t * iProd + t/2 * sqNormG;
		  
		  t = t * gamma[0];
		}
	      t = t / gamma[0];

	      check = 0;
	      
	      for(int j = 0; j < groupLen[i]; j++)
		{
		  check = check + fabs(theta[j + rangeGroupInd[i]] - U[j]);
		  for(int k = 0; k < nrow[0]; k++)
		    {
		      eta[k] = eta[k] - X[k + nrow[0] * (j + rangeGroupInd[i])]*beta[j + rangeGroupInd[i]];
		    }
		  beta[j + rangeGroupInd[i]] = U[j] + count%reset/(count%reset+3) * (U[j] - theta[j + rangeGroupInd[i]]);
		  theta[j + rangeGroupInd[i]] = U[j];

		  for(int k = 0; k < nrow[0]; k++)
		    {
		      eta[k] = eta[k] + X[k + nrow[0] * (j + rangeGroupInd[i])]*beta[j + rangeGroupInd[i]];
		    }
		}
	    }
	  delete [] z;
	  delete [] U;
	  delete [] G;
	  delete [] betaNew;
	}
      delete [] grad;
	}
    }
  delete [] etaNew;
  delete [] etaNull;
  delete [] theta;
  delete [] thetaNew;
}
  


  int coxSolver(double *X, int* index, int *nrow, int *ncol, int *numGroup, int *rangeGroupInd, int *groupLen, double *lambda1, double *lambda2, double *beta, int *innerIter, int *outerIter, double *thresh, double *outerThresh, int *riskSetInd, int *riskSet, int *numDeath, int *status, int *nDeath, double *eta, double *gamma, int *deathInd, int *totDeath, int *betaIsZero, double *step)
{
  double* nullBeta = NULL;
  nullBeta = new double[ncol[0]];
  double* V = NULL;
  V = new double[nDeath[0]];
  int n = nrow[0];
  int p = ncol[0];
  double *ldot = NULL;
  ldot = new double[n];
  int groupChange = 1;
  int* isActive = NULL;
  isActive = new int[numGroup[0]];
  int* useGroup = NULL;
  useGroup = new int[numGroup[0]];
  int* tempIsActive = NULL;
  tempIsActive = new int[numGroup[0]];
  
  for(int i = 0; i < numGroup[0]; i++)
    {
      isActive[i] = 0;
      useGroup[i] = 1;
    }

  // outer most loop creating response etc...
  int outermostCounter = 0;
  double outermostCheck = 100;
  double* outerOldBeta = NULL;
  outerOldBeta = new double[p];

   while(groupChange == 1)
     {
       groupChange = 0;

       Solver(X, index, nrow, ncol, numGroup, beta, rangeGroupInd, groupLen, lambda1, lambda2, innerIter, thresh, ldot, nullBeta, gamma, V, deathInd, totDeath, riskSetInd, riskSet, numDeath, status, nDeath, eta, betaIsZero, groupChange, isActive, useGroup, step);
 
  while(outermostCounter < outerIter[0] && outermostCheck > outerThresh[0])
    {
      outermostCounter ++;
      for(int i = 0; i < p; i++)
	{
	  outerOldBeta[i] = beta[i];
	}

      for(int i = 0; i < numGroup[0]; i++)
	{
	  tempIsActive[i] = isActive[i];
	}

      Solver(X, index, nrow, ncol, numGroup, beta, rangeGroupInd, groupLen, lambda1, lambda2, innerIter, thresh, ldot, nullBeta, gamma, V, deathInd, totDeath, riskSetInd, riskSet, numDeath, status, nDeath, eta, betaIsZero, groupChange, isActive, tempIsActive, step);

	outermostCheck = 0;
      for(int i = 0; i < p; i++)
	{
	  outermostCheck = outermostCheck + fabs(outerOldBeta[i] - beta[i]);
	}
    }
     }
  
  delete [] V;
  delete [] nullBeta;
  delete [] outerOldBeta;
  delete [] ldot;
  delete [] useGroup;
  delete [] isActive;
  delete [] tempIsActive;

  return 1;
}

void Cox(int *riskSetInd, int *riskSet, int *numDeath, int *status, int *nDeath, int *nrow, int *ncol, double *beta, double *eta, double *y, double *weights)
{

  int ndeath = nDeath[0];
  double* C = NULL;
  C = new double[ndeath];
  double* V = NULL;
  V = new double[ndeath];
  double w;
  double r;
  int n = nrow[0];
  double* expEta = NULL;
  expEta = new double[n];
  double* ldot = NULL;
  ldot = new double[n];

  for(int i = 0; i < n; i++)
    {
      expEta[i] = exp(eta[i]);
    }

  for(int i = 0; i < ndeath; i++)
    {
      C[i] = 0;
	  for(int j = (riskSetInd[i]-1); j < (riskSetInd[i+1]-1); j++)
	    {
	      C[i] = C[i] + expEta[j];
	    }
    }



  V[ndeath - 1] = C[ndeath - 1];
  
  for(int i = 2; i <= ndeath; i++) 
    {
      V[ndeath - i] = V[ndeath - i + 1] + C[ndeath - i];
    }

  for(int k = 0; k < n; k++)
    {
      w = 0;
      r = 0;
      for(int i = 0; i < riskSet[k]; i++)
	{
	  w = w + expEta[k] * numDeath[i]/V[i];
	  r = r + (expEta[k] * V[i] - pow(expEta[k],2)*numDeath[i])/(pow(V[i],2));
	}
      weights[k] = r;
      ldot[k] = -(status[k] - w);
    }

  for(int i = 0; i < n; i++)
    {
      y[i] = eta[i] - ldot[i]/weights[i];
    }
  delete [] C;
  delete [] V;
  delete [] ldot;
  delete [] expEta;
}





  //////////////////////////////////

  void linGradCalc(int *nrow, double *eta, double *y, double *ldot)
{
  for(int i = 0; i < nrow[0]; i++)
    {
      ldot[i] = (eta[i] - y[i])/nrow[0];  /* OR MAYBE NOT? */
    }
}


double linNegLogLikelihoodCalc(int *nrow, double *eta, double *y)
{
  double squareSum = 0;

  for(int i = 0; i < nrow[0]; i++)
    {
      squareSum = squareSum + pow(eta[i] - y[i], 2)/2; 
    }

  return squareSum/nrow[0];   /* OR MAYBE NOT? */
}

  void linSolver(double *X, double *y, int* index, int *nrow, int *ncol, int *numGroup, double *beta, int *rangeGroupInd, int *groupLen, double *lambda1, double *lambda2, int *innerIter, double *thresh, double *ldot, double *nullBeta, double *gamma, double *eta, int* betaIsZero, int& groupChange, int* isActive, int* useGroup, double *step, int *reset, double *groupW)
{
  double *theta = new double[ncol[0]];
  int startInd = 0;
  double zeroCheck = 0;
  double check = 0;
  int count = 0;
  double t = step[0];
  double diff = 1;
  double norm = 0;
  double uOp = 0;
  double Lnew = 0;
  double Lold = 0;
  double sqNormG = 0;
  double iProd = 0;
  double *etaNew = NULL;
  etaNew = new double[nrow[0]];
  double *etaNull = NULL;
  etaNull = new double[nrow[0]];
  //  int reset = 20;

  for(int i = 0; i < numGroup[0]; i++)
    {
      if(useGroup[i] == 1)
	{
      startInd = rangeGroupInd[i];
      
      // Setting up null gradient calc to check if group is 0
      for(int k = 0; k < nrow[0]; k++)
	{
	  etaNull[k] = eta[k];
	  for(int j = startInd; j < rangeGroupInd[i] + groupLen[i]; j++)
	    {
	      etaNull[k] = etaNull[k] - X[k + nrow[0] * j] * beta[j]; 
	    }
 	}

      // Calculating Null Gradient
      linGradCalc(nrow, etaNull, y, ldot);

      double *grad = NULL;
      grad = new double[groupLen[i]];

      for(int j = 0; j < groupLen[i]; j++)
	{
	  grad[j] = 0;
	  for(int k = 0; k < nrow[0]; k++)
	    {
	      grad[j] = grad[j] + X[k + nrow[0] * (j + rangeGroupInd[i])] * ldot[k];
	    }
	  if(grad[j] < lambda1[0] && grad[j] > -lambda1[0])
	    {
	      grad[j] = 0;
	    }
	  if(grad[j] > lambda1[0])
	    {
	      grad[j] = grad[j] - lambda1[0];
	    }
	  if(grad[j] < -lambda1[0])
	    {
	      grad[j] = grad[j] + lambda1[0];
	    }
	  if(pow(grad[j],2) == pow(lambda1[0],2))
	    {
	      grad[j] = 0;
	    }

	}

      zeroCheck = 0;
      for(int j = 0; j < groupLen[i]; j++)
	{
	  zeroCheck = zeroCheck + pow(grad[j],2);
	}

      if(zeroCheck <= pow(lambda2[0]*groupW[i], 2) )  //Or not?
	{
	  if(betaIsZero[i] == 0)
	    {
	      for(int k = 0; k < nrow[0]; k++)
		{
		  for(int j = rangeGroupInd[i]; j < rangeGroupInd[i] + groupLen[i]; j++)
		    {
		      eta[k] = eta[k] - X[k + nrow[0] * j] * beta[j];
		    }
		}
	    }
	  betaIsZero[i] = 1;
	  for(int j = 0; j < groupLen[i]; j++)
	    {
	      beta[j + rangeGroupInd[i]] = 0;
	    }
	}
      else
	{
	  if(isActive[i] == 0)
	    {
	      groupChange = 1;
	    }
	  isActive[i] = 1;

	  for(int k = 0; k < ncol[0]; k++)
	    {
	      theta[k] = beta[k];
	    }

	  betaIsZero[i] = 0;
	  double *z = NULL;
	  z = new double[groupLen[i]];
	  double *U = NULL;
	  U = new double[groupLen[i]];
	  double *G = NULL;
	  G = new double[groupLen[i]];
	  double *betaNew = NULL;
	  betaNew = new double[ncol[0]];

	  count = 0;
	  check = 100000;
	  


	  while(count <= innerIter[0] && check > thresh[0])
	    {

	      count++;

	      linGradCalc(nrow, eta, y ,ldot);

	      for(int j = 0; j < groupLen[i]; j++)
		{		  
		  grad[j] = 0;
		  for(int k = 0; k < nrow[0]; k++)
		    {
		      grad[j] = grad[j] + X[k + nrow[0] * (j + rangeGroupInd[i])] * ldot[k];
		    }

		}
	      
	      diff = -1;
	      //	      t = 0.5;
	      Lold = linNegLogLikelihoodCalc(nrow, eta, y);

	      // Back-tracking

	      while(diff < 0)
		{
		  for(int j = 0; j < groupLen[i]; j++)
		    {
		      z[j] = beta[j + rangeGroupInd[i]] - t * grad[j];
		      if(z[j] < lambda1[0] * t && z[j] > -lambda1[0] * t)
			{
			  z[j] = 0;
			}
		      if(z[j] > lambda1[0] * t)
			{
			  z[j] = z[j] - lambda1[0] * t;
			}
		      if(z[j] < -lambda1[0] * t)
			{
			  z[j] = z[j] + lambda1[0] * t;
			}
		    }
		  
		  norm = 0;
		  for(int j = 0; j < groupLen[i]; j++)
		    {
		      norm = norm + pow(z[j],2);
		    }
		  norm = sqrt(norm);

		  if(norm != 0){
		    uOp = (1 - lambda2[0]*groupW[i]*t/norm);   //Or not?
		  }
		  else{uOp = 0;}

		  if(uOp < 0)
		    {
		      uOp = 0;
		    }

		  for(int j = 0; j < groupLen[i]; j++)
		    {
		      U[j] = uOp*z[j];
		      G[j] = 1/t *(beta[j + rangeGroupInd[i]] - U[j]);
		      
		    }

		  // Setting up betaNew and etaNew in direction of Grad for descent step

		  for(int k = 0; k < nrow[0]; k++)
		    {
		      etaNew[k] = eta[k];
			for(int j = 0; j < groupLen[i]; j++)
			  {
			    etaNew[k] = etaNew[k] - t*G[j] * X[k + nrow[0]*(rangeGroupInd[i] + j)];
			  }
		    }

		  Lnew = linNegLogLikelihoodCalc(nrow, etaNew, y);
		    
		  sqNormG = 0;
		  iProd = 0;

		  for(int j = 0; j < groupLen[i]; j++)
		    {
		      sqNormG = sqNormG + pow(G[j],2);
		      iProd = iProd + grad[j] * G[j];
		    }
		  
		  diff = Lold - Lnew - t * iProd + t/2 * sqNormG;
		  
		  t = t * gamma[0];
		}
	      t = t / gamma[0];

	      check = 0;
	      
	      for(int j = 0; j < groupLen[i]; j++)
		{
		  check = check + fabs(theta[j + rangeGroupInd[i]] - U[j]);
		  for(int k = 0; k < nrow[0]; k++)
		    {
		      eta[k] = eta[k] - X[k + nrow[0] * (j + rangeGroupInd[i])]*beta[j + rangeGroupInd[i]];
		    }
		  beta[j + rangeGroupInd[i]] = U[j] + count%reset[0]/(count%reset[0]+3) * (U[j] - theta[j + rangeGroupInd[i]]);
		  theta[j + rangeGroupInd[i]] = U[j];

		  for(int k = 0; k < nrow[0]; k++)
		    {
		      eta[k] = eta[k] + X[k + nrow[0] * (j + rangeGroupInd[i])]*beta[j + rangeGroupInd[i]];
		    }
		}
	    }
	  delete [] z;
	  delete [] U;
	  delete [] G;
	  delete [] betaNew;
	}
      delete [] grad;
	}
    }
  delete [] etaNew;
  delete [] etaNull;
  delete [] theta;
}

  int linNest(double *X, double* y, int* index, int *nrow, int *ncol, int *numGroup, int *rangeGroupInd, int *groupLen, double *lambda1, double *lambda2, double *beta, int *innerIter, int *outerIter, double *thresh, double *outerThresh, double *eta, double *gamma, int *betaIsZero, double *step, int *reset, double *groupW)
{
  double* prob = NULL;
  prob = new double[nrow[0]];
  double* nullBeta = NULL;
  nullBeta = new double[ncol[0]];
  int n = nrow[0];
  int p = ncol[0];
  double *ldot = NULL;
  ldot = new double[n];
  int groupChange = 1;
  int* isActive = NULL;
  isActive = new int[numGroup[0]];
  int* useGroup = NULL;
  useGroup = new int[numGroup[0]];
  int* tempIsActive = NULL;
  tempIsActive = new int[numGroup[0]];
  
  for(int i = 0; i < numGroup[0]; i++)
    {
      isActive[i] = 0;
      useGroup[i] = 1;
    }

  // outer most loop creating response etc...
  int outermostCounter = 0;
  double outermostCheck = 100000;
  double* outerOldBeta = NULL;
  outerOldBeta = new double[p];

   while(groupChange == 1)
     {
       groupChange = 0;

       linSolver(X, y, index, nrow, ncol, numGroup, beta, rangeGroupInd, groupLen, lambda1, lambda2, innerIter, thresh, ldot, nullBeta, gamma, eta, betaIsZero, groupChange, isActive, useGroup, step, reset, groupW);
 
  while(outermostCounter < outerIter[0] && outermostCheck > outerThresh[0])
    {
      outermostCounter ++;
      for(int i = 0; i < p; i++)
	{
	  outerOldBeta[i] = beta[i];
	}

      for(int i = 0; i < numGroup[0]; i++)
	{
	  tempIsActive[i] = isActive[i];
	}

      linSolver(X, y, index, nrow, ncol, numGroup, beta, rangeGroupInd, groupLen, lambda1, lambda2, innerIter, thresh, ldot, nullBeta, gamma, eta, betaIsZero, groupChange, isActive, tempIsActive, step, reset, groupW);

	outermostCheck = 0;
      for(int i = 0; i < p; i++)
	{
	  outermostCheck = outermostCheck + fabs(outerOldBeta[i] - beta[i]);
	}
    }
     }

  delete [] nullBeta;
  delete [] outerOldBeta;
  delete [] ldot;
  delete [] isActive;
  delete [] useGroup;
  delete [] tempIsActive;
  delete [] prob;

  return 1;
}



  //////////////////////////////////



  void logitGradCalc(int *nrow, double *prob, int *y, double *ldot)
  {
    for(int i = 0; i < nrow[0]; i++)
      {
	ldot[i] = (-y[i] + prob[i])/nrow[0]; /* OR MAYBE NOT? */
      }
  }
  
  void pCalc(int *nrow, double *eta, double *prob)
  {
    for(int i = 0; i < nrow[0]; i++)
      {
	prob[i] = exp(eta[i]) / (1 + exp(eta[i]));
      }
  }
  
  
  double logitNegLogLikelihoodCalc(int *nrow, double *prob, int *y)
  {
    double logLik = 0;
    
    for(int i = 0; i < nrow[0]; i++)
      {
	logLik = logLik + y[i] * log(prob[i]) + (1 - y[i]) * log(1 - prob[i]);
    }
    
    return -logLik/nrow[0];  /* OR MAYBE NOT? */
  }
  
  void betaZeroSolve(int *nrow, double *betaZero, double *eta, double *prob, double *thresh, int *innerIter, int *y)
  {
    double diff = 10;
    double num = 0;
    double denom = 0;
    int count = 0;

    while(pow(diff,2) > pow(thresh[0],2) && count < innerIter[0])
      {
	pCalc(nrow, eta, prob);
	diff = 0;

	for(int i = 0; i < nrow[0]; i++)
	  {
	    num = num + y[i] - prob[i];
	    denom = denom + prob[i] * (1 - prob[i]);
	  }
	diff = num/denom;
	betaZero[0] = betaZero[0] - diff;
	
	for(int i = 0; i < nrow[0]; i++)
	  {
	    eta[i] = eta[i] + diff;
	  }
      }
  }


  void logitSolver(double *X, int *y, int* index, int *nrow, int *ncol, int *numGroup, double *beta, int *rangeGroupInd, int *groupLen, double *lambda1, double *lambda2, int *innerIter, double *thresh, double *ldot, double *nullBeta, double *gamma, double *eta, int* betaIsZero, int& groupChange, int* isActive, int* useGroup, double *prob, double *betaZero, double *step, double *groupW)
{
  double *theta = new double[ncol[0]];
  double *thetaNew = new double[ncol[0]];
  int startInd = 0;
  double zeroCheck = 0;
  double check = 0;
  int count = 0;
  double t = step[0];
  double diff = 1;
  double norm = 0;
  double uOp = 0;
  double Lnew = 0;
  double Lold = 0;
  double sqNormG = 0;
  double iProd = 0;
  double *etaNew = NULL;
  etaNew = new double[nrow[0]];
  double *etaNull = NULL;
  etaNull = new double[nrow[0]];

  for(int i = 0; i < numGroup[0]; i++)
    {
      if(useGroup[i] == 1)
	{
      startInd = rangeGroupInd[i];
      
      // Setting up null gradient calc to check if group is 0
      for(int k = 0; k < nrow[0]; k++)
	{
	  etaNull[k] = eta[k];
	  for(int j = startInd; j < rangeGroupInd[i] + groupLen[i]; j++)
	    {
	      etaNull[k] = etaNull[k] - X[k + nrow[0] * j] * beta[j]; 
	    }
 	}

      // Calculating Null Gradient
      pCalc(nrow, etaNull, prob);
      logitGradCalc(nrow, prob, y, ldot);

      double *grad = NULL;
      grad = new double[groupLen[i]];

      for(int j = 0; j < groupLen[i]; j++)
	{
	  grad[j] = 0;
	  for(int k = 0; k < nrow[0]; k++)
	    {
	      grad[j] = grad[j] + X[k + nrow[0] * (j + rangeGroupInd[i])] * ldot[k];
	    }
	  if(grad[j] < lambda1[0] && grad[j] > -lambda1[0])
	    {
	      grad[j] = 0;
	    }
	  if(grad[j] > lambda1[0])
	    {
	      grad[j] = grad[j] - lambda1[0];
	    }
	  if(grad[j] < -lambda1[0])
	    {
	      grad[j] = grad[j] + lambda1[0];
	    }
	  if(pow(grad[j],2) == pow(lambda1[0],2))
	    {
	      grad[j] = 0;
	    }

	}
      
      zeroCheck = 0;
      for(int j = 0; j < groupLen[i]; j++)
	{
	  zeroCheck = zeroCheck + pow(grad[j],2);
	}

      if(zeroCheck <= pow(lambda2[0] * groupW[i], 2) )  //Or not?
	{
	  if(betaIsZero[i] == 0)
	    {
	      for(int k = 0; k < nrow[0]; k++)
		{
		  for(int j = rangeGroupInd[i]; j < rangeGroupInd[i] + groupLen[i]; j++)
		    {
		      eta[k] = eta[k] - X[k + nrow[0] * j] * beta[j];
		    }
		}
	    }
	  betaIsZero[i] = 1;
	  for(int j = 0; j < groupLen[i]; j++)
	    {
	      beta[j + rangeGroupInd[i]] = 0;
	    }
	}
      else
	{
	  if(isActive[i] == 0)
	    {
	      groupChange = 1;
	    }
	  isActive[i] = 1;

	  for(int k = 0; k < ncol[0]; k++)
	    {
	      theta[k] = beta[k];
	    }

	  betaIsZero[i] = 0;
	  double *z = NULL;
	  z = new double[groupLen[i]];
	  double *U = NULL;
	  U = new double[groupLen[i]];
	  double *G = NULL;
	  G = new double[groupLen[i]];
	  double *betaNew = NULL;
	  betaNew = new double[ncol[0]];

	  count = 0;
	  check = 1000000;
	  


	  while(count <= innerIter[0] && check > thresh[0])
	    {

	      count++;

	      pCalc(nrow, eta, prob);
	      logitGradCalc(nrow, prob, y ,ldot);

	      for(int j = 0; j < groupLen[i]; j++)
		{		  
		  grad[j] = 0;
		  for(int k = 0; k < nrow[0]; k++)
		    {
		      grad[j] = grad[j] + X[k + nrow[0] * (j + rangeGroupInd[i])] * ldot[k];
		    }

		}
	      
	      diff = -1;
	      //	      t = 0.5;
	      pCalc(nrow, eta, prob);
	      Lold = logitNegLogLikelihoodCalc(nrow, prob, y);

	      // Back-tracking

	      while(diff < 0)
		{
		  for(int j = 0; j < groupLen[i]; j++)
		    {
		      z[j] = beta[j + rangeGroupInd[i]] - t * grad[j];
		      if(z[j] < lambda1[0] * t && z[j] > -lambda1[0] * t)
			{
			  z[j] = 0;
			}
		      if(z[j] > lambda1[0] * t)
			{
			  z[j] = z[j] - lambda1[0] * t;
			}
		      if(z[j] < -lambda1[0] * t)
			{
			  z[j] = z[j] + lambda1[0] * t;
			}
		    }
		  
		  norm = 0;
		  for(int j = 0; j < groupLen[i]; j++)
		    {
		      norm = norm + pow(z[j],2);
		    }
		  norm = sqrt(norm);

		  if(norm != 0){
		    uOp = (1 - lambda2[0]*groupW[i]*t/norm);   //Or not?
		  }
		  else{uOp = 0;}

		  if(uOp < 0)
		    {
		      uOp = 0;
		    }

		  for(int j = 0; j < groupLen[i]; j++)
		    {
		      U[j] = uOp*z[j];
		      G[j] = 1/t *(beta[j + rangeGroupInd[i]] - U[j]);
		      
		    }

		  // Setting up betaNew and etaNew in direction of Grad for descent step

		  for(int j = 0; j < rangeGroupInd[i]; j++)
		    {
		      thetaNew[j] = beta[j];
		    }
		  for(int j = rangeGroupInd[i]; j < rangeGroupInd[i] + groupLen[i]; j++)
		    {
		      thetaNew[j] = beta[j] - t * G[j];
		    }
		  
		  for(int j = rangeGroupInd[i] + groupLen[i]; j < ncol[0]; j++)
		    {
		      thetaNew[j] = beta[j];
		    }
		  for(int k = 0; k < nrow[0]; k++)
		    {
		      etaNew[k] = eta[k];
			for(int j = 0; j < groupLen[i]; j++)
			  {
			    etaNew[k] = etaNew[k] - t*G[j] * X[k + nrow[0]*(rangeGroupInd[i] + j)];
			  }
		    }

		  pCalc(nrow, etaNew, prob);
		  Lnew = logitNegLogLikelihoodCalc(nrow, prob, y);
		    
		  sqNormG = 0;
		  iProd = 0;

		  for(int j = 0; j < groupLen[i]; j++)
		    {
		      sqNormG = sqNormG + pow(G[j],2);
		      iProd = iProd + grad[j] * G[j];
		    }
		  
		  diff = Lold - Lnew - t * iProd + t/2 * sqNormG;
		  
		  t = t * gamma[0];
		}
	      t = t / gamma[0];

	      check = 0;
	      
	      for(int j = 0; j < groupLen[i]; j++)
		{
		  check = check + fabs(theta[j + rangeGroupInd[i]] - U[j]);
		  for(int k = 0; k < nrow[0]; k++)
		    {
		      eta[k] = eta[k] - X[k + nrow[0] * (j + rangeGroupInd[i])]*beta[j + rangeGroupInd[i]];
		    }
		  beta[j + rangeGroupInd[i]] = U[j] + count/(count+3) * (U[j] - theta[j + rangeGroupInd[i]]);
		  theta[j + rangeGroupInd[i]] = U[j];

		  for(int k = 0; k < nrow[0]; k++)
		    {
		      eta[k] = eta[k] + X[k + nrow[0] * (j + rangeGroupInd[i])]*beta[j + rangeGroupInd[i]];
		    }
		}
	    }
	  delete [] z;
	  delete [] U;
	  delete [] G;
	  delete [] betaNew;
	}
      delete [] grad;
	}
    }
  betaZeroSolve(nrow, betaZero, eta, prob, thresh, innerIter, y);
    
  delete [] etaNew;
  delete [] etaNull;
  delete [] theta;
  delete [] thetaNew;
}
  


  int logitNest(double *X, int* y, int* index, int *nrow, int *ncol, int *numGroup, int *rangeGroupInd, int *groupLen, double *lambda1, double *lambda2, double *beta, int *innerIter, int *outerIter, double *thresh, double *outerThresh, double *eta, double *gamma, int *betaIsZero, double* betaZero, double *step, double *groupW)
{
  double oldBetaZero = betaZero[0];
  double* prob = NULL;
  prob = new double[nrow[0]];
  double* nullBeta = NULL;
  nullBeta = new double[ncol[0]];
  int n = nrow[0];
  int p = ncol[0];
  double *ldot = NULL;
  ldot = new double[n];
  int groupChange = 1;
  int* isActive = NULL;
  isActive = new int[numGroup[0]];
  int* useGroup = NULL;
  useGroup = new int[numGroup[0]];
  int* tempIsActive = NULL;
  tempIsActive = new int[numGroup[0]];
  
  for(int i = 0; i < numGroup[0]; i++)
    {
      isActive[i] = 0;
      useGroup[i] = 1;
    }

  // outer most loop creating response etc...
  int outermostCounter = 0;
  double outermostCheck = 1000000;
  double* outerOldBeta = NULL;
  outerOldBeta = new double[p];

  while(groupChange == 1)
  {
      groupChange = 0;
      logitSolver(X, y, index, nrow, ncol, numGroup, beta, rangeGroupInd, groupLen, lambda1, lambda2, innerIter, thresh, ldot, nullBeta, gamma, eta, betaIsZero, groupChange, isActive, useGroup, prob, betaZero, step, groupW);
 
      while(outermostCounter < outerIter[0] && outermostCheck > outerThresh[0])
      {
          outermostCounter ++;
          for(int i = 0; i < p; i++)
	        {
	            outerOldBeta[i] = beta[i];
	        }
          oldBetaZero = betaZero[0];

          for(int i = 0; i < numGroup[0]; i++)
	        {
	            tempIsActive[i] = isActive[i];
	        }
          logitSolver(X, y, index, nrow, ncol, numGroup, beta, rangeGroupInd, groupLen, lambda1, lambda2, innerIter, thresh, ldot, nullBeta, gamma, eta, betaIsZero, groupChange, tempIsActive, isActive, prob, betaZero, step, groupW);
	        outermostCheck = 0;
          for(int i = 0; i < p; i++)
	        {
	            outermostCheck = outermostCheck + fabs(outerOldBeta[i] - beta[i]);
	        }
          outermostCheck = outermostCheck + fabs(oldBetaZero - betaZero[0]);
    }
  }

  delete [] nullBeta;
  delete [] outerOldBeta;
  delete [] ldot;
  delete [] isActive;
  delete [] useGroup;
  delete [] tempIsActive;
  delete [] prob;

  return 1;
}
}





