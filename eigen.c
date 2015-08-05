/*  RAxML-VI-HPC (version 2.2) a program for sequential and parallel estimation of phylogenetic trees 
 *  Copyright August 2006 by Alexandros Stamatakis
 *
 *  Partially derived from
 *  fastDNAml, a program for estimation of phylogenetic trees from sequences by Gary J. Olsen
 *  
 *  and 
 *
 *  Programs of the PHYLIP package by Joe Felsenstein.
 *
 *  This program is free software; you may redistribute it and/or modify its
 *  under the terms of the GNU General Public License as published by the Free
 *  Software Foundation; either version 2 of the License, or (at your option)
 *  any later version.
 *
 *  This program is distributed in the hope that it will be useful, but
 *  WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 *  or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
 *  for more details.
 * 
 *
 *  For any other enquiries send an Email to Alexandros Stamatakis
 *  Alexandros.Stamatakis@epfl.ch
 *
 *  When publishing work that is based on the results from RAxML-VI-HPC please cite:
 *
 *  Alexandros Stamatakis:"RAxML-VI-HPC: maximum likelihood-based phylogenetic analyses with thousands 
 *  of taxa and mixed models". 
 *  Bioinformatics 2006; doi: 10.1093/bioinformatics/btl446
 */

#include <math.h>
#include <assert.h>


static void mytred2(double **a, const int n, double *d, double *e)
{
  int     l, k, j, i;
  double  scale, hh, h, g, f; 
 
  for (i = n; i > 1; i--)
    {
      l = i - 1;
      h = 0.0;
      scale = 0.0;
      
      if (l > 1)
        {
	  for (k = 1; k <= l; k++)
	    scale += fabs(a[k - 1][i - 1]);
	  if (scale == 0.0)
	    e[i - 1] = a[l - 1][i - 1];
	  else
            {
	      for (k = 1; k <= l; k++)
                {
		  a[k - 1][i - 1] /= scale;
		  h += a[k - 1][i - 1] * a[k - 1][i - 1];
                }
	      f = a[l - 1][i - 1];
	      g = ((f > 0) ? -sqrt(h) : sqrt(h)); /* diff */
	      e[i - 1] = scale * g;
	      h -= f * g;
	      a[l - 1][i - 1] = f - g;
	      f = 0.0;
	      for (j = 1; j <= l; j++)
		{
		  a[i - 1][j - 1] = a[j - 1][i - 1] / h;
		  g = 0.0;
		  for (k = 1; k <= j; k++)
		    g += a[k - 1][j - 1] * a[k - 1][i - 1];
		  for (k = j + 1; k <= l; k++)
		    g += a[j - 1][k - 1] * a[k - 1][i - 1];
		  e[j - 1] = g / h;
		  f += e[j - 1] * a[j - 1][i - 1];
		}
	      hh = f / (h + h);
	      for (j = 1; j <= l; j++)
		{
		  f = a[j - 1][i - 1];
		  g = e[j - 1] - hh * f;
		  e[j - 1] = g;
		  for (k = 1; k <= j; k++)
		    a[k - 1][j - 1] -= (f * e[k - 1] + g * a[k - 1][i - 1]);
                }
	    }
	} 
      else
	e[i - 1] = a[l - 1][i - 1];
      d[i - 1] = h;
    }
  d[0] = 0.0;
  e[0] = 0.0;
  
  for (i = 1; i <= n; i++)
    {
      l = i - 1;
      if (d[i - 1] != 0.0)
	{
	  for (j = 1; j <= l; j++)
            {
                g = 0.0;
                for (k = 1; k <= l; k++)
		  g += a[k - 1][i - 1] * a[j - 1][k - 1];
                for(k = 1; k <= l; k++)
		  a[j - 1][k - 1] -= g * a[i - 1][k - 1];
            }
        }
      d[i - 1] = a[i - 1][i - 1];
      a[i - 1][i - 1] = 1.0;
      for (j = 1; j <= l; j++)
	a[i - 1][j - 1] = a[j - 1][i - 1] = 0.0;
    }
 
 
}
/*#define MYSIGN(a,b) ((b)<0 ? -fabs(a) : fabs(a))*/

static int mytqli(double *d, double *e, const int n, double **z)
{
  int     m, l, iter, i, k;
  double  s, r, p, g, f, dd, c, b;
   
  for (i = 2; i <= n; i++)
    e[i - 2] = e[i - 1];

  e[n - 1] = 0.0;

  for (l = 1; l <= n; l++)
    {
      iter = 0;
      do
	{
	  for (m = l; m <= n - 1; m++)
            {
	      dd = fabs(d[m - 1]) + fabs(d[m]);
	      if (fabs(e[m - 1]) + dd == dd)
		break;
	    }

	  if (m != l)
           {
	     assert(iter < 30);
	     
	     g = (d[l] - d[l - 1]) / (2.0 * e[l - 1]);
	     r = sqrt((g * g) + 1.0);
	     g = d[m - 1] - d[l - 1] + e[l - 1] / (g + ((g < 0)?-fabs(r):fabs(r)));/*MYSIGN(r, g));*/
	     s = c = 1.0;
	     p = 0.0;

	     for (i = m - 1; i >= l; i--)
               {
		 f = s * e[i - 1];
		 b = c * e[i - 1];
		 if (fabs(f) >= fabs(g))
		   {
		     c = g / f;
		     r = sqrt((c * c) + 1.0);
		     e[i] = f * r;
		     c *= (s = 1.0 / r);
		   } 
		 else
		   {
		     s = f / g;
		     r = sqrt((s * s) + 1.0);
		     e[i] = g * r;
		     s *= (c = 1.0 / r);
		   }
		 g = d[i] - p;
		 r = (d[i - 1] - g) * s + 2.0 * c * b;
		 p = s * r;
		 d[i] = g + p;
		 g = c * r - b;
		 for (k = 1; k <= n; k++)
		   {
		     f = z[i][k-1];
		     z[i][k-1] = s * z[i - 1][k - 1] + c * f;
		     z[i - 1][k - 1] = c * z[i - 1][k - 1] - s * f;
		   }
	       }

	     d[l - 1] = d[l - 1] - p;
	     e[l - 1] = g;
	     e[m - 1] = 0.0;
	   }
	} 
      while (m != l);
    }

    
 
    return (1);
 }


void makeEigen(double **_a, const int n, double *d, double *e)
{
  mytred2(_a, n, d, e);
  mytqli(d, e, n, _a);
}


