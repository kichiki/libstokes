/*************************************************************************/
/* Beenakker's Idea of Ewald summation technique for RP tensor           */
/*  ３次元配列 ３次元周期境界版                                          */
/*  無次元版 non-dimension formulation                                   */
/*  from cal_f_p3.c on '93/ 9/18                                         */
/*  Reviced on '93/ 9/20  Original (but all removed !!)                  */
/*                  9/21  rebuild (almost same above !!)                 */
/*                  9/22  ludcmp() で determinant を計算する             */
/*                 10/08  calc_mob_ewald_3fのミスチェックのためシンプル化*/
/*                 10/09  波数格子和のミスの修正                         */
/*             '94/06/09  calc_steady_vel_3f_fix の追加                  */
/*             '95/02/24  2D との比較、シンプル化                        */
/*                 02/27  debug の結果シンプル化                         */
/*                 03/02  calc_steady_vel() を削除( _fix に含まれる)     */
/*             '96/04/23  無次元化                                       */
/*                 05/08 式の演算順序を整える                            */
/*                       グローバル変数を整理                            */
/*************************************************************************/
#define CALF3

#include <math.h>
#include <stdio.h> /* for printf() */
#include <stdlib.h> /* for exit() */
#include "fun3.h" /* Wmat が定義されている */
#include "lub3.h" /* calc_lub_3f() */

#include "calf3.h" /* 本ファイルのヘッダー */

void calc_mob_ewald_3f(int n,double x[],double mat[]){
/* n は固定粒子も含んだ粒子数 */
  /* Matrix work area */
  extern int 
    Indx[];
  extern double 
    Wmat1[];

  extern int 
    pcellx,pcelly,pcellz,/* 周期セルの取る個数 */
    kmaxx,kmaxy,kmaxz; /* 波数ベクトルの取る個数 */

  extern double 
    zeta,zeta2,zaspi,za2,
    pi2,
    pi6vol,
    lx,ly,lz; /* cell size */

  double 
    d,
    m1i,m1n,
    m1xx,m1yy,m1zz,
    m1xy,m1yz,m1zx,
    ex,ey,ez,
    exx,eyy,ezz,
    exy,eyz,ezx,
    xx,yy,zz,rr,
    zr,zr2,
    s,s2,
    rlx,rly,rlz;

  int 
    n3,
    i,j,
    i1,i2,i3,j1,j2,j3,
    k,m1,m2,m3;

  double 
    k1,k2,k3,kk,k4z,
    cfactor,
    mob,mob11,mob22,mob33,
    mob12,mob23,mob31;

  n3=n*3;
/* 行列の消去 clear matrix */
  for(i=0;i<n3*n3;i++){
    Wmat1[i]=0.0;
  }
/* diagonal part ( self part ) 対角部分 */
  for(i=0;i<n3;i++){
    Wmat1[i*n3+i] = 1.0 - zaspi*(6.0 - 40.0/3.0*za2);
  }

/* ----- first Ewald part ( real space ) 実空間の格子和 -------*/
  for(i=0;i<n;i++){
    i1=i*3;
    i2=i*3+1;
    i3=i*3+2;
    for(j=0;j<n;j++){
      j1=j*3;
      j2=j*3+1;
      j3=j*3+2;

      for(m1=-pcellx;m1<=pcellx;m1++){
	rlx=lx*(double)m1;
	for(m2=-pcelly;m2<=pcelly;m2++){
	  rly=ly*(double)m2;
	  for(m3=-pcellz;m3<=pcellz;m3++){
	    rlz=lz*(double)m3;
  
	    xx = x[j1]-x[i1]+rlx;
	    yy = x[j2]-x[i2]+rly;
	    zz = x[j3]-x[i3]+rlz;
	    rr = sqrt(xx*xx+yy*yy+zz*zz);

	    if(rr>0.0){
	      zr = zeta*rr;  zr2 = zr*zr;
	      s  = rr; s2 = s*s;

	      ex = xx/rr; ey = yy/rr; ez = zz/rr;
	      exx=ex*ex; eyy=ey*ey; ezz=ez*ez;
	      exy=ex*ey; eyz=ey*ez; ezx=ez*ex;

	      m1i = ( 0.75 + 0.5/s2 )/s * erfc(zr)
		+ ( (1.0+zr2*(14.0+4.0*zr2*(-5.0+zr2)))/s2
		   -4.5+3.0*zr2)*zaspi * exp( -zr2 );
	      m1n = ( 0.75 - 1.5/s2 )/s * erfc(zr)
		+ ( (-3.0+zr2*(-2.0+4.0*zr2*(4.0-zr2)))/s2
		   +1.5-3.0*zr2)*zaspi * exp( -zr2 );

	      m1xx =( m1i + m1n * exx );m1yy =( m1i + m1n * eyy );
	      m1zz =( m1i + m1n * ezz );
	      m1xy = m1n * exy;m1yz = m1n * eyz;m1zx = m1n * ezx;

	      Wmat1[i1 *n3 +j1] += m1xx;
	      Wmat1[i2 *n3 +j2] += m1yy;
	      Wmat1[i3 *n3 +j3] += m1zz;

	      Wmat1[i1 *n3 +j2] += m1xy;Wmat1[i2 *n3 +j1] += m1xy;
	      Wmat1[i2 *n3 +j3] += m1yz;Wmat1[i3 *n3 +j2] += m1yz;
	      Wmat1[i3 *n3 +j1] += m1zx;Wmat1[i1 *n3 +j3] += m1zx;
	    }
	  }
	}
      }
    }
  }

/* Second Ewald part 逆空間での格子和 ( reciprocal space ) */
  for(m1=-kmaxx;m1<=kmaxx;m1++){
    k1=pi2*(double)m1/lx;
    for(m2=-kmaxy;m2<=kmaxy;m2++){
      k2=pi2*(double)m2/ly;
      for(m3=-kmaxz;m3<=kmaxz;m3++){
	k3=pi2*(double)m3/lz;
	if(m1!=0 || m2!=0 || m3!=0){
	  kk=k1*k1+k2*k2+k3*k3;
	  k4z = kk/4.0/zeta2;
	  mob = pi6vol*(1.0-kk/3.0)*(1.0+k4z*(1.0+2.0*k4z))/kk*exp(-k4z);

	  mob11= (1.0-k1*k1/kk)*mob;
	  mob22= (1.0-k2*k2/kk)*mob;
	  mob33= (1.0-k3*k3/kk)*mob;
	  mob12=-k1*k2/kk*mob;
	  mob23=-k2*k3/kk*mob;
	  mob31=-k3*k1/kk*mob;

	  for(i=0;i<n;i++){
	    i1=i*3;
	    i2=i*3+1;
	    i3=i*3+2;
	    for(j=0;j<n;j++){
	      j1=j*3;
	      j2=j*3+1;
	      j3=j*3+2;
	      xx=x[j1]-x[i1];
	      yy=x[j2]-x[i2];
	      zz=x[j3]-x[i3];

	      cfactor=cos(+k1*xx+k2*yy+k3*zz);

	      Wmat1[i1 *n3 +j1] += cfactor*mob11;
	      Wmat1[i2 *n3 +j2] += cfactor*mob22;
	      Wmat1[i3 *n3 +j3] += cfactor*mob33;

	      Wmat1[i1 *n3 +j2] += cfactor*mob12;
	      Wmat1[i2 *n3 +j1] += cfactor*mob12;
	      Wmat1[i2 *n3 +j3] += cfactor*mob23;
	      Wmat1[i3 *n3 +j2] += cfactor*mob23;
	      Wmat1[i3 *n3 +j1] += cfactor*mob31;
	      Wmat1[i1 *n3 +j3] += cfactor*mob31;
	    }
	  }
	}
      }
    }
  }

  ludcmp(Wmat1,n3,Indx,&d);
  if(d<0.0){
    fprintf(stderr,"Negative Matrix in F ewald(det = %e)\n",d);
    /*exit(0);*/
  }
  luinv(Wmat1,n3,Indx,mat);
}

void set_Uniform_3f(int n,double m[],double x[],double f[],
int neighbor[27][N][NB]){
  extern double 
    Uniform,
    Wmat1[],Wmat2[];
  extern int 
    Indx[];
  double 
    d,
    v,vtmp;
  int 
    i,j,
    n3;

  n3 = n*3;

  calc_lub_3f(n,n,x,Wmat2,neighbor);
  for(i=0;i<n3*n3;i++){
    Wmat2[i]+=m[i];
  }

  ludcmp(Wmat2,n3,Indx,&d);
  if(d<0.0){
    fprintf(stderr,"Negative Matrix in Calc. Uniform Flow Velocity!\n");
    printf("Negative Matrix in Calc. Uniform Flow Velocity!\n");
    /*exit(1);*/
  }
  luinv(Wmat2,n3,Indx,Wmat1);
  
  vtmp=0.0;
  for(i=0;i<n;i++){
    v=0.0;
    for(j=0;j<n3;j++){
      v+=Wmat1[(i*3+2)*n3+j]*f[j];
    }
    vtmp+=v;
  }
  Uniform=-vtmp/(double)n;
}

void calc_fluid_force_3f(int n,int nm,double mat[],
double x[],double f[],double v[],int neighbor[27][N][NB]){
  extern double 
    Uniform,
    Wmat1[];
  extern int 
    Indx[];
  int 
    i,j,n3,nm3;

  n3 = n*3;
  nm3= nm*3;

  calc_lub_3f(n,nm,x,Wmat1,neighbor);

  for(i=0;i<n3*n3;i++){
    Wmat1[i]+=mat[i];
  }
  for(i=0;i<nm3;i++){
    f[i]=0.0;
    for(j=0;j<nm;j++){
      f[i]+=Wmat1[i*n3+j*3  ]* v[j*3  ]
           +Wmat1[i*n3+j*3+1]* v[j*3+1]
           +Wmat1[i*n3+j*3+2]*(v[j*3+2]-Uniform);
    }
  }
}

void calc_steady_vel_3f_fix(int n,int nm,double m[],double x[],double f[],
double v[],int neighbor[27][N][NB]){
  extern double 
    Uniform;
  extern int 
    Indx[];
  extern double 
    Wmat1[],Wmat2[];
  double 
    d,ff[3*N];
  int 
    i,j,
    n3,nm3;

  n3 = n*3;
  nm3= nm*3;

  calc_lub_3f(n,n,x,Wmat2,neighbor);

  for(i=0;i<n3*n3;i++){
    Wmat2[i]+=m[i];
  }
  for(i=0;i<nm3;i++){
    for(j=0;j<nm3;j++){
      Wmat1[i*nm3+j]=Wmat2[i*n3+j];
    }
  }
  for(i=0;i<nm3;i++){
    ff[i]=f[i];
    for(j=0;j<(n-nm);j++){
      ff[i]+=Wmat2[i*n3+(j+nm)*3+2]*Uniform;
    }
  }

  ludcmp(Wmat1,nm3,Indx,&d);
  if(d<0.0){
    fprintf(stderr,"Negative Matrix in Calc. Steady Velocity!!\n");
    printf("Negative Matrix in Calc. Steady Velocity!!\n");
    /*exit(1);*/
  }
  luinv(Wmat1,nm3,Indx,Wmat2);

  for(i=0;i<nm3;i++){
    v[i]=0.0;
    for(j=0;j<nm3;j++){
      v[i]+=Wmat2[i*nm3+j]*ff[j];
    }
  }
  for(i=0;i<nm;i++){
    v[i*3+2]+=Uniform;
  }
}

/*=============================================================================
 from Numerical Recipes
 =============================================================================*/
void ludcmp(double a[],int n,int indx[],double *d){
/* modefied on '93/09/22 -- return 'd' as determinant */
  double tiny=1.0e-20;
  double vv[3*N],
  aamax,sum,dum,tmp;
  int i,j,k,
  imax;
  
  *d=1.0;
  for(i=1;i<=n;i++){
    aamax=0.0;
    for(j=1;j<=n;j++){
      if((tmp=fabs(a[(j-1)*n+i-1]))>aamax){
	aamax=tmp;
      }
    }
    if(aamax==0.0){
      printf("Singular matrix\n");
      exit(1);
    }
    vv[i-1]=1.0/aamax;
  }
  for(j=1;j<=n;j++){
    for(i=1;i<=j-1;i++){
      sum=a[(j-1)*n+i-1];
      for(k=1;k<=i-1;k++){
	sum-=a[(k-1)*n+i-1]*a[(j-1)*n+k-1];
      }
      a[(j-1)*n+i-1]=sum;
    }
    aamax=0.0;
    for(i=j;i<=n;i++){
      sum=a[(j-1)*n+i-1];
      for(k=1;k<=j-1;k++){
	sum-=a[(k-1)*n+i-1]*a[(j-1)*n+k-1];
      }
      a[(j-1)*n+i-1]=sum;
      dum=vv[i-1]*fabs(sum);
      if(dum>=aamax){
	imax=i;
	aamax=dum;
      }
    }
    if(j!=imax){
      for(k=1;k<=n;k++){
	dum=a[(k-1)*n+imax-1];
	a[(k-1)*n+imax-1]=a[(k-1)*n+j-1];
	a[(k-1)*n+j-1]=dum;
      }
      *d=-(*d);
      vv[imax-1]=vv[j-1];
    }
    indx[j-1]=imax;
    if(a[(j-1)*n+j-1]==0.0){
      a[(j-1)*n+j-1]=tiny;
    }
    if(j!=n){
      dum=1.0/a[(j-1)*n+j-1];
      for(i=j+1;i<=n;i++){
	a[(j-1)*n+i-1]=a[(j-1)*n+i-1]*dum;
      }
    }
  }
  /* calc. determinant */
  for(i=0;i<n;i++){
    *d *= a[i*n+i];
  }
}

void lubksb(double a[],int n,int indx[],double b[]){
  int i,j,
  ii,ll;
  double sum;
  
  ii=0;
  for(i=1;i<=n;i++){
    ll=indx[i-1];
    sum=b[ll-1];
    b[ll-1]=b[i-1];
    if(ii!=0){
      for(j=ii;j<=i-1;j++){
	sum-=a[(j-1)*n+i-1]*b[j-1];
      }
    }else if(sum!=0.0){
      ii=i;
    }
    b[i-1]=sum;
  }
  for(i=n;i>=1;i--){
    sum=b[i-1];
    if(i<n){
      for(j=i+1;j<=n;j++){
	sum-=a[(j-1)*n+i-1]*b[j-1];
      }
    }
    b[i-1]=sum/a[(i-1)*n+i-1];
  }
}

void luinv(double a[],int n,int indx[],double y[]){
/* a[] must be decomposed by subroutine 'ludcmp' */
  int i,j;
  
  for(i=0;i<n;i++){
    for(j=0;j<n;j++){
      y[i*n+j]=0.0;
    }
    y[i*n+i]=1.0;
  }
  
  for(j=0;j<n;j++){
    lubksb(a,n,indx,&y[j*n]);
  }
}
