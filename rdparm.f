c $Log$
c Revision 1.1  2003/08/13 01:24:20  dix043
c Initial revision
c
c Revision 1.1  1996/10/17  05:14:05  mrd
c Initial revision
c
c Revision 1.5  1994/12/07  23:53:36  mrd
c Explicitly declare all variables in include files.
c
c Revision 1.4  1992/12/09  14:42:34  ldr
c Replaced all n's with nl's for compatibility with ogcm.
c
c Revision 1.3  92/08/24  14:51:28  mrd
c Updated (12 band) solar radiation code
c 
c Revision 1.2  92/05/11  15:15:34  ldr
c Removed include newmpar.h, placing it in source files instead.
c 
c Revision 1.1  92/04/15  11:13:38  mrd
c Initial revision
c 
c     parameter settings for the longwave and shortwave radiation code: 
c          imax   =  no. points along the lat. circle used in calcs.
c          l      =  no. vertical levels (also layers) in model 
c***note: the user normally will modify only the imax and l parameters
c          nblw   =  no. freq. bands for approx computations. see 
c                      bandta for definition
c          nblx   =  no. freq bands for approx cts computations 
c          nbly   =  no. freq. bands for exact cts computations. see
c                      bdcomb for definition
c          inlte  =  no. levels used for nlte calcs.
c          nnlte  =  index no. of freq. band in nlte calcs. 
c          nb,ko2 are shortwave parameters; other quantities are derived
c                    from the above parameters. 
c   Source file must also include newmpar.h, before RDPARM.f
      
      integer l, imax, nblw, nblx, nbly, nblm, lp1, lp2, lp3,
     &        lm1, lm2, lm3, ll, llp1, llp2, llp3, llm1, llm2,
     &        llm3, lp1m, lp1m1, lp1v, lp121, ll3p, nb, inlte,
     &        inltep, nnlte, lp1i, llp1i, ll3pi, nb1, ko2, ko21, ko2m

      parameter (l=kl)
      parameter (imax=il)
      parameter (nblw=163,nblx=47,nbly=15)
      parameter (nblm=nbly-1) 
      parameter (lp1=l+1,lp2=l+2,lp3=l+3) 
      parameter (lm1=l-1,lm2=l-2,lm3=l-3) 
      parameter (ll=2*l,llp1=ll+1,llp2=ll+2,llp3=ll+3)
      parameter (llm1=ll-1,llm2=ll-2,llm3=ll-3) 
      parameter (lp1m=lp1*lp1,lp1m1=lp1m-1) 
      parameter (lp1v=lp1*(1+2*l/2))
      parameter (lp121=lp1*nbly)
      parameter (ll3p=3*l+2)
      parameter (nb=12)
      parameter (inlte=3,inltep=inlte+1,nnlte=56) 
      parameter (lp1i=imax*lp1,llp1i=imax*llp1,ll3pi=imax*ll3p) 
      parameter (nb1=nb-1)
      parameter (ko2=12)
      parameter (ko21=ko2+1,ko2m=ko2-1) 

