      parameter (ncomp=0)        ! to turn off compare3 diagnostics
!     parameter (ncomp=1)        ! to turn on  compare3 diagnostics
      parameter (ifull_ncomp=il**ncomp * jl**ncomp)
      common/comhtg/dqc(ifull_ncomp,kl**ncomp)
     .             ,dql(ifull_ncomp,kl**ncomp)
     .             ,dqv(ifull_ncomp,kl**ncomp)
     .             ,dtc(ifull_ncomp,kl**ncomp)
     .             ,dtl(ifull_ncomp,kl**ncomp)
     .             ,dtv(ifull_ncomp,kl**ncomp)
     .            ,dtlw(ifull_ncomp,kl**ncomp)
     .            ,dtsw(ifull_ncomp,kl**ncomp)
     .            , duv(ifull_ncomp,kl**ncomp)
     .            , dvv(ifull_ncomp,kl**ncomp)
     .            , egm(ifull_ncomp)
     .            , fgm(ifull_ncomp)
