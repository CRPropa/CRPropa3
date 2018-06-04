      subroutine channels
c
c +---------------------------------------------------------------------
c | Author: Marieke Duijvestijn and Arjan Koning
c | Date  : April 1, 2014
c | Task  : Exclusive reaction channels
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      character*1 apos
      integer     type,nen,Zend,Nend,Zix,Nix,Aix,inend,ipend,idend,
     +            itend,ihend,iaend,in,ip,id,it,ih,ia,Ztot,Ntot,npart,
     +            ident,i,i2,nexout,identorg(0:numpar),idd,idorg,Zcomp,
     +            Ncomp,nex,type2,NL,idc
      real        specexcl(0:numchantot,0:numpar,0:numex+1,0:numen),
     +            term1,term2,term3,term,fissum,emissum,Eaveragesum,frac
c
c ************************ Initialization ******************************
c
c channelsum: sum over exclusive channel cross sections
c
      channelsum=0.
c
c Initially, all the flux is in the initial compound state.
c
c flaginitpop: flag for initial population distribution
c xsexcl     : exclusive cross section per excitation energy
c xsinitpop  : initial population cross section
c gamexcl    : exclusive gamma cross section per excitation energy
c maxex      : maximum excitation energy bin for compound nucleus
c xsreacinc  : reaction cross section for incident channel
c
      if (flaginitpop) then
        xsexcl(0,maxex(0,0)+1)=xsinitpop
      else
        xsexcl(0,maxex(0,0)+1)=xsreacinc
      endif
      gamexcl(0,maxex(0,0)+1)=0.
c
c ********** Construction of exclusive channel cross sections **********
c
c 1. Loop over all residual nuclei, starting with the first residual
c    nucleus, and then according to decreasing Z and N.
c
c idnum   : counter for exclusive channel
c Zend    : maximal charge number
c Nend    : maximal neutron number
c numZchan: maximal number of outgoing proton units in individual
c           channel description
c numNchan: maximal number of outgoing neutron units in individual
c           channel description
c numZ    : maximal number of protons away from the initial compound
c           nucleus
c Zinit   : charge number of initial compound nucleus
c Ninit   : neutron number of initial compound nucleus
c numN    : maximal number of neutron away from the initial compound
c           nucleus
c Zix     : charge number index for residual nucleus
c Nix     : neutron number index for residual nucleus
c
c Each idnum represents a different exclusive channel.
c
      idnum=-1
      Zend=min(numZchan,numZ)
      Zend=min(Zend,Zinit)
      Nend=min(numNchan,numN)
      Nend=min(Nend,Ninit)
      do 110 Zix=0,Zend
        do 110 Nix=0,Nend
c
c 2. To minimize the loops, the maximal possible number of each particle
c    type is determined, given Zix and Nix. E.g. for Zix=1, Nix=2 there
c    can be at most one proton, 2 neutrons, one deuteron or 1 triton in
c    the exit channel. These maximal numbers are determined by simple
c    formulae.
c
c Aix       : mass number index for residual nucleus
c inend,..  : help variables
c parinclude: logical to include outgoing particle
c maxchannel: maximal number of outgoing particles in individual
c             channel description (e.g. this is 3 for (n,2np))
c numin,....: maximal number of ejectile in channel description
c
          Aix=Zix+Nix
          inend=0
          ipend=0
          idend=0
          itend=0
          ihend=0
          iaend=0
          if (Aix.ne.0) then
            if (parinclude(1)) inend=min(Nix,maxchannel)
            if (parinclude(2)) ipend=min(Zix,maxchannel)
            if (parinclude(3)) idend=min(2*Zix*Nix/Aix,maxchannel)
            if (parinclude(4)) itend=min(3*Zix*Nix/(2*Aix),maxchannel)
            if (parinclude(5)) ihend=min(3*Zix*Nix/(2*Aix),maxchannel)
            if (parinclude(6)) iaend=min(Zix*Nix/Aix,maxchannel)
          endif
          inend=min(inend,numin)
          ipend=min(ipend,numip)
          idend=min(idend,numid)
          itend=min(itend,numit)
          ihend=min(ihend,numih)
          iaend=min(iaend,numia)
c
c 3. Determine whether residual nucleus under consideration can be
c    reached by particle combination.
c    Increase running index idnum and identifier ident for particle
c    combination.
c
c chanopen : flag to open channel with first non-zero cross section
c idnumfull: flag to designate maximum number of exclusive channels
c Ztot,Ntot: number of nucleon units in exit channel
c npart    : number of particles in outgoing channel
c ident    : exclusive channel identifier
c idchannel: identifier for exclusive channel
c
          do 120 ih=0,ihend
          do 120 it=0,itend
          do 120 id=0,idend
          do 120 ia=0,iaend
          do 120 ip=0,ipend
          do 120 in=0,inend
            if (.not.chanopen(in,ip,id,it,ih,ia).and.idnumfull) goto 120
            if (idnum.eq.numchantot) goto 120
            Ztot=ip+id+it+2*ih+2*ia
            Ntot=in+id+2*it+ih+2*ia
            npart=in+ip+id+it+ih+ia
            if (npart.gt.maxchannel) goto 120
            if (Ztot.ne.Zix.or.Ntot.ne.Nix) goto 120
            ident=100000*in+10000*ip+1000*id+100*it+10*ih+ia
            idnum=idnum+1
            idchannel(idnum)=ident
c
c Initialization of arrays. Since the idnum counter may be reset at the
c end of loop 120, in the case that the exclusive cross section is
c below a threshold, the various exclusive channel arrays need to be
c initialized to zero here.
c
c xschannel     : channel cross section
c xsgamchannel  : gamma channel cross section
c xsfischannel  : fission channel cross section
c xschancheck   : integrated channel spectra
c xsfischancheck: integrated fission channel spectra
c xsratio       : ratio of exclusive cross section over residual
c                 production cross section (for exclusive gamma ray
c                 intensities)
c numlev        : maximum number of included discrete levels
c Qexcl         : Q-value for exclusive channel
c Ethrexcl      : threshold incident energy for exclusive channel
c xschaniso     : channel cross section per isomer
c exclyield     : exclusive channel yield per isomer
c xsgamdischan  : discrete gamma channel cross section
c S             : separation energy per particle
c k0            : index of incident particle
c targetE       : energy of target
c numex         : maximal number of excitation energies
c                 (set in talys.cmb)
c Eavchannel    : channel average energy
c xschannelsp   : channel cross section spectra
c xsfischannelsp: fission channel cross section spectra
c specexcl      : exclusive spectra per excitation energy
c
            xschannel(idnum)=0.
            xsgamchannel(idnum)=0.
            xsfischannel(idnum)=0.
            xschancheck(idnum)=0.
            xsfischancheck(idnum)=0.
            xsratio(idnum)=0.
            do 130 i=0,numlev
              Qexcl(idnum,i)=0.
              Ethrexcl(idnum,i)=0.
              xschaniso(idnum,i)=0.
              exclyield(idnum,i)=0.
              do 130 i2=0,numlev
                xsgamdischan(idnum,i,i2)=0.
  130       continue
            Qexcl(0,0)=S(0,0,k0)+targetE
            do 140 nexout=0,numex+1
              xsexcl(idnum,nexout)=0.
              gamexcl(idnum,nexout)=0.
  140       continue
            if (flagspec) then
              do 150 type=0,6
                Eavchannel(idnum,type)=0.
                do 160 nen=0,numen
                  xschannelsp(idnum,type,nen)=0.
                  xsfischannelsp(idnum,type,nen)=0.
                  do 170 nexout=0,numex+1
                    specexcl(idnum,type,nexout,nen)=0.
  170             continue
  160           continue
  150         continue
            endif
c
c 4. Determine source paths for exclusive channel.
c
c parskip       : logical to skip outgoing particle
c identorg,idorg: identifier for previous channel
c Zcomp         : charge number index for compound nucleus
c parZ          : charge number of particle
c Ncomp         : neutron number index for compound nucleus
c parN          : neutron number of particle
c Nlast         : last discrete level
c edis          : energy of level
c specmass      : specific mass for residual nucleus
c
c The exclusive channel under consideration may have been reached
c through different paths. Here the previous path is determined by
c looking at the type of the last emitted particle.
c
            do 210 type=0,6
              if (parskip(type)) goto 210
              identorg(type)=-1
              if (type.eq.0) idd=ident
              if (type.eq.1) idd=ident-100000
              if (type.eq.2) idd=ident-10000
              if (type.eq.3) idd=ident-1000
              if (type.eq.4) idd=ident-100
              if (type.eq.5) idd=ident-10
              if (type.eq.6) idd=ident-1
              if (idd.lt.0) goto 210
              do 220 idorg=0,idnum
                if (idchannel(idorg).eq.idd) then
                  identorg(type)=idorg
c
c The Q-value for exclusive channels is determined, both for the
c ground state and isomers.
c
                  Zcomp=Zix-parZ(type)
                  Ncomp=Nix-parN(type)
                  Qexcl(idnum,0)=Qexcl(idorg,0)-S(Zcomp,Ncomp,type)
                  do 230 nex=0,Nlast(Zix,Nix,0)
                    Qexcl(idnum,nex)=Qexcl(idnum,0)-edis(Zix,Nix,nex)
                    Ethrexcl(idnum,nex)=-(Qexcl(idnum,nex)/
     +                specmass(parZ(k0),parN(k0),k0))
                    Ethrexcl(idnum,nex)=max(Ethrexcl(idnum,nex),0.d0)
  230             continue
                  goto 210
                endif
  220         continue
  210       continue
c
c Initialize fission channel.
c
c flagfission: flag for fission
c fisfeedex  : fission contribution from excitation energy bin
c
            if (flagfission.and.Zix.eq.0.and.Nix.eq.0) then
              nex=maxex(0,0)+1
              xsfischannel(idnum)=xsfischannel(idnum)+fisfeedex(0,0,nex)
            endif
c
c 310: Loop over excitation energies of residual nucleus
c
            do 310 nexout=maxex(Zix,Nix),0,-1
              do 320 type=0,6
                if (parskip(type)) goto 320
                if (identorg(type).eq.-1) goto 320
c
c Determine the compound nucleus belonging to the emitted particle type.
c
                idorg=identorg(type)
                Zcomp=Zix-parZ(type)
                Ncomp=Nix-parN(type)
c
c 5. Exclusive cross sections per excitation energy
c
c feedexcl: feeding terms from compound excitation energy bin to
c           residual excitation energy bin
c ebegin  : first energy point of energy grid
c eend    : last energy point of energy grid
c binemis : emission spectra from initial compound nucleus
c
c The inclusive cross section per excitation energy S is
c tracked in xsexcl. The inclusive spectrum per excitation energy
c dS/dEk' is tracked in specexcl.
c
c Special case for binary channel. The binary emission spectrum has
c already been determined in subroutine binemission.
c
                if (Zcomp.eq.0.and.Ncomp.eq.0) then
                  nex=maxex(Zcomp,Ncomp)+1
                  xsexcl(idnum,nexout)=xsexcl(idnum,nexout)+
     +              feedexcl(Zcomp,Ncomp,type,nex,nexout)
                  if (type.eq.0) gamexcl(idnum,nexout)=
     +              gamexcl(idnum,nexout)+
     +              feedexcl(Zcomp,Ncomp,0,nex,nexout)
                  if (flagspec) then
                    do 330 nen=ebegin(type),eend(type)
                      specexcl(idnum,type,nexout,nen)=
     +                  specexcl(idnum,type,nexout,nen)+
     +                  binemis(type,nexout,nen)
  330               continue
                  endif
                endif
c
c 340: Loop over excitation energies of mother nucleus
c
c Calculation of inclusive cross section per excitation energy S
c
c term1,..: help variables
c popexcl : population cross section of bin just before decay
c
                do 340 nex=maxex(Zcomp,Ncomp),1,-1
                  if (feedexcl(Zcomp,Ncomp,type,nex,nexout).eq.0.)
     +              goto 340
                  term1=feedexcl(Zcomp,Ncomp,type,nex,nexout)/
     +              popexcl(Zcomp,Ncomp,nex)
                  term2=term1*xsexcl(idorg,nex)
                  xsexcl(idnum,nexout)=xsexcl(idnum,nexout)+term2
                  if (type.eq.0) then
                    gamexcl(idnum,nexout)=gamexcl(idnum,nexout)+term2
                    if (nex.le.Nlast(Zcomp,Ncomp,0))
     +                xsgamdischan(idnum,nex,nexout)=term2
                  endif
                    gamexcl(idnum,nexout)=gamexcl(idnum,nexout)+term1*
     +                gamexcl(idorg,nex)
c                 endif
c
c 6. Exclusive spectra per excitation energy
c
c Calculation of inclusive spectrum per excitation energy dS/dEk'
c
                  if (flagspec) then
                    do 350 type2=0,6
                      if (parskip(type2)) goto 350
c
c Calculation of first term of exclusive spectrum.
c
                      do 360 nen=ebegin(type2),eend(type2)
                        term3=term1*specexcl(idorg,type2,nex,nen)
                        specexcl(idnum,type2,nexout,nen)=
     +                    specexcl(idnum,type2,nexout,nen)+term3
  360                 continue
c
c Calculation of second term of exclusive spectrum. The spectrum of the
c last emitted particle is obtained by interpolating the feeding terms
c in subroutine specemission. The discrete gamma-ray production is
c excluded from the exclusive gamma spectra.
c
c NL,Nlast    : last discrete level
c specemission: subroutine for exclusive emission spectra
c specemis    : exclusive emission contribution
c
                      if (type2.eq.0.and.nex.le.Nlast(Zcomp,Ncomp,0))
     +                  goto 350
                      if (type2.eq.type) then
                        call specemission(Zcomp,Ncomp,nex,idorg,type,
     +                    nexout)
                        do 370 nen=ebegin(type),eend(type)
                          specexcl(idnum,type,nexout,nen)=
     +                      specexcl(idnum,type,nexout,nen)+
     +                      specemis(nen)
  370                   continue
                      endif
  350               continue
                  endif
  340           continue
  320         continue
c
c 7. Exclusive cross sections: total and per isomer
c
c xspopnuc: population cross section per nucleus
c
c At the end of the decay, the nucleus can end up at an isomer or the
c ground state. The exclusive cross sections for the particular channel
c can now be established. The same is true for the spectra.
c The spectra per isomer/ground state are always added.
c
              if (nexout.gt.Nlast(Zix,Nix,0)) goto 310
              if (tau(Zix,Nix,nexout).ne.0..or.nexout.eq.0) then
                xschaniso(idnum,nexout)=xsexcl(idnum,nexout)
                xschannel(idnum)=xschannel(idnum)+xsexcl(idnum,nexout)
                xsgamchannel(idnum)=xsgamchannel(idnum)+
     +            gamexcl(idnum,nexout)
                if (flagspec) then
                  do 380 type=0,6
                    if (parskip(type)) goto 380
                    do 390 nen=ebegin(type),eend(type)
                      xschannelsp(idnum,type,nen)=
     +                  xschannelsp(idnum,type,nen)+
     +                  specexcl(idnum,type,nexout,nen)
  390               continue
  380             continue
                endif
              endif
  310       continue
            if (xspopnuc(Zix,Nix).gt.0.) xsratio(idnum)=
     +        xschannel(idnum)/xspopnuc(Zix,Nix)
            channelsum=channelsum+xschannel(idnum)
c
c For non-threshold reactions (positive Q-value) we always assign a
c minimum value to the exclusive cross section. (The transmission
c coefficients for these reactions might have been zero (from ECIS),
c but non-threshold reactions theoretically have a non-zero cross
c section.)
c
c xseps: limit for cross sections
c
            if (Qexcl(idnum,0).gt.0..and.xschannel(idnum).le.xseps) then
              xschannel(idnum)=xseps
              xsgamchannel(idnum)=xseps
              xschaniso(idnum,0)=xseps
              do 400 i=1,Nlast(Zix,Nix,0)
                if (tau(Zix,Nix,i).ne.0.) then
                  xschaniso(idnum,0)=0.5*xseps
                  xschaniso(idnum,i)=0.5*xseps
                endif
  400         continue
            endif
c
c For each exclusive channel, the multiplicity or yield per
c isomer/ground state is determined.
c
            if (xschannel(idnum).ne.0.) then
              do 410 i=0,Nlast(Zix,Nix,0)
                exclyield(idnum,i)=xschaniso(idnum,i)/xschannel(idnum)
  410         continue
            endif
c
c 8. Exclusive fission cross sections
c
c All exclusive fission cross sections and spectra are determined.
c
            if (flagfission) then
              do 510 nex=maxex(Zix,Nix),1,-1
                if (popexcl(Zix,Nix,nex).ne.0.) then
                  term=fisfeedex(Zix,Nix,nex)/popexcl(Zix,Nix,nex)
                  xsfischannel(idnum)=xsfischannel(idnum)+term*
     +              xsexcl(idnum,nex)
                  if (flagspec) then
                    do 520 type=0,6
                      if (parskip(type)) goto 520
                      do 530 nen=ebegin(type),eend(type)
                        xsfischannelsp(idnum,type,nen)=
     +                    xsfischannelsp(idnum,type,nen)+
     +                    term*specexcl(idnum,type,nex,nen)
  530                 continue
  520               continue
                  endif
                endif
  510         continue
              channelsum=channelsum+xsfischannel(idnum)
            endif
c
c 9. Check.
c
c The summed exclusive cross sections should equal the particle
c production cross sections and analogously for the spectra.
c The result will appear in the output, if requested. Also the
c exclusive spectra are integrated so that they can be compared
c with the exclusive cross sections.
c
c flagcheck  : flag for output of numerical checks
c xsparcheck : total particle production cross section
c flagspec   : flag for output of spectra
c fissum     : help variable
c Especsum   : total emission energy
c emissum    : integrated emission spectrum
c Eaveragesum: help variable
c xsspeccheck: total particle production spectra
c deltaE     : energy bin around outgoing energies
c egrid      : outgoing energy grid
c eoutdis    : outgoing energy of discrete state reaction
c frac       : help variable
c Etop       : top of outgoing energy bin
c nendisc    : last discrete bin
c gmult      : continuum gamma multiplicity
c
            if (flagcheck) then
              xsparcheck(1)=xsparcheck(1)+in*xschannel(idnum)
              xsparcheck(2)=xsparcheck(2)+ip*xschannel(idnum)
              xsparcheck(3)=xsparcheck(3)+id*xschannel(idnum)
              xsparcheck(4)=xsparcheck(4)+it*xschannel(idnum)
              xsparcheck(5)=xsparcheck(5)+ih*xschannel(idnum)
              xsparcheck(6)=xsparcheck(6)+ia*xschannel(idnum)
              if (flagspec) then
                fissum=0.
                do 610 type=0,6
                  if (parskip(type)) goto 610
                  emissum=0.
                  Eaveragesum=0.
                  do 620 nen=ebegin(type),eend(type)
                    xsspeccheck(type,nen)=xsspeccheck(type,nen)+
     +                xschannelsp(idnum,type,nen)
                    emissum=emissum+xschannelsp(idnum,type,nen)*
     +                deltaE(nen)
                    Eaveragesum=Eaveragesum+egrid(nen)*
     +                xschannelsp(idnum,type,nen)*deltaE(nen)
                    if (flagfission) fissum=fissum+
     +                xsfischannelsp(idnum,type,nen)*deltaE(nen)
  620             continue
                  if (type.eq.0) then
                    do 630 i=1,numlev
                      do 630 i2=0,i
                        if (xsgamdischan(idnum,i,i2).gt.0.) then
                          emissum=emissum+xsgamdischan(idnum,i,i2)
                          Eaveragesum=Eaveragesum+
     +                      xsgamdischan(idnum,i,i2)*
     +                      (edis(Zcomp,Ncomp,i)-edis(Zcomp,Ncomp,i2))
                        endif
  630               continue
                  endif
                  if (npart.eq.1) then
                    NL=Nlast(Zix,Nix,0)
                    if (eoutdis(type,NL).gt.0.) then
                      frac=Etop(nendisc(type))-eoutdis(type,NL)
                      emissum=emissum-
     +                  xschannelsp(idnum,type,nendisc(type))*frac
                      Eaveragesum=Eaveragesum-egrid(nendisc(type))*
     +                  xschannelsp(idnum,type,nendisc(type))*frac
                      if (flagfission) fissum=fissum-
     +                  xsfischannelsp(idnum,type,nendisc(type))*frac
                    endif
                  endif
                  if (type.gt.0.or.npart.eq.0)
     +              xschancheck(idnum)=xschancheck(idnum)+emissum
                  if (emissum.gt.0.) Eavchannel(idnum,type)=
     +              Eaveragesum/emissum
  610           continue
                if (xschannel(idnum).gt.0.) then
                  gmult(idnum)=xsgamchannel(idnum)/xschannel(idnum)
                else
                  gmult(idnum)=0.
                endif
                Especsum(idnum)=gmult(idnum)*Eavchannel(idnum,0)+
     +            in*Eavchannel(idnum,1)+ip*Eavchannel(idnum,2)+
     +            id*Eavchannel(idnum,3)+it*Eavchannel(idnum,4)+
     +            ih*Eavchannel(idnum,5)+ia*Eavchannel(idnum,6)
                if (flagfission) xsfischancheck(idnum)=fissum
              endif
            endif
c
c 10. Create reaction string for output
c
c apos      : '
c reacstring: string for exclusive reaction channel
c fisstring : string for exclusive fission reaction channel
c parsym    : symbol of particle
c numchantot: maximal number of exclusive channels
c
            apos="'"
            reacstring(idnum)='( ,            '
            reacstring(idnum)(2:2)=parsym(k0)
            i=4
            if (npart.eq.0.) then
              reacstring(idnum)(4:4)='g'
              i=i+1
            endif
            if (in.eq.1) then
              write(reacstring(idnum)(i:i),'(a1)') parsym(1)
              i=i+1
            endif
            if (in.gt.1) then
              write(reacstring(idnum)(i:i+1),'(i1,a1)') in,parsym(1)
              i=i+2
            endif
            if (ip.eq.1) then
              write(reacstring(idnum)(i:i),'(a1)') parsym(2)
              i=i+1
            endif
            if (ip.gt.1) then
              write(reacstring(idnum)(i:i+1),'(i1,a1)') ip,parsym(2)
              i=i+2
            endif
            if (id.eq.1) then
              write(reacstring(idnum)(i:i),'(a1)') parsym(3)
              i=i+1
            endif
            if (id.gt.1) then
              write(reacstring(idnum)(i:i+1),'(i1,a1)') id,parsym(3)
              i=i+2
            endif
            if (it.eq.1) then
              write(reacstring(idnum)(i:i),'(a1)') parsym(4)
              i=i+1
            endif
            if (it.gt.1) then
              write(reacstring(idnum)(i:i+1),'(i1,a1)') it,parsym(4)
              i=i+2
            endif
            if (ih.eq.1) then
              write(reacstring(idnum)(i:i),'(a1)') parsym(5)
              i=i+1
            endif
            if (ih.gt.1) then
              write(reacstring(idnum)(i:i+1),'(i1,a1)') ih,parsym(5)
              i=i+2
            endif
            if (ia.eq.1) then
              write(reacstring(idnum)(i:i),'(a1)') parsym(6)
              i=i+1
            endif
            if (ia.gt.1) then
              write(reacstring(idnum)(i:i+1),'(i1,a1)') ia,parsym(6)
              i=i+2
            endif
            if (npart.eq.1) then
              if (k0.eq.1.and.in.eq.1) then
                reacstring(idnum)(i:i)=apos
                i=i+1
              endif
              if (k0.eq.2.and.ip.eq.1) then
                reacstring(idnum)(i:i)=apos
                i=i+1
              endif
              if (k0.eq.3.and.id.eq.1) then
                reacstring(idnum)(i:i)=apos
                i=i+1
              endif
              if (k0.eq.4.and.it.eq.1) then
                reacstring(idnum)(i:i)=apos
                i=i+1
              endif
              if (k0.eq.5.and.ih.eq.1) then
                reacstring(idnum)(i:i)=apos
                i=i+1
              endif
              if (k0.eq.6.and.ia.eq.1) then
                reacstring(idnum)(i:i)=apos
                i=i+1
              endif
            endif
            reacstring(idnum)(i:i)=')'
            if (flagfission)
     +        fisstring(idnum)=reacstring(idnum)(1:i-1)//'f)'
c
c Reset idnum counter in the case that the cross section is too small.
c
c opennum: total number of open channels
c
            if ((xschannel(idnum).ge.xseps.and..not.idnumfull).or.
     +        npart.eq.0) then
              if (.not.chanopen(in,ip,id,it,ih,ia)) opennum=opennum+1
              chanopen(in,ip,id,it,ih,ia)=.true.
            endif
            if (xschannel(idnum).lt.xseps.and.npart.gt.1.and.
     +        .not.chanopen(in,ip,id,it,ih,ia)) idnum=idnum-1
            if (opennum.eq.numchantot-10) idnumfull=.true.
            if (idnum.lt.0) goto 120
            if (xschannel(idnum).lt.0.) xschannel(idnum)=xseps
  120     continue
  110 continue
c
c Set threshold energy for inelastic scattering to that of first
c excited state.
c
      do 710 idc=0,idnum
        if (idchannel(idc).eq.100000) then
          Ethrexcl(idc,0)=Ethrexcl(idc,1)
          goto 720
        endif
  710 continue
  720 return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
