C     -*-fortran-*-
c     
c     
c     VERSION 2.2
c     
c     This library contains routines for B-spline interpolation in
c     one, two, and three dimensions. Part of the routines are based
c     on the book by Carl de Boor: A practical guide to Splines (Springer,
c     New-York 1978) and have the same calling sequence and names as
c     the corresponding routines from the IMSL library. For documen-
c     tation see the additional files. NOTE: The results in the demo
c     routines may vary slightly on different architectures.
c     
c     by W. Schadow 07/19/98
c     last changed by W. Schadow 07/28/2000

c     
c     
c     Wolfgang Schadow
c     TRIUMF
c     4004 Wesbrook Mall
c     Vancouver, B.C. V6T 2A3
c     Canada
c     
c     email: schadow@triumf.ca  or  schadow@physik.uni-bonn.de
c     
c     www  : http://www.triumf.ca/people/schadow
c     
c     
c     ------------------------------------------------------------------
c     
c     
c     Copyright (C) 2000 Wolfgang Schadow
c     
c     This library is free software; you can redistribute it and/or
c     modify it under the terms of the GNU Library General Public
c     License as published by the Free Software Foundation; either
c     version 2 of the License, or (at your option) any later version.
c     
c     This library is distributed in the hope that it will be useful,
c     but WITHOUT ANY WARRANTY; without even the implied warranty of
c     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
c     Library General Public License for more details.
c     
c     You should have received a copy of the GNU Library General Public
c     License along with this library; if not, write to the
c     Free Software Foundation, Inc., 59 Temple Place - Suite 330,
c     Boston, MA  02111-1307, USA.
c     
c     
c     ------------------------------------------------------------------
c     
c     
c     The following routines are included:
c     
c     dbsnak
c     
c     dbsint
c     dbsval
c     dbsder
c     dbs1gd
c     
c     dbs2in
c     dbs2dr
c     dbs2vl
c     dbs2gd
c     
c     dbs3in
c     dbs3vl
c     dbs3dr
c     dbs3gd
c     
c     ------------------------------------------------------------------
c     
c     NEW: corrected some error messages
c     some changes in the checks of dbs3dg to find a possible
c     error earlier.
c     
c     ------------------------------------------------------------------
c     
c     NEW: documentation included, changed some comments
c     
c     ------------------------------------------------------------------
c     

c     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


      subroutine dbsnak(nx,xvec,kxord,xknot)

c     
c     Compute the `not-a-knot' spline knot sequence.
c     (see de Boor p. 167)
c     
c     nx     - number of data points.  (input)
c     xvec   - array of length ndata containing the location of the
c     data points.  (input)
c     kxord  - order of the spline.  (input)
c     xknot  - array of length ndata+korder containing the knot
c     sequence.  (output)
c     

      implicit double precision (a-h,o-z), integer (i-n)

      dimension xvec(nx),xknot(nx+kxord)

      logical first

      save first,eps

      data first/.true./

      if (first) then
         first=.false.
         eps = dlamch('precision')
c     write(6,*) 'subroutine dbsnak: '
c     write(6,*) 'eps = ',eps
      endif

      if((kxord .lt. 0) .or. (kxord .gt. nx)) then
         write(6,*) 'subroutine dbsnak: error'
         write(6,*) '0 <= kxord <= nx is required.'
         write(6,*) 'kxord = ', kxord, ' and nx = ', nx,  ' is given.'
         stop
      endif

      do 30 i = 1, kxord
         xknot(i) = xvec(1)
 30   continue

      if(mod(kxord,2) .eq. 0) then
         do 40 i = kxord+1,nx
            xknot(i) = xvec(i-kxord/2)
 40      continue
      else
         do 50 i = kxord+1,nx
            xknot(i) = 0.5d0 * (xvec(i-kxord/2) + xvec(i-kxord/2-1))
 50      continue
      endif

      do 60 i = nx+1,nx+kxord
         xknot(i) = xvec(nx) * (1.0d0 + eps)
 60   continue

      return

      end


c     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


      subroutine dbsint(nx,xvec,xdata,kx,xknot,bcoef)

c     
c     Computes the spline interpolant, returning the B-spline coefficients.
c     (see de Boor p. 204)
c     
c     nx     - number of data points.  (input)
c     xvec   - array of length nx containing the data point
c     abscissas.  (input)
c     xdata  - array of length ndata containing the data point
c     ordinates.  (input)
c     kx     - order of the spline.  (input)
c     korder must be less than or equal to ndata.
c     xknot  - array of length nx+kx containing the knot
c     sequence.  (input)
c     xknot must be nondecreasing.
c     bscoef - array of length ndata containing the B-spline
c     coefficients.  (output)
c     

      implicit double precision (a-h,o-z), integer (i-n)

      parameter(kxmax=8,nxmax=5000)

      dimension xdata(nx),xvec(nx),xknot(nx+kx),bcoef(nx)
      dimension work((2*kxmax-1)*nxmax)

      if ((kx .gt. kxmax) .or. (nx .gt. nxmax)) then
         write(6,*) 'subroutine dbsint: error'
         write(6,*) 'kx > kxmax or nx > nxmax'
         write(6,*) 'kx = ',kx,'  kxmax = ',kxmax
         write(6,*) 'nx = ',nx,'  nxmax = ',nxmax
         stop
      endif

      nxp1  = nx + 1
      kxm1  = kx - 1
      kpkm2 = 2 * kxm1
      leftx = kx
      lenq  = nx * (kx + kxm1)

      do 10 i = 1, lenq
         work(i) = 0.d0
 10   continue

      do 20 ix = 1,nx
         xveci  = xvec(ix)
         ilp1mx = min0(ix+kx,nxp1)
         leftx   = max0(leftx,ix)
         if (xveci .lt. xknot(leftx)) goto 998
 30      if (xveci .lt. xknot(leftx+1)) go to 40
         leftx = leftx + 1
         if (leftx .lt. ilp1mx) go to 30
         leftx = leftx - 1
         if (xveci .gt. xknot(leftx+1)) goto 998
 40      call bsplvb (xknot,nx+kx,kx,1,xveci,leftx,bcoef)
         jj = ix - leftx + 1 + (leftx - kx) * (kx + kxm1)
         do 50 ik = 1,kx
            jj       = jj + kpkm2
            work(jj) = bcoef(ik)
 50      continue
 20   continue

      call banfac(work,kx+kxm1,nx,kxm1,kxm1,iflag)
      go to (60,999), iflag

 60   do 70 ix = 1,nx
         bcoef(ix) = xdata(ix)
 70   continue

      call banslv(work,kx+kxm1,nx,kxm1,kxm1,bcoef)

      return

 998  write(6,*) 'subroutine dbsint:'
      write(6,*) 'xknot(ix) <= xknot(ix+1) required.'
      write(6,*) ix,xknot(ix),xknot(ix+1)

      stop

 999  write(6,*) 'subroutine dbsint: error'
      write(6,*) 'no solution of linear equation system !!!'

      stop

      end


c     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


      function dbsval(x,kx,xknot,nx,bcoef)

c     
c     Evaluates a spline, given its B-spline representation.
c     
c     x      - point at which the spline is to be evaluated.  (input)
c     kx     - order of the spline.  (input)
c     xknot  - array of length nx+kx containing the knot
c     sequence.  (input)
c     xknot must be nondecreasing.
c     nx     - number of B-spline coefficients.  (input)
c     bcoef  - array of length nx containing the B-spline
c     coefficients.  (input)
c     dbsval - value of the spline at x.  (output)
c     

      implicit double precision (a-h,o-z), integer (i-n)

      parameter(kxmax=8)

      dimension xknot(nx+kx),bcoef(nx)
      dimension work(kxmax),dl(kxmax),dr(kxmax)


c     
c     check if kx <= kxmax
c     

      if(kx .gt. kxmax) then
         write(6,*) 'subroutine dbsval:'
         write(6,*) 'kx <= kxmax required.'
         stop
      endif

c     
c     check if xknot(i) <= xknot(i+1) and calculation of i so that
c     xknot(i) <= x < xknot(i+1)
c     

      leftx = 0

      do 10 ix = 1,nx+kx-1
         if (xknot(ix) .gt. xknot(ix+1)) then
            write(6,*) 'subroutine dbsval:'
            write(6,*) 'xknot(ix) <= xknot(ix+1) required.'
            write(6,*) ix,xknot(ix),xknot(ix+1)
            stop
         endif
         if((xknot(ix) .le. x) .and. (x .lt. xknot(ix+1))) leftx = ix
 10   continue

      if(leftx .eq. 0) then
         write(6,*) 'subroutine dbsval:'
         write(6,*) 'ix with xknot(ix) <= x < xknot(ix+1) required.'
         write(6,*) 'x = ', x
         stop
      endif

      do 20 ik = 1,kx-1
         work(ik) = bcoef(leftx+ik-kx)
         dl(ik)   = x - xknot(leftx+ik-kx)
         dr(ik)   = xknot(leftx+ik) - x
 20   continue

      work(kx)  = bcoef(leftx)
      dl(kx)    = x - xknot(leftx)

      do 30 ik = 1,kx-1
         save2 = work(ik)
         do 40  il = ik+1,kx
            save1 = work(il)
            work(il) = (dl(il) * work(il) + dr(il-ik) * save2)
     .           / (dl(il) + dr(il - ik))
            save2 = save1
 40      continue
 30   continue

      dbsval = work(kx)

      return

      end


c     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


      function dbsder(iderx,x,kx,xknot,nx,bcoef)

c     
c     Evaluates the derivative of a spline, given its B-spline representation.
c     
c     
c     iderx  - order of the derivative to be evaluated.  (input)
c     in particular, iderx = 0 returns the value of the
c     spline.
c     x      - point at which the spline is to be evaluated.  (input)
c     kx     - order of the spline.  (input)
c     xknot  - array of length nx+kx containing the knot
c     sequence.  (input)
c     xknot must be nondecreasing.
c     nx     - number of B-spline coefficients.  (input)
c     bcoef  - array of length nx containing the B-spline
c     coefficients.  (input)
c     dbsder - value of the iderx-th derivative of the spline at x.
c     (output)
c     

      implicit double precision (a-h,o-z), integer (i-n)

      parameter(kxmax=8)

      dimension xknot(nx+kx),bcoef(nx)
      dimension work(kxmax),dl(kxmax),dr(kxmax),bsp(kxmax)

c     
c     check if <= kxmax
c     

      if(kx .gt. kxmax) then
         write(6,*) 'subroutine dbsder:'
         write(6,*) 'kx <= kxmax required.'
         stop
      endif

c     
c     check if xknot(i) <= xknot(i+1) and calculation of i so that
c     xknot(i) <= x < xknot(i+1)
c     

      leftx = 0
      do 10 ix = 1,nx+kx-1
         if (xknot(ix) .gt. xknot(ix+1)) then
            write(6,*) 'subroutine dbsder:'
            write(6,*) 'xknot(ix) <= xknot(ix+1) required.'
            stop
         endif
         if ((xknot(ix) .le. x) .and. (x .lt. xknot(ix+1))) leftx = ix
 10   continue

      if (leftx .eq. 0) then
         write(6,*) 'subroutine dbsder:'
         write(6,*) 'ix with xknot(ix) <= x < xknot(ix+1) required.'
         write(6,*) 'xknot(1)     = ', xknot(1)
         write(6,*) 'xknot(nx+kx) = ', xknot(nx+kx)
         write(6,*) '         x   = ', x
         stop
      endif

      if (iderx .eq. 0) then

         do 20 ik = 1,kx-1
            work(ik) = bcoef(leftx+ik-kx)
            dl(ik)   = x - xknot(leftx+ik-kx)
            dr(ik)   = xknot(leftx+ik) - x
 20      continue

         work(kx)  = bcoef(leftx)
         dl(kx)    = x - xknot(leftx)

         do 30 ik = 1,kx-1
            save2 = work(ik)
            do 40  il = ik+1,kx
               save1 = work(il)
               work(il) = (dl(il) * work(il) + dr(il-ik) * save2)
     .              / (dl(il) + dr(il - ik))
               save2 = save1
 40         continue
 30      continue

         dbsder = work(kx)

      elseif ((iderx .ge. 1) .and. (iderx .lt. kx)) then
         bsp(1) = 1.0d0
         do 50 ik = 1,kx-iderx-1
            dr(ik) = xknot(leftx+ik) - x
            dl(ik) = x - xknot(leftx+1-ik)
            save   = bsp(1)
            bsp(1) = 0.0d0
            do 60 il = 1,ik
               y         = save / (dr(il) + dl(ik+1-il))
               bsp(il)   = bsp(il) + dr(il) * y
               save      = bsp(il+1)
               bsp(il+1) = dl(ik+1-il) * y
 60         continue
 50      continue

         do 70 ik = 1,kx
            work(ik) = bcoef(leftx+ik-kx)
            dr(ik)   = xknot(leftx+ik) - x
            dl(ik)   = x - xknot(leftx+ik-kx)
 70      continue

         do 80 ik = 1,iderx
            dik   = dble(kx - ik)
            save2 = work(ik)
            do 90  il = ik+1,kx
               save1    = work(il)
               work(il) = dik * (work(il) - save2) /(dl(il) + dr(il-ik))
               save2    = save1
 90         continue
 80      continue

         sum = 0.0d0

         do 100 i = 1,kx-iderx
            sum = sum + bsp(i) * work(iderx+i)
 100     continue

         dbsder = sum

      else
         dbsder = 0.0d0
      endif

      return

      end


c     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


      subroutine dbs1gd(iderx,nxvec,xvec,kx,xknot,nx,bcoef,val)

c     
c     Evaluates the derivative of a spline on a grid, given its B-spline
c     representation.
c     
c     iderx  - order of the derivative to be evaluated.  (input)
c     in particular, iderx = 0 returns the value of the
c     spline.
c     nxvec  - length of vector xvec.  (input)
c     xvec   - array of length nxvec containing the points at which the
c     spline is to be evaluated.  (input)
c     xvec should be strictly increasing.
c     kx     - order of the spline.  (input)
c     xknot  - array of length nx+kx containing the knot
c     sequence.  (input)
c     xknot must be nondecreasing.
c     nx     - number of B-spline coefficients.  (input)
c     bcoef  - array of length nx containing the B-spline
c     coefficients.  (input)
c     val    - array of length nxvec containing the values of the
c     iderx-th derivative of the spline at the points in
c     xvec.  (output)
c     

      implicit double precision (a-h,o-z), integer (i-n)

      parameter(kxmax=8,nxmax=4000)

      dimension xvec(nxvec),xknot(nx+kx),bcoef(nx),val(nxvec)
      dimension dl(nxmax,kxmax),dr(nxmax,kxmax),save1(nxmax)
      dimension biatx(nxmax,kxmax),term(nxmax),save2(nxmax)
      dimension leftx(nxmax),work(nxmax,kxmax)

      logical same,next

c     
c     check if kx <= kxmax
c     

      if(kx .gt. kxmax) then
         write(6,*) 'subroutine dbs1gd:'
         write(6,*) 'kx <= kxmax required.'
         stop
      endif

      leftx(1) = 0

      call huntn(xknot,nx+kx,kx,xvec(1),leftx(1))

      do 10 ix = 2,nxvec
         leftx(ix) = leftx(ix-1)
         same = (xknot(leftx(ix)) .le. xvec(ix))
     .        .and. (xvec(ix) .le. xknot(leftx(ix)+1))
         if(.not. same ) then
            leftx(ix) = leftx(ix) + 1
            next      = (xknot(leftx(ix)) .le. xvec(ix))
     .           .and. (xvec(ix) .le. xknot(leftx(ix)+1))
            if (.not. next)
     .           call huntn(xknot,nx+kx,kx,xvec(ix),leftx(ix))
         endif
 10   continue

      do 20 i = 1,nx+kx-1
         if (xknot(i) .gt. xknot(i+1)) then
            write(6,*) 'subroutine dbs1gd:'
            write(6,*) 'xknot(i) <= xknot(i+1) required.'
            write(6,*) i, xknot(i), xknot(i+1)
            write(6,*)
            write(6,*) xknot
            stop
         endif
 20   continue

      do 30 i = 1,nxvec
         if ((xvec(i).lt.xknot(1)).or.(xvec(i).gt.xknot(nx+kx))) then
            write(6,*) 'subroutine dbs1gd:'
            write(6,*) 'ix with xknot(ix) <= x < xknot(ix+1) required.'
            write(6,*) 'x = ', xvec(i)
            stop
         endif
 30   continue

      if (iderx .eq. 0) then

         do 40 ix = 1,nxvec
            biatx(ix,1) = 1.d0
            val(ix)     = 0.d0
 40      continue

         do 50 ik = 1,kx-1
            do 60 ix = 1,nxvec
               dr(ix,ik) = xknot(leftx(ix)+ik) - xvec(ix)
               dl(ix,ik) = xvec(ix) - xknot(leftx(ix)+1-ik)
               save1(ix) = 0.d0
 60         continue

            do 70 il = 1,ik
               do 80 ix = 1,nxvec
                  term(ix)     = biatx(ix,il)
     .                 / (dr(ix,il) + dl(ix,ik+1-il))
                  biatx(ix,il) = save1(ix) + dr(ix,il) * term(ix)
                  save1(ix)    = dl(ix,ik+1-il) * term(ix)
 80            continue
 70         continue

            do 90 ix = 1,nxvec
               biatx(ix,ik+1) = save1(ix)
 90         continue
 50      continue

         do 100 ik = 1,kx
            do 110 ix = 1,nxvec
               val(ix) = val(ix) + biatx(ix,ik) * bcoef(leftx(ix)-kx+ik)
 110        continue
 100     continue

      elseif ((iderx .ge. 1) .and. (iderx .lt. kx)) then

         do 120 ix = 1,nxvec
            biatx(ix,1) = 1.d0
            val(ix)     = 0.d0
 120     continue

         do 130 ik = 1,kx-iderx-1
            do 140 ix = 1,nxvec
               dr(ix,ik)   = xknot(leftx(ix)+ik) - xvec(ix)
               dl(ix,ik)   = xvec(ix) - xknot(leftx(ix)+1-ik)
               save1(ix)    = biatx(ix,1)
               biatx(ix,1) = 0.0d0
               do 150 il = 1,ik
                  term(ix)       = save1(ix)
     .                 / (dr(ix,il) + dl(ix,ik+1-il))
                  biatx(ix,il)   = biatx(ix,il) + dr(ix,il) * term(ix)
                  save1(ix)      = biatx(ix,il+1)
                  biatx(ix,il+1) = dl(ix,ik+1-il) * term(ix)
 150           continue
 140        continue
 130     continue

         do 160 ik = 1,kx
            do 170 ix = 1,nxvec
               work(ix,ik) = bcoef(leftx(ix)+ik-kx)
               dr(ix,ik)   = xknot(leftx(ix)+ik) - xvec(ix)
               dl(ix,ik)   = xvec(ix) - xknot(leftx(ix)+ik-kx)
 170        continue
 160     continue

         do 180 ik = 1,iderx
            dik   = dble(kx - ik)
            do 190 ix = 1,nxvec
               save2(ix) = work(ix,ik)
               do 200  il = ik+1,kx
                  save1(ix)   = work(ix,il)
                  work(ix,il) = dik * (work(ix,il) - save2(ix))
     .                 /(dl(ix,il) + dr(ix,il-ik))
                  save2(ix)   = save1(ix)
 200           continue
 190        continue
 180     continue

         do 210 i = 1,kx-iderx
            do 220 ix = 1,nxvec
               val(ix) = val(ix) + biatx(ix,i) * work(ix,iderx+i)
 220        continue
 210     continue

      else

         do 230 ix = 1,nxvec
            val(ix) = 0.0d0
 230     continue

      endif

      return

      end


c     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


      function dbsdca(iderx,x,kx,xknot,nx,bcoef,leftx)

c     
c     This routine is equivalent to the routine dbsder, but it does not
c     check the parameters!!c
c     
c     Evaluates the derivative of a spline, given its B-spline representation.
c     
c     
c     iderx  - order of the derivative to be evaluated.  (input)
c     in particular, iderx = 0 returns the value of the
c     spline.
c     x      - point at which the spline is to be evaluated.  (input)
c     kx     - order of the spline.  (input)
c     xknot  - array of length nx+kx containing the knot
c     sequence.  (input)
c     xknot must be nondecreasing.
c     nx     - number of B-spline coefficients.  (input)
c     bcoef  - array of length nx containing the B-spline
c     coefficients.  (input)
c     leftx  - number of the intervall of xknot that includes x
c     dbsdca - value of the ideriv-th derivative of the spline at x.
c     (output)
c     

      implicit double precision (a-h,o-z), integer (i-n)

      parameter(kxmax=8)

      dimension xknot(nx+kx),bcoef(nx)
      dimension work(kxmax),dl(kxmax),dr(kxmax),bsp(kxmax)


      if (iderx .eq. 0) then

         do 20 ik = 1,kx-1
            work(ik) = bcoef(leftx+ik-kx)
            dl(ik)   = x - xknot(leftx+ik-kx)
            dr(ik)   = xknot(leftx+ik) - x
 20      continue

         work(kx)  = bcoef(leftx)
         dl(kx)    = x - xknot(leftx)

         do 30 ik = 1,kx-1
            save2 = work(ik)
            do 40  il = ik+1,kx
               save1 = work(il)
               work(il) = (dl(il) * work(il) + dr(il-ik) * save2)
     .              / (dl(il) + dr(il - ik))
               save2 = save1
 40         continue
 30      continue

         dbsdca = work(kx)

      elseif ((iderx .ge. 1) .and. (iderx .lt. kx)) then
         bsp(1) = 1.0d0
         do 50 ik = 1,kx-iderx-1
            dr(ik) = xknot(leftx+ik) - x
            dl(ik) = x - xknot(leftx+1-ik)
            save   = bsp(1)
            bsp(1) = 0.0d0
            do 60 il = 1,ik
               y         = save / (dr(il) + dl(ik+1-il))
               bsp(il)   = bsp(il) + dr(il) * y
               save      = bsp(il+1)
               bsp(il+1) = dl(ik+1-il) * y
 60         continue
 50      continue

         do 70 ik = 1,kx
            work(ik) = bcoef(leftx+ik-kx)
            dr(ik)   = xknot(leftx+ik) - x
            dl(ik)   = x - xknot(leftx+ik-kx)
 70      continue

         do 80 ik = 1,iderx
            dik   = dble(kx - ik)
            save2 = work(ik)
            do 90  il = ik+1,kx
               save1    = work(il)
               work(il) = dik * (work(il) - save2) /(dl(il) + dr(il-ik))
               save2    = save1
 90         continue
 80      continue

         sum = 0.0d0

         do 100 i = 1,kx-iderx
            sum = sum + bsp(i) * work(iderx+i)
 100     continue

         dbsdca = sum

      else
         dbsdca = 0.0d0
      endif

      return

      end


c     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


      subroutine bsplvb(t,n,jhigh,index,x,left,biatx)

      implicit double precision (a-h,o-z), integer (i-n)

      parameter(kmax=8)

      dimension biatx(jhigh),t(n),dl(kmax),dr(kmax)

      data j/1/

      go to (10,20), index

 10   j = 1

      biatx(1) = 1.d0

      if (j .ge. jhigh) go to 99

 20   jp1 = j + 1

      dr(j) = t(left+j) - x
      dl(j) = x - t(left+1-j)
      saved = 0.d0

      do 30 i = 1, j
         term     = biatx(i) / (dr(i) + dl(jp1-i))
         biatx(i) = saved + dr(i) * term
         saved    = dl(jp1-i) * term
 30   continue

      biatx(jp1) = saved
      j          = jp1

      if (j .lt. jhigh) go to 20

 99   return

      end


c     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


      subroutine banfac(w,nroww,nrow,nbandl,nbandu,iflag)

      implicit double precision (a-h,o-z), integer (i-n)

      dimension w(nroww,nrow)

      iflag  = 1
      middle = nbandu + 1
      nrowm1 = nrow - 1

      if (nrowm1) 999,900,10

 10   if (nbandl .gt. 0) go to 30

      do 20 i = 1,nrowm1
         if (w(middle,i) .eq. 0.d0) go to 999
 20   continue

      go to 900

 30   if (nbandu .gt. 0) go to 60
      do 40 i = 1,nrowm1
         pivot = w(middle,i)
         if(pivot .eq. 0.d0) go to 999
         jmax = min0(nbandl, nrow - i)
         do 50 j = 1,jmax
            w(middle+j,i) = w(middle+j,i) / pivot
 50      continue
 40   continue

      return

 60   do 70 i = 1,nrowm1
         pivot = w(middle,i)
         if (pivot .eq. 0.d0) go to 999
         jmax = min0(nbandl,nrow - i)
         do 80 j = 1,jmax
            w(middle+j,i) = w(middle+j,i) / pivot
 80      continue

         kmax = min0(nbandu,nrow - i)

         do 90 k = 1,kmax
            ipk    = i + k
            midmk  = middle - k
            factor = w(midmk,ipk)
            do 100 j = 1,jmax
               w(midmk+j,ipk) = w(midmk+j,ipk) - w(middle+j,i)
     .              * factor
 100        continue
 90      continue
 70   continue

 900  if (w(middle,nrow) .ne. 0.d0) return
 999  iflag = 2

      return

      end


c     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


      subroutine banslv(w,nroww,nrow,nbandl,nbandu,b)

      implicit double precision (a-h,o-z), integer (i-n)

      dimension w(nroww,nrow),b(nrow)

      middle = nbandu + 1
      if (nrow .eq. 1) goto 99
      nrowm1 = nrow - 1
      if (nbandl .eq. 0) goto 30

      do 10 i = 1, nrowm1
         jmax = min0(nbandl, nrow - i)
         do 20 j = 1, jmax
            b(i+j) = b(i+j) - b(i) * w(middle+j,i)
 20      continue
 10   continue

 30   if (nbandu .gt. 0)  goto 50
      do 40 i = 1, nrow
         b(i) = b(i) / w(1,i)
 40   continue

      return

 50   do 60 i = nrow, 2, -1
         b(i) = b(i)/w(middle,i)
         jmax = min0(nbandu,i-1)
         do 70 j = 1,jmax
            b(i-j) = b(i-j) - b(i) * w(middle-j,i)
 70      continue
 60   continue

 99   b(1) = b(1) / w(middle,1)

      return

      end


c     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


      subroutine huntn(xx,n,kord,x,jlo)

      implicit double precision (a-h,o-z), integer (i-n)

      dimension xx(n)

c     
c     works only for B-Splines (order n)
c     

      max  = n - kord
      null = kord

      if (jlo.le.null.or.jlo.gt.max) then
         jlo = null
         jhi = max+1
         goto 30
      endif

      inc = 1

      if (x .ge. xx(jlo)) then
 10      jhi = jlo + inc
         if (jhi .gt. max) then
            jhi = max + 1
         else if (x .ge. xx(jhi)) then
            jlo = jhi
            inc = inc + inc
            goto 10
         endif
      else
         jhi = jlo
 20      jlo = jhi - inc
         if (jlo .le. null) then
            jlo = null
         else if (x .lt. xx(jlo)) then
            jhi = jlo
            inc = inc + inc
            goto 20
         endif
      endif

 30   if (jhi-jlo.eq.1) return

      jm = (jhi + jlo) / 2
      if (x .gt. xx(jm)) then
         jlo = jm
      else
         jhi = jm
      endif

      goto 30

      end
