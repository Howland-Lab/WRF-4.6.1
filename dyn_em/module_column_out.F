!WRF
!
MODULE module_column_out

CONTAINS

  SUBROUTINE compute_column_output( ucol, vcol, wcol, tcol, pcol,       &
                                    tkecol, qkecol, tkepblcol,          &
                                    exchmcol, qcol,                     &
                                    t11col, t22col, t33col,             &
				    t12col, t13col, t23col,             &
                                    h1col, h2col, h3col,                &
				    q1col, q2col, q3col,                &
                                    zcol, uscol, z0col, col_ij,         &
                                    col_i, col_j,                       &
                                    u, v, w, t, p, pb, ph, phb,         &
                                    tke, qke, tke_pbl, exch_m, q,       &
                                    t11, t22, t33, t12, t13, t23,       &
				    h1, h2, h3, q1, q2, q3,             &
                                    ust, z0, ht,                        &
                                    ids, ide, jds, jde, kds, kde,       &
                                    ims, ime, jms, jme, kms, kme,       &
                                    its, ite, jts, jte, kts, kte        )

!=======================================================================
!
!  This subroutine prepares a 1D column of parameters for an auxiliary 
!  output stream. All variables are at cell center horizontal locations.
!
!  The approach is to find the patch upon which the desired column
!  of data reside, and compute the variables (cell center, z, etc.).
!  Next, transfer the column information from the patch upon which it
!  resides onto the master processor where it will be written to the
!  output file. This is accomplished by setting each local patch
!  column value to a large negative number. Then, for the patch
!  upon which the column exists, we compute the desired quantities.
!  These values will be greater than the large negative values
!  originally prescribed. Then we use the built-in function
!  wrf_dm_max_real to select the maximum value from all the patches,
!  which recovers the desired column.
!    
!  Perhaps not the most elegant approach, but it seems to work.
!    
!  Wrtitten by Jeff Mirocha (LLNL) May, 2021
!
!======================================================================= 

#ifdef DM_PARALLEL
    USE module_dm
#endif
    
   IMPLICIT NONE
   
   INTEGER, INTENT(IN) :: col_ij

   REAL, DIMENSION(col_ij, kms:kme), INTENT( OUT ) ::    &
                     ucol, vcol, wcol, tcol, pcol,       &
                     tkecol, qkecol, tkepblcol,          &
                     exchmcol, qcol,                     &
                     t11col, t22col, t33col,             &
                     t12col, t13col, t23col,             &
                     h1col, h2col, h3col,                &
                     q1col, q2col, q3col, zcol 

   REAL, DIMENSION(col_ij), INTENT( OUT ) :: uscol, z0col

   REAL, DIMENSION(ims:ime,kms:kme,jms:jme), INTENT( IN ) :: u, v, w, t, p, pb, ph, phb,   &
                                                             tke, qke, tke_pbl, exch_m, q, &
                                                             t11, t22, t33, t12, t13, t23, & 
                                                             h1, h2, h3, q1, q2, q3 
                                                         
   REAL, DIMENSION(ims:ime,jms:jme), INTENT( IN ) :: ust, z0, ht

   INTEGER, DIMENSION(col_ij), INTENT( IN ) :: col_i, col_j

   INTEGER, INTENT( IN  ) :: ids, ide, jds, jde, kds, kde, &
                             ims, ime, jms, jme, kms, kme, &
                             its, ite, jts, jte, kts, kte

! local variable
   INTEGER :: n, k, nij

   nij = col_ij
!-----------------------------------------------------------------------
! Executable starts here
!
   DO n=1, nij
     DO k=kts,kte - 1
        ucol(n,k) = -999.0
        vcol(n,k) = -999.0
        wcol(n,k) = -999.0
        tcol(n,k) = -999.0
	exchmcol(n,k) = -999.0
	qcol(n,k) = -999.0
        pcol(n,k) = -999.0
        zcol(n,k) = -999.0
        tkecol(n,k) = -999.0
	qkecol(n,k) = -999.0
	tkepblcol(n,k) = -999.0
        t11col(n,k) = -999.0
        t22col(n,k) = -999.0
        t33col(n,k) = -999.0
        t12col(n,k) = -999.0
        t13col(n,k) = -999.0
        t23col(n,k) = -999.0
	h1col(n,k)  = -999.0
        h2col(n,k)  = -999.0 	
        h3col(n,k)  = -999.0
	q1col(n,k)  = -999.0
	q2col(n,k)  = -999.0
	q3col(n,k)  = -999.0
     ENDDO
!     wcol(n,kte) = -999.0
!     zcol(n,kte) = -999.0
!     t13col(n,kte) = -999.0
!     t23col(n,kte) = -999.0
!     h3col(n,kte)  = -999.0
!     q3col(n,kte)  = -999.0
     uscol(n) = -999.0
     z0col(n) = -999.0
   ENDDO
! If this is the processor containing the desired column, compute/assign the values.
   DO n = 1, nij
     IF ( ( col_i(n) .GE. its ) .AND. ( col_i(n) .LE. ite ) .AND. (col_j(n) .GE. jts ) .AND. ( col_j(n) .LE. jte ) ) THEN
        DO k = kts,kte - 1
          ucol(n,k) = 0.5*( u(col_i(n),k,col_j(n)) + u(col_i(n)+1,k,col_j(n)) )
          vcol(n,k) = 0.5*( v(col_i(n),k,col_j(n)) + v(col_i(n),k,col_j(n)+1) )
          wcol(n,k) = 0.5*( w(col_i(n),k,col_j(n)) + w(col_i(n),k+1,col_j(n)) )
	  tkepblcol(n,k) = 0.5*( tke_pbl(col_i(n),k,col_j(n)) + tke_pbl(col_i(n),k+1,col_j(n)) )
	  exchmcol(n,k) = 0.5*( exch_m(col_i(n),k,col_j(n)) + exch_m(col_i(n),k+1,col_j(n)) )
          tcol(n,k) = t(col_i(n),k,col_j(n))
	  qcol(n,k) = q(col_i(n),k,col_j(n))
          pcol(n,k) = pb(col_i(n),k,col_j(n)) + p(col_i(n),k,col_j(n))
          zcol(n,k) = ( 0.5*( ( phb(col_i(n),k,col_j(n)) + ph(col_i(n),k,col_j(n)) ) + &
                              ( phb(col_i(n),k+1,col_j(n)) + ph(col_i(n),k+1,col_j(n)) ) ) )/9.81 - ht(col_i(n),col_j(n))
          tkecol(n,k) = tke(col_i(n),k,col_j(n))
	  qkecol(n,k) = qke(col_i(n),k,col_j(n))
          t11col(n,k) = t11(col_i(n),k,col_j(n))
          t22col(n,k) = t22(col_i(n),k,col_j(n))      
          t33col(n,k) = t33(col_i(n),k,col_j(n))   
          t12col(n,k) = 0.25*( t12(col_i(n),k,col_j(n)) + t12(col_i(n)+1,k,col_j(n)) + t12(col_i(n)+1,k,col_j(n)+1) + t12(col_i(n),k,col_j(n)+1) )
          t13col(n,k) = 0.25*( t13(col_i(n),k,col_j(n)) + t13(col_i(n)+1,k,col_j(n)) + t13(col_i(n)+1,k+1,col_j(n)) + t13(col_i(n),k+1,col_j(n)) )
          t23col(n,k) = 0.25*( t23(col_i(n),k,col_j(n)) + t23(col_i(n),k,col_j(n)+1) + t23(col_i(n),k+1,col_j(n)+1) + t23(col_i(n),k+1,col_j(n)) )
	  h1col(n,k)  = 0.5*( h1(col_i(n),k,col_j(n)) + h1(col_i(n)+1,k,col_j(n)) )
	  h2col(n,k)  = 0.5*( h2(col_i(n),k,col_j(n)) + h2(col_i(n),k,col_j(n)+1) )
          h3col(n,k)  = 0.5*( h3(col_i(n),k,col_j(n)) + h3(col_i(n),k+1,col_j(n)) )
          q1col(n,k)  = 0.5*( q1(col_i(n),k,col_j(n)) + q1(col_i(n)+1,k,col_j(n)) )
          q2col(n,k)  = 0.5*( q2(col_i(n),k,col_j(n)) + q2(col_i(n),k,col_j(n)+1) )
	  q3col(n,k)  = 0.5*( q3(col_i(n),k,col_j(n)) + q3(col_i(n),k+1,col_j(n)) )
         ENDDO
!         wcol(n,kte) = w(col_i(n),kte,col_j(n))
!         zcol(n,kte) = ( phb(col_i(n),kte,col_j(n)) + ph(col_i(n),kte,col_j(n)) )/9.81 - ht(col_i(n),col_j(n))
!         t13col(n,kte) = 0.5*( t13(col_i(n),kte,col_j(n)) + t13(col_i(n)+1,kte,col_j(n)) )
!         t23col(n,kte) = 0.5*( t23(col_i(n),kte,col_j(n)) + t23(col_i(n),kte,col_j(n)+1) )
!         h3col(n,kte)  = h3(col_i(n),kte,col_j(n))
!         q3col(n,kte)  = q3(col_i(n),kte,col_j(n))
         uscol(n) = ust(col_i(n),col_j(n))  
         z0col(n) = z0(col_i(n),col_j(n))  
     ENDIF
   ENDDO
! Select maximum value from all processors
   DO n=1,nij
     DO k=kts,kte - 1
       ucol(n,k) = wrf_dm_max_real(ucol(n,k)) 
       vcol(n,k) = wrf_dm_max_real(vcol(n,k))
       wcol(n,k) = wrf_dm_max_real(wcol(n,k))
       tcol(n,k) = wrf_dm_max_real(tcol(n,k))
       qcol(n,k) = wrf_dm_max_real(qcol(n,k))
       pcol(n,k) = wrf_dm_max_real(pcol(n,k))
       zcol(n,k) = wrf_dm_max_real(zcol(n,k))
       tkecol(n,k) = wrf_dm_max_real(tkecol(n,k))
       qkecol(n,k) = wrf_dm_max_real(qkecol(n,k))
       tkepblcol(n,k) = wrf_dm_max_real(tkepblcol(n,k))
       exchmcol(n,k) = wrf_dm_max_real(exchmcol(n,k))
       t11col(n,k) = wrf_dm_max_real(t11col(n,k))
       t22col(n,k) = wrf_dm_max_real(t22col(n,k))
       t33col(n,k) = wrf_dm_max_real(t33col(n,k))
       t12col(n,k) = wrf_dm_max_real(t12col(n,k))
       t13col(n,k) = wrf_dm_max_real(t13col(n,k))
       t23col(n,k) = wrf_dm_max_real(t23col(n,k))
       h1col(n,k)  = wrf_dm_max_real(h1col(n,k))
       h2col(n,k)  = wrf_dm_max_real(h2col(n,k))
       h3col(n,k)  = wrf_dm_max_real(h3col(n,k))
       q1col(n,k)  = wrf_dm_max_real(q1col(n,k))
       q2col(n,k)  = wrf_dm_max_real(q2col(n,k))
       q3col(n,k)  = wrf_dm_max_real(q3col(n,k))
     ENDDO
!     wcol(n,kte) =  wrf_dm_max_real(wcol(n,kte))
!     zcol(n,kte) =  wrf_dm_max_real(zcol(n,kte))
!     t13col(n,kte) =  wrf_dm_max_real(t13col(n,kte))
!     t23col(n,kte) =  wrf_dm_max_real(t23col(n,kte))
!     h3col(n,kte)  =  wrf_dm_max_real(h3col(n,kte))
!     q3col(n,kte)  =  wrf_dm_max_real(q3col(n,kte))
     uscol(n) =  wrf_dm_max_real(uscol(n)) 
     z0col(n) =  wrf_dm_max_real(z0col(n))
   ENDDO

! Set out-of-range values to zero using WHERE
!   WHERE (ucol .LE. -999.0 .OR. ucol .GE. 9999999.0) ucol = 0.0
!   WHERE (vcol .LE. -999.0 .OR. vcol .GE. 9999999.0) vcol = 0.0
!   WHERE (wcol .LE. -999.0 .OR. wcol .GE. 9999999.0) wcol = 0.0
!   WHERE (tcol .LE. -999.0 .OR. tcol .GE. 9999999.0) tcol = 0.0
!   WHERE (qcol .LE. -999.0 .OR. qcol .GE. 9999999.0) qcol = 0.0
!   WHERE (pcol .LE. -999.0 .OR. pcol .GE. 9999999.0) pcol = 0.0
!   WHERE (zcol .LE. -999.0 .OR. zcol .GE. 9999999.0) zcol = 0.0
!   WHERE (tkecol .LE. -999.0 .OR. tkecol .GE. 9999999.0) tkecol = 0.0
!   WHERE (qkecol .LE. -999.0 .OR. qkecol .GE. 9999999.0) qkecol = 0.0
!   WHERE (tkepblcol .LE. -999.0 .OR. tkepblcol .GE. 9999999.0) tkepblcol = 0.0
!   WHERE (exchmcol .LE. -999.0 .OR. exchmcol .GE. 9999999.0) exchmcol = 0.0
!   WHERE (t11col .LE. -999.0 .OR. t11col .GE. 9999999.0) t11col = 0.0
!   WHERE (t22col .LE. -999.0 .OR. t22col .GE. 9999999.0) t22col = 0.0
!   WHERE (t33col .LE. -999.0 .OR. t33col .GE. 9999999.0) t33col = 0.0
!   WHERE (t12col .LE. -999.0 .OR. t12col .GE. 9999999.0) t12col = 0.0
!   WHERE (t13col .LE. -999.0 .OR. t13col .GE. 9999999.0) t13col = 0.0
!   WHERE (t23col .LE. -999.0 .OR. t23col .GE. 9999999.0) t23col = 0.0
!   WHERE (h1col .LE. -999.0 .OR. h1col .GE. 9999999.0) h1col = 0.0
!   WHERE (h2col .LE. -999.0 .OR. h2col .GE. 9999999.0) h2col = 0.0
!   WHERE (h3col .LE. -999.0 .OR. h3col .GE. 9999999.0) h3col = 0.0
!   WHERE (q1col .LE. -999.0 .OR. q1col .GE. 9999999.0) q1col = 0.0
!   WHERE (q2col .LE. -999.0 .OR. q2col .GE. 9999999.0) q2col = 0.0
!   WHERE (q3col .LE. -999.0 .OR. q3col .GE. 9999999.0) q3col = 0.0 
!   WHERE (uscol .LE. -999.0 .OR. uscol .GE. 9999999.0) uscol = 0.0 
!   WHERE (z0col .LE. -999.0 .OR. z0col .GE. 9999999.0) z0col = 0.0

   RETURN

END SUBROUTINE compute_column_output

END MODULE module_column_out

