subroutine GET_FREQUENCY(ndata,dt,frq)
implicit none
integer:: ndata,nfrq,i
real*8 :: dt,frq(ndata),fNyq
fNyq = 0.5d0/dt
if(mod(ndata-1,2)==0) then
	nfrq=(ndata-1)/2
else
	nfrq=ndata/2
endif
frq(1)=0d0
do i=2,nfrq+1
frq(i) = (i-1)*fNyq/nfrq
enddo
do i=nfrq+2,ndata
frq(i) = frq(ndata+2-i)
enddo
return
end subroutine GET_FREQUENCY
